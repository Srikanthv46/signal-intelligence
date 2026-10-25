#!/usr/bin/env python3
"""
det_analog_vtx.py -- find analog FPV video transmitters (FM-modulated PAL / NTSC
composite video) anywhere in a wideband complex capture.

    detect(x, fs, fc, ctx) -> list[Emission]     family "analog_vtx", label "pal" | "ntsc" | "analog_fm_video"

Why this works far below the FM threshold
-----------------------------------------
An FM video carrier is constant-envelope, so |x| carries nothing, and below
~10 dB CNR a discriminator outputs mostly clicks. But the composite signal parks
the carrier at ONE frequency (the sync-tip frequency) for 4.7 us at the start of
every line, 15625 (PAL) / 15734.26 (NTSC) times a second, and nowhere near it
for the other ~93 % of the line (video swings the carrier away from it). In a
time-frequency picture with ~1 us x ~1 MHz cells, all of the signal's power sits
in a single ~1 MHz cell during the sync tip, i.e. the per-cell SNR there is the
in-band SNR times B/1 MHz (~+13 dB for a 20 MHz wide VTX). Folding that picture
modulo the line period over the whole capture (32 ms = ~500 lines) stacks those
cells: a VTX at -15 dB in-band SNR still stands out at >15 sigma.

  1. STFT, ~1 us Hann windows, 50 % overlap -> P[frame, bin]. Each bin is normalised by its
     own time-median and clipped at 8x -> bursty Wi-Fi / DJI frames cannot dominate a fold.
  2. For each line-rate hypothesis (PAL and NTSC, +-500 ppm, step small enough that the
     residual drift over the capture is < 0.3 us), fold the frames by line phase (128
     phase bins of 0.5 us) -> F[phase, bin]. Per bin, robust z-score over phase
     (median/MAD: a periodic emitter fills few phase bins, so MAD stays honest),
     then a 4.7 us box along phase (the sync-pulse matched filter) and 1-or-2 adjacent
     bins along frequency (sync tip between two bins). Score = max.
     Covers every carrier frequency at once: the 'channel scan step' is the STFT bin
     (fs/nfft, 0.96 MHz at 61.44 MSps), well inside the 4-10 MHz deviation swing.
  3. Score >= THRESHOLD -> emitter. Its columns (contiguous band with line-locked
     structure near the same line phase) are masked and the scan repeated, so several
     VTXs on different channels are reported separately.
  4. Per emitter (refine_timing / line_locked_spectrum / demod_features):
     - sync-tip bin = cleanest single pulse in the (lightly clipped) fold;
     - line rate: the period that makes the folded sync-band envelope sharpest (coarse
       4 ppm then 0.25 ppm grid), then -- if the SNR allows -- regression of the
       block-to-block drift of the sync leading edge (leave-one-out templates) -> ~0.2 ppm
       at high SNR, ~2 ppm at -5 dB, grid-limited ~5-20 ppm at -10..-15 dB;
     - precise sync-tip frequency from the spectrum of the sync-tip cores of all lines;
     - vertical sync (broad-pulse runs, pulse-vs-serration contrast) -> field timing and
       lines per field when two field syncs fall in the capture;
     - LINE-LOCKED spectrum: mean over all frames minus mean over sync-core frames. Only
       what is locked to this VTX's line clock survives -- co-channel Wi-Fi / DJI cancel --
       giving centre (PSD centroid of a constant-envelope FM signal = mean instantaneous
       frequency = the AC-coupled VTX's carrier), occupied band, in-band SNR, and from the
       active-video phases the peak-white edge -> deviation;
     - above ~4 dB SNR an FM discriminator adds the averaged line waveform (sync / blank /
       white levels), colour burst (4.43 / 3.58 MHz) and audio subcarrier (5.5..7 MHz).

Why other periodic emitters do not alias in
  Wi-Fi beacons (102.4 ms), Herelink 5 ms frames, DJI 15 kHz-OFDM symbols (71.4 us) and
  TDD frames are not multiples of 64 / 63.56 us, so folding spreads them across phase.
  Exact-rational cases (1 ms subframes = 15.625 lines, Wi-Fi 4 us symbols = 1/16 line)
  fold onto 8 or 16 evenly spaced phases, never into ONE 4.7 us pulse per line, and the
  median/clip normalisation caps a burst's weight. THRESHOLD was set on noise and on the
  five real captures without a VTX, as recorded and time-reversed and with +3 dB noise
  (max score 8.5 on the busy Wi-Fi capture, 6.6 on noise; threshold 10).

Known shared-code issue worked around here: Channelizer.slice (spectrum_survey) maps
bins by raw FFT index, so bands at/below DC come out scrambled unless L is a multiple
of the output length -> slice_band() below maps by signed frequency.

Channel tables (MHz, CH1..CH8)
  5.8 GHz A/B/E/F/R and the less common D/U/O/L/H bands: Oscar Liang, "5.8GHz FPV
  Channels & Frequency Chart" (oscarliang.com/fpv-channels). LowRace (L) = RaceBand - 296 MHz
  (= Oscar Liang's "D" row; Betaflight vtxtable LOWRACE). 1.2/1.3 GHz: the 9-channel
  table 1080..1360 step 40 + 1258 (Oscar Liang 1.2 GHz guide; same as 8-ch Lawmate/Boscam
  1.2G). 2.4 GHz: 4-ch 2414/2432/2450/2468 (Fat Shark/Lawmate) and the 8-ch/intl set
  2370..2510 step 20 (getfpv / flitetest charts).
"""
from __future__ import annotations
import numpy as np
from scipy import sparse, ndimage, signal

from spectrum_intel.core import common as C
from spectrum_intel.core.emission import Emission
from spectrum_intel.core.spectral import find_spurs
# ============================================================== constants
NOMINAL = {"pal": dict(fl=15625.0, lines=625, fsc=4.43361875e6),
           "ntsc": dict(fl=15734.264, lines=525, fsc=3.579545e6)}
SYNC_S = 4.7e-6
THRESHOLD = 10.0       # detection score; 15 real-capture nulls (5 captures x orig/reversed/+3 dB noise) max 8.5, noise max ~6.6
PPM_MAX = 500.0        # line-rate search span: crystal-clocked cameras are within +-100 ppm
NB = 128               # phase bins per line (0.5 us)
CLIP = 8.0             # per-bin normalised power clip (burst suppression)
MAX_EMITTERS = 4
DEV_CORR_HZ = 1.0e6      # STFT smearing of the peak-white edge (calibrated, see test)
MAX_SAMPLES = 4_000_000  # detection uses at most this many samples (~65 ms at 61.44 MSps)

CHANNELS = {
    "A": [5865, 5845, 5825, 5805, 5785, 5765, 5745, 5725],
    "B": [5733, 5752, 5771, 5790, 5809, 5828, 5847, 5866],
    "E": [5705, 5685, 5665, 5645, 5885, 5905, 5925, 5945],
    "F": [5740, 5760, 5780, 5800, 5820, 5840, 5860, 5880],      # Fatshark / ImmersionRC
    "R": [5658, 5695, 5732, 5769, 5806, 5843, 5880, 5917],      # RaceBand
    "L": [5362, 5399, 5436, 5473, 5510, 5547, 5584, 5621],      # LowRace (RaceBand - 296; "D" / Diatone)
    "LB": [5333, 5373, 5413, 5453, 5493, 5533, 5573, 5613],     # Lowband, 40 MHz step (research: 01_analog_video.md)
    "U": [5325, 5348, 5366, 5384, 5402, 5420, 5438, 5456],
    "O": [5474, 5492, 5510, 5528, 5546, 5564, 5582, 5600],
    "H": [5653, 5693, 5733, 5773, 5813, 5853, 5893, 5933],
    "1G3": [1080, 1120, 1160, 1200, 1240, 1280, 1320, 1360, 1258],
    "2G4": [2414, 2432, 2450, 2468],                            # Fat Shark / Lawmate 4-ch
    "2G4X": [2410, 2430, 2450, 2470, 2490, 2510, 2390, 2370],   # 8-ch / international
    # 3.3 GHz analog VTX / VRX, 16 channels in two groups (PLAN G9; vendor spec of one 3.3 GHz VRX, file 14 section 2)
    "3G3A": [3330, 3350, 3370, 3390, 3410, 3430, 3450, 3470],
    "3G3B": [3170, 3190, 3210, 3230, 3250, 3270, 3290, 3310],
}
# commonest first: when two tables put a channel within tolerance, prefer these names
BAND_PRIORITY = ["R", "F", "A", "B", "E", "L", "1G3", "2G4", "2G4X", "H", "U", "O", "LB", "3G3A", "3G3B"]


def match_channel(f_mhz, tol_mhz=3.0):
    """All standard channels within tol of f, nearest first (ties by band priority)."""
    out = []
    for band in BAND_PRIORITY:
        for i, fc in enumerate(CHANNELS[band]):
            d = f_mhz - fc
            if abs(d) <= tol_mhz:
                name = f"{band}{i + 1}" if band != "1G3" or i < 8 else "1G3-1258"
                out.append((abs(d), BAND_PRIORITY.index(band), name, fc))
    out.sort()
    return [(n, fc) for _, _, n, fc in out]


# ============================================================== time-frequency fold
def tf_power(x, fs):
    """~1 us Hann STFT, 50 % overlap. Returns q (median-normalised, clipped power),
    raw P, frame times (s), bin offsets (Hz)."""
    nfft = int(2 ** round(np.log2(fs * 1e-6)))
    nfft = max(nfft, 8)
    hop = nfft // 2
    n = (len(x) - nfft) // hop + 1
    idx = np.arange(nfft)[None, :] + hop * np.arange(n)[:, None]
    win = np.hanning(nfft).astype(np.float32)
    P = np.abs(np.fft.fftshift(np.fft.fft(x[idx] * win, axis=1), axes=1)).astype(np.float32) ** 2
    t = (np.arange(n) * hop + nfft / 2) / fs
    f = (np.arange(nfft) - nfft // 2) / nfft * fs
    med = np.median(P, 0) + 1e-30
    q = np.minimum(P / med, CLIP)
    return q, P, t, f


def fold(q, t, rate, nb=NB):
    """Mean of q per line-phase bin -> F[nb, bins]. Empty bins (rate commensurate with
    the hop) get the column median."""
    b = (np.mod(t * rate, 1.0) * nb).astype(np.int64)
    M = sparse.csr_matrix((np.ones(len(t), np.float32), (b, np.arange(len(t)))), shape=(nb, len(t)))
    cnt = np.bincount(b, minlength=nb).astype(float)
    F = np.asarray(M @ q) / np.maximum(cnt, 1)[:, None]
    if (cnt == 0).any():
        F[cnt == 0] = np.median(F[cnt > 0], 0)
    return F, cnt


def zscore(F, q_sd, n_per_bin, nb_box):
    """Robust z over phase per frequency bin, then the sync-pulse box along phase and
    the best of 1 or 2 adjacent bins along frequency."""
    med = np.median(F, 0)
    mad = 1.4826 * np.median(np.abs(F - med), 0)
    # MAD of 128 values can be small by chance; never trust less than 0.7x the iid value
    z = (F - med) / np.maximum(mad, 0.7 * q_sd / np.sqrt(n_per_bin) + 1e-12)
    zb = ndimage.convolve1d(z, np.ones(nb_box) / np.sqrt(nb_box), axis=0, mode="wrap")
    z2 = np.full_like(zb, -np.inf)
    z2[:, :-1] = (zb[:, :-1] + zb[:, 1:]) / np.sqrt(2)
    return zb, np.maximum(zb, z2)


def scan(q, t, f, col_mask=None, ppm_max=PPM_MAX):
    """Best (score, std, ppm, phase_bin, bin) over all line-rate hypotheses; also the
    per-hypothesis maxima (for the ppm profile)."""
    D = t[-1] - t[0]
    step = 0.6e-6 / D * 1e6                     # residual drift < 0.3 us over the capture
    ppms = np.arange(-ppm_max, ppm_max + step / 2, step)
    q_sd = q.std(0)
    best = None
    prof = {}
    for std, nom in NOMINAL.items():
        nb_box = int(round(SYNC_S * nom["fl"] * NB))
        prof[std] = []
        for p in ppms:
            r = nom["fl"] * (1 + p * 1e-6)
            F, cnt = fold(q, t, r)
            zb, zs = zscore(F, q_sd, cnt.mean(), nb_box)
            if col_mask is not None:
                zs[:, col_mask] = -np.inf
            i = np.unravel_index(np.argmax(zs), zs.shape)
            s = float(zs[i])
            prof[std].append(s)
            if best is None or s > best["score"]:
                best = dict(score=s, std=std, ppm=float(p), rate=r, b0=int(i[0]), k0=int(i[1]), zb=zb, F=F)
    best["ppm_profile"] = {k: (ppms, np.array(v)) for k, v in prof.items()}
    return best


def emitter_columns(zb, b0, k0, f, nbins_win, thr=3.5):
    """Frequency bins belonging to the emitter found at (b0, k0): contiguous run around
    k0 whose fold shows line-locked structure (|z| peak) within the blanking window."""
    nb = zb.shape[0]
    ph = (np.arange(-nbins_win, 2 * nbins_win) + b0) % nb
    s = np.max(np.abs(zb[ph]), 0)
    lo = hi = k0
    gap = 0
    while lo > 0 and gap < 3:
        lo -= 1; gap = 0 if s[lo] > thr else gap + 1
    lo += gap
    gap = 0
    while hi < len(s) - 1 and gap < 3:
        hi += 1; gap = 0 if s[hi] > thr else gap + 1
    hi -= gap
    return lo, hi


# ============================================================== refinement
def _xcorr_shift(a, b):
    """Circular shift (fractional bins) that best aligns profile a to template b."""
    A = np.fft.fft(a - a.mean()); B = np.fft.fft(b - b.mean())
    c = np.real(np.fft.ifft(A * np.conj(B)))
    k = int(np.argmax(c)); n = len(c)
    y0, y1, y2 = c[(k - 1) % n], c[k], c[(k + 1) % n]
    d = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2 + 1e-30)
    s = k + d
    return (s + n / 2) % n - n / 2, c[k] / (np.linalg.norm(a - a.mean()) * np.linalg.norm(b - b.mean()) + 1e-30)


def _lead_edge(p, smooth=1):
    """Positive part of the circular derivative of a folded sync-band profile. Every
    line type -- normal, equalising, broad pulse -- starts its sync at the same instant,
    but they END at different times, so aligning on leading edges keeps the field-
    blanking lines from biasing the block timing."""
    if smooth > 1:
        p = ndimage.uniform_filter1d(p, smooth, mode="wrap")
    d = p - np.roll(p, 1)
    return np.maximum(ndimage.uniform_filter1d(d, 3, mode="wrap"), 0)


def _fold_prof(b, e, nbf):
    return np.bincount(b, e, nbf) / np.maximum(np.bincount(b, minlength=nbf), 1)


def _sync_start(prof, rate, nbf):
    """Time (s, modulo the line) of the sync leading edge: start of the 4.7 us run with
    the most sync-band power in the folded envelope."""
    w = max(int(round(SYNC_S * rate * nbf)), 1)
    sm = ndimage.convolve1d(prof, np.ones(w), mode="wrap")
    c = int(np.argmax(sm))
    return ((c - w / 2) % nbf) / nbf / rate


def grid_line_rate(e, efs, rate, span_ppm=160.0):
    """Line rate that makes the folded sync leading edge sharpest (energy of the
    leading-edge profile): coarse 4 ppm / 0.25 us bins, then fine 0.25 ppm / 62 ns bins.
    Robust at low SNR (it integrates the whole capture), quantised at ~0.3 ppm."""
    tt = np.arange(len(e)) / efs
    def sharp(r, nb):
        # energy of the folded profile = sum of |line-rate harmonics|^2 (classic period
        # estimator); a short box tames the per-bin noise of the fine fold
        p = _fold_prof((np.mod(tt * r, 1.0) * nb).astype(np.int64), e, nb)
        p = ndimage.uniform_filter1d(p, max(nb // 256, 1), mode="wrap")
        return np.sum((p - p.mean()) ** 2)
    for span, step, nb in ((span_ppm, 4.0, 256), (6.0, 0.25, 1024)):
        ppms = np.arange(-span, span + step / 2, step)
        met = np.array([sharp(rate * (1 + p * 1e-6), nb) for p in ppms])
        if step < 1:
            met = ndimage.uniform_filter1d(met, 3)            # bin-quantisation ripple
        k = int(np.argmax(met)); d = 0.0
        if 0 < k < len(met) - 1:
            y0, y1, y2 = met[k - 1], met[k], met[k + 1]
            d = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2 + 1e-30)
        rate = rate * (1 + (ppms[k] + d * step) * 1e-6)
    return rate


def drift_line_rate(e, efs, rate, nblk, nbf=1024, iters=2, smooth=8):
    """Refine the line rate by regressing the sync position of each block of the capture
    (leading-edge cross-correlation against a leave-one-out template, so a block cannot
    lock onto its own noise) against time. Returns (rate, sigma_rate)."""
    n = len(e); tt = np.arange(n) / efs
    edges = np.linspace(0, n, nblk + 1).astype(int)
    sig = np.inf
    for _ in range(iters):
        b = (np.mod(tt * rate, 1.0) * nbf).astype(np.int64)
        S = np.bincount(b, e, nbf); N = np.bincount(b, minlength=nbf)
        sh, tc = [], []
        for a, z in zip(edges[:-1], edges[1:]):
            Sb = np.bincount(b[a:z], e[a:z], nbf); Nb = np.bincount(b[a:z], minlength=nbf)
            tm = _lead_edge((S - Sb) / np.maximum(N - Nb, 1), smooth)
            s_, _ = _xcorr_shift(_lead_edge(Sb / np.maximum(Nb, 1), smooth), tm)
            sh.append(s_ / nbf / rate); tc.append(tt[(a + z) // 2])
        sh, tc = np.array(sh), np.array(tc)
        A = np.vstack([tc - tc.mean(), np.ones_like(tc)]).T
        coef, *_ = np.linalg.lstsq(A, sh, rcond=None)
        res = sh - A @ coef
        sig_slope = np.sqrt(np.sum(res ** 2) / max(len(sh) - 2, 1) / np.sum((tc - tc.mean()) ** 2))
        # a leave-one-out template sits at minus the mean of the other blocks' offsets,
        # which inflates every measured offset by nblk/(nblk-1) -> undo that
        slope = coef[0] * (nblk - 1) / nblk
        sig = rate * sig_slope * (nblk - 1) / nblk
        rate = rate / (1 + slope)            # sync arriving later each block -> longer period
    return rate, sig


def sync_tip_freq(ch, fs, f_guess, rate, t0, half=2.0e6):
    """Precise sync-tip frequency: power spectrum of the middle 3.3 us of every sync tip
    (the carrier is a pure tone there), averaged over all lines."""
    y, ofs = slice_band(ch, f_guess, half, min(fs, 2.5 * half))
    T = 1 / rate
    nseg = int(3.3e-6 * ofs)
    starts = ((t0 + 0.7e-6 + np.arange(int(len(y) / ofs / T) - 1) * T) * ofs).astype(int)
    starts = starts[starts + nseg < len(y)]
    nfft = 512
    seg = y[starts[:, None] + np.arange(nseg)[None, :]] * np.hanning(nseg)
    p = np.mean(np.abs(np.fft.fftshift(np.fft.fft(seg, nfft, axis=1), axes=1)) ** 2, 0)
    fr = (np.arange(nfft) - nfft // 2) / nfft * ofs
    k = int(np.argmax(p)); d = 0.0
    if 0 < k < nfft - 1:
        y0, y1, y2 = np.log(p[k - 1:k + 2] + 1e-30)
        d = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2 + 1e-30)
    prom = p[k] / np.median(p)
    return f_guess + fr[k] + d * ofs / nfft, prom


def refine_timing(ch, fs, f_sync, rate, score, eclip=24.0):
    """Sync band envelope -> grid line rate -> sync phase -> precise sync-tip frequency
    -> narrower envelope -> drift regression. Returns dict."""
    def env(fq, hb):
        y, efs = slice_band(ch, fq, hb, min(fs, 8e6))
        e = np.abs(y) ** 2
        return np.minimum(e / np.median(e), eclip), efs
    e, efs = env(f_sync, 1.5e6)                             # tolerant of the ~0.5 MHz STFT error
    rate = grid_line_rate(e, efs, rate)
    t0 = _sync_start(_fold_prof((np.mod(np.arange(len(e)) / efs * rate, 1.0) * 1024).astype(np.int64), e, 1024), rate, 1024)
    f2, prom = sync_tip_freq(ch, fs, f_sync, rate, t0)
    if prom > 4 and abs(f2 - f_sync) < 1.5e6:
        f_sync = f2
        e, efs = env(f_sync, 1.0e6)
        rate = grid_line_rate(e, efs, rate, span_ppm=8.0)
    # few long blocks: the slope precision does not depend on the block count as long
    # as every block's sync is well above noise, and 4 blocks keep that true to ~-10 dB
    nblk = 8 if score > 100 else 4
    r2, sig = drift_line_rate(e, efs, rate, nblk, smooth=16)
    sig_ppm = sig / rate * 1e6
    if np.isfinite(sig) and sig_ppm < 3.0 and abs(r2 / rate - 1) * 1e6 < 8.0:
        rate = r2
    else:
        # regression failed (low SNR): keep the grid estimate; its spread on synthetic
        # VTXs at -10..-15 dB is ~5..20 ppm
        sig = rate * (10e-6 if score > 30 else 40e-6)
    tt = np.arange(len(e)) / efs
    prof = _fold_prof((np.mod(tt * rate, 1.0) * 1024).astype(np.int64), e, 1024)
    t0 = _sync_start(prof, rate, 1024)
    return dict(rate=rate, sig=sig, t0=t0, f_sync=f_sync, e=e, efs=efs, sync_prom=prom)


def field_timing(e, efs, rate, t0, std):
    """Vertical sync: broad pulses keep the carrier at the sync-tip frequency for ~27 us
    of every half line for 2.5 (PAL) / 3 (NTSC) lines. Score each half line by the mean
    sync-band power 6..26 us after its start (zero in every other kind of half line),
    find runs of 5 / 6, return their times and the half-line spacing between them."""
    T = 1 / rate
    nh = int((len(e) / efs - t0) / (T / 2))
    if nh < 20:
        return []
    hs = t0 + np.arange(nh) * T / 2
    a = (hs * efs).astype(int)
    on0, on1 = int(6e-6 * efs), int(26e-6 * efs)          # inside the broad pulse
    g0, g1 = int(28.2e-6 * efs), int(31.3e-6 * efs)       # inside its 4.7 us serration
    ok = a + g1 < len(e)
    a = a[ok]
    cs = np.concatenate([[0], np.cumsum(e)])
    hp = (cs[a + on1] - cs[a + on0]) / (on1 - on0) - (cs[a + g1] - cs[a + g0]) / (g1 - g0)
    # (a long co-channel burst fills pulse AND serration and cancels; only the broad-pulse
    #  comb gives pulse-high / gap-low half line after half line)
    run = 5 if std == "pal" else 6
    sm = np.convolve(hp, np.ones(run) / run, "same")
    med = np.median(sm); mad = 1.4826 * np.median(np.abs(sm - med)) + 1e-30
    zz = (sm - med) / mad
    pk = []
    for k in np.argsort(zz)[::-1]:
        if zz[k] < 8: break
        if all(abs(k - j) > 100 for j in pk):
            pk.append(int(k))
    if not pk:
        return []
    # keep VBIs whose spacing from the strongest is a whole number of fields (+-2 lines)
    lpf2 = NOMINAL[std]["lines"]                            # half lines per field
    k0 = pk[0]
    keep = [k for k in pk if abs(((k - k0) / lpf2) - round((k - k0) / lpf2)) * lpf2 <= 4]
    return sorted((k, float(zz[k]), t0 + k * T / 2) for k in keep)


def spectrum_features(P, f, k_sync, max_half_hz=25e6):
    """Occupied band, centre, in-band SNR and RMS deviation from the time-averaged PSD.
    Frames with a strong burst in the emitter's neighbourhood (power > median + 6 MAD)
    are dropped -- a constant-envelope FM carrier has constant in-band power, so this
    removes co-channel Wi-Fi/DJI bursts without biasing the VTX at any SNR. The floor is
    the receiver floor (25th pct per bin) interpolated across the emitter from the bins
    just outside it, so a sloping analog front-end response does not bias low-SNR
    estimates. For a constant-envelope signal the PSD centroid is the mean instantaneous
    frequency and the PSD second moment the mean-square frequency deviation."""
    nb = len(f); df = f[1] - f[0]
    W = int(np.ceil(max_half_hz / df))
    a0, b0 = max(k_sync - W, 0), min(k_sync + W, nb - 1)
    bp = P[:, a0:b0 + 1].sum(1)
    med = np.median(bp); mad = 1.4826 * np.median(np.abs(bp - med))
    keep = bp <= med + 6 * mad
    m = P[keep].mean(0)
    flo = np.percentile(P[keep], 25, axis=0) / 0.2877
    n0 = np.full(nb, np.percentile(flo, 20))
    sig = 1 / np.sqrt(keep.sum() / 2)                     # relative std of a bin mean (50 % overlap)
    lo = hi = k_sync
    for _ in range(2):
        ex = m / n0 - 1
        exs = ndimage.uniform_filter1d(ex, 3)
        thr = max(4 * sig / np.sqrt(3), 0.03 * np.percentile(exs[a0:b0 + 1], 90))
        on = exs > thr
        lo = hi = k_sync
        # grow across gaps of up to 3 bins (sync spike and video hump can be separated)
        while lo > a0 and on[max(lo - 4, a0):lo].any(): lo -= 1
        while hi < b0 and on[hi + 1:min(hi + 5, b0 + 1)].any(): hi += 1
        while lo < k_sync and not on[lo]: lo += 1
        while hi > k_sync and not on[hi]: hi -= 1
        # local floor: straight line through the receiver floor just outside the band
        side = [k for k in list(range(lo - 4, lo - 1)) + list(range(hi + 2, hi + 5)) if 0 <= k < nb]
        if len(side) >= 2:
            c = np.polyfit(f[side], flo[side], 1 if len(set(np.sign(np.array(side) - k_sync))) > 1 else 0)
            n0 = np.maximum(np.polyval(c, f), 0.5 * np.percentile(flo, 20))
    exa = np.maximum(m - n0, 0)[lo:hi + 1]
    fb = f[lo:hi + 1]
    tot = exa.sum() + 1e-30
    centre = float(np.sum(exa * fb) / tot)
    cdf = np.cumsum(exa) / tot
    e_lo = float(fb[min(np.searchsorted(cdf, 0.005), len(fb) - 1)] - df / 2)
    e_hi = float(fb[min(np.searchsorted(cdf, 0.995), len(fb) - 1)] + df / 2)
    rms = float(np.sqrt(max(np.sum(exa * (fb - centre) ** 2) / tot - df ** 2 / 12, 0)))
    snr_db = float(10 * np.log10(max(tot / n0[lo:hi + 1].sum(), 1e-6)))
    return dict(centre=centre, e_lo=e_lo, e_hi=e_hi, rms=rms, snr_db=snr_db, lo=lo, hi=hi,
                edge_cut=bool(lo <= 1 or hi >= nb - 2))


def line_locked_spectrum(P, t, f, rate, t0, k_sync, max_half_hz=25e6, clip=50.0):
    """VTX power spectrum with everything that is NOT locked to its line clock removed.
    During the core of the sync tip the VTX is in the sync bin only, so frames from that
    phase measure, in every other bin, noise + Wi-Fi + DJI + ... exactly as they are
    during the rest of the line (those emitters do not know the line phase). Mean over
    all frames minus mean over sync-core frames = the VTX's own spectrum, even when a
    busy Wi-Fi channel sits on top of it. Video-phase frames minus the same baseline =
    the distribution of the active-video carrier frequency -> peak-white edge."""
    T = 1 / rate
    nb = len(f); df = f[1] - f[0]
    ph = np.mod(t - t0, T)
    med = np.median(P, 0) + 1e-30
    Pn = np.minimum(P / med, clip)       # weak VTX: low clip, so Wi-Fi bursts cannot swamp it
    core = (ph > 1.2e-6) & (ph < 3.5e-6)
    act = (ph > 12.5e-6) & (ph < T - 2.5e-6)
    mid = (ph > 20e-6) & (ph < T - 10e-6)
    # baseline: sync-core frames, except in the bins the sync tip itself occupies (its
    # carrier spreads over ~2-4 bins), where mid-line frames are the VTX-free ones. The
    # minimum of the two picks the right one per bin without knowing the spread.
    b = np.minimum(Pn[core].mean(0), Pn[mid].mean(0))
    S = (Pn.mean(0) - b) * med
    Sv = (Pn[act].mean(0) - b) * med
    W = int(np.ceil(max_half_hz / df))
    a0, b0 = max(k_sync - W, 0), min(k_sync + W, nb - 1)
    rel = S / med
    # per-bin noise of S (frames ~half independent at 50 % overlap): bins under bursty
    # emitters get their own, larger, sigma instead of inflating everyone's
    nc = max(core.sum(), 2)
    sig = np.sqrt(Pn[core].var(0) / (nc / 2) + Pn.var(0) / (len(Pn) / 2)) + 1e-9
    rs = ndimage.uniform_filter1d(rel, 3)
    thr = np.maximum(4 * sig / np.sqrt(3), 0.03 * np.percentile(rs[a0:b0 + 1], 90))
    on = rs > thr
    lo = hi = k_sync
    while lo > a0 and on[max(lo - 4, a0):lo].any(): lo -= 1
    while hi < b0 and on[hi + 1:min(hi + 5, b0 + 1)].any(): hi += 1
    while lo < k_sync and not on[lo]: lo += 1
    while hi > k_sync and not on[hi]: hi -= 1
    if hi - lo < 3:
        # nothing passes the per-bin threshold (weak VTX under a busy Wi-Fi channel, whose
        # bins carry a large sigma): take the one-sided window on the video side of the
        # sync tip, where the VTX must be, and rely on the line-locked difference being
        # zero-mean for everything else
        side = np.sign(np.sum(S[a0:b0 + 1] * (f[a0:b0 + 1] - f[k_sync]))) or 1
        span = int(np.ceil(2.6 * 10e6 / df))                 # up to 2 x 13 MHz swing
        lo, hi = (max(k_sync - 2, 0), min(k_sync + span, nb - 1)) if side > 0 else \
                 (max(k_sync - span, 0), min(k_sync + 2, nb - 1))
        cs = np.cumsum(np.maximum(ndimage.uniform_filter1d(S[lo:hi + 1], 3), 0))
        cs /= cs[-1] + 1e-30
        a, z = int(np.searchsorted(cs, 0.01)), int(np.searchsorted(cs, 0.99))
        lo, hi = lo + a, lo + min(z, hi - lo)
    seg = np.maximum(S[lo:hi + 1], 0); fb = f[lo:hi + 1]
    tot = seg.sum() + 1e-30
    centre = float(np.sum(seg * fb) / tot)
    cdf = np.cumsum(seg) / tot
    e_lo = float(fb[min(np.searchsorted(cdf, 0.005), len(fb) - 1)] - df / 2)
    e_hi = float(fb[min(np.searchsorted(cdf, 0.995), len(fb) - 1)] + df / 2)
    rms = float(np.sqrt(max(np.sum(seg * (fb - centre) ** 2) / tot - df ** 2 / 12, 0)))
    # peak white: far edge (99.5 % from the sync side) of the active-video distribution
    v = np.maximum(Sv[lo:hi + 1], 0)
    if f[k_sync] > centre:
        v = v[::-1]; fv = fb[::-1]
    else:
        fv = fb
    cv = np.cumsum(v) / (v.sum() + 1e-30)
    f_white = float(fv[min(np.searchsorted(cv, 0.995), len(fv) - 1)])
    flo = np.percentile(P, 25, axis=0) / 0.2877
    n0 = np.percentile(flo, 20)
    snr_db = float(10 * np.log10(max(tot / (n0 * len(seg)), 1e-6)))
    return dict(centre=centre, e_lo=e_lo, e_hi=e_hi, rms=rms, snr_db=snr_db, lo=lo, hi=hi,
                f_white=f_white, edge_cut=bool(lo <= 1 or hi >= nb - 2))


def slice_band(ch, f_off, half_bw, out_fs):
    """Like Channelizer.slice, but maps bins by SIGNED frequency. (Channelizer.slice
    indexes the output by the raw FFT index, which is only right when L is a multiple
    of the output length or the band lies wholly at positive frequency; bands at or
    below DC come out scrambled.)"""
    L, fs = ch.L, ch.fs
    M = int(round(L * out_fs / fs))
    sel = np.where(np.abs(ch.f - f_off) <= half_bw)[0]
    kk = np.round(ch.f[sel] / fs * L).astype(np.int64)          # signed bin index
    k0 = int(round(f_off / fs * L))
    Y = np.zeros(M, complex)
    Y[(kk - k0) % M] = ch.X[sel]
    return np.fft.ifft(Y) * (M / L), M * fs / L          # actual output rate (M is rounded)


def fm_discriminate(ch, f_off, half_bw, fs):
    y, ofs = slice_band(ch, f_off, half_bw, min(fs, max(2.4 * half_bw, 4e6)))
    d = np.angle(y[1:] * np.conj(y[:-1])) * ofs / (2 * np.pi)
    return d, ofs


def demod_features(d, ofs, rate, t0, std, polarity):
    """Averaged line waveform (instantaneous frequency vs line phase), colour burst and
    audio subcarrier from the discriminator output. Meaningful above the FM threshold."""
    out = {}
    T = 1 / rate
    tt = np.arange(len(d)) / ofs
    ph = np.mod(tt - t0, T)
    nbf = 256
    b = (ph / T * nbf).astype(int)
    line = np.bincount(b, d, nbf) / np.maximum(np.bincount(b, minlength=nbf), 1)
    us = np.arange(nbf) / nbf * T * 1e6
    sync = np.median(line[(us > 1.0) & (us < 3.7)])
    blank = np.median(line[(us > 8.5) & (us < 10.0)] if std == "pal" else line[(us > 8.0) & (us < 9.2)])
    out["sync_level_mhz"] = sync / 1e6
    out["blank_level_mhz"] = blank / 1e6
    # peak white: 99.5th percentile of smoothed active-video frequency (polarity-aware)
    # peak white: brightest point of each line (0.3 us smoothing ~ luma bandwidth, kills
    # discriminator clicks), then the 90th percentile over lines -- robust to lines hit
    # by a co-channel burst, which only ever give a few wild maxima
    act = (ph > 12e-6) & (ph < T - 3e-6)
    sm = ndimage.uniform_filter1d(d, max(int(ofs * 0.3e-6), 1)) * polarity
    ln = np.floor((tt - t0) / T).astype(int)
    la, va = ln[act], sm[act]
    mx = np.full(la.max() - la.min() + 1, -np.inf)
    np.maximum.at(mx, la - la.min(), va)
    mx = mx[np.isfinite(mx)]
    white = polarity * np.percentile(mx, 90)
    out["white_level_mhz"] = white / 1e6
    out["deviation_demod_mhz"] = abs(white - sync) / 2e6
    # colour burst: gated back porch, energy at 4.43 / 3.58 MHz vs neighbours
    bw = (ph > 5.6e-6) & (ph < 7.6e-6)
    seg = d[bw]
    tb = tt[bw]
    def tone(fr):
        # per-line energy of the gated burst at fr, line-to-line phase ignored
        z = seg * np.exp(-2j * np.pi * fr * tb)
        ln = np.floor((tb - t0) / T).astype(int); ln -= ln.min()
        acc = np.bincount(ln, z.real) + 1j * np.bincount(ln, z.imag)
        return np.mean(np.abs(acc) ** 2)
    if ofs > 9.5e6:
        e_pal, e_ntsc = tone(4.43361875e6), tone(3.579545e6)
        ref = np.median([tone(fr) for fr in (2.6e6, 3.1e6, 4.0e6, 5.0e6, 5.4e6)])
        r_pal, r_ntsc = 10 * np.log10(e_pal / ref), 10 * np.log10(e_ntsc / ref)
        out["burst_db"] = dict(f4_43=round(r_pal, 1), f3_58=round(r_ntsc, 1))
        if max(r_pal, r_ntsc) > 6:
            out["color_burst_mhz"] = 4.434 if r_pal > r_ntsc else 3.580
    # colour standard: the burst demodulated line by line at each standard's subcarrier (a fixed multiple of
    # the measured line rate: the camera derives both from one crystal)
    if ofs > 9.5e6:
        # the sync time can sit on the trailing edge of H-sync (4.7 us late) at low SNR: try both edges and keep
        # the reading in which a colour standard passes, else the stronger coherence
        rs = [chroma_standard(d, ofs, rate, t0 + sh, std) for sh in (0.0, -SYNC_S)]
        rs = [r for r in rs if r]
        if rs:
            r = max(rs, key=lambda r: (bool(r.get("color")), r.get("chroma_coherence", 0.0)))
            if r is rs[-1] and len(rs) == 2 and r.get("color"):
                r["chroma_gate_note"] = "burst found 4.7 us before the measured sync time (trailing-edge lock)"
            out.update(r)
    # audio subcarrier: narrow FM bump at 5.5 / 6.0 / 6.5 / 7.0 MHz
    if ofs > 15e6:
        fr, p = signal.welch(d, ofs, nperseg=8192)
        cands = []
        for fa in (5.5e6, 6.0e6, 6.5e6, 7.0e6):
            if fa > ofs / 2 * 0.9: continue
            inn = np.abs(fr - fa) < 100e3
            ring = (np.abs(fr - fa) > 300e3) & (np.abs(fr - fa) < 700e3)
            cands.append((10 * np.log10(p[inn].max() / np.median(p[ring])), fa))
        if cands:
            prom, fa = max(cands)
            out["audio_prominence_db"] = round(prom, 1)
            if prom > 10:
                out["audio_subcarrier_mhz"] = fa / 1e6
                if prom >= AUDIO_PROM_DB:
                    out.update(audio_features(d, ofs, fa))
    return out


AUDIO_PROM_DB = 15.0              # audio content only above this subcarrier prominence (else video noise dominates)


def audio_features(d, ofs, fa, half=100e3, afs=200e3):
    """Audio on the FM subcarrier (research, PLAN A5): the subcarrier sliced +-100 kHz around fa (to ~200 kS/s),
    FM-discriminated -> audio (Hz of deviation). Returns the subcarrier's exact frequency, the RMS audio deviation,
    the dominant audio tone (50 Hz-8 kHz, with its prominence; on a drone the propeller blade-pass tone) and the
    spectral flatness (tone ~0, noise ~1)."""
    n = len(d)
    X = np.fft.fft(d - np.mean(d))
    k0 = int(round(fa / ofs * n))
    kk = int(half / ofs * n)
    M = 1 << int(np.ceil(np.log2(2 * kk + 1)))
    Y = np.zeros(M, complex)
    m = np.arange(-kk, kk + 1)
    Y[m % M] = X[(k0 + m) % n]
    fs2 = M * ofs / n
    z = np.fft.ifft(Y)
    a = np.angle(z[1:] * np.conj(z[:-1])) * fs2 / (2 * np.pi)          # Hz: carrier offset + audio
    off = float(np.median(a))
    a = a - np.mean(a)
    nps = min(len(a), 1 << int(np.floor(np.log2(len(a)))))
    f, P = signal.welch(a, fs2, nperseg=nps)
    band = (f >= 50) & (f <= 8000)
    if not band.any():
        return {}
    flat = float(np.exp(np.mean(np.log(P[band] + 1e-30))) / (np.mean(P[band]) + 1e-30))
    # the tone: whole-block Hann periodogram, 8x zero padding, parabolic peak (~29 Hz bins -> a few Hz)
    nf = 8 * (1 << int(np.ceil(np.log2(len(a)))))
    S2 = np.abs(np.fft.rfft(a * np.hanning(len(a)), nf)) ** 2
    f2 = np.arange(len(S2)) * fs2 / nf
    b2 = (f2 >= 50) & (f2 <= 8000)
    j = int(np.argmax(np.where(b2, S2, -1)))
    l0, l1, l2 = np.log(S2[j - 1] + 1e-30), np.log(S2[j] + 1e-30), np.log(S2[j + 1] + 1e-30)
    dj = 0.5 * (l0 - l2) / (l0 - 2 * l1 + l2) if (l0 - 2 * l1 + l2) < 0 else 0.0
    fpk = (j + dj) * fs2 / nf
    jp = int(np.argmin(np.abs(f - fpk)))
    lp = signal.sosfiltfilt(signal.butter(4, 15e3, fs=fs2, output="sos"), a)
    return dict(audio_subcarrier_hz=round(fa + off + (k0 / n * ofs - fa), 0),
                audio_dev_rms_khz=round(float(np.std(lp)) / 1e3, 2),
                audio_peak_hz=round(float(fpk), 1), audio_peak_db=round(float(10 * np.log10(P[jp] / np.median(P[band]))), 1),
                audio_flatness=round(flat, 3))


# ============================================================== colour standard (research: ITU-R BT.470 / SMPTE 170M)
# (name, fsc / fH, offset Hz at fH = 15625 Hz scaled with the line rate)
CHROMA = (("ntsc", 227.5, 0.0), ("pal_m", 227.25, 0.0), ("pal", 283.75, 25.0), ("pal_n", 229.25, 25.0))
CHROMA_WIN_S = (5.7e-6, 7.5e-6)   # burst gate after the sync leading edge (NTSC 5.3-7.8 us, PAL 5.6-7.85 us)
CHROMA_SEARCH = 0.06              # cycles/line searched around each candidate (+-940 Hz): line-rate error;
                                  # the nearest wrong standard (NTSC vs PAL-M) sits 0.25 cycles/line away
CHROMA_COH = 0.35                 # coherence of the line-to-line burst phasors for "colour present"
CHROMA_PROM_DB = 6.0              # burst energy at the subcarrier vs off-subcarrier references
CHROMA_REF_STEP, CHROMA_REF_GAP = 0.25e6, 1.0e6   # reference grid from 1.5 MHz to 0.45 ofs, >= 1 MHz from every fsc


SECAM_WIN_S = (5.8e-6, 10.0e-6)   # SECAM: rest-frequency subcarrier over the back porch (no burst)
SECAM_F = (4.25e6, 4.40625e6)     # Db / Dr rest frequencies (ITU-R BT.470)
SECAM_ALT = 0.25                  # line-to-line alternation of (e_Db - e_Dr) / (e_Db + e_Dr) (SECAM ~0.7, PAL ~0)
SECAM_T = 6.0                     # its significance (t statistic over the lines)


def secam_test(d, ofs, rate, t0):
    """SECAM (research, ITU-R BT.470): no burst; the back porch carries the FM subcarrier at its rest frequency,
    Db 4.250 MHz and Dr 4.40625 MHz on alternate lines. Per line, the energies at the two rest frequencies over a
    4.2 us gate (resolution ~240 kHz) give r = (e_Db - e_Dr) / (e_Db + e_Dr); SECAM makes r alternate line by line
    (|mean((-1)^n r)| ~ 0.7), a PAL burst (constant 4.43 MHz) or noise does not. The two rest frequencies are then
    measured per line parity (5 kHz grid): the camera's subcarrier generators."""
    T = 1.0 / rate
    k0 = int(np.ceil((-t0 - SECAM_WIN_S[0]) / T))
    k1 = int(np.floor((len(d) / ofs - t0 - SECAM_WIN_S[1]) / T))
    m = int((SECAM_WIN_S[1] - SECAM_WIN_S[0]) * ofs)
    if k1 - k0 + 1 < 32 or m < 8:
        return {}
    i0 = np.ceil((t0 + np.arange(k0, k1 + 1) * T + SECAM_WIN_S[0]) * ofs).astype(np.int64)
    keep = (i0 >= 0) & (i0 + m <= len(d))
    i0, lines = i0[keep], np.arange(k0, k1 + 1)[keep]
    if len(i0) < 32:
        return {}
    gi = i0[:, None] + np.arange(m)[None, :]
    seg = d[gi]
    seg = seg - seg.mean(axis=1, keepdims=True)
    tn = np.arange(m) / ofs                                  # |phasor| per line: the line's start phase drops out

    def energies(frs):
        ph = 2 * np.pi * np.outer(tn, frs)                   # real products (dgemm): lines x len(frs)
        return (seg @ np.cos(ph)) ** 2 + (seg @ np.sin(ph)) ** 2
    REF = (2.6e6, 3.1e6, 5.0e6, 5.4e6)
    E0 = energies(SECAM_F + REF)
    eb, er = E0[:, 0], E0[:, 1]
    r = (eb - er) / (eb + er + 1e-30)
    alt = np.where(lines % 2 == 0, 1.0, -1.0)
    a = r * alt
    mu = float(np.mean(a))
    tstat = abs(mu) / (float(np.std(a)) / np.sqrt(len(a)) + 1e-12)
    ref = float(np.median(E0[:, 2:].mean(axis=0))) + 1e-30
    prom = 10 * np.log10(float(np.mean(eb + er)) / 2 / ref)
    out = dict(alternation=round(abs(mu), 3), t_stat=round(float(tstat), 1), prominence_db=round(float(prom), 1))
    if abs(mu) >= SECAM_ALT:
        grid = np.arange(4.10e6, 4.56e6 + 1, 5e3)
        par_b = (alt > 0) if mu > 0 else (alt < 0)          # the lines where Db (lower) dominates
        E = energies(grid)
        out["db_hz"] = float(grid[int(np.argmax(E[par_b].mean(axis=0)))])
        out["dr_hz"] = float(grid[int(np.argmax(E[~par_b].mean(axis=0)))])
    return out


def chroma_standard(d, ofs, rate, t0, std):
    """Colour burst -> colour present, colour standard (NTSC / PAL-M on 525 lines, PAL / PAL-N on 625), subcarrier
    and burst phase swing. Per line, the gated burst is demodulated at each candidate fsc of the line family.
    The phasors z_n are
    - constant over lines for NTSC (180 deg burst);
    - alternating +-45 deg for the PAL family (V switch): mean and (-1)^n-mean both ~0.71 of |z|;
    - rotating at a wrong candidate.
    A residual rotation from the line-rate error is searched (|nu| <= CHROMA_SEARCH cycles/line) by an FFT over
    lines; coherence = (|Z(nu)|^2 + |Z(nu + 1/2)|^2) / (N sum |z|^2): 1 for a clean burst, ~2/N for noise.
    The line-to-line test only sees fsc modulo fH, and its score is symmetric under a half-line shift, so:
    - the candidates are restricted to the measured line family;
    - PAL (4.43 MHz) vs PAL-N (3.58 MHz) is decided by the burst energy at each candidate (the 1.8 us gate
      resolves ~0.6 MHz).
    Without colour (black-and-white camera, OSD-only field) only the line standard is known."""
    T = 1.0 / rate
    # the burst gate of every line, indexed directly (lines x gate samples)
    k0 = int(np.ceil((-t0 - CHROMA_WIN_S[0]) / T))
    k1 = int(np.floor((len(d) / ofs - t0 - CHROMA_WIN_S[1]) / T))
    m = int((CHROMA_WIN_S[1] - CHROMA_WIN_S[0]) * ofs)
    if k1 - k0 + 1 < 32 or m < 4:
        return {}
    i0 = np.ceil((t0 + np.arange(k0, k1 + 1) * T + CHROMA_WIN_S[0]) * ofs).astype(np.int64)
    i0 = i0[(i0 >= 0) & (i0 + m <= len(d))]
    if len(i0) < 32:
        return {}
    gi = i0[:, None] + np.arange(m)[None, :]
    seg, tb = d[gi], gi / ofs
    # per line, remove a straight line (the porch level and the tail of the sync edge): in a 1.8 us gate its
    # leakage otherwise swamps every frequency (a real PAL burst read -5 dB below the references,
    # 20261022 F1 capture; the sync-shifted reading scored coherence 0.6-0.9 from the edge alone)
    u = np.arange(m) - (m - 1) / 2
    seg = seg - seg.mean(axis=1, keepdims=True) - np.outer((seg * u).sum(axis=1) / (u * u).sum(), u)
    nl = seg.shape[0]
    nfft = 1 << int(np.ceil(np.log2(nl * 8)))
    kmax = int(np.ceil(CHROMA_SEARCH * nfft))
    idx = np.r_[0:kmax + 1, nfft - kmax:nfft]
    alt = (-1.0) ** np.arange(nl)

    def phasors(fr):
        return np.sum(seg * np.exp(-2j * np.pi * fr * tb), axis=1)
    # references below 0.45 ofs (5.4 MHz aliased at ~10 MSps) and >= 1 MHz from every candidate subcarrier
    fscs = [k * rate + off * rate / 15625.0 for _, k, off in CHROMA]
    grid = [fr for fr in np.arange(1.5e6, 0.45 * ofs, CHROMA_REF_STEP)
            if min(abs(fr - f0) for f0 in fscs) >= CHROMA_REF_GAP] or [2.0e6]
    ref = float(np.median([np.mean(np.abs(phasors(fr)) ** 2) for fr in grid])) + 1e-30
    fam = ("ntsc", "pal_m") if std == "ntsc" else ("pal", "pal_n")
    cands = {}
    for name, k, off in CHROMA:
        if name not in fam:
            continue
        fsc = k * rate + off * rate / 15625.0
        zl = phasors(fsc)
        e = np.sum(np.abs(zl) ** 2) + 1e-30
        Z, Zh = np.fft.fft(zl, nfft), np.fft.fft(zl * alt, nfft)
        c = (np.abs(Z[idx]) ** 2 + np.abs(Zh[idx]) ** 2) / (nl * e)
        j = int(np.argmax(c))
        nu = idx[j] / nfft if idx[j] <= nfft // 2 else idx[j] / nfft - 1
        swing = float(np.degrees(np.arctan2(np.abs(Zh[idx[j]]), np.abs(Z[idx[j]]))))
        cands[name] = dict(coherence=round(float(c[j]), 3), prominence_db=round(float(10 * np.log10(e / nl / ref)), 1),
                           swing_deg=round(swing, 1), fsc_hz=round(fsc + nu * rate, 1))
    out = {"chroma_candidates": cands}
    ok = {n: v for n, v in cands.items() if v["coherence"] >= CHROMA_COH and v["prominence_db"] >= CHROMA_PROM_DB
          and ((v["swing_deg"] > 25.0) if n != "ntsc" else (v["swing_deg"] < 20.0))}
    if not ok and std != "ntsc":
        sec = secam_test(d, ofs, rate, t0)
        if sec:
            out["secam_test"] = sec
            if sec["alternation"] >= SECAM_ALT and sec["t_stat"] >= SECAM_T and sec["prominence_db"] >= CHROMA_PROM_DB:
                out.update(color=True, color_standard="secam", secam_db_hz=sec["db_hz"], secam_dr_hz=sec["dr_hz"])
                return out
    if ok:
        name = max(ok, key=lambda n: (ok[n]["prominence_db"], ok[n]["coherence"]))
        v = ok[name]
        out.update(color=True, color_standard=name, chroma_fsc_hz=v["fsc_hz"], burst_swing_deg=v["swing_deg"],
                   chroma_coherence=v["coherence"])
    else:
        out["color"] = False
    return out


# ============================================================== main entry
def detect(x, fs, fc, ctx=None):
    ctx = ctx or {}
    x = np.asarray(x[:MAX_SAMPLES], np.complex64)
    dur_us = len(x) / fs * 1e6
    q, P, t, f = tf_power(x, fs)
    ch = ctx.get("channelizer")
    if ch is None or getattr(ch, "L", None) != len(x):
        # notch CW spurs / LO leakage: a pure tone in the sync band would otherwise be
        # taken for the (equally pure) sync-tip carrier by sync_tip_freq
        notch = [(sp["offset_hz"], 50e3) for sp in find_spurs(x, fs, fc)]
        ch = C.Channelizer(x.astype(np.complex128), fs, notch)
    mask = np.zeros(len(f), bool)
    found = []
    for _ in range(MAX_EMITTERS):
        best = scan(q, t, f, mask if mask.any() else None)
        if best["score"] < THRESHOLD:
            break
        nbw = int(round(10e-6 * best["rate"] * NB))
        cols = emitter_columns(best["zb"], best["b0"], best["k0"], f, nbw)
        mask[max(cols[0] - 1, 0):cols[1] + 2] = True
        if _is_part_of(best, f, fc, found):
            continue
        em = _characterise(x, fs, fc, ch, q, P, t, f, best, cols, dur_us)
        if em is None or _capture_locked(em, fs, fc):
            continue
        # the whole occupied band (+ audio sidebands) belongs to this VTX
        df = f[1] - f[0]
        lo_k = int(np.searchsorted(f, (em.f_lo_mhz * 1e6 - fc) - 1e6))
        hi_k = int(np.searchsorted(f, (em.f_hi_mhz * 1e6 - fc) + 1e6))
        mask[max(lo_k, 0):hi_k + 1] = True
        # same line clock and same sync phase as an emitter already reported -> a piece
        # of that transmitter (audio sideband, blanking level), not a second VTX
        dup = [o for o in found if o.features["standard"] == em.features["standard"]
               and abs(o.features["line_rate_ppm"] - em.features["line_rate_ppm"]) < 5
               and abs(((o.features["sync_phase_us"] - em.features["sync_phase_us"]) + 32) % (1e6 / o.features["line_rate_hz"]) - 32) < 3]
        if dup:
            continue
        found.append(em)
    return found


LOCK_PPM = 2.0          # line rate within this of fs / N: generated by the capture chain
LOCK_DC_HZ = 100e3     # ... with the "sync tip" on the LO


def _capture_locked(em, fs, fc):
    """LO leakage modulated by a periodic process of the capture chain (live Pluto, 8 MSps: a 'PAL' line clock of
    exactly fs / 512 with its sync tip 31-43 kHz from the LO, no field sync): its line rate is an integer fraction of
    the sample rate, its 'sync tip' sits on the LO, and no vertical sync is seen. A real VTX meets all three only by
    coincidence (its crystal error would have to cancel the fs / N offset within 2 ppm)."""
    F = em.features
    r = F["line_rate_hz"]
    lock_ppm = (r * round(fs / r) / fs - 1) * 1e6
    return (abs(lock_ppm) < LOCK_PPM and abs(F["sync_tip_mhz"] * 1e6 - fc) < LOCK_DC_HZ
            and F.get("field_rate_source") != "measured_vbi")


def _is_part_of(best, f, fc, found):
    """A further scan peak that has the same line standard, a line rate within the scan's
    resolution, a line phase inside an already reported VTX's horizontal blanking and a
    frequency inside its band is a leftover piece of that VTX (blanking level, audio
    sideband, band edge), not a second transmitter."""
    T = 1 / best["rate"]
    ph_us = best["b0"] / NB * T * 1e6
    fk = (fc + f[best["k0"]]) / 1e6
    for o in found:
        F = o.features
        if F["standard"] != best["std"] or abs(F["line_rate_ppm"] - best["ppm"]) > 150:
            continue
        d = (ph_us - F["sync_phase_us"]) % (T * 1e6)
        fsm = (f[1] - f[0]) * len(f) / 1e6                  # capture bandwidth: VTX skirts wrap
        dsync = abs(((fk - F["sync_tip_mhz"]) + fsm / 2) % fsm - fsm / 2)
        near = (o.f_lo_mhz - 5 <= fk <= o.f_hi_mhz + 5) or dsync < 25
        if (d < 14 or d > T * 1e6 - 4) and near:
            return True
    return False


def _characterise(x, fs, fc, ch, q, P, t, f, best, cols, dur_us):
    std, rate = best["std"], best["rate"]
    df = f[1] - f[0]
    # sync-tip bin: the scan's best cell can be the blanking level or, at high SNR, an
    # audio sideband of the sync tip (clipping equalises them). In the UNCLIPPED fold the
    # sync tip is unambiguous: it holds all of the carrier power for a full 4.7 us.
    Pn = np.minimum(P / (np.median(P, 0) + 1e-30), 50.0)
    Fr, _ = fold(Pn, t, rate)
    nb_box = int(round(SYNC_S * rate * NB))
    Fb = ndimage.convolve1d(Fr, np.ones(nb_box) / nb_box, axis=0, mode="wrap")
    exc = Fb - np.median(Fr, 0)
    # per-bin robust z: strong bursty emitters (DJI TDD) also give big unclipped excess,
    # but with a big spread over phase; the sync tip is one clean pulse on a flat floor
    zz = exc / (1.4826 * np.median(np.abs(Fb - np.median(Fb, 0)), 0) + 1e-12)
    k_lo, k_hi = max(cols[0] - 3, 0), min(cols[1] + 3, len(f) - 1)
    # Candidates: bins with a strong line-locked pulse. Among them the sync tip is
    #  (a) one of the LOUDEST at its pulse phase -- the carrier sits wholly in that bin for
    #      4.7 us (audio / chroma sidebands of it are >= 20 dB weaker), and
    #  (b) the EARLIEST of the loud ones -- the blanking-level bin, equally loud during the
    #      back porch, peaks ~5 us after the sync tip.
    zmax = zz[:, k_lo:k_hi + 1].max(0)
    cand = [k_lo + i for i in np.where(zmax >= 0.4 * zmax.max())[0]]
    bph = (np.mod(t * rate, 1.0) * NB).astype(int)
    half = max(nb_box // 2 - 1, 1)
    def circ(d):
        return (d + NB // 2) % NB - NB // 2
    cinfo = []
    for k in cand:
        bk = int(np.argmax(exc[:, k]))
        ins = np.abs(circ(bph - bk)) <= half
        cinfo.append((float(np.median(P[ins, k])), k, bk))
    lmax, kref, bref = max(cinfo)
    loud = [(circ(bk - bref), k, bk) for lv, k, bk in cinfo
            if lv >= 0.5 * lmax and abs(circ(bk - bref)) <= int(8e-6 * rate * NB)]
    _, k0, b0 = min(loud)
    kk = np.arange(max(k0 - 1, 0), min(k0 + 2, len(f)))
    w = np.maximum(exc[b0, kk], 0)
    f_sync = float(np.sum(w * f[kk]) / (w.sum() + 1e-12)) if w.sum() > 0 else float(f[k0])
    # VTX sync-tip level in units of the bin's median power (median over sync frames:
    # robust to bursts). Sets the clip levels below: high enough not to flatten the VTX,
    # low enough that Wi-Fi / DJI bursts cannot dominate the averages.
    insync = np.abs(((bph - b0) + NB // 2) % NB - NB // 2) <= max(nb_box // 2 - 1, 1)
    lvl = float(np.median(P[insync, k0]) / (np.median(P[:, k0]) + 1e-30))
    clip_v = max(6.0, 2.0 * lvl)
    # line rate, sync phase, precise sync-tip frequency from the sync-band envelope
    tm = refine_timing(ch, fs, f_sync, rate, best["score"], clip_v)
    rate, sig_rate, t0, f_sync, e, efs = (tm[k] for k in ("rate", "sig", "t0", "f_sync", "e", "efs"))
    nom = NOMINAL[std]
    ppm = (rate / nom["fl"] - 1) * 1e6
    # spectrum: centre, band, SNR
    k_s = int(np.argmin(np.abs(f - f_sync)))
    sp = line_locked_spectrum(P, t, f, rate, t0, k_s, clip=clip_v)
    sp["source"] = "line_locked"
    if sp["hi"] - sp["lo"] < 3:
        # too weak for the line-locked spectrum (few sync-core frames): plain PSD excess,
        # good on a clean floor, biased by co-channel emitters
        sp2 = spectrum_features(P, f, k_s)
        sp2["f_white"] = sp2["e_hi"] if f_sync < sp2["centre"] else sp2["e_lo"]
        sp2["source"] = "psd_excess"
        sp = sp2
    spec_ok = sp["hi"] - sp["lo"] >= 3
    centre, e_lo, e_hi, rms_dev, snr_db, edge_cut = (sp[k] for k in ("centre", "e_lo", "e_hi", "rms", "snr_db", "edge_cut"))
    polarity = 1 if f_sync < centre else -1
    # peak deviation = half the sync-to-white swing; the 1 us STFT smears the white edge
    # outward by ~DEV_CORR_HZ (calibrated on synthetic VTXs)
    dev = max(abs(sp["f_white"] - f_sync) - DEV_CORR_HZ, 0) / 2
    feats = dict(
        standard=std, line_rate_hz=round(rate, 4), line_rate_ppm=round(ppm, 2),
        line_rate_ppm_sigma=round(sig_rate / nom["fl"] * 1e6, 3) if np.isfinite(sig_rate) else None,
        line_period_us=round(1e6 / rate, 6), sync_phase_us=round(t0 * 1e6, 3),
        sync_tip_mhz=round((fc + f_sync) / 1e6, 3), sync_offset_mhz=round((f_sync - centre) / 1e6, 3),
        polarity="sync_low" if polarity > 0 else "sync_high",
        rms_dev_mhz=round(rms_dev / 1e6, 3), deviation_mhz=round(dev / 1e6, 2),
        bw_mhz=round((e_hi - e_lo) / 1e6, 2), snr_db=round(snr_db, 1), score=round(best["score"], 1),
        band_edge_truncated=bool(edge_cut), scan_bin_mhz=round(df / 1e6, 3),
        ppm_search_best=round(best["ppm"], 1))
    # field timing
    vbi = field_timing(e, efs, rate, t0, std)
    feats["vbi_us"] = [round(v[2] * 1e6, 1) for v in vbi]
    lpf_nom = nom["lines"] / 2
    if len(vbi) >= 2:
        nf = max(round((vbi[1][0] - vbi[0][0]) / nom["lines"]), 1)
        lpf = (vbi[1][0] - vbi[0][0]) / 2 / nf
        feats["lines_per_field"] = lpf
        feats["scan"] = "interlaced" if abs(lpf - round(lpf)) > 0.25 else "progressive"
        feats["field_rate_source"] = "measured_vbi"
    else:
        lpf = lpf_nom
        feats["field_rate_source"] = "standard"
    feats["field_rate_hz"] = round(rate / lpf, 5)
    # FM discriminator features when above threshold
    if snr_db > 4:
        hb = max(abs(e_hi - centre), abs(centre - e_lo)) + 1.0e6
        hb = min(hb, fs / 2 - abs(centre))
        d, ofs = fm_discriminate(ch, centre, hb, fs)
        dm = demod_features(d, ofs, rate, t0, std, polarity)
        feats.update({k: (round(v, 3) if isinstance(v, float) else v) for k, v in dm.items()})
        if dm.get("color") is False and snr_db < 10:
            # no burst found, but below ~10 dB a burst is lost in FM-threshold noise: not measurable, not B/W
            dm["color"] = None
            dm["color_note"] = f"colour burst not measurable at {snr_db:.1f} dB SNR"
            feats["color"] = None
            feats["color_note"] = dm["color_note"]
        if snr_db > 6 and dm.get("deviation_demod_mhz"):
            feats["deviation_mhz"] = round(dm["deviation_demod_mhz"], 2)
            feats["deviation_source"] = "discriminator"
    feats.setdefault("deviation_source", "spectrum")
    feats.setdefault("audio_subcarrier_mhz", None)
    # label: line standard; analog_fm_video if the timing disagrees with it
    label = std
    if feats.get("lines_per_field") and abs(feats["lines_per_field"] - lpf_nom) > 1.0:
        label = "analog_fm_video"
    if abs(ppm) > 0.9 * PPM_MAX:
        label = "analog_fm_video"
    # colour standard (research): PAL-M rides 525-line timing, PAL-N 625-line; a colour standard of the other
    # line family is reported, not used
    cs = feats.get("color_standard")
    if cs:
        fam = "ntsc" if cs in ("ntsc", "pal_m") else "pal"
        if fam != std:
            feats["color_note"] = f"colour burst of {cs} on {std} line timing: colour standard not used"
            feats.pop("color_standard")
        elif label == std and cs in ("pal_m", "pal_n", "secam"):
            label = cs
    cen_mhz = (fc + centre) / 1e6
    cands = match_channel(cen_mhz, 3.0) if spec_ok else []
    feats["spectrum_source"] = sp["source"] if spec_ok else "unmeasurable"
    if not spec_ok:
        # band/centre not measurable (very low SNR under strong co-channel emitters):
        # report the sync tip as the only reliable frequency
        feats["deviation_mhz"] = None; feats["bw_mhz"] = None; feats["rms_dev_mhz"] = None
    if cands:
        feats["channel"], ch_mhz = cands[0]
        feats["channel_candidates"] = [c for c, _ in cands]
        if any(c.startswith("LB") for c, _ in cands):
            feats["channel_note"] = "LB = Lowband 40 MHz-step table (vendor tables, research)"
        feats["nominal_mhz"] = ch_mhz
        feats["cfo_hz"] = round((cen_mhz - ch_mhz) * 1e6, -3)
        feats["cfo_ppm"] = round(feats["cfo_hz"] / (ch_mhz * 1e6) * 1e6, 1)
    else:
        feats["channel"] = None
    feats["center_mhz"] = round(cen_mhz, 3)
    conf = float(np.clip(1 - THRESHOLD / (2 * best["score"]), 0.3, 0.99))
    return Emission("analog_vtx", label, 0.0, round(dur_us, 1), round((fc + e_lo) / 1e6, 3),
                    round((fc + e_hi) / 1e6, 3), round(cen_mhz, 3), round(snr_db, 1), conf, None, feats)


def video_frame(x, fs, fc, em, ncols=720, max_lines=1250, ch=None):
    """FM-discriminate a detected VTX and stack its active lines into a picture
    (lines in time order, both fields). Returns uint8 [lines, ncols] (0 = black),
    or None when the emission has no line timing. Same method as the test's
    demod_picture()."""
    F = em.features
    if not F.get("line_rate_hz") or F.get("sync_phase_us") is None:
        return None
    if ch is None:
        ch = C.Channelizer(np.asarray(x, np.complex128), fs)
    hb = max(em.f_hi_mhz - em.center_mhz, em.center_mhz - em.f_lo_mhz) * 1e6 + 1e6
    hb = min(hb, fs / 2 - abs(em.center_mhz * 1e6 - fc))
    d, ofs = fm_discriminate(ch, em.center_mhz * 1e6 - fc, hb, fs)
    T = 1 / F["line_rate_hz"]
    t0 = F["sync_phase_us"] * 1e-6
    act0, act1 = (10.5e-6, T - 1.65e-6) if F.get("standard") == "pal" else (9.4e-6, T - 1.5e-6)
    nl = min(int((len(d) / ofs - t0) / T) - 1, max_lines)
    if nl < 2:
        return None
    ts = t0 + np.arange(nl)[:, None] * T + np.linspace(act0, act1, ncols)[None, :]
    img = np.interp(ts * ofs, np.arange(len(d)), d)
    if F.get("polarity") != "sync_low":
        img = -img
    lo, hi = np.percentile(img, [1, 99])
    return (np.clip((img - lo) / max(hi - lo, 1e-12), 0, 1) * 255).astype(np.uint8)


if __name__ == "__main__":
    import sys, json
    x, fs, fc, meta = C.load_capture(sys.argv[1])
    for em in detect(x, fs, fc, {}):
        print(json.dumps(em.to_dict(), default=float, indent=1))
