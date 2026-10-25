#!/usr/bin/env python3
"""
analog_vtx_synth.py -- synthetic analog FPV video transmitter (FM-modulated
PAL / NTSC composite video) for testing det_analog_vtx.py.

There is no real analog-VTX capture yet, so this generates one from the
standards, as an FPV VTX (RTC6705-class: direct FM of the VCO by the camera's
composite output, audio on 6.0 / 6.5 MHz FM subcarriers) would emit it:

  composite  : ITU-R BT.470 line timing -- 4.7 us sync tip at -0.3 V, blanking 0 V,
               colour burst on the back porch (PAL 4.43361875 MHz 10 cycles with
               +-45 deg swinging phase; NTSC 3.579545 MHz 9 cycles), active video
               from a synthetic image (gradients, sky/ground, bars, blocks, noise
               texture + chroma), and the full field-blanking structure (pre/post
               equalising pulses, broad pulses, blank VBI lines, 2:1 interlace).
               The whole video clock -- line rate, field rate, subcarrier -- is
               scaled by (1 + clock_ppm*1e-6): one crystal drives a camera's timing.
  VTX        : AC-coupled input (the VCO's mean frequency is the channel
               frequency whatever the picture brightness), optional first-order
               pre-emphasis, optional FM audio subcarrier, FM with peak deviation
               `dev_hz` defined as half the sync-tip-to-peak-white frequency swing.
               Polarity (sync tip low or high in frequency) is a parameter: VTX
               vendors differ and the receiver only cares through its discriminator.

Deviation: analog FPV deviation is not standardised and RTC6705-class datasheets
give only "1 Vpp video in". Published SDR receivers default to ~5 MHz peak
(5G8atv-rf-hackrf-decoder --dev 5e6), and the ~30 MHz analog channel occupancy
quoted for 5.8 GHz FPV implies, by Carson's rule 2*(dev + ~6.5 MHz), dev ~ 8 MHz.
Tests therefore draw dev_hz from 4..10 MHz.

embed(): scale a clean VTX waveform so that its in-band SNR -- measured with
drone_common.measure_snr, the same yardstick as other detector evaluations --
equals a target, and add it to a real capture.
"""
from __future__ import annotations
import numpy as np
from scipy import signal

from spectrum_intel.core import common as C
# (line rate, lines/frame, fsc, sync, active start, active end-from-line-end, burst start, burst cycles,
#  setup (black) level, first/last blanked lines of each field)
STANDARDS = {
    "pal": dict(fl=15625.0, lines=625, fsc=4.43361875e6, sync=4.7e-6, act0=10.5e-6, fp=1.65e-6,
                burst0=5.6e-6, burst_cyc=10, burst_amp=0.15, black=0.0, broad=27.3e-6, eq=2.35e-6),
    "ntsc": dict(fl=15734.264, lines=525, fsc=3.579545e6, sync=4.7e-6, act0=9.4e-6, fp=1.5e-6,
                 burst0=5.3e-6, burst_cyc=9, burst_amp=0.143, black=0.054, broad=27.1e-6, eq=2.3e-6),
    # PAL-M (Brazil): 525/59.94 timing, PAL colour at 227.25 fH; PAL-N (Argentina): 625/50 timing, 229.25 fH + 25 Hz
    "pal_m": dict(fl=15734.264, lines=525, fsc=227.25 * 15734.264, sync=4.7e-6, act0=9.4e-6, fp=1.5e-6,
                  burst0=5.3e-6, burst_cyc=9, burst_amp=0.143, black=0.054, broad=27.1e-6, eq=2.3e-6),
    "pal_n": dict(fl=15625.0, lines=625, fsc=229.25 * 15625.0 + 25.0, sync=4.7e-6, act0=10.5e-6, fp=1.65e-6,
                  burst0=5.6e-6, burst_cyc=9, burst_amp=0.15, black=0.0, broad=27.3e-6, eq=2.35e-6),
    # SECAM (ITU-R BT.470): 625/50 timing, no burst; FM chroma, Db on one line (rest 4.250 MHz, 230 kHz / unit) and
    # Dr on the next (4.40625 MHz, 280 kHz / unit), unblanked from the back porch, HF pre-emphasis "bell" at 4.286 MHz
    "secam": dict(fl=15625.0, lines=625, fsc=4.286e6, sync=4.7e-6, act0=10.5e-6, fp=1.65e-6,
                  burst0=5.6e-6, burst_cyc=0, burst_amp=0.0, black=0.0, broad=27.3e-6, eq=2.35e-6),
}
SWING = ("pal", "pal_m", "pal_n")       # PAL colour: V switch, burst +-45 deg line by line
SECAM_F0 = (4.25e6, 4.40625e6)          # Db, Dr rest frequencies
SECAM_DEV = (230e3, 280e3)              # Hz per unit of Db / Dr
SECAM_AMP = 0.08                        # ~23 % of the luma range peak-to-peak at 4.286 MHz


# ============================================================== picture
def make_image(rng, h=576, w=720, kind="mixed"):
    """Luma (0..1) + chroma (U, V in -0.5..0.5) of a plausible FPV picture:
    sky/ground horizon with gradient, a few objects, texture -- so active-video
    lines differ from each other (a constant picture would make the FM spectrum
    a comb and flatter the detector)."""
    y = np.zeros((h, w)); u = np.zeros((h, w)); v = np.zeros((h, w))
    rr, cc = np.mgrid[0:h, 0:w] / np.array([h, w])[:, None, None]
    if kind in ("mixed", "fpv"):
        horizon = 0.35 + 0.3 * rng.random() + 0.08 * np.sin(2 * np.pi * (cc * rng.uniform(0.5, 2) + rng.random()))
        sky = rr < horizon
        y = np.where(sky, 0.75 - 0.3 * rr, 0.35 + 0.2 * rng.random() * cc)
        u = np.where(sky, 0.25, -0.1); v = np.where(sky, -0.1, 0.05)
        for _ in range(rng.integers(4, 12)):                     # objects / buildings / trees
            r0, c0 = rng.random(2); hh, ww = rng.uniform(0.03, 0.3, 2)
            m = (np.abs(rr - r0) < hh / 2) & (np.abs(cc - c0) < ww / 2)
            y[m] = rng.uniform(0.05, 1.0); u[m] = rng.uniform(-0.4, 0.4); v[m] = rng.uniform(-0.4, 0.4)
        tex = ndimage_smooth(rng.standard_normal((h, w)), rng.uniform(0.5, 3))
        y = y + 0.12 * tex / (tex.std() + 1e-9)
    if kind in ("mixed", "bars") and rng.random() < 0.3:          # colour bars on part of the picture
        bars = np.array([[1, 0, 0], [0.89, -0.44, 0.1], [0.7, 0.15, -0.62], [0.59, -0.29, -0.52],
                         [0.41, 0.29, 0.52], [0.3, -0.15, 0.62], [0.11, 0.44, -0.1], [0, 0, 0]])
        top = rr < rng.uniform(0.1, 0.5)
        k = np.minimum((cc * 8).astype(int), 7)
        y = np.where(top, bars[k, 0], y); u = np.where(top, bars[k, 1], u); v = np.where(top, bars[k, 2], v)
    if kind == "noise":
        y = rng.random((h, w)); u = rng.uniform(-.3, .3, (h, w)); v = rng.uniform(-.3, .3, (h, w))
    return np.clip(y, 0, 1), np.clip(u, -.5, .5), np.clip(v, -.5, .5)


def ndimage_uniform(a, n):
    from scipy import ndimage
    return ndimage.uniform_filter1d(a, n)


def ndimage_smooth(a, s):
    from scipy import ndimage
    return ndimage.gaussian_filter(a, s)


# ============================================================== composite
def composite(n, fs, std="pal", clock_ppm=0.0, rng=None, image=None, t0_lines=None, color=True,
              progressive=False):
    """Composite video (volts: sync -0.3, blank 0, white 0.7) sampled at fs, n samples.
    Time is converted to the VTX's own video-clock time (scaled by clock_ppm), so
    every timing -- line, field, burst -- carries the same clock error.
    t0_lines: frame position (in lines) at capture start; random if None.
    progressive: 'non-interlaced 240p/288p' camera -- every field has an integer
    number of lines (262/312), common on cheap CMOS FPV cameras."""
    rng = rng or np.random.default_rng()
    S = STANDARDS[std]
    T = 1 / S["fl"]; L = S["lines"]
    six25 = L == 625                                          # 625/50 field structure (PAL, PAL-N)
    if image is None:
        image = make_image(rng, h=L - 49 if six25 else L - 39)
    Y, U, V = image
    if t0_lines is None:
        t0_lines = rng.uniform(0, L)
    # video-clock time in units of lines since frame start
    u = t0_lines + np.arange(n) / fs * S["fl"] * (1 + clock_ppm * 1e-6)
    if progressive:                                           # 2 fields of (L-1)/2 lines, no half line
        Lf = (L - 1) // 2
        fieldpos = np.mod(u, Lf)
        line_in_field = np.floor(fieldpos).astype(int)
        field = (np.floor(u / Lf).astype(int)) % 2
        tau = (fieldpos - line_in_field) * T
        # map onto an interlaced-like half-line index with field 1's VBI layout
        h = 2 * line_in_field + (tau >= T / 2)
        hl_frame = h                                          # both fields use field-1 VBI pattern
        Lh = 2 * Lf
    else:
        fpos = np.mod(u, L)                                   # line position in frame
        line = np.floor(fpos).astype(int)                     # 0-based line number (line 1 = 0)
        tau = (fpos - line) * T                               # time since line sync (s)
        hl_frame = 2 * line + (tau >= T / 2)                  # half-line index within frame
        Lh = 2 * L
    th = np.where(tau >= T / 2, tau - T / 2, tau)             # time within half-line
    # ---- VBI pulse type per half-line (ITU-R BT.470 / SMPTE 170M)
    if six25:
        # field 1: pre-eq h=1245..1249, broad 0..4, post-eq 5..9 ; field 2: pre-eq 620..624,
        # broad 625..629, post-eq 630..634 (312.5 lines later -> 2:1 interlace)
        pre = [(1245, 1250), (620, 625)]; brd = [(0, 5), (625, 630)]; post = [(5, 10), (630, 635)]
        blank_lines = [(0, 22), (310, 335), (622, 625)]
    else:
        pre = [(0, 6), (525, 531)]; brd = [(6, 12), (531, 537)]; post = [(12, 18), (537, 543)]
        blank_lines = [(0, 21), (262, 284)]
    if progressive:                                           # VBI of field 1 only, period Lh
        pre = [(Lh - 5, Lh)] if six25 else [(0, 6)]
        brd = [(0, 5)] if six25 else [(6, 12)]
        post = [(5, 10)] if six25 else [(12, 18)]
        blank_lines = [(0, 22), ((Lh // 2) - 1, Lh // 2)] if six25 else [(0, 21)]
    def inset(hh, rngs):
        m = np.zeros(hh.shape, bool)
        for a, b in rngs:
            m |= (hh >= a) & (hh < b)
        return m
    is_eq = inset(hl_frame, pre + post)
    is_brd = inset(hl_frame, brd)
    is_vbi = is_eq | is_brd
    line_no = hl_frame // 2
    blank = inset(line_no, blank_lines) | is_vbi
    v = np.zeros(n)
    # sync pulses
    norm_sync = (~is_vbi) & (tau < S["sync"])
    eq_sync = is_eq & (th < S["eq"])
    br_sync = is_brd & (th < S["broad"])
    v[norm_sync | eq_sync | br_sync] = -0.3
    # active video
    act = (~blank) & (tau >= S["act0"]) & (tau < T - S["fp"])
    if act.any():
        if progressive:
            row = line_no[act] - 22
        else:
            # interlace: field 1 lines -> even rows, field 2 -> odd rows
            f2 = line_no[act] >= (L // 2)
            a1, a2 = (22, 335) if six25 else (21, 284)       # first active line of each field
            fl_line = np.where(f2, line_no[act] - a2, line_no[act] - a1)
            row = 2 * np.maximum(fl_line, 0) + f2
        row = np.clip(row, 0, Y.shape[0] - 1)
        col = np.clip(((tau[act] - S["act0"]) / (T - S["fp"] - S["act0"]) * Y.shape[1]).astype(int), 0, Y.shape[1] - 1)
        luma = S["black"] + (0.7 - S["black"]) * Y[row, col]
        luma += 0.01 * rng.standard_normal(luma.size)          # camera noise
        v[act] = luma
    # chroma + burst (subcarrier phase-continuous in video-clock time)
    if color and std == "secam":
        ln = np.floor(u).astype(int)
        dr_line = (ln % 2) == 1                                # line-sequential Db / Dr
        on = (~is_vbi) & (tau >= 5.6e-6) & (tau < T - S["fp"])
        D = np.zeros(n)
        if act.any():
            db = np.clip(3.06 * U[row, col], -1.2, 1.2)        # Db = 1.505 (B-Y), U = 0.492 (B-Y)
            dr = np.clip(-2.17 * V[row, col], -1.2, 1.2)       # Dr = -1.902 (R-Y), V = 0.877 (R-Y)
            D[act] = np.where(dr_line[act], dr, db)
        f = np.where(dr_line, SECAM_F0[1] + SECAM_DEV[1] * D, SECAM_F0[0] + SECAM_DEV[0] * D)
        ph = 2 * np.pi * np.cumsum(f) / fs * (1 + clock_ppm * 1e-6)
        F = f / 4.286e6 - 4.286e6 / f
        bell = np.sqrt((1 + (16 * F) ** 2) / (1 + (1.26 * F) ** 2))
        v[on] += SECAM_AMP * bell[on] * np.cos(ph[on])
    elif color:
        tv = u * T                                             # video-clock seconds
        ph = 2 * np.pi * S["fsc"] * tv
        ln = np.floor(u).astype(int)
        if std in SWING:
            sw = np.where(ln % 2 == 0, 1.0, -1.0)              # V switch
            burst_ph = np.pi + sw * np.pi / 4
        else:
            sw = np.ones(n); burst_ph = np.pi * np.ones(n)
        bmask = (~is_vbi) & (tau >= S["burst0"]) & (tau < S["burst0"] + S["burst_cyc"] / S["fsc"])
        bmask &= ~inset(line_no, [(0, 9), (L // 2 - 3, L // 2 + 6)])  # no burst in field sync lines
        v[bmask] += S["burst_amp"] * np.sin(ph[bmask] + burst_ph[bmask])
        if act.any():
            c = U[row, col] * np.sin(ph[act]) + sw[act] * V[row, col] * np.cos(ph[act])
            v[act] += 0.6 * c
    # band-limit (5.5 MHz video bandwidth) -> realistic ~150 ns edges
    taps = signal.firwin(129, 5.5e6, fs=fs)
    v = signal.fftconvolve(v, taps, "same")
    return v


# ============================================================== FM VTX
def preemphasis(v, fs, f1=0.5e6, f2=2.5e6):
    """First-order lead network (zero f1, pole f2), unity DC gain: raises edges
    and chroma like broadcast FM-video pre-emphasis (CCIR 405 class)."""
    w1, w2 = 2 * np.pi * f1, 2 * np.pi * f2
    b, a = signal.bilinear([1 / w1, 1], [1 / w2, 1], fs)
    return signal.lfilter(b, a, v)


def synth_vtx(n, fs, f_off_hz=0.0, std="pal", dev_hz=8e6, clock_ppm=0.0, polarity=+1,
              audio_hz=6.5e6, audio_dev_hz=50e3, audio_dbc=-25.0, preemph=False, rng=None,
              image_kind="mixed", progressive=False, return_info=False, color=True, audio_tone_hz=None,
              t0_lines=None):
    """Complex baseband FM analog video at offset f_off_hz from the capture centre.
    Unit amplitude (constant envelope). dev_hz = half of the sync-tip-to-peak-white
    frequency swing (the composite is nominally 1 Vpp -> kf = 2*dev_hz per volt).
    audio_hz None -> no audio subcarrier. audio_dbc: level of the first audio
    sidebands vs carrier (RTC6705 datasheet VTAA: -25 dBc typ)."""
    rng = rng or np.random.default_rng()
    S = STANDARDS[std]
    img = make_image(rng, h=S["lines"] - 49 if S["lines"] == 625 else S["lines"] - 39, kind=image_kind)
    v = composite(n, fs, std, clock_ppm, rng, img, t0_lines=t0_lines, progressive=progressive, color=color)
    if preemph:
        v0pp = np.ptp(v)
        v = preemphasis(v, fs)
        v *= v0pp / np.ptp(v)                                   # keep the set swing
    kf = 2 * dev_hz / 1.0                                       # Hz per volt (1 Vpp composite)
    # deviation actually exercised by THIS picture (sync tip to its brightest 0.1 %):
    # the nominal dev_hz is only reached if the picture contains peak white
    vs = ndimage_uniform(v, max(int(fs * 0.3e-6), 1))
    dev_obs = kf * (np.percentile(vs, 99.9) - np.percentile(vs, 1.0)) / 2
    v = v - v.mean()                                            # AC coupling at the VTX input
    finst = f_off_hz + polarity * kf * v
    if audio_hz:
        beta = 2 * 10 ** (audio_dbc / 20)                       # J1/J0 ~ beta/2
        t = np.arange(n) / fs
        nz = signal.sosfilt(signal.butter(2, 8e3, fs=fs, output="sos"), rng.standard_normal(n))  # voice-band noise
        ft = audio_tone_hz if audio_tone_hz else rng.uniform(300, 3000)
        tone = 0.6 * np.sin(2 * np.pi * ft * t) + 0.4 * nz / (nz.std() + 1e-12)
        tone /= np.max(np.abs(tone)) + 1e-12
        sub_ph = 2 * np.pi * audio_hz * t + 2 * np.pi * audio_dev_hz * np.cumsum(tone) / fs
        finst = finst + beta * audio_hz * np.cos(sub_ph)
    phase = 2 * np.pi * np.cumsum(finst) / fs + rng.uniform(0, 2 * np.pi)
    x = np.exp(1j * phase).astype(np.complex64)
    if not return_info:
        return x
    info = dict(std=std, f_off_hz=f_off_hz, dev_hz=dev_hz, dev_obs_hz=float(dev_obs), clock_ppm=clock_ppm, polarity=polarity,
                audio_hz=audio_hz, line_rate_hz=S["fl"] * (1 + clock_ppm * 1e-6),
                field_rate_hz=S["fl"] * (1 + clock_ppm * 1e-6) / (S["lines"] / 2 if not progressive else (S["lines"] - 1) // 2),
                mean_hz=float(np.mean(finst)), kf=kf,
                sync_hz=float(np.percentile(finst, 0.5 if polarity > 0 else 99.5)),
                bw99_hz=occupied_bw(x, fs))
    return x, info


def occupied_bw(x, fs, frac=0.99, nfft=4096):
    """99 % power bandwidth (Hz) of a clean waveform."""
    f, p = signal.welch(x, fs, nperseg=nfft, return_onesided=False)
    o = np.argsort(f); f, p = f[o], p[o]
    c = np.cumsum(p) / p.sum()
    return float(f[np.searchsorted(c, 1 - (1 - frac) / 2)] - f[np.searchsorted(c, (1 - frac) / 2)])


def band_edges(x, fs, frac=0.99, nfft=4096):
    f, p = signal.welch(x, fs, nperseg=nfft, return_onesided=False)
    o = np.argsort(f); f, p = f[o], p[o]
    c = np.cumsum(p) / p.sum()
    return float(f[np.searchsorted(c, (1 - frac) / 2)]), float(f[np.searchsorted(c, 1 - (1 - frac) / 2)])


# ============================================================== embedding
GATE = 0.12   # fraction of the capture used for the gated measure_snr probe


def inband_snr(bg, sig, fs, fc, f_lo_mhz, f_hi_mhz, nfft=1024):
    """In-band SNR (dB) of `sig` over the floor of `bg`, with exactly C.measure_snr's
    definition -- mean spectrogram power in the box over the per-bin noise floor
    (C.noise_floor_psd) in the same bins -- but with the signal term taken from the clean
    waveform. Why not call measure_snr on bg+sig: its floor is a 25th-percentile-over-time
    estimate that a CONTINUOUS emitter raises itself (it would read ~-inf dB), and
    subtracting a busy background's own in-box excess (Wi-Fi, DJI) is unstable at low SNR."""
    Ps, t, f = C.spectrogram(sig, fs, fc, nfft)
    Ps /= np.sum(np.hanning(nfft) ** 2)
    fl = C.noise_floor_psd(bg, nfft)["floor_bin"]
    fi = (f >= f_lo_mhz) & (f <= f_hi_mhz)
    return float(10 * np.log10(Ps[:, fi].mean() / fl[fi].mean()))


def measure_snr_gated(y_bg, sig, fs, fc, f_lo_mhz, f_hi_mhz):
    """Cross-check with C.measure_snr itself: add sig over only the first GATE of the
    capture (so measure_snr's floor comes from signal-free frames) and measure that box.
    Valid on a quiet background; ~0.5 dB high-floor bias at high SNR."""
    g = int(len(y_bg) * GATE)
    probe = y_bg.copy(); probe[:g] += sig[:g]
    return C.measure_snr(probe, fs, fc, 0, g / fs * 1e6, f_lo_mhz, f_hi_mhz)


def embed(bg, sig, fs, fc, snr_db, f_lo_mhz, f_hi_mhz):
    """Scale sig so its in-band SNR in [f_lo, f_hi] over bg's floor is snr_db (inband_snr),
    return (bg + g*sig, measured snr, g)."""
    m0 = inband_snr(bg, sig, fs, fc, f_lo_mhz, f_hi_mhz)
    g = 10 ** ((snr_db - m0) / 20)
    m = inband_snr(bg, g * sig, fs, fc, f_lo_mhz, f_hi_mhz)
    return (bg + g * sig).astype(np.complex64), m, g


def awgn(n, rng, power=1.0):
    return ((rng.standard_normal(n) + 1j * rng.standard_normal(n)) * np.sqrt(power / 2)).astype(np.complex64)


if __name__ == "__main__":
    rng = np.random.default_rng(1)
    fs = 61.44e6
    x, info = synth_vtx(2_000_000, fs, 0.0, "pal", 8e6, 35.0, rng=rng, return_info=True)
    print(info)
