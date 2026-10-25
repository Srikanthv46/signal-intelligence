#!/usr/bin/env python3
"""
det_dji.py -- DJI emission detector (OcuSync video / RC uplink / DroneID).

    detect(x, fs, fc, ctx) -> list[Emission]          (drone_common interface)

What DJI links look like (measured on the bundled captures, see test_det_dji.py):

  OcuSync 2 / O3 video  15 kHz-SCS OFDM, LTE numerology: Tu 66.67 us, normal CP
                        4.69 us (5.21 us on the first symbol of each 0.5 ms slot),
                        symbol period 71.43 us averaged over the slot. 10 MHz mode:
                        600 used subcarriers (+ empty DC), 20 MHz mode: 1200.
                        Every 0.5 ms slot starts with a Zadoff-Chu symbol whose
                        length is used+1 (601 / 1201, DJI drops the middle ZC value
                        and maps the rest around DC) -- the DJI-specific part: an
                        LTE-like emitter that is not DJI (e.g. Herelink, also
                        15 kHz 20 MHz OFDM) has no such symbols. Bursts are TDD
                        and start on a 1 ms grid; the symbol grid is continuous
                        across bursts, so the CP correlation can be folded over
                        the WHOLE capture (hundreds of symbols) at one centre.
                        Not every link: a 10 MHz link at 2466.5 MHz (2484 MHz
                        capture, test_signal-20351019) moves its grid ~10-20
                        samples per burst and sends one ZC symbol (root 116) at
                        the start of each ~5 ms burst, 1225 samples before the
                        next symbol -> fallback: re-time burst by burst and test
                        the strongest individual symbols (_retimed_link).
  O4 video (Mini 4 Pro) 30 kHz-SCS OFDM: Tu 33.33 us, CP 2.34 us, symbol 35.71 us;
                        ~1.3-1.8 ms bursts every 2.5 ms. No ZC symbols found; DJI
                        evidence is a reference symbol repeated every 4 symbols
                        (correlation peak at a 142.9 us lag).
  RC uplink             15 kHz OFDM, 72 used subcarriers (1.08 MHz), 7-symbol
                        0.5 ms slots; symbols 0 and 6 are length-73 ZC (root is
                        link-specific: 50 and 6 seen). Hops every 4/6 ms.
  DroneID               15 kHz OFDM, 600 subcarriers (9 MHz), 9 symbols (CP 80/72
                        samples at 15.36 MSps), ZC roots 600 / 147 in symbols 4 and
                        6; decoded with droneid-decoder-reference/src/droneid_api.

Detectors are coherent (CP correlation folded over many symbols, ZC structure,
known DroneID CP schedule), not energy based, so they keep working below 0 dB
in-band SNR where a spectrogram shows nothing.

Fingerprint features (Emission.features): carrier_mhz, cfo_hz (vs the 100 kHz
raster, an assumption: DJI's channel raster is not public), cfo_ppm, symbol_us,
cp_us, scs_khz, bw_mhz, burst_us, n_symbols, frame_period_ms, frame_phase_us,
clock_ppm (link-level, from the symbol period over the capture), zc_root,
serial (DroneID), link_id (bursts sharing one continuous symbol grid).
"""
from __future__ import annotations
import os, sys
from collections import Counter
import numpy as np
from scipy import ndimage
from scipy.signal import find_peaks

from spectrum_intel.core import common as C
from spectrum_intel.core.emission import Emission
from spectrum_intel.core import spectral as SP, ofdm15k as O15, bursts as BU
from ._paths import droneid_src as _droneid_src   # external DroneID decoder, configurable
DRONEID_SRC = _droneid_src()

__all__ = ["detect", "VIDEO_PROFILES", "zc_root_search", "droneid_decode"]


# ============================================================== small DSP helpers
def _movsum(v, k):
    c = np.concatenate(([0], np.cumsum(v)))
    return c[k:] - c[:-k]


def _cp_metric(y, N, cp):
    """Normalised CP correlation (same as spectrum_survey._cp_metric, cumsum-based
    so long slices are cheap). Returns |corr|/energy and the complex corr sum."""
    p = y[:-N] * np.conj(y[N:])
    e = (np.abs(y[:-N]) ** 2 + np.abs(y[N:]) ** 2) / 2
    cs = _movsum(p, cp)
    return np.abs(cs) / (_movsum(e, cp) + 1e-20), cs


_FOLD_IDX = {}


def _fold(v, P):
    """Mean of v folded modulo a (fractional) period P samples."""
    key = (len(v), P)
    if key not in _FOLD_IDX:
        idx = np.floor(np.mod(np.arange(len(v)), P)).astype(np.int64)
        _FOLD_IDX[key] = (idx, np.maximum(np.bincount(idx, None, int(np.ceil(P))), 1))
    idx, cnt = _FOLD_IDX[key]
    return np.bincount(idx, v, len(cnt)) / cnt


class _FastChannelizer:
    """Same output as drone_common.Channelizer.slice (brick-wall FFT-domain
    channelizer) but selects bins by index arithmetic instead of a mask over the
    whole spectrum, and keeps a few recent slices: the scans below cut hundreds
    of slices from one 2 M-point FFT."""

    def __init__(self, ch, keep=6):
        self.fs, self.L, self.X = ch.fs, ch.L, ch.X
        self.keep, self.cache = keep, {}

    def slice(self, f_off, half_bw, out_fs):
        key = (round(f_off, 3), round(half_bw, 3), out_fs)
        if key in self.cache:
            return self.cache[key]
        L, fs = self.L, self.fs
        M = int(round(L * out_fs / fs))
        df = fs / L
        ks = np.arange(int(np.ceil((f_off - half_bw) / df - 1e-9)), int(np.floor((f_off + half_bw) / df + 1e-9)) + 1)
        k0 = int(round(f_off / fs * L))
        Y = np.zeros(M, complex)
        Y[(ks - k0) % M] = self.X[ks % L]
        y = np.fft.ifft(Y) * (M / L)
        if len(self.cache) >= self.keep:
            self.cache.pop(next(iter(self.cache)))
        self.cache[key] = y
        return y


def _robust_z(prof):
    med = np.median(prof)
    sig = 1.4826 * np.median(np.abs(prof - med)) + 1e-12
    return (prof - med) / sig


def zc_freq(root, L):
    """ZC as DJI maps it into the used subcarriers: L values, middle one dropped,
    the rest placed around an empty DC bin (droneid-decoder zcsequence_f)."""
    n = np.arange(L)
    z = np.exp(-1j * np.pi * root * n * (n + 1) / L)
    h = L // 2
    out = np.zeros(L, complex)
    out[:h], out[h + 1:] = z[:h], z[h + 1:]
    return out


def _zc_dd(X, L):
    """Second difference of the used subcarriers. For a ZC of root u,
    X[k+1]X*[k] advances by exp(-j2*pi*u/L) per subcarrier, so DD is CONSTANT
    whatever the timing offset (linear phase) and channel phase -> summing DD is
    a coherent ZC test (noise averages as 1/sqrt(L)) that also yields the root.
    X: fftshifted symbol spectrum with the transmitter DC in the centre."""
    N = len(X); h = (L - 1) // 2
    acc, mag = 0j, 0.0
    for seg in (X[N // 2 - h:N // 2], X[N // 2 + 1:N // 2 + 1 + h]):   # skip the DC bin
        D = seg[1:] * np.conj(seg[:-1])
        DD = D[1:] * np.conj(D[:-1])
        acc += DD.sum(); mag += np.abs(DD).sum()
    return acc, mag


def _root_from_dd(acc, L):
    return int(np.round(-np.angle(acc) * L / (2 * np.pi))) % L


def zc_root_search(X, L, roots=None):
    """Best ZC root by direct correlation over cyclic shifts (timing offset of a
    ZC symbol = cyclic shift of its spectrum). X: L used bins incl. DC. Returns
    (corr, root); corr ~1 for a clean ZC, ~3/sqrt(L) for noise."""
    roots = range(1, L) if roots is None else roots
    Xn = X / (np.linalg.norm(X) + 1e-20)
    F = np.fft.fft(Xn, 2 * L)
    best = (0.0, None)
    for r in roots:
        z = zc_freq(r, L)
        cc = np.abs(np.fft.ifft(F * np.conj(np.fft.fft(z, 2 * L)))).max() / np.linalg.norm(z)
        if cc > best[0]:
            best = (float(cc), int(r))
    return best


class _SNR:
    """One spectrogram + per-bin floor per capture; in-band SNR of any box
    (same definition as drone_common.measure_snr, computed once)."""

    def __init__(self, x, fs, fc, nfft=1024):
        P, self.t, self.f = C.spectrogram(x, fs, fc, nfft)
        self.P = P / np.sum(np.hanning(nfft) ** 2)
        # receiver (thermal) floor, flat across the band. The per-bin floor of
        # noise_floor_psd is a time percentile, which a link busy > 50 % of the time
        # (OcuSync video) lifts to its own level -> in-band SNR would read ~0 dB
        self.fl = np.full(len(self.f), C.noise_floor_psd(x, nfft)["floor_power_per_sample"])

    def __call__(self, t0, t1, f0, f1):
        ti = (self.t >= t0) & (self.t <= t1); fi = (self.f >= f0) & (self.f <= f1)
        if not ti.any() or not fi.any():
            return None
        sig = self.P[np.ix_(ti, fi)].mean(0) - self.fl[fi]
        return round(float(10 * np.log10(max(sig.mean(), 1e-12) / self.fl[fi].mean())), 1)


def _cfo_feats(carrier_hz, raster=100e3):
    nom = np.round(carrier_hz / raster) * raster
    cfo = carrier_hz - nom
    return dict(carrier_mhz=round(carrier_hz / 1e6, 5), cfo_hz=round(float(cfo), 1),
                cfo_ppm=round(float(cfo / carrier_hz * 1e6), 3))


# ============================================================== video (continuous symbol grid)
# (name, SCS kHz, output rate, FFT N, normal CP @rate, occupied BW, scan step, ZC length)
VIDEO_PROFILES = [
    ("15k_10MHz", 15, 15.36e6, 1024, 72, 9.0e6, 1.0e6, 601),
    ("15k_20MHz", 15, 30.72e6, 2048, 144, 18.0e6, 2.0e6, 1201),
    ("30k_10MHz", 30, 15.36e6, 512, 36, 9.0e6, 1.0e6, None),
    ("30k_20MHz", 30, 30.72e6, 1024, 72, 18.0e6, 2.0e6, None),
]
SLOT_S = 0.5e-3                      # 7 (15 kHz) or 14 (30 kHz) symbols
ZC_MATCH = 0.25                      # mean best-symbol DJI-ZC correlation needed (noise 0.10-0.13, Herelink <= 0.2)


def _scan_video(ch, usable_half, z_thr):
    """Capture-level search: fold the CP metric of each (profile, centre) modulo
    the nominal symbol period over the whole capture. A link whose symbol grid is
    continuous across its TDD bursts stacks every symbol into one bin.
    Profiles with the same output rate / width share one slice (15 and 30 kHz)."""
    groups = {}
    for prof in VIDEO_PROFILES:
        groups.setdefault((prof[2], prof[5], prof[6]), []).append(prof)
    cands = []
    for (ofs, occ, step), profs in groups.items():
        offs = np.arange(-usable_half + occ / 2, usable_half - occ / 2 + 1, step)
        zs = {p[0]: [] for p in profs}
        for off in offs:
            y = ch.slice(off, occ / 2 * 1.05, ofs)
            for name, scs, _, N, cp, _, _, L in profs:
                P = SLOT_S * ofs / (7 if scs == 15 else 14)
                z = _robust_z(_fold(_cp_metric(y, N, cp)[0], P))
                k = int(np.argmax(z))
                zs[name].append((float(z[k]), k))
        for prof in profs:
            zv = np.array([z for z, _ in zs[prof[0]]])
            P = SLOT_S * ofs / (7 if prof[1] == 15 else 14)
            for i in range(len(offs)):                    # local maxima over centre
                if zv[i] < z_thr or zv[i] < zv[max(i - 1, 0)] or zv[i] < zv[min(i + 1, len(zv) - 1)]:
                    continue
                cands.append(dict(profile=prof, f_off=float(offs[i]), z=zv[i], phase=zs[prof[0]][i][1], P=P))
    return cands


def _segment(m, phase, P, n_avg=4, k_sigma=4.0):
    """Burst segmentation on a KNOWN symbol grid: the CP metric sampled at the
    grid phase vs half a symbol away. Knowing the phase is what makes this far
    more sensitive than a blind per-window search."""
    K = int((len(m) - phase - P) // P)
    t = phase + P * np.arange(K)
    samp = lambda dt: np.array([m[max(int(round(a + dt)) - 2, 0):int(round(a + dt)) + 3].max() for a in t])
    on, off = samp(0), samp(P / 2)
    # noise statistics from two OFF-grid phases (pure noise in every symbol, whatever
    # the link duty cycle -- on-off itself is half signal for a busy link)
    nref = ndimage.uniform_filter1d(off - samp(P / 4), n_avg)
    sig = 1.4826 * np.median(np.abs(nref - np.median(nref))) + 1e-9
    d = ndimage.uniform_filter1d(on - off, n_avg)
    mask = d > k_sigma * sig
    mask = ndimage.binary_closing(mask, np.ones(3))
    lab, n = ndimage.label(mask)
    runs = []
    for sl in ndimage.find_objects(lab):
        a, b = sl[0].start, sl[0].stop
        # the smoothing widens runs by ~n_avg/2 each side; trim to symbols that are on
        while a < b and on[a] - off[a] < sig: a += 1
        while b > a and on[b - 1] - off[b - 1] < sig: b -= 1
        if b - a >= 3:
            runs.append((a, b))
    return runs, t, on, off


def _segment_multi(ms, phase, P, n_avg=4, k_sigma=4.0, min_len=3, exclude=()):
    """_segment over several CP metrics with the same timing (full band + sub-bands): a
    symbol is on when ANY of them shows the CP structure, so a strong signal overlapping
    part of the band (Wi-Fi over a DJI burst) does not cut the burst. Edges are trimmed and
    extended at 3 per-symbol sigma on the max (the smoothed sigma is about half the
    per-symbol one; 1 sigma on the max of 5 metrics passed 6-11 % of noise symbols).
    exclude: sample intervals of another known emission (an in-band RC uplink hop), never
    part of a run. Returns runs, the grid, and the smoothed and per-symbol z of the first
    (full-band) metric (the per-symbol one -inf where excluded; ~N(0, 1) under noise, unlike
    the max over the metrics, whose noise mean is ~1.2)."""
    K = int((len(ms[0]) - phase - P) // P)
    t = phase + P * np.arange(K)
    zs, zr = [], []
    for mm in ms:
        samp = lambda dt: np.array([mm[max(int(round(a + dt)) - 2, 0):int(round(a + dt)) + 3].max() for a in t])
        on, off = samp(0), samp(P / 2)
        q = off - samp(P / 4)
        nref = ndimage.uniform_filter1d(q, n_avg)
        sig = 1.4826 * np.median(np.abs(nref - np.median(nref))) + 1e-9
        zs.append(ndimage.uniform_filter1d(on - off, n_avg) / sig)
        zr.append((on - off) / (1.4826 * np.median(np.abs(q - np.median(q))) + 1e-9))
    z, zraw = np.max(zs, axis=0), np.max(zr, axis=0)
    z1 = zr[0].copy()
    if len(exclude):
        ex = np.array([any(a < e and a + P > s_ for s_, e in exclude) for a in t], bool)
        z[ex] = -np.inf; zraw[ex] = -np.inf; z1[ex] = -np.inf
    trim = 3.0
    mask = ndimage.binary_closing(z > k_sigma, np.ones(3))
    lab, n = ndimage.label(mask)
    runs = []
    for sl in ndimage.find_objects(lab):
        a, b = sl[0].start, sl[0].stop
        while a < b and zraw[a] < trim: a += 1
        while b > a and zraw[b - 1] < trim: b -= 1
        if b > a:
            # the smoothed mask misses a strong first/last symbol next to a masked stretch
            while a > 0 and zraw[a - 1] >= trim: a -= 1
            while b < len(zraw) and zraw[b] >= trim: b += 1
        if b - a >= min_len:
            if runs and a <= runs[-1][1]:
                runs[-1] = (runs[-1][0], max(b, runs[-1][1]))
            else:
                runs.append((a, b))
    return runs, t, zs[0], z1


N_SUB = 4                            # sub-bands of the occupied band for _segment


def _sub_slices(ch, car, ofs, occ):
    """N_SUB equal sub-bands of the occupied band (same output rate and timing as the
    full-band slice)."""
    w = occ / N_SUB
    return [ch.slice(car + (i - (N_SUB - 1) / 2) * w, w / 2, ofs) for i in range(N_SUB)]


MAX_BRIDGE_S = 0.6e-3                # longest masked gap bridged inside one burst
FRAME_S = (2.5e-3, 5e-3)             # TDD frame periods tried by _frame_complete


def _frame_complete(runs, t, zraw, P, ofs):
    """Low-SNR burst completion from the TDD frame structure. The per-symbol evidence is
    folded over the frame period (tried: FRAME_S) within the link's active span; a clear
    template gives the burst region of the frame. In every frame where the evidence over that
    region is significant (z > 4), the stretch between the first and last symbol with
    smoothed evidence >= 1 sigma is a burst. D_20260609-151042 ch0 (O4 video at about -5 dB):
    burst 20.13-21.84 ms was cut to 21.20-21.67 ms by the per-symbol tests alone."""
    if len(runs) < 2:
        return runs
    e = np.where(np.isfinite(zraw), np.clip(zraw, -3.0, 6.0), 0.0)
    K = len(e)
    best = None
    for T in FRAME_S:
        F = int(round(T / (P / ofs)))
        if abs(F * P / ofs - T) > 0.02 * T or 3 * F > K:
            continue
        a0, a1 = runs[0][0] // F, min(runs[-1][1] // F, K // F - 1)
        if a1 - a0 + 1 < 3:
            continue
        M = e[a0 * F:(a1 + 1) * F].reshape(-1, F)
        zt = M.mean(0) * np.sqrt(len(M))                      # each e ~ N(0, 1) under noise
        con = float(np.sort(zt)[-max(F // 4, 1):].mean())
        if best is None or con > best[0]:
            best = (con, F, zt, a0, a1)
    if best is None or best[0] < 6.0:
        return runs
    con, F, zt, a0, a1 = best
    # rotate the frame so that it starts at the template minimum (no burst straddles it)
    off = int(np.argmin(ndimage.uniform_filter1d(zt, 5, mode="wrap")))
    ztr = np.roll(zt, -off)
    lab, n = ndimage.label(ztr > 3.0)
    regions = [(sl[0].start, sl[0].stop) for sl in ndimage.find_objects(lab) if sl[0].stop - sl[0].start >= 3]
    add, present = [], []
    for f in range(a0, a1 + 1):
        for ra, rb in regions:
            lo, hi = f * F + off + ra, min(f * F + off + rb, K)
            lo = max(lo, 0)
            if hi - lo < 3:
                continue
            seg = e[lo:hi]
            if seg.mean() * np.sqrt(len(seg)) < 4.0:
                continue
            present.append((lo, hi))
            # the contiguous stretch (dips of up to 3 symbols closed) with the most evidence: a
            # burst shorter than the template is not stretched to it by a stray noise symbol
            on_ = ndimage.uniform_filter1d(seg, 4) >= 1.0
            idx = np.where(on_)[0]
            for i0, i1 in zip(idx[:-1], idx[1:]):
                if 1 < i1 - i0 <= 4:
                    on_[i0:i1] = True
            lb, nb = ndimage.label(on_)
            st = [(float(seg[sl[0]].sum()), sl[0].start, sl[0].stop) for sl in ndimage.find_objects(lb)
                  if sl[0].stop - sl[0].start >= 3]
            if st:
                _, sa, sb = max(st)
                # the 4-symbol smoothing spills ~2 symbols past each edge: trim to symbols that
                # are on by their own evidence
                while sa < sb and seg[sa] < 1.0: sa += 1
                while sb > sa and seg[sb - 1] < 1.0: sb -= 1
                if sb - sa >= 3:
                    add.append((lo + sa, lo + sb))
    out = []
    for a, b in sorted(runs + add):
        # inside a frame region where the burst is present, a dip of up to 6 symbols (0.2 ms at
        # 30 kHz; gaps between DJI bursts are >= 0.8 ms) is the weak signal, not a gap
        same = out and any(lo <= out[-1][1] and a < hi for lo, hi in present) and a - out[-1][1] <= 6
        if out and (a <= out[-1][1] + 1 or same):     # a 1-symbol gap is never a TDD gap
            out[-1] = (out[-1][0], max(b, out[-1][1]))
        else:
            out.append((a, b))
    return out


def _bridge(runs, t, P, subs, max_gap):
    """Join runs separated by a short gap in which the link cannot be observed: every
    sub-band carries at least the power of the neighbouring link symbols (a stronger
    signal covers the whole band). A gap where any sub-band drops below that level is
    a real gap and stays. Both runs are on the same grid (one capture-level phase)."""
    if len(runs) < 2:
        return runs
    cums = [np.concatenate(([0.0], np.cumsum(np.abs(y) ** 2))) for y in subs]
    Pi = int(round(P))

    def pw(c, ks):
        return np.array([(c[min(int(round(t[k])) + Pi, len(c) - 1)] - c[min(int(round(t[k])), len(c) - 1)]) / Pi
                         for k in ks])
    out = [runs[0]]
    for a, b in runs[1:]:
        a0, b0 = out[-1]
        gap = range(b0, a)
        if 0 < len(gap) <= max_gap:
            near = list(range(max(a0, b0 - 8), b0)) + list(range(a, min(b, a + 8)))
            if all(np.median(pw(c, gap)) >= np.median(pw(c, near)) for c in cums):
                out[-1] = (a0, b)
                continue
        out.append((a, b))
    return out


def _avg_spectrum(y, N, cp, t_syms, shift_hz, ofs):
    """Symbol-synchronous power spectrum (fftshifted) averaged over t_syms."""
    n = np.arange(len(y))
    yc = y * np.exp(-2j * np.pi * shift_hz * n / ofs)
    S = np.zeros(N)
    for t in t_syms:
        s = int(round(t)) + cp
        if 0 <= s and s + N <= len(yc):
            S += np.abs(np.fft.fftshift(np.fft.fft(yc[s:s + N]))) ** 2
    return S


def _med_spectrum(y, N, cp, t_syms, shift_hz, ofs):
    """Symbol-synchronous power spectrum (fftshifted), per-bin median over t_syms."""
    n = np.arange(len(y))
    yc = y * np.exp(-2j * np.pi * shift_hz * n / ofs)
    rows = [np.abs(np.fft.fftshift(np.fft.fft(yc[int(round(t)) + cp:int(round(t)) + cp + N]))) ** 2
            for t in t_syms if 0 <= int(round(t)) + cp and int(round(t)) + cp + N <= len(yc)]
    return np.median(rows, axis=0) if rows else np.zeros(N)


def _carrier(ch, f_off, ofs, N, cp, t_syms, L_used):
    """Transmitter carrier (DC subcarrier) relative to fc, three steps:
    1. fractional part (+-SCS/2) from the CP correlation phase;
    2. coarse band position on a slice as wide as the output rate allows (the scan
       grid can be several MHz off, and a narrow slice would cut the band edge);
    3. integer subcarrier from DJI's empty DC bin, the weakest bin near the centre
       of the symbol-synchronous averaged spectrum."""
    scs = ofs / N
    y = ch.slice(f_off, ofs / 2 * 0.98, ofs)
    m, cs = _cp_metric(y, N, cp)
    idx = [int(round(t)) for t in t_syms if int(round(t)) < len(cs)]
    frac = -np.angle(np.sum(cs[idx])) / (2 * np.pi) * scs
    # per-bin MEDIAN over the symbols, not the mean: an interferer in fewer than half of the
    # symbols (Wi-Fi 14 dB above the video, D_20260515-192835) would otherwise pull the
    # band position onto itself (carrier 5 MHz off -> ZC test failed)
    S = _med_spectrum(y, N, cp, idx, frac, ofs)
    w = np.convolve(S, np.ones(L_used), "valid")
    c0 = int(np.argmax(w)) + L_used // 2
    loc = ndimage.uniform_filter1d(S, 7)
    cand = range(max(c0 - 3, 1), min(c0 + 4, N - 1))
    dip = min(cand, key=lambda k: S[k] / (loc[k] + 1e-20))
    k_int = dip - N // 2
    return f_off + k_int * scs + frac, S, k_int


def _symbol_spectra(y, ofs, N, cp, t_syms, shift_hz):
    n = np.arange(len(y))
    yc = y * np.exp(-2j * np.pi * shift_hz * n / ofs)
    out = []
    for t in t_syms:
        s = int(round(t)) + cp
        out.append(np.fft.fftshift(np.fft.fft(yc[s:s + N])) if s + N <= len(yc) else None)
    return out


def _clock_ppm(m, t_syms, P, ncls):
    """Symbol-clock offset: refine every on-grid CP peak (centroid) and fit
    t_k = a_(k mod ncls) + k * P_true  (per-class offsets absorb the longer CP of
    the first symbol in each slot)."""
    ks, ts = [], []
    med = np.median(m)
    for k, t in t_syms:
        a = int(round(t))
        seg = m[max(a - 6, 0):a + 7]
        if len(seg) < 13:
            continue
        w = np.maximum(seg - med, 0)
        if w.sum() <= 0:
            continue
        ks.append(k); ts.append(a - 6 + float(np.sum(w * np.arange(13)) / w.sum()))
    if len(ks) < 20 or max(ks) - min(ks) < 50:
        return None, None
    ks = np.array(ks); ts = np.array(ts)
    A = np.zeros((len(ks), ncls + 1)); A[:, 0] = ks
    A[np.arange(len(ks)), 1 + ks % ncls] = 1
    sol, *_ = np.linalg.lstsq(A, ts, rcond=None)
    res = ts - A @ sol
    sk = np.sqrt(np.sum(res ** 2) / max(len(ks) - ncls - 1, 1) / np.sum((ks - ks.mean()) ** 2))
    return (sol[0] / P - 1) * 1e6, sk / P * 1e6


def _frame_period(starts_ms):
    """TDD frame cadence of a burst train: constant / alternating spacing
    (spectrum_survey._cycle), else the largest multiple of 0.5 ms that >= 80 % of
    the starts are commensurate with (bursts start on the slot grid)."""
    s = np.sort(np.asarray(starts_ms))
    if len(s) < 3:
        return None
    T = O15._cycle(s)
    if T is not None:
        return T
    iv = np.diff(s)
    best = None
    for Tc in np.arange(1.0, 10.01, 0.5):
        r = np.mod(s - s[0] + Tc / 2, Tc) - Tc / 2
        if np.mean(np.abs(r) < 0.1) >= 0.8 and Tc <= np.median(iv) + 0.1:
            best = Tc
    return best


def _video_link(ch, c, fc, hops=()):
    """Segment one capture-level candidate into bursts and measure the link."""
    name, scs, ofs, N, cp, occ, step, L = c["profile"]
    y = ch.slice(c["f_off"], occ / 2 * 1.05, ofs)
    m, _ = _cp_metric(y, N, cp)
    runs, t, _, _ = _segment(m, c["phase"], c["P"])
    if not runs:
        return None
    syms = [k for a, b in runs for k in range(a, b)]
    L_used = int(round(occ / (ofs / N))) + 1
    car, S, _ = _carrier(ch, c["f_off"], ofs, N, cp, t[syms], L_used)
    S_off = c["f_off"]
    # The strongest scan centre can be several MHz off the band centre; its slice then cuts
    # the band edge and biases the band position (D_20260515-192148: 20 MHz Mini 2 video
    # behind 5 GHz Wi-Fi, scan peak at +8.4 MHz, carrier estimated 2.35 MHz high, so the
    # +-2-bin ZC test failed). Re-centre the slice on the estimate until it stops moving.
    for _ in range(3):
        car2, S2, _ = _carrier(ch, car, ofs, N, cp, t[syms], L_used)
        moved = abs(car2 - car) > 1.5 * ofs / N
        car, S, S_off = car2, S2, car
        if not moved:
            break
    # re-slice centred on the transmitter's DC and redo the grid there: the scan
    # centre may be off by up to step/2 and the slice edge cuts into the band
    f_off = car
    y = ch.slice(f_off, occ / 2 * 1.05, ofs)
    m, _ = _cp_metric(y, N, cp)
    # 2026-10-10, D_20260515-192835: Wi-Fi 14 dB above the video over 3/4 of its band cut
    # the 0-3.05 ms burst at 1.73 ms; the clean top sub-band still carries the CP structure.
    # Where a wider signal covers the whole band for a moment, _bridge joins the two halves.
    subs = _sub_slices(ch, car, ofs, occ)
    # short runs are kept until bridging: the first symbols of a burst can stand alone before
    # a masked stretch (4.90 ms on D_20260515-192835)
    ex = [(u.start_us * 1e-6 * ofs, u.end_us * 1e-6 * ofs) for u in hops
          if abs((u.center_mhz * 1e6 - fc) - car) <= occ / 2]
    runs, t, z_full, z1 = _segment_multi([m] + [_cp_metric(u, N, cp)[0] for u in subs], c["phase"], c["P"],
                                     min_len=1, exclude=ex)
    # Sub-band evidence only extends or bridges a burst: a run must contain a symbol where the
    # full-band test alone fires. A narrowband signal of the same numerology in one sub-band
    # (a weak RC uplink hop, tv_2p4_122710 18.79 ms) is not a video burst.
    runs = [r for r in _bridge(runs, t, c["P"], subs, int(MAX_BRIDGE_S / (c["P"] / ofs)))
            if r[1] - r[0] >= 3 and np.any(z_full[r[0]:r[1]] > 4.0)]
    runs = _frame_complete(runs, t, z1, c["P"], ofs)
    runs = _bridge(runs, t, c["P"], subs, int(MAX_BRIDGE_S / (c["P"] / ofs)))
    if not runs:
        return None
    syms = [k for a, b in runs for k in range(a, b)]
    # DJI evidence and the clock fit use only the symbols the full-band test finds on its own:
    # bridged / sub-band-only symbols are masked or weak and dilute the slot statistics (Mini 2
    # + 5 dB noise: z 6.8 < 8 with them)
    ev_runs = _segment(m, c["phase"], c["P"])[0]
    ev_syms = [k for a, b in ev_runs for k in range(a, b)]
    return dict(c=c, f_off=f_off, carrier=car, y=y, m=m, runs=runs, t=t, syms=syms, ev_syms=ev_syms or syms,
                S=S, S_off=S_off, n_on=len(syms))


def _zc_exact(best, acc, L, N):
    """Exact DJI-ZC test on the strongest symbols. best: [(DD coherence, fftshifted
    spectrum)], strongest first; acc: their coherent DD sum (root estimate).
    Returns (zc_root, zc_corr, ref, bin_shift, top_coh, ok)."""
    h = (L - 1) // 2
    # Per symbol (8 strongest of that slot position): exact DJI-ZC correlation, max
    # over cyclic shifts and over roots within +-L/60 of its DD root estimate (slots
    # use different roots -- Mini 2: 192 / 1 / 137 -- and the DD estimate scatters at
    # low SNR). Noise level of one such max ~0.13 (L=601) / 0.10 (L=1201).
    # Exact-sequence test, coherent over many symbols: for a candidate root, the
    # per-symbol correlation X . conj(Z) brought back to time gives a peak at the
    # symbol's timing error, which the grid bounds to tens of samples (first-symbol CP
    # is longer) -> only +-24 lags
    # are searched (small noise bias, unlike a search over all L cyclic shifts),
    # and the per-symbol peaks are averaged over every symbol of that slot position.
    sel = [X for _, X in best[:64]]
    # root candidates: the ZC correlation is sharp in root (+-1 already loses it),
    # while per-symbol DD estimates scatter by ~+-10 at low SNR and slots carry
    # different roots -> modes of a smoothed histogram of the per-symbol estimates
    r0 = _root_from_dd(acc, L)
    hist = np.zeros(L)
    for p_, X in best:
        hist[_root_from_dd(_zc_dd(X, L)[0], L)] += p_
    hs = np.convolve(np.r_[hist[-8:], hist, hist[:8]], np.ones(17), "valid")
    modes = []
    for i in np.argsort(hs)[::-1]:
        if all(min((i - m_) % L, (m_ - i) % L) > 12 for m_ in modes):
            modes.append(int(i))
        if len(modes) == 3:
            break
    cand = set()
    for m_ in modes + [r0]:
        cand |= {(m_ + d) % L for d in range(-10, 11)}
    lags = np.r_[0:25, N - 24:N]
    # the DC-null carrier estimate can be a few subcarriers off at low SNR: also
    # search integer bin shifts (the winning shift corrects the link carrier)
    Xfull = np.array(sel)
    Xfull /= np.linalg.norm(Xfull[:, N // 2 - h:N // 2 + h + 1], axis=1, keepdims=True) + 1e-20
    scores, shifts = {}, {}
    for r in cand - {0}:
        Z = np.zeros(N, complex); Z[N // 2 - h:N // 2 + h + 1] = zc_freq(r, L) / np.sqrt(L - 1)
        for dsh in range(-2, 3):
            q = np.abs(np.fft.ifft(np.fft.ifftshift(np.roll(Xfull, -dsh, axis=1) * np.conj(Z), axes=1),
                                   axis=1)[:, lags]) * N
            v = float(np.mean(q.max(1))) * (1.0 if dsh == 0 else 0.92)   # prefer the DC-null estimate
            if v > scores.get(r, 0):
                scores[r], shifts[r] = v, dsh
    zc_root = max(scores, key=scores.get)
    zc_corr = scores[zc_root]
    # reference: the same statistic for roots far from the best (noise + any
    # generic quadratic-phase similarity)
    others = [v for r, v in scores.items() if min((r - zc_root) % L, (zc_root - r) % L) > 3]
    ref = float(np.median(others)) if others else 0.0
    top_coh = float(np.mean([p_ for p_, _ in best[:64]]))
    # A quadratic-phase sequence of another length/mapping (Herelink's reference
    # symbols) has high DD coherence but matches no DJI ZC well: the match must be
    # commensurate with sqrt(coherence) of the same symbols (DJI ~0.6-1, Herelink 0.2)
    ok = bool(zc_corr - ref > 0.08 and zc_corr > 0.5 * np.sqrt(max(top_coh, 1e-6)))
    return zc_root, zc_corr, ref, shifts[zc_root], top_coh, ok


def _dji_evidence(lk):
    """15 kHz: coherent ZC test of each slot position (k mod 7), accumulated over
    every on-grid symbol of the link. 30 kHz (O4): correlation of each symbol
    with the one 4 symbols later (repeated reference) vs 3 and 5 later."""
    name, scs, ofs, N, cp, occ, step, L = lk["c"]["profile"]
    y, t, syms = lk["y"], lk["t"], lk["ev_syms"]
    ncls = 7 if scs == 15 else 14
    if L:
        spec = _symbol_spectra(y, ofs, N, cp, t[syms], 0.0)
        acc = np.zeros(ncls, complex); mag = np.zeros(ncls); per = []
        for k, X in zip(syms, spec):
            if X is None:
                per.append(0.0); continue
            a, mg = _zc_dd(X, L)
            acc[k % ncls] += a; mag[k % ncls] += mg
            per.append(abs(a) / (mg + 1e-20))
        # Slot position of the ZC symbol: highest MEAN per-symbol coherence. Not the
        # coherent sum over symbols: roots differ between slots (Mini 2: 192 / 1 / 137),
        # so a coherent sum favours any weak constant pattern over the real ZC slot.
        coh = np.abs(acc) / (mag + 1e-20)
        pk = np.array([k % ncls for k in syms]); pv = np.array(per)
        means = np.array([pv[pk == c_].mean() if np.any(pk == c_) else 0 for c_ in range(ncls)])
        cls = int(np.argmax(means))
        others = pv[pk != cls]
        n_cls = int(np.sum(pk == cls))
        sd = others.std() + 1e-9
        z = float((means[cls] - np.median(others)) / (sd / np.sqrt(max(n_cls, 1))))
        # The DD test fires on ANY quadratic-phase sequence (e.g. LTE-style DMRS of
        # another length -- Herelink has one). DJI-specific: a length-L ZC in DJI's
        # mapping. Root from the accumulated DD, then the exact sequence: per-symbol
        # correlation at every cyclic shift, |.|^2 accumulated over the symbols of
        # that slot position (the shift is set by the grid timing, common to all).
        best = sorted([(p_, X) for k, X, p_ in zip(syms, spec, per) if X is not None and k % ncls == cls],
                      key=lambda v: -v[0])[:48]
        zc_root, zc_corr, ref, bshift, top_coh, ok = _zc_exact(best, acc[cls], L, N)
        return dict(kind="zc_slot_symbol", z=z, coh=round(top_coh, 3), slot_pos=cls,
                    zc_corr=round(zc_corr, 3), zc_ref=round(ref, 3), zc_root=zc_root if ok else None, ok=ok,
                    bin_shift=bshift if ok else 0)


    # O4: some symbols repeat exactly 4 symbols later (reference symbols; which
    # positions varies between the captures seen, so no position is assumed).
    # Per symbol: |corr| with the symbol 4 later minus the mean for 3 and 5 later
    # (unrelated symbols, same spacing statistics). Null mean 0.
    on = set(syms)

    def corr(k, lag):
        s = int(round(t[k])) + cp
        s2 = int(round(t[k + lag])) + cp
        if s2 + N + 3 > len(y):
            return None
        a_ = y[s:s + N]
        # grid positions are good to ~1 sample and a 9 MHz signal decorrelates within
        # ~2 samples: best of +-3 samples (same bias for every lag)
        return max(abs(np.vdot(a_, y[s2 + d:s2 + d + N])) /
                   np.sqrt(np.vdot(a_, a_).real * np.vdot(y[s2 + d:s2 + d + N], y[s2 + d:s2 + d + N]).real + 1e-20)
                   for d in range(-3, 4))
    ks, d4, d35 = [], [], []
    for k in syms:
        if not all(k + l in on for l in (3, 4, 5)):
            continue
        v = [corr(k, l) for l in (3, 4, 5)]
        if None in v:
            continue
        ks.append(k); d4.append(v[1] - (v[0] + v[2]) / 2); d35.append(v[0] - v[2])
    if len(d4) < 30:
        return dict(kind="repeated_ref_symbol_4", z=0.0)
    ks = np.array(ks); d4 = np.array(d4); d35 = np.array(d35)
    # var(v) from the 3-vs-5 difference; var(d4) = 1.5 var(v)
    sd = np.sqrt(1.5 * np.var(d35) / 2 / len(d4)) + 1e-9
    z_all = float(np.mean(d4) / sd)
    # The repeat sits at two slot positions 8 apart (a comb-4 pilot pair, like NR's TRS):
    # 4 and 12 of the 14-symbol slot on D_20260609-151042, 5 and 13 on the Mini 4 Pro; the
    # other pairs only add noise. Best of the 14 rotations of {r, r+8} (the grid's slot
    # phase is arbitrary): about twice the all-pairs z on every O4 link seen.
    sd1 = np.sqrt(1.5 * np.var(d35) / 2) + 1e-9
    pos = (0.0, -1, 0)
    for r in range(ncls):
        sel = np.isin(ks % ncls, [r, (r + 8) % ncls])
        if sel.sum() >= 4:
            zr_ = float(np.mean(d4[sel]) / (sd1 / np.sqrt(sel.sum())))
            if zr_ > pos[0]:
                pos = (zr_, r, int(sel.sum()))
    return dict(kind="repeated_ref_symbol_4", z=max(z_all, pos[0]), z_all=round(z_all, 2), slot_pair=pos[1],
                n_pair=pos[2], excess=round(float(np.mean(d4)), 4), n=len(d4))


def _local_bursts(m, P):
    """Symbol starts found one by one (CP-metric peaks spaced ~P), grouped into bursts of
    >= 3 symbols. For a link whose symbol grid is NOT continuous across bursts, where the
    capture-level grid keeps only the bursts that happen to sit on it."""
    thr = 4 * np.median(m)
    Pi = int(round(P))
    starts, i = [], 0
    while i + Pi < len(m):
        k = int(np.argmax(m[i:i + Pi]))
        if m[i + k] > thr:
            starts.append(i + k); i += k + int(0.9 * P)
        else:
            i += Pi
    groups, cur = [], starts[:1]
    for s in starts[1:]:
        # P is the slot average (the first CP of a slot is 8 samples longer, and this link's grid
        # also jumps 15-20 samples inside a burst, well within the 72-sample CP); up to 3 symbols
        # in a row may be missing inside a burst (Wi-Fi over the band on the 2484 MHz capture)
        n = int(round((s - cur[-1]) / P))
        if 1 <= n <= 4 and abs(s - cur[-1] - n * P) <= 24:
            cur.append(s)
        else:
            groups.append(cur); cur = [s]
    groups.append(cur)
    # a burst's first symbol can stand further off: the 2484 MHz link's ZC symbol leads each
    # burst 1225 samples (P + 128) before the next symbol
    merged = []
    for g in groups:
        if merged and len(merged[-1]) == 1 and P - 12 <= g[0] - merged[-1][0] <= P + 256:
            merged[-1] = merged[-1] + g
        else:
            merged.append(g)
    return [g for g in merged if len(g) >= 3]


def _retimed_link(ch, lk):
    """The link re-timed burst by burst: every burst on its own symbol grid, carrier
    re-estimated from those symbols. 2484 MHz capture (test_signal-20351019, 10 MHz
    OcuSync at 2466.5 MHz): the grid moves ~10 samples per burst gap, ~400 samples over
    20 ms, so the capture-level grid kept 138 of ~330 symbols and none of the 3 ZC symbols."""
    name, scs, ofs, N, cp, occ, step, L = lk["c"]["profile"]
    groups = _local_bursts(lk["m"], lk["c"]["P"])
    if not groups:
        return None
    t = np.array([s for g in groups for s in g], float)
    runs, a = [], 0
    for g in groups:
        runs.append((a, a + len(g))); a += len(g)
    car, S, _ = _carrier(ch, lk["carrier"], ofs, N, cp, t, int(round(occ / (ofs / N))) + 1)
    y = ch.slice(car, occ / 2 * 1.05, ofs)
    m, _ = _cp_metric(y, N, cp)
    return dict(lk, f_off=car, carrier=car, y=y, m=m, runs=runs, t=t, syms=list(range(len(t))),
                ev_syms=list(range(len(t))), S=S,
                S_off=lk["carrier"], n_on=len(t), retimed=True)


def _sparse_zc_evidence(lk):
    """DJI-ZC test without a slot position: the symbols whose DD coherence stands > 8 sigma
    above the link's own symbols, then the exact-sequence test on them. For re-timed links,
    where slot positions are not known across bursts and Wi-Fi in the band leaves only a
    few clean ZC symbols (3 in the 2484 MHz capture)."""
    name, scs, ofs, N, cp, occ, step, L = lk["c"]["profile"]
    spec = [X for X in _symbol_spectra(lk["y"], ofs, N, cp, lk["t"], 0.0) if X is not None]
    if len(spec) < 20:
        return dict(kind="zc_any_symbol", z=0.0, ok=False)
    dd = [_zc_dd(X, L) for X in spec]
    per = np.array([abs(a) / (mg + 1e-20) for a, mg in dd])
    med = float(np.median(per))
    sd = 1.4826 * float(np.median(np.abs(per - med))) + 1e-9
    top = [i for i in np.argsort(per)[::-1][:48] if per[i] > med + 8 * sd]
    if len(top) < 2:
        return dict(kind="zc_any_symbol", z=0.0, ok=False, n_zc=len(top))
    z = float((per[top].mean() - med) / (sd / np.sqrt(len(top))))
    best = [(float(per[i]), spec[i]) for i in top]
    zc_root, zc_corr, ref, bshift, top_coh, ok = _zc_exact(best, sum(dd[i][0] for i in top), L, N)
    return dict(kind="zc_any_symbol", z=z, coh=round(top_coh, 3), n_zc=len(top), zc_corr=round(zc_corr, 3),
                zc_ref=round(ref, 3), zc_root=zc_root if ok else None, ok=ok, bin_shift=bshift if ok else 0)


EV_LIKE = 5.0                        # evidence z of a "-like" DJI video link (confirmed: > 8)
MAX_GRIDS = 2                        # further symbol grids searched per confirmed link (_regrid)
MIN_GRID_SYMS = 20                   # symbols a further grid must add to count


def _regrid(ch, lk, fc, hops, z_thr):
    """Bursts of a confirmed link on another symbol grid. A link can re-time mid-capture
    (D_20260609-151042: the O4 grid jumps by 316 samples, 10.3 us, at ~16.8 ms, and its TDD frame
    phase with it); the capture-level fold then keeps one grid and the bursts before the jump are
    lost. Further peaks of the link's own CP fold (same carrier and profile, more than a CP away
    from the grids already used) are segmented on their own phase; their runs are kept where they
    do not overlap the link in time. Returns segment links (they share the parent's identity)."""
    P, cp = lk["c"]["P"], lk["c"]["profile"][4]
    f = _robust_z(_fold(lk["m"], P))
    nb = len(f)
    taken = [lk["c"]["phase"]]
    busy = [(lk["t"][a], lk["t"][b - 1] + P) for a, b in lk["runs"]]
    out = []
    for _ in range(MAX_GRIDS):
        far = [k for k in range(nb) if all(min(abs(k - p), nb - abs(k - p)) > cp + 16 for p in taken)]
        if not far:
            break
        k = max(far, key=lambda i: (f[i], -i))
        if f[k] < z_thr:
            break
        taken.append(k)
        sg = _video_link(ch, dict(lk["c"], phase=k, z=float(f[k]), f_off=lk["carrier"]), fc, hops)
        if sg is None:
            continue
        t = sg["t"]
        runs = [(a, b) for a, b in sg["runs"] if not any(t[a] < e and t[b - 1] + P > s for s, e in busy)]
        if sum(b - a for a, b in runs) < MIN_GRID_SYMS:
            continue
        sg["grid_ev"] = _dji_evidence(sg)
        sg.update(runs=runs, parent=lk, ev=lk["ev"], carrier=lk["carrier"], S=lk["S"], S_off=lk["S_off"],
                  image_of=None)
        busy += [(t[a], t[b - 1] + P) for a, b in runs]
        out.append(sg)
    return out


def detect_video(ch, x, fs, fc, snr, usable_half, z_thr=8.0, ev_thr=8.0, log=None, hops=(), tf=None):
    """OcuSync video: capture-level CP fold -> segmentation on the grid ->
    per-link DJI test -> per-burst features."""
    cands = sorted(_scan_video(ch, usable_half, z_thr), key=lambda c: -c["z"])
    links = []

    def blocks(l, c, f):
        # A link blocks a candidate of the same SCS near it only if it passed the DJI test or
        # has the same profile: the 15k_20MHz hypothesis of a 10 MHz link folds as strongly as
        # the right one, fails the ZC test (wrong length), and must not hide the 10 MHz profile
        # (D_20260916-171218: 10 MHz Mini 2 video lost behind its own 20 MHz candidate).
        # A link that failed the DJI test hides only another hypothesis of itself: the same
        # profile at the same refined carrier. Hiding every same-profile candidate within half
        # the band lost D_20260515-192148, whose right candidate (5769 MHz) sat 4 MHz from a
        # failed one whose carrier estimate was wrong.
        same = l["c"]["profile"][1] == c["profile"][1]
        if l["ev"].get("ok", True) and l["ev"]["z"] > 8.0:
            return same and abs(l[f] - c["f_off"]) < max(l["c"]["profile"][5], c["profile"][5]) / 2
        return l["c"]["profile"][0] == c["profile"][0] and abs(l["carrier"] - c["f_off"]) < 0.5e6

    for c in cands:
        if len(links) >= 8:
            break
        # cheap pre-dedup on the scan centre before the expensive refinement
        if any(blocks(l, c, "f_off") for l in links):
            continue
        lk = _video_link(ch, c, fc, hops)
        if lk is None:
            continue
        if any(blocks(l, dict(c, f_off=lk["carrier"]), "carrier") for l in links):
            continue
        lk["ev"] = _dji_evidence(lk)
        if c["profile"][7] and not lk["ev"].get("ok", True):
            # no ZC on the capture-level grid: the grid may move between bursts
            rt = _retimed_link(ch, lk)
            if rt is not None:
                ev = _sparse_zc_evidence(rt)
                if log:
                    log(f"  retimed {c['profile'][0]} carrier {(fc + rt['carrier']) / 1e6:.4f} MHz "
                        f"bursts {len(rt['runs'])} symbols {rt['n_on']} evidence {ev}")
                if ev["ok"] and ev["z"] > ev_thr:
                    lk = rt; lk["ev"] = ev
        links.append(lk)
    links = [l for l in links if not any(o is not l and o["ev"].get("ok", True) and o["ev"]["z"] > 8.0
                                         and not (l["ev"].get("ok", True) and l["ev"]["z"] > 8.0)
                                         and o["c"]["profile"][1] == l["c"]["profile"][1]
                                         and abs(o["carrier"] - l["carrier"]) < 9e6 for o in links)]
    # IQ-imbalance image: a weak copy mirrored about fc (ZC root L-u), same timing
    for lk in links:
        lk["image_of"] = None
        for o in links:
            if o is lk or o["c"]["profile"][1] != lk["c"]["profile"][1]:
                continue
            if abs(o["carrier"] + lk["carrier"]) < 0.3e6 and o["c"]["z"] > lk["c"]["z"] * 1.5:
                lk["image_of"] = (fc + o["carrier"]) / 1e6
    # Tiers: "dji" -- evidence z > ev_thr, or an overwhelming exact-sequence match (noise
    # ~0.1-0.13, Herelink <= 0.2; a saturated Mini 2 at 1 m, tv_5p8_125115, matches 0.51 vs 0.05
    # at z 6.2); "like" -- z in [EV_LIKE, ev_thr) with the DJI structure present (15 kHz: the
    # exact ZC test passed), reported as an inferred "-like" label at confidence <= 0.6 instead
    # of being dropped; one per band, never next to a confirmed link.
    tier = {}
    for lk in links:
        ev = lk["ev"]
        if ev.get("ok", True) and lk["image_of"] is None and (
                ev["z"] > ev_thr or ev.get("zc_corr", 0.0) - ev.get("zc_ref", 0.0) > 0.3):
            tier[id(lk)] = "dji"
    for lk in sorted(links, key=lambda l: -l["ev"]["z"]):
        ev = lk["ev"]
        if (id(lk) not in tier and ev.get("ok", True) and lk["image_of"] is None and ev["z"] >= EV_LIKE
                and not any(id(o) in tier and o["c"]["profile"][1] == lk["c"]["profile"][1]
                            and abs(o["carrier"] - lk["carrier"]) < 9e6 for o in links)):
            tier[id(lk)] = "like"
    # further symbol grids of every reported link (a link re-timed mid-capture)
    lid = {id(lk): li for li, lk in enumerate(links)}
    for lk in list(links):
        if id(lk) in tier and not lk.get("retimed"):
            for sg in _regrid(ch, lk, fc, hops, z_thr):
                tier[id(sg)] = tier[id(lk)]; lid[id(sg)] = lid[id(lk)]
                links.append(sg)
    groups = {}
    for lk in links:
        li = lid[id(lk)]
        name, scs, ofs, N, cp, occ, step, L = lk["c"]["profile"]
        ev = lk["ev"]
        is_dji = id(lk) in tier
        if log:
            log(f"video link {name} carrier {(fc + lk['carrier']) / 1e6:.4f} MHz z_cp {lk['c']['z']:.1f} "
                f"bursts {len(lk['runs'])} evidence {ev} image_of {lk['image_of']} -> "
                f"{({'dji': 'DJI', 'like': 'DJI-like'}[tier[id(lk)]]) if is_dji else 'rejected'}"
                + (f" (grid {lk['c']['phase']} of the link above, its own evidence {lk['grid_ev']['z']:.1f})"
                   if "parent" in lk else ""))
        if not is_dji:
            continue
        if "parent" in lk:
            lk["carrier"] = lk["parent"]["carrier"]          # a further grid: its link's (corrected) carrier
        elif ev.get("bin_shift"):
            # ZC match says the DC bin is elsewhere: correct the carrier (integer SCS)
            lk["carrier"] += ev["bin_shift"] * ofs / N
        t, runs = lk["t"], lk["runs"]
        if hops:
            # The RC uplink is TDD with the video and has the same numerology: a video symbol
            # overlapping an in-band hop is the hop (D_20260515-192835: burst started inside
            # the 23.17 ms hop; tv_2p4_122710 + 15 dB: a 3-symbol "burst" inside a hop).
            lo_, hi_ = (fc + lk["carrier"] - occ / 2) / 1e6, (fc + lk["carrier"] + occ / 2) / 1e6
            hs = [(u.start_us, u.end_us) for u in hops if lo_ <= u.center_mhz <= hi_]
            Pus = lk["c"]["P"] / ofs * 1e6
            nr = []
            for a, b in runs:
                cur = None
                for k in range(a, b):
                    tk = t[k] / ofs * 1e6
                    if not any(tk < e_ and tk + Pus > s_ for s_, e_ in hs):
                        cur = (cur[0], k + 1) if cur else (k, k + 1)
                        continue
                    if cur and cur[1] - cur[0] >= 3:
                        nr.append(cur)
                    cur = None
                if cur and cur[1] - cur[0] >= 3:
                    nr.append(cur)
            runs = nr
        lk["runs_out"] = runs
        groups.setdefault(li, []).append(lk)
    out = []
    for li, grp in groups.items():
        lk = grp[0]                                          # the link (further grids follow it)
        name, scs, ofs, N, cp, occ, step, L = lk["c"]["profile"]
        ev = lk["ev"]
        P = lk["c"]["P"]
        S, m = lk["S"], lk["m"]
        ncls = 7 if scs == 15 else 14
        # symbol-clock fit needs one continuous grid; a re-timed link has one per burst
        cppm, cppm_sd = (None, None) if lk.get("retimed") else \
            _clock_ppm(m, [(k, lk["t"][k]) for k in lk["ev_syms"]], P, ncls)
        # occupied width from the symbol-synchronous averaged spectrum (slice wider than
        # the band): edges where the power falls to midway between the in-band level
        # and the out-of-band floor -- a fixed -10 dB rule fails at low SNR, where
        # the whole slice is within 10 dB of the in-band level
        ps = ndimage.uniform_filter1d(S, 9)
        k0 = int(round(N // 2 + (lk["carrier"] - lk["S_off"]) / (ofs / N)))
        k0 = int(np.clip(k0, 0, N - 1))
        inb = np.abs(np.arange(N) - k0) < 0.25 * occ / (ofs / N)
        occm = ps > (np.median(ps[inb]) + np.percentile(ps, 5)) / 2
        # dips inside the band (O4's pilot/null pattern, a faded stretch) must not end the walk:
        # D_20260609-151042, an 18 MHz O4 link, measured 7-11 MHz and was labelled 10 MHz
        occm = ndimage.binary_closing(occm, np.ones(max(3, int(0.6e6 / (ofs / N)))))
        lo = hi = k0
        while lo > 0 and occm[lo - 1]: lo -= 1
        while hi < N - 1 and occm[hi + 1]: hi += 1
        f_lo = float(lk["S_off"] + (lo - N // 2) * ofs / N)
        f_hi = float(lk["S_off"] + (hi + 1 - N // 2) * ofs / N)
        # bursts on every grid of the link, in seconds
        seeds = sorted((float(g["t"][a] / ofs), float((g["t"][b - 1] + P) / ofs)) for g in grp for a, b in g["runs_out"])
        zc_bw = bool(L and ev.get("zc_root") is not None)
        band_note = "symbol-synchronous spectrum edges"
        if tf is not None and seeds:
            # 2026-10-13, D_20260609-151042 ch0 (O4 at ~0 dB): the band from in-band power while the link
            # is on vs off (multipath dips no longer end it: 14.1 -> 18.0 MHz) and the bursts from the
            # in-band power seeded by the CP evidence (no more split / late / short bursts)
            span = f_hi - f_lo if f_hi > f_lo else occ
            if not zc_bw:
                fb = BU.fit_band(tf, seeds, f_lo, f_hi, max(f_lo - 0.5 * span, -usable_half),
                                 min(f_hi + 0.5 * span, usable_half))
                if fb is not None:
                    f_lo, f_hi = fb["lo"], fb["hi"]; band_note = "in-band power on vs off (rectangle fit)"
            elo, ehi = (f_lo, f_hi) if not zc_bw else (lk["carrier"] - (L - 1) * ofs / N / 2,
                                                       lk["carrier"] + (L - 1) * ofs / N / 2)
            seeds, einfo = BU.energy_bursts(tf, seeds, elo, ehi, usable_half, min_sym_s=P / ofs)
            if log:
                log(f"  video link {li}: band {(fc + elo) / 1e6:.3f}-{(fc + ehi) / 1e6:.3f} MHz, energy bursts {einfo}")
            if einfo.get("used") and not zc_bw:
                fb = BU.fit_band(tf, seeds, f_lo, f_hi, max(f_lo - 0.5 * (f_hi - f_lo), -usable_half),
                                 min(f_hi + 0.5 * (f_hi - f_lo), usable_half))
                if fb is not None:
                    f_lo, f_hi = fb["lo"], fb["hi"]; band_note = "in-band power on vs off (rectangle fit)"
        bw_meas = f_hi - f_lo
        bw = bw_meas
        cen = (fc + lk["carrier"]) / 1e6
        if zc_bw:
            # a verified DJI ZC of length L occupies exactly L-1 subcarriers: that is the mode.
            # The spectrum edge walk is a measurement only -- overlapping Wi-Fi (5 GHz ch 157 on
            # D_20260916-171137) or a weak far link cuts it short (4.95 MHz for a 20 MHz link)
            bw = (L - 1) * ofs / N
            b_lo, b_hi = cen - bw / 2e6, cen + bw / 2e6
        else:
            b_lo, b_hi = (fc + f_lo) / 1e6, (fc + f_hi) / 1e6
        bw_label = "10MHz" if bw < 13e6 else ("20MHz" if bw < 26e6 else "40MHz")
        starts_ms = [a * 1e3 for a, b in seeds]
        T = _frame_period(starts_ms)
        for t0s, t1s in seeds:
            # per-burst residual carrier from the CP phase, on the grid of the link that covers the burst
            g = max(grp, key=lambda g: sum(max(0.0, min(t1s, (g["t"][b - 1] + P) / ofs) - max(t0s, g["t"][a] / ofs))
                                           for a, b in g["runs_out"]))
            _, csum = _cp_metric(g["y"], N, cp)
            idx = [int(round(v)) for v in g["t"] if t0s <= v / ofs < t1s and int(round(v)) < len(csum)]
            t0, t1 = t0s * 1e6, t1s * 1e6
            f = dict(scs_khz=scs, symbol_us=round(P / ofs * 1e6, 3), cp_us=round(cp / ofs * 1e6, 2),
                     tu_us=round(N / ofs * 1e6, 2), bw_mhz=round(bw / 1e6, 2), bw_measured_mhz=round(bw_meas / 1e6, 2),
                     burst_us=round(t1 - t0, 1), n_symbols=int(round((t1s - t0s) / (P / ofs))), link_id=f"video{li}",
                     cp_fold_z=round(lk["c"]["z"], 1), evidence=ev["kind"], evidence_z=round(ev["z"], 1),
                     evidence_tier=tier[id(lk)], band_from=band_note)
            f.update(_cfo_feats(fc + lk["carrier"]))
            if idx:
                f["cfo_burst_hz"] = round(f["cfo_hz"] - float(np.angle(np.sum(csum[idx]))) / (2 * np.pi) * (ofs / N), 1)
            if T:
                f["frame_period_ms"] = round(T, 3)
                f["frame_phase_us"] = round(float(np.mod(t0, T * 1e3)), 1)
            if cppm is not None:
                f["clock_ppm"] = round(float(cppm), 2); f["clock_ppm_sd"] = round(float(cppm_sd), 2)
            if ev.get("zc_root") is not None:
                f["zc_root"] = ev["zc_root"]
            lab = ("video_" if scs == 15 else "video_o4_") + bw_label
            out.append(Emission("dji_ocusync", lab, round(t0, 1), round(t1, 1), round(b_lo, 3), round(b_hi, 3),
                                round(cen if zc_bw else (b_lo + b_hi) / 2, 4),
                                snr(t0, t1, (b_lo + b_hi) / 2 - (b_hi - b_lo) / 2.4, (b_lo + b_hi) / 2 + (b_hi - b_lo) / 2.4),
                                confidence=float(min(0.99 if tier[id(lk)] == "dji" else 0.6, 0.5 + ev["z"] / 60)),
                                features=f))
    return out


# ============================================================== RC uplink (hopping, 1.4 MHz)
UP_L = 73                             # ZC length: 72 used subcarriers + DC


def detect_uplink(ch, x, fs, fc, snr, usable_half, k_sigma=5.0, log=None):
    """Blind 1.4 MHz CP scan + hop-train tracking (spectrum_survey), then DJI
    test: symbols 0 and 6 of the 7-symbol slot are length-73 ZC."""
    prof = [O15.PROFILES[0]]
    hits, stats = O15.ofdm15k_scan(ch, fc, profiles=prof, k_sigma=k_sigma, usable_half=usable_half)
    bs = O15.merge_ofdm15k([b for b in (O15.refine_ofdm15k(ch, h) for h in hits) if b is not None])
    bs = O15.merge_ofdm15k(bs + O15.track_hops(ch, fc, bs, stats, len(x) / fs, usable_half=usable_half))
    name, ofs, N, cp, occ, step = O15.PROFILES[0]
    out = []
    for b in bs:
        dur = (b["end_s"] - b["start_s"]) * 1e6
        if b["bw_hz"] > 2.5e6 or not (250 <= dur <= 1200):
            continue
        y = ch.slice(b["f_off"], 0.9e6, ofs)
        s0, e0 = int(b["start_s"] * ofs), int(b["end_s"] * ofs)
        a, e = max(0, s0 - N), min(len(y), e0 + N)
        m, cs = _cp_metric(y[a:e], N, cp)
        pk, _ = find_peaks(m, height=0.3, distance=int(0.8 * (N + cp)))
        zcs = []
        for p in pk:
            s = a + p + cp
            if s + N > len(y):
                continue
            X = np.fft.fftshift(np.fft.fft(y[s:s + N]))
            cc, r = zc_root_search(X[N // 2 - 36:N // 2 + 37], UP_L)
            zcs.append((cc, r, (a + p) / ofs * 1e6, X))
        zcs.sort(key=lambda v: -v[0])
        # DJI: two ZC symbols with the SAME root, 6 symbols (428.6 us) apart (first and
        # last of the slot), above the best-of-72-roots noise level (~0.33 for L=73)
        ok = (len(zcs) >= 2 and zcs[0][0] > 0.45 and zcs[1][0] > 0.33 and zcs[0][1] == zcs[1][1]
              and abs(abs(zcs[0][2] - zcs[1][2]) - 428.6) < 30)
        if log:
            log(f"uplink hop {(fc + b['f_off']) / 1e6:.3f} MHz t {b['start_s'] * 1e3:.3f} ms "
                f"zc {[(round(z[0], 2), z[1], round(z[2], 1)) for z in zcs[:2]]} -> {'DJI' if ok else '?'}")
        if not ok:
            continue
        # carrier: the ZC symbols are constant-modulus except DJI's empty DC bin, so the
        # weakest bin near the centre of |X0|^2+|X6|^2 is the transmitter DC (integer
        # subcarrier); the fraction from the CP phase over the slot
        P2 = np.abs(zcs[0][3]) ** 2 + np.abs(zcs[1][3]) ** 2
        k_dc = min(range(N // 2 - 4, N // 2 + 5), key=lambda k: P2[k]) - N // 2
        strong = m > np.percentile(m, 85)
        frac = -np.angle(np.sum(cs[strong])) / (2 * np.pi) * (ofs / N)
        car = fc + b["f_off"] + k_dc * ofs / N + frac
        cen = car / 1e6
        t0, t1 = b["start_s"] * 1e6, b["end_s"] * 1e6
        f = dict(scs_khz=15, symbol_us=b["period_us"], cp_us=round(cp / ofs * 1e6, 2), tu_us=66.67,
                 bw_mhz=round(b["bw_hz"] / 1e6, 2), burst_us=round(dur, 1), n_symbols=b["nsym"],
                 zc_root=zcs[0][1], zc_corr=round(zcs[0][0], 2), tracked=bool(b.get("tracked", False)))
        f.update(_cfo_feats(car))
        out.append(Emission("dji_ocusync", "uplink", round(t0, 1), round(t1, 1),
                            round(cen - b["bw_hz"] / 2e6, 3), round(cen + b["bw_hz"] / 2e6, 3), round(cen, 4),
                            snr(t0, t1, cen - 0.5, cen + 0.5), confidence=0.8, features=f))
    # IQ-imbalance images: a much weaker copy mirrored about fc at the same time
    # (its ZC root comes out as 73-u)
    fcm = fc / 1e6
    out = [e for e in out if not any(
        o is not e and abs(o.center_mhz + e.center_mhz - 2 * fcm) < 0.15 and abs(o.start_us - e.start_us) < 50
        and (o.snr_db or -99) > (e.snr_db or -99) + 10 for o in out)]
    # hop-train timing: frame period / phase over all hops of the link
    if len(out) >= 3:
        T = _frame_period([e.start_us / 1e3 for e in out])
        if T:
            for e in out:
                e.features["frame_period_ms"] = round(T, 3)
                e.features["frame_phase_us"] = round(float(np.mod(e.start_us, T * 1e3)), 1)
    return out


# ============================================================== DroneID
DID_CP = np.array([80, 72, 72, 72, 72, 72, 72, 72, 80])     # at 15.36 MSps
DID_N = 1024
DID_LEN = 9 * DID_N + int(DID_CP.sum())


def _droneid_cp_score(y):
    """DroneID's complete 9-symbol CP schedule (packetizer._droneid_cp_candidates)
    evaluated on a 15.36 MSps slice: in-band, so ~8 dB more SNR than at the full
    61.44 MHz rate."""
    N = DID_N
    p = y[N:] * np.conj(y[:-N])
    ea = np.abs(y[:-N]) ** 2; eb = np.abs(y[N:]) ** 2
    nmax = len(y) - DID_LEN
    if nmax <= 0:
        return np.zeros(0)
    score = np.zeros(nmax + 1)
    cache = {}
    off = 0
    for i, cpl in enumerate(DID_CP):
        if i:
            if cpl not in cache:
                cache[cpl] = np.abs(_movsum(p, cpl)) / np.sqrt(np.maximum(_movsum(ea, cpl) * _movsum(eb, cpl), 1e-30))
            score += cache[cpl][off:off + nmax + 1]
        off += N + cpl
    return score / 8


DID8_LEN = DID_LEN - DID_N - int(DID_CP[1])      # 8-symbol frame: the first symbol dropped (CP 80 kept first [U])
DID_SYM0_MIN, DID_SYM0_FRAC = 0.3, 0.5


def _droneid_sym_cp(y, start, cp):
    """Normalised CP correlation of the OFDM symbol whose CP starts at `start` (0 where there is no symbol)."""
    N = DID_N
    if start < 0 or start + cp + N > len(y):
        return 0.0
    a, b = y[start:start + cp], y[start + N:start + N + cp]
    return float(np.abs(np.sum(b * np.conj(a))) / np.sqrt(np.sum(np.abs(a) ** 2) * np.sum(np.abs(b) ** 2) + 1e-30))


def _droneid_zc(y, start):
    """Correlation of DroneID symbols 4 and 6 (index 3, 5) with roots 600 / 147,
    plus the integer-bin offset of the occupied band."""
    offs = np.concatenate(([0], np.cumsum(DID_N + DID_CP[:-1])))
    spectra = []
    for i in (3, 5):
        s = start + int(offs[i]) + int(DID_CP[i])
        if s + DID_N > len(y) or s < 0:
            return 0.0, 0.0, 0
        spectra.append(np.fft.fftshift(np.fft.fft(y[s:s + DID_N])))
    L = 601
    w = np.convolve(np.abs(spectra[0]) ** 2 + np.abs(spectra[1]) ** 2, np.ones(L), "valid")
    b0 = int(np.argmax(w))
    sc = []
    for X, r in zip(spectra, (600, 147)):
        cc, _ = zc_root_search(X[b0:b0 + L], L, [r])
        sc.append(cc)
    # integer carrier bin: every DroneID symbol leaves the DC bin empty -> the weakest bin within +-4 of the
    # energy-window centre, with the power summed over all the frame's symbols (the two ZC symbols alone were 4 bins
    # off at 3 dB on synthetic bursts, PLAN D4)
    P2 = np.zeros(DID_N)
    for i in range(len(DID_CP)):
        s = start + int(offs[i]) + int(DID_CP[i])
        if 0 <= s and s + DID_N <= len(y):
            P2 += np.abs(np.fft.fftshift(np.fft.fft(y[s:s + DID_N]))) ** 2
    c0 = b0 + L // 2
    dc = min(range(max(c0 - 4, 0), min(c0 + 5, DID_N)), key=lambda k: P2[k])
    return sc[0], sc[1], dc - DID_N // 2


class _CodeFallback(dict):
    """The DroneID decoder's product-type table with a numeric fallback: an unmapped code comes back as its number."""

    def get(self, k, default=None):
        return super().get(k, str(k))


def droneid_decode(x, fs, fc, center_hz, t0_s, t1_s, pad_s=50e-6, src=DRONEID_SRC):
    """Shift the burst to DC and hand a crop with ~50 us of real capture either
    side to droneid_api.decode_ddc_packet (its readme: needed for decoding)."""
    if src not in sys.path:
        sys.path.insert(0, src)
    try:
        from droneid_api import decode_ddc_packet
        import droneid_packet as _dp
    except Exception as e:                                           # pragma: no cover
        return dict(status=f"import_failed: {e}", success=False)
    if not isinstance(_dp.DRONE_TYPES, _CodeFallback):
        # the decoder's own name table stops at code 70 and reports newer types (Mini 3 Pro 73, Inspire 3 76, ...) as
        # "Unknown": return the code itself instead, so the plugin maps it with core.drone_models.dji_product
        _dp.DRONE_TYPES = _CodeFallback(_dp.DRONE_TYPES)
    s = max(0, int((t0_s - pad_s) * fs)); e = min(len(x), int((t1_s + pad_s) * fs))
    n = np.arange(s, e)
    iq = (x[s:e] * np.exp(-2j * np.pi * (center_hz - fc) / fs * n)).astype(np.complex64)
    r = decode_ddc_packet(iq, sample_rate=fs, strip_pre_buffer=False)
    pk = r.get("packet") or {}
    return dict(status=r.get("status"), success=r.get("success"),
                payload=pk.get("payload"), error=(r.get("error") or {}).get("code"),
                crc24_ok=pk.get("crc24_ok"), crc16_ok=pk.get("crc16_ok"),
                residual_cfo_hz=(r.get("signal") or {}).get("residual_cfo_hz"),
                zc_correlation=(r.get("signal") or {}).get("zc_correlation"))


def detect_droneid(ch, x, fs, fc, snr, usable_half, cp_thr=0.2, step=3e6, decode=True, log=None):
    ofs = 15.36e6
    offs = np.arange(-usable_half + 4.5e6, usable_half - 4.5e6 + 1, step)
    cands = []
    for off in offs:
        y = ch.slice(off, 5.0e6, ofs)
        sc = _droneid_cp_score(y)
        if not sc.size:
            continue
        pk, _ = find_peaks(sc, height=cp_thr, distance=DID_LEN // 2)
        for p in pk:
            # the CP schedule alone also fits OcuSync video (same numerology); the two
            # DroneID ZC symbols are what make it DroneID. The CP score is ambiguous by one symbol (a frame scores
            # almost as well one symbol early), so the ZC pair is tested at p and one symbol either side
            best = None
            for q in (int(p), int(p) + DID_N + int(DID_CP[0]), int(p) - DID_N - int(DID_CP[1])):
                z600, z147, kb = _droneid_zc(y, q)
                if best is None or z600 + z147 > best[0]:
                    best = (z600 + z147, q, z600, z147, kb)
            _, q, z600, z147, kb = best
            # fixed-root correlation noise over 601 bins ~ 3/sqrt(601) = 0.12
            if z600 > 0.25 and z147 > 0.2:
                cands.append((z600 + z147, float(sc[p]), float(off), q, z600, z147, kb))
    cands.sort(reverse=True)
    out, seen = [], []
    for _, score, off, p, z600, z147, kb in cands:
        t0 = p / ofs
        if any(abs(t0 - s) < 0.5e-3 for s in seen):
            continue
        y = ch.slice(off, 5.0e6, ofs)
        if log:
            log(f"droneid cand t {t0 * 1e3:.3f} ms off {off / 1e6:.1f} MHz cp {score:.2f} zc {z600:.2f},{z147:.2f}")
        seen.append(t0)
        # carrier: integer bin from the ZC band, fraction from the CP phase
        pp = y[DID_N:] * np.conj(y[:-DID_N])
        offs_s = np.concatenate(([0], np.cumsum(DID_N + DID_CP[:-1])))
        cs = sum(np.sum(pp[p + int(o):p + int(o) + int(c_)]) for o, c_ in zip(offs_s[1:], DID_CP[1:]))
        frac = np.angle(cs) / (2 * np.pi) * (ofs / DID_N)
        car = fc + off + kb * ofs / DID_N + frac
        # the search slice can cut the 9 MHz band (slices every 3 MHz): the ZC energy window then shifts and the DC bin
        # search (+-4 bins) misses (35 bins on a burst 1 MHz from a slice centre). Re-slice on the coarse carrier and
        # take the integer bin again (synthetic DroneID, PLAN D4)
        for _ in range(3):                               # until the band sits in the slice centre (|bin| <= 1)
            off2 = float(np.round((car - fc) / 1e3) * 1e3)
            kb2 = _droneid_zc(ch.slice(off2, 5.0e6, ofs), p)[2]
            car = fc + off2 + kb2 * ofs / DID_N + frac
            if abs(kb2) <= 1:
                break
        # 9 symbols (OcuSync 2/3) or 8 (legacy: Mavic 2 and older; the first symbol dropped, ZC in symbols 3 / 5):
        # is there an OFDM symbol in the 9-symbol layout's first slot? (research, PLAN D4)
        c0 = _droneid_sym_cp(y, p, int(DID_CP[0]))
        nsym = 9 if c0 >= max(DID_SYM0_MIN, DID_SYM0_FRAC * score) else 8
        if nsym == 8:
            t0 = (p + DID_N + int(DID_CP[0])) / ofs
        t1 = t0 + (DID_LEN if nsym == 9 else DID8_LEN) / ofs
        dec = droneid_decode(x, fs, fc, car, t0, t1) if decode and nsym == 9 else None
        cen = car / 1e6
        f = dict(scs_khz=15, symbol_us=round((DID_N + 72) / ofs * 1e6, 2), cp_us=round(72 / ofs * 1e6, 2),
                 tu_us=66.67, bw_mhz=9.0, burst_us=round((t1 - t0) * 1e6, 1), n_symbols=nsym,
                 zc_root=[600, 147], zc_corr=[round(z600, 2), round(z147, 2)], cp_score=round(score, 2),
                 droneid_frame=("9-symbol (OcuSync 2/3 generation: Mini 2, Air 2 and later)" if nsym == 9 else
                                "8-symbol (first symbol absent: Mavic Pro / Mavic 2 family, and measured on the DJI Mini 4 Pro (O4); payload not decoded)"),
                 symbol0_cp=round(c0, 2))
        f.update(_cfo_feats(car))
        payload = None
        if dec:
            f["decode_status"] = dec.get("status")
            # research (NDSS 2023 / AntSDR reports): newer firmware may encrypt the flight information; a valid
            # outer block (CRC-24) with a failing inner CRC-16 is "DroneID present, payload opaque"
            if dec.get("crc24_ok") and dec.get("crc16_ok"):
                f["payload_state"] = "decoded"
            elif dec.get("crc24_ok"):
                f["payload_state"] = ("opaque: CRC-24 valid, flight-information CRC-16 failed "
                                      "(encrypted or unknown payload layout)")
            if dec.get("success") and dec.get("payload"):
                # CRC-verified flight information: the adapter reports such a burst as "success"
                payload = dict(dec["payload"], crc_ok=bool(dec.get("crc16_ok") or dec.get("crc24_ok")))
                f["serial"] = payload.get("serial_number")
                f["device_type"] = payload.get("device_type")
        out.append(Emission("dji_droneid", "droneid", round(t0 * 1e6, 1), round(t1 * 1e6, 1),
                            round(cen - 4.5, 3), round(cen + 4.5, 3), round(cen, 4),
                            snr(t0 * 1e6, t1 * 1e6, cen - 4, cen + 4),
                            confidence=0.99 if payload else 0.85, decoded=payload, features=f))
    return out


# ============================================================== entry point
def detect(x, fs, fc, ctx=None, log=None):
    """All DJI emissions in a complex-baseband capture.
    ctx (optional): {"channelizer": drone_common.Channelizer over x, "decode_droneid": bool,
                     "z_thr": capture-level CP-fold threshold, "verbose": bool}."""
    ctx = ctx or {}
    if log is None and ctx.get("verbose"):
        log = print
    x = np.asarray(x)
    if ctx.get("channelizer") is None:
        # notch CW tones / LO leakage first: a tone correlates at every lag and masks
        # hops next to it (the far Mini 2 capture has a hop 0.6 MHz from the LO spur)
        notch = [(s["offset_hz"], 30e3) for s in SP.find_spurs(x, fs, fc)]
        ctx["channelizer"] = C.Channelizer(x.astype(np.complex128), fs, notch)
    ch = _FastChannelizer(ctx["channelizer"])
    usable_half = min(fs, 56e6) / 2
    snr = _SNR(x, fs, fc)
    em = []
    up = detect_uplink(ch, x, fs, fc, snr, usable_half, log=log)
    em += detect_video(ch, x, fs, fc, snr, usable_half, z_thr=ctx.get("z_thr", 8.0), log=log, hops=up,
                       tf=ctx.get("tfpower") or BU.TFPower(x, fs))
    em += up
    em += detect_droneid(ch, x, fs, fc, snr, usable_half, decode=ctx.get("decode_droneid", True), log=log)
    em.sort(key=lambda e: (e.start_us, e.center_mhz))
    return em


if __name__ == "__main__":
    import argparse, time
    ap = argparse.ArgumentParser(description="DJI emission detector")
    ap.add_argument("meta"); ap.add_argument("--out"); ap.add_argument("-v", action="store_true")
    a = ap.parse_args()
    x, fs, fc, meta = C.load_capture(a.meta)
    t = time.time()
    em = detect(x, fs, fc, {"verbose": a.v})
    print(f"{len(em)} DJI emissions in {time.time() - t:.1f} s: {dict(Counter(e.label for e in em))}")
    for e in em:
        print(f"  {e.label:16s} {e.start_us:9.1f}-{e.end_us:9.1f} us  {e.center_mhz:9.3f} MHz  "
              f"snr {e.snr_db}  {({k: v for k, v in e.features.items() if k in ('cfo_hz', 'zc_root', 'n_symbols', 'serial', 'clock_ppm')})}")
    if a.out:
        C.save_emissions(a.out, em)
