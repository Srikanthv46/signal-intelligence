"""Barrage-noise and noise-FM jammers (research, signal gaps S12; research/signal_gaps/05_aviation_radar_jammers.md
section 5.4): a wide, continuous energy box tested for what a noise jammer is and no service is.

    analyse(ch, f0, f1) -> dict      measurements + "class": "barrage" | "noise_fm" | None

1. Plateau (block FFT, the box's inner 90 % in 64 sub-bands): >= 90 % of sub-bands within +-3 dB of their median;
   each measurable outside strip (5-15 % of the width beyond an edge) >= 6 dB below it; no narrow line (16-bin
   averages: max / median < 6 dB).
2. Statistics on the box at 1.25 x its width (inner 90 %; the FFT-based tests on its first 2^17 samples):
   envelope k = E|x|^4 / E^2|x|^2 (Gaussian 2, constant envelope 1), I and Q kurtosis (Gaussian 3).
3. No structure:
   - autocorrelation of the Hann-tapered spectrum, |R(tau)| / R(0) over tau = max(1 us, 8 / BW) ... 1000 us, below
     max(0.025, 5 / sqrt(N)) (an OFDM cyclic prefix gives its CP fraction: LTE at 0 dB ~0.035);
   - no cyclic line in the FFT of |x|^2, of x^2, and (constant envelope) of the instantaneous frequency: each bin
     over (ln(N) + 9) x its 512-bin block's mean (P_fa ~ 1e-4 per spectrum);
   - power per 1 ms within +-1 dB over the whole block (no TDD gaps, no gating).
4. Barrage: Gaussian (k 1.85-2.15, I / Q kurtosis 2.7-3.3) and every test of 1 and 3, and >= 10 MHz wide unless its
   edges match an allocation (one WCDMA carrier is just as Gaussian). Noise-FM: k <= 1.2, instantaneous frequency
   spread >= 15 % of the width with a Gaussian-like kurtosis (2.2-4.5: FSK is two-level, a chirp uniform) and the
   tests of 3. covers_band: an allocation (core.bands_india) whose edges are within 2 MHz of the box's.
Depends on numpy and spectrum_intel.core only.
"""
from __future__ import annotations

import numpy as np

from spectrum_intel.core.bands_india import BANDS

MIN_BW = 1e6
N_SUB = 64
FLAT_DB = 3.0
FLAT_FRAC = 0.9
EDGE_DB = 6.0
LINE_DB = 6.0
MAX_N = 2 ** 17
CP_MIN = 0.025
CP_LAG_BW = 8.0
CYC_MARGIN = 9.0
CYC_BLOCK = 512
STAT_DB = 1.0
MIN_JAM_BW = 10e6
MATCH_HZ = 2e6


def _db(x: float) -> float:
    return float(10 * np.log10(x + 1e-30))


def plateau(ch, f0: float, f1: float) -> dict:
    bw = f1 - f0
    sel = ch._bins((f0 + f1) / 2 - ch.fc, 0.45 * bw)
    P = np.abs(ch.X[sel]).astype(np.float64) ** 2
    sub = np.array([float(np.mean(b)) for b in np.array_split(P, N_SUB)])
    med = float(np.median(sub))
    flat = float(np.mean(np.abs(10 * np.log10(sub / med + 1e-30)) <= FLAT_DB))
    edges = []
    for a, b in ((f0 - 0.15 * bw, f0 - 0.05 * bw), (f1 + 0.05 * bw, f1 + 0.15 * bw)):
        if a < ch.fc - 0.45 * ch.fs or b > ch.fc + 0.45 * ch.fs:
            continue
        k = ch._bins((a + b) / 2 - ch.fc, (b - a) / 2)
        edges.append(round(_db(med / float(np.mean(np.abs(ch.X[k]).astype(np.float64) ** 2))), 1))
    m = len(P) // 16
    sm = P[:m * 16].reshape(m, 16).mean(axis=1) if m else P
    line = _db(float(sm.max()) / float(np.median(sm)))
    return {"flatness": round(flat, 3), "edge_db": edges, "line_db": round(line, 1),
            "plateau_ok": bool(flat >= FLAT_FRAC and edges and min(edges) >= EDGE_DB and line < LINE_DB)}


def _cyclic_db(v: np.ndarray, fs: float, f_lo: float, f_hi: float) -> float:
    """Peak of |FFT(v - mean)|^2 over |f| in [f_lo, f_hi], each bin against its 512-bin block's mean (the block
    median / ln 2 of exponential bins), in units of ln(N) + 9: > 0 dB = a line (P_fa ~ 1e-4 per spectrum). The local
    reference keeps a shaped spectrum (a filtered IF, the triangle of |x|^2) from reading as a line."""
    S = np.fft.fftshift(np.abs(np.fft.fft(v - np.mean(v))) ** 2)
    fr = np.abs(np.fft.fftshift(np.fft.fftfreq(len(v), 1 / fs)))
    s = S[(fr >= f_lo) & (fr <= f_hi)]
    m = len(s) // CYC_BLOCK
    if m < 1:
        return -99.0
    blk = s[:m * CYC_BLOCK].reshape(m, CYC_BLOCK)
    ratio = blk.max(axis=1) / (np.median(blk, axis=1) / np.log(2))
    return _db(float(ratio.max()) / (np.log(m * CYC_BLOCK) + CYC_MARGIN))


def covers_band(f0: float, f1: float) -> list[str]:
    return [key for lo, hi, key, *_ in BANDS if abs(f0 - lo) <= MATCH_HZ and abs(f1 - hi) <= MATCH_HZ]


def analyse(ch, f0: float, f1: float) -> dict | None:
    bw = f1 - f0
    if bw < MIN_BW:
        return None
    out = {"plateau_bw_mhz": round(bw / 1e6, 2), **plateau(ch, f0, f1)}
    yall, fs = ch.slice((f0 + f1) / 2, 0.45 * bw, min(1.25 * bw, ch.fs))
    if len(yall) < 4096:
        return None
    pall = np.abs(yall).astype(np.float64) ** 2
    y = yall[:MAX_N].astype(np.complex128)          # the FFT tests on the first 2^17 samples (latency)
    n = len(y)
    p = np.abs(y) ** 2
    k = float(np.mean(p ** 2) / np.mean(p) ** 2)
    ki = float(np.mean(y.real ** 4) / np.mean(y.real ** 2) ** 2)
    kq = float(np.mean(y.imag ** 4) / np.mean(y.imag ** 2) ** 2)
    out.update(envelope_k=round(k, 3), kurtosis_i=round(ki, 3), kurtosis_q=round(kq, 3))
    Y = np.fft.fft(y, 2 * n)
    fy = np.fft.fftfreq(2 * n, 1 / fs)
    w = np.where(np.abs(fy) < 0.45 * bw, np.cos(np.pi * fy / (0.9 * bw)) ** 2, 0.0)   # Hann taper: no sinc sidelobes
    R = np.fft.ifft(np.abs(Y) ** 2 * w)[:n]
    r = np.abs(R) / np.abs(R[0])
    l0 = int(np.ceil(max(1e-6, CP_LAG_BW / bw) * fs))   # past the tapered main lobe
    l1 = min(int(1000e-6 * fs), n // 4)
    if l1 > l0:
        j = l0 + int(np.argmax(r[l0:l1]))
        out.update(cp_corr=round(float(r[j]), 4), cp_lag_us=round(j / fs * 1e6, 2))
    else:
        out.update(cp_corr=0.0, cp_lag_us=None)
    cp_thr = max(CP_MIN, 5 / np.sqrt(n))
    out["cyclic_env_db"] = round(_cyclic_db(p, fs, 1e3, bw), 1)
    out["cyclic_conj_db"] = round(_cyclic_db(y ** 2, fs, 0.0, fs / 2), 1)
    ms = int(1e-3 * fs)
    nb = len(pall) // ms
    if nb >= 4:                                     # stationarity over the whole block
        pw = pall[:nb * ms].reshape(nb, ms).mean(axis=1)
        out["stationarity_db"] = round(float(np.max(np.abs(10 * np.log10(pw / np.mean(pw))))), 2)
    else:
        out["stationarity_db"] = 0.0
    structure_free = (out["cp_corr"] < cp_thr and out["cyclic_env_db"] < 0 and out["cyclic_conj_db"] < 0
                      and out["stationarity_db"] <= STAT_DB)
    cb = covers_band(f0, f1)
    out["covers_band"] = cb
    cls = None
    if 1.85 <= k <= 2.15 and 2.7 <= ki <= 3.3 and 2.7 <= kq <= 3.3 and structure_free and out["plateau_ok"] \
            and (bw >= MIN_JAM_BW or cb):
        cls = "barrage"
    elif k <= 1.2:
        d = np.angle(y[1:] * np.conj(y[:-1])) * (fs / (2 * np.pi))
        sd = float(np.std(d))
        kd = float(np.mean((d - d.mean()) ** 4) / np.var(d) ** 2)
        out.update(if_std_mhz=round(sd / 1e6, 3), if_kurtosis=round(kd, 3),
                   cyclic_if_db=round(_cyclic_db(d, fs, 1e3, fs / 2), 1))
        if sd >= 0.15 * bw and 2.2 <= kd <= 4.5 and out["cyclic_if_db"] < 0 and out["cp_corr"] < cp_thr \
                and out["stationarity_db"] <= STAT_DB:
            cls = "noise_fm"
    out["class"] = cls
    return out
