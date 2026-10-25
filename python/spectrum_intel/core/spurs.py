"""Receiver-generated lines and LO-harmonic images (shared by the pipeline and receivers).

Receiver spurs. A line generated inside the receiver from its own clocks is displayed at its exact nominal
frequency, because the LO and the sample clock come from the same reference:
  - harmonics of the sample clock, k x fs (Pluto: 430.08 MHz = 7 x 61.44 MHz);
  - lines locked to the reference at an integer MHz (Pluto: 536.000, 552.000 MHz). An external transmitter on an
    integer MHz is displayed off by the receiver's frequency error (+7.4 ppm -> ~4 kHz at 550 MHz), so this test
    needs the receiver ppm and only applies where that offset is many times the tolerance.

LO-harmonic images. A switching mixer also converts RF near n x LO (n odd) to baseband: a strong emitter at F
shows up at baseband d = 3 LO - F (inverted) or d = F - 5 LO. One block cannot prove it (a second block at
LO + delta can: the image moves by (1 -+ n) delta); harmonic_image_sources() only lists where a signal would come
from, for a caution.
Depends on numpy only.
"""
from __future__ import annotations

import numpy as np

LINE_TOL_HZ = 150.0                      # refined line vs nominal
REF_GRID_HZ = 1e6
STRONG_CLASSES = {"lte_dl", "nr_dl", "gsm_dl", "umts_dl", "lte_tdd", "nr_tdd"}


def refine_line(X: np.ndarray, fs: float, fc: float, f_abs: float, half_hz: float) -> tuple[float, float]:
    """Peak of a narrow line on the full-length FFT X of the capture (bins fs/L), parabolic on the log power.
    Returns (frequency in Hz, peak over the window median in dB)."""
    L = len(X)
    df = fs / L
    k0 = int(round((f_abs - fc) / df))
    h = max(2, int(np.ceil(half_hz / df)))
    k = np.arange(k0 - h, k0 + h + 1)
    k = k[(k > -L // 2) & (k < L // 2)]
    if len(k) < 5:
        return f_abs, 0.0
    p = np.abs(X[k % L]).astype(np.float64) ** 2
    i = int(np.argmax(p))
    off = 0.0
    if 0 < i < len(p) - 1:
        a, b, c = np.log(p[i - 1] + 1e-30), np.log(p[i] + 1e-30), np.log(p[i + 1] + 1e-30)
        den = a - 2 * b + c
        off = 0.5 * (a - c) / den if den < 0 else 0.0
    return fc + (k[i] + off) * df, float(10 * np.log10(p[i] / (np.median(p) + 1e-30)))


def sample_clock_harmonic(f_abs: float, fs: float, tol_hz: float = LINE_TOL_HZ) -> int | None:
    """k if the line is k x fs (k >= 2) within tol_hz."""
    k = int(round(f_abs / fs))
    return k if k >= 2 and abs(f_abs - k * fs) <= tol_hz else None


def reference_locked(f_abs: float, ppm: float | None, tol_hz: float = LINE_TOL_HZ) -> bool:
    """True if the line sits on the integer-MHz grid in the receiver's own (uncorrected) frequency scale while an
    external line there would be displayed >= 8 x tol_hz away because of the known receiver error."""
    if ppm is None or abs(ppm) * 1e-6 * f_abs < 8 * tol_hz:
        return False
    return abs(f_abs - round(f_abs / REF_GRID_HZ) * REF_GRID_HZ) <= tol_hz


LINE_MIN_DB = 10.0                       # a receiver line must stand this far over its +-window median


def receiver_lines(X: np.ndarray, fs: float, fc: float, ppm: float | None, band: tuple[float, float] | None = None,
                   win_hz: float = 5e3) -> list[dict]:
    """The receiver's own lines present in this capture: the sample-clock harmonic(s) k x fs in band and, with the
    receiver ppm known, lines locked to the reference on an integer MHz. X: full-length FFT of the capture.
    band: (lo, hi) analysed band (default fc +- 0.45 fs)."""
    lo, hi = band or (fc - 0.45 * fs, fc + 0.45 * fs)
    nominal = [(k * fs, "sample_clock_harmonic", k) for k in range(max(2, int(np.ceil(lo / fs))), int(hi // fs) + 1)]
    if ppm is not None:
        nominal += [(n * REF_GRID_HZ, "reference_locked", None) for n in range(int(np.ceil(lo / REF_GRID_HZ)),
                                                                              int(hi // REF_GRID_HZ) + 1)]
    out = []
    for f_nom, kind, k in nominal:
        if not lo <= f_nom <= hi or abs(f_nom - fc) < 3 * win_hz:          # the DC bin is LO leakage
            continue
        f, pk = refine_line(X, fs, fc, f_nom, win_hz)
        if pk < LINE_MIN_DB:
            continue
        ok = sample_clock_harmonic(f, fs) == k if kind == "sample_clock_harmonic" else reference_locked(f, ppm)
        if ok and abs(f - f_nom) <= LINE_TOL_HZ:
            out.append({"freq_hz": round(f, 1), "kind": kind, "harmonic": k, "peak_db": round(pk, 1)})
    return out


def harmonic_image_sources(f0: float, f1: float, lo: float, bands, harmonics=(3, 5)) -> list[dict]:
    """Strong-emitter bands an LO-harmonic image at [f0, f1] (tuned to lo) would come from.
    bands: [(f_lo, f_hi, key, description, classes), ...] (core.bands_india.BANDS)."""
    out = []
    for n in harmonics:
        d0, d1 = f0 - lo, f1 - lo
        src = sorted((n * lo - d0, n * lo - d1)) if n % 4 == 3 else sorted((n * lo + d0, n * lo + d1))
        for b_lo, b_hi, key, desc, cls in bands:
            if STRONG_CLASSES & set(cls) and src[0] >= b_lo and src[1] <= b_hi:
                out.append({"n": n, "inverted": n % 4 == 3, "source_low_hz": round(src[0]),
                            "source_high_hz": round(src[1]), "band": key, "band_description": desc})
                break
    return out


def in_strong_band(f0: float, f1: float, bands) -> bool:
    return any(STRONG_CLASSES & set(cls) and f0 < b_hi and f1 > b_lo for b_lo, b_hi, _, _, cls in bands)
