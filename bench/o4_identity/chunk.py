"""Bench: crystal ppm of a DJI O4 link from one short chunk of an extracted packet.

Input is what signal-vault's ``extract_packet`` produces: baseband IQ of one burst at the packet rate, plus the
absolute mixing frequency (capture centre + the DDC's integer offset (lo_n + hi_n) // 2).

1. Resample the chunk to an O4-native rate R (15.36 / 30.72 / 61.44 MSps, the first >= packet rate), so a symbol is
   an integer FFT: N = R / 30 kHz, CP = 144 R / 61.44 MHz.
2. Symbols one by one from the CP metric peaks (no slot fold: works on 3 symbols and across recorder jumps).
3. Fine CFO from those CPs (+-15 kHz).
4. Integer subcarrier n: correlate the chunk's comb-2 differential spectra with a stored O4 template (reference
   symbols of a known link, n_ref) over integer shifts; 10 MHz mode carries the slice centred 150 bins below the
   20 MHz centre.
5. carrier = mix + CFO + 30 kHz x n; ppm against the 0.5 MHz raster.
"""
from __future__ import annotations

from fractions import Fraction

import numpy as np
from scipy.signal import find_peaks, resample_poly

import o4

SCS = 30e3
NATIVE = (15.36e6, 30.72e6, 61.44e6)


def native_rate(fs):
    return next((r for r in NATIVE if r >= fs - 1), NATIVE[-1])


def to_native(iq, fs):
    R = native_rate(fs)
    if abs(R - fs) < 1:
        return iq, R
    q = Fraction(R / fs).limit_denominator(20000)
    return resample_poly(iq, q.numerator, q.denominator), fs * q.numerator / q.denominator


def diff_spec(X, N, K):
    """d(k) = X(k+2) X*(k), k in [-K, K), unit modulus; X is fftshifted (bin 0 at N//2)."""
    v = X[N // 2 - K: N // 2 + K + 2]
    d = v[2:] * np.conj(v[:-2])
    return d / (np.abs(d) + 1e-30)


def build_template(x, fc, mix, bw, a, b, K=400):
    """Template from a long, clean stretch of a known 20 MHz link at 61.44 MSps (n_ref = 0 relative to mix + cfo)."""
    y = o4.baseband(x, fc, mix, bw + 2e6)
    g = o4.slot_grid(y, a, b)
    X, ok = o4.symbol_spectra(y, g["starts"], g["cfo"])
    sel = g["on"] & ok
    rows = np.array([diff_spec(X[i], o4.N, K) for i in np.flatnonzero(sel)])
    return dict(rows=rows, K=K, n_ref=0, mix=mix, cfo=g["cfo"])


def chunk_symbols(y, R, rho_min=0.3):
    """Symbol spectra of a native-rate chunk: CP-metric peaks, fine CFO, FFT inside each CP."""
    r = R / 61.44e6
    N, CP = int(round(2048 * r)), int(round(144 * r)); SYM = int(round(2192 * r))
    if len(y) < N + CP + 2:
        return None
    p = y[:-N] * np.conj(y[N:]); c = np.cumsum(np.concatenate([[0], p])); P = c[CP:] - c[:-CP]
    e = np.abs(y[:-N]) ** 2 + np.abs(y[N:]) ** 2; ce = np.cumsum(np.concatenate([[0], e])); E = (ce[CP:] - ce[:-CP]) / 2
    rho = np.abs(P) / (E + 1e-30)
    back = CP // 4
    pk, _ = find_peaks(rho, height=rho_min, distance=int(0.85 * SYM))
    pk = pk[pk + CP - back + N <= len(y)]
    if len(pk) == 0:
        return None
    cfo = float(-np.angle(np.sum(P[pk])) * R / (2 * np.pi * N))
    yc = y * np.exp(-2j * np.pi * cfo * np.arange(len(y)) / R)
    X = np.array([np.fft.fftshift(np.fft.fft(yc[s + CP - back: s + CP - back + N])) for s in pk])
    return dict(X=X, N=N, cfo=cfo, n_sym=len(pk), rho=rho[pk])


def match_template(D, tmpl, Kc):
    """Score S(s) = sum over chunk symbols of the best |corr| with any template row at shift s (chunk bin k <->
    template bin k + s). Returns best s, its score and a robust z against the other shifts."""
    T = tmpl["rows"]; Kt = tmpl["K"]
    L = 1 << int(np.ceil(np.log2(2 * Kt + 2 * Kc)))
    FT = np.fft.fft(T, L, axis=1)
    S = None
    for d in D:
        r = np.abs(np.fft.ifft(FT * np.conj(np.fft.fft(d, L)), axis=1))[:, :2 * Kt - 2 * Kc + 1] / (2 * Kc)
        best = r.max(0)
        S = best if S is None else S + best
    s_vals = np.arange(len(S)) + Kc - Kt
    i = int(np.argmax(S))
    med = np.median(S); mad = np.median(np.abs(S - med)) + 1e-12
    return int(s_vals[i]), float(S[i] / len(D)), float((S[i] - med) / (1.4826 * mad))


def chunk_ppm(iq, fs, mix, bw_hz, tmpl, raster=0.5e6, z_min=8.0):
    """Crystal ppm from one chunk. bw_hz: the packet's box width (selects 10/20 MHz mode and the usable bins)."""
    y, R = to_native(np.asarray(iq, complex), fs)
    sy = chunk_symbols(y, R)
    if sy is None:
        return dict(ok=False, why="no symbols")
    mode10 = bw_hz < 12e6
    Kc = min(142 if mode10 else 292, int(min(fs, bw_hz) / 2 / SCS) - 4, sy["N"] // 2 - 4)
    D = [diff_spec(X, sy["N"], Kc) for X in sy["X"]]
    s, score, z = match_template(D, tmpl, Kc)
    n = tmpl["n_ref"] - (150 if mode10 else 0) - s
    car = mix + sy["cfo"] + SCS * n
    nom = raster * round(car / raster)
    return dict(ok=z >= z_min, why="" if z >= z_min else f"weak template match (z {z:.1f})", ppm=(car - nom) / nom * 1e6,
                carrier=car, nominal=nom, n=n, shift=s, z=z, score=score, cfo=sy["cfo"], n_sym=sy["n_sym"], rate=R)
