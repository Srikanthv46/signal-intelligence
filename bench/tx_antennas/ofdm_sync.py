"""Blind CP-OFDM symbol extraction for one link (bench tool).

Mixes the link to baseband, resamples to N*SCS, finds the symbol grid of every burst from the CP correlation
and the fractional CFO from the CP phase, and returns the FFT of every symbol (fftshifted, bin N//2 = DC).
"""
import numpy as np
from fractions import Fraction
from scipy.signal import resample_poly


def load(path):
    import json
    meta = json.load(open(path.replace('.sigmf-data', '.sigmf-meta')))
    g = meta['global']
    dt = {'cf32_le': np.complex64, 'ci16_le': np.int16}[g['core:datatype']]
    x = np.fromfile(path.replace('.sigmf-meta', '.sigmf-data'), dtype=dt)
    if dt == np.int16:
        x = (x[0::2] + 1j * x[1::2]).astype(np.complex64)
    return x, g['core:sample_rate'], meta['captures'][0]['core:frequency']


def baseband(x, fs, f_off, fs2):
    n = np.arange(len(x))
    y = x * np.exp(-2j * np.pi * f_off / fs * n).astype(np.complex64)
    r = Fraction(fs2 / fs).limit_denominator(2000)
    return resample_poly(y, r.numerator, r.denominator).astype(np.complex128)


def cp_metric(y, N, cp):
    c = y[:-N] * np.conj(y[N:])
    cs = np.concatenate([[0], np.cumsum(c)])
    p = np.abs(y) ** 2
    ps = np.concatenate([[0], np.cumsum(p)])
    num = cs[cp:] - cs[:-cp]
    e1 = ps[cp:len(c) + 1] - ps[:len(c) + 1 - cp]
    e2 = ps[N + cp:] - ps[N:len(ps) - cp]
    m = len(num)
    den = np.sqrt(e1[:m] * e2[:m]) + 1e-12
    return num / den   # complex, |.| ~ 1 at CP starts of an OFDM symbol


def symbols(y, N, cp, thr=0.6):
    """Every symbol whose CP correlation is >= thr, following the grid; returns (starts, FFTs, cfo_cycles/sample)."""
    g = cp_metric(y, N, cp)
    a = np.abs(g)
    L = N + cp
    starts = []
    i = 0
    while i < len(a) - L:
        w = a[i:i + L]
        j = i + int(np.argmax(w))
        if a[j] >= thr:
            # follow the grid while the symbols continue
            k = j
            while k < len(a) and a[k] >= thr * 0.8:
                lo, hi = max(0, k - 3), min(len(a), k + 4)
                k = lo + int(np.argmax(a[lo:hi]))
                starts.append(k)
                k += L
            i = k
        else:
            i += L
    starts = np.array(sorted(set(starts)))
    ph = np.angle(np.sum(g[starts]))
    eps = -ph / (2 * np.pi * N)        # cycles per sample (fractional CFO)
    n = np.arange(len(y))
    yc = y * np.exp(-2j * np.pi * eps * n)
    X = np.array([np.fft.fftshift(np.fft.fft(yc[s + cp:s + L])) for s in starts if s + L <= len(y)])
    return starts[:len(X)], X, eps, a
