"""
synth.py -- synthetic DJI DroneID bursts for tests (research, PLAN D4; research/drone_transmitters/
09_ble_coded_droneid_legacy.md). OFDM 15.36 MSps, FFT 1024, 600 used subcarriers + empty DC, 15 kHz.
  9 symbols (OcuSync 2/3: Mini 2, Air 2 and later): CP 80, 72 x 7, 80; ZC roots 600 / 147 in symbols 4 / 6;
  8 symbols (legacy: Mavic 2 and older): the first symbol dropped (CP 80, 72 x 6, 80 [U]); ZC in symbols 3 / 5.
Payload symbols carry random QPSK (no turbo code / scrambler: detection tests only).
"""
from __future__ import annotations

import numpy as np

from .detector import zc_freq

N, L = 1024, 601
FS = 15.36e6
LAYOUT = {9: ([80, 72, 72, 72, 72, 72, 72, 72, 80], (3, 5)), 8: ([80, 72, 72, 72, 72, 72, 72, 80], (2, 4))}


def droneid_burst(n_sym=9, rng=None):
    rng = rng or np.random.default_rng()
    cps, (z1, z2) = LAYOUT[n_sym]
    out = []
    for i, cp in enumerate(cps):
        X = np.zeros(N, complex)
        if i == z1:
            used = zc_freq(600, L)
        elif i == z2:
            used = zc_freq(147, L)
        else:
            used = (rng.choice([-1, 1], L) + 1j * rng.choice([-1, 1], L)) / np.sqrt(2)
            used[L // 2] = 0
        X[N // 2 - L // 2:N // 2 + L // 2 + 1] = used
        s = np.fft.ifft(np.fft.ifftshift(X)) * np.sqrt(N)
        out.append(np.r_[s[-cp:], s])
    return np.concatenate(out).astype(np.complex64)


def upsample(x, k):
    """FFT interpolation by an integer factor k (15.36 -> 61.44 MSps: k = 4)."""
    n = len(x)
    X = np.fft.fft(x)
    Y = np.zeros(n * k, complex)
    h = n // 2
    Y[:h], Y[-(n - h):] = X[:h], X[h:]
    return (np.fft.ifft(Y) * k).astype(np.complex64)
