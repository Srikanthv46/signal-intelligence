"""Synthetic AIS bursts (tests): payload -> HDLC frame -> NRZI -> GMSK (BT 0.4, h 0.5) at 9600 bit/s."""
from __future__ import annotations

import numpy as np

from . import ais as A

RATE = 9600.0


def gmsk(levels: list[int], fs: float, bt: float = 0.4) -> np.ndarray:
    """Phase-continuous GMSK baseband for NRZ levels (0/1 -> -1/+1), frequency deviation +-RATE/4."""
    sps = fs / RATE
    n = int(len(levels) * sps)
    nrz = (2.0 * np.array(levels, float) - 1)[(np.arange(n) / sps).astype(int)]
    sig = bt * RATE
    tt = np.arange(-int(1.5 * sps), int(1.5 * sps) + 1) / fs
    g = np.exp(-(2 * np.pi * sig * tt) ** 2 / (2 * np.log(2)))
    f = np.convolve(nrz, g / g.sum(), "same") * (RATE / 4)
    return np.exp(1j * 2 * np.pi * np.cumsum(f) / fs)


def burst(payload: str, fs: float, fill: int = 0) -> np.ndarray:
    """One AIS transmission (8 ramp bits of carrier, then the frame) from an NMEA armoured payload."""
    lv = A.nrzi(A.frame_bits(A.armor_to_bits(payload, fill)))
    lv = [lv[0]] * 8 + lv + [lv[-1]] * 8
    return gmsk(lv, fs)


def capture(fs: float, dur: float, fc: float, bursts: list[tuple[float, float, str]], snr_db: float = 20.0,
            seed: int = 0) -> np.ndarray:
    """bursts: (t0 s, carrier Hz, payload); SNR in a 25 kHz channel."""
    rng = np.random.default_rng(seed)
    n = int(dur * fs)
    x = np.zeros(n, np.complex128)
    for t0, f, p in bursts:
        b = burst(p, fs)
        i = int(t0 * fs)
        m = min(len(b), n - i)
        t = (i + np.arange(m)) / fs
        x[i:i + m] += b[:m] * np.exp(1j * (2 * np.pi * (f - fc) * t + rng.uniform(0, 6)))
    npow = 10 ** (-snr_db / 10) * fs / 25e3
    x += np.sqrt(npow / 2) * (rng.standard_normal(n) + 1j * rng.standard_normal(n))
    return x.astype(np.complex64)
