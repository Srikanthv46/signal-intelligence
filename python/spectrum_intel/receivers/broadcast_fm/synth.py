"""Synthetic FM broadcast transmitters (tests): MPX = audio + 19 kHz pilot + RDS on 57 kHz (3 x pilot), FM +-75 kHz."""
from __future__ import annotations

import numpy as np
from scipy import signal as ss

from . import rds as R


def rds_wave(groups: list[tuple[int, int, int, int]], fs: float, n: int, phase0: float = 0.0) -> np.ndarray:
    """RDS subcarrier (unit peak) for n samples: groups repeated, differential + biphase, low-passed, x cos(57 kHz)."""
    bits = []
    while len(bits) * fs / 1187.5 < n + 2 * fs / 1187.5:
        for g in groups:
            bits += R.encode_group(g)
    enc, prev = [], 0
    for b in bits:
        prev ^= b
        enc.append(prev)
    spb = fs / 1187.5
    t = np.arange(n)
    k = (t / spb).astype(int)
    first = (t - k * spb) < spb / 2
    sym = (2.0 * np.array(enc)[k] - 1) * np.where(first, 1.0, -1.0)
    sym = ss.sosfiltfilt(ss.butter(4, 2400, fs=fs, output="sos"), sym)
    return sym * np.cos(2 * np.pi * 57e3 * t / fs + phase0)


def fm_station(fs: float, dur: float, f_off: float, groups=None, stereo: bool = True, snr_db: float = 30.0,
               seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    fm = 480e3                                       # MPX rate
    n = int(dur * fm)
    t = np.arange(n) / fm
    a = ss.sosfiltfilt(ss.butter(4, [100, 15e3], "bandpass", fs=fm, output="sos"), rng.standard_normal(n))
    mpx = 0.8 * a / (np.percentile(np.abs(a), 99.5) + 1e-12)
    if stereo:
        mpx += 0.09 * np.sin(2 * np.pi * 19e3 * t)
    if groups:
        mpx += 0.05 * rds_wave(groups, fm, n)
    dev = 75e3 * mpx
    ni = int(dur * fs)
    di = np.interp(np.arange(ni) / fs, t, dev)
    x = np.exp(1j * (2 * np.pi * np.cumsum(f_off + di) / fs + rng.uniform(0, 6)))
    npow = 10 ** (-snr_db / 10) * fs / 200e3
    x = x + np.sqrt(npow / 2) * (rng.standard_normal(ni) + 1j * rng.standard_normal(ni))
    return x.astype(np.complex64)


def ps_groups(pi: int, pty: int, name: str) -> list[tuple[int, int, int, int]]:
    """Four 0A groups carrying an 8-character programme service name."""
    name = (name + " " * 8)[:8]
    out = []
    for seg in range(4):
        b = (0 << 12) | (0 << 11) | (pty << 5) | seg
        d = (ord(name[2 * seg]) << 8) | ord(name[2 * seg + 1])
        out.append((pi, b, 0xE0CD, d))
    return out
