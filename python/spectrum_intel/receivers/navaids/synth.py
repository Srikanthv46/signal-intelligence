"""Synthetic ILS (glide slope / localizer) and VOR carriers for tests (research file 05 section 2.4, 06 section 9)."""
from __future__ import annotations

import numpy as np


def _noise(x, snr_db, fs, rng):
    if snr_db is None:
        return x
    npow = 10 ** (-snr_db / 10) * fs / 25e3                 # SNR of the carrier in a 25 kHz channel
    return x + np.sqrt(npow / 2) * (rng.standard_normal(len(x)) + 1j * rng.standard_normal(len(x)))


def ils(fs: float, dur: float, f_off: float, m90: float = 0.4, m150: float = 0.4, scale: float = 1.0,
        phi: float = 0.3, snr_db: float | None = 30.0, two_freq_hz: float | None = None, seed: int = 0) -> np.ndarray:
    """Carrier AM by 90 s and 150 s Hz (phase-locked to 30 Hz) at f_off; two_freq_hz adds a clearance carrier
    (-6 dB, same tones) that far above."""
    rng = np.random.default_rng(seed)
    t = np.arange(int(dur * fs)) / fs
    env = 1 + m90 * np.cos(2 * np.pi * 90 * scale * t + phi) + m150 * np.cos(2 * np.pi * 150 * scale * t + 5 / 3 * phi)
    x = env * np.exp(1j * (2 * np.pi * f_off * t + rng.uniform(0, 6)))
    if two_freq_hz:
        x = x + 0.5 * env * np.exp(1j * (2 * np.pi * (f_off + two_freq_hz) * t + rng.uniform(0, 6)))
    return _noise(x, snr_db, fs, rng).astype(np.complex64)


def vor(fs: float, dur: float, f_off: float, radial_deg: float, am: float = 0.3, sub: float = 0.3,
        dev_hz: float = 480.0, snr_db: float | None = 30.0, seed: int = 0) -> np.ndarray:
    """Conventional VOR: 30 Hz AM lagging the 30 Hz FM (9960 Hz subcarrier, +-dev_hz) by the radial."""
    rng = np.random.default_rng(seed)
    t = np.arange(int(dur * fs)) / fs
    th = np.radians(radial_deg)
    env = (1 + am * np.cos(2 * np.pi * 30 * t - th)
           + sub * np.cos(2 * np.pi * 9960 * t + dev_hz / 30 * np.sin(2 * np.pi * 30 * t)))
    x = env * np.exp(1j * (2 * np.pi * f_off * t + rng.uniform(0, 6)))
    return _noise(x, snr_db, fs, rng).astype(np.complex64)
