"""Synthetic GPS L1 C/A skies (research file 02 section 10): satellites with C/N0, Doppler, code phase and 50 bps
nav bits, a receiver clock error eps (LO and sample clock from one oscillator), complex noise of unit variance."""
from __future__ import annotations

import numpy as np

from .codes import CHIPS, F_L1, RATE, ca_code


def sky(fs: float, dur: float, fc: float, sats: list[dict], eps_ppm: float = 0.0, seed: int = 0,
        noise: bool = True) -> np.ndarray:
    """sats: [{prn, cn0, doppler, code_phase (chips), [phase]}]. Unit noise variance per sample: N0 = 1 / fs."""
    rng = np.random.default_rng(seed)
    n = int(dur * fs)
    eps = eps_ppm * 1e-6
    t = np.arange(n) / fs
    tt = t / (1 + eps)                                # true time of each sample (the sample clock runs fast)
    x = np.zeros(n, np.complex128)
    for s in sats:
        c = 1.0 - 2.0 * ca_code(s["prn"]).astype(np.float64)
        cp = s.get("code_phase", 0.0) + RATE * (1 + s["doppler"] / F_L1) * tt
        chips = c[np.floor(cp).astype(np.int64) % CHIPS]
        bit0 = rng.uniform(0, 0.02)
        bits = rng.choice([-1.0, 1.0], int(dur / 0.02) + 3)
        d = bits[np.floor((tt + bit0) / 0.02).astype(np.int64)]
        a = np.sqrt(10 ** (s["cn0"] / 10) / fs)
        f_bb = (F_L1 + s["doppler"]) - fc * (1 + eps)     # true carrier minus the (fast) LO, in true Hz
        ph = 2 * np.pi * f_bb * tt + s.get("phase", rng.uniform(0, 6))
        x += a * d * chips * np.exp(1j * ph)
    if noise:
        x += np.sqrt(0.5) * (rng.standard_normal(n) + 1j * rng.standard_normal(n))
    return x.astype(np.complex64)


def ppd(fs: float, dur: float, centre_hz: float, span_hz: float, period_s: float, amp: float = 3.0,
        ripple: float = 0.0, seed: int = 0) -> np.ndarray:
    """A sawtooth swept tone (GNSS personal privacy device) of the given span and period around centre_hz
    (baseband), with an optional relative period ripple; unit noise added."""
    rng = np.random.default_rng(seed)
    n = int(dur * fs)
    t = np.arange(n) / fs
    per = period_s * (1 + ripple * np.sin(2 * np.pi * t / 1e-3))
    ph = np.cumsum(1 / (per * fs)) % 1.0
    f = centre_hz + span_hz * (ph - 0.5)
    x = amp * np.exp(1j * 2 * np.pi * np.cumsum(f) / fs)
    x += np.sqrt(0.5) * (rng.standard_normal(n) + 1j * rng.standard_normal(n))
    return x.astype(np.complex64)
