"""Synthetic FMCW radio altimeters (tests): a triangle (or sawtooth) sweep of bandwidth B at f_rep around a centre,
seen through a capture at fc; unit noise."""
from __future__ import annotations

import numpy as np


def fmcw(fs: float, dur: float, fc: float, centre: float = 4300e6, bw: float = 104e6, f_rep: float = 50.0,
         shape: str = "triangle", amp: float = 1.0, phase0: float = 0.1, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n = int(dur * fs)
    t = np.arange(n) / fs
    u = (f_rep * t + phase0) % 1.0
    tri = 2 * np.abs(u - 0.5) if shape == "triangle" else u          # 0..1
    f = centre - bw / 2 + bw * tri - fc
    x = amp * np.exp(1j * 2 * np.pi * np.cumsum(f) / fs)
    x *= np.abs(f) < 0.45 * fs                 # the receiver's analog filter: no aliased copies of the sweep
    x += np.sqrt(0.5) * (rng.standard_normal(n) + 1j * rng.standard_normal(n))
    return x.astype(np.complex64)
