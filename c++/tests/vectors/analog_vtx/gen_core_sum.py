#!/usr/bin/env python3
"""Vectors for c++/tests/unit/test_core_npcompat_analog.cpp: np.sum / np.mean of float32,
float64 and complex64 arrays longer than numpy's 8192-element reduction buffer."""
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
SIZES = [7, 100, 8192, 8193, 30000]

r = np.random.default_rng(11)
with open(OUT / "core_sum.bin", "wb") as f:
    for n in SIZES:
        x = r.standard_normal(n) * 1000 + 3000
        z = ((r.standard_normal(n) + 1j * r.standard_normal(n)) * 1000).astype(np.complex64)
        x32 = x.astype(np.float32)
        x32.tofile(f); x.tofile(f); z.tofile(f)
        np.array([np.sum(x32), np.mean(x32)], np.float32).tofile(f)
        np.array([np.sum(x), np.mean(x)], np.float64).tofile(f)
        np.array([np.sum(z)], np.complex64).tofile(f)
