"""Reference FFT vectors for tests/unit/test_core_fft.cpp (numpy 2.x / scipy semantics).

Run: .venv/bin/python c++/tests/vectors/core/gen_fft.py
Files are raw little-endian arrays; names are lower-case only (macOS file systems are
case-insensitive)."""
from pathlib import Path

import numpy as np
import scipy.fft as sf

out = Path(__file__).parent / "fft"
out.mkdir(exist_ok=True)
rng = np.random.default_rng(7)
for n in (1000, 1967):
    x = (rng.standard_normal(n) + 1j * rng.standard_normal(n)).astype(np.complex64)
    r = rng.standard_normal(n).astype(np.float32)
    d = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    arrays = {
        "in_c64": x, "in_f32": r, "in_c128": d,
        "sp_fft": sf.fft(x), "sp_ifft": sf.ifft(x), "sp_rfft": sf.rfft(r),
        "np_fft": np.fft.fft(x), "np_ifft": np.fft.ifft(x), "np_rfft": np.fft.rfft(r),
        "np_fft_pad": np.fft.fft(x, 2 * n), "fft_c128": np.fft.fft(d), "ifft_c128": np.fft.ifft(d),
    }
    for k, v in arrays.items():
        v.tofile(out / f"{k}_{n}.bin")
