#!/usr/bin/env python3
"""Bench: time profile (per ms) and on-minus-off spectrum (4 MHz steps) of broadband energy in a stored LO snapshot.

    .venv/bin/python bench/live_lo/floor_profile.py 1192
"""
import json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
cap = json.loads((HERE / "out" / sys.argv[1] / "capture.json").read_text()); fs, fc = cap["rate"], cap["freq_mhz"] * 1e6
for k, sn in enumerate(cap["snapshots"]):
    z = np.load(sn["file"]); x = (z[:, 0] + 1j * z[:, 1]).astype(np.complex64); del z
    n = 1024; m = len(x) // n
    X = np.abs(np.fft.fftshift(np.fft.fft(x[:m * n].reshape(m, n) * np.hanning(n), axis=1), axes=1)) ** 2; del x
    f = fc + np.fft.fftshift(np.fft.fftfreq(n, 1 / fs))
    band = 10 * np.log10(X[:, 64:960].mean(1)); ms = band[:m // 60 * 60].reshape(-1, 60).mean(1)
    print(f"snap{k} per-ms dB:", " ".join(f"{v - ms.min():.0f}" for v in ms))
    s = np.argsort(band)[::-1][:m // 10]; q = np.argsort(band)[:m // 2]
    on = 10 * np.log10(X[s].mean(0)) - 10 * np.log10(X[q].mean(0))
    print("   on-off by 4 MHz:", " ".join(f"{f[i] / 1e6:.0f}:{on[i:i + 64].mean():.0f}" for i in range(64, 960, 64)))
