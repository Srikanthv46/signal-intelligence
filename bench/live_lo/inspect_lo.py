#!/usr/bin/env python3
"""Bench: look at one LO's stored snapshots. Full-band spectrogram + zooms, and per-tone stats.

    .venv/bin/python bench/live_lo/inspect_lo.py 328 [--zoom 352:359 --zoom 306.9:307.5] [--snap 1]
"""
import argparse, json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import signal as ss

HERE = Path(__file__).resolve().parent

def load(tag, k):
    cap = json.loads((HERE / "out" / tag / "capture.json").read_text())
    a = np.load(cap["snapshots"][k]["file"])
    return a[:, 0].astype(np.float32) + 1j * a[:, 1].astype(np.float32), cap

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("freq"); ap.add_argument("--snap", type=int, default=0)
    ap.add_argument("--zoom", action="append", default=[]); ap.add_argument("--out")
    a = ap.parse_args()
    x, cap = load(a.freq, a.snap); fs, fc = cap["rate"], cap["freq_mhz"] * 1e6
    out = Path(a.out or HERE / "out" / a.freq)
    panels = 1 + len(a.zoom)
    fig, axs = plt.subplots(panels, 1, figsize=(14, 4 * panels), squeeze=False)
    f, t, S = ss.spectrogram(x, fs, nperseg=1024, noverlap=0, return_onesided=False)
    f = np.fft.fftshift(f); S = np.fft.fftshift(S, axes=0)
    ax = axs[0][0]; ax.pcolormesh(t * 1e3, (fc + f) / 1e6, 10 * np.log10(S + 1e-9), shading="auto", vmin=np.percentile(10*np.log10(S+1e-9), 5))
    ax.set_title(f"{a.freq} MHz snap{a.snap}"); ax.set_ylabel("MHz")
    for i, z in enumerate(a.zoom, 1):
        lo, hi = (float(v) * 1e6 for v in z.split(":")); fz = (lo + hi) / 2 - fc; bw = hi - lo
        n = np.arange(len(x)); y = x * np.exp(-2j * np.pi * fz / fs * n).astype(np.complex64)
        dec = max(1, int(fs // (bw * 1.25))); y = ss.decimate(y, dec, ftype="fir") if dec > 1 else y
        fz2, tz, Sz = ss.spectrogram(y, fs / dec, nperseg=512, noverlap=256, return_onesided=False)
        fz2 = np.fft.fftshift(fz2); Sz = 10 * np.log10(np.fft.fftshift(Sz, axes=0) + 1e-9)
        ax = axs[i][0]; ax.pcolormesh(tz * 1e3, (fc + fz + fz2) / 1e6, Sz, shading="auto", vmin=np.percentile(Sz, 10))
        ax.set_ylim(lo / 1e6, hi / 1e6); ax.set_ylabel("MHz")
    axs[-1][0].set_xlabel("ms"); fig.tight_layout(); p = out / f"inspect_snap{a.snap}.png"; fig.savefig(p, dpi=70); print(p)

if __name__ == "__main__":
    main()
