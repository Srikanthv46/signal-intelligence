#!/usr/bin/env python3
"""Bench: compare two LOs' mean spectra on their overlap, aligned by absolute frequency and by LO offset.
Air signals agree in absolute frequency; instrument artefacts (LO/clock/IQ image/filter) agree in offset.

    .venv/bin/python bench/live_lo/overlap.py 328 364
"""
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent

def mean_psd(tag, m=4096):
    c = json.loads((HERE / "out" / tag / "capture.json").read_text()); acc = 0
    for s in c["snapshots"]:
        a = np.load(s["file"]); x = a[:, 0].astype(np.float32) + 1j * a[:, 1]
        acc = acc + np.mean(np.abs(np.fft.fftshift(np.fft.fft(x[:len(x)//m*m].reshape(-1, m) * np.hanning(m), axis=1), axes=1)) ** 2, 0)
    fs, fc = c["rate"], c["freq_mhz"] * 1e6
    return (np.arange(m) - m // 2) * fs / m, fc, 10 * np.log10(acc / len(c["snapshots"]))

def main():
    (fa, ca, pa), (fb, cb, pb) = mean_psd(sys.argv[1]), mean_psd(sys.argv[2])
    pa -= np.median(pa); pb -= np.median(pb)
    lo, hi = max(ca + fa[0], cb + fb[0]) + 3e6, min(ca + fa[-1], cb + fb[-1]) - 3e6
    grid = np.arange(lo, hi, 0.5e6)
    ia = np.interp(grid, ca + fa, pa); ib = np.interp(grid, cb + fb, pb)
    print("absolute-aligned corr %.2f   (ca-%d offset grid)" % (np.corrcoef(ia, ib)[0, 1], 0))
    off = np.arange(-25e6, 25e6, 0.5e6)
    oa = np.interp(off, fa, pa); ob = np.interp(off, fb, pb)
    print("offset-aligned  corr %.2f" % np.corrcoef(oa, ob)[0, 1])
    print("offset MHz  : " + " ".join("%5.0f" % (o / 1e6) for o in off[::4]))
    print(f"{sys.argv[1]:>10} : " + " ".join("%5.1f" % v for v in oa[::4]))
    print(f"{sys.argv[2]:>10} : " + " ".join("%5.1f" % v for v in ob[::4]))

if __name__ == "__main__":
    main()
