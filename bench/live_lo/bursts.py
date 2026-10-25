#!/usr/bin/env python3
"""Bench: find bursts in a frequency window of a stored LO snapshot and measure each one.

    .venv/bin/python bench/live_lo/bursts.py 328 351.5 356 [--snap 1] [--thr 10]
Per burst: start/duration, centre, -10 dB width, and the instantaneous-frequency spread (FSK test).
"""
import argparse, json
from pathlib import Path
import numpy as np
from scipy import signal as ss, ndimage

HERE = Path(__file__).resolve().parent

def window(tag, k, lo, hi):
    cap = json.loads((HERE / "out" / tag / "capture.json").read_text())
    a = np.load(cap["snapshots"][k]["file"]); x = a[:, 0].astype(np.float32) + 1j * a[:, 1].astype(np.float32)
    fs, fc = cap["rate"], cap["freq_mhz"] * 1e6
    f0 = (lo + hi) / 2 - fc; n = np.arange(len(x))
    y = x * np.exp(-2j * np.pi * f0 / fs * n).astype(np.complex64)
    dec = max(1, int(fs // ((hi - lo) * 1.25)))
    y = ss.resample_poly(y, 1, dec).astype(np.complex64) if dec > 1 else y
    return y, fs / dec, fc + f0

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("freq"); ap.add_argument("lo", type=float); ap.add_argument("hi", type=float)
    ap.add_argument("--snap", type=int, default=1); ap.add_argument("--thr", type=float, default=10); ap.add_argument("--nfft", type=int, default=256)
    a = ap.parse_args()
    y, fs, f0 = window(a.freq, a.snap, a.lo * 1e6, a.hi * 1e6)
    f, t, S = ss.stft(y, fs, nperseg=a.nfft, noverlap=a.nfft // 2, return_onesided=False)
    S = np.abs(np.fft.fftshift(S, axes=0)) ** 2; f = np.fft.fftshift(f)
    nf = np.median(S, axis=1, keepdims=True)
    m = 10 * np.log10(S / nf) > a.thr
    m = ndimage.binary_closing(m, np.ones((3, 3))); lab, k = ndimage.label(m)
    rows = []
    for sl in ndimage.find_objects(lab):
        fi, ti = sl
        if (ti.stop - ti.start) < 3: continue
        s0, s1 = int(t[ti.start] * fs), int(t[ti.stop - 1] * fs)
        seg_f = f[fi]; P = S[fi, ti].sum(1); fcen = float(np.sum(seg_f * P) / P.sum())
        w = seg_f[P > P.max() / 10]; bw = float(w.max() - w.min() + fs / a.nfft)
        z = y[s0:s1] * np.exp(-2j * np.pi * fcen / fs * np.arange(s1 - s0))
        z = ss.lfilter(*ss.butter(4, min(0.99, bw / fs * 1.5)), z)
        inst = np.diff(np.unwrap(np.angle(z))) * fs / (2 * np.pi); inst = inst[len(inst)//10: -len(inst)//10 or None]
        h, e = np.histogram(inst, bins=60, range=(-bw, bw))
        rows.append((s0 / fs * 1e3, (s1 - s0) / fs * 1e3, (f0 + fcen) / 1e6, bw / 1e3, np.std(inst) / 1e3,
                     np.percentile(inst, 10) / 1e3, np.percentile(inst, 90) / 1e3))
    print(f"fs={fs/1e6:.3f} MHz  {len(rows)} bursts\n t_ms   dur_ms  f_MHz      bw_kHz  instF_std  p10   p90 (kHz)")
    for r in sorted(rows):
        print("%7.2f %6.3f %10.4f %7.1f %8.1f %7.1f %7.1f" % r)

if __name__ == "__main__":
    main()
