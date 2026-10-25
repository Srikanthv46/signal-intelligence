#!/usr/bin/env python3
"""Bench: per-carrier checks on a stored LO. For each frequency (MHz):
  raw peak frequency (no ppm correction; our own clock harmonics land exactly on n x 40 MHz / n x 61.44 MHz
  here, external signals are off by the receiver ppm), -10 dB width, power in the first 12 ms vs after 13 ms
  (capture-locked self noise), and FM: frequency deviation p5/p95 and audio-band energy.

    .venv/bin/python bench/live_lo/carrier_check.py 436 420.003 430.08 440.0 456.055
"""
import sys
from pathlib import Path
import numpy as np
from scipy import signal as ss
sys.path.insert(0, str(Path(__file__).resolve().parent))
from bursts import window

def main():
    tag = sys.argv[1]
    for f in map(float, sys.argv[2:]):
        rows = []
        for k in range(3):
            try:
                y, fs, f0 = window(tag, k, (f - 0.05) * 1e6, (f + 0.05) * 1e6)
            except (IndexError, FileNotFoundError):
                break
            fr, P = ss.welch(y, fs, nperseg=2048, return_onesided=False); fr = np.fft.fftshift(fr); P = 10 * np.log10(np.fft.fftshift(P))
            fl = np.median(P); i = np.argmax(P); w = fr[P > P[i] - 10]
            # fine peak with long FFT
            Y = np.abs(np.fft.fft(y * np.hanning(len(y)), 1 << 20)); fy = np.fft.fftfreq(1 << 20, 1 / fs); j = np.argmax(Y)
            pk = f0 + fy[j]
            z = y * np.exp(-2j * np.pi * fy[j] / fs * np.arange(len(y)))
            z = ss.lfilter(*ss.butter(4, 12.5e3 / (fs / 2)), z)
            n12, n13 = int(12e-3 * fs), int(13e-3 * fs)
            pa, pb = 10 * np.log10(np.mean(np.abs(z[:n12]) ** 2)), 10 * np.log10(np.mean(np.abs(z[n13:]) ** 2))
            inst = np.diff(np.unwrap(np.angle(z[n13:]))) * fs / (2 * np.pi)
            rows.append((pk, P[i] - fl, w.max() - w.min(), pa - pb, np.percentile(inst, 5), np.percentile(inst, 95)))
        print(f"{f}:")
        for pk, snr, wd, lock, p5, p95 in rows:
            print(f"   raw peak {pk/1e6:.6f} MHz  peak/floor {snr:5.1f} dB  -10dB width {wd/1e3:5.1f} kHz  "
                  f"first12ms-minus-after {lock:+5.1f} dB  inst-freq p5/p95 {p5/1e3:+6.2f}/{p95/1e3:+6.2f} kHz")

if __name__ == "__main__":
    main()
