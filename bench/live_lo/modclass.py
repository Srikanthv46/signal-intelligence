#!/usr/bin/env python3
"""Bench: IQ-only characterisation of one signal in a stored LO snapshot (no images).

Channelises the signal, then measures from the samples:
  occupied bandwidth (99 % power), carrier/CW component, envelope statistics (constant-envelope or not),
  duty cycle, instantaneous-frequency spread, cyclic features: symbol-rate lines of |x|^2 and of the
  differential, the 2nd/4th-power carrier lines (BPSK / QPSK / FSK), and OFDM cyclic-prefix lags.

    .venv/bin/python bench/live_lo/modclass.py 508 533.335 --bw 0.6 [--snap 1]
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
from scipy import signal as ss
sys.path.insert(0, str(Path(__file__).resolve().parent))
from bursts import window

def lines(v, fs, lo, hi, k=4):
    V = np.abs(np.fft.fft(v - v.mean())); f = np.fft.fftfreq(len(v), 1 / fs)
    s = (np.abs(f) >= lo) & (np.abs(f) <= hi); Vs = V[s] / np.median(V[s]); i = np.argsort(Vs)[-k:][::-1]
    return [(round(float(f[s][j]), 1), round(float(Vs[j]), 1)) for j in i]

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("tag"); ap.add_argument("f", type=float)
    ap.add_argument("--bw", type=float, required=True, help="MHz analysis span"); ap.add_argument("--snap", type=int, default=1)
    a = ap.parse_args()
    y, fs, f0 = window(a.tag, a.snap, (a.f - a.bw / 2) * 1e6, (a.f + a.bw / 2) * 1e6)
    y = y[int(13e-3 * fs):]                                   # after the capture self-noise step
    f, P = ss.welch(y, fs, nperseg=4096, return_onesided=False); f = np.fft.fftshift(f); P = np.fft.fftshift(P)
    nf = np.median(P); Ps = np.clip(P - nf, 0, None); c = np.cumsum(Ps) / Ps.sum()
    occ = f[np.searchsorted(c, 0.995)] - f[np.searchsorted(c, 0.005)]
    pk = np.argmax(P); cw = 10 * np.log10(P[pk] * fs / 4096 / (Ps.sum() * fs / 4096 + 1e-12))
    snr = 10 * np.log10(Ps.sum() / (nf * np.sum(Ps > 0) + 1e-12))
    z = y * np.exp(-2j * np.pi * f[pk] / fs * np.arange(len(y)))
    z = ss.lfilter(*ss.butter(5, min(0.95, max(occ, 5e3) / fs)), z)
    env = np.abs(z); papr = 10 * np.log10(env.max() ** 2 / np.mean(env ** 2)); kurt = np.mean(env ** 4) / np.mean(env ** 2) ** 2
    b = int(1e-4 * fs); eb = env[:len(env) // b * b].reshape(-1, b).mean(1); duty = np.mean(eb > 2 * np.median(eb))
    inst = np.diff(np.unwrap(np.angle(z))) * fs / (2 * np.pi)
    print(f"{a.f} MHz (fs {fs/1e3:.0f} kHz): peak at {(f0 + f[pk])/1e6:.6f} MHz, 99% occ {occ/1e3:.1f} kHz, "
          f"SNR in occ {snr:.1f} dB, carrier share {cw:.1f} dB")
    print(f"  envelope: PAPR {papr:.1f} dB, kurtosis {kurt:.2f} (1 = const env, 2 = Gaussian/OFDM), high-duty {duty:.2f}")
    print(f"  inst freq: std {np.std(inst)/1e3:.2f} kHz, p5/p95 {np.percentile(inst,5)/1e3:+.2f}/{np.percentile(inst,95)/1e3:+.2f} kHz")
    print(f"  |x|^2 lines (symbol/frame rate): {lines(env ** 2, fs, 50, fs / 2)}")
    for p in (2, 4, 8):
        print(f"  x^{p} lines (carrier x{p}): {lines(z ** p / (np.abs(z) ** (p - 1) + 1e-9), fs, 0, fs / 2, 3)}")
    print(f"  FM-baseband lines: {lines(inst, fs, 20, min(fs / 2, 50e3))}")
    for lag_us in (3.2, 4.0, 16, 64, 66.67, 133.3, 224, 896):
        L = int(round(lag_us * 1e-6 * fs))
        if 0 < L < len(z) // 4:
            r = np.abs(np.vdot(z[:-L], z[L:])) / np.vdot(z, z).real
            if r > 0.05: print(f"  correlation at lag {lag_us} us: {r:.3f}")

if __name__ == "__main__":
    main()
