#!/usr/bin/env python3
"""Bench: is a narrowband carrier TETRA (pi/4-DQPSK, 18 ksym/s, RRC 0.35, 25 kHz)? Uses a dwell capture.

    .venv/bin/python bench/live_lo/tetra_probe.py tetra390 390.0297 [391.527 ...]
Prints the PSD profile, the symbol-rate line, and the differential-phase histogram at the recovered timing.
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy import signal as ss
sys.path.insert(0, str(Path(__file__).resolve().parent))
from lo_step import SCRATCH

SR = 18000.0

def rrc(beta, sps, span=8):
    t = np.arange(-span * sps, span * sps + 1) / sps
    h = np.where(np.abs(1 - (4 * beta * t) ** 2) < 1e-8, np.pi / 4 * np.sinc(1 / (4 * beta)),
                 np.sinc(t) * np.cos(np.pi * beta * t) / (1 - (4 * beta * t) ** 2 + 1e-12))
    return h / np.sqrt(np.sum(h ** 2))

def load(name):
    m = json.loads((SCRATCH / "dwell" / f"{name}.json").read_text()); a = np.load(SCRATCH / "dwell" / f"{name}.npy")
    return a[:, 0].astype(np.float32) + 1j * a[:, 1].astype(np.float32), m

def baseband(x, m, f, sps=8):
    fs = m["rate"]; n = np.arange(len(x))
    y = x * np.exp(-2j * np.pi * (f * 1e6 - m["freq_mhz"] * 1e6) / fs * n).astype(np.complex64)
    from fractions import Fraction
    q = Fraction(int(SR * sps), int(fs)).limit_denominator(1000)
    return ss.resample_poly(y, q.numerator, q.denominator).astype(np.complex64), SR * sps

def main():
    name, freqs = sys.argv[1], [float(v) for v in sys.argv[2:]]
    x, m = load(name)
    for f in freqs:
        y, fs = baseband(x, m, f)
        fr, P = ss.welch(y, fs, nperseg=4096, return_onesided=False); fr = np.fft.fftshift(fr); P = 10 * np.log10(np.fft.fftshift(P))
        P -= np.median(P[np.abs(fr) > 40e3])
        prof = " ".join(f"{int(o/1e3):+d}:{P[np.argmin(np.abs(fr - o))]:.0f}" for o in np.arange(-25e3, 25.1e3, 2.5e3))
        # centre by spectral centroid of the >6 dB region, then matched filter
        sel = (P > 6) & (np.abs(fr) < 20e3); df = np.sum(fr[sel] * 10 ** (P[sel] / 10)) / np.sum(10 ** (P[sel] / 10)) if sel.any() else 0
        y = y * np.exp(-2j * np.pi * df / fs * np.arange(len(y)))
        z = np.convolve(y, rrc(0.35, 8), "same")
        e = np.abs(z) ** 2; E = np.abs(np.fft.rfft(e - e.mean())); fe = np.fft.rfftfreq(len(e), 1 / fs)
        s = (fe > 10e3) & (fe < 30e3); line = fe[s][np.argmax(E[s])]; snr_line = E[s].max() / np.median(E[s])
        # timing phase: maximise symbol-rate envelope, then differential phase
        best = None
        for ph in range(8):
            v = z[ph::8]; d = v[1:] * np.conj(v[:-1]); w = np.abs(d)
            q = np.angle(d) ; err = np.mean(w * np.abs(np.angle(np.exp(1j * 4 * (q - np.pi / 4))))) / np.mean(w)
            if best is None or err < best[0]: best = (err, ph, q, w)
        err, ph, q, w = best
        strong = w > np.percentile(w, 50)
        h, _ = np.histogram(np.degrees(q[strong]), bins=16, range=(-180, 180))
        print(f"\n{f} MHz: centre corr {df:+.0f} Hz; PSD dB rel floor (kHz:dB): {prof}")
        print(f"  envelope line {line:.1f} Hz ({snr_line:.0f}x median) ; pi/4 lattice err {err:.2f} rad (0 = perfect, ~0.79 = random)")
        print("  dphase hist (22.5 deg bins from -180):", h.tolist())

if __name__ == "__main__":
    main()
