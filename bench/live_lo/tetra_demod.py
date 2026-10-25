#!/usr/bin/env python3
"""Bench: TETRA downlink pi/4-DQPSK demodulator + training-sequence search (ETSI EN 300 392-2 cl. 5, 9).

    .venv/bin/python bench/live_lo/tetra_demod.py tetra390 390.0297
"""
import sys
from pathlib import Path
import numpy as np
from scipy import signal as ss
sys.path.insert(0, str(Path(__file__).resolve().parent))
from tetra_probe import load, rrc, SR
from fractions import Fraction

# EN 300 392-2 9.4.4.3: training sequences (bits)
Y_SYNC = np.array([1,1,0,0,0,0,0,1,1,0,0,1,1,1,0,0,1,1,1,0,1,0,0,1,1,1,0,0,0,0,0,1,1,0,0,1,1,1], np.int8)
N_TS = np.array([1,1,0,1,0,0,0,0,1,1,1,0,1,0,0,1,1,1,0,1,0,0], np.int8)
P_TS = np.array([0,1,1,1,1,0,1,0,0,1,0,0,0,0,1,1,0,1,1,1,1,0], np.int8)
SPS = 8

def symbol_rate(y, fs):
    e = np.abs(np.convolve(y, rrc(0.35, SPS), "same")) ** 2
    nf = 1 << 24; E = np.abs(np.fft.rfft(e - e.mean(), nf)); fe = np.fft.rfftfreq(nf, 1 / fs)
    s = (fe > SR - 50) & (fe < SR + 50); return float(fe[s][np.argmax(E[s])])

def demod(x, m, f):
    fs0 = m["rate"]; n = np.arange(len(x))
    y = x * np.exp(-2j * np.pi * (f * 1e6 - m["freq_mhz"] * 1e6) / fs0 * n).astype(np.complex64)
    q = Fraction(int(SR * SPS), int(fs0)).limit_denominator(1000)
    y = ss.resample_poly(y, q.numerator, q.denominator).astype(np.complex64); fs = SR * SPS
    sr = symbol_rate(y, fs)                              # true symbol rate in our (ppm-off) sample clock
    q2 = Fraction(sr / SR).limit_denominator(20000)      # resample so that exactly SPS samples/symbol
    y = ss.resample_poly(y, q2.denominator, q2.numerator).astype(np.complex64)
    z = np.convolve(y, rrc(0.35, SPS), "same")
    # coarse frequency from the 4th power of the 1-symbol differential
    d = z[SPS:] * np.conj(z[:-SPS]); rot = (np.angle(np.mean(d ** 4)) - np.pi) / 4
    z = z * np.exp(-1j * rot / SPS * np.arange(len(z)))
    out = []
    blk = 2000 * SPS
    for b0 in range(0, len(z) - blk // 4, blk):           # per block: timing phase + residual rotation
        seg = z[b0:b0 + blk + SPS]
        best = max(range(SPS), key=lambda p: np.mean(np.abs(seg[p::SPS]) ** 2))
        v = seg[best::SPS][:blk // SPS + 1]
        dd = v[1:] * np.conj(v[:-1]); r = (np.angle(np.mean(dd ** 4)) - np.pi) / 4
        out.append(dd * np.exp(-1j * r))
    dd = np.concatenate(out)
    ph = np.angle(dd)
    # dibit map (EN 300 392-2 5.5.2): 00:+pi/4  01:+3pi/4  10:-pi/4  11:-3pi/4
    b1 = (ph < 0).astype(np.int8); b2 = (np.abs(ph) > np.pi / 2).astype(np.int8)
    bits = np.empty(2 * len(ph), np.int8); bits[0::2], bits[1::2] = b1, b2
    evm = np.mean(np.abs(np.angle(np.exp(1j * 4 * (ph - np.pi / 4))))) / 4
    return bits, dd, sr, rot * SR / (2 * np.pi), evm

def search(bits, seq, maxerr):
    L = len(seq); s = 2 * seq - 1; b = 2 * bits.astype(np.float32) - 1
    c = np.correlate(b, s.astype(np.float32), "valid"); err = (L - c) / 2
    hits = np.where(err <= maxerr)[0]
    return [(int(h), int(err[h])) for h in hits]

def main():
    name, f = sys.argv[1], float(sys.argv[2])
    x, m = load(name)
    bits, dd, sr, foff, evm = demod(x, m, f)
    print(f"{f} MHz: symbol rate {sr:.2f} (ppm {((sr / SR) - 1) * 1e6:+.1f}), residual freq {foff:+.0f} Hz, "
          f"phase err {np.degrees(evm):.1f} deg rms-ish, {len(bits)} bits")
    for nm, seq, me in (("SYNC y (38)", Y_SYNC, 4), ("normal n (22)", N_TS, 2), ("normal p (22)", P_TS, 2)):
        h = search(bits, seq, me)
        print(f"  {nm}: {len(h)} hits", [(p, e, p % 510) for p, e in h[:12]])

if __name__ == "__main__":
    main()
