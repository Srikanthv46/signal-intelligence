#!/usr/bin/env python3
"""Bench: 4FSK (DMR, ETSI TS 102 361-1) demodulator + sync search on a dwell capture.

    .venv/bin/python bench/live_lo/dmr_probe.py d860 859.3455
Symbol rate from the FM-baseband transition line, 10 samples/symbol at 48 kHz, levels from the symbol histogram,
dibits: +3 -> 01, +1 -> 00, -1 -> 10, -3 -> 11 (TS 102 361-1 10.1.2), then the 48-bit SYNC patterns (9.1.1).
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy import signal as ss
sys.path.insert(0, str(Path(__file__).resolve().parent))
from lo_step import SCRATCH

SYNC = {"BS voice": 0x755FD7DF75F7, "BS data": 0xDFF57D75DF5D, "MS voice": 0x7F7D5DD57DFD, "MS data": 0xD5D7F77FD757,
        "RC": 0x77D55F7DFD77, "TDMA1 voice": 0x5D577F7757FF, "TDMA1 data": 0xF7FDD5DDFD55,
        "TDMA2 voice": 0x7DFFD5F55D5F, "TDMA2 data": 0xD7557F5FF7F5}

def bits48(v):
    return np.array([(v >> (47 - i)) & 1 for i in range(48)], np.int8)

def fm48k(name, f):
    m = json.loads((SCRATCH / "dwell" / f"{name}.json").read_text()); s = np.load(SCRATCH / "dwell" / f"{name}.npy")
    x = s[:, 0].astype(np.float32) + 1j * s[:, 1].astype(np.float32); fs = m["rate"]; fc = m["freq_mhz"] * 1e6
    y = x * np.exp(-2j * np.pi * (f * 1e6 - fc) / fs * np.arange(len(x))).astype(np.complex64); del x
    y = ss.resample_poly(y, 1, 50); fo = fs / 50                 # 48 kHz
    y = ss.lfilter(*ss.butter(5, 6.25e3 / (fo / 2)), y)
    d = np.angle(y[1:] * np.conj(y[:-1])) * fo / (2 * np.pi)
    return d - np.median(d), fo

def demod(d, fo, sr=4800.0):
    sps = fo / sr
    h = np.ones(int(round(sps * 0.6))) / int(round(sps * 0.6))   # integrate-and-dump-ish
    z = np.convolve(d, h, "same")
    # symbol rate line from |d'|^2
    e = np.diff(z) ** 2; E = np.abs(np.fft.rfft(e - e.mean(), 1 << 22)); fe = np.fft.rfftfreq(1 << 22, 1 / fo)
    s = (fe > sr - 30) & (fe < sr + 30); sr_m = fe[s][np.argmax(E[s])]
    sps = fo / sr_m
    best = None
    for ph in np.arange(0, sps, sps / 10):
        idx = (ph + sps * np.arange(int((len(z) - ph) / sps) - 1)).astype(int); v = z[idx]
        # spread around the 4 nearest levels after scaling by the outer level
        a = np.percentile(np.abs(v), 75); q = np.clip(np.round((v / a * 3 + 3) / 2) * 2 - 3, -3, 3)
        err = np.mean((v / a * 3 - q) ** 2)
        if best is None or err < best[0]: best = (err, v, a)
    err, v, a = best
    q = np.clip(np.round((v / a * 3 + 3) / 2) * 2 - 3, -3, 3).astype(int)
    m = {3: (0, 1), 1: (0, 0), -1: (1, 0), -3: (1, 1)}
    bits = np.array([b for s in q for b in m[s]], np.int8)
    return bits, sr_m, err

def main():
    name, f = sys.argv[1], float(sys.argv[2])
    d, fo = fm48k(name, f)
    bits, sr, err = demod(d, fo)
    print(f"{f} MHz: symbol rate {sr:.1f} sym/s, level-fit mse {err:.2f} (4 levels at ±1, ±3), {len(bits)} bits")
    b = 2 * bits.astype(np.float32) - 1
    for nm, v in SYNC.items():
        for sgn, lab in ((1, ""), (-1, " (inverted)")):
            c = np.correlate(b, sgn * (2 * bits48(v).astype(np.float32) - 1), "valid"); errs = (48 - c) / 2
            hits = np.where(errs <= 4)[0]
            if len(hits):
                print(f"  {nm}{lab}: {len(hits)} hits, spacing {np.diff(hits)[:8].tolist()} bits, errors {errs[hits][:8].astype(int).tolist()}")

if __name__ == "__main__":
    main()
