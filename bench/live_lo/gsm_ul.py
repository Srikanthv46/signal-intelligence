#!/usr/bin/env python3
"""Bench: GSM uplink normal bursts in a stored LO snapshot -> GMSK demod + training sequence (TSC = cell BCC).

    .venv/bin/python bench/live_lo/gsm_ul.py 904 899.6 [--snaps 0 1 2]
Per burst: start, the best of the 8 normal-burst TSCs (3GPP TS 45.002 5.2.3) and its bit errors, and the
burst-to-burst spacing in TDMA frames (4.615 ms). Differential decoding: sign of the GMSK frequency at symbol i
is d_i xor d_(i-1) (TS 45.004 2.4), so the TSC is matched in that form.
"""
import argparse, json
from pathlib import Path
import numpy as np
from scipy import signal as ss

HERE = Path(__file__).resolve().parent
RS = 1625e3 / 6                                    # 270.833 ksym/s
TSC = ["00100101110000100010010111", "00101101110111100010110111", "01000011101110100100001110",
       "01000111101101000100011110", "00011010111001000001101011", "01001110101100000100111010",
       "10100111110110001010011111", "11101111000100101110111100"]
FRAME = 4.615e-3


def diff(bits):
    b = np.array([int(c) for c in bits])
    return 1 - 2 * (b[1:] ^ b[:-1])               # +-1 per symbol, 25 values


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("tag"); ap.add_argument("f", type=float)
    ap.add_argument("--snaps", type=int, nargs="+", default=[0, 1, 2]); ap.add_argument("--ppm", type=float, default=7.4)
    a = ap.parse_args()
    cap = json.loads((HERE / "out" / a.tag / "capture.json").read_text())
    fs, fc = cap["rate"], cap["freq_mhz"] * 1e6
    ref = [diff(t) for t in TSC]
    for k in a.snaps:
        z = np.load(cap["snapshots"][k]["file"]); x = z[:, 0].astype(np.float32) + 1j * z[:, 1].astype(np.float32)
        f0 = a.f * 1e6 * (1 + a.ppm * 1e-6) - fc
        y = x * np.exp(-2j * np.pi * f0 / fs * np.arange(len(x))).astype(np.complex64)
        y = ss.resample_poly(y, 52, 2949)          # 61.44 MHz * 52/2949 = 1.0834 MHz = 4 sps
        sps = 4; fsr = fs * 52 / 2949
        p = np.abs(y) ** 2
        env = np.convolve(p, np.ones(sps * 140) / (sps * 140), "same")
        thr = 10 * np.median(p)
        on = env > thr
        edges = np.flatnonzero(np.diff(on.astype(int)))
        starts = edges[::2] + 1 if on[0] == 0 else np.r_[0, edges[1::2] + 1]
        print(f"== snap{k}  {len(starts)} burst candidates (fs {fsr / 1e3:.1f} kHz)")
        ph = np.angle(y[1:] * np.conj(y[:-1]))      # frequency per sample
        prev = None
        for s in starts:
            seg = slice(max(0, s - sps * 20), s + sps * 170)
            fr = ph[seg]
            if len(fr) < sps * 160: continue
            best = (0, None, None, None)
            for off in range(sps):                     # symbol timing
                sym = np.sign(np.convolve(fr, np.ones(sps), "same")[off::sps])
                for t, r in enumerate(ref):
                    c = np.correlate(sym, r, "valid")
                    i = int(np.argmax(c))                 # polarity is fixed by TS 45.004 (TSC0 inverted = TSC3 shifted)
                    if c[i] > best[0]: best = (c[i], t, i, off)
            c, t, i, off = best
            errs = int((25 - c) / 2)
            t_ms = (seg.start + i * sps + off) / fsr * 1e3
            gap = "" if prev is None else f"  +{(t_ms - prev) / (FRAME * 1e3):.2f} frames"
            if errs <= 4:
                prev = t_ms
            pos = (seg.start + i * sps + off - s) / sps   # symbols after the energy onset (normal burst: ~61)
            print(f"  t {t_ms:7.3f} ms  TSC {t}  errors {errs:2d}/25  pol {'+' if c > 0 else '-'}  at symbol {pos:5.1f}{gap}")


if __name__ == "__main__":
    main()
