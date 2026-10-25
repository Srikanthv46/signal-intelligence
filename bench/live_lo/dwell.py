#!/usr/bin/env python3
"""Bench: one long contiguous narrowband capture (a single Pluto buffer) at a chosen centre, for decoding work.
Same flush as lo_step: buffer destroyed before the retune, --discard buffers thrown away.

    .venv/bin/python bench/live_lo/dwell.py --freq 390.78 --rate 2.4e6 --bw 2e6 --n 4194304 --gain 50 --name tetra
-> <scratch>/dwell/<name>.npy (ci16) + .json
"""
import argparse, json, sys
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from lo_step import SCRATCH, FULL_SCALE, open_pluto

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--freq", type=float, required=True); ap.add_argument("--rate", type=float, default=2.4e6)
    ap.add_argument("--bw", type=float, default=2e6); ap.add_argument("--n", type=int, default=1 << 22)
    ap.add_argument("--gain", type=float, default=50); ap.add_argument("--discard", type=int, default=3)
    ap.add_argument("--name", required=True); ap.add_argument("--uri", default="usb:1.2.5")
    a = ap.parse_args()
    sdr = open_pluto(a.uri, a.rate, a.bw, a.n)
    try:
        sdr.rx_destroy_buffer(); sdr.rx_lo = int(a.freq * 1e6); sdr.rx_hardwaregain_chan0 = a.gain
        for _ in range(a.discard):
            sdr.rx()
        t0 = datetime.now(timezone.utc); x = sdr.rx()
    finally:
        sdr.rx_destroy_buffer(); del sdr
    clip = float(np.mean((np.abs(x.real) >= FULL_SCALE - 1) | (np.abs(x.imag) >= FULL_SCALE - 1)))
    d = SCRATCH / "dwell"; d.mkdir(parents=True, exist_ok=True)
    iq = np.empty((len(x), 2), np.int16); iq[:, 0], iq[:, 1] = x.real, x.imag; del x
    np.save(d / f"{a.name}.npy", iq)
    (d / f"{a.name}.json").write_text(json.dumps(dict(freq_mhz=a.freq, rate=a.rate, bw=a.bw, n=a.n, gain_db=a.gain,
                                                      utc=t0.isoformat(), clip=clip)))
    print(f"{d / a.name}.npy  {a.n / a.rate:.2f} s  clip {clip:.2e}  rms {np.sqrt(np.mean(iq.astype(np.float32) ** 2)) * 1.414:.1f}")

if __name__ == "__main__":
    main()
