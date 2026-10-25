#!/usr/bin/env python3
"""Bench: run the Engine on testdata inputs and print the signals (optionally only some receivers).

    .venv/bin/python bench/run_engine.py testdata/live/LIVE_868MHz/LIVE_868MHz.sigmf-meta [--enable lte] [--ppm 7.4] [--all]
"""
import argparse, json, time
from spectrum_intel.api import Engine, IQBlock
from spectrum_intel.io.sigmf import load

ap = argparse.ArgumentParser(); ap.add_argument("paths", nargs="+"); ap.add_argument("--enable"); ap.add_argument("--ppm", type=float, default=7.4)
ap.add_argument("--all", action="store_true"); ap.add_argument("--bw", type=float)
a = ap.parse_args()
eng = Engine(enable=a.enable.split(",") if a.enable else None, detail="compact")
for p in a.paths:
    cap = load(p); t = time.time()
    r = eng.analyze(IQBlock(cap.iq, sample_rate_hz=cap.fs, center_frequency_hz=cap.fc, bandwidth_hz=a.bw or min(56e6, cap.fs), receiver_ppm=a.ppm))
    print(f"== {p.split('/')[-1]}  {time.time()-t:.1f}s  timing {json.dumps(r.get('timing_ms'))[:300]}")
    for s in r["signals"]:
        if a.all or s["decode"]["status"] != "not_possible" or s["technology"] not in ("unknown",):
            f = s["frequency"]
            print(f"   {f['center_hz']/1e6:10.4f} {(f.get('bandwidth_hz') or 0)/1e3:8.1f}k {s['technology']:14s} {s['label'][:70]:70s} [{s['decode']['status']}] {s.get('flags') or ''}")
