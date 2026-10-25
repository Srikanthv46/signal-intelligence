#!/usr/bin/env python3
"""Bench: does the LTE plugin find the cells if bursty wideband candidates are also searched?
Monkeypatches LTEDecoder.capture_signals in memory only (the library is not changed).

    .venv/bin/python bench/live_lo/lte_fix_test.py 868 [--snap 1] [--ppm 7.4]
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
from spectrum_intel.receivers.lte import plugin as P
from spectrum_intel.api import Engine, IQBlock

HERE = Path(__file__).resolve().parent

def patched(self, ctx):
    cand = [s for s in ctx.det.signals
            if not s.flags and s.snr_db >= 6 and s.bw >= 1.0e6 and self.can_decode(s, ctx)]   # no 'continuous' test
    if not cand:
        return []
    cal = self._ppm(ctx, cand); ppm = cal["ppm"] if cal else 0.0
    out, seen = [], set()
    for sig in cand:
        for c in P.search_cells(ctx, sig, ppm):
            key = (round(c["raster_hz"] / 100e3), c["nid2"])
            if key in seen:
                continue
            seen.add(key); out.append(self._signal(ctx, sig, c, c["sss"], cal))
    return out

def run(tag, k, ppm, patch):
    cls = next(v for v in vars(P).values() if isinstance(v, type) and hasattr(v, "capture_signals") and v.__module__ == P.__name__)
    orig = cls.capture_signals
    if patch:
        cls.capture_signals = patched
    try:
        if tag.endswith(".sigmf-meta"):
            from spectrum_intel.io.sigmf import load
            cap = load(tag); x, fs, fc, bw = cap.iq, cap.fs, cap.fc, 56e6
        else:
            c = json.loads((HERE / "out" / tag / "capture.json").read_text()); s = np.load(c["snapshots"][k]["file"])
            x = s[:, 0].astype(np.float32) + 1j * s[:, 1].astype(np.float32); fs, fc, bw = c["rate"], c["freq_mhz"] * 1e6, c["bw"]
        r = Engine(enable=["lte"], detail="compact").analyze(IQBlock(x, sample_rate_hz=fs, center_frequency_hz=fc,
                                                                     bandwidth_hz=bw, receiver_ppm=ppm))
    finally:
        cls.capture_signals = orig
    return [(s["label"], s["decode"]["status"], s.get("identifiers")) for s in r["signals"] if s["technology"] == "lte"]

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("tag"); ap.add_argument("--snap", type=int, default=1); ap.add_argument("--ppm", type=float, default=7.4)
    a = ap.parse_args()
    print("as is   :", run(a.tag, a.snap, a.ppm, False))
    print("patched :", run(a.tag, a.snap, a.ppm, True))

if __name__ == "__main__":
    main()
