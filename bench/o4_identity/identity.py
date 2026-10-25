#!/usr/bin/env python3
"""Bench: crystal fingerprint and reference-sequence alignment of DJI O4 links in real recordings.

    cd bench/o4_identity && nice -n 19 ../../.venv/bin/python identity.py links_2026-10-16.json

The spec lists links: {"name", "path", "fc", "f0" (approximate centre), "bw" (18e6 or 9e6), "a", "b" (a contiguous
stretch, no recorder discontinuity inside)}. The first link is the alignment reference. Per link: clock ppm
(timing drift of recurring reference symbols), the integer bin n and nominal centre that make the carrier ppm equal
the clock ppm on the 0.5 MHz raster, the carrier ppm, and, independently, the n implied by aligning its pilot
sequence on the reference link's (20 MHz: n_ref - shift; 10 MHz: n_ref - 150 - shift).
"""
import json
import sys
from pathlib import Path

import numpy as np

import o4

SH = 200


def main(spec_path):
    spec = json.loads(Path(spec_path).read_text())
    rows = []
    for s in spec["links"]:
        x = o4.load_cf32(s["path"])
        L = o4.analyse_link(x, s["fc"], s["f0"], s["bw"], s["a"], s["b"], SH)
        del x
        sco = L["sco"]
        rc = o4.resolve_carrier(L["f0"], L["cfo"], sco["ppm"]) if sco["ppm"] is not None else None
        car, ppm = o4.carrier_ppm(L["f0"], L["cfo"], rc["n"], rc["f_nom"]) if rc else (None, None)
        rows.append((s, L, rc, ppm))
        print(f"{s['name']:34s} grid {L['grid_quality']:4.1f}  on {L['n_on']:3d}  CFO {L['cfo']:+8.0f} Hz  "
              + (f"clock {sco['ppm']:+.3f} +- {sco['sigma']:.3f} ppm (n={sco['n']})  -> n {rc['n']:+d}, nominal "
                 f"{rc['f_nom'] / 1e6:.1f} MHz (residual {rc['raster_residual_hz']:+5.0f} Hz), carrier {ppm:+.3f} ppm"
                 if rc else "clock: too few recurring reference symbols"), flush=True)
    s0, L0, rc0, _ = rows[0]
    print(f"\nsequence alignment on {s0['name']} (n_ref = {rc0['n'] if rc0 else '?'}):")
    for s, L, rc, ppm in rows[1:]:
        b, sh = o4.match(L["rows"], L["K"], L0["rows"], L0["K"], SH)
        st = o4.match_stats(b, sh, SH)
        slice_off = 150 if (s["bw"] < 12e6 and s0["bw"] > 12e6) else 0
        n_seq = (rc0["n"] - slice_off - st["shift"]) if rc0 else None
        agree = "" if (n_seq is None or rc is None or st["p"] > 1e-6) else \
            ("  n agrees" if n_seq == rc["n"] else f"  n DISAGREES (clock says {rc['n']:+d})")
        print(f"  {s['name']:34s} shift {st['shift']:+4d}  x{st['hits']}/{st['n']}  p {st['p']:.1e}"
              + (f"  -> n {n_seq:+d}{agree}" if st["p"] < 1e-6 and n_seq is not None else "  (no common sequence)"))


if __name__ == "__main__":
    main(sys.argv[1])
