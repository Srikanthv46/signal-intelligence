"""Reference for tests/unit/test_core_ofdm15k.cpp: core.ofdm15k on a DJI RC-uplink capture
(testdata/D_20260916-171343), the way the emitters receiver calls it (spur notches, all
profiles, refine, merge, hop tracking).

Run: nice -n 19 .venv/bin/python c++/tests/vectors/core/gen_ofdm15k.py"""
import json
from pathlib import Path

import numpy as np

from spectrum_intel.core import ofdm15k as O15
from spectrum_intel.core import spectral as SP
from spectrum_intel.io import sigmf

REPO = Path(__file__).resolve().parents[4]
cap = sigmf.load(REPO / "testdata/D_20260916-171343/D_20260916-171343.sigmf-meta")
x, fs, fc = cap.iq, cap.fs, cap.fc
usable_half = min(fs, 56e6) / 2
spurs = SP.find_spurs(x, fs, fc)
ch = SP.Channelizer(x.astype(np.complex128), fs, [(s["offset_hz"], 30e3) for s in spurs])
hits, stats = O15.ofdm15k_scan(ch, fc, usable_half=usable_half)
refined = [O15.refine_ofdm15k(ch, h) for h in hits]
ob = O15.merge_ofdm15k([b for b in refined if b is not None])
tracked = O15.track_hops(ch, fc, ob, stats, len(x) / fs, usable_half=usable_half)
final = O15.merge_ofdm15k(ob + tracked)


def clean(d):
    return None if d is None else {k: (float(v) if isinstance(v, (np.floating, float)) else
                                       int(v) if isinstance(v, (np.integer, int)) and not isinstance(v, bool) else v)
                                   for k, v in d.items()}


out = dict(spurs=[s["offset_hz"] for s in spurs], stats=stats, hits=[clean(h) for h in hits],
           refined=[clean(b) for b in refined], merged=[clean(b) for b in ob],
           tracked=[clean(b) for b in tracked], final=[clean(b) for b in final],
           labels=[list(O15.label_ofdm15k(b)) for b in final])
p = Path(__file__).parent / "ofdm15k_D_20260916-171343.json"
p.write_text(json.dumps(out, indent=1))
print(p, len(hits), "hits,", len(ob), "merged,", len(tracked), "tracked,", len(final), "final")
