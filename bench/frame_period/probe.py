"""D7: blind frame-period estimate on the Walksnail-like capture (testdata/live/LIVE_OFDM_VTX_5805MHz).
Records frame_period's inputs inside digital_vtx and compares estimators."""
import sys
import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "c++" / "tools"))
import parity                                                                  # noqa: E402
from spectrum_intel.receivers.digital_vtx import blind as BL                   # noqa: E402

rec = []
orig_fp, orig_seg = BL.frame_period, BL.segment


def fp(tk, zr, ts, **kw):
    r = orig_fp(tk, zr, ts, **kw)
    rec.append(dict(tk=tk.copy(), zr=zr.copy(), ts=ts, T=r[0], con=r[1]))
    return r


BL.frame_period = fp
name = sys.argv[1] if len(sys.argv) > 1 else "LIVE_OFDM_VTX_5805MHz"
inp = [i for i in parity.inputs() if i["name"] == name][0]
rep = parity.python_report(inp, ["digital_vtx"])
for s in rep["signals"]:
    f = (s.get("decoded") or {}).get("fields") or {}
    print(s.get("label"), s.get("center_mhz"), {k: f.get(k) for k in ("frame_period_ms", "frame_ms", "tu_us")})
np.savez(Path(__file__).with_name(f"rec_{name}.npz"), **{f"{k}{i}": np.asarray(v) for i, r in enumerate(rec)
                                                        for k, v in r.items()})
for r in rec:
    print("frame_period ->", r["T"], r["con"], "symbols", len(r["tk"]), "ts", r["ts"])
