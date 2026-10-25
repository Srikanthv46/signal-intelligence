#!/usr/bin/env python3
"""Synthetic reference for Remote ID over Wi-Fi NAN (wifi/droneid.from_nan_action via mac80211.parse; research,
PLAN T3): the frame bytes and Python's remote_id / drone output.

    nice -n 19 .venv/bin/python -I c++/tests/vectors/wifi_nan/gen_wifi_nan.py
"""
import json
import sys
from pathlib import Path

OUT = Path(__file__).resolve().parent
for p in ("python", "python/tests/receivers"):
    sys.path.insert(0, str(OUT.parents[3] / p))
from test_det_droneid import nan_action_frame                     # noqa: E402
from spectrum_intel.receivers.wifi import mac80211 as M          # noqa: E402

f = nan_action_frame()
r = M.parse(f)
(OUT / "nan.json").write_text(json.dumps({"hex": f.hex(), "remote_id": r["remote_id"], "drone": r["drone"]}, indent=1))
print(r["remote_id"]["uas_id"], r["drone"]["source"])
