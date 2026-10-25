#!/usr/bin/env python3
"""Reference for rc_link/models (research, PLAN M5): model blocks for hand-built emission sets.

    .venv/bin/python -I c++/tests/vectors/rc_models/gen_rc_models.py
"""
import json
import sys
from pathlib import Path

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parents[3] / "python"))
from spectrum_intel.core.emission import Emission                 # noqa: E402
from spectrum_intel.receivers.rc_link.models import model         # noqa: E402

CASES = [
    [("frsky_d16", 2440, dict(variant="D16 v2", crc_ok=True), {}), ("frsky_d16", 2441, dict(variant="D16", crc_ok=False), {})],
    [("frsky_d8", 2440, dict(variant="D8", crc_ok=False), {})],
    [("afhds2a_id", 2440, dict(crc_ok=None), {})],
    [("afhds2a", 2440, dict(crc_ok=True), {})],
    [("crossfire_sync", 868.3, dict(grid=dict(band="868", channel=3, link="uplink")), {})],
    [("crossfire_sync", 868.3, dict(sync="0x3E6B3E"), {})],
    [("sik_telemetry", 915, dict(air_rate_bps=64000), {})],
    [("sik_telemetry", 915, dict(air_rate_bps=125000), {})],
    [("elrs_lora_sf5_bw500k", 915.5, None, dict(modulation="lora", sf=5, bw_hz=500e3, elrs_channel=["FCC915", 7]))],
    [("elrs_lora_sf7_bw500k", 866.0, None, dict(modulation="lora", sf=7, bw_hz=500e3, elrs_channel=["IN866", 1]))],
    [("elrs_like_lora_sf6_bw812.5k", 2440.4, None, dict(modulation="lora", sf=6, bw_hz=812.5e3))],
    [("elrs_flrc_like", 2440, None, dict(modulation="gfsk"))],
    [("elrs_gfsk_like", 2440, None, dict(modulation="gfsk"))],
    [("elrs_gfsk_like", 915, None, dict(modulation="gfsk"))],
    [("nb_burst", 2440, None, {})],
]
ref = []
for c in CASES:
    em = [Emission("rc_fhss", l, 0, 100, f - 0.1, f + 0.1, f, 20, decoded=d, features=ft) for l, f, d, ft in c]
    ref.append(dict(em=[dict(label=l, center_mhz=f, decoded=d, features=ft) for l, f, d, ft in c], out=model(em)))
(OUT / "cases.json").write_text(json.dumps(ref, indent=1))
print(len(ref), sum(r["out"] is not None for r in ref))
