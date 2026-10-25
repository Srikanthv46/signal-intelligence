#!/usr/bin/env python3
"""Synthetic reference for the DroneID frame variants (dji detect_droneid; research, PLAN D4): two 15.36 MSps
captures (2^18 samples, fc 2440 MHz) with a 9-symbol and an 8-symbol DroneID burst at 4 dB, and the Python DroneID
emissions on them (decode off).

    nice -n 19 .venv/bin/python -I c++/tests/vectors/droneid_frames/gen_droneid_frames.py
"""
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parents[3] / "python"))
from spectrum_intel.receivers.dji import detector as D, synth as DS          # noqa: E402
from spectrum_intel.receivers.rc_link import synth as RS                     # noqa: E402

FS, FC = 15.36e6, 2440e6
ref = {"fs": FS, "fc": FC, "cases": []}
for ns, f_hz in ((9, 2441.3e6), (8, 2438.2e6)):
    rng = np.random.default_rng(20 + ns)
    x = RS.noise_capture(1 << 18, rng)
    y, _ = RS.embed(x, FS, FC, [dict(t0_us=4000, f_hz=f_hz, iq=DS.droneid_burst(ns, rng), bw_hz=9e6)], 4)
    y = y.astype(np.complex64)
    y.tofile(OUT / f"droneid_{ns}.cf32")
    em = [e for e in D.detect(y, FS, FC, ctx={"decode_droneid": False}) if e.family == "dji_droneid"]
    ref["cases"].append(dict(file=f"droneid_{ns}.cf32", n_symbols=ns,
                             emissions=[dict(start_us=e.start_us, end_us=e.end_us, center_mhz=float(e.center_mhz),
                                             features=e.features) for e in em]))
    print(ns, [(e.start_us, float(e.center_mhz), e.features["n_symbols"]) for e in em])
(OUT / "cases.json").write_text(json.dumps(ref, indent=1, default=float))
