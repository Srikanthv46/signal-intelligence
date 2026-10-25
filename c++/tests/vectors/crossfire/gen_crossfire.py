#!/usr/bin/env python3
"""Synthetic reference for TBS Crossfire 150 Hz (rc_link/crossfire; research, PLAN R9): a 7.68 MSps capture
(2^18 samples, fc 864 MHz) with three uplink packets on the 868 grid at 14 dB, rc_link detect's Crossfire emissions
and the plugin's crossfire field; plus grid_channel cases.

    nice -n 19 .venv/bin/python -I c++/tests/vectors/crossfire/gen_crossfire.py
"""
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parents[3] / "python"))
from spectrum_intel.receivers.rc_link import crossfire as CF, detector as D, synth as S   # noqa: E402
from spectrum_intel.receivers.rc_link.plugin import PLUGINS                             # noqa: E402

FS, FC = 7.68e6, 864e6
rng = np.random.default_rng(31)
x = S.noise_capture(1 << 18, rng)
data = bytes(rng.integers(0, 256, 23).tolist())
bursts = []
for k, n in enumerate((3, 9, 15)):
    iq, _ = S.crossfire_burst(data, cfo_hz=2e3, fs=FS, rng=rng)
    bursts.append(dict(t0_us=1500 + 6666 * k, f_hz=860.165e6 + n * 260e3, iq=iq, bw_hz=170e3))
y, _ = S.embed(x, FS, FC, bursts, 14)
y = y.astype(np.complex64)
y.tofile(OUT / "capture.cf32")
em = [e for e in D.detect(y, FS, FC) if e.label == "crossfire_sync"]
ref = dict(fs=FS, fc=FC, decoded=[e.decoded for e in em], field=PLUGINS[0].summarise(em, None)["fields"].get("crossfire"),
           grid=[[f, CF.grid_channel(f)] for f in (860.165e6 + 7 * 260e3, 860.165e6 + 57 * 260e3 + 5e3, 902.165e6 + 3 * 260e3,
                                                    863.093e6 + 10 * 112e3, 861.0e6)])
print(len(em), ref["field"])
(OUT / "cases.json").write_text(json.dumps(ref, indent=1, default=float))
