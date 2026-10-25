#!/usr/bin/env python3
"""Synthetic reference for FrSky D8 / D16 (rc_link/cc2500; research, PLAN R6a): a 7.68 MSps capture (2^18 samples,
fc 2440 MHz) with one D8, D16 v1, D16 v2 and D16 LBT packet at 12 dB, rc_link detect's FrSky emissions (label,
decoded) and the plugin's frsky field; plus CRC cases.

    nice -n 19 .venv/bin/python -I c++/tests/vectors/frsky/gen_frsky.py
"""
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parents[3] / "python"))
from spectrum_intel.receivers.rc_link import cc2500 as CC, detector as D, synth as S   # noqa: E402
from spectrum_intel.receivers.rc_link.plugin import PLUGINS                          # noqa: E402

FS, FC = 7.68e6, 2440e6
LEN = {"D8": 17, "D16 v1": 27, "D16 v2": 27, "D16 LBT": 30}
rng = np.random.default_rng(41)
x = S.noise_capture(1 << 18, rng)
bursts = []
for k, (v, f) in enumerate(zip(LEN, (2437.9e6, 2439.4e6, 2441.6e6, 2442.7e6))):
    pl = bytes([0x34, 0x12]) + bytes(rng.integers(0, 256, LEN[v] - 2).tolist())
    iq, _ = S.frsky_burst(v, pl, cfo_hz=3e3, fs=FS, rng=rng)
    bursts.append(dict(t0_us=1000 + 8000 * k, f_hz=f, iq=iq, bw_hz=250e3))
y, _ = S.embed(x, FS, FC, bursts, 12)
y = y.astype(np.complex64)
y.tofile(OUT / "capture.cf32")
em = [e for e in D.detect(y, FS, FC) if e.label.startswith("frsky_")]
data = [bytes(rng.integers(0, 256, n).tolist()) for n in (0, 1, 18, 29, 33)]
ref = dict(fs=FS, fc=FC, labels=[e.label for e in em], decoded=[e.decoded for e in em],
           field=PLUGINS[0].summarise([e for e in em if e.decoded], None)["fields"].get("frsky"),
           crc=[[d.hex(), CC.crc16_cc2500(d), CC.crc16_d16(d)] for d in data])
print(ref["labels"], ref["field"])
(OUT / "cases.json").write_text(json.dumps(ref, indent=1, default=float))
