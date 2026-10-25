#!/usr/bin/env python3
"""Synthetic reference for mLRS (rc_link/mlrs; research, PLAN R3 / T6): two 7.68 MSps captures (2^19 samples): 915 MHz
with the decodable SX127x SF6 frame, 2440 MHz with an SF5 812.5 kHz frame (identified by air mode only); rc_link
detect's mLRS emissions; plus CRC / parse cases.

    nice -n 19 .venv/bin/python -I c++/tests/vectors/mlrs/gen_mlrs.py
"""
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parents[3] / "python"))
from spectrum_intel.receivers.rc_link import detector as D, mlrs as ML, synth as S   # noqa: E402

rng = np.random.default_rng(111)
caps = []
for name, fc, f, sf, bw in (("c915", 915e6, 913.8e6, 6, 500e3), ("c2440", 2440e6, 2441.0e6, 5, 812.5e3)):
    x = S.noise_capture(1 << 19, rng)
    fr = ML.build_frame(0xBEEF if sf == 6 else 0x1234, bytes(rng.integers(0, 256, 87).tolist()))
    iq, _ = S.lora_packet(sf, bw, n_pre=12, payload=S.lora_implicit_symbols(fr, sf, 1), sync_word=0x12, cfo_hz=3e3,
                          fs=7.68e6, rng=rng)
    y, _ = S.embed(x, 7.68e6, fc, [dict(t0_us=2000, f_hz=f, iq=iq, bw_hz=bw)], 10)
    y = y.astype(np.complex64)
    y.tofile(OUT / f"{name}.cf32")
    em = [e for e in D.detect(y, 7.68e6, fc) if e.label.startswith("mlrs")]
    caps.append(dict(file=f"{name}.cf32", fs=7.68e6, fc=fc, labels=[e.label for e in em], decoded=[e.decoded for e in em],
                     modes=[e.features.get("mlrs_modes") for e in em]))
frames = [ML.build_frame(0xBEEF, bytes(range(87))), bytes(91), bytes(rng.integers(0, 256, 91).tolist())]
ref = dict(captures=caps, parse=[[f.hex(), ML.parse(f)] for f in frames],
           crc=[[d.hex(), ML.crc16_mcrf4xx(d)] for d in (b"", b"123456789", bytes(range(89)))])
print([(c["labels"], c["modes"]) for c in caps])
(OUT / "cases.json").write_text(json.dumps(ref, default=float))
