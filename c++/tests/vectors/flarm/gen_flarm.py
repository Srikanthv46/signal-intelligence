#!/usr/bin/env python3
"""Synthetic reference for FLARM (rc_link/flarm; research, PLAN T5): a 7.68 MSps capture (2^18 samples, fc 868 MHz)
with two beacons (FLARM ID, ICAO address) at 14 dB; rc_link detect's FLARM emissions; plus standalone cases.

    nice -n 19 .venv/bin/python -I c++/tests/vectors/flarm/gen_flarm.py
"""
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parents[3] / "python"))
from spectrum_intel.receivers.rc_link import detector as D, flarm as F, synth as S   # noqa: E402

FS, FC = 7.68e6, 868.0e6
rng = np.random.default_rng(101)
x = S.noise_capture(1 << 18, rng)
bursts = []
for k, (f, rid, at) in enumerate(((868.2e6, 0x123456, 0x20), (868.4e6, 0x4CA1B2, 0x10))):
    iq, _ = S.flarm_burst(rid, at, cfo_hz=3e3, fs=FS, rng=rng)
    bursts.append(dict(t0_us=2000 + 12000 * k, f_hz=f, iq=iq, bw_hz=200e3))
y, _ = S.embed(x, FS, FC, bursts, 14)
y = y.astype(np.complex64)
y.tofile(OUT / "capture.cf32")
em = [e for e in D.detect(y, FS, FC) if e.label == "flarm_like"]
sa = []
for snr in (20, 6, 3):
    iq, _ = S.flarm_burst(0xABCDEF, fs=8e5, rng=rng)
    w = (iq + (rng.normal(size=len(iq)) + 1j * rng.normal(size=len(iq))) * np.sqrt(10 ** (-snr / 10) / 2)).astype(np.complex64)
    sa.append(dict(iq=[w.real.tolist(), w.imag.tolist()], out=F.try_identify(w, 8e5)))
ref = dict(fs=FS, fc=FC, labels=[e.label for e in em], decoded=[e.decoded for e in em], standalone=sa)
print(ref["labels"], [s["out"] and s["out"]["radio_id"] for s in sa])
(OUT / "cases.json").write_text(json.dumps(ref, default=float))
