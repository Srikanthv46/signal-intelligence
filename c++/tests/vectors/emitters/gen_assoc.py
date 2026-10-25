"""Reference for tests/unit/test_rx_emitters.cpp: the standalone emitters path on a DJI
RC-uplink capture (testdata/D_20260916-171343): extract_bursts -> measure_fingerprints
-> associate (python receivers/emitters analyse_capture).

Run: nice -n 19 .venv/bin/python c++/tests/vectors/emitters/gen_assoc.py"""
import json
from pathlib import Path

import numpy as np

from spectrum_intel.io import sigmf
from spectrum_intel.receivers.emitters import association as A, bursts as B

REPO = Path(__file__).resolve().parents[4]
cap = sigmf.load(REPO / "testdata/D_20260916-171343/D_20260916-171343.sigmf-meta")
em = B.extract_bursts(cap.iq, cap.fs, cap.fc)
A.measure_fingerprints(cap.iq, cap.fs, cap.fc, em)
emitters = A.associate(em, cap.fc)
llr01, ev01 = A.same_emitter_score(em[0], em[1])


def py(v):
    if isinstance(v, dict):
        return {k: py(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [py(x) for x in v]
    if isinstance(v, np.generic):
        return v.item()
    return v


out = dict(emissions=[dict(family=e.family, label=e.label, start_us=e.start_us, end_us=e.end_us,
                           center_mhz=e.center_mhz, snr_db=e.snr_db, features=py(e.features)) for e in em],
           emitters=[py(m.to_dict()) for m in emitters], pair01=dict(llr=llr01, evidence=py(ev01)))
Path(__file__).with_name("assoc_D_20260916-171343.json").write_text(json.dumps(out, indent=1))
print(len(em), "emissions,", len(emitters), "emitters")
