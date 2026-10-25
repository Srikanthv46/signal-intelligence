"""Reference vectors for c++/tests/unit/test_rx_gsm.cpp, written by the Python receiver
(spectrum_intel.receivers.gsm.sch) on a synthetic FCCH + SCH (MSK, 4 samples/symbol).

    .venv/bin/python c++/tests/vectors/gsm/gen_gsm.py   -> gsm.json next to this script
"""
import json
from pathlib import Path

import numpy as np

from spectrum_intel.receivers.gsm import sch as S

rng = np.random.default_rng(3)
cx = lambda a: [[float(v.real), float(v.imag)] for v in np.asarray(a).ravel()]


def crc10_par(b25):
    poly = 0b10101110101
    reg = list(b25) + [0] * 10
    for i in range(25):
        if reg[i]:
            for j in range(11):
                reg[i + j] ^= (poly >> (10 - j)) & 1
    return [r ^ 1 for r in reg[25:]]


def conv(u):
    s, out = 0, []
    for b in u:
        out += [int(v) for v in S._OUT[s, b]]
        s = S._NXT[s, b]
    return out


def msk(bits, sps=4):
    """phase steps of +-pi/2 per symbol (differential), linearly interpolated at sps."""
    d = [0] + [bits[k] ^ bits[k - 1] for k in range(1, len(bits))]
    th = np.cumsum([0.0] + [np.pi / 2 * (1 - 2 * x) for x in d[1:]])
    t = np.arange((len(bits) - 1) * sps + 1) / sps
    return np.exp(1j * np.interp(t, np.arange(len(bits)), th))


out = {"cases": []}
for case in range(3):
    info = [int(v) for v in rng.integers(0, 2, 25)]
    coded = conv(info + crc10_par(info) + [0, 0, 0, 0])
    burst = [0, 0, 0] + coded[:39] + [int(v) for v in S.TS_SCH] + coded[39:] + [0, 0, 0, 0, 0]
    n = 8000
    y = (rng.normal(size=n) + 1j * rng.normal(size=n)) * (0.1 + 0.15 * case)
    s0 = 900 + 37 * case
    y[s0:s0 + 148 * 4 + 1] += msk([0] * 149)
    st = s0 + 5000 + 3 * case
    seg = msk(burst)
    y[st:st + len(seg)] += seg
    cfo = (-400.0, 150.0, 1200.0)[case]
    y = (y * np.exp(2j * np.pi * cfo * np.arange(n) / S.OFS)).astype(np.complex64)
    fc = S.find_fcch(y)
    c = {"y": cx(y)}
    if fc is not None:
        a0, f, pk = fc
        y2 = y * np.exp(-2j * np.pi * f * np.arange(len(y)) / S.OFS)
        c.update(fcch=[int(a0), f, pk], y2=cx(y2))
        d = S.demod_sch(y2, a0 + int(round(1250 * S.SPS)))
        if d is not None:
            soft, inf = d
            bits = S.viterbi(soft)
            ok = S.crc10_ok(bits[:25], bits[25:35])
            c.update(soft=soft.tolist(), ts_match=inf["ts_match"], sample=inf["sample"], bits=[int(v) for v in bits],
                     crc_ok=ok, fields=S.sch_fields(bits[:25]) if ok else None)
    out["cases"].append(c)
soft = rng.normal(size=78)
out["vit_in"] = soft.tolist()
out["vit_out"] = [int(v) for v in S.viterbi(soft)]
b = [int(v) for v in rng.integers(0, 2, 25)]
out["fields_in"] = b
out["fields"] = S.sch_fields(b)
Path(__file__).with_name("gsm.json").write_text(json.dumps(out))
print("wrote gsm.json", [(c.get("fcch", [None])[0], c.get("crc_ok"), c.get("fields")) for c in out["cases"]])
