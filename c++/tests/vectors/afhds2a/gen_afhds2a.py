#!/usr/bin/env python3
"""Synthetic reference for FlySky AFHDS 2A (rc_link/afhds2a; research, PLAN R7): a 7.68 MSps capture (2^18 samples,
fc 2440 MHz) with two packets in the default candidate convention, one in another (classic Hamming, low nibble first,
CRC 0x8005 / 0xFFFF) and one ID-code-only packet at 14 dB; rc_link detect's AFHDS 2A emissions (label, decoded) and
the plugin's afhds2a field; plus CRC cases.

    nice -n 19 .venv/bin/python -I c++/tests/vectors/afhds2a/gen_afhds2a.py
"""
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parents[3] / "python"))
from spectrum_intel.receivers.rc_link import afhds2a as A, detector as D, synth as S   # noqa: E402
from spectrum_intel.receivers.rc_link.plugin import PLUGINS                          # noqa: E402

FS, FC = 7.68e6, 2440e6
rng = np.random.default_rng(61)
x = S.noise_capture(1 << 18, rng)


def pl(typ):
    return bytes([typ, 1, 2, 3, 4, 5, 6, 7, 8]) + bytes(rng.integers(0, 256, 29).tolist())


bursts = []
for k, (f, typ, conv) in enumerate(((2437.5e6, 0x58, {}), (2439.0e6, 0xAA, {}),
                                    (2441.5e6, 0x58, dict(layout="classic", high_first=False, crc=(0x8005, 0xFFFF))))):
    iq, _ = S.afhds2a_burst(pl(typ), cfo_hz=6e3, fs=FS, rng=rng, **conv)
    bursts.append(dict(t0_us=1000 + 3850 * k, f_hz=f, iq=iq, bw_hz=1.2e6))
bits = np.r_[np.arange(32) % 2, [(A.ID_CODE >> (31 - i)) & 1 for i in range(32)], rng.integers(0, 2, 560),
             np.zeros(8)].astype(np.uint8)
iq, _ = S.gfsk_burst(500e3, 186e3, bits=bits, bt=2.0, fs=FS, rng=rng)
bursts.append(dict(t0_us=20000, f_hz=2442.5e6, iq=iq, bw_hz=1.2e6))
y, _ = S.embed(x, FS, FC, bursts, 14)
y = y.astype(np.complex64)
y.tofile(OUT / "capture.cf32")
em = [e for e in D.detect(y, FS, FC) if e.label.startswith("afhds2a")]
data = [bytes(rng.integers(0, 256, n).tolist()) for n in (0, 1, 38)]
ref = dict(fs=FS, fc=FC, labels=[e.label for e in em], decoded=[e.decoded for e in em],
           field=PLUGINS[0].summarise(em, None)["fields"].get("afhds2a"),
           crc=[[d.hex(), p, i, A.crc16(d, p, i)] for d in data for p, i in A.CRC_VARIANTS])
print(ref["labels"], ref["field"])
(OUT / "cases.json").write_text(json.dumps(ref, indent=1, default=float))
