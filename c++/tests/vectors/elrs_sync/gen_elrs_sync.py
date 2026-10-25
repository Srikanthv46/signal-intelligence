#!/usr/bin/env python3
"""Synthetic reference for ExpressLRS SYNC packets (rc_link/elrs + lora_phy.decode_implicit; research, PLAN R2c):
parse() on built packets, decode_implicit on encoded symbols, and rc_link detect on a 7.68 MSps capture (fc 915 MHz)
with SF8 and SF9 SYNC packets at 4 dB.

    nice -n 19 .venv/bin/python -I c++/tests/vectors/elrs_sync/gen_elrs_sync.py
"""
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parents[3] / "python"))
from spectrum_intel.receivers.rc_link import detector as D, elrs as E, lora_phy as P, synth as S   # noqa: E402

UID = bytes([0x9A, 0x01, 0x5C, 0x7E, 0x33, 0xC4])
ref = {"parse": [], "implicit": []}
for kw in (dict(ver=4), dict(ver=3), dict(ver=4, binding=True), dict(ver=4, ota8=True)):
    p = E.build_sync(12, 77, 3, 2, UID, **kw)
    ref["parse"].append(dict(hex=p.hex(), out=E.parse(p)))
bad = bytearray(E.build_sync(12, 77, 3, 2, UID)); bad[3] ^= 0x10
ref["parse"].append(dict(hex=bytes(bad).hex(), out=E.parse(bytes(bad))))
rng = np.random.default_rng(3)
for sf, cr in ((6, 3), (7, 3), (8, 3), (9, 4)):
    syms = S.lora_implicit_symbols(E.build_sync(5, 6, 7, 1, UID), sf, cr)
    syms = [int(v) for v in syms]
    syms[3] = (syms[3] + 1) % 2 ** sf                     # one symbol error (Hamming corrects)
    d = P.decode_implicit(syms, sf, cr, 8)
    ref["implicit"].append(dict(sf=sf, cr=cr, syms=syms, payload=d["payload"].hex(), corr=d["fec_corrected"],
                                bad=d["fec_detected_errors"], used=d["symbols_used"]))

fs, fc = 7.68e6, 915e6
x = S.noise_capture(1 << 18, rng)
bursts = []
for k, (mode, t0) in enumerate((("900_200hz", 1000), ("900_50hz", 7000))):     # both packets inside the 34 ms capture
    band, sf, bw, cr_den, pre, pl, hop, ival = D.ELRS_LORA[mode]
    iq, _ = S.lora_packet(sf, bw, n_pre=pre, payload=S.lora_implicit_symbols(E.build_sync(40 + k, 9, 2, 3, UID), sf, cr_den - 4),
                          sync_word=0x12, cfo_hz=4e3, fs=fs, rng=rng)
    bursts.append(dict(t0_us=t0, f_hz=903.5e6 + (926.9e6 - 903.5e6) * (20 + k) / 39, iq=iq, bw_hz=bw))
y, _ = S.embed(x, fs, fc, bursts, 4)
y = y.astype(np.complex64)
y.tofile(OUT / "capture.cf32")
ref["capture"] = dict(fs=fs, fc=fc, decoded=[e.decoded for e in D.detect(y, fs, fc) if (e.decoded or {}).get("protocol") == "ExpressLRS"])
print([(d["sync"]["fhss_index"], d["cr"]) for d in ref["capture"]["decoded"]])
(OUT / "cases.json").write_text(json.dumps(ref, indent=1, default=float))
