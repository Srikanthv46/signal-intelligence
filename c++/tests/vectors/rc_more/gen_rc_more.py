#!/usr/bin/env python3
"""Synthetic reference for the PLAN G3-G6 RC decoders (rc_link radiolink / afhds2a AFHDS 1 / hubsan / toys / elrs GFSK;
research): a 7.68 MSps capture (2^18 samples, fc 2440 MHz) with a RadioLink air packet, two AFHDS 1 packets, a Hubsan
H107 packet, an MJXQ (XN297) and a V2x2 (nRF24) toy packet and an ExpressLRS GFSK SYNC packet at 15 dB; rc_link detect's
emissions of those labels (label, decoded), the plugin fields and model blocks, plus hop-table cases.

    nice -n 19 .venv/bin/python -I c++/tests/vectors/rc_more/gen_rc_more.py
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parents[3] / "python"))
from spectrum_intel.receivers.rc_link import afhds2a as A, detector as D, elrs as E, hubsan as H   # noqa: E402
from spectrum_intel.receivers.rc_link import radiolink as RL, synth as S, toys as T               # noqa: E402
from spectrum_intel.receivers.rc_link.plugin import PLUGINS                                     # noqa: E402

FS, FC = 7.68e6, 2440e6
LABELS = ("radiolink", "afhds1", "hubsan", "toy_drone", "elrs_gfsk")
rng = np.random.default_rng(71)
x = S.noise_capture(1 << 18, rng)
bursts = []
tx = bytes.fromhex("fc110d20")
hf, _ = RL.hop_table(tx)
ps = 0
while True:                                             # a hop slot whose channel is within +-2.5 MHz of fc
    ps = (ps * 0xAA + 3) % 0x7673
    if abs(RL.F0 + hf[ps & 15] * RL.STEP - FC) < 2.5e6:
        break
iq, _ = S.radiolink_burst(RL.data_packet(0x23, tx, list(range(900, 1060, 10)), ps), cfo_hz=4e3, fs=FS, rng=rng)
bursts.append(dict(t0_us=500, f_hz=RL.F0 + hf[ps & 15] * RL.STEP, iq=iq, bw_hz=150e3))
pl = A.afhds1_packet(0x5A6B7C3D, [1500, 1480, 1100, 1520, 1000, 2000, 1500, 1500])
hops = [h for h in A.afhds1_hops(pl[1]) if abs(A.A7105_F0 + h * 0.5e6 - FC) < 2.5e6]
for k in range(2):
    iq, _ = S.a7105_burst(A.ID_CODE, pl, crc=(0x1021, 0x0000), cfo_hz=6e3, fs=FS, rng=rng)
    bursts.append(dict(t0_us=11000 + 1510 * k, f_hz=A.A7105_F0 + hops[k % len(hops)] * 0.5e6, iq=iq, bw_hz=1.0e6))
pk = H.packet(0x20, bytes([0, 0x80, 0, 0x7F, 0, 0x80, 0, 0x80, 0x02, 0x64]) + bytes.fromhex("db042679"))
iq, _ = S.a7105_burst(0x3C2D1E0F, pk, rate=100e3, cfo_hz=5e3, fs=FS, rng=rng)
bursts.append(dict(t0_us=14500, f_hz=H.F0 + 0x50 * 0.5e6, iq=iq, bw_hz=0.6e6))
mj = bytearray([0x80, 0x00, 0x80, 0x80, 0x40, 0x40, 0x40, 0xA1, 0xB2, 0x00, 0, 0, 0, 0, 0x04])
mj = bytes(mj) + bytes([sum(mj) & 0xFF])
iq, _ = S.xn297_burst(bytes.fromhex("6d6a777777"), mj, cfo_hz=8e3, fs=FS, rng=rng)
bursts.append(dict(t0_us=17500, f_hz=2438e6, iq=iq, bw_hz=1.2e6))
v2 = bytearray([0x80, 0x00, 0x80, 0x80, 0x40, 0x40, 0x40, 0x12, 0x34, 0x56, 0, 0, 0, 0, 0x00])
v2 = bytes(v2) + bytes([sum(v2) & 0xFF])
bits = np.r_[T.nrf_packet_bits(bytes.fromhex("6688686868"), v2), np.zeros(8, np.uint8)]
iq, _ = S.gfsk_burst(1e6, 160e3, bits=bits, bt=0.5, pre_bits=0, cfo_hz=-8e3, fs=FS, rng=rng)
bursts.append(dict(t0_us=19500, f_hz=2442e6, iq=iq, bw_hz=1.2e6))
uid = hashlib.md5(b'-DMY_BINDING_PHRASE="myquad"').digest()[:6]
bits = np.r_[E.gfsk_frame_bits(uid, E.build_sync(5, 17, 9, 3, uid), "pn9_msb"), np.zeros(8, np.uint8)]
iq, _ = S.gfsk_burst(300e3, 100e3, bits=bits, bt=2.0, pre_bits=0, cfo_hz=6e3, fs=FS, rng=rng)
bursts.append(dict(t0_us=22000, f_hz=2440.4e6, iq=iq, bw_hz=0.5e6))
y, _ = S.embed(x, FS, FC, bursts, 15)
y = y.astype(np.complex64)
y.tofile(OUT / "capture.cf32")
em = [e for e in D.detect(y, FS, FC) if e.label.startswith(LABELS)]
fields = {}
for lab in LABELS:
    g = [e for e in em if e.label.startswith(lab)]
    if g:
        f = PLUGINS[0].summarise(g, None)["fields"]
        fields[lab] = {k: f.get(k) for k in ("radiolink", "afhds1", "hubsan", "toy_drone", "elrs_gfsk", "model")
                       if f.get(k) is not None}
ref = dict(fs=FS, fc=FC, labels=[e.label for e in em], decoded=[e.decoded for e in em], fields=fields,
           radiolink_hops=[[t.hex(), *RL.hop_table(t)] for t in (tx, bytes.fromhex("3a99223a"), bytes([1, 2, 3, 4]))],
           afhds1_hops=[[v, A.afhds1_hops(v)] for v in (0x78, 0x2B, 0x55, 0x00, 0x9F)])
print(ref["labels"])
(OUT / "cases.json").write_text(json.dumps(ref, indent=1, default=float))
