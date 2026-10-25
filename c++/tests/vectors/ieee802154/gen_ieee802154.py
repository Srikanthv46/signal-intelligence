#!/usr/bin/env python3
"""Synthetic reference for IEEE 802.15.4 O-QPSK (rc_link/ieee802154; research, PLAN T4): frames at 8 MHz
(4 samples / chip) with noise, demod output; CRC and MAC-parse cases.

    nice -n 19 .venv/bin/python -I c++/tests/vectors/ieee802154/gen_ieee802154.py
"""
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parents[3] / "python"))
from spectrum_intel.receivers.rc_link import ieee802154 as Z   # noqa: E402

rng = np.random.default_rng(81)


def mpdu(fcf, body):
    m = bytes([fcf & 0xFF, fcf >> 8]) + bytes(body)
    c = Z.crc16_kermit(m)
    return m + bytes([c & 0xFF, c >> 8])


frames = [
    mpdu(1 | 1 << 5 | 1 << 6 | 2 << 10 | 3 << 14, bytes([0x42, 0x34, 0x12, 0xFF, 0xFF]) +
         bytes.fromhex("0013a20040a1b2c3")[::-1] + bytes(rng.integers(0, 256, 20).tolist())),   # XBee data
    mpdu(0 | 2 << 14, bytes([7, 0xCD, 0xAB, 0x01, 0x00]) + bytes(10)),                          # beacon
    mpdu(3 | 2 << 10 | 2 << 14 | 1 << 6, bytes([9, 0x34, 0x12, 0x00, 0x00, 0x01, 0x00, 0x04])),  # MAC command
]
iqs, outs = [], []
for k, (m, snr) in enumerate(zip(frames + frames[:1], (20, 8, 4, 30))):
    iq = Z.modulate(Z.chips_for(m))
    iq = np.r_[np.zeros(100), iq, np.zeros(100)] * np.exp(1j * rng.uniform(0, 6.28))
    w = iq + (rng.normal(size=len(iq)) + 1j * rng.normal(size=len(iq))) * np.sqrt(10 ** (-snr / 10) / 2)
    if k == 3:
        w = w[:len(w) // 2]                                                  # truncated: no frame
    w = w.astype(np.complex64)
    iqs.append(w)
    outs.append(Z.demod(w.astype(np.complex128), 8e6))
np.concatenate(iqs).tofile(OUT / "frames.cf32")
ref = dict(lengths=[len(w) for w in iqs], out=outs, channels=[[f, Z.channel(f)] for f in (2405e6, 2480e6, 2440.2e6, 2441e6, 2400e6)],
           crc=[[m.hex(), Z.crc16_kermit(m)] for m in frames], mac=[[m[:-2].hex(), Z.parse_mac(m[:-2])] for m in frames])
# whole rc_link pipeline: 7.68 MSps, fc 2440 MHz, channel 18: two XBee frames and one with a corrupted FCS
from spectrum_intel.receivers.rc_link import detector as D, synth as S   # noqa: E402
from spectrum_intel.receivers.rc_link.plugin import PLUGINS            # noqa: E402
FS, FC = 7.68e6, 2440e6
x = S.noise_capture(1 << 18, rng)
bad = bytearray(frames[0]); bad[-1] ^= 0x55
bursts = []
for k, m in enumerate((frames[0], frames[1], bytes(bad))):
    iq, _ = S.oqpsk_burst(m, cfo_hz=12e3, fs=FS, rng=rng)
    bursts.append(dict(t0_us=1500 + 9000 * k, f_hz=2440e6, iq=iq, bw_hz=2.6e6))
y, _ = S.embed(x, FS, FC, bursts, 14)
y = y.astype(np.complex64)
y.tofile(OUT / "capture.cf32")
em = [e for e in D.detect(y, FS, FC) if e.label.startswith("ieee802154")]
ref["fs"], ref["fc"] = FS, FC
ref["labels"], ref["decoded"] = [e.label for e in em], [e.decoded for e in em]
fl = PLUGINS[0].summarise(em, None)["fields"]
ref["field"], ref["model"] = fl.get("ieee802154"), fl.get("model")
print([o and (o["crc_ok"], o.get("frame_type")) for o in outs], ref["labels"], ref["field"])
(OUT / "cases.json").write_text(json.dumps(ref, indent=1, default=float))
