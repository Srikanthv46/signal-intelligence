#!/usr/bin/env python3
"""Synthetic reference for MAVLink (core/mavlink) and SiK telemetry (rc_link/sik); research, no capture.
  mavlink.json   a byte stream of v1 / v2 frames + junk, and core.mavlink.scan / summary on it;
  sik.cf32       one SiK packet (NETID 0x19, 64 kb/s, PN9 LSB-first, inverted polarity) at 512 kS/s, 15 dB SNR,
                 and rc_link.sik.sik_try_decode on it (sik.json).

    nice -n 19 .venv/bin/python -I c++/tests/vectors/sik/gen_sik.py
"""
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
REPO = OUT.parents[3]
sys.path.insert(0, str(REPO / "python"))
sys.path.insert(0, str(REPO / "python" / "tests"))
from test_sik_telemetry import frames                                     # noqa: E402
from spectrum_intel.core import mavlink as MAV                             # noqa: E402
from spectrum_intel.receivers.rc_link import sik as SIK, synth as S        # noqa: E402

stream = b"\x00\x11" + frames() + b"\x07\xfe\x03"
fr = MAV.scan(stream)
(OUT / "mavlink.json").write_text(json.dumps({"hex": stream.hex(), "scan": fr, "summary": MAV.summary(fr)}, indent=1))

fs, rate = 512e3, 64e3
rng = np.random.default_rng(21)
bits = S.sik_packet_bits(0x0019, frames() + b"\x12\x34", whitening="pn9_lsb_first") ^ 1
iq, _ = S.gfsk_burst(rate, 63.75e3, bits=bits, bt=0.5, cfo_hz=1.5e3, fs=fs, rng=rng)
pad = np.zeros(300, np.complex64)
y = np.r_[pad, iq, pad]
y = (y + (rng.standard_normal(len(y)) + 1j * rng.standard_normal(len(y))) * np.sqrt(fs / 200e3 * 10 ** (-15 / 10) / 2)).astype(np.complex64)
y.tofile(OUT / "sik.cf32")
d = SIK.sik_try_decode(y, fs, rate)
(OUT / "sik.json").write_text(json.dumps({"fs": fs, "rate": rate, "out": d}, indent=1, default=float))
print("netid", d.get("netid"), "mavlink ok", d.get("mavlink_crc_ok"), "inverted", d.get("polarity_inverted"), d.get("whitening"))
