#!/usr/bin/env python3
"""Synthetic reference for Remote ID over Bluetooth LE (core/remote_id + rc_link ble_try_decode; research:
research/drone_transmitters/04_telemetry_remoteid.md). No capture is used:
  ad.json        advertising-data bytes (flags + F3411 service data + name) and core.remote_id.parse_ble_ad output;
  ble_rid.cf32   one ADV_NONCONN_IND on channel 37 at 4 MSps (complex baseband, 25 dB SNR in 1 MHz), and
                 rc_link ble_try_decode's output on it (ble.json).

    nice -n 19 .venv/bin/python -I c++/tests/vectors/remote_id/gen_remote_id.py
"""
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
REPO = OUT.parents[3]
sys.path.insert(0, str(REPO / "python"))
sys.path.insert(0, str(REPO / "python" / "tests" / "receivers"))
from test_det_droneid import odid_ie_body                                  # noqa: E402
from spectrum_intel.core import remote_id as RID                           # noqa: E402
from spectrum_intel.receivers.rc_link import detector as D, synth as S     # noqa: E402

ad = bytes([2, 0x01, 0x06]) + S.remote_id_adv_data(odid_ie_body()) + bytes([5, 0x09]) + b"DRON" + \
    bytes([5, 0xFF, 0x22, 0x08, 0x01, 0x02])
(OUT / "ad.json").write_text(json.dumps({"adv_data_hex": ad.hex(), "out": RID.parse_ble_ad(ad)}, indent=1))

fs = 4e6
rng = np.random.default_rng(9)
iq, _ = S.ble_adv_burst(bytes([0x11, 0x22, 0x33, 0x44, 0x55, 0xC6]), S.remote_id_adv_data(odid_ie_body()), ch=37,
                        cfo_hz=15e3, fs=fs, rng=rng)
pad = np.zeros(200, np.complex64)
y = np.r_[pad, iq, pad]
noise = (rng.standard_normal(len(y)) + 1j * rng.standard_normal(len(y))) * np.sqrt(fs / 1e6 * 10 ** (-25 / 10) / 2)
y = (y + noise).astype(np.complex64)
y.tofile(OUT / "ble_rid.cf32")
d = D.ble_try_decode(y.astype(np.complex128), fs, 2402e6)
(OUT / "ble.json").write_text(json.dumps({"fs": fs, "f_hz": 2402e6, "out": d}, indent=1, default=float))
print("crc_ok", d.get("crc_ok"), "uas_id", d.get("ad", {}).get("remote_id", {}).get("uas_id"))
