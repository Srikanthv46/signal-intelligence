"""Remote ID (ASTM F3411) over Bluetooth LE legacy advertising (research: research/drone_transmitters/04_telemetry_remoteid.md):
synthetic ADV_NONCONN_IND on channel 37 carrying a message pack (Basic ID, Location, System, Operator ID), through
rc_link detect + the plugin summary, from high to low SNR."""
import sys
import time
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent / "receivers"))
from test_det_droneid import LAT, LON, OP_LAT, OP_LON, UAS, odid_ie_body   # noqa: E402

from spectrum_intel.core import remote_id as RID                          # noqa: E402
from spectrum_intel.receivers.rc_link import detector as D, synth as S    # noqa: E402

FC, FS = 2420e6, 61.44e6


def test_ad_parser():
    ad = RID.parse_ble_ad(bytes([2, 0x01, 0x06]) + S.remote_id_adv_data(odid_ie_body()) + bytes([5, 0x09]) + b"DRON")
    assert ad["flags"] == 6 and ad["local_name"] == "DRON"
    r = ad["remote_id"]
    assert r["transport"].startswith("Bluetooth") and r["uas_id"] == UAS
    assert abs(r["lat"] - LAT) < 1e-5 and abs(r["operator_lat"] - OP_LAT) < 1e-5


def capture(snr_db, seed=3, n_adv=3):
    rng = np.random.default_rng(seed)
    n = 1 << 21
    x = ((rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)).astype(np.complex64)
    addr = bytes([0x11, 0x22, 0x33, 0x44, 0x55, 0xC6])
    bursts = []
    for k in range(n_adv):                                   # the same advert every ~10 ms (3 per block)
        iq, info = S.ble_adv_burst(addr, S.remote_id_adv_data(odid_ie_body()), ch=37, cfo_hz=12e3, rng=rng)
        bursts.append(dict(t0_us=2000 + 10000 * k, f_hz=2402e6, iq=iq, bw_hz=1.0e6))
    y, _ = S.embed(x, FS, FC, bursts, snr_db)
    return y


@pytest.mark.parametrize("snr_db", [20, 12, 8])   # 8 dB: coherent detector + Chase (PLAN L1; was 12 dB)
def test_remote_id_ble(snr_db):
    y = capture(snr_db)
    t = time.perf_counter()
    em = D.detect(y, FS, FC)
    dt = time.perf_counter() - t
    rid = [e.decoded["ad"]["remote_id"] for e in em if (e.decoded or {}).get("ad", {}).get("remote_id")]
    assert rid, f"no Remote ID at {snr_db} dB ({len(em)} emissions, {dt:.2f} s)"
    r = rid[0]
    assert r["uas_id"] == UAS and r["status"] == "airborne"
    assert abs(r["lat"] - LAT) < 1e-5 and abs(r["lon"] - LON) < 1e-5
    assert abs(r["operator_lat"] - OP_LAT) < 1e-5 and abs(r["operator_lon"] - OP_LON) < 1e-5
