"""Remote ID over Bluetooth 5 Long Range (LE Coded PHY; research, PLAN T2b; bench/ble_coded): ADV_EXT_IND on a
primary channel pointing to an AUX_ADV_IND that carries the F3411 message pack, through rc_link detect."""
import sys
import time
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent / "receivers"))
from test_det_droneid import LAT, LON, UAS, odid_ie_body                    # noqa: E402

from spectrum_intel.receivers.rc_link import ble_coded as BC, detector as D, synth as S   # noqa: E402

FC, FS = 2420e6, 61.44e6
ADDR = bytes([0x21, 0x43, 0x65, 0x87, 0xA9, 0xCB])


def test_viterbi_roundtrip():
    rng = np.random.default_rng(0)
    b = list(rng.integers(0, 2, 200)) + [0, 0, 0]
    c = 2.0 * np.array(BC.conv_encode(b)) - 1
    c = c + 0.8 * rng.standard_normal(len(c))
    assert np.array_equal(BC.viterbi(c.reshape(-1, 2)), np.array(b, np.uint8))


def test_parse_extended():
    p = S.ext_adv_payload(ADDR, b"", did=0x456, sid=3, aux=(8, 166, 2))
    e = BC.parse_extended(p)
    assert e["adv_addr"] == "CB:A9:87:65:43:21" and e["adi"] == {"did": 0x456, "sid": 3}
    assert e["aux_ptr"] == {"channel": 8, "ca": 0, "offset_us": 166 * 30, "phy": "Coded"}


def capture(snr_db, seed=4, ci=0):
    rng = np.random.default_rng(seed)
    x = S.noise_capture(1 << 21, rng)
    body = bytes([0]) + odid_ie_body()[1:]                       # counter + message pack (F3411 over BT5)
    aux = S.ext_adv_payload(ADDR, S.remote_id_adv_data(body), did=0x456, sid=3)
    prim = S.ext_adv_payload(ADDR, b"", did=0x456, sid=3, aux=(8, 166, 2))
    ip, _ = S.ble_coded_burst(prim, ch=37, cfo_hz=12e3, rng=rng)
    ia, _ = S.ble_coded_burst(aux, ch=8, ci=ci, cfo_hz=12e3, rng=rng)
    bursts = [dict(t0_us=2000, f_hz=2402e6, iq=ip, bw_hz=1.0e6), dict(t0_us=8000, f_hz=2420e6, iq=ia, bw_hz=1.0e6)]
    return S.embed(x, FS, FC, bursts, snr_db)[0]


@pytest.mark.parametrize("snr_db,ci", [(12, 0), (1, 0), (8, 1)])
def test_remote_id_long_range(snr_db, ci):
    t = time.perf_counter()
    em = D.detect(capture(snr_db, ci=ci), FS, FC)
    dt = time.perf_counter() - t
    dec = [e.decoded for e in em if (e.decoded or {}).get("phy", "").startswith("LE Coded")]
    prim = [d for d in dec if d["channel_index_rf"] == 37]
    aux = [d for d in dec if d["channel_index_rf"] == 8]
    assert aux, f"no AUX_ADV_IND at {snr_db} dB ({dt:.1f} s)"
    a = aux[0]
    assert a["crc_ok"] and a["adv_addr"] == "CB:A9:87:65:43:21" and a["phy"] == ("LE Coded S=8" if ci == 0 else "LE Coded S=2")
    r = a["ad"]["remote_id"]
    assert r["uas_id"] == UAS and abs(r["lat"] - LAT) < 1e-5 and abs(r["lon"] - LON) < 1e-5
    if ci == 0:
        assert prim and prim[0]["extended"]["aux_ptr"]["channel"] == 8
