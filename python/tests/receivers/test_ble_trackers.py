"""BLE item trackers and notable devices (research, signal gaps S7) from advertising data."""
import struct

import numpy as np

from spectrum_intel.core.ble_trackers import recognize
from spectrum_intel.core.remote_id import parse_ble_ad


def mfg(cid, payload):
    v = struct.pack("<H", cid) + payload
    return bytes([len(v) + 1, 0xFF]) + v


def svc(uuid, payload):
    v = struct.pack("<H", uuid) + payload
    return bytes([len(v) + 1, 0x16]) + v


AIRTAG = mfg(0x004C, bytes([0x12, 0x19, 0x10 | 0x04 | 0xC0]) + bytes(22) + bytes([0x00, 0x00]))
SMARTTAG = svc(0xFD5A, bytes.fromhex("156bfa00c84062b28f00e260c3000000ad018b47"))   # Yu et al. example, state 0


def test_apple_find_my_and_ibeacon_and_airpods():
    d = parse_ble_ad(AIRTAG)["device"]
    assert d["label"] == "Apple Find My / AirTag" and d["tracker"] and d["battery"] == "full" and d["offline_finding"]
    ib = parse_ble_ad(mfg(0x004C, bytes([0x02, 0x15]) + bytes(range(16)) + bytes([0, 1, 0, 2, 0xC5])))["device"]
    assert ib["kind"] == "iBeacon" and ib["major"] == 1 and ib["minor"] == 2 and ib["measured_power_dbm"] == -59
    assert not ib["tracker"]
    ap = parse_ble_ad(mfg(0x004C, bytes([0x07, 0x05, 0x01, 0x0E, 0x20, 0x55, 0x66])))["device"]
    assert ap["model"] == "AirPods Pro" and not ap["tracker"]


def test_samsung_google_tile_dult():
    # the published example has tag state 0 (outside 1-6): not claimed
    assert parse_ble_ad(SMARTTAG).get("device") is None
    st = svc(0xFD5A, bytes([0x15 | (2 << 5)]) + bytes(11) + bytes([0xE0]) + bytes(7))
    d = parse_ble_ad(st)["device"]
    assert d["label"] == "Samsung SmartTag" and d["state"] == "lost" and d["uwb"] and d["battery"] == "full"
    g = parse_ble_ad(svc(0xFEAA, bytes([0x41]) + bytes(20) + bytes([0x20])))["device"]
    assert g["kind"] == "Find My Device (FMDN)" and g["unwanted_tracking_protection"] and g["tracker"]
    e = parse_ble_ad(svc(0xFEAA, bytes([0x10, 0xEB, 0x03]) + b"example"))["device"]
    assert e["label"] == "Eddystone-URL" and not e["tracker"]
    assert parse_ble_ad(svc(0xFEED, b"\x02\x00\x11\x22"))["device"]["label"] == "Tile-like"
    du = parse_ble_ad(svc(0xFCB2, bytes([0x05, 0x01, 0xAA])))["device"]
    assert du["near_owner"] and du["tracker"]
    assert parse_ble_ad(svc(0xFCB2, bytes([0x05, 0x03]))).get("device") is None      # reserved bits set


def test_negatives_and_swift_pair():
    assert parse_ble_ad(bytes([2, 0x01, 0x06])).get("device") is None
    assert recognize({"company_id": "0x004C", "manufacturer_data_hex": "12"}) is None     # truncated TLV
    assert parse_ble_ad(mfg(0x0006, bytes([0x03, 0x05, 0x80]))).get("device") is None    # bad sub-scenario
    sp = parse_ble_ad(mfg(0x0006, bytes([0x03, 0x00, 0x80]) + b"Mouse"))["device"]
    assert sp["kind"] == "Swift Pair" and sp["display_name"] == "Mouse"


def test_airtag_through_rc_link():
    from spectrum_intel.receivers.rc_link import detector as D
    from spectrum_intel.receivers.rc_link import synth as S
    from spectrum_intel.receivers.rc_link.plugin import PLUGINS
    FC, FS = 2420e6, 61.44e6
    rng = np.random.default_rng(3)
    n = 1 << 21
    x = ((rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)).astype(np.complex64)
    addr = bytes([0x11, 0x22, 0x33, 0x44, 0x55, 0xC6])
    bursts = []
    for k in range(3):
        iq, _ = S.ble_adv_burst(addr, AIRTAG, ch=37, cfo_hz=8e3, rng=rng)
        bursts.append(dict(t0_us=2000 + 10000 * k, f_hz=2402e6, iq=iq, bw_hz=1.0e6))
    y, _ = S.embed(x, FS, FC, bursts, 20)
    em = [e for e in D.detect(y, FS, FC) if (e.decoded or {}).get("crc_ok")]
    assert em and em[0].decoded["ad"]["device"]["label"] == "Apple Find My / AirTag"
    rc = next(p for p in PLUGINS if p.name == "rc_link")
    assert rc.pretty("ble", em) == "Bluetooth LE tracker (Apple Find My / AirTag)"
    out = rc.summarise(em, None)
    assert out["fields"]["ble_devices"][0]["label"] == "Apple Find My / AirTag"
    assert "item tracker" in out["evidence"]
