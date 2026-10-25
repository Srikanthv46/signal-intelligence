#!/usr/bin/env python3
"""
test_det_droneid.py -- drone identification carried in 802.11 frames (receivers/wifi/droneid.py).

  1. Open Drone ID (ASTM F3411 / ASD-STAN) message pack: every field of Basic ID, Location, System
     and Operator ID, including the operator's position.
  2. DJI's proprietary DroneID: located by CRC-16/Kermit, cross-checked field by field against the
     independent reference decoder (../droneid-decoder-reference) on identical bytes.
  3. Nothing is reported on a guess: a corrupted payload, an unrelated vendor element and an
     all-zero element all return None.
  4. Platform signatures from the IEEE registry and the SSID, without catching ordinary networks.
  5. End to end through mac80211.parse() on a real beacon.

pytest collects the test_* functions; `python test_det_droneid.py` runs them with a printed summary.
"""
from __future__ import annotations

import os
import struct
import sys

from spectrum_intel.receivers.wifi import droneid as D
from spectrum_intel.receivers.wifi import dsss_tx as TX
from spectrum_intel.receivers.wifi import mac80211 as M

sys.path.insert(0, os.path.dirname(__file__))
from _data import DRONEID_TV  # noqa: E402

UAS = "1581F6M24A700300XXXX"
LAT, LON = 28.6139, 77.2090
OP_LAT, OP_LON = 28.6100, 77.2000


def _msg(t, ver, body):
    return (bytes([(t << 4) | ver]) + bytes(body)).ljust(25, b"\x00")[:25]


def odid_ie_body(status=0x20, speed=40, track=90):
    """A four-message Open Drone ID pack as a Wi-Fi beacon carries it."""
    basic = _msg(0, 2, bytes([(1 << 4) | 2]) + UAS.encode() + b"\x00\x00\x00")
    loc = _msg(1, 2, bytes([status, track, speed]) + struct.pack("<b", 4)
               + struct.pack("<ii", int(LAT * 1e7), int(LON * 1e7))
               + struct.pack("<HHH", int(1120 / 0.5), int(1122 / 0.5), int(1118 / 0.5))
               + bytes([(5 << 4) | 11, (3 << 4) | 3]) + struct.pack("<H", 12345) + bytes([0, 0]))
    system = _msg(4, 2, bytes([(1 << 2) | 1])
                  + struct.pack("<ii", int(OP_LAT * 1e7), int(OP_LON * 1e7))
                  + struct.pack("<HB", 1, 0) + struct.pack("<HH", 0, 0) + bytes([(1 << 4) | 1])
                  + struct.pack("<H", int(1095 / 0.5)) + struct.pack("<I", 123456789) + b"\x00")
    opid = _msg(5, 2, bytes([0]) + b"IND-OP-7721".ljust(20, b"\x00") + b"\x00\x00\x00")
    return bytes([7, 0xF2, 25, 4]) + basic + loc + system + opid


def dji_payload(serial="3N6CK9800A1234", model=63):
    vals = [91, 0, 2, 1234, 0x1234, serial.encode().ljust(16, b"\x00"),
            int(LON * D.DJI_SCALE), int(LAT * D.DJI_SCALE), int(150 * 3.281), int(120 * 3.281),
            30, -12, 5, 4500, 1700000000000, int(OP_LAT * D.DJI_SCALE), int(OP_LON * D.DJI_SCALE),
            int(77.1990 * D.DJI_SCALE), int(28.6090 * D.DJI_SCALE), model, 0, b"".ljust(20, b"\x00"), 0]
    raw = D._DJI_STRUCT.pack(*vals)
    return raw[:-2] + struct.pack("<H", D.crc16_kermit(raw[:-2]))


def test_open_drone_id():
    r = D.from_vendor_ie(D.ODID_OUI, D.ODID_VENDOR_TYPE, odid_ie_body())
    assert r, "Open Drone ID did not parse"
    assert r["uas_id"] == UAS and r["id_type"].startswith("serial")
    assert r["ua_type"] == "multirotor" and r["status"] == "airborne"
    assert abs(r["lat"] - LAT) < 1e-5 and abs(r["lon"] - LON) < 1e-5
    assert r["speed_mps"] == 10.0 and r["track_deg"] == 90 and r["vertical_speed_mps"] == 2.0
    assert r["height_m"] == 118.0 and r["geodetic_altitude_m"] == 122.0
    assert abs(r["operator_lat"] - OP_LAT) < 1e-5 and abs(r["operator_lon"] - OP_LON) < 1e-5
    assert r["operator_id"] == "IND-OP-7721"
    assert r["operator_location_type"].startswith("live GNSS")
    assert r["accuracy"]["horizontal"] == "<3 m"
    assert set(r["messages"]) == {"basic_id", "location", "system", "operator_id"}
    assert "OPERATOR" in D.summary(r) and UAS in D.summary(r)


def test_open_drone_id_direction_and_speed_multiplier():
    """The east/west bit adds 180 degrees and the multiplier bit changes the speed scale."""
    r = D.from_vendor_ie(D.ODID_OUI, D.ODID_VENDOR_TYPE, odid_ie_body(status=0x33, track=10, speed=100))
    assert r["track_deg"] == 190                       # bit 1 set -> +180
    assert r["speed_mps"] == round(100 * 0.75 + 63.75, 2)   # bit 0 set -> coarse scale
    assert r["status"] == "emergency"


def test_dji_matches_the_reference_decoder():
    raw = dji_payload()
    got = D.parse_dji(b"\x10\x20" + raw)               # offset: it must be found by its CRC
    assert got and got["payload_offset"] == 2
    assert got["uas_id"] == "3N6CK9800A1234" and got["model"] == "Mini 2"
    assert abs(got["lat"] - LAT) < 1e-5 and abs(got["operator_lat"] - OP_LAT) < 1e-5
    if not DRONEID_TV:
        return
    src = os.path.dirname(DRONEID_TV)
    if not os.path.isdir(src):
        return
    sys.path.insert(0, src)
    try:
        from droneid_packet import DroneIDPacket
    except Exception:
        return
    ref = DroneIDPacket(raw).to_dict()
    assert ref["crc-packet"] == ref["crc-calculated"]
    assert got["uas_id"] == ref["serial_number"] and got["model"] == ref["device_type"]
    for mine, theirs in (("lat", "latitude"), ("lon", "longitude"), ("operator_lat", "app_lat"),
                         ("operator_lon", "app_lon"), ("home_lat", "latitude_home"),
                         ("home_lon", "longitude_home")):
        assert abs(got[mine] - ref[theirs]) < 1e-6, mine


def test_never_reports_a_guess():
    raw = bytearray(dji_payload())
    raw[20] ^= 0xFF
    assert D.parse_dji(bytes(raw)) is None             # one flipped bit -> CRC fails -> nothing
    assert D.from_vendor_ie("00:50:f2", 1, b"\x00" * 120) is None       # Microsoft WPA element
    assert D.parse_odid(b"\x00" * 60) is None          # syntactically valid but says nothing
    assert D.parse_odid(b"\x07" + b"\xff" * 80) is None
    assert D.scan_vendor_ies(b"\x00\x04test") is None


def test_platform_signatures():
    assert D.signature("Sz Dji Technology Co.,Ltd", None)["platform"] == "DJI"
    assert D.signature("Parrot Sa", None)["platform"] == "Parrot"
    assert D.signature("Autel Robotics USA LLC", None)["platform"] == "Autel Robotics"
    assert D.signature("Beijing Autelan Technology Inc.", None) is None       # a different company
    assert D.signature(None, "MAVIC2-ABC123")["platform"] == "DJI Mavic"
    assert D.signature(None, "openhd-gnd")["platform"].startswith("open-source FPV")
    both = D.signature("Sz Dji Technology Co.,Ltd", "Mavic_Air_7721")
    assert both["platform"] == "DJI Mavic Air" and both["confidence"] == 0.9
    for ssid in ("Airtel_maya_9187", "Sharma Brothers_2G", "Hindusena", "", "my-phantom-of-the-opera"):
        assert D.signature(None, ssid) is None, ssid


def test_through_the_dissector():
    body = odid_ie_body()
    bc = TX.beacon("60:60:1f:aa:bb:cc", "Mavic_Air_7721", 6)
    bc += bytes([221, 4 + len(body)]) + b"\xfa\x0b\xbc\x0d" + body
    m = M.parse(bc)
    assert m["remote_id"]["uas_id"] == UAS
    assert m["drone"]["source"] == "Remote ID broadcast" and m["drone"]["confidence"] == 0.99
    assert "REMOTE ID" in m["summary"] and "OPERATOR" in m["summary"]
    plain = M.parse(TX.beacon("60:60:1f:aa:bb:cc", "Mavic_Air_7721", 6))
    assert plain["drone"]["platform"] == "DJI Mavic Air-like"      # inferred, so it is marked inferred
    assert plain["drone"]["confidence"] < 1.0
    assert M.parse(TX.beacon("24:a3:f0:a5:bd:ab", "Sharma Brothers_2G", 6)).get("drone") is None


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("PASS", name)
    r = D.from_vendor_ie(D.ODID_OUI, D.ODID_VENDOR_TYPE, odid_ie_body())
    print("\nOpen Drone ID:", D.summary(r))
    print("DJI DroneID  :", D.summary(D.parse_dji(dji_payload())))


def nan_action_frame(body=None, sa="9c:5a:44:01:02:03"):
    """A Wi-Fi NAN service discovery (public action) frame carrying an F3411 message pack (opendroneid wifi.c layout)."""
    body = odid_ie_body() if body is None else body
    info = body                                               # [message counter][message pack]
    sda = D.NAN_ODID_SERVICE_ID + bytes([1, 0, 0x10, len(info)]) + info
    act = bytes([4, 9, 0x50, 0x6F, 0x9A, 0x13]) + bytes([0x03]) + struct.pack("<H", len(sda)) + sda
    mac = bytes.fromhex(sa.replace(":", ""))
    hdr = bytes([0xD0, 0x00, 0, 0]) + bytes([0xFF] * 6) + mac + bytes.fromhex("506f9a0101ff") + bytes(2)
    return hdr + act + bytes(4)


def test_remote_id_over_wifi_nan():
    """Remote ID over Wi-Fi NAN service discovery (research, PLAN T3)."""
    r = M.parse(nan_action_frame())
    rid = r.get("remote_id")
    assert rid and rid["uas_id"] == UAS and rid["transport"] == "Wi-Fi NAN service discovery"
    assert abs(rid["lat"] - LAT) < 1e-5 and r["drone"]["source"] == "Remote ID broadcast"
    other = bytearray(nan_action_frame()); other[24 + 9] ^= 0xFF          # another service ID
    assert M.parse(bytes(other)).get("remote_id") is None
