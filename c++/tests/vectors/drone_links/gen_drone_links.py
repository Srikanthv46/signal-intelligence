#!/usr/bin/env python3
"""Reference for tests/unit/test_drone_links.cpp (PLAN G1 / G2 / G7): 802.11 MPDUs and LoRa payloads carrying
MAVLink (UDP 14550 / TCP 5760 over Wi-Fi, transparent LoRa), parsed by python receivers.wifi.mac80211.parse and
receivers.rc_link.lora_phy.mavlink_in. Synthetic only.

    nice -n 19 .venv/bin/python -I c++/tests/vectors/drone_links/gen_drone_links.py
"""
import json
import struct
import sys
from pathlib import Path

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parents[3] / "python"))
from spectrum_intel.core import mavlink                                   # noqa: E402
from spectrum_intel.receivers.rc_link import lora_phy as P                # noqa: E402
from spectrum_intel.receivers.wifi import mac80211 as M                   # noqa: E402


def mav_v1(msgid, extra, sysid, payload, seq=7):
    hdr = bytes([0xFE, len(payload), seq, sysid, 1, msgid])
    c = mavlink.x25(bytes([extra]), mavlink.x25(hdr[1:] + payload))
    return hdr + payload + bytes([c & 0xFF, c >> 8])


HB = mav_v1(0, 50, 1, bytes([0, 0, 0, 0, 2, 3, 0x81, 4, 3]))             # quad, ArduPilot, armed, STABILIZE
GPI = mav_v1(33, 104, 1, struct.pack("<IiiiihhhH", 1000, 283874406, 768928889, 250000, 12000, 0, 0, 0, 9000))
GCS_HB = mav_v1(0, 50, 255, bytes([0, 0, 0, 0, 6, 8, 0, 0, 3]))


def ip_udp(payload, sport=14550, dport=14550):
    u = struct.pack(">HHHH", sport, dport, 8 + len(payload), 0) + payload
    return bytes([0x45, 0, 0, 20 + len(u), 0, 0, 0, 0, 64, 17, 0, 0, 192, 168, 4, 1, 192, 168, 4, 2]) + u


def ip_tcp(payload, sport=5760, dport=40000):
    t = struct.pack(">HHIIBBHHH", sport, dport, 1, 0, 5 << 4, 0x18, 1024, 0, 0) + payload
    return bytes([0x45, 0, 0, 20 + len(t), 0, 0, 0, 0, 64, 6, 0, 0, 192, 168, 4, 1, 192, 168, 4, 2]) + t


def data_frame(ip, fc1=0x02):
    return (bytes([0x08, fc1, 0, 0]) + bytes.fromhex("ffffffffffff") + bytes.fromhex("5ccf7f010203") * 2 + b"\x10\x00"
            + b"\xaa\xaa\x03\x00\x00\x00\x08\x00" + ip)


frames = {"udp_hb_gpi": data_frame(ip_udp(HB + GPI)), "tcp_hb": data_frame(ip_tcp(b"\x00\x01" + HB)),
          "udp_gcs_only": data_frame(ip_udp(GCS_HB)), "udp_noise": data_frame(ip_udp(bytes(range(40)))),
          "udp_bad_crc": data_frame(ip_udp(HB[:-1] + bytes([HB[-1] ^ 1])))}
import hashlib                                                             # noqa: E402
_a = bytes([0x57, 0x42]) + ((0x7505D6 << 8) | 0x10).to_bytes(4, "big")
frames["wfb_data"] = bytes([0x08, 0x01, 0, 0]) + b"\xff" * 6 + _a + _a + b"\x10\x00" + bytes([1]) + \
    ((77 << 8) | 2).to_bytes(8, "big") + bytes(30)
_b = bytes([0x57, 0x42, 0x12, 0x34, 0x56, 0xA0])
frames["wfb_session"] = bytes([0x08, 0x01, 0, 0]) + b"\xff" * 6 + _b + _b + b"\x20\x00" + bytes([2]) + bytes(60)
_uid = bytearray(hashlib.md5(b'-DMY_BINDING_PHRASE="myquad"').digest()[:6])
_uid[0] &= 0xFE
_uid = bytes(_uid)


def _msp(func, pay):
    fr = b"\x00" + func.to_bytes(2, "little") + len(pay).to_bytes(2, "little") + pay
    return b"$X<" + fr + bytes([M.crc8_dvb_s2(fr)])


def _espnow(msp, dst=None):
    body = bytes([127]) + bytes.fromhex("18fe34") + b"\x01\x02\x03\x04" + bytes([221, 5 + len(msp)]) + \
        bytes.fromhex("18fe34") + bytes([4, 1]) + msp
    return bytes([0xD0, 0, 0, 0]) + (dst or _uid) + _uid + b"\xff" * 6 + b"\x20\x00" + body


frames["espnow_vtx"] = _espnow(_msp(0x59, (4 * 8).to_bytes(2, "little") + b"\0\0"))
frames["espnow_freq"] = _espnow(_msp(0x59, (5800).to_bytes(2, "little") + b"\0\0"))
frames["espnow_chidx"] = _espnow(_msp(0x0301, bytes([10])))
frames["espnow_bind"] = _espnow(_msp(0x09, _uid), dst=b"\xff" * 6)
_bad = bytearray(_msp(0x59, b"\x20\x00\x00\x00"))
_bad[-1] ^= 1
frames["espnow_badcrc"] = _espnow(bytes(_bad))
ref = {"mpdu": {k: dict(hex=v.hex(), out=M.parse(v)) for k, v in frames.items()},
       "lora": {k: dict(hex=v.hex(), out=P.mavlink_in(v)) for k, v in
                {"hb": b"\x11" + HB + b"\x22", "hb_gpi": HB + GPI, "split": HB[:10], "none": bytes(32)}.items()}}
(OUT / "ref.json").write_text(json.dumps(ref, indent=1))
for k, v in ref["mpdu"].items():
    print(k, v["out"].get("drone"), (v["out"].get("ipv4") or {}).get("mavlink"))
for k, v in ref["lora"].items():
    print("lora", k, v["out"])
