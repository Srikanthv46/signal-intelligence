"""Drone links inside 802.11 / LoRa frames (research, PLAN G1 / G2 / G7; research/drone_transmitters/14): wfb-ng /
OpenIPC raw injection, ESP-NOW + ExpressLRS Backpack, MAVLink over UDP / TCP and transparent LoRa. Synthetic frames."""
import hashlib
import struct

from spectrum_intel.core import mavlink
from spectrum_intel.receivers.rc_link import lora_phy as P
from spectrum_intel.receivers.wifi import mac80211 as M


def _hb(sysid=1):
    pl = bytes([0, 0, 0, 0, 2, 3, 0x81, 4, 3])
    h = bytes([0xFE, len(pl), 7, sysid, 1, 0])
    c = mavlink.x25(bytes([50]), mavlink.x25(h[1:] + pl))
    return h + pl + bytes([c & 0xFF, c >> 8])


def test_wfb_ng():
    a = bytes([0x57, 0x42]) + ((0x7505D6 << 8) | 0x10).to_bytes(4, "big")
    i = M.parse(bytes([0x08, 0x01, 0, 0]) + b"\xff" * 6 + a + a + b"\x10\x00" + bytes([1]) +
                ((77 << 8) | 2).to_bytes(8, "big") + bytes(30))
    w = i["wfb_ng"]
    assert w["link_id"] == 0x7505D6 and w["default_link_domain"] and w["radio_port"] == 0x10
    assert w["stream"] == "MAVLink" and (w["block_idx"], w["fragment_idx"]) == (77, 2)
    assert i["drone"]["platform"].startswith("wfb-ng")
    a2 = bytes([0x57, 0x42, 0, 0, 0, 0x01])                     # not broadcast: not wfb-ng
    assert "wfb_ng" not in M.parse(bytes([0x08, 0x01, 0, 0]) + a2 + a2 + a2 + b"\x10\x00" + bytes(10))


def _espnow(uid, msp, dst=None):
    body = bytes([127]) + bytes.fromhex("18fe34") + b"\x01\x02\x03\x04" + bytes([221, 5 + len(msp)]) + \
        bytes.fromhex("18fe34") + bytes([4, 1]) + msp
    return bytes([0xD0, 0, 0, 0]) + (dst or uid) + uid + b"\xff" * 6 + b"\x20\x00" + body


def _msp(func, pay):
    fr = b"\x00" + func.to_bytes(2, "little") + len(pay).to_bytes(2, "little") + pay
    return b"$X<" + fr + bytes([M.crc8_dvb_s2(fr)])


def test_espnow_backpack():
    uid = bytearray(hashlib.md5(b'-DMY_BINDING_PHRASE="myquad"').digest()[:6])
    uid[0] &= 0xFE
    uid = bytes(uid)
    i = M.parse(_espnow(uid, _msp(0x59, (4 * 8 + 0).to_bytes(2, "little") + b"\0\0")))        # R1
    assert i["esp_now"]["msp"]["crc_ok"] and i["elrs_backpack"]["bind_uid"] == M._mac(uid)
    assert (i["elrs_backpack"]["vtx_band"], i["elrs_backpack"]["vtx_channel"]) == ("R", 1)
    assert i["drone"]["platform"].startswith("ExpressLRS Backpack")
    bad = bytearray(_msp(0x59, b"\x20\x00\x00\x00"))
    bad[-1] ^= 1
    i = M.parse(_espnow(uid, bytes(bad)))
    assert "elrs_backpack" not in i and i["esp_now"]["msp"]["crc_ok"] is False
    i = M.parse(_espnow(uid, _msp(0x09, uid), dst=b"\xff" * 6))                           # bind to broadcast
    assert i["elrs_backpack"]["msp"] == "ELRS bind"


def test_mavlink_udp_tcp_lora():
    u = struct.pack(">HHHH", 14550, 14550, 8 + len(_hb()), 0) + _hb()
    ip = bytes([0x45, 0, 0, 20 + len(u), 0, 0, 0, 0, 64, 17, 0, 0, 192, 168, 4, 1, 192, 168, 4, 2]) + u
    i = M.parse(bytes([0x08, 0x02, 0, 0]) + b"\xff" * 6 + bytes(6) * 2 + b"\x10\x00" +
                b"\xaa\xaa\x03\x00\x00\x00\x08\x00" + ip)
    assert i["ipv4"]["mavlink"]["vehicle_type"] == "quadrotor" and i["drone"]["confidence"] == 0.9
    assert P.mavlink_in(b"\x00" + _hb(3))["sysid"] == 3 and P.mavlink_in(bytes(20)) is None
    i = M.parse(bytes([0x08, 0x02, 0, 0]) + b"\xff" * 6 + bytes(6) * 2 + b"\x10\x00" +
                b"\xaa\xaa\x03\x00\x00\x00\x08\x00" + ip[:-1] + bytes([ip[-1] ^ 1]))
    assert "mavlink" not in i["ipv4"] and "drone" not in i
