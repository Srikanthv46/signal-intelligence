#!/usr/bin/env python3
"""Reference vectors for tests/unit/test_rx_wifi*.cpp, written by the Python Wi-Fi receiver.

    nice -n 19 .venv/bin/python c++/tests/vectors/wifi/gen_wifi.py

  mac_frames.json   real FCS-valid MPDUs that dsss.scan decodes on the two 2.4 GHz testdata captures
                    (Phantom 4 and D_20251119-183758) plus synthetic frames that reach every branch of
                    mac80211.parse (data/LLC/EAPOL/ARP/IP, control, auth, assoc, action, WPS, P2P,
                    Open Drone ID and DJI DroneID vendor elements), each with mac80211.parse(bytes)
  droneid.json      parse_odid / parse_dji / signature / summary cases, and every IEEE vendor name
                    whose droneid.signature() matches (the C++ must match exactly the same names)
  links.json        links.attribute() on synthetic PPDU records (input, output rows, updated records)
  dsss_blocks.json  descramble, CRC-16, bytes, psdu_bytes, CCK codebooks, dbpsk_bits / find_sfds /
                    diff_bits / rebuild_scrambled on synthetic symbols
  scene_<x>.cf32 + scene_<x>.json   synthetic 802.11b captures (dsss_tx.scene, 30.72 MSps at 2437 MHz:
                    every rate and preamble, CFO, a weak and an excised frame) with candidates() of each
                    channel, the Barker MF of channel 6 at a few points, and dsss.scan()'s records
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
sys.path.insert(0, str(REPO / "python"))

from spectrum_intel.core.bands_india import WIFI_24  # noqa: E402
from spectrum_intel.core.dsp import IQContext  # noqa: E402
from spectrum_intel.io.sigmf import Capture  # noqa: E402
from spectrum_intel.receivers.wifi import droneid, dsss, links, mac80211  # noqa: E402
from spectrum_intel.receivers.wifi import dsss_tx as TX  # noqa: E402


def clean(o):
    if isinstance(o, dict):
        return {str(k) if not isinstance(k, str) else k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, np.ndarray):
        return clean(o.tolist())
    return o


def dump(name, obj):
    (HERE / name).write_text(json.dumps(clean(obj), indent=None, separators=(",", ":")))
    print(f"{name}: {(HERE / name).stat().st_size} bytes")


def mac(s):
    return bytes(int(x, 16) for x in s.split(":"))


# ----------------------------------------------------------------------------- MAC frames
def ie(tag, body):
    return bytes([tag, len(body)]) + bytes(body)


def mgmt(sub, da, sa, bssid, body, flags=0, seq=5):
    return bytes([sub << 4, flags, 0x3A, 0x01]) + mac(da) + mac(sa) + mac(bssid) + (seq << 4).to_bytes(2, "little") + body


def odid_pack():
    def msg(t, body):
        return bytes([(t << 4) | 2]) + bytes(body).ljust(24, b"\0")
    basic = msg(0, bytes([0x12]) + b"1581F45TB2211Z0A".ljust(20, b"\0"))
    lat, lon = int(28.6139 * 1e7), int(77.2090 * 1e7)
    loc = msg(1, bytes([0x20 | 0x04 | 0x02, 45, 40, 6]) + lat.to_bytes(4, "little", signed=True)
              + lon.to_bytes(4, "little", signed=True) + (2400).to_bytes(2, "little") + (2460).to_bytes(2, "little")
              + (2120).to_bytes(2, "little") + bytes([0x4B, 0x23]) + (1234).to_bytes(2, "little"))
    selfid = msg(3, bytes([0]) + b"survey flight")
    system = msg(4, bytes([0x05]) + int(28.61 * 1e7).to_bytes(4, "little", signed=True)
                 + int(77.20 * 1e7).to_bytes(4, "little", signed=True) + (1).to_bytes(2, "little") + bytes([0]) * 5
                 + bytes([0x12]) + (2300).to_bytes(2, "little") + (123456).to_bytes(4, "little"))
    opid = msg(5, bytes([0]) + b"IND-OP-0042")
    auth = msg(2, bytes([0x10, 2, 40]) + (99).to_bytes(4, "little"))
    msgs = [basic, loc, auth, selfid, system, opid]
    return bytes([7, 0xF2, 25, len(msgs)]) + b"".join(msgs)          # counter, pack header


def dji_payload(lat=28.5, lon=77.1):
    import struct
    s = struct.Struct("<BBBHH16siihhhhhhQiiiiBB20sH")
    sc = droneid.DJI_SCALE
    v = [0x58, 0x62, 2, 1234, 0x1f33, b"3NZCH8800A0XYZ", int(lon * sc), int(lat * sc), 300, 120, 5, -3, 1, 9000,
         1700000000000, int(28.49 * sc), int(77.09 * sc), int(77.08 * sc), int(28.48 * sc), 63, 6, b"uuid-1234", 0]
    raw = bytearray(s.pack(*v))
    crc = droneid.crc16_kermit(bytes(raw[:-2]))
    raw[-2:] = crc.to_bytes(2, "little")
    return bytes(raw)


def synthetic_frames():
    ap, sta, bc = "40:bb:56:00:6b:85", "dc:f5:05:e5:3d:39", "ff:ff:ff:ff:ff:ff"
    tsf = (123 * 86400 + 4 * 3600 + 5 * 60 + 6) * 1_000_000
    fixed = tsf.to_bytes(8, "little") + (100).to_bytes(2, "little") + (0x0431).to_bytes(2, "little")
    rsn = bytes([1, 0]) + bytes([0, 0x0f, 0xac, 4]) + (2).to_bytes(2, "little") + bytes([0, 0x0f, 0xac, 4, 0, 0x0f, 0xac, 2]) \
        + (3).to_bytes(2, "little") + bytes([0, 0x0f, 0xac, 2, 0, 0x0f, 0xac, 8, 0x00, 0x10, 0x18, 7]) + (0xC0).to_bytes(2, "little") \
        + (0).to_bytes(2, "little") + bytes([0, 0x0f, 0xac, 6])
    wps = b"".join(t.to_bytes(2, "big") + len(v).to_bytes(2, "big") + v for t, v in (
        (0x1044, b"\x02"), (0x1021, b"ZTE\x00"), (0x1023, b"ZXHN "), (0x1011, b"GateWay\xff\xfe"), (0x1047, bytes(range(16))),
        (0x1054, bytes([0, 6, 0, 0x50, 0xf2, 4, 0, 1])), (0x103B, b"\x03"), (0x1041, b"\x01"), (0x1012, b"\x00\x04")))
    p2p_dev = mac("02:00:00:00:01:02") + (0x0188).to_bytes(2, "big") + bytes([0, 10, 0, 0x50, 0xf2, 4, 0, 5]) + bytes([0]) \
        + (0x1011).to_bytes(2, "big") + (6).to_bytes(2, "big") + "Pixelé".encode()[:6]
    p2p = bytes([13]) + len(p2p_dev).to_bytes(2, "little") + p2p_dev + bytes([3, 6, 0]) + mac("02:00:00:00:09:09")
    elems = (ie(0, "Café-Ü 2.4G".encode()) + ie(1, [0x82, 0x84, 0x8b, 0x96, 0x0c, 0x12, 0x18, 0x24]) + ie(3, [6])
             + ie(5, [0, 3, 3, 0x1e, 0x80]) + ie(7, b"IN\x20" + bytes([1, 13, 20, 36, 4, 23, 250, 1, 1]))
             + ie(11, [3, 0, 88, 0x10, 0x27]) + ie(42, [5]) + ie(50, [0x30, 0x48, 0x60, 0x6c, 0xff])
             + ie(45, bytes([0x6f, 0x11, 0x17, 0xff, 0xff, 0, 0]) + bytes(19)) + ie(61, bytes([6, 0x07]) + bytes(20))
             + ie(48, rsn) + ie(54, [0x34, 0x12, 1]) + ie(70, [0] * 5) + ie(107, [0]) + ie(127, bytes([0, 0, 8, 0x80, 0, 0, 0, 0x40, 0xc0, 0x60]))
             + ie(191, bytes([0x32, 0x08, 0x88, 0x01, 0xfa, 0xff, 0, 0, 0xfa, 0xff, 0, 0]))
             + ie(192, [1, 42, 0]) + ie(221, bytes([0, 0x50, 0xf2, 1]) + bytes([1, 0, 0, 0x50, 0xf2, 2, 1, 0, 0, 0x50, 0xf2, 2, 1, 0, 0, 0x50, 0xf2, 2]))
             + ie(221, bytes([0, 0x50, 0xf2, 2, 1, 1])) + ie(221, bytes([0, 0x50, 0xf2, 4]) + wps)
             + ie(221, bytes([0x50, 0x6f, 0x9a, 9]) + p2p) + ie(221, bytes([0x50, 0x6f, 0x9a, 0x1c]) + mac(ap) + b"\x03abc")
             + ie(221, bytes([0x50, 0x6f, 0x9a, 0x10, 0])) + ie(221, bytes([0x50, 0x6f, 0x9a, 0x77]))
             + ie(221, bytes([0x00, 0x17, 0xf2, 0x0a, 0])) + ie(221, bytes([0x12, 0x34, 0x56, 1]))
             + ie(255, [35, 0]) + ie(255, [36, 0, 0, 0, 0x2a, 0, 0]) + ie(255, [108]) + ie(200, [1, 2]))
    frames = {
        "beacon_rich": mgmt(8, bc, ap, ap, fixed + elems),
        "probe_resp_retry": mgmt(5, sta, ap, ap, fixed + ie(0, b"") + ie(1, [2, 4, 11, 22]) + ie(3, [1]), flags=0x08),
        "beacon_hidden": mgmt(8, bc, ap, ap, fixed[:10] + (0x0011).to_bytes(2, "little") + ie(0, b"\0\0\0") + ie(1, [0x82])),
        "beacon_wep": mgmt(8, bc, ap, ap, fixed[:10] + (0x0011).to_bytes(2, "little") + ie(0, b"wep") + ie(1, [0x82])),
        "beacon_truncated_ie": mgmt(8, bc, ap, ap, fixed + ie(0, b"x") + bytes([1, 10, 2])),
        "beacon_malformed_rsn": mgmt(8, bc, ap, ap, fixed + ie(48, [1]) + ie(48, [1, 0, 0, 0x0f, 0xac, 4, 5, 0, 0, 0x0f]) +
                                     ie(221, bytes([0x50, 0x6f, 0x9a, 9, 13, 20, 0]) + bytes(10))),
        "beacon_sae_owe": mgmt(8, bc, ap, ap, fixed + ie(0, b"wpa3") + ie(48, bytes([1, 0, 0, 0x0f, 0xac, 4, 1, 0, 0, 0x0f, 0xac, 4, 2, 0, 0, 0x0f, 0xac, 8, 0, 0x0f, 0xac, 2, 0x40, 0]))),
        "beacon_owe": mgmt(8, bc, ap, ap, fixed + ie(48, bytes([1, 0, 0, 0x0f, 0xac, 4, 1, 0, 0, 0x0f, 0xac, 4, 1, 0, 0, 0x0f, 0xac, 18]))),
        "beacon_enterprise": mgmt(8, bc, ap, ap, fixed + ie(48, bytes([1, 0, 0, 0x0f, 0xac, 4, 1, 0, 0, 0x0f, 0xac, 4, 1, 0, 0, 0x0f, 0xac, 3]))),
        "beacon_odid": mgmt(8, bc, "fa:0b:bc:12:34:56", "fa:0b:bc:12:34:56", fixed + ie(0, b"") + ie(221, bytes([0xfa, 0x0b, 0xbc, 0x0d]) + odid_pack())),
        "beacon_dji": mgmt(8, bc, "60:60:1f:aa:bb:cc", "60:60:1f:aa:bb:cc", fixed + ie(0, b"DJI-MAVIC3") + ie(221, bytes([0x26, 0x37, 0x12]) + b"\x00\x00" + dji_payload())),
        "beacon_drone_ssid": mgmt(8, bc, "02:11:22:33:44:55", "02:11:22:33:44:55", fixed + ie(0, b"PHANTOM4_abc")),
        "beacon_dji_vendor": mgmt(8, bc, "60:60:1f:aa:bb:cc", "60:60:1f:aa:bb:cc", fixed + ie(0, b"Mavic-Air_123")),
        "probe_req": mgmt(4, bc, "da:a1:19:00:00:01", bc, ie(0, b"HomeNet") + ie(1, [2, 4, 11, 22, 12, 18]) + ie(45, bytes(26))
                          + ie(221, bytes([0x50, 0x6f, 0x9a, 9]) + p2p)),
        "probe_req_bcast": mgmt(4, bc, sta, bc, ie(0, b"") + ie(1, [2, 4])),
        "assoc_req": mgmt(0, ap, sta, ap, (0x0431).to_bytes(2, "little") + (10).to_bytes(2, "little") + ie(0, b"FTTH") + ie(1, [0x82, 0x84])),
        "reassoc_req": mgmt(2, ap, sta, ap, (0x0431).to_bytes(2, "little") + (10).to_bytes(2, "little") + mac(ap) + ie(0, b"FTTH")),
        "assoc_resp": mgmt(1, sta, ap, ap, (0x0431).to_bytes(2, "little") + (0).to_bytes(2, "little") + (0xC005).to_bytes(2, "little") + ie(1, [0x82])),
        "reassoc_resp_fail": mgmt(3, sta, ap, ap, (0x0431).to_bytes(2, "little") + (17).to_bytes(2, "little") + (0).to_bytes(2, "little")),
        "auth_open": mgmt(11, ap, sta, ap, bytes([0, 0, 1, 0, 0, 0])),
        "auth_sae": mgmt(11, ap, sta, ap, bytes([3, 0, 1, 0, 0, 0, 19, 0]) + bytes(32)),
        "auth_unknown_alg": mgmt(11, ap, sta, ap, bytes([9, 0, 2, 0, 99, 0])),
        "deauth": mgmt(12, sta, ap, ap, (7).to_bytes(2, "little")),
        "disassoc_unknown": mgmt(10, sta, ap, ap, (200).to_bytes(2, "little")),
        "action_p2p": mgmt(13, bc, sta, ap, bytes([4, 9, 0x50, 0x6f, 0x9a, 9, 1])),
        "action_vendor": mgmt(14, bc, sta, ap, bytes([127, 3, 0])),
        "action_unknown_cat": mgmt(13, bc, sta, ap, bytes([0x95, 3])),
        "mgmt_protected": mgmt(13, ap, sta, ap, bytes([1, 0, 0, 0x20, 0, 0, 0, 0]) + bytes(20), flags=0x40),
        "timing_adv": mgmt(6, bc, ap, ap, bytes(4)),
        "rts": bytes([0xB4, 0, 0x10, 0]) + mac(ap) + mac(sta),
        "cts": bytes([0xC4, 0, 0, 0]) + mac(sta),
        "ack": bytes([0xD4, 0x10, 0, 0]) + mac(sta),
        "ps_poll": bytes([0xA4, 0, 0x05, 0xC0]) + mac(ap) + mac(sta),
        "block_ack": bytes([0x94, 0, 0, 0]) + mac(sta) + mac(ap) + (0x5004).to_bytes(2, "little") + (0x0120).to_bytes(2, "little") + bytes([0xff] * 8),
        "block_ack_req": bytes([0x84, 0, 0, 0]) + mac(sta) + mac(ap) + (0x3004).to_bytes(2, "little") + (0x0230).to_bytes(2, "little"),
        "ctrl_short": bytes([0xD4, 0, 0, 0x80]),
        "too_short": bytes([0x08]),
        "ext_frame": bytes([0x1C, 0, 0, 0]) + mac(ap) + mac(sta) + bytes(10),
        "qos_null": bytes([0xC8, 0x11, 0x2C, 0]) + mac(ap) + mac(sta) + mac(ap) + (9 << 4).to_bytes(2, "little") + bytes([0x07, 0]),
        "data_ccmp": TX.data_frame(bssid=ap, sta=sta, n_payload=40, seq=3),
        "data_tkip": bytes([0x08, 0x41, 0, 0]) + mac(ap) + mac(sta) + mac("01:00:5e:00:00:fb") + (4 << 4).to_bytes(2, "little")
        + bytes([0x12, (0x12 | 0x20) & 0x7f, 0x34, 0x60, 1, 2, 3, 4]) + bytes(20),
        "data_wep": bytes([0x08, 0x41, 0, 0]) + mac(ap) + mac(sta) + mac(ap) + bytes(2) + bytes([1, 2, 3, 0x80]) + bytes(8),
        "data_enc_short": bytes([0x08, 0x41, 0, 0]) + mac(ap) + mac(sta) + mac(ap) + bytes(2) + bytes([1, 2, 3, 0x20, 5]),
        "data_eapol_m1": bytes([0x08, 0x02, 0, 0]) + mac(sta) + mac(ap) + mac(ap) + bytes(2) + bytes([0xaa, 0xaa, 3, 0, 0, 0, 0x88, 0x8e])
        + bytes([2, 3, 0, 95, 2, 0x00, 0x8a, 0, 16]) + (1).to_bytes(8, "big") + bytes(95 - 12),
        "data_eapol_m4": bytes([0x08, 0x01, 0, 0]) + mac(ap) + mac(sta) + mac(ap) + bytes(2) + bytes([0xaa, 0xaa, 3, 0, 0, 0, 0x88, 0x8e])
        + bytes([2, 3, 0, 95, 2, 0x03, 0x0a, 0, 16]) + (2).to_bytes(8, "big") + bytes(95 - 12),
        "data_eap_identity": bytes([0x08, 0x01, 0, 0]) + mac(ap) + mac(sta) + mac(ap) + bytes(2) + bytes([0xaa, 0xaa, 3, 0, 0, 0, 0x88, 0x8e])
        + bytes([1, 0, 0, 14, 2, 1, 0, 14, 1]) + b"user@x",
        "data_arp": bytes([0x08, 0x02, 0, 0]) + mac(bc) + mac(ap) + mac(sta) + bytes(2) + bytes([0xaa, 0xaa, 3, 0, 0, 0, 0x08, 0x06])
        + bytes([0, 1, 8, 0, 6, 4, 0, 1]) + mac(sta) + bytes([192, 168, 1, 5]) + bytes(6) + bytes([192, 168, 1, 1]),
        "data_ipv4_udp": bytes([0x88, 0x01, 0, 0]) + mac(ap) + mac(sta) + mac("01:00:5e:7f:ff:fa") + bytes(2) + bytes([0x80, 0])
        + bytes([0xaa, 0xaa, 3, 0, 0, 0, 0x08, 0x00]) + bytes([0x45, 0, 0, 40, 0, 0, 0, 0, 1, 17, 0, 0, 192, 168, 1, 5, 239, 255, 255, 250])
        + (50000).to_bytes(2, "big") + (1900).to_bytes(2, "big") + bytes(4),
        "data_ipv4_proto0": bytes([0x08, 0x01, 0, 0]) + mac(ap) + mac(sta) + mac(ap) + bytes(2)
        + bytes([0xaa, 0xaa, 3, 0, 0, 0, 0x08, 0x00]) + bytes([0x45, 0, 0, 40, 0, 0, 0, 0, 64, 0, 0, 0, 10, 0, 0, 1, 10, 0, 0, 2]),
        "data_ipv6": bytes([0x08, 0x02, 0, 0]) + mac("33:33:00:00:00:01") + mac(ap) + mac(sta) + bytes(2)
        + bytes([0xaa, 0xaa, 3, 0, 0, 0, 0x86, 0xdd]) + bytes([0x60, 0, 0, 0, 0, 8, 17, 255]) + bytes([0xfe, 0x80] + [0] * 13 + [1])
        + bytes([0xff, 2] + [0] * 13 + [0xfb]) + (5353).to_bytes(2, "big") + (5353).to_bytes(2, "big"),
        "data_wds": bytes([0x88, 0x03, 0, 0]) + mac(ap) + mac(sta) + mac("00:11:22:33:44:55") + bytes(2) + mac("00:aa:bb:cc:dd:ee")
        + bytes([0x80, 0]) + bytes(10),
        "data_ibss_order": bytes([0x88, 0x80, 0, 0]) + mac(sta) + mac(ap) + mac("00:11:22:33:44:55") + bytes(2) + bytes([0, 0]) + bytes(4)
        + bytes([0xaa, 0xaa, 3, 0, 0, 0, 0x88, 0xcc, 1, 2]),
        "data_other_llc": bytes([0x08, 0x02, 0, 0]) + mac(sta) + mac(ap) + mac(ap) + bytes(2) + bytes([0xaa, 0xaa, 3, 0, 0, 0, 0x12, 0x34, 9]),
        "null_frame": bytes([0x48, 0x11, 0, 0]) + mac(ap) + mac(sta) + mac(ap) + bytes(2),
    }
    return frames


def real_frames():
    out = []
    for meta, fc in (("DJI_Phantom_4_20251119-182450/D_20251119-182450.sigmf-meta", 2450e6),
                     ("D_20251119-183758/D_20251119-183758.sigmf-meta", 2420e6)):
        cap = Capture.load(REPO / "testdata" / meta) if hasattr(Capture, "load") else None
        x = cap.iq if cap is not None else np.fromfile(REPO / "testdata" / meta.replace("-meta", "-data"), np.complex64)
        for r in dsss.scan(IQContext(x, 61.44e6, fc)):
            if r.get("fcs_ok"):
                out.append(bytes.fromhex(r["mpdu_hex"]))
    return out


def gen_mac():
    cases = []
    for b in real_frames():
        cases.append({"name": "real", "hex": b.hex(), "parsed": mac80211.parse(b)})
    for name, b in synthetic_frames().items():
        cases.append({"name": name, "hex": bytes(b).hex(), "parsed": mac80211.parse(b)})
    vend = [{"mac": m, "vendor": mac80211.vendor(m)} for m in
            ("ff:ff:ff:ff:ff:ff", "01:00:5e:00:00:01", "da:a1:19:00:00:01", "40:bb:56:00:6b:85", "70:b3:d5:00:10:01",
             "8c:1f:64:00:10:01", "00:00:00:00:00:01", "fc:ff:ff:00:00:01", "60:60:1f:00:00:01")]
    names = [{"oui": o, "name": mac80211.oui_name(o)} for o in ("00:50:f2", "00:0F:AC", "60:60:1f", "12:34:56", "fa:0b:bc")]
    crc = [{"hex": b.hex(), "crc": __import__("zlib").crc32(b)} for b in (b"", b"123456789", bytes(range(256)))]
    dump("mac_frames.json", {"frames": cases, "vendor": vend, "oui_name": names, "crc32": crc,
                             "oui_db_size": len(mac80211._oui_db())})


def gen_droneid():
    pack = odid_pack()
    dji = dji_payload()
    cases = {
        "odid_pack": {"hex": pack.hex(), "out": droneid.parse_odid(pack)},
        "odid_single": {"hex": pack[4:29].hex(), "out": droneid.parse_odid(pack[4:29])},
        "odid_counter_single": {"hex": (bytes([3]) + pack[29:54]).hex(), "out": droneid.parse_odid(bytes([3]) + pack[29:54])},
        "odid_zero": {"hex": bytes(26).hex(), "out": droneid.parse_odid(bytes(26))},
        "odid_bad_version": {"hex": bytes([0x0F] + [1] * 24).hex(), "out": droneid.parse_odid(bytes([0x0F] + [1] * 24))},
        "dji": {"hex": (b"\x01\x02\x03" + dji).hex(), "out": droneid.parse_dji(b"\x01\x02\x03" + dji)},
        "dji_badcrc": {"hex": dji[:-1].hex() + "00", "out": droneid.parse_dji(dji[:-1] + b"\0")},
        "dji_short": {"hex": dji[:50].hex(), "out": droneid.parse_dji(dji[:50])},
    }
    for k, c in cases.items():
        c["summary"] = droneid.summary(c["out"]) if c["out"] else None
    sigs = [{"vendor": v, "ssid": s, "out": droneid.signature(v, s)} for v, s in (
        ("SZ DJI Technology Co.,Ltd", None), (None, "Mavic-Air_12"), ("SZ DJI Technology Co.,Ltd", "mavic_pro-1"),
        ("SZ DJI Technology Co.,Ltd", "PHANTOM"), ("Parrot SA", "anafi-1"), ("Parrot SA", "Bebop_2"),
        (None, "DJI-RC"), (None, "openhd"), (None, "tello-AB"), (None, "mini2_x"), (None, "Home"), ("Apple, Inc.", None),
        ("Autel Robotics Co., Ltd.", "autel-evo"), ("Hubsan", "H501S"), ("dji", "SKYDIO_2"))]
    names = sorted(set(mac80211._oui_db().values()))
    hits = [{"vendor": v, "out": droneid.signature(v, None)} for v in names if droneid.signature(v, None)]
    dump("droneid.json", {"cases": cases, "signature": sigs, "vendor_hits": hits, "n_vendor_names": len(names)})


def gen_links():
    rng = np.random.default_rng(3)
    A, B, C = "00:11:22:33:44:01", "00:11:22:33:44:02", "00:11:22:33:44:03"
    recs, t, seed_a = [], 1e-3, 40
    for i in range(30):
        kind = i % 5
        r = {"t0": t, "t1": t + 200e-6, "fmt": "HT" if i % 2 else "NON_HT", "power_db": float(-50 + rng.normal(0, 0.5)),
             "cfo_khz": float(12.0 + rng.normal(0, 0.1))}
        if kind in (0, 1, 2):
            r["scrambler_seed"], r["seed_ok"] = seed_a, True
            seed_a = seed_a % 127 + 1
            if kind == 0:
                r["mac"] = {"subtype": "qos-data", "ta": A, "ra": B, "ta_vendor": "Vendor A"}
        elif kind == 3:
            r["mac"] = {"subtype": "ack", "ra": A}
            r["power_db"] = float(-70 + rng.normal(0, 0.5))
        else:
            r["power_db"] = None if i == 9 else float(-61 + rng.normal(0, 0.5))
            r["scrambler_seed"], r["seed_ok"] = int(rng.integers(1, 128)), True
            if i == 14:
                r["mac"] = {"subtype": "beacon", "ta": C, "ra": "ff:ff:ff:ff:ff:ff", "addr2_vendor": "Vendor C"}
        recs.append(r)
        t = r["t1"] + (30e-6 if kind in (2,) else 400e-6)
    recs.append({"t0": t, "t1": t + 50e-6, "fmt": "NON_HT", "mac": {"subtype": "cts", "ra": B}, "power_db": -49.9})
    recs.append({"t0": t + 66e-6, "t1": t + 300e-6, "fmt": "VHT", "power_db": -60.8})
    shuffled = [recs[i] for i in rng.permutation(len(recs))]
    inp = json.loads(json.dumps(shuffled))
    rows = links.attribute(shuffled, sifs_us=16.0)
    dump("links.json", {"input": inp, "rows": rows, "records": shuffled})


def gen_blocks():
    rng = np.random.default_rng(7)
    bits = rng.integers(0, 2, 200).astype(np.uint8)
    hdr = TX.plcp_header(100, 11.0)
    hdr_bad = list(hdr)
    hdr_bad[5] ^= 1
    syms = np.exp(1j * (np.pi * np.cumsum(rng.integers(0, 2, 300)) + 0.01 * np.arange(301)[1:])) * (1 + 0.1 * rng.normal(size=300))
    syms = syms + 0.2 * (rng.normal(size=300) + 1j * rng.normal(size=300))
    # a DBPSK SYNC + SFD run, scrambled, as Barker symbols
    pre = [1] * 128 + [(0xF3A0 >> i) & 1 for i in range(16)] + list(rng.integers(0, 2, 60))
    sc = TX._scramble(pre, 0x6C)
    ph = np.pi * np.cumsum(sc)
    sfd_syms = np.exp(1j * (ph + 0.03 * np.arange(len(ph)))) + 0.1 * (rng.normal(size=len(ph)) + 1j * rng.normal(size=len(ph)))
    raw = dsss._dbpsk_bits(sfd_syms)
    desc = np.concatenate([np.full(128, 1, np.uint8), dsss.SFD_LONG])
    reb, agree = dsss._rebuild_scrambled(raw[:144], desc)
    out = {
        "descramble_in": bits, "descramble_out": dsss._descramble(bits),
        "crc_hdr": hdr, "crc_ok": dsss._crc16_ok(hdr[:32], hdr[32:48]), "crc_bad": hdr_bad,
        "crc_bad_ok": dsss._crc16_ok(hdr_bad[:32], hdr_bad[32:48]),
        "bytes_in": bits[:203 // 8 * 8], "bytes_out": list(dsss._bytes(bits)),
        "psdu_bytes": [[r, L, s, dsss.psdu_bytes(r, L, s)] for r in (1.0, 2.0, 5.5, 11.0) for L, s in ((800, 0x04), (1101, 0x84), (73, 0x80))],
        "cck11_re": dsss.CCK11.real, "cck11_im": dsss.CCK11.imag, "cck11_bits": dsss.CCK11_BITS,
        "cck55_re": dsss.CCK55.real, "cck55_im": dsss.CCK55.imag, "cck55_bits": dsss.CCK55_BITS,
        "syms_re": syms.real, "syms_im": syms.imag, "dbpsk": dsss._dbpsk_bits(syms),
        "diff1": dsss._diff_bits(syms, 0.013, 1), "diff2": dsss._diff_bits(syms, -0.2, 2),
        "sfd_re": sfd_syms.real, "sfd_im": sfd_syms.imag, "sfd_raw": raw,
        "sfds": [list(x) for x in dsss._find_sfds(dsss._descramble(raw))],
        "rebuilt": reb, "agree": agree,
    }
    dump("dsss_blocks.json", out)


def gen_scene(name, ppdus, dur_s, seed, fs=30.72e6, fc=2437e6, tone=None):
    n = int(dur_s * fs)
    x = TX.scene(ppdus, n, fs, fc, seed=seed)
    if tone:                                    # (t0, t1, f_hz, amplitude): a narrowband interferer
        k = np.arange(n)
        m = (k > int(tone[0] * fs)) & (k < int(tone[1] * fs))
        x = (x + m * tone[3] * np.exp(2j * np.pi * (tone[2] - fc) * k / fs)).astype(np.complex64)
    x.astype(np.complex64).tofile(HERE / f"scene_{name}.cf32")
    ctx = IQContext(x, fs, fc)
    cands = {}
    for ch, f in WIFI_24.items():
        lo, hi = fc - 0.47 * fs, fc + 0.47 * fs
        if f - 9e6 < lo or f + 9e6 > hi:
            continue
        y2, ofs2 = ctx.ch.slice(f, 11e6, dsss.FS)
        cands[str(ch)] = [list(c) for c in dsss._candidates(y2)]
    y4, ofs4 = ctx.ch.slice(WIFI_24[6], 11e6, dsss.FS_DEC)
    chan = dsss._Chan(y4, ofs4, WIFI_24[6])
    pts = list(range(1000, len(chan.c) - 10, 7919))
    tau, A = 4000.25, float(np.max(chan.ac[3990:4010]))
    s, taus, T = chan.track(tau, 60, A=A)
    out = dsss.scan(ctx)
    dump(f"scene_{name}.json", {
        "fs": fs, "fc": fc, "n": n, "candidates": cands,
        "mf_idx": pts, "c_re": chan.c[pts].real, "c_im": chan.c[pts].imag, "ym_re": chan.ym[pts].real,
        "ym_im": chan.ym[pts].imag, "track_tau": tau, "track_A": A, "track_s_re": s.real, "track_s_im": s.imag,
        "track_taus": taus, "track_T": T,
        "ppdus": out, "n_candidates": ctx.cache["wifi_dsss_candidates"], "unexplained": ctx.cache["wifi_dsss_unexplained"]})
    print(f"  scene {name}: {sum(1 for r in out if r.get('fcs_ok'))} FCS ok of {len(out)}: "
          + "; ".join(f"{r.get('rate_mbps')}/{r['preamble']}/{r.get('psdu_demod')}" for r in out))


def gen_scenes():
    f6 = WIFI_24[6]
    ap = "40:bb:56:00:6b:85"
    a = [dict(mpdu=TX.beacon(bssid=ap, ssid="test-1", channel=6), rate=1.0, preamble="long", f_hz=f6, t_s=0.1e-3, snr_db=14, cfo_hz=30e3),
         dict(mpdu=TX.data_frame(bssid=ap, n_payload=60, seq=1), rate=2.0, preamble="short", f_hz=f6, t_s=1.6e-3, snr_db=14, cfo_hz=-25e3),
         dict(mpdu=TX.data_frame(bssid=ap, n_payload=80, seq=2), rate=5.5, preamble="long", f_hz=f6, t_s=2.4e-3, snr_db=16, cfo_hz=40e3),
         dict(mpdu=TX.data_frame(bssid=ap, n_payload=120, seq=3), rate=11.0, preamble="short", f_hz=f6, t_s=3.0e-3, snr_db=18, cfo_hz=40e3),
         dict(mpdu=TX.data_frame(bssid=ap, n_payload=60, seq=4), rate=11.0, preamble="long", f_hz=f6, t_s=3.5e-3, snr_db=18, cfo_hz=-10e3),
         dict(mpdu=TX.data_frame(bssid=ap, n_payload=40, seq=5), rate=2.0, preamble="long", f_hz=f6, t_s=4.0e-3, snr_db=14, cfo_hz=5e3)]
    gen_scene("rates", a, 4.9e-3, seed=11)
    # weak frames: FCS-guided correction on the 11 Mb/s one, STFT excision of a tone on the 2 Mb/s one,
    # and a channel-7 beacon also seen (weaker) on channels 4..10
    b = [dict(mpdu=TX.beacon(bssid=ap, ssid="weak", channel=6), rate=1.0, preamble="long", f_hz=f6, t_s=0.1e-3, snr_db=-2.0),
         dict(mpdu=TX.data_frame(bssid=ap, n_payload=100, seq=6), rate=11.0, preamble="long", f_hz=f6, t_s=1.6e-3, snr_db=3.0,
              taps=((50.0, 0.4j), (120.0, -0.2))),
         dict(mpdu=TX.data_frame(bssid=ap, n_payload=100, seq=7), rate=5.5, preamble="short", f_hz=f6, t_s=2.3e-3, snr_db=1.0),
         dict(mpdu=TX.data_frame(bssid=ap, n_payload=150, seq=8), rate=2.0, preamble="long", f_hz=f6, t_s=2.9e-3, snr_db=-0.5),
         dict(mpdu=TX.beacon(bssid=ap, ssid="adj", channel=7), rate=1.0, preamble="long", f_hz=WIFI_24[7], t_s=4.0e-3, snr_db=12.0)]
    gen_scene("weak", b, 5.6e-3, seed=5, tone=(3.1e-3, 3.6e-3, f6 + 2.5e6, 3.0))


if __name__ == "__main__":
    which = set(sys.argv[1:]) or {"mac", "droneid", "links", "blocks", "scenes"}
    if "blocks" in which:
        gen_blocks()
    if "links" in which:
        gen_links()
    if "droneid" in which:
        gen_droneid()
    if "scenes" in which:
        gen_scenes()
    if "mac" in which:
        gen_mac()
