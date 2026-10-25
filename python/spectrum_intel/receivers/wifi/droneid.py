"""Drone identification broadcast inside 802.11 frames.

Two independent things live here, both read straight out of a vendor-specific element:

  parse_odid(body)   ASTM F3411-22a / ASD-STAN prEN 4709-002 "Open Drone ID" (OUI fa:0b:bc, type 0x0d).
                     The open regulatory standard: UAS serial, position, altitude, speed, and the
                     OPERATOR's position. Plain-text, no authentication, so it is simply read.
  parse_dji(body)    DJI's proprietary 91-byte DroneID payload (same structure DJI sends over OcuSync;
                     layout from ../droneid-decoder-reference). CRC-16/Kermit protected, so it is
                     located by checksum rather than by an assumed offset -- a parse means the CRC passed.

  signature(vendor, ssid)   platform guess from the IEEE OUI name and the SSID, "-like" when inferred.

Nothing is reported unless it is checked: an Open Drone ID message must carry a known type and a valid
protocol version, and a DJI payload must pass its CRC. Pure stdlib (no numpy), so mac80211 can call it.
"""
from __future__ import annotations

import re
import struct

from spectrum_intel.core.remote_id import ODID_OUI, ODID_VENDOR_TYPE, _text, parse_odid  # moved to core (shared with BLE)

DJI_OUIS = ("26:37:12",)                       # DJI DroneID vendor element seen on Wi-Fi


# --------------------------------------------------------------------------- DJI DroneID (proprietary)
_DJI_STRUCT = struct.Struct("<BBBHH16siihhhhhhQiiiiBB20sH")
DJI_LEN = _DJI_STRUCT.size                                      # 91
DJI_CRC_INIT = 0x3692
DJI_SCALE = 174533.0                                            # DJI's radian-ish lat/lon fixed point

DJI_MODELS = {1: "Inspire 1", 2: "Phantom 3 Series", 3: "Phantom 3 Series", 4: "Phantom 3 Std", 5: "M100",
              6: "ACEONE", 7: "WKM", 8: "NAZA", 9: "A2", 10: "A3", 11: "Phantom 4", 12: "MG1", 14: "M600",
              15: "Phantom 3 4k", 16: "Mavic Pro", 17: "Inspire 2", 18: "Phantom 4 Pro", 20: "N2",
              21: "Spark", 23: "M600 Pro", 24: "Mavic Air", 25: "M200", 26: "Phantom 4 Series",
              27: "Phantom 4 Adv", 28: "M210", 30: "M210RTK", 31: "A3_AG", 32: "MG2", 34: "MG1A",
              35: "Phantom 4 RTK", 36: "Phantom 4 Pro V2.0", 38: "MG1P", 40: "MG1P-RTK", 41: "Mavic 2",
              44: "M200 V2 Series", 51: "Mavic 2 Enterprise", 53: "Mavic Mini", 58: "Mavic Air 2",
              59: "P4M", 60: "M300 RTK", 61: "DJI FPV", 63: "Mini 2", 64: "AGRAS T10", 65: "AGRAS T30",
              66: "Air 2S", 67: "M30", 68: "DJI Mavic 3", 69: "Mavic 2 Enterprise Advanced", 70: "Mini SE"}


def crc16_kermit(data: bytes, crc: int = DJI_CRC_INIT) -> int:
    """CRC-16/Kermit as DJI uses it (verified against ../droneid-decoder-reference fast_kernels)."""
    for b in data:
        crc ^= b
        for _ in range(8):
            crc = (crc >> 1) ^ (0x8408 if crc & 1 else 0)
    return crc


def _dji_at(raw: bytes) -> dict | None:
    v = _DJI_STRUCT.unpack(raw)
    if crc16_kermit(raw[:-2]) != v[22]:
        return None
    lat, lon = v[7] / DJI_SCALE, v[6] / DJI_SCALE
    out = {"standard": "DJI DroneID (proprietary)", "transport": "Wi-Fi vendor element",
           "crc_ok": True, "version": v[2], "sequence": v[3], "state_info": f"0x{v[4]:04x}",
           "uas_id": _text(v[5]), "model": DJI_MODELS.get(v[19], f"unknown type {v[19]}"),
           "ua_type": "multirotor"}
    if abs(lat) <= 90 and abs(lon) <= 180 and (lat or lon):
        out["lat"], out["lon"] = round(lat, 7), round(lon, 7)
    out["geodetic_altitude_m"] = round(v[8] / 3.281, 2)
    out["height_m"] = round(v[9] / 3.281, 2)
    out["velocity_ne_up_mps"] = [v[10], v[11], v[12]]
    out["yaw_deg"] = round(v[13] / 100.0, 2) if abs(v[13]) <= 36000 else None
    if v[14]:
        out["gps_time_ms"] = v[14]
    olat, olon = v[15] / DJI_SCALE, v[16] / DJI_SCALE
    if abs(olat) <= 90 and abs(olon) <= 180 and (olat or olon):
        out["operator_lat"], out["operator_lon"] = round(olat, 7), round(olon, 7)
        out["operator_location_type"] = "live GNSS (operator moving)"
    hlat, hlon = v[18] / DJI_SCALE, v[17] / DJI_SCALE
    if abs(hlat) <= 90 and abs(hlon) <= 180 and (hlat or hlon):
        out["home_lat"], out["home_lon"] = round(hlat, 7), round(hlon, 7)
    uuid = _text(v[21])
    if uuid:
        out["uuid"] = uuid
    return {k: v for k, v in out.items() if v is not None}


def parse_dji(body: bytes) -> dict | None:
    """Find DJI's 91-byte DroneID payload inside a vendor element and parse it.

    The payload is located by its CRC-16/Kermit, not by an assumed offset, so a result is always a
    frame that checked out. Returns None when no window passes.
    """
    body = bytes(body)
    for off in range(0, max(0, len(body) - DJI_LEN) + 1):
        try:
            got = _dji_at(body[off:off + DJI_LEN])
        except (struct.error, ValueError):
            continue
        if got:
            got["payload_offset"] = off
            return got
    return None


# --------------------------------------------------------------------------- from a vendor element
def from_vendor_ie(oui: str, vendor_type: int, body: bytes) -> dict | None:
    """Vendor-specific element -> Remote ID, or None. `body` is everything after OUI + vendor type."""
    o = (oui or "").lower()
    if o == ODID_OUI and vendor_type == ODID_VENDOR_TYPE:
        return parse_odid(body)
    if o in DJI_OUIS:
        return parse_dji(bytes([vendor_type]) + bytes(body))
    return None


# ASTM F3411 over Wi-Fi NAN (research, PLAN T3; F3411-22a / opendroneid-core-c wifi.c): a public action frame (category
# 4, vendor specific 9, OUI 50:6f:9a, type 0x13 = NAN service discovery) whose Service Descriptor Attribute (0x03) has
# the service ID SHA-256("org.opendroneid.remoteid")[:6] and carries [message counter][message pack] as service info.
NAN_ODID_SERVICE_ID = bytes.fromhex("8869199d9209")


def from_nan_action(b: bytes) -> dict | None:
    """802.11 action frame body (from the category byte) -> Remote ID, or None."""
    b = bytes(b)
    if len(b) < 7 or b[0] != 4 or b[1] != 9 or b[2:5] != bytes([0x50, 0x6F, 0x9A]) or b[5] != 0x13:
        return None
    i = 6
    while i + 3 <= len(b):
        aid, ln = b[i], b[i + 1] | (b[i + 2] << 8)
        v = b[i + 3:i + 3 + ln]
        i += 3 + ln
        if aid != 0x03 or len(v) < 9 or v[:6] != NAN_ODID_SERVICE_ID:
            continue
        ctrl, j = v[8], 9
        if ctrl & 0x40:                                      # binding bitmap
            j += 2
        if ctrl & 0x04 and j < len(v):                       # matching filter
            j += 1 + v[j]
        if ctrl & 0x08 and j < len(v):                       # service response filter
            j += 1 + v[j]
        if ctrl & 0x10 and j < len(v):                       # service info
            info = v[j + 1:j + 1 + v[j]]
            rid = parse_odid(info, transport="Wi-Fi NAN service discovery")
            if rid:
                return rid
    return None


def scan_vendor_ies(tlv: bytes) -> dict | None:
    """Last resort: walk raw elements looking for Remote ID (used on frames that are not beacons)."""
    i = 0
    tlv = bytes(tlv)
    while i + 2 <= len(tlv):
        ln = tlv[i + 1]
        v = tlv[i + 2:i + 2 + ln]
        if tlv[i] == 221 and len(v) >= 4:
            oui = ":".join(f"{x:02x}" for x in v[0:3])
            got = from_vendor_ie(oui, v[3], v[4:])
            if got:
                return got
        i += 2 + ln
    return None


# --------------------------------------------------------------------------- platform signatures
# Registry names (IEEE OUI database) of companies that build drones. Matched on the vendor string so the
# prefixes stay in one place (data/oui.tsv.gz) instead of being duplicated here.
VENDOR_PATTERNS = (
    (r"\bdji\b|\bsz dji\b|dji (osmo|baiwang|ronin)", "DJI"),
    (r"\bparrot\b", "Parrot"),
    (r"autel (robotics|lntelligent|intelligent)", "Autel Robotics"),
    (r"\bskydio\b", "Skydio"),
    (r"\byuneec\b", "Yuneec"),
    (r"\bhubsan\b", "Hubsan"),
    (r"powervision tech", "PowerVision"),
    (r"\bwalkera\b", "Walkera"),
    (r"\behang\b", "EHang"),
)
# SSID patterns drones and their controllers broadcast. Anchored so ordinary networks do not match.
SSID_PATTERNS = (
    (r"^mavic[_-]?air[_-]", "DJI Mavic Air"),
    (r"^mavic[_-]?mini[_-]", "DJI Mavic Mini"),
    (r"^phantom3[_-]", "DJI Phantom 3 Standard / 4K (Wi-Fi video models)"),
    (r"^(mavic|mavic[_-]?(air|mini|pro|2|3))[_-]", "DJI Mavic"),
    (r"^phantom[_-]?\d*", "DJI Phantom"),
    (r"^spark[_-]rc", "DJI Spark remote controller"),
    (r"^spark[_-]", "DJI Spark"),
    (r"^tello[_-]", "DJI Tello"),
    (r"^dji[_-]", "DJI"),
    (r"^(mini)[_-]?[23][_-]", "DJI Mini"),
    (r"^anafi[_-]", "Parrot Anafi"),
    (r"^(bebop|anafi|disco|parrot)[_-]", "Parrot"),
    (r"^autel[_-]", "Autel Robotics"),
    (r"^cgo3p?_", "Yuneec CGO3 camera video link (Typhoon Q500 4K / Typhoon H)"),
    (r"^(typhoon|yuneec)[_-]", "Yuneec"),
    (r"^(hubsan|h501|h107)", "Hubsan"),
    (r"^skydio[_-]", "Skydio"),
    (r"^(fimi|xiaomi[_-]fimi)[_-]", "FIMI"),
    (r"^(openhd|wfb|wifibroadcast)", "open-source FPV link (OpenHD / wfb-ng)"),
    # PLAN W3 (research file 18): vendor-documented default hotspot names
    (r"^rmtt-", "DJI RoboMaster TT (Tello Talent)"),
    (r"^holystone ?fpv", "Holy Stone (Wi-Fi FPV: HS720 / HS100 series)"),
    (r"^ruko[_-]f11", "Ruko F11"),
    (r"^ruko[_-]", "Ruko"),
    (r"^wifi-720p-", "toy Wi-Fi FPV drone (Eachine E58 class)"),
    (r"^expresslrs (tx|rx)$", "ExpressLRS radio in Wi-Fi setup mode"),
    (r"^(ardupilot|pixracer)$", "ArduPilot Wi-Fi telemetry bridge (MAVESP8266)"),
    (r"^cuavwlink", "CUAV PW-Link Wi-Fi telemetry"),
    (r"^dronebridge", "DroneBridge ESP32 telemetry"),
)


def signature(vendor_name: str | None = None, ssid: str | None = None) -> dict | None:
    """Platform guess from the IEEE vendor name and/or the SSID. Returns None when nothing matches.

    `confidence` stays below 1: a drone-vendor OUI or a drone-shaped SSID is strong evidence of the
    platform, never proof that the thing is airborne. Remote ID (above) is the authoritative source.
    """
    out: dict = {}
    if vendor_name:
        v = vendor_name.lower()
        for pat, name in VENDOR_PATTERNS:
            if re.search(pat, v):
                out = {"platform": name, "evidence": f"IEEE OUI registered to {vendor_name}",
                       "confidence": 0.75}
                break
    if ssid:
        s = ssid.lower()
        for pat, name in SSID_PATTERNS:
            if re.search(pat, s):
                if out and out["platform"].split()[0] == name.split()[0]:
                    out["platform"] = name if len(name) > len(out["platform"]) else out["platform"]
                    out["evidence"] += f"; SSID '{ssid}' matches {name}"
                    out["confidence"] = 0.9
                elif not out:
                    out = {"platform": name, "evidence": f"SSID '{ssid}' matches {name}",
                           "confidence": 0.7}
                break
    return out or None


# --------------------------------------------------------------------------- summary
def summary(rid: dict) -> str:
    """One line for an operator: what the aircraft said about itself."""
    s = ["REMOTE ID"]
    if rid.get("uas_id"):
        s.append(f"UAS {rid['uas_id']}")
    if rid.get("model"):
        s.append(f"({rid['model']})")
    elif rid.get("ua_type") and rid["ua_type"] != "undeclared":
        s.append(f"({rid['ua_type']})")
    if rid.get("status"):
        s.append(rid["status"].upper())
    if rid.get("lat") is not None:
        s.append(f"at {rid['lat']:.6f},{rid['lon']:.6f}")
    for k, lab in (("height_m", "AGL"), ("geodetic_altitude_m", "alt")):
        if rid.get(k) is not None:
            s.append(f"{rid[k]:.0f} m {lab}")
            break
    if rid.get("speed_mps") is not None:
        s.append(f"{rid['speed_mps']:.1f} m/s")
    if rid.get("track_deg") is not None:
        s.append(f"hdg {rid['track_deg']}")
    if rid.get("operator_lat") is not None:
        s.append(f"OPERATOR at {rid['operator_lat']:.6f},{rid['operator_lon']:.6f}")
    if rid.get("operator_id"):
        s.append(f"operator ID {rid['operator_id']}")
    if rid.get("description"):
        s.append(f"'{rid['description']}'")
    return " ".join(s)
