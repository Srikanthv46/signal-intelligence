"""ASTM F3411-22a / ASD-STAN prEN 4709-002 "Open Drone ID" (broadcast Remote ID): message parsing shared by every
transport (Wi-Fi Beacon vendor element in receivers.wifi; Bluetooth LE advertising in receivers.rc_link), and the
Bluetooth LE advertising-data (AD structure) parser.

Plain text, unauthenticated: the values are what the aircraft declares. Nothing is reported unless a message has a
known type, a valid protocol version and actual content. Pure stdlib.

Bluetooth transport (F3411-22a Annex, opendroneid-core-c, Apache-2.0, used as reference): AD type 0x16 (Service
Data, 16-bit UUID) with UUID 0xFFFA (little endian FA FF), application code 0x0D, then [message counter]
[one 25-byte message | message pack]: the same body as the Wi-Fi Beacon element after its OUI and type.
"""
from __future__ import annotations

import struct

from .drone_models import serial_manufacturer

ODID_OUI, ODID_VENDOR_TYPE = "fa:0b:bc", 0x0D

MSG_TYPE = {0: "basic_id", 1: "location", 2: "authentication", 3: "self_id", 4: "system",
            5: "operator_id", 0xF: "message_pack"}
ID_TYPE = {0: "none", 1: "serial number (CTA-2063-A)", 2: "CAA registration ID", 3: "UTM assigned UUID",
           4: "specific session ID"}
UA_TYPE = {0: "undeclared", 1: "aeroplane", 2: "multirotor", 3: "gyroplane", 4: "hybrid lift",
           5: "ornithopter", 6: "glider", 7: "kite", 8: "free balloon", 9: "captive balloon",
           10: "airship", 11: "free fall / parachute", 12: "rocket", 13: "tethered powered aircraft",
           14: "ground obstacle", 15: "other"}
OP_STATUS = {0: "undeclared", 1: "on ground", 2: "airborne", 3: "emergency", 4: "remote ID system failure"}
OP_LOC_TYPE = {0: "takeoff point", 1: "live GNSS (operator moving)", 2: "fixed location"}
DESC_TYPE = {0: "text", 1: "emergency", 2: "extended status"}
HORIZ_ACC = {0: None, 1: "<10 NM", 2: "<4 NM", 3: "<2 NM", 4: "<1 NM", 5: "<0.5 NM", 6: "<0.3 NM",
             7: "<0.1 NM", 8: "<0.05 NM", 9: "<30 m", 10: "<10 m", 11: "<3 m", 12: "<1 m"}
VERT_ACC = {0: None, 1: "<150 m", 2: "<45 m", 3: "<25 m", 4: "<10 m", 5: "<3 m", 6: "<1 m"}
SPEED_ACC = {0: None, 1: "<10 m/s", 2: "<3 m/s", 3: "<1 m/s", 4: "<0.3 m/s"}
UA_CLASS_EU = {0: "undeclared", 1: "open", 2: "specific", 3: "certified"}

_INVALID_LATLON = 0                            # 0,0 means "not reported" in F3411
MSG_LEN = 25


# --------------------------------------------------------------------------- Open Drone ID (ASTM F3411)
def _s32(b, i):
    return struct.unpack_from("<i", b, i)[0]


def _u16(b, i):
    return struct.unpack_from("<H", b, i)[0]


def _text(b):
    return b.split(b"\x00", 1)[0].decode("utf-8", "replace").strip()


def _alt(raw):
    """uint16 -> metres (0.5 m steps, -1000 m offset); the all-zero code means not reported."""
    return None if raw == 0 else round(raw * 0.5 - 1000, 1)


def _latlon(v):
    return None if v == _INVALID_LATLON else round(v * 1e-7, 7)


def _basic_id(b):
    out = {"id_type": ID_TYPE.get(b[0] >> 4, b[0] >> 4), "ua_type": UA_TYPE.get(b[0] & 0xF, b[0] & 0xF)}
    uas = _text(b[1:21])
    if uas:
        out["uas_id"] = uas
        mf = serial_manufacturer(uas) if b[0] >> 4 == 1 else None
        if mf:                                    # CTA-2063-A manufacturer code (research, PLAN M2)
            out["manufacturer"] = mf
    return out


def _location(b):
    st = b[0]
    mult, ew = st & 1, (st >> 1) & 1
    out = {"status": OP_STATUS.get(st >> 4, st >> 4),
           "height_type": "above ground" if (st >> 2) & 1 else "above takeoff"}
    if b[1] != 181:                                             # 181 = unknown direction
        out["track_deg"] = b[1] + 180 * ew
    if b[2] != 255:
        out["speed_mps"] = round(b[2] * 0.75 + 63.75 if mult else b[2] * 0.25, 2)
    vs = struct.unpack_from("<b", b, 3)[0]
    if vs != 63:
        out["vertical_speed_mps"] = round(vs * 0.5, 2)
    lat, lon = _latlon(_s32(b, 4)), _latlon(_s32(b, 8))
    if lat is not None and lon is not None and abs(lat) <= 90 and abs(lon) <= 180:
        out["lat"], out["lon"] = lat, lon
    for k, i in (("pressure_altitude_m", 12), ("geodetic_altitude_m", 14), ("height_m", 16)):
        v = _alt(_u16(b, i))
        if v is not None:
            out[k] = v
    acc = {"horizontal": HORIZ_ACC.get(b[18] & 0xF), "vertical": VERT_ACC.get(b[18] >> 4),
           "speed": SPEED_ACC.get(b[19] & 0xF), "barometric": VERT_ACC.get(b[19] >> 4)}
    acc = {k: v for k, v in acc.items() if v}
    if acc:
        out["accuracy"] = acc
    ts = _u16(b, 20)
    if ts != 0xFFFF:
        out["timestamp_s_into_hour"] = round(ts * 0.1, 1)
    return out


def _auth(b):
    out = {"auth_type": b[0] >> 4, "page": b[0] & 0xF}
    if out["page"] == 0 and len(b) >= 7:
        out["page_count"] = b[1]
        out["length"] = b[2]
        out["timestamp"] = struct.unpack_from("<I", b, 3)[0]
    return out


def _self_id(b):
    out = {"description_type": DESC_TYPE.get(b[0], b[0])}
    d = _text(b[1:24])
    if d:
        out["description"] = d
    return out


def _system(b):
    out = {"operator_location_type": OP_LOC_TYPE.get(b[0] & 3, b[0] & 3),
           "classification_type": "EU" if (b[0] >> 2) & 7 == 1 else (b[0] >> 2) & 7}
    lat, lon = _latlon(_s32(b, 1)), _latlon(_s32(b, 5))
    if lat is not None and lon is not None and abs(lat) <= 90 and abs(lon) <= 180:
        out["operator_lat"], out["operator_lon"] = lat, lon
    n = _u16(b, 9)
    if n > 1:
        out["area_count"] = n
        out["area_radius_m"] = b[11] * 10
    if out["classification_type"] == "EU":
        out["ua_class_eu"] = UA_CLASS_EU.get(b[16] >> 4, b[16] >> 4)
        out["ua_category_eu"] = b[16] & 0xF
    v = _alt(_u16(b, 17))
    if v is not None:
        out["operator_altitude_m"] = v
    ts = struct.unpack_from("<I", b, 19)[0]
    if ts:
        out["timestamp_s_since_2019"] = ts
    return out


def _operator_id(b):
    out = {"operator_id_type": "CAA-issued" if b[0] == 0 else b[0]}
    v = _text(b[1:21])
    if v:
        out["operator_id"] = v
    return out


_BODY = {0: _basic_id, 1: _location, 2: _auth, 3: _self_id, 4: _system, 5: _operator_id}


def _message(m: bytes) -> dict | None:
    """One 25-byte Open Drone ID message -> {type, protocol_version, ...fields}."""
    if len(m) < MSG_LEN:
        return None
    t, ver = m[0] >> 4, m[0] & 0xF
    if t not in _BODY or ver > 3:                               # F3411 versions 0-2 (+1 headroom)
        return None
    try:
        out = _BODY[t](m[1:MSG_LEN])
    except (IndexError, struct.error, ValueError):
        return None
    out["message_type"] = MSG_TYPE[t]
    out["protocol_version"] = ver
    return out


# A message only counts as Remote ID if it actually says something. An all-zero element parses as a
# syntactically valid but empty Basic ID, and reporting that as a drone would be a guess.
_CONTENT = ("uas_id", "lat", "operator_lat", "operator_id", "description", "speed_mps", "height_m",
            "geodetic_altitude_m", "track_deg")


def _informative(msgs: list[dict]) -> bool:
    for m in msgs:
        if any(m.get(k) is not None for k in _CONTENT):
            return True
        if m.get("status") not in (None, "undeclared") or m.get("ua_type") not in (None, "undeclared"):
            return True
    return False


# --------------------------------------------------------------------------- China GB 46750-2025
# (research, PLAN G14; research/drone_transmitters/16_uat_adsl_fanet_ogn.md section E): one variable-length packet,
# type byte 0xFF, version (bits 7-5 = 001), content length, flag bytes (7 presence bits each, bit 0 = another flag byte
# follows), then the present items in ID order at fixed lengths; little-endian; positions longitude first.
GB_ITEMS = [  # (id, bytes, key)
    (1, 20, "uas_id"), (2, 8, "registration_mark"), (3, 1, "operation_category"), (4, 1, "ua_class"),
    (5, 1, "operator_location_type"), (6, 8, "operator_position"), (7, 2, "operator_altitude_m"),
    (8, 8, "position"), (9, 2, "track_deg"), (10, 2, "speed_mps"), (11, 2, "height_m"), (12, 1, "vertical_speed_mps"),
    (13, 2, "geodetic_altitude_m"), (14, 2, "pressure_altitude_m"), (15, 1, "status"), (16, 1, "coordinate_system"),
    (17, 1, "horizontal_accuracy"), (18, 1, "vertical_accuracy"), (19, 1, "speed_accuracy"), (20, 6, "timestamp_ms"),
    (21, 1, "timestamp_accuracy")]
GB_CATEGORY = ("undefined", "open", "specific", "certified")
GB_CLASS = ("micro", "light", "small", "medium", "large")
GB_STATUS = ("not reported", "ground", "airborne", "emergency", "RID failure", "RID failure + emergency")


def _gb_lonlat(b):
    lon, lat = _s32(b, 0), _s32(b, 4)
    if lon == -1 or lat == -1:
        return None
    lon, lat = lon * 1e-7, lat * 1e-7
    return (round(lat, 7), round(lon, 7)) if abs(lat) <= 90 and abs(lon) <= 180 else None


def parse_gb46750(b: bytes, transport: str = "Wi-Fi Beacon vendor element") -> dict | None:
    """GB 46750-2025 packet (starting at the 0xFF type byte) -> fields, or None when it is not well formed (version,
    length, flag chain and item lengths must agree, and the mandatory serial must be present)."""
    b = bytes(b)
    if len(b) < 6 or b[0] != 0xFF or b[1] >> 5 != 1:
        return None
    ln = b[2]
    flags, i = [], 3
    while i < len(b):
        flags.append(b[i])
        i += 1
        if not b[i - 1] & 1:
            break
    present = [k for k in range(1, 8 * len(flags)) if k - 1 < 7 * len(flags)
               and flags[(k - 1) // 7] & (0x80 >> ((k - 1) % 7))]
    items = [it for it in GB_ITEMS if it[0] in present]
    if sum(n for _, n, _ in items) != ln or len(b) < i + ln or 1 not in present:
        return None
    out: dict = {"standard": "China GB 46750-2025", "transport": transport, "version": f"{b[1] >> 5}.{b[1] & 0x1F}",
                 "messages": ["GB 46750"]}
    for iid, n, key in items:
        v = b[i:i + n]
        i += n
        u = int.from_bytes(v, "little")
        if key in ("uas_id", "registration_mark"):
            t = _text(v)
            if t and t.strip("0"):                             # an all-'0' registration mark is a placeholder
                out[key] = t
        elif key == "operation_category":
            out[key] = GB_CATEGORY[u] if u < 4 else u
        elif key == "ua_class":
            out[key] = GB_CLASS[u] if u < 5 else u
        elif key == "operator_location_type":
            out[key] = "take-off point" if u == 0 else "remote-station position"
        elif key in ("operator_position", "position"):
            ll = _gb_lonlat(v)
            if ll:
                if key == "position":
                    out["lat"], out["lon"] = ll
                else:
                    out["operator_lat"], out["operator_lon"] = ll
        elif key in ("operator_altitude_m", "geodetic_altitude_m", "pressure_altitude_m"):
            if u:
                out[key] = round(u / 2 - 1000, 1)
        elif key == "height_m":
            if u:
                out[key] = round(u / 2 - 9000, 1)
        elif key == "track_deg":
            if u != 0xFFFF:
                out[key] = u / 10
        elif key == "speed_mps":
            if u != 0xFFFF:
                out[key] = u / 10
        elif key == "vertical_speed_mps":
            if u != 0xFF:
                out[key] = round((u & 0x7F) * 0.5 * (-1 if u & 0x80 else 1), 1)
        elif key == "status":
            out[key] = GB_STATUS[u] if u < len(GB_STATUS) else u
        elif key == "coordinate_system":
            out[key] = "CGCS2000" if u == 1 else "WGS-84"
        else:
            out[key] = u
    if out.get("uas_id"):
        mf = serial_manufacturer(out["uas_id"])
        if mf:
            out["manufacturer"] = mf
    return out


def parse_odid(body: bytes, transport: str = "Wi-Fi Beacon vendor element") -> dict | None:
    """Open Drone ID vendor-element body (after OUI fa:0b:bc and vendor type 0x0d) -> merged fields.

    Wi-Fi Beacon transport (F3411-22a Annex E): [message counter][message pack | single message].
    Returns None unless at least one message parses AND carries content, so neither a coincidental
    OUI match nor an all-zero element is ever reported as an aircraft.
    """
    body = bytes(body)
    for start in (1, 0):                                        # China GB 46750-2025 in the same element (0xFF first)
        g = parse_gb46750(body[start:], transport)
        if g:
            if start:
                g["message_counter"] = body[0]
            return g
    msgs, counter = [], None
    for start in (1, 0):                                        # with and without the message counter
        b = body[start:]
        if len(b) < MSG_LEN:
            continue
        if b[0] >> 4 == 0xF and len(b) >= 3 and b[1] == MSG_LEN:   # message pack
            n = b[2]
            got = [_message(b[3 + MSG_LEN * k:3 + MSG_LEN * (k + 1)]) for k in range(n)]
            got = [g for g in got if g]
            if got and _informative(got):
                msgs, counter = got, (body[0] if start else None)
                break
        got = [_message(b[MSG_LEN * k:MSG_LEN * (k + 1)]) for k in range(len(b) // MSG_LEN)]
        got = [g for g in got if g]
        if got and _informative(got):
            msgs, counter = got, (body[0] if start else None)
            break
    if not msgs:
        return None
    out: dict = {"standard": "ASTM F3411 / ASD-STAN prEN 4709-002 (Open Drone ID)",
                 "transport": transport, "messages": []}
    if counter is not None:
        out["message_counter"] = counter
    for m in msgs:
        out["messages"].append(m["message_type"])
        for k, v in m.items():
            if k not in ("message_type", "protocol_version"):
                out.setdefault(k, v)
        out["protocol_version"] = m["protocol_version"]
    return out



# --------------------------------------------------------------------------- Bluetooth LE advertising data
ODID_BLE_UUID = 0xFFFA
ODID_BLE_APP_CODE = 0x0D
AD_TYPES = {0x01: "flags", 0x02: "uuid16_incomplete", 0x03: "uuid16_complete", 0x06: "uuid128_incomplete",
            0x07: "uuid128_complete", 0x08: "short_name", 0x09: "local_name", 0x0A: "tx_power_dbm",
            0x16: "service_data_uuid16", 0x19: "appearance", 0xFF: "manufacturer_data"}
COMPANY_IDS = {0x004C: "Apple", 0x0006: "Microsoft", 0x0075: "Samsung", 0x00E0: "Google", 0x0059: "Nordic Semiconductor",
               0x02E5: "Espressif", 0x000D: "Texas Instruments", 0x0131: "Cypress", 0x0157: "Huami", 0x038F: "Xiaomi",
               0x0822: "DJI"}


def parse_ble_ad(adv_data: bytes) -> dict:
    """Advertising data (after AdvA) -> AD fields; a Remote ID service-data element is parsed with parse_odid.
    Malformed / truncated structures stop the walk (what was parsed before is kept)."""
    out: dict = {}
    b = bytes(adv_data)
    i = 0
    while i < len(b):
        ln = b[i]
        if ln == 0 or i + 1 + ln > len(b):
            break
        t, v = b[i + 1], b[i + 2:i + 1 + ln]
        i += 1 + ln
        if t == 0x01 and v:
            out["flags"] = v[0]
        elif t in (0x02, 0x03) and len(v) >= 2:
            out.setdefault("service_uuids16", []).extend(f"{struct.unpack_from('<H', v, k)[0]:04X}" for k in range(0, len(v) - 1, 2))
        elif t in (0x08, 0x09):
            out["local_name"] = v.decode("utf-8", "replace")
        elif t == 0x0A and v:
            out["tx_power_dbm"] = struct.unpack("<b", v[:1])[0]
        elif t == 0x19 and len(v) >= 2:
            out["appearance"] = struct.unpack_from("<H", v, 0)[0]
        elif t == 0xFF and len(v) >= 2:
            cid = struct.unpack_from("<H", v, 0)[0]
            out["company_id"] = f"0x{cid:04X}"
            if cid in COMPANY_IDS:
                out["company"] = COMPANY_IDS[cid]
            out["manufacturer_data_hex"] = v[2:].hex()
        elif t == 0x16 and len(v) >= 2:
            uuid = struct.unpack_from("<H", v, 0)[0]
            out.setdefault("service_data", []).append({"uuid16": f"{uuid:04X}", "data_hex": v[2:].hex()})
            if uuid == ODID_BLE_UUID and len(v) >= 3 and v[2] == ODID_BLE_APP_CODE:
                rid = parse_odid(v[3:], transport="Bluetooth LE advertising (legacy)")
                if rid:
                    out["remote_id"] = rid
    from .ble_trackers import recognize           # item trackers / notable devices (research, signal gaps S7)
    dev = recognize(out)
    if dev:
        out["device"] = dev
    return out
