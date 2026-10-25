"""Recognise item trackers and notable devices from a decoded BLE advertisement (research, signal gaps S7;
research/signal_gaps/04_short_range.md section 3). Runs on the AD structures that core.remote_id.parse_ble_ad already
decodes (manufacturer data, 16-bit service data); each recognizer has strict length / constant checks so noise is not
labelled. Nothing here recovers a stable identity: the trackers broadcast rotating anonymous tokens. Byte parsing
only (no numpy); mirrored in C++ core/ble_trackers.

    recognize(ad) -> dict | None     {vendor, kind, label, ...fields}
"""
from __future__ import annotations

APPLE, MICROSOFT, SAMSUNG, GOOGLE = 0x004C, 0x0006, 0x0075, 0x00E0
BATT = {0: "critical", 1: "low", 2: "medium", 3: "full"}
APPLE_MODELS = {0x0E20: "AirPods Pro", 0x0220: "AirPods", 0x0F20: "AirPods Pro (2nd gen)", 0x0A20: "AirPods Max",
                0x0E30: "AirPods (3rd gen)", 0x1420: "AirPods Pro (2nd gen, USB-C)"}


def _mfg(ad: dict) -> tuple[int, bytes] | None:
    """(company id, the bytes after it) from manufacturer data, or None."""
    h = ad.get("manufacturer_data_hex")
    cid = ad.get("company_id")
    if h is None or cid is None:
        return None
    try:
        return int(cid, 16), bytes.fromhex(h)
    except ValueError:
        return None


def _svc(ad: dict, uuid: int) -> bytes | None:
    for s in ad.get("service_data", []):
        if int(s["uuid16"], 16) == uuid:
            return bytes.fromhex(s["data_hex"])
    return None


def _apple(cid: int, m: bytes) -> dict | None:
    """Apple Continuity TLVs (type, length, payload). Find My (0x12), iBeacon (0x02), proximity pairing (0x07)."""
    i = 0
    tlvs = []
    while i + 2 <= len(m):
        t, ln = m[i], m[i + 1]
        if i + 2 + ln > len(m):
            break
        tlvs.append((t, m[i + 2:i + 2 + ln]))
        i += 2 + ln
    for t, v in tlvs:
        if t == 0x02 and len(v) == 0x15:                      # iBeacon: 16-byte UUID, major, minor, tx power
            return {"vendor": "Apple", "kind": "iBeacon", "label": "iBeacon",
                    "proximity_uuid": v[:16].hex(), "major": (v[16] << 8) | v[17],
                    "minor": (v[18] << 8) | v[19], "measured_power_dbm": v[20] - 256 if v[20] >= 128 else v[20]}
    for t, v in tlvs:
        if t == 0x12 and len(v) >= 1:                          # Find My / AirTag
            status = v[0]
            dev = {0: "Apple device (Find My)", 1: "AirTag", 3: "AirPods"}.get((status & 0x30) >> 4, "Apple device")
            o = {"vendor": "Apple", "kind": "Find My", "label": f"Apple Find My / {dev}", "device_type": dev,
                 "offline_finding": len(v) >= 0x18, "maintained": bool(status & 0x04)}
            if status & 0x04:
                o["battery"] = BATT[(status >> 6) & 3]
            return o
        if t == 0x07 and len(v) >= 3:                          # proximity pairing (AirPods and similar)
            model = (v[1] << 8) | v[2]
            return {"vendor": "Apple", "kind": "proximity pairing",
                    "label": f"Apple {APPLE_MODELS.get(model, 'device')}".rstrip(),
                    "model_id": f"0x{model:04X}", "model": APPLE_MODELS.get(model)}
    if tlvs:
        return {"vendor": "Apple", "kind": "Continuity", "label": "Apple device (Continuity)",
                "continuity_types": [f"0x{t:02X}" for t, _ in tlvs]}
    return None


# kinds that are item trackers (finding networks); the rest are beacons / accessories / phones
TRACKER_KINDS = {"Find My", "SmartTag", "Find My Device (FMDN)", "Tile", "Chipolo", "tracker"}


def recognize(ad: dict) -> dict | None:
    """The strongest tracker / device match for one advert, or None. Adds tracker=True for item trackers."""
    r = _recognize(ad)
    if r is not None:
        r["tracker"] = r["kind"] in TRACKER_KINDS
    return r


def _recognize(ad: dict) -> dict | None:
    mf = _mfg(ad)
    if mf:
        cid, m = mf
        if cid == APPLE and m:
            r = _apple(cid, m)
            if r:
                return r
        if cid == MICROSOFT and len(m) >= 3 and m[2] == 0x80 and m[1] in (0, 1, 2):
            scen = {0: "LE", 1: "BR/EDR", 2: "LE + BR/EDR"}[m[1]]
            name = m[3:].decode("utf-8", "replace") if len(m) > 3 else None
            return {"vendor": "Microsoft", "kind": "Swift Pair", "label": "Microsoft Swift Pair",
                    "pairing": scen, "display_name": name or None}
    # Samsung SmartTag: service data under 0xFD5A (registered) / 0xFD59
    for uuid, reg in ((0xFD5A, True), (0xFD59, False)):
        d = _svc(ad, uuid)
        if d is not None and len(d) >= 20:
            state = (d[0] >> 5) & 7
            if 1 <= state <= 6:
                states = {1: "premature-lost", 2: "lost", 3: "overmature-lost", 4: "paired", 5: "connected",
                          6: "connected (two)"}
                return {"vendor": "Samsung", "kind": "SmartTag", "label": "Samsung SmartTag",
                        "registered": reg, "state": states[state], "battery": BATT[(d[12] >> 6) & 3],
                        "uwb": bool(d[12] & 0x20)}
    d = _svc(ad, 0xFEAA)
    if d is not None and len(d) >= 1:
        fr = d[0]
        # Google Find My Device Network: frame (1) + ephemeral ID (20, or 32 on the 256-bit curve) + optional
        # hashed flags (1), after the UUID
        if fr in (0x40, 0x41) and len(d) in (21, 22, 33, 34):
            return {"vendor": "Google", "kind": "Find My Device (FMDN)", "label": "Google Find My (FMDN)",
                    "unwanted_tracking_protection": fr == 0x41, "eid_bytes": 20 if len(d) < 33 else 32}
        names = {0x00: "UID", 0x10: "URL", 0x20: "TLM", 0x30: "EID"}
        if fr in names:                                           # Eddystone beacon
            return {"vendor": "Eddystone", "kind": names[fr], "label": f"Eddystone-{names[fr]}"}
    d = _svc(ad, 0xFE2C)
    if d is not None:                                             # Google Fast Pair
        if len(d) == 3:
            return {"vendor": "Google", "kind": "Fast Pair", "label": "Google Fast Pair",
                    "model_id": d.hex().upper()}
        if len(d) >= 1 and d[0] == 0x00:
            return {"vendor": "Google", "kind": "Fast Pair", "label": "Google Fast Pair (paired)"}
    for uuid, name in ((0xFEED, "Tile"), (0xFEEC, "Tile"), (0xFD84, "Tile"), (0xFE33, "Chipolo"),
                       (0xFD69, "Samsung Find My Mobile")):
        d = _svc(ad, uuid)
        if d is not None:
            return {"vendor": name.split()[0], "kind": name, "label": f"{name}-like"}   # service UUID only
    d = _svc(ad, 0xFCB2)                                          # IETF DULT (cross-vendor tracker)
    if d is not None and len(d) >= 2 and (d[1] & 0xFE) == 0:
        near = bool(d[1] & 0x01)
        return {"vendor": "DULT", "kind": "tracker", "label": "Tracker (DULT, near-owner)" if near
                else "Tracker (DULT, separated)", "near_owner": near}
    return None
