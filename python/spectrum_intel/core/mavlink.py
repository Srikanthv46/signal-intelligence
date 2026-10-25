"""MAVLink v1 / v2 frame scanning and message decoding (research: research/drone_transmitters/05_sik_mavlink_details.md).

Any byte stream (a SiK radio packet, an mLRS / LoRa telemetry payload, a Wi-Fi UDP payload) is scanned for frames:
v1 `FE len seq sys comp msgid payload crc`, v2 `FD len incompat compat seq sys comp msgid(3, LE) payload crc [sig 13]`.
The checksum is CRC-16/MCRF4XX (X.25: poly 0x1021 reflected, init 0xFFFF) over everything after the start byte up to
the end of the payload, then the message's CRC_EXTRA byte (from the MAVLink c_library_v2 generated headers). A frame
whose message id has a known CRC_EXTRA is reported only when its CRC passes (crc_ok); unknown ids are listed by id only.
v2 payloads are zero-trimmed on air: they are zero-filled to the message's full length before decoding.
Pure stdlib.
"""
from __future__ import annotations

import struct

# id: (name, CRC_EXTRA, full length (v2 max), [(field, struct fmt, offset)])
_F = struct.calcsize
MSGS = {
    0: ("HEARTBEAT", 50, 9, [("custom_mode", "<I", 0), ("type", "B", 4), ("autopilot", "B", 5), ("base_mode", "B", 6),
                              ("system_status", "B", 7), ("mavlink_version", "B", 8)]),
    1: ("SYS_STATUS", 124, 43, [("load", "<H", 12), ("voltage_battery", "<H", 14), ("current_battery", "<h", 16),
                                ("drop_rate_comm", "<H", 18), ("errors_comm", "<H", 20), ("battery_remaining", "b", 30)]),
    2: ("SYSTEM_TIME", 137, 12, [("time_unix_usec", "<Q", 0), ("time_boot_ms", "<I", 8)]),
    22: ("PARAM_VALUE", 220, 25, [("param_value", "<f", 0), ("param_count", "<H", 4), ("param_index", "<H", 6),
                                  ("param_id", "16s", 8), ("param_type", "B", 24)]),
    24: ("GPS_RAW_INT", 24, 52, [("time_usec", "<Q", 0), ("lat", "<i", 8), ("lon", "<i", 12), ("alt", "<i", 16),
                                 ("eph", "<H", 20), ("epv", "<H", 22), ("vel", "<H", 24), ("cog", "<H", 26),
                                 ("fix_type", "B", 28), ("satellites_visible", "B", 29)]),
    27: ("RAW_IMU", 144, 29, [("time_usec", "<Q", 0)]),
    29: ("SCALED_PRESSURE", 115, 16, [("time_boot_ms", "<I", 0), ("press_abs", "<f", 4), ("temperature", "<h", 12)]),
    30: ("ATTITUDE", 39, 28, [("time_boot_ms", "<I", 0), ("roll", "<f", 4), ("pitch", "<f", 8), ("yaw", "<f", 12)]),
    33: ("GLOBAL_POSITION_INT", 104, 28, [("time_boot_ms", "<I", 0), ("lat", "<i", 4), ("lon", "<i", 8), ("alt", "<i", 12),
                                          ("relative_alt", "<i", 16), ("vx", "<h", 20), ("vy", "<h", 22), ("vz", "<h", 24),
                                          ("hdg", "<H", 26)]),
    35: ("RC_CHANNELS_RAW", 244, 22, [("time_boot_ms", "<I", 0), ("rssi", "B", 21)]),
    36: ("SERVO_OUTPUT_RAW", 222, 37, [("time_usec", "<I", 0)]),
    42: ("MISSION_CURRENT", 28, 18, [("seq", "<H", 0)]),
    62: ("NAV_CONTROLLER_OUTPUT", 183, 26, [("nav_bearing", "<h", 20), ("target_bearing", "<h", 22), ("wp_dist", "<H", 24)]),
    65: ("RC_CHANNELS", 118, 42, [("time_boot_ms", "<I", 0), ("chancount", "B", 40), ("rssi", "B", 41)]),
    74: ("VFR_HUD", 20, 20, [("airspeed", "<f", 0), ("groundspeed", "<f", 4), ("alt", "<f", 8), ("climb", "<f", 12),
                            ("heading", "<h", 16), ("throttle", "<H", 18)]),
    76: ("COMMAND_LONG", 152, 33, [("command", "<H", 28), ("target_system", "B", 30), ("target_component", "B", 31)]),
    77: ("COMMAND_ACK", 143, 10, [("command", "<H", 0), ("result", "B", 2)]),
    109: ("RADIO_STATUS", 185, 9, [("rxerrors", "<H", 0), ("fixed", "<H", 2), ("rssi", "B", 4), ("remrssi", "B", 5),
                                   ("txbuf", "B", 6), ("noise", "B", 7), ("remnoise", "B", 8)]),
    111: ("TIMESYNC", 34, 18, [("tc1", "<q", 0), ("ts1", "<q", 8)]),
    147: ("BATTERY_STATUS", 154, 54, [("current_consumed", "<i", 0), ("energy_consumed", "<i", 4), ("temperature", "<h", 8),
                                      ("voltage_cell1", "<H", 10), ("current_battery", "<h", 30), ("battery_remaining", "b", 35)]),
    148: ("AUTOPILOT_VERSION", 178, 78, [("capabilities", "<Q", 0), ("uid", "<Q", 8), ("flight_sw_version", "<I", 16),
                                         ("board_version", "<I", 28), ("vendor_id", "<H", 32), ("product_id", "<H", 34)]),
    166: ("RADIO", 21, 9, [("rxerrors", "<H", 0), ("fixed", "<H", 2), ("rssi", "B", 4), ("remrssi", "B", 5),
                           ("txbuf", "B", 6), ("noise", "B", 7), ("remnoise", "B", 8)]),
    242: ("HOME_POSITION", 104, 60, [("latitude", "<i", 0), ("longitude", "<i", 4), ("altitude", "<i", 8)]),
    245: ("EXTENDED_SYS_STATE", 130, 2, [("vtol_state", "B", 0), ("landed_state", "B", 1)]),
    253: ("STATUSTEXT", 83, 54, [("severity", "B", 0), ("text", "50s", 1)]),
}
MAV_TYPE = {0: "generic", 1: "fixed wing", 2: "quadrotor", 3: "coaxial", 4: "helicopter", 5: "antenna tracker", 6: "GCS",
            7: "airship", 8: "free balloon", 9: "rocket", 10: "ground rover", 11: "surface boat", 12: "submarine",
            13: "hexarotor", 14: "octorotor", 15: "tricopter", 16: "flapping wing", 17: "kite", 18: "onboard controller",
            19: "VTOL tailsitter duorotor", 20: "VTOL tailsitter quadrotor", 21: "VTOL tiltrotor", 22: "VTOL fixedrotor",
            23: "VTOL tailsitter", 24: "VTOL tiltwing", 26: "gimbal", 27: "ADS-B", 28: "parafoil", 29: "dodecarotor",
            30: "camera", 32: "FLARM", 34: "ODID", 35: "decarotor"}
MAV_AUTOPILOT = {0: "generic", 3: "ArduPilot", 4: "OpenPilot", 8: "invalid (GCS / peripheral)", 9: "Paparazzi", 12: "PX4",
                 13: "SMACCMPilot", 14: "AutoQuad", 18: "SmartAP", 20: "Reflex", 21: "Flix"}
MAV_STATE = {0: "uninit", 1: "boot", 2: "calibrating", 3: "standby", 4: "active", 5: "critical", 6: "emergency",
             7: "poweroff", 8: "flight termination"}
COPTER = {0: "STABILIZE", 1: "ACRO", 2: "ALT_HOLD", 3: "AUTO", 4: "GUIDED", 5: "LOITER", 6: "RTL", 7: "CIRCLE", 9: "LAND",
          11: "DRIFT", 13: "SPORT", 14: "FLIP", 15: "AUTOTUNE", 16: "POSHOLD", 17: "BRAKE", 18: "THROW", 19: "AVOID_ADSB",
          20: "GUIDED_NOGPS", 21: "SMART_RTL", 22: "FLOWHOLD", 23: "FOLLOW", 24: "ZIGZAG", 25: "SYSTEMID",
          26: "AUTOROTATE", 27: "AUTO_RTL", 28: "TURTLE"}
PLANE = {0: "MANUAL", 1: "CIRCLE", 2: "STABILIZE", 3: "TRAINING", 4: "ACRO", 5: "FBWA", 6: "FBWB", 7: "CRUISE",
         8: "AUTOTUNE", 10: "AUTO", 11: "RTL", 12: "LOITER", 13: "TAKEOFF", 14: "AVOID_ADSB", 15: "GUIDED",
         16: "INITIALIZING", 17: "QSTABILIZE", 18: "QHOVER", 19: "QLOITER", 20: "QLAND", 21: "QRTL", 22: "QAUTOTUNE",
         23: "QACRO", 24: "THERMAL", 25: "LOITER_ALT_QLAND", 26: "AUTOLAND"}
ROVER = {0: "MANUAL", 1: "ACRO", 3: "STEERING", 4: "HOLD", 5: "LOITER", 6: "FOLLOW", 7: "SIMPLE", 8: "DOCK", 9: "CIRCLE",
         10: "AUTO", 11: "RTL", 12: "SMART_RTL", 15: "GUIDED", 16: "INITIALIZING"}
PX4_MAIN = {1: "MANUAL", 2: "ALTCTL", 3: "POSCTL", 4: "AUTO", 5: "ACRO", 6: "OFFBOARD", 7: "STABILIZED",
            8: "RATTITUDE_LEGACY", 9: "SIMPLE", 10: "TERMINATION", 11: "ALTITUDE_CRUISE", 12: "MANUAL_PARKING"}
PX4_AUTO = {1: "READY", 2: "TAKEOFF", 3: "LOITER", 4: "MISSION", 5: "RTL", 6: "LAND", 8: "FOLLOW_TARGET", 9: "PRECLAND",
            10: "VTOL_TAKEOFF", 19: "GUIDED_COURSE", 20: "AUTO_DESCEND"}


def x25(data: bytes, crc: int = 0xFFFF) -> int:
    for b in data:
        t = (b ^ (crc & 0xFF)) & 0xFF
        t = (t ^ (t << 4)) & 0xFF
        crc = ((crc >> 8) ^ (t << 8) ^ (t << 3) ^ (t >> 4)) & 0xFFFF
    return crc


def flight_mode(autopilot: int, mav_type: int, custom_mode: int) -> str | None:
    if autopilot == 3:
        tab = PLANE if mav_type in (1, 16, 19, 20, 21, 22, 23, 24) else ROVER if mav_type in (10, 11) else COPTER
        return tab.get(custom_mode, f"mode {custom_mode}")
    if autopilot == 12:
        main, sub = (custom_mode >> 16) & 0xFF, (custom_mode >> 24) & 0xFF
        name = PX4_MAIN.get(main, f"main {main}")
        return f"{name}/{PX4_AUTO.get(sub, sub)}" if main == 4 else name
    return None


def _decode(msgid: int, payload: bytes) -> dict:
    name, _, full, fields = MSGS[msgid]
    p = payload + bytes(max(0, full - len(payload)))          # v2 zero-trimmed payloads
    out = {"msg": name}
    for f, fmt, off in fields:
        if off + _F(fmt) > len(p):
            continue
        v = struct.unpack_from(fmt, p, off)[0]
        if isinstance(v, bytes):
            v = v.split(b"\x00", 1)[0].decode("utf-8", "replace")
        elif isinstance(v, float):
            v = round(v, 6)
        out[f] = v
    if name == "HEARTBEAT":
        out["vehicle_type"] = MAV_TYPE.get(out["type"], f"type {out['type']}")
        out["autopilot_name"] = MAV_AUTOPILOT.get(out["autopilot"], f"autopilot {out['autopilot']}")
        out["armed"] = bool(out["base_mode"] & 0x80)
        out["state"] = MAV_STATE.get(out["system_status"], out["system_status"])
        fm = flight_mode(out["autopilot"], out["type"], out["custom_mode"]) if out["base_mode"] & 0x01 else None
        if fm:
            out["flight_mode"] = fm
    for k in ("lat", "lon", "latitude", "longitude"):
        if k in out:
            out[k + "_deg"] = round(out[k] * 1e-7, 7)
    return out


def scan(data: bytes) -> list[dict]:
    """All MAVLink frames in `data` (frames may start anywhere). Each: version, seq, sysid, compid, msgid, len,
    crc_ok (None when the id's CRC_EXTRA is unknown), and the decoded fields when crc_ok."""
    b = bytes(data)
    out, i = [], 0
    while i < len(b):
        stx = b[i]
        if stx not in (0xFE, 0xFD) or i + 2 > len(b):
            i += 1
            continue
        ln = b[i + 1]
        if stx == 0xFE:
            hdr, end = 6, i + 6 + ln + 2
            if end > len(b):
                i += 1
                continue
            seq, sysid, compid, msgid = b[i + 2], b[i + 3], b[i + 4], b[i + 5]
            sig = 0
        else:
            hdr = 10
            if i + hdr > len(b):
                i += 1
                continue
            incompat = b[i + 2]
            if incompat & ~0x01:
                i += 1
                continue
            sig = 13 if incompat & 0x01 else 0
            end = i + hdr + ln + 2
            if end + sig > len(b):
                i += 1
                continue
            seq, sysid, compid = b[i + 4], b[i + 5], b[i + 6]
            msgid = b[i + 7] | (b[i + 8] << 8) | (b[i + 9] << 16)
        payload = b[i + hdr:i + hdr + ln]
        crc_rx = b[end - 2] | (b[end - 1] << 8)
        info = MSGS.get(msgid)
        if info is not None:
            ok = x25(bytes([info[1]]), x25(b[i + 1:i + hdr + ln])) == crc_rx
        else:
            ok = None
        fr = dict(version=1 if stx == 0xFE else 2, seq=seq, sysid=sysid, compid=compid, msgid=msgid, len=ln, crc_ok=ok)
        if ok:
            fr.update(_decode(msgid, payload))
            out.append(fr)
            i = end + sig
            continue
        if ok is None and stx == 0xFD:                           # unknown id: keep only well-formed v2 frames
            out.append(fr)
            i = end + sig
            continue
        i += 1
    return out


def summary(frames: list[dict]) -> dict:
    """Latest value of the identity / position / state fields over CRC-valid frames (vehicle sysid first seen)."""
    s: dict = {}
    for f in frames:
        if not f.get("crc_ok"):
            continue
        s.setdefault("sysid", f["sysid"])
        s.setdefault("messages", [])
        if f["msg"] not in s["messages"]:
            s["messages"].append(f["msg"])
        m = f["msg"]
        if m == "HEARTBEAT" and f["sysid"] != 255 and f["type"] != 6:
            s.update(vehicle_type=f["vehicle_type"], autopilot=f["autopilot_name"], armed=f["armed"], state=f["state"])
            if "flight_mode" in f:
                s["flight_mode"] = f["flight_mode"]
        elif m == "GLOBAL_POSITION_INT":
            s.update(lat=f["lat_deg"], lon=f["lon_deg"], alt_msl_m=f["alt"] / 1000, alt_rel_m=f["relative_alt"] / 1000,
                     heading_deg=None if f["hdg"] == 65535 else f["hdg"] / 100)
        elif m == "GPS_RAW_INT":
            s.update(gps_fix_type=f["fix_type"], satellites=f["satellites_visible"])
            s.setdefault("lat", f["lat_deg"])
            s.setdefault("lon", f["lon_deg"])
        elif m == "HOME_POSITION":
            s.update(home_lat=f["latitude_deg"], home_lon=f["longitude_deg"], home_alt_m=f["altitude"] / 1000)
        elif m == "VFR_HUD":
            s.update(groundspeed_mps=f["groundspeed"], airspeed_mps=f["airspeed"], climb_mps=f["climb"], throttle_pct=f["throttle"])
        elif m == "SYS_STATUS":
            s.update(battery_v=f["voltage_battery"] / 1000)
            if f.get("battery_remaining", -1) >= 0:
                s["battery_pct"] = f["battery_remaining"]
        elif m == "ATTITUDE":
            s.update(roll_rad=f["roll"], pitch_rad=f["pitch"], yaw_rad=f["yaw"])
        elif m == "STATUSTEXT":
            s.setdefault("statustext", []).append(f["text"])
        elif m in ("RADIO", "RADIO_STATUS"):
            s.update(radio_rssi=f["rssi"], radio_remrssi=f["remrssi"], radio_noise=f["noise"])
    return s
