"""Drone Wi-Fi links among the 802.11 traffic of one signal (research, PLAN W1; research file 18).

    frames = [dict(t0, t1, dur_us, bytes, ta, macs, ctrl, video), ...]   # one per PPDU (plugin adapters)
    links(frames, block_s)      -> one profile per BSS: AP, clients, airtime / PPDUs / bytes per direction, drone
                                   evidence and a drone-link verdict (identity, else behaviour)

Device tracks, alerts (deauthentication floods, duplicate SSIDs, beacon floods) and randomised-MAC linking are
report-level, in receivers/devices; this module works per signal because it also counts the PPDUs whose payload is
out of reach (attributed by links.attribute), which never reach the device tracker.

Drone over Wi-Fi, by behaviour (Alipour-Fanid et al., IEEE TIFS 2020; Alsoliman et al., ICDCN 2021): the drone is the
AP (or the controller is), it has one client (the pilot's phone or controller), it sends a steady heavy downlink
(video) and receives small frequent uplink frames (control). Encrypted (WPA2) or undecodable PPDUs still count: the
PHY gives every PPDU's airtime and links.attribute its transmitter. One 34 ms block shows the shape, not the video
frame rate; an IP camera or a phone hotspot streaming video looks the same, so a behaviour-only verdict stays
"-like" at 0.55 and says so.
"""
from __future__ import annotations

import numpy as np

CTRL_MAX_BYTES = 40              # ACK 14, CTS 14, RTS 20, BlockAckReq 24, BlockAck 32
BIG_BYTES = 600                  # a video / bulk data frame


def _group(a):
    return a is None or bool(int(a[:2], 16) & 1)


def _mac_list(f):
    return f.get("macs") or []


def _bss_of(m):
    b = m.get("bssid")
    return None if _group(b) else b


def _side():
    return {"ppdus": 0, "data_ppdus": 0, "airtime_us": 0.0, "data_airtime_us": 0.0, "bytes": 0, "sizes": [],
            "big_ppdus": 0}


def links(frames, block_s):
    """Profiles of the BSSs heard in these frames (see the module docstring)."""
    aps, clients, info = {}, {}, {}
    for f in frames:                                     # who is an AP, and who talks to it
        for m in _mac_list(f):
            b = _bss_of(m)
            if not b:
                continue
            aps.setdefault(b, {"bssid": b})
            if m.get("subtype") in ("beacon", "probe-resp"):
                for k in ("ssid", "security", "phy"):
                    if m.get(k) is not None:
                        aps[b].setdefault(k, m[k])
            # wireless peers only: the transmitter of a to-DS frame, the receiver of a from-DS frame (SA / DA can be
            # wired hosts behind the AP), either end of an association / authentication exchange
            if m.get("type") == "data" and not (m.get("to_ds") and m.get("from_ds")):
                peers = [m.get("ta")] if m.get("to_ds") else [m.get("ra")] if m.get("from_ds") else \
                    [m.get("ta"), m.get("ra")]
            elif m.get("subtype") in ("assoc-req", "reassoc-req", "assoc-resp", "reassoc-resp", "auth"):
                peers = [m.get("ta"), m.get("ra")]
            else:
                peers = []
            for a in peers:
                if a and a != b and not _group(a):
                    clients.setdefault(b, set()).add(a)
    if not aps:
        return []
    owner = {}                                           # station -> its BSS (only when it belongs to one)
    for b, cs in clients.items():
        for c in cs:
            owner[c] = None if c in owner else b
    for b in aps:
        info[b] = {"down": _side(), "up": _side(), "beacons": 0, "drone": [], "video": 0, "t": []}
    for f in frames:
        ta = f.get("ta")
        b = side = None
        if ta in aps:
            b, side = ta, "down"
        elif ta and owner.get(ta):
            b, side = owner[ta], "up"
        if b is None:
            continue
        r = info[b]
        ms = _mac_list(f)
        if any(m.get("subtype") == "beacon" for m in ms):
            r["beacons"] += 1
            continue
        s = r[side]
        s["ppdus"] += 1
        s["airtime_us"] += f["dur_us"]
        ctrl = f.get("ctrl")
        if not ctrl:
            s["data_ppdus"] += 1
            s["data_airtime_us"] += f["dur_us"]
            if f.get("bytes") is not None:
                s["bytes"] += f["bytes"]
                s["sizes"].append(f["bytes"])
                s["big_ppdus"] += f["bytes"] >= BIG_BYTES
            if side == "down":
                r["t"].append(f["t0"])
        r["video"] += bool(f.get("video"))
        for m in ms:
            if m.get("drone"):
                r["drone"].append(m["drone"])
    out = []
    for b, ap in aps.items():
        r = info[b]
        if not (r["down"]["ppdus"] or r["up"]["ppdus"]):
            continue
        cs = sorted(clients.get(b, ()))
        row = {**ap, "clients": cs, "n_clients": len(cs), "beacons": r["beacons"]}
        for side in ("down", "up"):
            s = r[side]
            row[side + "link"] = {"ppdus": s["ppdus"], "data_ppdus": s["data_ppdus"],
                                  "airtime_ms": round(s["airtime_us"] / 1e3, 3),
                                  "data_airtime_ms": round(s["data_airtime_us"] / 1e3, 3),
                                  "bytes": s["bytes"], "big_ppdus": s["big_ppdus"],
                                  "median_bytes": int(np.median(s["sizes"])) if s["sizes"] else None}
        row["data_airtime_pct"] = round(100 * (r["down"]["data_airtime_us"] + r["up"]["data_airtime_us"])
                                        / (block_s * 1e6), 2)
        row["downlink_mbps"] = round(r["down"]["bytes"] * 8 / block_s / 1e6, 3)
        if len(r["t"]) >= 3:
            g = np.diff(sorted(r["t"])) * 1e3
            row["downlink_gap_ms"] = {"median": round(float(np.median(g)), 3), "max": round(float(g.max()), 3)}
        row["video_frames"] = r["video"]
        v = verdict(row, r["drone"], block_s)
        if v:
            row["drone_link"] = v
        out.append(row)
    return sorted(out, key=lambda x: -(x["downlink"]["airtime_ms"] + x["uplink"]["airtime_ms"]))


def verdict(row, drone, block_s):
    """Identity (Remote ID, Tello SDK, MAVLink, wfb-ng, ESP-NOW, SSID / OUI signature), else behaviour."""
    d, u = row["downlink"], row["uplink"]
    beh = {"one_client": row["n_clients"] == 1,
           "downlink_dominant": d["data_ppdus"] >= 5 and d["data_airtime_ms"] >= 3 * u["data_airtime_ms"],
           "sustained_downlink": d["data_airtime_ms"] >= 0.05 * block_s * 1e3 or row["downlink_mbps"] >= 1.0,
           "uplink_small": u["median_bytes"] is None or u["median_bytes"] < 200,
           "video_payload": row["video_frames"] > 0}
    if drone:
        best = max(drone, key=lambda x: x.get("confidence", 0))
        return {"verdict": "drone link", "platform": best["platform"], "basis": "identity",
                "source": best.get("source"), "confidence": best.get("confidence"), "evidence": best.get("evidence"),
                "behaviour": beh}
    if beh["one_client"] and beh["downlink_dominant"] and beh["sustained_downlink"] and beh["uplink_small"]:
        return {"verdict": "drone / FPV video link-like", "platform": None, "basis": "behaviour",
                "confidence": 0.65 if beh["video_payload"] else 0.55, "behaviour": beh,
                "caveat": "one block shows a single client receiving a steady heavy downlink; an IP camera or a "
                          "phone hotspot streaming video looks the same: confirm over ~1 s of blocks or by identity"}
    return None
