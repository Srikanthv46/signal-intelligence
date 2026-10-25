"""Per-aircraft summary from decoded Mode S frames (messages.fields): identity, CPR positions,
kinematics, intent, integrity, capabilities, Comm-B air data, derived wind / temperature, the
interrogation environment and plausibility flags.

    tracks(msgs, ref=None) -> list[dict]
      msgs: frames with 't' (absolute or capture-relative seconds), from messages.fields
      ref:  (lat, lon[, alt_m]) of the receiver, optional; adds range/bearing/elevation and allows
            single-frame (local) CPR decoding
"""
from __future__ import annotations

import math
from collections import Counter

from .tables import NACP_EPU_M, OPERATORS, country

PAIR_S = 10.0               # max even/odd spacing for a global CPR decode
LOCAL_S = 30.0              # a fix is a valid local-decode reference for this long
MAX_KT = 1000.0             # faster than this between fixes -> implausible position


def _haversine_km(a, b):
    la1, lo1, la2, lo2 = map(math.radians, (*a, *b))
    h = math.sin((la2 - la1) / 2) ** 2 + math.cos(la1) * math.cos(la2) * math.sin((lo2 - lo1) / 2) ** 2
    return 2 * 6371.0 * math.asin(math.sqrt(h))


def _bearing(a, b):
    la1, lo1, la2, lo2 = map(math.radians, (*a, *b))
    y = math.sin(lo2 - lo1) * math.cos(la2)
    x = math.cos(la1) * math.sin(la2) - math.sin(la1) * math.cos(la2) * math.cos(lo2 - lo1)
    return (math.degrees(math.atan2(y, x)) + 360) % 360


def _positions(pms, frames, ref):
    """Global CPR from an even/odd pair, then local decoding against the last fix (or the receiver)."""
    fixes, last, cpr = [], None, {}
    for m in frames:
        tc = m.get("typecode") or 0
        if "lat_cpr" not in m or not (9 <= tc <= 18 or 20 <= tc <= 22 or (5 <= tc <= 8 and ref)):
            continue
        t, h = m["t"], m["hex"]
        cpr[m["odd"]] = (h, t)
        pos = method = None
        if last and t - last[0] <= LOCAL_S:
            fn = pms.adsb.surface_position_with_ref if tc <= 8 else pms.adsb.airborne_position_with_ref
            pos, method = fn(h, last[1], last[2]), "local (previous fix)"
        elif 0 in cpr and 1 in cpr and abs(cpr[0][1] - cpr[1][1]) <= PAIR_S and tc >= 9:
            pos, method = pms.adsb.airborne_position(cpr[0][0], cpr[1][0], cpr[0][1], cpr[1][1]), "global (even/odd)"
        elif ref:
            fn = pms.adsb.surface_position_with_ref if tc <= 8 else pms.adsb.airborne_position_with_ref
            pos, method = fn(h, ref[0], ref[1]), "local (receiver)"
        if not pos:
            continue
        fix = {"t": round(t, 3), "lat": round(pos[0], 5), "lon": round(pos[1], 5), "method": method,
               "alt_ft": m.get("altitude_ft", m.get("gnss_altitude_ft"))}
        if last and t > last[0]:
            kt = _haversine_km(last[1:], pos) / 1.852 / ((t - last[0]) / 3600)
            if kt > MAX_KT:
                fix["implausible_speed_kt"] = round(kt)
                fixes.append(fix)
                continue                                 # do not use it as the next reference
        fixes.append(fix)
        last = (t, pos[0], pos[1])
    return fixes


def _nearest(items, t, max_dt):
    best = min(items, key=lambda x: abs(x[0] - t), default=None)
    return best[1] if best and abs(best[0] - t) <= max_dt else None


def _resolve_commb(pms, frames, vel, alts):
    """Pick one register for ambiguous DF20/21 replies (typically BDS50 vs BDS60) using the
    aircraft's own ADS-B ground speed, track and altitude as the reference (pyModeS is50or60)."""
    out = []
    for m in frames:
        cb = m.get("commb")
        if not cb or not cb["candidates"]:
            continue
        regs = list(cb["candidates"])
        reg = regs[0] if len(regs) == 1 else None
        if reg is None and set(regs) == {"BDS50", "BDS60"}:
            v, a = _nearest(vel, m["t"], 10), _nearest(alts, m["t"], 10)
            if v and a is not None:
                reg = pms.bds.is50or60(m["hex"], v[0], v[1], a)
        m["commb"]["register"] = reg or "ambiguous"
        if reg:
            out.append((m["t"], reg, cb["candidates"][reg]))
    return out


def _met(commb, vel, alts):
    """Wind from (TAS, heading) vs (GS, track); static air temperature from TAS and Mach."""
    tas = [(t, d["tas_kt"]) for t, r, d in commb if r == "BDS50" and d.get("tas_kt")]
    hdg = [(t, d["magnetic_heading_deg"]) for t, r, d in commb if r == "BDS60" and d.get("magnetic_heading_deg") is not None]
    mach = [(t, d["mach"]) for t, r, d in commb if r == "BDS60" and d.get("mach")]
    met = {}
    for t, h in hdg:
        a, g = _nearest(tas, t, 5), _nearest(vel, t, 5)
        if a and g:
            ax, ay = a * math.sin(math.radians(h)), a * math.cos(math.radians(h))
            gx, gy = g[0] * math.sin(math.radians(g[1])), g[0] * math.cos(math.radians(g[1]))
            wx, wy = gx - ax, gy - ay
            met.update(wind_kt=round(math.hypot(wx, wy)), wind_from_deg=round((math.degrees(math.atan2(wx, wy)) + 180) % 360),
                       wind_method="ground vector (ADS-B GS/track) - air vector (BDS50 TAS, BDS60 magnetic heading; "
                                   "declination not applied, ~0-2 deg over India)")
    for t, mc in mach:
        a, alt = _nearest(tas, t, 5), _nearest(alts, t, 10)
        if a and mc > 0.1:
            tk = (a / (mc * 38.967854)) ** 2
            sat = tk - 273.15
            # half-LSB quantisation of TAS (2 kt) and Mach (0.004), doubled by the square
            err = tk * 2 * math.hypot(1 / a, 0.002 / mc)
            met.update(static_air_temp_c=round(sat, 1), temp_uncertainty_c=round(err, 1),
                       temp_method="SAT = (TAS / (Mach * 38.97))^2, BDS50 TAS + BDS60 Mach")
            if alt is not None and alt < 36089:
                met["isa_deviation_c"] = round(sat - (15 - 1.98 * alt / 1000), 1)
    for t, r, d in commb:                                # MRAR (BDS44): reported directly
        if r == "BDS44":
            if d.get("wind"):
                met.update(wind_kt=d["wind"][0], wind_from_deg=d["wind"][1], wind_method="BDS44 (reported)")
            if d.get("temperature_c") is not None:
                met.update(static_air_temp_c=d["temperature_c"], temp_method="BDS44 (reported)")
    return met


def _last(frames, key):
    for m in reversed(frames):
        if m.get(key) is not None:
            return m[key]
    return None


def _one(pms, icao, frames, ref):
    frames = sorted(frames, key=lambda m: m["t"])
    vel = [(m["t"], (m["ground_speed_kt"], m["track_deg"])) for m in frames
           if m.get("typecode") == 19 and m.get("ground_speed_kt") is not None]
    alts = [(m["t"], m["altitude_ft"]) for m in frames if m.get("altitude_ft") is not None]
    a = {"icao": icao, "country": country(icao), "first_t": round(frames[0]["t"], 3), "last_t": round(frames[-1]["t"], 3),
         "n_frames": len(frames), "dfs": dict(Counter(m["df"] for m in frames)),
         "typecodes": dict(Counter(m["typecode"] for m in frames if m.get("typecode")))}
    flags = []

    cs = [m["callsign"] for m in frames if m.get("callsign")]
    if cs:
        a["callsign"] = cs[-1]
        op = OPERATORS.get(cs[-1][:3])
        if op and cs[-1][3:].isdigit():
            a["operator"], a["flight"] = op[0], op[1] + cs[-1][3:]
        if len(set(cs)) > 1:
            flags.append(f"callsign changed: {sorted(set(cs))}")
    for k in ("category", "capability", "control_field", "adsb_version", "squawk", "emergency"):
        v = _last(frames, k)
        if v is not None:
            a[k] = v
    sq = {m["squawk"] for m in frames if m.get("squawk")}
    if len(sq) > 1:
        flags.append(f"squawk changed: {sorted(sq)}")
    if sq & {"7500", "7600", "7700"}:
        flags.append(f"emergency squawk {sorted(sq & {'7500', '7600', '7700'})}")
    if a.get("emergency") not in (None, "no emergency"):
        flags.append(f"emergency: {a['emergency']}")
    if a.get("category") == "unmanned aerial vehicle":       # emitter category B6 (TC 3, CA 6): a drone on ADS-B
        a["uav"] = True
        flags.append("UAV (ADS-B emitter category B6)")

    # kinematics
    if alts:
        al = [x for _, x in alts]
        a["altitude_ft"] = {"first": al[0], "last": al[-1], "min": min(al), "max": max(al)}
        if len(alts) > 1 and alts[-1][0] > alts[0][0]:
            a["altitude_rate_fpm"] = round((al[-1] - al[0]) / (alts[-1][0] - alts[0][0]) * 60)
        vr = a.get("altitude_rate_fpm", 0)
        a["vertical_trend"] = "climbing" if vr > 300 else "descending" if vr < -300 else "level"
    for k in ("gnss_altitude_ft", "gnss_minus_baro_ft", "ground_speed_kt", "track_deg", "vrate_fpm", "vrate_source",
              "ias_kt", "tas_kt", "heading_deg", "heading_source"):
        v = _last(frames, k)
        if v is not None:
            a[k] = v
    if vel:
        trk = [v[1] for _, v in vel]
        a["ground_speed_range_kt"] = [min(v[0] for _, v in vel), max(v[0] for _, v in vel)]
        turn = ((trk[-1] - trk[0] + 540) % 360) - 180
        if abs(turn) >= 5:
            a["turn_deg"] = round(turn, 1)               # + right, - left, over the observation
    alt_last = alts[-1][1] if alts else 0
    if abs(a.get("gnss_minus_baro_ft") or 0) > max(1500, 0.08 * alt_last):   # beyond what temperature explains
        flags.append(f"GNSS and baro altitude differ by {a['gnss_minus_baro_ft']} ft")

    # position
    fixes = _positions(pms, frames, ref)
    good = [f for f in fixes if "implausible_speed_kt" not in f]
    if good:
        a["position"] = good[-1]
        a["n_positions"] = len(good)
        a["track"] = [[f["t"], f["lat"], f["lon"], f["alt_ft"]] for f in good]
        if ref:
            p = (good[-1]["lat"], good[-1]["lon"])
            d = _haversine_km(ref[:2], p)
            a["range_km"], a["bearing_deg"] = round(d, 2), round(_bearing(ref[:2], p), 1)
            if good[-1]["alt_ft"] is not None:
                dz = good[-1]["alt_ft"] * 0.3048 - (ref[2] if len(ref) > 2 else 0)
                a["elevation_deg"] = round(math.degrees(math.atan2(dz, d * 1000)), 1)
    if len(fixes) > len(good):
        flags.append(f"{len(fixes) - len(good)} position(s) implausible (> {MAX_KT:.0f} kt jump): possible spoofing or bad CPR")

    # intent (TC29, BDS40)
    for k in ("selected_altitude_ft", "selected_altitude_source", "baro_setting_hpa", "selected_heading_deg", "modes",
              "tcas_operational"):
        v = _last(frames, k)
        if v is not None:
            a.setdefault("intent", {})[k] = v

    # integrity and capabilities (TC29 / TC31 / DF0)
    integ = {k: _last(frames, k) for k in ("nac_p", "nic_a", "nic_b", "nic_baro", "sil", "sil_per", "gva", "nac_v", "sda")}
    integ = {k: v for k, v in integ.items() if v is not None}
    if "nac_p" in integ:
        integ["position_accuracy_m"] = f"< {NACP_EPU_M.get(integ['nac_p'], 'unknown')}"
    if integ:
        a["integrity"] = integ
    cap = dict(_last(frames, "capabilities") or {})
    for k in ("single_antenna", "heading_reference", "acas_reply_info", "acas_sensitivity", "crosslink_capable"):
        v = _last(frames, k)
        if v is not None:
            cap[k] = v
    if cap:
        a["capabilities"] = cap

    # Comm-B (Mode S ground radar EHS / MRAR) and derived met
    commb = _resolve_commb(pms, frames, vel, alts)
    if commb:
        cb = {"registers": dict(Counter(r for _, r, _ in commb))}
        for _, r, d in commb:
            cb.update({k: v for k, v in d.items() if v is not None})
        a["commb"] = cb
        met = _met(commb, vel, alts)
        if met:
            a["met"] = met
        bdscs = [d["callsign"].strip("_ ") for _, r, d in commb if r == "BDS20" and d.get("callsign")]
        if cs and bdscs and bdscs[-1] != cs[-1]:
            flags.append(f"Comm-B callsign {bdscs[-1]} != ADS-B callsign {cs[-1]}")

    # who is interrogating this aircraft
    acas = sum(1 for m in frames if m["df"] in (0, 16))
    radar = [m for m in frames if m["df"] in (4, 5, 20, 21)]
    dwells = []                                          # radar beam dwells: replies within 0.2 s
    for m in radar:
        if not dwells or m["t"] - dwells[-1][-1] > 0.2:
            dwells.append([m["t"]])
        else:
            dwells[-1].append(m["t"])
    a["interrogation"] = {
        "acas_replies": acas, "mode_s_radar_replies": len(radar), "radar_dwells": len(dwells),
        "radar_iis": sorted({m["iis"] for m in radar if m.get("ids") == 1}),
        "df11_interrogators": dict(Counter(m["interrogator"] for m in frames if m.get("interrogator"))),
        "squitters_df17_18": sum(1 for m in frames if m["df"] in (17, 18)),
    }
    if len(dwells) >= 2:
        gaps = sorted(dwells[i + 1][0] - dwells[i][0] for i in range(len(dwells) - 1))
        a["interrogation"]["radar_revisit_s"] = round(gaps[len(gaps) // 2], 2)   # antenna rotation period (or a multiple)
    if any(m.get("acas_ra_report") or m.get("tcas_ra_active") for m in frames):
        flags.append("ACAS resolution advisory reported")

    snr = sorted(m["snr_db"] for m in frames if m.get("snr_db") is not None)
    if snr:
        a["signal"] = {"snr_db_min": snr[0], "snr_db_median": snr[len(snr) // 2], "snr_db_max": snr[-1]}
    span = frames[-1]["t"] - frames[0]["t"]
    if span > 1:
        a["frame_rate_hz"] = round(len(frames) / span, 2)
    if flags:
        a["flags"] = flags
    return a


def tracks(msgs, ref=None):
    try:
        import pyModeS as pms
    except Exception:                                    # noqa: BLE001
        return []
    by = {}
    for m in msgs:
        if m.get("icao"):
            by.setdefault(m["icao"], []).append(m)
    out = [_one(pms, k, v, ref) for k, v in by.items()]
    return sorted(out, key=lambda a: -a["n_frames"])
