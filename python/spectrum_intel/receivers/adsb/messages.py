"""Per-frame Mode S / ADS-B field decoding (pyModeS for the encodings, DO-260B / ICAO Annex 10
bit layouts for the fields pyModeS leaves out).

    fields(msgs) -> msgs      adds the decoded fields to each frame dict from modes.demod

Fields that need more than one frame (CPR position, resolving ambiguous Comm-B registers, wind,
temperature) are left to tracks.tracks().
"""
from __future__ import annotations

from .tables import (ACAS_RI, ACAS_SL, CAPABILITY, CATEGORY, CF, EMERGENCY, FLIGHT_STATUS,
                     SURVEILLANCE_STATUS)

# Comm-B register -> {field: pyModeS.commb function}
COMMB = {
    "BDS20": {"callsign": "cs20"},
    "BDS40": {"selected_alt_mcp_ft": "selalt40mcp", "selected_alt_fms_ft": "selalt40fms",
              "baro_setting_hpa": "p40baro"},
    "BDS44": {"wind": "wind44", "temperature_c": "temp44", "pressure_hpa": "p44", "humidity_pct": "hum44",
              "turbulence": "turb44"},
    "BDS45": {"turbulence": "turb45", "wind_shear": "ws45", "microburst": "mb45", "icing": "ic45",
              "wake_vortex": "wv45", "temperature_c": "temp45", "pressure_hpa": "p45", "radio_height_ft": "rh45"},
    "BDS50": {"roll_deg": "roll50", "true_track_deg": "trk50", "ground_speed_kt": "gs50",
              "track_rate_dps": "rtrk50", "tas_kt": "tas50"},
    "BDS60": {"magnetic_heading_deg": "hdg60", "ias_kt": "ias60", "mach": "mach60",
              "baro_vrate_fpm": "vr60baro", "inertial_vrate_fpm": "vr60ins"},
}


def _bits(h):
    return bin(int(h, 16))[2:].zfill(len(h) * 4)


def _u(b, i, j):
    return int(b[i:j], 2)


def _r(v, n=3):
    return round(v, n) if isinstance(v, float) else v


def _commb(pms, h):
    """Infer the Comm-B register(s) of a DF20/21 MB field and decode each candidate."""
    inferred = pms.bds.infer(h, mrar=True)
    out = {"bds": inferred, "candidates": {}}
    for reg in (inferred or "").split(","):
        if reg == "BDS17":
            out["candidates"][reg] = {"supported_registers": pms.commb.cap17(h)}
        elif reg in COMMB:
            dec = {}
            for k, fn in COMMB[reg].items():
                try:
                    v = getattr(pms.commb, fn)(h)
                except Exception:                        # noqa: BLE001 — field not present in this frame
                    v = None
                if v is not None:
                    dec[k] = [_r(x) for x in v] if isinstance(v, tuple) else _r(v)
            out["candidates"][reg] = dec
        elif reg:
            out["candidates"][reg] = {}
    return out


def _me_position(b, m, tc):
    me = b[32:88]
    m["odd"] = int(me[21])
    m["lat_cpr"], m["lon_cpr"] = _u(me, 22, 39), _u(me, 39, 56)
    m["time_sync"] = int(me[20])
    if tc >= 9:
        m["surveillance_status"] = SURVEILLANCE_STATUS[_u(me, 5, 7)]
        m["nic_b"] = int(me[7])


def _me_target_state(b, m):
    """TC29 subtype 1 (ADS-B v2), DO-260B 2.2.3.2.7.1."""
    me = b[32:88]
    m["tss_subtype"] = _u(me, 5, 7)
    if m["tss_subtype"] != 1:
        return
    alt = _u(me, 9, 20)
    if alt:
        m["selected_altitude_ft"] = (alt - 1) * 32
        m["selected_altitude_source"] = "FMS" if me[8] == "1" else "MCP/FCU"
    baro = _u(me, 20, 29)
    if baro:
        m["baro_setting_hpa"] = round((baro - 1) * 0.8 + 800, 1)
    if me[29] == "1":
        m["selected_heading_deg"] = round(_u(me, 30, 39) * 180 / 256, 1)
    m["nac_p"], m["nic_baro"], m["sil"] = _u(me, 39, 43), int(me[43]), _u(me, 44, 46)
    if me[46] == "1":                                    # mode bits valid
        m["modes"] = [n for n, i in (("autopilot", 47), ("VNAV", 48), ("altitude hold", 49),
                                     ("approach", 51), ("LNAV", 53)) if me[i] == "1"]
    m["tcas_operational"] = me[52] == "1"


def _me_status(b, m):
    """TC31 aircraft operational status, DO-260B 2.2.3.2.7.2."""
    me = b[32:88]
    st = _u(me, 5, 8)
    m["status_subtype"] = "airborne" if st == 0 else "surface" if st == 1 else st
    m["adsb_version"] = _u(me, 40, 43)
    m["nic_a"], m["nac_p"], m["sil"] = int(me[43]), _u(me, 44, 48), _u(me, 50, 52)
    m["heading_reference"] = "magnetic north" if me[53] == "1" else "true north"
    m["sil_per"] = "sample" if me[54] == "1" else "hour"
    if m["adsb_version"] < 2:
        return
    if st == 0:
        m["gva"], m["nic_baro"] = _u(me, 48, 50), int(me[52])
        m["capabilities"] = {"tcas_operational": me[10] == "1", "1090es_in": me[11] == "1",
                             "air_referenced_velocity": me[14] == "1", "target_state_reports": me[15] == "1",
                             "trajectory_change_reports": _u(me, 16, 18), "uat_in": me[18] == "1"}
    else:
        m["length_width_code"] = _u(me, 20, 24)
        m["capabilities"] = {"1090es_in": me[11] == "1", "uat_in": me[15] == "1"}
    m["tcas_ra_active"] = me[26] == "1"
    m["ident"] = me[27] == "1"
    m["single_antenna"] = me[29] == "1"
    m["sda"] = _u(me, 30, 32)


def _extended_squitter(pms, h, b, m):
    df = m["df"]
    m["icao"] = pms.icao(h)
    if df == 17:
        m["capability"] = CAPABILITY.get(_u(b, 5, 8), f"CA {_u(b, 5, 8)}")
    else:
        m["control_field"] = CF[_u(b, 5, 8)]
    tc = m["typecode"] = pms.adsb.typecode(h)
    me = b[32:88]
    if 1 <= tc <= 4:
        m["callsign"] = pms.adsb.callsign(h).strip("_ ")
        m["category"] = CATEGORY.get((tc, _u(me, 5, 8)), "reserved")
    elif 5 <= tc <= 8:
        mov = pms.adsb.surface_velocity(h, source=True)
        m["ground_speed_kt"], m["track_deg"] = _r(mov[0]), _r(mov[1])
        _me_position(b, m, tc)
    elif 9 <= tc <= 18 or 20 <= tc <= 22:
        alt = pms.adsb.altitude(h)
        m["altitude_ft" if tc <= 18 else "gnss_altitude_ft"] = alt
        _me_position(b, m, tc)
    elif tc == 19:
        v = pms.adsb.velocity(h, source=True)
        m["velocity_subtype"] = _u(me, 5, 8)
        m["nac_v"] = _u(me, 10, 13)
        m["intent_change"] = me[8] == "1"
        if v:
            spd, ang, vr, kind, ang_src, vr_src = v
            if kind == "GS":
                m["ground_speed_kt"], m["track_deg"] = _r(spd), _r(ang)
            else:
                m[f"{kind.lower()}_kt"] = spd
                if ang is not None:
                    m["heading_deg"] = _r(ang)
                m["heading_source"] = ang_src
            m["vrate_fpm"], m["vrate_source"] = vr, vr_src
        m["gnss_minus_baro_ft"] = pms.adsb.altitude_diff(h)
    elif tc == 28:
        sub = m["status_subtype"] = _u(me, 5, 8)
        if sub == 1:
            m["emergency"] = EMERGENCY[_u(me, 8, 11)]
            m["squawk"] = pms.adsb.emergency_squawk(h)
        elif sub == 2:
            m["tcas_ra_broadcast"] = me[8:56]
    elif tc == 29:
        _me_target_state(b, m)
    elif tc == 31:
        _me_status(b, m)


def _decode(pms, m):
    h, df = m["hex"], m["df"]
    b = _bits(h)
    if df in (17, 18):
        _extended_squitter(pms, h, b, m)
        return
    m["icao"] = m["aa"]
    if df == 11:
        m["capability"] = CAPABILITY.get(_u(b, 5, 8), f"CA {_u(b, 5, 8)}")
        res = int(m["parity_residual"], 16)
        cl, ic = res >> 4, res & 0xF                     # code label + interrogator code
        if cl == 0:
            m["interrogator"] = f"II {ic}" + (" (or acquisition squitter)" if ic == 0 else "")
        else:
            m["interrogator"] = f"SI {ic + 16 * (cl - 1)}"
        m["interrogator_code"] = res
    elif df in (0, 16):                                  # ACAS (TCAS) air-air replies
        m["vertical_status"] = "on ground" if b[5] == "1" else "airborne"
        m["acas_sensitivity"] = ACAS_SL[_u(b, 8, 11)]
        m["acas_reply_info"] = ACAS_RI.get(_u(b, 13, 17), f"RI {_u(b, 13, 17)}")
        if df == 0:
            m["crosslink_capable"] = b[6] == "1"
        m["altitude_ft"] = pms.common.altcode(h)
        if df == 16:
            m["mv"] = h[8:22]
            if h[8:10] == "30":
                m["acas_ra_report"] = True
    elif df in (4, 5, 20, 21):                           # replies to Mode S ground radar
        m["flight_status"] = FLIGHT_STATUS[_u(b, 5, 8)]
        m["downlink_request"] = _u(b, 8, 13)
        um = _u(b, 13, 19)
        m["iis"], m["ids"] = um >> 2, um & 3            # interrogator identifier of the radar
        if df in (4, 20):
            m["altitude_ft"] = pms.common.altcode(h)
        else:
            m["squawk"] = pms.common.idcode(h)
        if df in (20, 21):
            m["commb"] = _commb(pms, h)


def fields(msgs):
    try:
        import pyModeS as pms
    except Exception:                                    # noqa: BLE001
        return msgs
    for m in msgs:
        try:
            _decode(pms, m)
        except Exception as e:                           # noqa: BLE001 — keep the raw frame
            m["field_error"] = str(e)[:80]
    return msgs
