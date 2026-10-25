"""ADS-B / Mode S field decoding and per-aircraft summaries, on frames from a live 1090 MHz Pluto
capture (2026-10-09, testdata/adsb/): IndiGo 6E311 climbing out of Hyderabad, and 6E6294."""
import json
from pathlib import Path

import numpy as np
import pytest

from spectrum_intel.receivers.adsb import demod, fields, tracks
from spectrum_intel.receivers.adsb.modes import SPU

pytest.importorskip("pyModeS")
D = Path(__file__).resolve().parents[2] / "testdata/adsb"


def _frame(h, df=None, res="000000", aa=None):
    df = int(h[:2], 16) >> 3 if df is None else df
    return fields([{"t": 0.0, "df": df, "hex": h, "parity_residual": res, "aa": aa}])[0]


def test_identity_and_category():
    m = _frame("8D801648232473F3C718206F859E")
    assert (m["icao"], m["callsign"], m["category"]) == ("801648", "IGO311", "large (75,000-300,000 lb)")
    assert m["capability"] == "level 2+ transponder, airborne"


def test_target_state_status_and_emergency():
    ts = _frame("8D801648EA2EF858013C08BD9323")
    assert (ts["selected_altitude_ft"], ts["selected_altitude_source"], ts["baro_setting_hpa"]) == (24000, "MCP/FCU", 1012.8)
    assert (ts["nac_p"], ts["sil"], ts["tcas_operational"]) == (9, 3, True)
    st = _frame("8D801648F82300030049B8B88048")
    assert (st["adsb_version"], st["nac_p"], st["gva"], st["sda"], st["single_antenna"]) == (2, 9, 2, 3, False)
    assert st["capabilities"]["tcas_operational"] and st["capabilities"]["target_state_reports"]
    em = _frame("8D801648E1063F00000000FD1401")
    assert (em["emergency"], em["squawk"]) == ("no emergency", "2727")


def test_velocity_and_reply_formats():
    v = _frame("8D8016489914741D50D00F879E03")
    assert (v["ground_speed_kt"], v["vrate_fpm"], v["vrate_source"], v["nac_v"]) == (259, 3264, "BARO", 2)
    acas = _frame("02A1871CA6D80D", aa="801648")
    assert acas["altitude_ft"] == 10500 and acas["acas_reply_info"] == "ACAS, vertical-only RAs"
    df11 = _frame("5D80164849E5B3", res="000000")
    assert df11["interrogator"].startswith("II 0")
    cb = _frame("A800063FEEBD9721FE3C85488827", aa="801648")
    assert cb["squawk"] == "2727" and cb["commb"]["bds"] == "BDS50"
    assert cb["commb"]["candidates"]["BDS50"]["tas_kt"] == 266


def test_tracks_on_live_frames():
    msgs = fields(json.load(open(D / "real_live_frames.json")))
    ac = {a["icao"]: a for a in tracks(msgs)}
    a = ac["801648"]
    assert (a["flight"], a["operator"], a["country"], a["adsb_version"]) == ("6E311", "IndiGo", "India", 2)
    assert a["n_positions"] >= 30 and abs(a["position"]["lat"] - 17.3645) < 0.01 and abs(a["position"]["lon"] - 78.492) < 0.01
    assert a["vertical_trend"] == "climbing" and a["turn_deg"] < -30            # left turn out of Hyderabad
    assert a["commb"]["registers"] == {"BDS40": 2, "BDS50": 2, "BDS60": 1}     # BDS50/60 ambiguity resolved
    assert a["met"]["wind_kt"] < 20 and -5 < a["met"]["static_air_temp_c"] < 20 and a["met"]["temp_uncertainty_c"] < 5
    assert a["interrogation"]["acas_replies"] == 18 and a["interrogation"]["mode_s_radar_replies"] == 5
    assert "flags" not in a
    b = ac["801903"]
    assert b["flight"] == "6E6294" and b["n_positions"] == 3 and "flags" not in b


def test_receiver_reference_gives_range_and_local_decode():
    msgs = fields(json.load(open(D / "real_live_frames.json")))
    one = [m for m in msgs if m["icao"] == "801648" and m.get("typecode") == 11][:1]  # a single frame: no pair
    a = tracks(one, ref=(17.40, 78.50, 500))[0]
    assert a["position"]["method"] == "local (receiver)"
    assert abs(a["position"]["lat"] - 17.339) < 0.01
    assert 5 < a["range_km"] < 10 and 130 < a["bearing_deg"] < 200 and a["elevation_deg"] > 10


def test_implausible_jump_is_flagged():
    msgs = fields(json.load(open(D / "real_live_frames.json")))
    pos = [m for m in msgs if m["icao"] == "801648" and m.get("typecode") == 11]
    late = dict(pos[-1], t=pos[-2]["t"] + 0.001)                                # 1 ms later, CPR latitude altered
    late["hex"] = late["hex"][:14] + "%05X" % ((int(late["hex"][14:19], 16) + 0x1000) & 0xFFFFF) + late["hex"][19:]
    a = tracks(pos[:-1] + [late])[0]
    assert any("implausible" in f for f in a.get("flags", []))


def test_end_to_end_on_real_iq():
    meta = json.load(open(D / "real_1090_commb.json"))
    x = np.fromfile(D / "real_1090_commb.cf32", np.complex64)
    a = np.abs(x)
    msgs = fields(demod(a, 6 * np.median(a) * SPU))
    assert [m["hex"] for m in msgs] == [f["hex"] for f in meta["frames"]]
    assert all(m["snr_db"] > 10 for m in msgs)

    from spectrum_intel.core.dsp import IQContext
    from spectrum_intel.receivers.adsb.plugin import AdsbDecoder
    sig = AdsbDecoder().capture_signals(IQContext(x, meta["fs"], meta["fc"]))[0]
    f = sig.decoded["fields"]
    assert f["n_frames"] == 6 and f["aircraft"][0]["icao"] == "801648"
    assert f["aircraft"][0]["interrogation"]["radar_dwells"] == 1


def _with_parity(h):
    """Mode S CRC-24 (generator 0xFFF409) appended to the first 88 bits of an extended squitter."""
    v = int(h[:22], 16) << 24
    for i in range(87, -1, -1):
        if v >> (i + 24) & 1:
            v ^= 0x1FFF409 << i
    return h[:22] + "%06X" % (v & 0xFFFFFF)


def test_uav_emitter_category():
    """Emitter category B6 (TC 3, CA 6) = unmanned aerial vehicle -> track flag 'uav' (research PLAN T7)."""
    h = _with_parity("8D801648" + "1E" + "2473F3C71820" + "000000")
    assert _with_parity("8D801648232473F3C718206F859E") == "8D801648232473F3C718206F859E"
    m = _frame(h)
    assert m["category"] == "unmanned aerial vehicle" and m["callsign"] == "IGO311"
    a = tracks([m])[0]
    assert a["uav"] is True and "UAV (ADS-B emitter category B6)" in a["flags"]
