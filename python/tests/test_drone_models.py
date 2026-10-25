"""Model-level identification (research, PLAN M0-M2; research/drone_transmitters/11_model_level_identification.md):
the common model block, DJI product types, CTA-2063-A manufacturer codes, Remote ID and DroneID model blocks."""
from spectrum_intel.core import drone_models as M
from spectrum_intel.core import remote_id as RID


def test_model_block_never_picks_one_of_several():
    b = M.model_block("Walksnail Avatar HD", ["Mini 1S", "Nano V3"], ["tx_antennas 1"])
    assert b["model"] is None and b["candidates"] == ["Mini 1S", "Nano V3"] and "unresolved_reason" in b
    assert M.model_block("Walksnail Avatar HD", ["Moonlight"], [])["model"] == "Moonlight-like"   # inferred: -like
    assert M.model_block("DJI", ["Mini 2"], [], "decoded")["model"] == "Mini 2"


def test_dji_product_types():
    assert M.dji_product(63) == "Mini 2" and M.dji_product(66) == "Air 2S" and M.dji_product("68") == "Mavic 3"
    assert M.dji_product(62) is None and M.dji_product(98) is None          # unknown: never the nearest name


def test_serial_manufacturer():
    assert M.serial_manufacturer("1581F45T0000000000AB") == dict(code="1581", maker="DJI", kind="aircraft")
    assert M.serial_manufacturer("1596A1234567890")["kind"] == "module"     # Dronetag add-on module
    assert M.serial_manufacturer("1581F45T") is None                        # length character disagrees
    assert M.serial_manufacturer("1581FO5T0000000000AB") is None            # O is not in the character set
    assert M.serial_manufacturer("9999312A") is None                        # unknown code


def test_remote_id_basic_id_manufacturer():
    body = bytes([0x12]) + b"1748C123456789ABC".ljust(20, b"\x00") + bytes(3)   # serial, multirotor
    bid = RID._basic_id(body)
    assert bid["manufacturer"]["maker"] == "Autel Robotics"
    b = M.remote_id_model(bid)
    assert b["family"] == "Autel Robotics" and b["model"] is None and b["basis"] == "decoded"
    bid = RID._basic_id(bytes([0x22]) + b"1748C123456789ABC".ljust(20, b"\x00") + bytes(3))  # CAA reg., not serial
    assert "manufacturer" not in bid


def test_dji_plugin_model():
    from spectrum_intel.core.emission import Emission
    from spectrum_intel.receivers.dji.plugin import PLUGINS
    e = Emission("dji_droneid", "droneid", 0, 600, 2440, 2449, 2444.5, 20, features=dict(
        droneid_frame="9-symbol (OcuSync 2/3 generation: Mini 2, Air 2 and later)", device_type=63))
    m = PLUGINS[0].summarise([e], None)["fields"]["model"]
    assert m["model"] == "Mini 2" and m["basis"] == "decoded"
    e.features["device_type"] = 62
    m = PLUGINS[0].summarise([e], None)["fields"]["model"]
    assert m["model"] is None and "62" in m["evidence"][0]
    del e.features["device_type"]
    m = PLUGINS[0].summarise([e], None)["fields"]["model"]
    assert m["model"] is None and m["basis"] == "inferred"


def test_fpv_model_candidates():
    """digital_vtx catalog -> model block (PLAN M3 / M4): antenna count, DJI O3 / O4 centres, HDZero 17 MHz."""
    from spectrum_intel.receivers.digital_vtx import catalog as C

    def mb(tu, cp, bw, f, tx, fr=None):
        return C.model_candidates(C.candidates(tu, cp, bw, f, tx, fr), f, bw, tx)

    w2 = mb(91.43, 0.125, 18.9, 5805.0, 2, 7.163)                   # 2-TX STBC Walksnail
    assert w2["family"] == "Walksnail Avatar HD" and "Moonlight" in w2["candidates"]
    assert "Mini 1S" not in w2["candidates"] and w2["model"] is None
    w = mb(91.43, 0.125, 18.9, 5805.0, None)                        # antennas unknown: both kits stay
    assert {"Mini 1S", "Moonlight"} <= set(w["candidates"])
    o4 = mb(33.33, 0.0703, 18.0, 5789.5, 2, 2.5)
    assert o4["family"] == "DJI O4" and "O4 Air Unit (Lite)" not in o4["candidates"]
    assert any("O4 centre" in e for e in o4["evidence"])
    o3 = mb(66.67, 0.0703, 18.0, 5804.5, 2)
    assert o3["family"] == "DJI O3" and "O3 Air Unit" in o3["candidates"]
    hz = mb(37.926, 0.03125, 17.0, 5806.0, None)                     # measured 27 MHz-mode numerology at 17 MHz
    assert any("Nano 90" in m for m in hz["candidates"])
    assert C.model_candidates([], 5800.0, 18.0, None) is None


def test_rc_link_models():
    """rc_link model rules (PLAN M5): FrSky / AFHDS / Crossfire / SiK / ELRS hardware class."""
    from spectrum_intel.core.emission import Emission
    from spectrum_intel.receivers.rc_link.models import model

    def em(label, mhz, decoded=None, **feats):
        return [Emission("rc_fhss", label, 0, 100, mhz - 0.1, mhz + 0.1, mhz, 20, decoded=decoded, features=feats)]

    m = model(em("frsky_d16", 2440, dict(variant="D16 v2", crc_ok=True)))
    assert m["family"] == "FrSky D16 v2" and m["model"] is None and "v2" in m["evidence"][0]
    assert model(em("afhds2a", 2440, dict(crc_ok=True)))["family"] == "FlySky AFHDS 2A"
    m = model(em("crossfire_sync", 868.3, dict(grid=dict(band="868", channel=3, link="uplink"))))
    assert m["family"] == "TBS Crossfire" and "Nano TX" in m["candidates"]
    assert model(em("sik_telemetry", 915, dict(air_rate_bps=64000)))["family"] == "SiK telemetry"
    assert "RFD900x" in model(em("sik_telemetry", 915, dict(air_rate_bps=125000)))["candidates"]
    m = model(em("elrs_lora_sf5_bw500k", 915.5, modulation="lora", sf=5, bw_hz=500e3, elrs_channel=("FCC915", 7)))
    assert m["model"] == "LR1121-based module-like" and any("FCC915" in e for e in m["evidence"])
    m = model(em("elrs_lora_sf6_bw500k", 2440.4, modulation="lora", sf=6, bw_hz=812.5e3))
    assert m["model"] is None and m["candidates"] == ["SX1280-based module", "LR1121-based module"]
    assert model(em("elrs_flrc_like", 2440, modulation="gfsk"))["model"] == "SX1280-based module-like"
    assert model(em("nb_burst", 2440)) is None


def test_analog_vtx_model():
    from spectrum_intel.receivers.analog_vtx.plugin import _model
    m = _model({"channel": "R1", "audio_subcarrier_mhz": 6.5}, 5658.0)
    assert m["model"] is None and m["candidates"] == [] and "microphone" in m["evidence"][-1]
    assert "pit mode" in _model({}, 5584.0)["model"]                  # single candidate: "-like"
    assert len(_model({}, 1280.0)["candidates"]) == 2


def test_dji_hop_model():
    """Uplink hop timing -> model candidates (PLAN M6 / D3; RFUAV, one unit per model: low confidence)."""
    from spectrum_intel.receivers.digital_vtx.catalog import hop_model
    m = hop_model(6.6, 410.0, 2.48)
    assert m["family"] == "DJI (hop timing)" and m["model"] is None
    assert m["candidates"] == ["DJI Mini 4 Pro", "DJI O4 aircraft (model unresolved)"]
    assert hop_model(7.0, 410.0, 2.0)["candidates"] == ["DJI Avata 2", "DJI O4 aircraft (model unresolved)"]
    assert hop_model(4.8, 940.8, 7.17) is None                 # the Walksnail-like uplink of LIVE_OFDM_VTX_5805MHz
    assert hop_model(1.14, 513.0, 4.0)["model"] == "DJI OcuSync 2 aircraft (Mini 2-class; model unresolved)-like"
    assert hop_model(2.96, 520.0, 5.16)["family"] == "Herelink (hop timing)"
    assert hop_model(6.06, 250.0, 2.97)["candidates"] == ["Skydroid H12", "Skydroid T10"]   # within tolerance of each other
    assert hop_model(3.32, 1550.0, 4.99)["candidates"] == ["Siyi MK15"]


def test_faa_doc_tables():
    """File 13 (PLAN M7a-c): DOC manufacturer codes, serial prefix -> model, DJI product codes > 70."""
    assert M.MFR_CODES["1809"][0].startswith("Aerobits") and M.MFR_CODES["1971"] == ("FIMI", "aircraft")
    assert M.serial_model("1581F6Z9A00000000001")["models"] == ["Mini 4 Pro"]
    assert M.serial_model("1581F446A00000000001")["models"] == ["Agras T10", "Agras T30"]     # shared DOC prefix
    assert M.serial_model("1869AF11P1234")["verified"] is False
    assert M.serial_model("9999") is None
    b = M.remote_id_model({"uas_id": "1581F6Z9A00000000001", "manufacturer": M.serial_manufacturer("1581F6Z9A00000000001")})
    assert b["model"] == "Mini 4 Pro" and b["basis"] == "decoded"
    assert M.dji_product(93) == "Mini 4 Pro" and "not yet confirmed" in M.dji_product_source(93)
    assert M.dji_product(71) is None


def test_dji_legacy_frame_candidates():
    from spectrum_intel.core.emission import Emission
    from spectrum_intel.receivers.dji.plugin import PLUGINS
    e = Emission("dji_droneid", "droneid", 0, 600, 2440, 2449, 2444.5, 20,
                 features=dict(droneid_frame="8-symbol (first symbol absent: ...)"))
    m = PLUGINS[0].summarise([e], None)["fields"]["model"]
    assert "Mavic 2 Pro" in m["candidates"] and "Mini 4 Pro" in m["candidates"] and m["model"] is None
    e.features["device_type"] = "Mini 2"                    # the external decoder's name (real Mini 2 decodes)
    m = PLUGINS[0].summarise([e], None)["fields"]["model"]
    assert m["model"] == "Mini 2" and m["basis"] == "decoded"
