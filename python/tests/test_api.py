"""Production API (spectrum_intel.api): the response contract, input handling and streaming,
on short real captures from testdata/ (ADS-B at 1090 MHz, a PMR446 walkie-talkie)."""
import json
from pathlib import Path

import jsonschema
import numpy as np
import pytest

from spectrum_intel.api import Engine, IQBlock, StreamAnalyzer, contract, technology

DATA = Path(__file__).resolve().parents[2] / "testdata"


def _vec(name):
    meta = json.loads((DATA / f"{name}.json").read_text())
    meta.setdefault("fs", meta.get("fs_hz"))
    meta.setdefault("fc", meta.get("fc_hz"))
    return np.fromfile(DATA / f"{name}.cf32", np.complex64), meta


@pytest.fixture(scope="module")
def engine():
    return Engine()


@pytest.fixture(scope="module")
def schema():
    return contract.load("response.schema.json")


@pytest.fixture(scope="module")
def adsb(engine):
    x, m = _vec("adsb/real_1090_commb")
    return x, m, engine.analyze(IQBlock(x, m["fs"], m["fc"], start_time="2026-10-09T08:00:00Z",
                                        source_id="balcony", request_id="r1", gain_db=m["gain_db"]))


def test_response_matches_schema(adsb, schema):
    _, _, r = adsb
    jsonschema.validate(r, schema)
    assert r["status"] == "ok" and r["request_id"] == "r1" and r["input"]["source_id"] == "balcony"
    assert r["engine"]["implementation"] == "python"


def test_adsb_decoded_with_identity(adsb):
    _, m, r = adsb
    s = [s for s in r["signals"] if s["technology"] == "mode_s_adsb"]
    assert s and s[0]["decode"]["status"] == "success" and s[0]["family"] == "aviation"
    assert s[0]["identifiers"].get("icao")
    assert s[0]["time"]["start_utc"].startswith("2026-10-09T08:00:00")
    assert r["summary"]["by_status"].get("success", 0) >= 1


def test_coverage_and_warnings(adsb):
    _, _, r = adsb
    cov = {c["receiver"]: c for c in r["coverage"]}
    assert cov["adsb"]["in_band"] and cov["adsb"]["duration_ok"] and not cov["adsb"]["full_duration"]
    assert not cov["wifi"]["in_band"]
    assert not any(w["code"] == "SHORT_BLOCK" for w in r["warnings"])     # 272k samples, 34 ms


def test_short_block_is_refused_not_crashed(engine, schema):
    x, m = _vec("adsb/real_1090_commb")
    for n in (2048, 4096, 65535):
        r = engine.analyze(IQBlock(x[:n], m["fs"], m["fc"]))
        jsonschema.validate(r, schema)
        assert r["status"] == "insufficient_input" and r["errors"][0]["code"] == "INVALID_INPUT"
        assert r["input"]["n_samples"] == n


def test_bad_metadata_is_refused(engine):
    x, m = _vec("adsb/real_1090_commb")
    assert engine.analyze(IQBlock(x, -1.0, m["fc"]))["status"] == "insufficient_input"
    assert engine.analyze(IQBlock(x, m["fs"], float("nan")))["status"] == "insufficient_input"
    y = x.copy()
    y[1000] = np.nan
    assert engine.analyze(IQBlock(y, m["fs"], m["fc"]))["status"] == "insufficient_input"


def test_interleaved_int16_equals_complex(engine):
    x, m = _vec("voice/real_pmr446_ch8_dcs")
    xi = np.round(x).astype(np.complex64)                  # the capture holds raw ADC counts
    inter = np.empty(2 * len(xi), np.int16)
    inter[0::2], inter[1::2] = xi.real, xi.imag
    a = Engine(enable=["voice_fm"]).analyze(IQBlock(xi, m["fs"], m["fc"]))
    b = Engine(enable=["voice_fm"]).analyze(IQBlock(inter, m["fs"], m["fc"], iq_format="ci16"))
    strip = lambda r: [(s["label"], s["frequency"]["center_hz"], s["decode"]["status"]) for s in r["signals"]]
    assert a["status"] == b["status"] == "ok" and strip(a) == strip(b)
    assert any(s["technology"] == "fm_voice" for s in a["signals"])


def test_bandwidth_limits_analysed_band(engine):
    x, m = _vec("adsb/real_1090_commb")
    r = Engine(enable=["adsb"]).analyze(IQBlock(x, m["fs"], m["fc"], bandwidth_hz=4e6))
    lo, hi = r["input"]["analysed_band_hz"]
    assert hi - lo == pytest.approx(4e6)


def test_stream_accumulates_small_blocks():
    x, m = _vec("adsb/real_1090_commb")
    st = StreamAnalyzer(Engine(enable=["adsb"]), target_samples=131072)
    out = []
    for i in range(0, 2 * 131072 + 4096 * 2, 4096):         # 4096-sample SDR buffers (x: 272000)
        out += st.push(x[i:i + 4096], m["fs"], m["fc"], start_time=1_760_000_000 + i / m["fs"])
    assert len(out) == 2 and all(r["input"]["n_samples"] == 131072 for r in out)
    t0 = [r["input"]["start_time_utc"] for r in out]
    assert t0[0] != t0[1]
    assert st.buffered == 2 * 4096
    assert st.flush() is None and st.buffered == 0          # 8192 samples: below the minimum, dropped


def test_stream_closes_block_on_retune():
    x, m = _vec("adsb/real_1090_commb")
    st = StreamAnalyzer(Engine(enable=["adsb"]), target_samples=200000)
    assert st.push(x[:100000], m["fs"], m["fc"]) == []
    out = st.push(x[100000:110000], m["fs"], m["fc"] + 1e6)   # retune: the first 100000 are analysed alone
    assert len(out) == 1 and out[0]["input"]["n_samples"] == 100000 and st.buffered == 10000


def test_receiver_ppm_reaches_the_receivers():
    """The caller's ppm (the Pluto's +7.44) reaches the receivers; the GSM receiver also finds the BCCH carriers of
    this 928 MHz snapshot (testdata/cellular) without it: the FCCH tone is searched over the carrier offset and every
    FCCH of the block measures the same receiver error (~+7.8 ppm)."""
    from spectrum_intel.io.sigmf import load
    cap = load(DATA / "cellular/SWEEP_20261004_141028_317_928.08MHz")
    eng = Engine(enable=["lte", "gsm"])
    blind = eng.analyze(IQBlock(cap.iq, cap.fs, cap.fc))
    cal = eng.analyze(IQBlock(cap.iq, cap.fs, cap.fc, receiver_ppm=7.44))
    for r in (cal, blind):
        gsm = [s for s in r["signals"] if s["technology"] == "gsm" and s["decode"]["status"] != "not_possible"]
        assert len(gsm) >= 5 and all("arfcn" in s["identifiers"] for s in gsm)
        assert all(abs(s["decode"]["fields"]["receiver_ppm_fcch"] - 7.8) < 0.5 for s in gsm)
    assert cal["input"]["receiver_ppm"] == 7.44 and blind["input"]["receiver_ppm"] is None


@pytest.mark.parametrize("label,decoder,tech", [
    ("LTE B3 FDD downlink PCI 12", "lte", "lte"),
    ("LTE/NR-like TDD B40/n40", None, "lte_or_nr"),
    ("5G NR n78 PCI 82 (SSB 3361.44 MHz)", "nr_ssb", "nr"),
    ("Wi-Fi 802.11b DSSS ch6", "wifi_dsss", "wifi_dsss"),
    ("Wi-Fi 802.11n-type link, 11 MHz channel at 2069.02 MHz", "wifi_ofdm_scaled", "wifi_nonstandard_clock"),
    ("Bluetooth LE advertising (ADV_IND)", "rc_link", "ble"),
    ("Drone Remote ID (Bluetooth LE, ADV_NONCONN_IND)", "rc_link", "remote_id_ble"),
    ("SiK telemetry radio (MAVLink)", "rc_link", "sik_telemetry"),
    ("IEEE 802.15.4 frames (Digi International (XBee))", "rc_link", "ieee802154"),
    ("IEEE 802.15.4-like (O-QPSK 2.4 GHz)", "rc_link", "ieee802154"),
    ("TBS Crossfire-like (150 Hz FSK, sync 0x3E6B3E)", "rc_link", "crossfire"),
    ("FLARM-like (868 MHz collision-avoidance beacon)", "rc_link", "flarm"),
    ("FrSky D16 v2 (ACCST)", "rc_link", "frsky"),
    ("FrSky D8-like (ACCST)", "rc_link", "frsky"),
    ("FlySky AFHDS 2A (A7105)", "rc_link", "afhds2a"),
    ("Bayang (XN297)", "rc_link", "xn297_toy"),
    ("XN297-like (toy drone)", "rc_link", "xn297_toy"),
    ("FlySky AFHDS 2A-like (A7105 ID code)", "rc_link", "afhds2a"),
    ("Drone Remote ID (Bluetooth LE Long Range, AUX_ADV_IND)", "rc_link", "remote_id_ble"),
    ("COFDM video link (DVB-T 2K, 8 MHz, QPSK, HP 1/2, GI 1/8)", "cofdm", "cofdm_video_link"),
    ("COFDM video link-like (DVB-T 8K numerology, 8 MHz, GI 1/32)", "cofdm", "cofdm_video_link"),
    ("Analog FPV video (PAL-M)", "analog_vtx", "analog_fpv_video"),
    ("ExpressLRS lora sf7 bw500k", "rc_link", "expresslrs"),
    ("mLRS LoRa sf6 bw500k", "rc_link", "mlrs"),
    ("mLRS like LoRa sf5 bw812.5k", "rc_link", "mlrs"),
    ("LoRa sf9 bw125k", "rc_link", "lora"),
    ("2.4 GHz GFSK link (ESB-like, 2 Mb/s)", "rc_link", "esb_like_link"),
    ("gfsk 250k hop", "rc_link", "fsk_rc_link"),
    ("DJI DroneID", "dji", "dji_droneid"),
    ("HDZero digital FPV video link", "digital_vtx", "hdzero"),
    ("RadioLink (air)", "rc_link", "radiolink"),
    ("Hubsan (A7105)", "rc_link", "hubsan"),
    ("Toy drone RC link: MJXQ (WLH08 / E010 class) (XN297)", "rc_link", "toy_drone"),
    ("Herelink (CubePilot) drone video/control link", "digital_vtx", "herelink"),
    ("unknown OFDM video/data link (44.7 kHz SCS)", "digital_vtx", "ofdm_data_link"),
    ("LO leakage (receiver artifact)", None, "lo_leakage"),
    ("Jammer (swept carrier)", "jammer", "swept_jammer"),
    ("Swept carrier (FMCW-like)", "jammer", "swept_carrier"),
    ("unknown", None, "unknown"),
])
def test_technology_codes(label, decoder, tech):
    assert technology(label, decoder)[0] == tech
