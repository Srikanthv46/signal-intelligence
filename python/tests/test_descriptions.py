"""Plain-language signal descriptions (schema/descriptions.json, api.contract.description): every technology has
one (drones named only where the signal is drone-specific), label variants pick the right text, artifacts are described as artifacts, and reports / responses carry it."""
import json
import re

from spectrum_intel.api import contract as K


def test_every_technology_described():
    techs = {r["technology"] for r in K.load("technologies.json")["rules"]}
    entries = K.load("descriptions.json")["entries"]
    assert techs <= {e["technology"] for e in entries}
    for e in entries:
        t = e["text"]
        assert 60 <= len(t) <= 420 and t.endswith("."), e
        assert len(re.findall(r"[.!?](\s|$)", t)) <= 3, e                  # two or three plain sentences
        assert "not a drone" not in t.lower(), e                       # drones only where the signal is drone-specific


def test_variants():
    d = K.description
    assert d("OFDM fixed-wireless link-like (PMP 450-like numerology, 43.75 kHz SCS)", "digital_vtx").startswith(
        "Fixed-wireless internet link")
    assert "Walksnail" in d("OFDM drone video link-like (Walksnail / WiMAX-like numerology, 10.94 kHz SCS)",
                            "digital_vtx")
    assert "Enhanced Wi-Fi" in d("Wi-Fi 802.11a/g-type link, 5 MHz channel at 2432.03 MHz (Sz Dji Technology "
                                 "Co.,Ltd), drone link-like (DJI-like)", "wifi_ofdm_scaled")
    assert d("Wi-Fi 802.11 OFDM (NON_HT) 2g6, drone link-like (by behaviour)", "wifi_ofdm").startswith(
        "Wi-Fi traffic that includes a drone link")
    assert d("Wi-Fi 802.11 OFDM (NON_HT) 2g6", "wifi_ofdm").startswith("Wi-Fi (802.11")
    assert d("Wi-Fi 802.11 OFDM (NON_HT) 2g6", "wifi_ofdm", ["iq_image"]).startswith("Not a real transmitter")
    assert "PMR446" in d("Analog FM voice (PMR446 ch 8, D053N)", "voice_fm")
    assert "not decoded" in d("LTE-like TDD B40/n40", None)
    assert d("never seen before", None) == K.description("unknown", None)


def test_report_and_response_carry_it():
    from pathlib import Path
    sch = json.loads((Path(K.schema_dir()) / "response.schema.json").read_text())
    assert "description" in json.dumps(sch)
    assert K.SCHEMA_VERSION == "1.4.0"
