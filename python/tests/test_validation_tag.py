"""schema/capabilities.json: the proven / research validation tag carried by every signal."""
from spectrum_intel.api import contract as K


def test_levels():
    a = K.validation("analog_fpv_video", "analog_vtx", {"standard": "pal"})
    assert a["level"] == "proven" and "proof" in a and "research_fields" not in a
    d = K.validation("ofdm_data_link", "digital_vtx", {"tx_antennas": {"verdict": "2tx_stbc_like"}})
    assert d["level"] == "proven" and d["research_fields"] == ["tx_antennas"]
    assert d["research_sources"] == ["bench/tx_antennas/FINDINGS.md"]
    assert "research_fields" not in K.validation("ofdm_data_link", "dji", {"tx_antennas": 2})
    assert K.validation("newthing", None, {})["level"] == "research"


def test_registry_complete():
    reg = K.load("capabilities.json")
    for r in K.load("technologies.json")["rules"]:
        assert r["technology"] in reg["technologies"], r["technology"]
    for t in reg["technologies"].values():
        assert t["level"] in reg["levels"]
        assert t.get("proof") if t["level"] == "proven" else t.get("sources")
    for c in reg["field_capabilities"]:
        assert c["level"] in reg["levels"] and c["fields"] and (c.get("proof") or c.get("sources"))
