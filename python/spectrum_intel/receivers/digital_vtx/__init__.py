"""Digital video links that are not DJI: Herelink (LTE Rel-8 DL PHY in 2.4 GHz ISM,
PSS/SSS -> PCI) and blind unknown-OFDM video/telemetry links.

    detect(x, fs, fc, ctx=None) -> list[Emission]      family herelink | digital_vtx | wifi_ofdm
        ctx["report_wifi"]=False drops Wi-Fi-numerology bursts

Depends on spectrum_intel.core only. See README.md.
"""
from .detector import HERELINK, KNOWN_NUMEROLOGIES, classify_numerology, detect, generic_detect, herelink_detect  # noqa: F401
