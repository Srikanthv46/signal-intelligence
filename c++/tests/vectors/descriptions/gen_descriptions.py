#!/usr/bin/env python3
"""Reference for tests/unit/test_descriptions.cpp: api.contract.description(label, decoder, flags) on every distinct
(label, decoder, flags) of the latest regression reports (out/regress/wbatch) plus one sample label per
schema/descriptions.json variant.

    .venv/bin/python -I c++/tests/vectors/descriptions/gen_descriptions.py
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
sys.path.insert(0, str(REPO / "python"))
from spectrum_intel.api import contract as K                     # noqa: E402

SAMPLES = [("Wi-Fi 802.11a/g-type link, 5 MHz channel at 2432.03 MHz (Sz Dji Technology Co.,Ltd), drone link-like (DJI-like)",
            "wifi_ofdm_scaled", []),
           ("Wi-Fi 802.11n-type link, 11 MHz channel at 2069.02 MHz, spectrum-inverted (Microhard Systems Inc.)",
            "wifi_ofdm_scaled", []),
           ("Wi-Fi 802.11b DSSS ch6, drone link-like (DJI Tello)", "wifi_dsss", []),
           ("Wi-Fi 802.11 OFDM (NON_HT) 2g6, drone link-like (by behaviour)", "wifi_ofdm", []),
           ("OFDM fixed-wireless link-like (PMP 450-like numerology, 43.75 kHz SCS)", "digital_vtx", []),
           ("OFDM RC uplink-like (DJI Lightbridge / WiMAX-like numerology, 10.94 kHz SCS)", "digital_vtx", []),
           ("OFDM hopping RC uplink-like (DJI O4-like numerology, 30 kHz SCS)", "digital_vtx", []),
           ("OFDM drone video link-like (DJI O4-like numerology, 30 kHz SCS)", "digital_vtx", []),
           ("unknown OFDM TDD link (25 kHz SCS)", "digital_vtx", []),
           ("Analog FM voice (PMR446 ch 8, D053N)", "voice_fm", []),
           ("Pulsed radar-like", None, []), ("Pulsed radar (S-band ATC)-like", None, []),
           ("Wi-Fi 802.11 OFDM (NON_HT) 5g149", "wifi_ofdm", ["iq_image"]),
           ("COFDM video link-like (DVB-T 2K numerology, 8 MHz, GI 1/4)", "cofdm", []),
           ("DVB-T2 broadcast", None, []), ("Radiosonde-like telemetry", None, []), ("Iridium-like TDMA burst", None, []),
           ("Jammer (swept carrier)", "jammer", []), ("Swept carrier (FMCW-like)", "jammer", []),
           ("Drone Remote ID (Bluetooth LE, ADV_NONCONN_IND)", "rc_link", []), ("DJI DroneID", "dji", []),
           ("TBS Crossfire 150 Hz (FSK)", "rc_link", []), ("FrSky D16 v2 (ACCST)", "rc_link", []),
           ("Toy drone (nRF24 / XN297)", "rc_link", []), ("SiK telemetry radio 64 kb/s", "rc_link", []),
           ("UAT 978 MHz ADS-B", "adsb", []), ("LO-harmonic image", None, []), ("something new", "rc_link", []),
           # signal gaps batch S (research)
           ("SSR Mode A/C replies (squawk 7700 emergency)", "adsb", []), ("SSR Mode A/C replies-like", "adsb", []),
           ("Analog AM voice (VHF air band, 118.25 MHz)", "voice_fm", []),
           ("Analog AM voice (aeronautical emergency 121.5 MHz)", "voice_fm", []),
           ("ILS glide slope (330.95 MHz, paired LOC 111.95 MHz)", "navaids", []),
           ("ILS localizer-like (109.90 MHz, paired GS 333.80 MHz, two-frequency)", "navaids", []),
           ("VOR beacon (113.10 MHz)", "navaids", []), ("GNSS: GPS L1 C/A (5 satellites)", "gnss", []),
           ("GNSS spoofing-like (GPS L1 C/A, 5 satellites)", "gnss", []), ("GNSS jamming-like (GPS L1)", "gnss", []),
           ("Bluetooth LE tracker (Apple Find My / AirTag)", "rc_link", []),
           ("FM broadcast (98.3 MHz, stereo, RDS 'AIR FM')", "broadcast_fm", []),
           ("FM broadcast (97.6 MHz)", "broadcast_fm", []), ("AIS (3 stations)", "marine", []),
           ("APRS (VU2ABC-9)", "voice_fm", []), ("AX.25 packet radio (VU2BCN)", "voice_fm", []),
           ("Pulsed radar (S-band ATC, solid-state)-like", None, []), ("Pulsed radar (weather, S-band)-like", None, []),
           ("Pulsed radar (weather, C-band)-like", None, []), ("Pulsed radar (marine S-band)-like", None, []),
           ("Pulsed radar (C-band)-like", None, []), ("Pulsed radar (L-band surveillance)-like", None, []),
           ("Wind profiler radar-like", None, [])]


def main():
    seen, cases = set(), []
    for f in sorted((REPO / "out/regress/wbatch").glob("*.json")):
        for s in json.loads(f.read_text()).get("signals", []):
            k = (s["label"], s.get("features", {}).get("decoder"), tuple(s.get("flags") or []))
            if k not in seen:
                seen.add(k)
    for lab, dec, fl in SAMPLES:
        seen.add((lab, dec, tuple(fl)))
    for lab, dec, fl in sorted(seen, key=lambda k: (k[0], str(k[1]), k[2])):
        cases.append({"label": lab, "decoder": dec, "flags": list(fl), "text": K.description(lab, dec, list(fl))})
    (HERE / "cases.json").write_text(json.dumps(cases, indent=0, ensure_ascii=False))
    print(len(cases), "cases;", sum(not c["text"] for c in cases), "empty")


main()
