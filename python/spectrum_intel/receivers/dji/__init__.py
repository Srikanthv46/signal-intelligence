"""DJI OcuSync / Lightbridge links: video downlink (15/30 kHz CP-OFDM, ZC roots), RC
uplink (72-subcarrier ZC hops), DroneID (decoded by the external droneid_api, see _paths).

    detect(x, fs, fc, ctx=None) -> list[Emission]      ctx["decode_droneid"]=True to decode DroneID
    droneid_decode(x, fs, fc, center_hz, t0_s, t1_s) -> dict

Depends on spectrum_intel.core only (+ optional droneid-decoder-reference). See README.md.
"""
from .detector import VIDEO_PROFILES, detect, droneid_decode, zc_root_search  # noqa: F401
