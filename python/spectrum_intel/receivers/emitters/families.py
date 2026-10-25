"""Which emission families are drone links, and the timing/oscillator model for each
(detector label -> drone_assoc.MODELS key). From the research drone_survey module."""
from __future__ import annotations


# what an emission family means for the operator
DRONE_FAMILIES = {"dji_ocusync": "DJI drone link", "dji_lightbridge": "DJI drone link",
                  "dji_droneid": "DJI DroneID", "herelink": "Herelink drone link",
                  "analog_vtx": "analog FPV video", "digital_vtx": "unknown digital video/datalink",
                  "rc_fhss": "RC control link", "telemetry": "telemetry link"}

# detector labels -> keys of drone_assoc.MODELS (timing/oscillator models)
LABEL_ALIASES = {("dji_ocusync", "uplink"): ("dji_ocusync", "ocusync_uplink")}

def _assoc_models():
    from . import association as A
    m = dict(A.MODELS)
    for k, v in LABEL_ALIASES.items():
        m.setdefault(k, A.MODELS[v])
    return m
