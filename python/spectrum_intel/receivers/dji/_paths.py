"""Where optional external pieces live."""
from __future__ import annotations

import os
from pathlib import Path

DRONEID_ENV = "SPECTRUM_INTEL_DRONEID_SRC"


def droneid_src() -> str:
    """The DJI DroneID decoder (``droneid-decoder-reference/src``, which holds
    ``droneid_api.py``; a separate repo, not vendored). ``$SPECTRUM_INTEL_DRONEID_SRC``
    wins; otherwise the first ``droneid-decoder-reference/src`` beside any parent of
    this file. Absent is fine: DroneID bursts are still detected, just not decoded."""
    env = os.environ.get(DRONEID_ENV)
    if env:
        return env
    for parent in Path(__file__).resolve().parents:
        cand = parent / "droneid-decoder-reference" / "src"
        if (cand / "droneid_api.py").is_file():
            return str(cand)
    return ""
