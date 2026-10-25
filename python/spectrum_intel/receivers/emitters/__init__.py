"""Emitter association: which bursts (Emission records from any receiver) come from one
radio, which radios form one drone system (frame-clock lattice, oscillator ppm).

    measure_fingerprints(x, fs, fc, emissions)          refines timing/carrier features in place
    associate(emissions, fc, models=None) -> list[Emitter]
    extract_bursts(x, fs, fc) -> list[Emission]         15 kHz-OFDM bursts for standalone use
    DRONE_FAMILIES, _assoc_models()                     which families are drone links + their models

Depends on spectrum_intel.core only. See README.md.
"""
from .association import MODELS, Emitter, associate, measure_fingerprints, same_emitter_score  # noqa: F401
from .bursts import extract_bursts  # noqa: F401
from .families import DRONE_FAMILIES, LABEL_ALIASES, _assoc_models  # noqa: F401
