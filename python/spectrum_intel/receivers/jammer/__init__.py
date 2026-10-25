"""Swept-carrier jammer receiver.

    detect(x, fs, fc, ctx=None) -> list[Emission]      family "jammer", label "swept_carrier"
    measure(x, fs, fc) -> (sweep | None, diagnostics)   the detector's internals (tests, bench)
    synth                                               synthetic sweeps and look-alikes (tests)

Depends on spectrum_intel.core only. See README.md.
"""
from .detector import detect, measure  # noqa: F401
