"""Radio altimeter receiver (research, signal gaps S8): FMCW sweeps in 4200-4400 MHz measured from the crossings a
capture sees (rate, direction, turnarounds -> sweep bandwidth, period). See README.md.

    detect(x, fs, fc) -> list[Emission]     family "radalt", label "radio_altimeter"
    synth: synthetic FMCW altimeters (tests)
"""
from .detector import detect  # noqa: F401
