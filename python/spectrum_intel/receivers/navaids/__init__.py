"""Aeronautical navigation aids receiver (research, signal gaps S9 / S21): ILS glide slope and localizer (90 / 150 Hz
DDM / SDM, channel pairing) and VOR (30 Hz AM / 9960 Hz FM bearing). See README.md.

    detect(ch, ppm=0.0) -> list[Emission]     family "navaids", label ils_gs | ils_loc | vor
    synth: synthetic ILS / VOR carriers (tests)
"""
from .detector import ILS_PAIRS, channels, detect  # noqa: F401
