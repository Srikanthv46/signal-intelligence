"""FM broadcast receiver (research, signal gaps S20): carriers on the 88-108 MHz raster, deviation, stereo pilot and
RDS (PI, PTY, programme name, RadioText, clock time). See README.md.

    detect(ch) -> list[Emission]   family "broadcast_fm", label "fm_broadcast"
    rds: the RDS block code and group decoder;  synth: synthetic stations (tests)
"""
from .detector import detect  # noqa: F401
