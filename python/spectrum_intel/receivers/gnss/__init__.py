"""GNSS receiver (research, signal gaps S3): GPS L1 C/A acquisition from one block, C/N0, jamming and spoofing
indicators. See README.md.

    analyse(ch, ppm=None, exclude_hz=()) -> dict | None
    codes.ca_code(prn) -> the 1023-chip Gold code;  synth.sky(...) -> synthetic skies (tests)
"""
from .acquire import analyse  # noqa: F401
