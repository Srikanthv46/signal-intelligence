"""Marine VHF receiver (research, signal gaps S15): AIS (GMSK 9600 bit/s, HDLC, CRC-16/X.25) on 161.975 /
162.025 MHz, message types 1/2/3, 4, 5, 18, 21, 24. See README.md.

    detect(ch, ppm=0.0) -> list[Emission]   family "marine", label "ais"
    ais: the codec;  synth: synthetic AIS bursts (tests)
"""
from .detector import detect  # noqa: F401
