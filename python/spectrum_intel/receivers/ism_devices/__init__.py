"""315 / 433 MHz short-range devices (research, signal gaps S6): OOK fixed-code remotes (EV1527 / PT2262), KeeLoq
remotes (clear serial / buttons only), FSK TPMS (Continental / VDO, Toyota). See README.md.

    detect(ch) -> list[Emission]     family "ism_devices"
    synth: synthetic OOK / FSK devices (tests)
"""
from .detector import BANDS, detect  # noqa: F401
