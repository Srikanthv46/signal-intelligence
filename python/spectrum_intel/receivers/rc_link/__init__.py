"""RC / telemetry links and BLE: LoRa / ExpressLRS CSS, GFSK/FSK FHSS links, BLE
advertising (CRC-24 verified), and the capture frequency-label check from BLE adverts.

    detect(x, fs, fc, ctx=None) -> list[Emission]
    check_frequency_label(emissions) -> dict | None     label error from CRC-valid BLE adverts
    synth                                               synthetic LoRa / GFSK link generator (tests)

Depends on spectrum_intel.core only. See README.md.
"""
from .detector import detect, detect_fsk, detect_lora  # noqa: F401
from .label_check import BLE_ADV_MHZ, check_frequency_label  # noqa: F401
