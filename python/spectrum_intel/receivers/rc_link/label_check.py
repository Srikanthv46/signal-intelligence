"""Capture centre-frequency label check from CRC-verified BLE advertising packets."""
from __future__ import annotations

import numpy as np


# anything not in DRONE_FAMILIES (Wi-Fi, Bluetooth, generic narrowband, plain LoRa -- often
# IoT) is reported separately and never forms a drone system
BLE_ADV_MHZ = {37: 2402.0, 38: 2426.0, 39: 2480.0}

def check_frequency_label(emissions):
    """A CRC-verified BLE advert can only pass when de-whitened with its true channel
    index, so its labelled frequency minus the true advertising-channel frequency is
    the capture's frequency-label error. Returns None without such evidence."""
    offs = [e.center_mhz - BLE_ADV_MHZ[e.decoded["whitening_channel"]] for e in emissions
            if e.label == "ble_adv" and (e.decoded or {}).get("crc_ok")
            and (e.decoded or {}).get("whitening_channel") in BLE_ADV_MHZ]
    if not offs:
        return None
    off = float(np.median(offs))
    return dict(offset_mhz=round(off, 3), n_adverts=len(offs), spread_mhz=round(float(np.ptp(offs)), 3),
                suspect=abs(off) > 0.5,
                note=("capture frequency label looks wrong by %+.1f MHz: rerun with --fc <labelled - %.1f MHz>"
                      % (off, off)) if abs(off) > 0.5 else "frequency label consistent with BLE adverts")
