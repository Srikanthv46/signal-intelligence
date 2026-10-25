"""
flarm.py -- FLARM collision-avoidance beacons, detection and the clear header only (research, PLAN T5;
research/drone_transmitters/10_cc2500_a7105_xn297_flarm_xbee.md section 5). Synthetic validation only.

PHY: 868.2 / 868.4 MHz (Europe / Africa / Asia; 866.0 MHz in the SoftRF India plan), 2-FSK +-50 kHz, 100 kchip/s,
Manchester coded (chip pair 01 -> bit 1, 10 -> bit 0) = 50 kb/s. After the preamble, the 32-bit sync F5 31 FA B6
(0x31FAB6 is the static start of the header), then the 24-bit radio ID (little-endian) and the address type (0x20 =
FLARM / manufacturer ID, 0x10 = ICAO 24-bit address). The rest of the packet is obfuscated: it is not decoded here.
Nothing is transmitted.
"""
from __future__ import annotations

import numpy as np

from . import sik as SIK

CHIP_RATE, DEV = 100e3, 50e3
SYNC = 0xF531FAB6
BANDS = ((865.5e6, 866.5e6), (867.9e6, 868.7e6))
ADDR_TYPES = {0x20: "FLARM ID", 0x10: "ICAO 24-bit address"}


def manchester(bits):
    """Data bits -> chips (1 -> 01, 0 -> 10)."""
    return np.array([c for b in bits for c in ((0, 1) if b else (1, 0))], np.uint8)


def _decode(chips):
    """Chips -> (data bits, number of invalid pairs 00 / 11)."""
    c = np.asarray(chips[:len(chips) // 2 * 2], np.int8).reshape(-1, 2)
    return (c[:, 1] > c[:, 0]).astype(np.uint8), int(np.count_nonzero(c[:, 0] == c[:, 1]))


def try_identify(y, ofs, chip_rate=CHIP_RATE, max_err=2):
    """One burst at ~8 samples / chip centred on the carrier -> dict (radio ID, address type) or None. Every timing
    phase, polarity and chip-pair alignment is tried; the hit with the fewest header / total Manchester violations
    (then sync errors) wins."""
    pat = np.array([(SYNC >> (31 - i)) & 1 for i in range(32)], np.int8) * 2 - 1
    best = None
    for ph, chips in SIK.demod_bits(y, ofs, chip_rate):
        for inv in (False, True):
            ch = 1 - chips if inv else chips
            for odd in (0, 1):                                 # chip-pair alignment
                bits, bad = _decode(ch[odd:])
                if len(bits) < 32 + 32:
                    continue
                s = bits.astype(np.int8) * 2 - 1
                c = np.correlate(s, pat, "valid")
                for p in np.flatnonzero(c >= 32 - 2 * max_err):
                    hb = bits[p + 32:p + 64]
                    if len(hb) < 32:
                        continue
                    by = [int("".join(map(str, hb[8 * i:8 * i + 8])), 2) for i in range(4)]
                    pairs = ch[odd + 2 * p:odd + 2 * (p + 64)].reshape(-1, 2)
                    viol = int(np.count_nonzero(pairs[:, 0] == pairs[:, 1]))
                    if viol > 4:                               # not Manchester: chance match
                        continue
                    rid = by[0] | by[1] << 8 | by[2] << 16
                    # no CRC verifies the header: a chip error inside a Manchester pair shows as 00 / 11, so the ID
                    # is reported only when its 32 header bits have no violation
                    hv = int(np.count_nonzero(pairs[32:, 0] == pairs[32:, 1]))
                    out = dict(protocol="FLARM-like (868 MHz beacon)", sync=f"0x{SYNC:08X}",
                               sync_bit_errors=int((32 - c[p]) // 2), manchester_violations=viol,
                               radio_id=f"{rid:06X}" if hv == 0 else None,
                               address_type=ADDR_TYPES.get(by[3], f"0x{by[3]:02X}") if hv == 0 else None,
                               polarity_inverted=inv, payload="obfuscated (not decoded)")
                    if hv:
                        out["note"] = f"{hv} Manchester violation(s) in the header: ID not reliable"
                    key = (hv, viol, out["sync_bit_errors"])
                    if best is None or key < best[0]:
                        best = (key, out)
    return None if best is None else best[1]
