"""
crossfire.py -- TBS Crossfire 150 Hz FSK air packets (research, PLAN R9; research/drone_transmitters/08_crossfire_dsm.md,
from g3gg0's analysis and MIT-licensed ESP32 sniffer). Synthetic validation only.

PHY (SX1272/6 FSK packet mode): 2-FSK 85.1 kb/s, deviation +-42.48 kHz (h ~ 1), preamble 0xAA..., sync word
0x3E 0x6B 0x3E (whether every unit uses it is unverified), whitening on (SX127x CCITT PN9, seed 0x1FF), fixed length:
uplink 23 bytes (CRC-16 0x1021 reflected over bytes 0-20 in bytes 21-22, per-hop init unknown), downlink 13 bytes
(CRC-8 0x07 over 0-11, per-hop init unknown). Grid: 860.165 / 902.165 MHz + n x 260 kHz, uplink n 0-49, downlink n+50;
CE 863.093 + n x 112 kHz.
Without the per-hop CRC init the packet cannot be verified: the 24-bit sync word after an alternating preamble is the
evidence; the packet bytes (raw and de-whitened) and the grid channel are reported.
"""
from __future__ import annotations

import numpy as np

from . import sik as SIK

RATE, DEV = 85.1e3, 42.48e3
SYNC = 0x3E6B3E
GRIDS = (("868", 860.165e6, 260e3, 100), ("915", 902.165e6, 260e3, 100), ("868 CE", 863.093e6, 112e3, 50))


def grid_channel(f_hz, tol_hz=20e3):
    """(band, n, link) of the nearest Crossfire grid channel within tol_hz, or None."""
    best = None
    for band, f0, step, nmax in GRIDS:
        n = int(round((f_hz - f0) / step))
        if 0 <= n < nmax:
            err = abs(f_hz - (f0 + n * step))
            if err <= tol_hz and (best is None or err < best[0]):
                link = ("uplink" if n < 50 else "downlink") if band != "868 CE" else "unknown"
                best = (err, band, n if n < 50 or band == "868 CE" else n - 50, link)
    return None if best is None else dict(band=best[1], channel=best[2], link=best[3])


def find_sync(bits, max_err=2):
    """Positions after the 24-bit sync preceded by >= 8 alternating bits, either polarity: [(pos, inverted)]."""
    pat = np.array([(SYNC >> (23 - i)) & 1 for i in range(24)], np.int8) * 2 - 1
    s = bits.astype(np.int8) * 2 - 1
    if len(s) < 48:
        return []
    c = np.correlate(s, pat, "valid")
    hits = []
    for inv, cc in ((False, c), (True, -c)):
        for p in np.flatnonzero(cc >= 24 - 2 * max_err):
            if p >= 8:
                q = s[p - 8:p] * (-1 if inv else 1)
                if int(np.sum(np.abs(np.diff(q)))) >= 2 * 6:          # >= 6 of 7 transitions: 1010 preamble
                    hits.append((int(p) + 24, inv))
    return hits


def _bytes(bits):
    n = len(bits) // 8
    return bytes(int("".join(map(str, bits[8 * i:8 * i + 8])), 2) for i in range(n))      # MSB first


def try_decode(y, ofs, f_hz):
    """One burst at ~8 samples / bit centred on the FSK carrier -> dict (sync, polarity, link, packet bytes raw and
    de-whitened) or None."""
    ch = grid_channel(f_hz)
    nbytes = 13 if (ch or {}).get("link") == "downlink" else 23
    for ph, bits in SIK.demod_bits(y, ofs, RATE):
        for pos, inv in find_sync(bits):
            b = bits[pos:pos + 8 * nbytes]
            if len(b) < 8 * nbytes:
                continue
            if inv:
                b = 1 - b
            raw = _bytes(b)
            dw = bytes(v ^ w for v, w in zip(raw, SIK.WHITEN["pn9_msb_first"][:nbytes]))
            out = dict(sync=f"0x{SYNC:06X}", polarity_inverted=bool(inv), packet_bytes=nbytes, packet_hex=raw.hex(),
                       dewhitened_hex=dw.hex(), crc_ok=None,
                       crc_note="per-hop CRC init unknown: packet not verifiable (research)")
            if ch:
                out.update(grid=ch)
            return out
    return None
