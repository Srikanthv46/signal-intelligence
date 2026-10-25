"""
cc2500.py -- passive identification of FrSky D8 / D16 (ACCST) drone control links (research, PLAN R6;
research/drone_transmitters/10_cc2500_a7105_xn297_flarm_xbee.md). Synthetic validation only.

A received GFSK burst in the 2.4 GHz band is checked for the CC2500 sync word after the alternating preamble:
- D8: sync 0xD391, length 0x11, CC2500 hardware CRC-16 (0x8005, init 0xFFFF, MSB first) over length + payload.
- D16: sync 0xD391D391, length 0x1D (FCC v1 / v2) or 0x20 (LBT v1); a software CRC-16 in the last two packet bytes
  (reflected-0x1021 table applied MSB first, seed 0) over bytes 3 .. L-3; v2 appends the hardware CRC after it.
A CRC pass confirms the link; the transmitter ID (bytes 1-2, in the clear) identifies it. An exact 32-bit D16 sync with
an unexpected length byte (e.g. the receiver's telemetry frame, whose length is unverified) is reported as a sync-only
match. Nothing is transmitted or controlled.
"""
from __future__ import annotations

import numpy as np

from . import sik as SIK

D8_SYNC = [0xD3, 0x91]
D16_SYNC = [0xD3, 0x91, 0xD3, 0x91]
D8_LEN = 0x11
D16_LENS = {0x1D: "FCC", 0x20: "LBT v1"}
RATES = {"D8": 31.04e3, "D16 v1": 70.0e3, "D16 v2": 76.97e3, "D16 LBT": 99.98e3}
_T = []
for _v in range(256):
    _c = _v
    for _ in range(8):
        _c = (_c >> 1) ^ 0x8408 if _c & 1 else _c >> 1
    _T.append(_c)                                       # KERMIT table: T[1] = 0x1189


def crc16_cc2500(data: bytes) -> int:
    """CC2500 hardware CRC-16 (0x8005, init 0xFFFF, MSB first), over the given bytes (DN502)."""
    c = 0xFFFF
    for b in data:
        c ^= b << 8
        for _ in range(8):
            c = ((c << 1) ^ 0x8005) & 0xFFFF if c & 0x8000 else (c << 1) & 0xFFFF
    return c


def crc16_d16(data: bytes) -> int:
    """FrSky D16 software CRC: KERMIT table applied MSB first, seed 0."""
    c = 0
    for b in data:
        c = ((c << 8) & 0xFFFF) ^ _T[(c >> 8) ^ b]
    return c


def d16_packet(payload: bytes, v2=False) -> bytes:
    """Length byte + payload (bytes 1 .. L-3) + software CRC (+ v2 hardware CRC): as the CC2500 sends after the sync."""
    p = bytes([len(payload) + 2]) + bytes(payload)
    c = crc16_d16(p[3:])
    p += bytes([c >> 8, c & 0xFF])
    if v2:
        h = crc16_cc2500(p)
        p += bytes([h >> 8, h & 0xFF])
    return p


def d8_packet(payload: bytes) -> bytes:
    """Length byte + payload + hardware CRC."""
    p = bytes([len(payload)]) + bytes(payload)
    h = crc16_cc2500(p)
    return p + bytes([h >> 8, h & 0xFF])


def _find_sync(bits: np.ndarray, sync_bytes):
    """[(pos just after the sync, inverted, sync bit errors)] where sync_bytes (MSB first) match with <= 2 bit errors,
    either polarity, after >= 6 of 7 alternating preamble transitions."""
    pat = np.array([(v >> (7 - i)) & 1 for v in sync_bytes for i in range(8)], np.int8) * 2 - 1
    s = bits.astype(np.int8) * 2 - 1
    L = len(pat)
    if len(s) < L + 8:
        return []
    c = np.correlate(s, pat, "valid")
    hits = []
    for inv, cc in ((False, c), (True, -c)):
        for p in np.flatnonzero(cc >= L - 4):
            if p >= 8 and int(np.sum(np.abs(np.diff(s[p - 8:p] * (-1 if inv else 1))))) >= 2 * 6:
                hits.append((int(p) + L, inv, int((L - cc[p]) // 2)))
    return hits


def _bytes(bits) -> bytes:
    n = len(bits) // 8
    return bytes(int("".join(map(str, bits[8 * i:8 * i + 8])), 2) for i in range(n))      # MSB first


def _check_d8(raw):
    if len(raw) < D8_LEN + 3 or raw[0] != D8_LEN:
        return None
    n = D8_LEN + 1
    ok = crc16_cc2500(raw[:n]) == (raw[n] << 8 | raw[n + 1])
    return dict(variant="D8", crc_ok=bool(ok), crc="hardware 0x8005" if ok else None, packet=raw[:n + 2])


def _check_d16(raw):
    if not raw or raw[0] not in D16_LENS:
        return None
    L = raw[0] + 1                                      # bytes incl. the length byte; software CRC in the last two
    if len(raw) < L:
        return None
    ok = crc16_d16(raw[3:L - 2]) == (raw[L - 2] << 8 | raw[L - 1])
    hw = len(raw) >= L + 2 and crc16_cc2500(raw[:L]) == (raw[L] << 8 | raw[L + 1])
    var = "D16 v2" if ok and hw else ("D16 v1" if raw[0] == 0x1D else "D16 LBT v1") if ok else "D16"
    return dict(variant=var, crc_ok=bool(ok), crc=("software + hardware" if hw else "software") if ok else None,
                packet=raw[:L + 2] if hw else raw[:L])


def try_identify(y, ofs, rate):
    """One 2.4 GHz GFSK burst at ~8 samples/bit centred on the carrier -> dict (variant, TX ID, CRC state) or None.
    Every timing phase and sync hit is tried: the first CRC pass wins, else the hit with the fewest sync bit errors."""
    best = None
    for name, sync, check in (("D16", D16_SYNC, _check_d16), ("D8", D8_SYNC, _check_d8)):
        sync_only = None
        for ph, bits in SIK.demod_bits(y, ofs, rate):
            for pos, inv, err in _find_sync(bits, sync):
                b = bits[pos:]
                raw = _bytes(1 - b if inv else b)
                r = check(raw)
                if r is None:
                    if name == "D16" and err == 0 and len(raw) >= 1 and sync_only is None:
                        sync_only = dict(protocol="FrSky D16-like (ACCST)", variant="D16", sync="0xD391D391",
                                         length=raw[0], tx_id=None, crc_ok=None, crc=None,
                                         packet_hex=raw[:min(len(raw), raw[0] + 3)].hex(),
                                         polarity_inverted=bool(inv), note="32-bit sync only; length not a data frame")
                    continue
                out = dict(protocol=f"FrSky {r['variant']} (ACCST)" if r["crc_ok"] else f"FrSky {name}-like (ACCST)",
                           variant=r["variant"], sync=f"0x{''.join(f'{v:02X}' for v in sync)}", length=raw[0],
                           tx_id=f"0x{raw[2]:02X}{raw[1]:02X}", crc_ok=r["crc_ok"], crc=r["crc"],
                           packet_hex=r["packet"].hex(), polarity_inverted=bool(inv))
                if r["crc_ok"]:
                    return out
                if best is None or err < best[0]:
                    best = (err, out)
        if best is not None:
            return best[1]                                      # D16 sync contains D8's: do not fall through
        if sync_only is not None:
            return sync_only
    return None
