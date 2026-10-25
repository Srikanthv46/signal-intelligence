"""TPMS frames on a PCM symbol stream (research, signal gaps S6; research/signal_gaps/04_short_range.md section 1.4):
Continental / VDO (Hyundai-Kia family) and Toyota / Pacific PMV-C210. Bit conventions follow the rtl_433 decoders
(used as a fact source only): VDO = inverted stream, preamble aa aa aa a9 (32 symbols), Manchester (the second
symbol of each pair), 80 bits, CRC-8 poly 0x07 init 0xaa over bytes 0-8; Toyota = preamble a9 e (12 symbols),
differential Manchester from its last symbol (no mid-bit transition = 1), 72 bits, CRC-8 poly 0x07 init 0x80 over
bytes 0-7 and the inverted pressure copy. Pure Python; mirrored in C++ receivers/ism_devices/tpms.
"""
from __future__ import annotations

VDO_PRE = [int(c) for c in "10101010101010101010101010101001"]
TOYOTA_PRE = [int(c) for c in "101010011110"]


def crc8(data: bytes, poly: int, init: int) -> int:
    c = init
    for b in data:
        c ^= b
        for _ in range(8):
            c = ((c << 1) ^ poly) & 0xFF if c & 0x80 else (c << 1) & 0xFF
    return c


def _find(s: list[int], pat: list[int], start: int) -> int:
    n, m = len(s), len(pat)
    for i in range(start, n - m + 1):
        if s[i:i + m] == pat:
            return i
    return -1


def _bytes(bits: list[int]) -> bytes:
    return bytes(sum(bits[k + q] << (7 - q) for q in range(8)) for k in range(0, len(bits) - 7, 8))


def manchester(s: list[int], start: int, max_bits: int) -> list[int]:
    out = []
    i = start
    end = min(len(s), start + 2 * max_bits)
    while i + 1 < end:
        if s[i] == s[i + 1]:
            break
        out.append(s[i + 1])
        i += 2
    return out


def diff_manchester(s: list[int], start: int, max_bits: int) -> tuple[list[int], int]:
    """(bits, symbols consumed): the first long pulse sets the clock; a missing boundary transition ends it."""
    out = []
    i = start
    end = min(len(s), start + 2 * max_bits)
    b2 = 0
    while i + 2 < end:
        b1, b2, b3 = s[i], s[i + 1], s[i + 2]
        i += 2
        if b1 != b2:
            if b2 != b3:
                out.append(0)
            else:
                b2 = b1
                i -= 1
                break
        else:
            b2 = 1 - b1
            i -= 2
            break
    while i + 1 < end:
        b1 = s[i]
        i += 1
        if b1 == b2:
            break
        b2 = s[i]
        i += 1
        out.append(1 if b1 == b2 else 0)
    return out, i - start


def vdo(s: list[int]) -> dict | None:
    inv = [1 - v for v in s]
    pos = 0
    while True:
        pos = _find(inv, VDO_PRE, pos)
        if pos < 0 or pos + 32 + 160 > len(inv):
            return None
        bits = manchester(inv, pos + 32, 80)
        if len(bits) >= 80:
            b = _bytes(bits[:80])
            if crc8(b[:9], 0x07, 0xAA) == b[9]:
                return {"protocol": "Continental / VDO", "id": f"{(b[1] << 24) | (b[2] << 16) | (b[3] << 8) | b[4]:08x}",
                        "state": b[0], "flags": b[5] >> 4, "repeat": b[5] & 0x0F,
                        "pressure_kpa": round(b[6] * 1.375, 3), "temperature_c": b[7] - 50, "battery_raw": b[8],
                        "crc_ok": True, "frame_hex": b.hex()}
        pos += 2


def toyota(s: list[int]) -> dict | None:
    pos = 0
    while True:
        pos = _find(s, TOYOTA_PRE, pos)
        if pos < 0 or pos + 156 > len(s):
            return None
        bits, used = diff_manchester(s, pos + 11, 80)
        if used >= 144 and len(bits) >= 72:
            b = _bytes(bits[:72])
            p1 = ((b[4] & 0x7F) << 1) | (b[5] >> 7)
            if crc8(b[:8], 0x07, 0x80) == b[8] and p1 == b[7] ^ 0xFF:
                temp = ((b[5] & 0x7F) << 1) | (b[6] >> 7)
                return {"protocol": "Toyota / Pacific PMV-C210",
                        "id": f"{(b[0] << 24) | (b[1] << 16) | (b[2] << 8) | b[3]:08x}",
                        "status": (b[4] & 0x80) | (b[6] & 0x7F), "pressure_psi": p1 * 0.25 - 7.0,
                        "pressure_kpa": round((p1 * 0.25 - 7.0) * 6.894757, 3),
                        "temperature_c": temp - 40, "crc_ok": True, "frame_hex": b.hex()}
        pos += 1


DECODERS = (vdo, toyota)


def decode(s: list[int]) -> dict | None:
    """The first CRC-valid TPMS frame in the symbol stream, either FSK polarity."""
    for pol in (0, 1):
        t = s if pol == 0 else [1 - v for v in s]
        for dec in DECODERS:
            r = dec(t)
            if r is not None:
                return r
    return None
