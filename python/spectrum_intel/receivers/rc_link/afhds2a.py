"""
afhds2a.py -- FlySky AFHDS 2A (A7105) passive identification (research, PLAN R7;
research/drone_transmitters/10_cc2500_a7105_xn297_flarm_xbee.md section 3). Synthetic validation only.

PHY: 2-FSK 500 kb/s, +-186 kHz, 2400.5-2482 MHz on a 500 kHz grid, 4-byte preamble, 32-bit ID code 0x5475C52A
(MSB first; the A7105 accepts <= 3 bit errors), then 38 payload bytes + CRC-16, every nibble (7,4) Hamming coded.
- Detection: the ID code after the alternating preamble (<= 3 bit errors, either polarity) -> "AFHDS 2A-like".
- Decode: the datasheet leaves the Hamming bit layout and the CRC polynomial / init unverified, so a small set of
  candidate conventions is tried (FEC_LAYOUTS x nibble order x CRC_VARIANTS); a decode is reported only when the CRC
  passes and the packet type is a known one, and it names the convention that passed. Fields in clear: packet type,
  TX ID (bytes 1-4), RX ID (bytes 5-8; FF FF FF FF while binding).
FlySky AFHDS 1 (FS-CT6B, FS-T6, 9X modules; research, PLAN G5) shares the PHY and the ID code: 21 payload bytes, no
FEC: 0x55 data / 0xAA bind, TX ID (4 bytes, little-endian), 8 channels (uint16 LE, 860-2140 us); A7105 CRC-16 (CCITT,
init unverified). 16 hop channels from TX ID byte 0 (an 8 x 4 nibble table); f = 2400.001 + 0.5 ch MHz; 1510 us.
It is reported when the header and all 8 channel values are plausible and the measured channel is in the ID's hop
table, or when a CRC convention passes.
Nothing is transmitted or controlled.
"""
from __future__ import annotations

import numpy as np

from . import sik as SIK

RATE, DEV = 500e3, 186e3
ID_CODE = 0x5475C52A
ID_MAX_ERR = 3
PAYLOAD = 38
TYPES = {0x58: "sticks", 0x56: "failsafe", 0xAA: "settings / telemetry", 0xAC: "telemetry", 0xBB: "bind",
         0xBC: "bind"}
# (7,4) Hamming layouts, codeword bits sent MSB first; d3 is the nibble's MSB
FEC_LAYOUTS = ("systematic", "classic")
CRC_VARIANTS = ((0x1021, 0x0000), (0x1021, 0xFFFF), (0x8005, 0x0000), (0x8005, 0xFFFF))


def _encode7(nib: int, layout: str) -> int:
    d3, d2, d1, d0 = (nib >> 3) & 1, (nib >> 2) & 1, (nib >> 1) & 1, nib & 1
    if layout == "systematic":           # d3 d2 d1 d0 p0 p1 p2
        bits = (d3, d2, d1, d0, d3 ^ d2 ^ d0, d3 ^ d1 ^ d0, d2 ^ d1 ^ d0)
    else:                                # classic positions 1..7: p1 p2 d3 p4 d2 d1 d0
        bits = (d3 ^ d2 ^ d0, d3 ^ d1 ^ d0, d3, d2 ^ d1 ^ d0, d2, d1, d0)
    v = 0
    for b in bits:
        v = (v << 1) | b
    return v


_ENC = {lay: [_encode7(n, lay) for n in range(16)] for lay in FEC_LAYOUTS}
_DEC = {}
for _lay in FEC_LAYOUTS:                 # nearest codeword (Hamming(7,4) is perfect: one per 7-bit word)
    _d = np.zeros(128, np.uint8)
    for _w in range(128):
        _d[_w] = min(range(16), key=lambda n: bin(_ENC[_lay][n] ^ _w).count("1"))
    _DEC[_lay] = _d


def crc16(data: bytes, poly: int, init: int) -> int:
    """MSB-first CRC-16, no reflection, no final XOR."""
    c = init
    for b in data:
        c ^= b << 8
        for _ in range(8):
            c = ((c << 1) ^ poly) & 0xFFFF if c & 0x8000 else (c << 1) & 0xFFFF
    return c


def encode_packet(payload: bytes, layout="systematic", high_first=True, crc=(0x1021, 0x0000)) -> np.ndarray:
    """payload (38 bytes) + CRC-16 (big-endian) -> FEC-coded bits after the ID code (synthesis / tests)."""
    c = crc16(payload, *crc)
    out = []
    for b in bytes(payload) + bytes([c >> 8, c & 0xFF]):
        for nib in ((b >> 4, b & 15) if high_first else (b & 15, b >> 4)):
            w = _ENC[layout][nib]
            out += [(w >> (6 - i)) & 1 for i in range(7)]
    return np.array(out, np.uint8)


AFHDS1_LEN = 21
AFHDS1_HEADERS = {0x55: "data", 0xAA: "bind"}
AFHDS1_TABLE = ((0x12, 0x34, 0x56, 0x78), (0x18, 0x27, 0x36, 0x45), (0x41, 0x82, 0x36, 0x57), (0x84, 0x13, 0x65, 0x72),
                (0x87, 0x64, 0x15, 0x32), (0x76, 0x84, 0x13, 0x52), (0x71, 0x62, 0x84, 0x35), (0x71, 0x86, 0x43, 0x52))
A7105_F0, A7105_STEP = 2400.000977e6, 0.5e6


def afhds1_hops(id3: int):
    """The 16 hop channels of an AFHDS 1 transmitter from TX ID byte 0 as sent (packet byte 1, already clamped)."""
    row, off = id3 & 0x0F, id3 >> 4
    hf = [0] * 16
    for i in range(16):
        b = AFHDS1_TABLE[row >> 1][i >> 2]
        n = b & 0x0F if i & 2 else b >> 4
        t = 10 * n + (80 if i & 1 else 0)
        hf[15 - i if row & 1 else i] = t - off
    return hf


def afhds1_packet(tx_id: int, channels_us, bind=False) -> bytes:
    """Synthesis / tests: the 21-byte AFHDS 1 payload (header, TX ID LE with the hop clamp, 8 channels LE)."""
    t = tx_id.to_bytes(4, "little")
    b0 = t[0] - 0x70 if (t[0] & 0xF0) > 0x90 else t[0]
    p = bytes([0xAA if bind else 0x55, b0, t[1], t[2], t[3]])
    return p + b"".join(int(v).to_bytes(2, "little") for v in list(channels_us)[:8])


def _afhds1(raw, f_hz):
    """AFHDS 1 fields of 23 raw bytes after the ID code (payload + CRC), or None."""
    p = raw[:AFHDS1_LEN]
    if len(raw) < AFHDS1_LEN + 2 or p[0] not in AFHDS1_HEADERS:
        return None
    ch = [p[5 + 2 * i] | p[6 + 2 * i] << 8 for i in range(8)]
    if not all(800 <= v <= 2200 for v in ch):
        return None
    hops = afhds1_hops(p[1])
    out = dict(protocol="FlySky AFHDS 1 (A7105)", packet_type=f"0x{p[0]:02X}", packet_name=AFHDS1_HEADERS[p[0]],
               tx_id=p[1:5][::-1].hex(), channels_us=ch, hop_channels=hops, packet_hex=p.hex(), crc_ok=None)
    for poly, init in CRC_VARIANTS:
        if crc16(p, poly, init) == (raw[AFHDS1_LEN] << 8 | raw[AFHDS1_LEN + 1]):
            out.update(crc_ok=True, convention=dict(crc_poly=f"0x{poly:04X}", crc_init=f"0x{init:04X}"))
            break
    if f_hz is not None:
        n = int(round((f_hz - A7105_F0) / A7105_STEP))
        out["channel"] = n
        out["hop_check"] = "match" if n in hops or (p[0] == 0xAA and n == 1) else "mismatch"
    if not out["crc_ok"] and out.get("hop_check") != "match":
        return None
    if not out["crc_ok"]:
        out["protocol"] = "FlySky AFHDS 1-like (A7105; header, channels and hop table verified, CRC convention not)"
    return out


def _find_id(bits):
    """[(pos after the ID code, inverted, bit errors)] with <= ID_MAX_ERR errors after >= 6 of 7 alternating bits."""
    pat = np.array([(ID_CODE >> (31 - i)) & 1 for i in range(32)], np.int8) * 2 - 1
    s = bits.astype(np.int8) * 2 - 1
    if len(s) < 40:
        return []
    c = np.correlate(s, pat, "valid")
    hits = []
    for inv, cc in ((False, c), (True, -c)):
        for p in np.flatnonzero(cc >= 32 - 2 * ID_MAX_ERR):
            if p >= 8 and int(np.sum(np.abs(np.diff(s[p - 8:p] * (-1 if inv else 1))))) >= 2 * 6:
                hits.append((int(p) + 32, inv, int((32 - cc[p]) // 2)))
    return hits


def _bytes(bits) -> bytes:
    n = len(bits) // 8
    return bytes(int("".join(map(str, bits[8 * i:8 * i + 8])), 2) for i in range(n))      # MSB first


def _fec_decode(cw_bits, layout, high_first):
    """7-bit codewords -> bytes and the number of corrected codewords."""
    n = len(cw_bits) // 7
    w = cw_bits[:7 * n].reshape(n, 7).astype(np.int64) @ (1 << np.arange(6, -1, -1))
    nib = _DEC[layout][w]
    fixes = int(np.sum(np.array(_ENC[layout])[nib] != w))
    nib = nib[:n - n % 2].reshape(-1, 2)
    hi, lo = (nib[:, 0], nib[:, 1]) if high_first else (nib[:, 1], nib[:, 0])
    return bytes((hi.astype(np.uint8) << 4) | lo.astype(np.uint8)), fixes


def try_identify(y, ofs, rate=RATE, f_hz=None):
    """One burst at ~8 samples/bit centred on the FSK carrier -> dict or None. f_hz: the burst's carrier (AFHDS 1 hop
    table check)."""
    nbits = (PAYLOAD + 2) * 2 * 7
    best = None
    for ph, bits in SIK.demod_bits(y, ofs, rate):
        for pos, inv, err in _find_id(bits):
            b1 = bits[pos:pos + 8 * (AFHDS1_LEN + 2)]
            if len(b1) == 8 * (AFHDS1_LEN + 2):
                d1 = _afhds1(_bytes(1 - b1 if inv else b1), f_hz)
                if d1:
                    d1.update(id_code=f"0x{ID_CODE:08X}", id_bit_errors=err, polarity_inverted=bool(inv))
                    return d1
            b = bits[pos:pos + nbits]
            if inv:
                b = 1 - b
            if len(b) == nbits:
                for layout in FEC_LAYOUTS:
                    for high_first in (True, False):
                        raw, fixes = _fec_decode(b, layout, high_first)
                        if raw[0] not in TYPES:
                            continue
                        for poly, init in CRC_VARIANTS:
                            if crc16(raw[:PAYLOAD], poly, init) == (raw[PAYLOAD] << 8 | raw[PAYLOAD + 1]):
                                return dict(protocol="FlySky AFHDS 2A (A7105)", id_code=f"0x{ID_CODE:08X}",
                                            id_bit_errors=err, polarity_inverted=bool(inv), crc_ok=True,
                                            packet_type=f"0x{raw[0]:02X}", packet_name=TYPES[raw[0]],
                                            tx_id=raw[1:5].hex(), rx_id=raw[5:9].hex(), binding=raw[5:9] == b"\xff" * 4,
                                            fec_corrections=fixes, packet_hex=raw[:PAYLOAD].hex(),
                                            convention=dict(fec=layout, nibble="high first" if high_first else
                                                            "low first", crc_poly=f"0x{poly:04X}",
                                                            crc_init=f"0x{init:04X}"))
            if best is None or err < best["id_bit_errors"]:
                best = dict(protocol="FlySky AFHDS 2A-like (A7105 ID code)", id_code=f"0x{ID_CODE:08X}",
                            id_bit_errors=err, polarity_inverted=bool(inv), crc_ok=None,
                            note="ID code only: FEC / CRC convention not verified (research)")
    return best
