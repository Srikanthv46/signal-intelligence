"""
toys.py -- toy-drone 2.4 GHz protocols on nRF24L01 / XN297 radios (research, PLAN G6;
research/drone_transmitters/14_market_gap_survey.md section 3). Synthetic validation only.

1 Mb/s GFSK at 2400 + n MHz. Addresses are given in register order (the air sends them last byte first).
- XN297-scrambled (the xn297 module demodulates address + payload; this module names the protocol from the address,
  the payload length and the payload's own checksum):
  MJXQ WLH08 / E010-class (address 6D 6A 77 77 77, 16 bytes, byte 15 = sum of 0-14, TX ID bytes 7-9),
  MJXQ X600 / X800 (6D 6A 73 73 73), CX10 family (CC x 5, 15 bytes, byte 0 = 0x55 data / 0xAA bind, TX ID 1-4),
  JJRC345 / SkyTumbler (CC x 5, 16 bytes, byte 0 = 0, 16-bit sum split over bytes 13 and 3, TX ID 14-15),
  H8_3D / CG023 / H20 (C4 57 09 65 21, 20 bytes, byte 8 = sum of the TX ID 1-4, byte 19 = sum of 9-18),
  E016H (data address A5 00 + TX ID, 10 bytes, byte 9 = sum of 0-8),
  Potensic A20 (data address = TX ID, 10 bytes, byte 0 = 0x64, low nibble of byte 8 = sum of 1-7).
- Plain nRF24 (classic ShockBurst, CRC-16 CCITT init 0xFFFF over address + payload; verified by the CRC):
  V2x2 (WLtoys V202 / V222 / V262, 66 88 68 68 68, 16 bytes, byte 15 = sum, TX ID 7-9),
  SymaX X5C (6D 6A 73 73 73, 16 bytes, byte 15 = sum), SymaX (bind AB AC AD AE AF / data A2 + TX ID, 10 bytes,
  byte 9 = XOR of 0-8 + 0x55; a 9-bit ESB control field may follow the address: both layouts are tried).
Not covered here: the 250 kb/s XN297 protocols (V911S, KF606, MJXQ E010 / Phoenix) and the non-nRF chips (E012 / E015
HS6200, E129 / E010 r5, E016H v2). Nothing is transmitted or controlled.
"""
from __future__ import annotations

import numpy as np

RATE = 1e6
_REV = [int(f"{v:08b}"[::-1], 2) for v in range(256)]


def _sum(b):
    return sum(b) & 0xFF


def _xn_rules():
    def mjxq(p):
        return len(p) >= 16 and _sum(p[:15]) == p[15]

    def cx10(p):
        return len(p) >= 15 and p[0] in (0x55, 0xAA)

    def jjrc(p):
        if len(p) < 16 or p[0] != 0:
            return False
        q = bytearray(p[:13])
        q[3] = 0
        s = 2 + sum(q)
        return p[13] == s & 0xFF and p[3] == ((((s >> 8) & 0xFF) << 2) + 2) & 0xFF

    def h83d(p):
        return len(p) >= 20 and p[0] in (0x13, 0x14) and _sum(p[1:5]) == p[8] and _sum(p[9:19]) == p[19]

    def e016h(p):
        return len(p) >= 10 and _sum(p[:9]) == p[9]

    def potensic(p):
        return len(p) >= 10 and p[0] == 0x64 and (sum(p[1:8]) & 0x0F) == (p[8] & 0x0F)

    return [  # (name, address register order or None = any, payload length, rule, TX ID slice)
        ("MJXQ (WLH08 / E010 class)", bytes.fromhex("6d6a777777"), 16, mjxq, (7, 10)),
        ("MJXQ (X600 / X800)", bytes.fromhex("6d6a737373"), 16, mjxq, (7, 10)),
        ("JJRC345 / SkyTumbler", bytes.fromhex("cccccccccc"), 16, jjrc, (14, 16)),
        ("CX10 family (CX10 / DM007 / JC3015 / MK33041)", bytes.fromhex("cccccccccc"), 15, cx10, (1, 5)),
        ("H8_3D / CG023 / H20", bytes.fromhex("c457096521"), 20, h83d, (1, 5)),
        ("H8_3D (H20H)", bytes.fromhex("eeddccbb11"), 20, h83d, (1, 5)),
        ("E016H", None, 10, lambda p: e016h(p), None),
        ("Potensic A20", None, 10, potensic, None),
    ]


XN_RULES = _xn_rules()


def match_xn297(addr: bytes, unscrambled: bytes):
    """XN297 address (register order) + de-scrambled, bit-reversed payload bytes -> dict or None. Protocols with a
    fixed address need it; ID-addressed ones (E016H, Potensic) need their stricter payload checks."""
    for name, a, n, rule, tid in XN_RULES:
        if a is not None and addr != a:
            continue
        p = unscrambled[:n]
        if len(p) < n or not rule(p):
            continue
        if a is None and name == "E016H" and addr[:2] != b"\xa5\x00":
            continue                                    # E016H data address: A5 00 + TX ID
        out = dict(protocol=f"{name} (XN297)", address=addr.hex(), payload_len=n, payload_hex=p.hex(),
                   checksum_ok=True)
        out["tx_id"] = p[tid[0]:tid[1]].hex() if tid else addr.hex()
        if name.startswith("CX10"):
            out.update(checksum_ok=None, packet_name="bind" if p[0] == 0xAA else "data",
                       note="no payload checksum: fixed address + header byte only")
        return out
    return None


# ---------------------------------------------------------------- plain nRF24 (classic ShockBurst)
def crc16_nrf(data: bytes, init=0xFFFF, bits=None) -> int:
    """nRF24L01 CRC-16 (CCITT 0x1021, MSB first) over whole bytes, or over a bit list when bits is given."""
    c = init
    seq = bits if bits is not None else [(v >> (7 - i)) & 1 for v in data for i in range(8)]
    for b in seq:
        c = ((c << 1) ^ 0x1021) & 0xFFFF if ((c >> 15) & 1) ^ b else (c << 1) & 0xFFFF
    return c


NRF_PROTOS = [  # (name, address register order or None, payload length, check, TX ID slice)
    ("V2x2 (WLtoys V202 / V222 / V262)", bytes.fromhex("6688686868"), 16, lambda p: _sum(p[:15]) == p[15], (7, 10)),
    ("SymaX X5C", bytes.fromhex("6d6a737373"), 16, lambda p: _sum(p[:15]) == p[15], None),
    ("SymaX (bind)", bytes.fromhex("abacadaeaf"), 10, lambda p: ((_xor(p[:9]) + 0x55) & 0xFF) == p[9], (0, 5)),
    ("SymaX", None, 10, lambda p: ((_xor(p[:9]) + 0x55) & 0xFF) == p[9], None),
]


def _xor(b):
    v = 0
    for x in b:
        v ^= x
    return v


def nrf_packet_bits(addr_reg: bytes, payload: bytes, pcf=None) -> np.ndarray:
    """Synthesis / tests: preamble, address (air order), [9-bit PCF], payload, CRC-16 as bits."""
    a = bytes(addr_reg[::-1])
    bits = [(v >> (7 - i)) & 1 for v in a for i in range(8)]
    if pcf is not None:
        bits += [(pcf >> (8 - i)) & 1 for i in range(9)]
    bits += [(v >> (7 - i)) & 1 for v in payload for i in range(8)]
    c = crc16_nrf(b"", bits=bits)
    bits += [(c >> (15 - i)) & 1 for i in range(16)]
    pre = [(i + bits[0]) % 2 for i in range(8)]              # 0x55 / 0xAA: alternating into the first address bit
    return np.array(pre + bits, np.uint8)


def _fullbit_streams(y, ofs, rate, phases=8):
    """Hard decisions on the phase change over a whole bit at every timing phase (as xn297 L1c: at h ~0.3 the
    half-bit discriminator misses an isolated bit), threshold midway between the two tones."""
    sps = ofs / rate
    ph = np.unwrap(np.angle(np.asarray(y, np.complex128)))
    out = []
    for p in np.arange(0, sps, max(sps / phases, 1)):
        n = int((len(y) - p) / sps) - 1
        if n < 1:
            continue
        t0 = np.round(p + np.arange(n) * sps).astype(int)
        t1 = np.minimum(np.round(p + (np.arange(n) + 1) * sps).astype(int), len(ph) - 1)
        d = ph[t1] - ph[t0]
        from .sik import two_level_threshold
        thr = two_level_threshold(d)
        out.append((float(p), (d > thr).astype(np.uint8)))
    return out


def _nrf_check(b, s0, name, n, check, tid, inv):
    """CRC (classic layout, then with a 9-bit ESB control field) and the protocol's checksum at address start s0."""
    air = bytes(int("".join(map(str, b[s0 + 8 * i:s0 + 8 * i + 8])), 2) for i in range(5))
    reg = air[::-1]
    for pcf in (0, 9):
        seg = b[s0:s0 + 40 + pcf + 8 * n + 16]
        if len(seg) < 40 + pcf + 8 * n + 16:
            continue
        body = seg[:40 + pcf + 8 * n]
        if crc16_nrf(b"", bits=list(body)) != int("".join(map(str, seg[-16:])), 2):
            continue
        pb = seg[40 + pcf:40 + pcf + 8 * n]
        p = bytes(int("".join(map(str, pb[8 * i:8 * i + 8])), 2) for i in range(n))
        if not check(p):
            continue
        return dict(protocol=f"{name} (nRF24)", address=reg.hex(), payload_len=n, payload_hex=p.hex(),
                    esb_control_field=bool(pcf), crc_ok=True, checksum_ok=True, polarity_inverted=bool(inv),
                    tx_id=p[tid[0]:tid[1]].hex() if tid else reg.hex())
    return None


def _pattern(addr_reg: bytes) -> np.ndarray:
    return np.array([(v >> (7 - i)) & 1 for v in addr_reg[::-1] for i in range(8)], np.int8) * 2 - 1


def nrf24_try_decode(y, ofs, rate=RATE):
    """One 1 Mb/s burst at ~8 samples/bit -> dict (CRC-verified plain nRF24 toy protocol) or None. Fixed addresses are
    found by correlation (<= 1 bit error); the ID-based SymaX data address (A2 + TX ID) after the preamble run."""
    from . import sik as SIK
    fixed = [(pr, _pattern(pr[1])) for pr in NRF_PROTOS if pr[1] is not None]
    for ph, bits in _fullbit_streams(y, ofs, rate) + SIK.demod_bits_mid(y, ofs, rate):
        for inv in (False, True):
            b = (1 - bits) if inv else bits
            if len(b) < 40:
                continue
            sb = b.astype(np.int8) * 2 - 1
            for (name, a, n, check, tid), pat in fixed:
                for s0 in np.flatnonzero(np.correlate(sb, pat, "valid") >= 38):
                    d = _nrf_check(b, int(s0), name, n, check, tid, inv)
                    if d:
                        return d
            name, _, n, check, tid = NRF_PROTOS[-1]                       # SymaX data: A2 + TX ID
            dd = np.flatnonzero(np.diff(b.astype(np.int8)) == 0)
            edges = np.r_[-1, dd, len(b) - 1]
            for a0, a1 in zip(edges[:-1], edges[1:]):
                if a1 - a0 < 6:
                    continue
                for s0 in range(max(int(a0) + 1, int(a1) - 12), int(a1) + 3):
                    if s0 < 0 or s0 + 40 > len(b):
                        continue
                    if int("".join(map(str, b[s0 + 32:s0 + 40])), 2) != 0xA2:     # register byte 0 is sent last
                        continue
                    d = _nrf_check(b, s0, name, n, check, tid, inv)
                    if d:
                        return d
    return None
