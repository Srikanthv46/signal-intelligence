"""
xn297.py -- XN297 toy-drone links and the Bayang protocol (research, PLAN R8;
research/drone_transmitters/10_cc2500_a7105_xn297_flarm_xbee.md section 4). Synthetic validation only.

PHY: GFSK 1 Mb/s at 2.4 GHz. On air after the preamble: the XN297 prefix 0x710F55, the address (sent in reverse byte
order, each byte XORed with the scramble table), the payload (each byte bit-reversed, then scrambled), CRC-16.
- Scramble table: a 7-bit Fibonacci LFSR (feedback = state bits 6 xor 2, output bit 6, MSB first) seeded 0x71. This
  generator was found independently and reproduces the published table start E3 B1 4B EA 85 BC E5 66 exactly.
- CRC-16 0x1021, init 0xB5D2 over the scrambled on-air bytes, XORed with a length-dependent constant whose table is
  published only partly (GPL source): not verified here. The residual (computed ^ received) is reported; it is the
  same for every packet of one length, which a real capture can confirm.
- Bayang (5-byte address, 15-byte payload) verifies itself: header byte A1-A6 / 0x53 and byte 14 = 8-bit sum of
  bytes 0-13. Bind packets carry the TX address and the 4 hop channels; data packets carry sticks and flags.
Demodulation: the full-bit phase change with the decision level from the known prefix first (PLAN L1c, ~2.5 dB over
the half-bit discriminator), then the discriminator slicer.
Nothing is transmitted or controlled.
"""
from __future__ import annotations

import numpy as np

from . import sik as SIK
from . import toys as TOYS

RATE = 1e6
PREFIX = 0x710F55
ADDR_LEN, BAYANG_LEN = 5, 15
MIN_PAYLOAD, MAX_PAYLOAD = 10, 21                 # toy protocols (toys.py): 10-21 payload bytes
BAYANG_BIND = {0xA1, 0xA2, 0xA3, 0xA4, 0x53}
BAYANG_DATA = {0xA5, 0xA6}


def scramble_table(n=40):
    st, out = 0x71, []
    for _ in range(n):
        b = 0
        for _ in range(8):
            bit = (st >> 6) & 1
            st = ((st << 1) | (((st >> 6) ^ (st >> 2)) & 1)) & 0x7F
            b = (b << 1) | bit
        out.append(b)
    return out


SCRAMBLE = scramble_table()
_REV = [int(f"{v:08b}"[::-1], 2) for v in range(256)]


def crc16(data: bytes, init=0xB5D2) -> int:
    c = init
    for b in data:
        c ^= b << 8
        for _ in range(8):
            c = ((c << 1) ^ 0x1021) & 0xFFFF if c & 0x8000 else (c << 1) & 0xFFFF
    return c


def encode(addr: bytes, payload: bytes, crc_xor=0) -> bytes:
    """On-air bytes after the prefix (synthesis / tests): scrambled reversed address, scrambled bit-reversed payload,
    CRC (with the given length constant)."""
    a = bytes(addr[len(addr) - 1 - j] ^ SCRAMBLE[j] for j in range(len(addr)))
    p = bytes(_REV[v] ^ SCRAMBLE[len(addr) + i] for i, v in enumerate(payload))
    c = crc16(a + p) ^ crc_xor
    return a + p + bytes([c >> 8, c & 0xFF])


def _find_prefix(bits, max_err=1):
    pat = np.array([(PREFIX >> (23 - i)) & 1 for i in range(24)], np.int8) * 2 - 1
    s = bits.astype(np.int8) * 2 - 1
    if len(s) < 24 + 8 * (ADDR_LEN + MIN_PAYLOAD + 2):
        return []
    c = np.correlate(s, pat, "valid")
    hits = [(int(p) + 24, False, int((24 - c[p]) // 2)) for p in np.flatnonzero(c >= 24 - 2 * max_err)]
    hits += [(int(p) + 24, True, int((24 + c[p]) // 2)) for p in np.flatnonzero(-c >= 24 - 2 * max_err)]
    return hits


def _bytes(bits):
    n = len(bits) // 8
    return bytes(int("".join(map(str, bits[8 * i:8 * i + 8])), 2) for i in range(n))


def _bayang(p):
    t = p[0]
    out = dict(packet_type=f"0x{t:02X}")
    if t in BAYANG_BIND:
        out.update(packet_name="bind", tx_id=p[1:6].hex(), hop_channels=list(p[6:10]))
    else:
        st = [((p[4 + 2 * k] & 0x03) << 8) | p[5 + 2 * k] for k in range(4)]
        out.update(packet_name="data", sticks=dict(zip(("aileron", "elevator", "throttle", "rudder"), st)),
                   flags=[f"0x{p[2]:02X}", f"0x{p[3]:02X}"], tx_id_byte2=f"0x{p[12]:02X}")
    return out


def _soft_streams(y, ofs, rate):
    """Full-bit phase change at every timing phase (PLAN L1c): for CPFSK at h ~ 0.3 the phase change over a whole bit
    is a better statistic than the half-bit-smoothed discriminator (synthetic: ~2.5 dB). Each stream is centred on
    its median; the prefix fixes the decision level later. Returns [(phase, soft)]."""
    sps = ofs / rate
    ph = np.unwrap(np.angle(np.asarray(y, np.complex128)))
    out = []
    for p in np.arange(0, sps, max(sps / 8, 1)):
        n = int((len(y) - p) / sps) - 1
        if n < 1:
            continue
        t0 = np.round(p + np.arange(n) * sps).astype(int)
        t1 = np.minimum(np.round(p + (np.arange(n) + 1) * sps).astype(int), len(ph) - 1)
        d = ph[t1] - ph[t0]
        out.append((float(p), d - np.median(d)))
    return out


_PAT = np.array([(PREFIX >> (23 - i)) & 1 for i in range(24)], np.int8) * 2 - 1


def _packets(y, ofs, rate):
    """Candidate on-air bit streams after the prefix: [(bits after the prefix, inverted, prefix errors)]; full-bit
    soft decisions with the level from the 24 known prefix bits first, then the hard discriminator slicer."""
    nb, nb_min = 8 * (ADDR_LEN + MAX_PAYLOAD + 2), 8 * (ADDR_LEN + MIN_PAYLOAD + 2)
    for ph, s in _soft_streams(y, ofs, rate):
        if len(s) < 24 + nb_min:
            continue
        c = np.correlate(np.sign(s).astype(np.int8), _PAT, "valid")
        for p in np.flatnonzero(np.abs(c) >= 22):
            inv = bool(c[p] < 0)
            v = s * (-1 if inv else 1)
            seg = v[p:p + 24]
            thr = 0.5 * (seg[_PAT > 0].mean() + seg[_PAT < 0].mean())
            b = v[p + 24:p + 24 + nb] - thr
            if len(b) >= nb_min:
                yield (b > 0).astype(np.uint8), inv, int((24 - abs(int(c[p]))) // 2)
    for ph, bits in SIK.demod_bits(y, ofs, rate):
        for pos, inv, err in _find_prefix(bits):
            b = bits[pos:pos + nb]
            if len(b) >= nb_min:
                yield (1 - b if inv else b), inv, err


def try_decode(y, ofs, rate=RATE):
    """One 2.4 GHz burst at ~8 samples/bit centred on the carrier -> dict or None. A Bayang packet whose header and
    checksum verify wins; otherwise an XN297 prefix with no bit error gives the 'XN297-like' address."""
    like = None
    for b, inv, err in _packets(y, ofs, rate):
        air = _bytes(b)
        addr = bytes(air[ADDR_LEN - 1 - j] ^ SCRAMBLE[ADDR_LEN - 1 - j] for j in range(ADDR_LEN))
        unscr = bytes(_REV[air[ADDR_LEN + i] ^ SCRAMBLE[ADDR_LEN + i]] for i in range(min(MAX_PAYLOAD, len(air) - ADDR_LEN)))
        p = unscr[:BAYANG_LEN]
        if len(air) >= ADDR_LEN + BAYANG_LEN + 2:
            rx = air[ADDR_LEN + BAYANG_LEN] << 8 | air[ADDR_LEN + BAYANG_LEN + 1]
            resid = f"0x{crc16(air[:ADDR_LEN + BAYANG_LEN]) ^ rx:04X}"
        else:
            resid = None
        if (len(p) == BAYANG_LEN and (p[0] in BAYANG_BIND or p[0] in BAYANG_DATA)
                and p[14] == sum(p[:14]) & 0xFF and resid is not None):
            return dict(protocol="Bayang (XN297)", prefix=f"0x{PREFIX:06X}", prefix_bit_errors=err,
                        polarity_inverted=bool(inv), address=addr.hex(), checksum_ok=True, crc_ok=None,
                        crc_residual=resid, payload_hex=p.hex(), **_bayang(p))
        toy = TOYS.match_xn297(addr, unscr)                  # other toy protocols (research, PLAN G6)
        if toy:
            toy.update(prefix=f"0x{PREFIX:06X}", prefix_bit_errors=err, polarity_inverted=bool(inv))
            return toy
        if like is None and err == 0:
            like = dict(protocol="XN297-like (toy drone)", prefix=f"0x{PREFIX:06X}", prefix_bit_errors=0,
                        polarity_inverted=bool(inv), address=addr.hex(), checksum_ok=None, crc_ok=None,
                        note="5-byte address assumed; payload protocol unknown")
    return like
