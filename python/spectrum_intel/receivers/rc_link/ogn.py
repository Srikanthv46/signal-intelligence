"""
ogn.py -- OGN tracker radio protocol (OGNTP v1) on the FLARM PHY: detection, clear header and de-whitened position
(research, PLAN G10; research/drone_transmitters/16_uat_adsl_fanet_ogn.md section D). Synthetic validation only.

PHY = FLARM (868.2 / 868.4 MHz, 2-FSK +-50 kHz, 100 kchip/s Manchester). After the preamble the decoded stream reads
0A F3 65 6C, then 20 data bytes (5 little-endian words) and 6 LDPC (208,160) parity bytes. The LDPC matrix exists only
in unlicensed / GPL code and is not used: a packet is verified instead by zero Manchester violations over all 208
coded bits (noise: 2^-208), the header's even-parity bit and a plausible position (latitude / longitude in range, time
< 60, fix quality <= 2). Label "-like": no FEC check.
Header word: address (24 bits), address type (random / ICAO / FLARM / OGN), non-position flag, parity, relay count,
encrypted, emergency. Position words 1-4 are whitened (TEA, zero key, 8 loops on word pairs) unless encrypted
(owner key: only the header is reported). Aircraft type 13 = UAV. Nothing is transmitted.
"""
from __future__ import annotations

import numpy as np

from . import flarm as FL
from . import sik as SIK

SYNC = 0x0AF3656C
DELTA = 0x9E3779B9
ADDR_TYPES = ("random", "ICAO", "FLARM", "OGN")
AIRCRAFT = ("unknown", "glider", "tow plane", "helicopter", "skydiver", "drop plane", "hang-glider", "paraglider",
            "powered aircraft", "jet", "UFO", "balloon", "airship", "UAV", "ground vehicle", "static object")
M = 0xFFFFFFFF


def tea_key0(v0, v1, loops=8, decrypt=False):
    if not decrypt:
        s = 0
        for _ in range(loops):
            s = (s + DELTA) & M
            v0 = (v0 + (((v1 << 4) & M) ^ ((v1 + s) & M) ^ (v1 >> 5))) & M
            v1 = (v1 + (((v0 << 4) & M) ^ ((v0 + s) & M) ^ (v0 >> 5))) & M
        return v0, v1
    s = (DELTA * loops) & M
    for _ in range(loops):
        v1 = (v1 - (((v0 << 4) & M) ^ ((v0 + s) & M) ^ (v0 >> 5))) & M
        v0 = (v0 - (((v1 << 4) & M) ^ ((v1 + s) & M) ^ (v1 >> 5))) & M
        s = (s - DELTA) & M
    return v0, v1


def _uvr(v, n):
    r, m = v >> n, v & ((1 << n) - 1)
    return (m, (1 << n) + 1 + 2 * m, 3 * (1 << n) + 2 + 4 * m, 7 * (1 << n) + 4 + 8 * m)[r]


def _s24(v):
    return v - (1 << 24) if v & 0x800000 else v


def parse(data20: bytes):
    """20 data bytes -> dict, or None when the header parity or the position plausibility fails."""
    w = list(np.frombuffer(bytes(data20), "<u4").astype(np.int64))
    h = int(w[0])
    if bin(h & 0x0FFFFFFF).count("1") % 2:
        return None
    out = dict(protocol="OGN tracker-like (OGNTP)", address=f"{h & 0xFFFFFF:06X}", address_type=ADDR_TYPES[(h >> 24) & 3],
               position_packet=not (h >> 26) & 1, relay_count=(h >> 28) & 3, encrypted=bool((h >> 30) & 1),
               emergency=bool(h >> 31), header_parity_ok=True, fec="LDPC not checked (matrix not available)")
    if out["encrypted"] or not out["position_packet"]:
        return out
    w[1], w[2] = tea_key0(int(w[1]), int(w[2]), decrypt=True)
    w[3], w[4] = tea_key0(int(w[3]), int(w[4]), decrypt=True)
    lat = ((_s24(w[1] & 0xFFFFFF) << 3) + 4) / 600000
    lon = ((_s24(w[2] & 0xFFFFFF) << 4) + 8) / 600000
    t, fixq = (w[1] >> 24) & 0x3F, (w[1] >> 30) & 3
    if abs(lat) > 90 or abs(lon) > 180 or t >= 60 or fixq > 2:
        return None
    at = (w[4] >> 20) & 0xF
    out.update(lat=round(lat, 5), lon=round(lon, 5), time_s=t, fix_quality=fixq, alt_m=_uvr(w[3] & 0x3FFF, 12),
               speed_ms=_uvr((w[3] >> 14) & 0x3FF, 8) / 10, heading_deg=round(((w[4] & 0x3FF) * 3600 + 512 >> 10) / 10, 1),
               aircraft_type=AIRCRAFT[at], uav=at == 13, stealth=bool((w[4] >> 19) & 1))
    return out


def packet(address: int, addr_type=3, lat=0.0, lon=0.0, alt_m=100, aircraft=13, t=10) -> bytes:
    """Synthesis / tests: 20 data bytes + 6 placeholder parity bytes (the LDPC is not computed: never checked)."""
    h = (address & 0xFFFFFF) | (addr_type & 3) << 24
    if bin(h & 0x0FFFFFFF).count("1") % 2:
        h |= 1 << 27
    w1 = (int(round(lat * 600000 / 8)) & 0xFFFFFF) | t << 24 | 1 << 30
    w2 = (int(round(lon * 600000 / 16)) & 0xFFFFFF) | 1 << 31
    w3 = min(alt_m, 4095) & 0x3FFF
    w4 = (aircraft & 0xF) << 20
    w1, w2 = tea_key0(w1, w2)
    w3, w4 = tea_key0(w3, w4)
    return np.array([h, w1, w2, w3, w4], "<u4").tobytes() + bytes(6)


def air_chips(pkt26: bytes) -> np.ndarray:
    by = SYNC.to_bytes(4, "big") + bytes(pkt26)
    return FL.manchester([(v >> (7 - i)) & 1 for v in by for i in range(8)])


def try_decode(y, ofs, chip_rate=FL.CHIP_RATE):
    """One burst at ~8 samples / chip centred on the carrier -> dict or None."""
    pat = np.array([(SYNC >> (31 - i)) & 1 for i in range(32)], np.int8) * 2 - 1
    for ph, chips in SIK.demod_bits(y, ofs, chip_rate):
        for inv in (False, True):
            ch = 1 - chips if inv else chips
            for odd in (0, 1):
                bits, _ = FL._decode(ch[odd:])
                if len(bits) < 32 + 208:
                    continue
                c = np.correlate(bits.astype(np.int8) * 2 - 1, pat, "valid")
                for p in np.flatnonzero(c >= 30):
                    pairs = ch[odd + 2 * (p + 32):odd + 2 * (p + 32 + 208)]
                    if len(pairs) < 416:
                        continue
                    pairs = pairs.reshape(-1, 2)
                    if np.any(pairs[:, 0] == pairs[:, 1]):
                        continue                                    # every coded bit must be valid Manchester
                    b = bits[p + 32:p + 32 + 160]
                    d = parse(bytes(int("".join(map(str, b[8 * i:8 * i + 8])), 2) for i in range(20)))
                    if d:
                        d.update(sync_bit_errors=int((32 - c[p]) // 2), polarity_inverted=inv)
                        return d
    return None
