"""
adsl.py -- ADS-L 4 SRD-860 M-band (EASA, Issue 2): passive decode of iConspicuity Traffic / Status packets (research,
PLAN G10; research/drone_transmitters/16_uat_adsl_fanet_ogn.md section B). Synthetic validation only.

PHY = the FLARM PHY (868.2 / 868.4 MHz, 2-GFSK +-50 kHz, 100 kchip/s Manchester -> 50 kb/s); after the preamble the
decoded stream reads F5 72 4B, then the packet length (0x18 = 24 for Traffic), the network header byte (version,
signature flag, key index: 0 = XXTEA with the zero key, 6 rounds; 3 = not scrambled; error-control mode), the ADS-L
Data (4k bytes, little-endian words) and the CRC-24 (Mode-S polynomial 0xFFF409, init 0, big-endian) over header +
data. A packet is reported only when the CRC passes. ADS-L Data: payload type (2 Traffic, 3 Status, 6 Remote ID
(reserved in Issue 2)), address-mapping table, 24-bit address, relay; Traffic: timestamp, flight state, aircraft
category (11-13 = UAS open / specific / certified), emergency, position (1/93206, 1/46603 deg), exponential speed /
altitude / climb, track, integrity fields. Nothing is transmitted.
"""
from __future__ import annotations

import numpy as np

from . import flarm as FL
from . import sik as SIK

SYNC24 = 0xF5724B
DELTA = 0x9E3779B9
TYPES = {2: "traffic", 3: "status", 4: "traffic uplink", 5: "FIS-B uplink", 6: "remote identification", 66: "OGN diagnostics"}
AMT = {0: "random / privacy", 5: "ICAO", 6: "FLARM", 7: "OGN tracker", 8: "FANET"}
CATEGORY = ("unknown", "light fixed-wing", "small-heavy fixed-wing", "light rotorcraft", "glider", "lighter-than-air",
            "ultralight / motorised hang-glider", "paraglider", "parachutist", "eVTOL / UAM", "gyrocopter",
            "UAS open category", "UAS specific category", "UAS certified category", "model plane", "heavy rotorcraft",
            "hang-glider", "paramotor")
FLIGHT = ("undefined", "ground", "airborne", "reserved")
EMERGENCY = ("undefined", "none", "general", "medical", "no communications", "unlawful interference", "downed",
             "reserved")


def crc24(data: bytes) -> int:
    """Mode-S CRC-24 (0xFFF409, init 0, MSB first)."""
    c = 0
    for v in data:
        c ^= v << 16
        for _ in range(8):
            c = ((c << 1) ^ 0xFFF409) & 0xFFFFFF if c & 0x800000 else (c << 1) & 0xFFFFFF
    return c


def _mx(z, y, s):
    return ((((z >> 5) ^ (y << 2)) + ((y >> 3) ^ (z << 4))) ^ ((s ^ y) + z)) & 0xFFFFFFFF   # key 0


def xxtea_encrypt(data: bytes, rounds=6) -> bytes:
    v = list(np.frombuffer(bytes(data), "<u4").astype(np.int64))
    n, s, z = len(v), 0, int(v[-1])
    for _ in range(rounds):
        s = (s + DELTA) & 0xFFFFFFFF
        for p in range(n):
            y = int(v[(p + 1) % n])
            v[p] = (int(v[p]) + _mx(z, y, s)) & 0xFFFFFFFF
            z = v[p]
    return np.array(v, "<u4").tobytes()


def xxtea_decrypt(data: bytes, rounds=6) -> bytes:
    v = list(np.frombuffer(bytes(data), "<u4").astype(np.int64))
    n, s, y = len(v), (rounds * DELTA) & 0xFFFFFFFF, int(v[0])
    for _ in range(rounds):
        for p in range(n - 1, -1, -1):
            z = int(v[(p - 1) % n])
            v[p] = (int(v[p]) - _mx(z, y, s)) & 0xFFFFFFFF
            y = v[p]
        s = (s - DELTA) & 0xFFFFFFFF
    return np.array(v, "<u4").tobytes()


def _exp_u(raw, n):
    e, base = raw >> n, raw & ((1 << n) - 1)
    return (1 << e) * ((1 << n) + base) - (1 << n)


def _exp_s(raw, n):
    sign = raw >> (n + 2)
    v = _exp_u(raw & ((1 << (n + 2)) - 1), n)
    return -v if sign else v


def _s24(v):
    return v - (1 << 24) if v & 0x800000 else v


def traffic(pl: bytes) -> dict:
    """Traffic payload (15 bytes, LSB-first bit fields) -> fields."""
    v = int.from_bytes(bytes(pl), "little")

    def f(off, n):
        return (v >> off) & ((1 << n) - 1)
    lat, lon = f(16, 24), f(40, 24)
    spd, alt, vr = f(64, 8), f(72, 14), f(86, 9)
    cat = f(8, 5)
    out = dict(timestamp_qs=f(0, 6), flight_state=FLIGHT[f(6, 2)], aircraft_category=CATEGORY[cat]
               if cat < len(CATEGORY) else f"reserved {cat}", uas=cat in (11, 12, 13), emergency=EMERGENCY[f(13, 3)],
               ground_track_deg=round(f(95, 9) * 360 / 512, 1), sil=f(104, 2), design_assurance=f(106, 2),
               nic=f(108, 4), hfom=f(112, 3), vfom=f(115, 2), velocity_accuracy=f(117, 2))
    if lat != 0x800000:
        out.update(lat=round(_s24(lat) / 93206, 5), lon=round(_s24(lon) / 46603, 5))
    out["ground_speed_ms"] = None if spd == 0xFF else _exp_u(spd, 6) / 4
    out["alt_m"] = None if alt == 0x3FFF else _exp_u(alt, 12) - 320
    out["climb_ms"] = None if vr == 0x1FF else _exp_s(vr, 6) / 8
    return out


def parse(pkt: bytes):
    """Length byte + network header + ADS-L Data + CRC-24 -> dict or None (CRC-verified)."""
    L = pkt[0]
    if L < 8 or len(pkt) < L + 1 or (L - 4) % 4:
        return None
    body = pkt[1:L + 1]
    if crc24(body[:-3]) != int.from_bytes(body[-3:], "big"):
        return None
    hdr, data = body[0], body[1:-3]
    key = (hdr >> 5) & 3
    if key == 0:
        data = xxtea_decrypt(data)
    elif key != 3:
        return None                                              # reserved key: not decodable
    typ = data[0]
    w = int.from_bytes(data[1:5], "little")
    amt = w & 0x3F
    out = dict(protocol="ADS-L (SRD-860 M-band)", crc_ok=True, length=L, version=hdr & 0x0F, signed=bool(hdr & 0x10),
               scrambled=key == 0, error_control="RS FEC" if hdr & 0x80 else "CRC-24",
               payload_type=TYPES.get(typ, f"0x{typ:02X}"),
               address_table=AMT.get(amt, f"manufacturer page {amt - 9}" if 9 <= amt <= 63 else f"reserved {amt}"),
               address=f"{(w >> 6) & 0xFFFFFF:06X}", relay=bool(w >> 31))
    if typ == 2 and len(data) >= 20:
        out.update(traffic(data[5:20]))
    return out


def packet(address: int, amt: int, payload15: bytes, typ=2, key=0, version=1) -> bytes:
    """Synthesis / tests: length + header + (scrambled) data + CRC-24."""
    data = bytes([typ]) + ((address & 0xFFFFFF) << 6 | amt & 0x3F).to_bytes(4, "little") + bytes(payload15)
    if key == 0:
        data = xxtea_encrypt(data)
    body = bytes([version & 0x0F | (key & 3) << 5]) + data
    c = crc24(body)
    body += c.to_bytes(3, "big")
    return bytes([len(body)]) + body


def traffic_payload(lat, lon, alt_m, category=11, speed_ms=10.0, track_deg=90.0, state=2) -> bytes:
    """Synthesis / tests: a Traffic payload (exponential fields by search over the code space)."""
    def enc_u(val, n):
        return min(range(4 << n), key=lambda r: abs(_exp_u(r, n) - val))
    v = 0
    v |= 2 << 0 | state << 6 | category << 8 | 1 << 13
    v |= (int(round(lat * 93206)) & 0xFFFFFF) << 16 | (int(round(lon * 46603)) & 0xFFFFFF) << 40
    v |= enc_u(speed_ms * 4, 6) << 64 | enc_u(alt_m + 320, 12) << 72 | 0 << 86
    v |= (int(round(track_deg * 512 / 360)) & 0x1FF) << 95
    return v.to_bytes(15, "little")


def air_bits(pkt: bytes) -> np.ndarray:
    """Synthesis / tests: chips after the FLARM-style preamble: Manchester(F5 72 4B + packet)."""
    by = bytes([0xF5, 0x72, 0x4B]) + bytes(pkt)
    return FL.manchester([(v >> (7 - i)) & 1 for v in by for i in range(8)])


def try_decode(y, ofs, chip_rate=FL.CHIP_RATE):
    """One burst at ~8 samples / chip centred on the carrier -> dict or None (CRC-24 verified)."""
    pat = np.array([(SYNC24 >> (23 - i)) & 1 for i in range(24)], np.int8) * 2 - 1
    for ph, chips in SIK.demod_bits(y, ofs, chip_rate):
        for inv in (False, True):
            ch = 1 - chips if inv else chips
            for odd in (0, 1):
                bits, _ = FL._decode(ch[odd:])
                if len(bits) < 24 + 8 * 12:
                    continue
                s = bits.astype(np.int8) * 2 - 1
                c = np.correlate(s, pat, "valid")
                for p in np.flatnonzero(c >= 22):
                    b = bits[p + 24:]
                    if len(b) < 8:
                        continue
                    L = int("".join(map(str, b[:8])), 2)
                    if len(b) < 8 * (L + 1):
                        continue
                    pkt = bytes(int("".join(map(str, b[8 * i:8 * i + 8])), 2) for i in range(L + 1))
                    d = parse(pkt)
                    if d:
                        d.update(sync_bit_errors=int((24 - c[p]) // 2), polarity_inverted=inv)
                        return d
    return None
