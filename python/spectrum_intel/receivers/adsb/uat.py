"""
uat.py -- UAT 978 MHz ADS-B downlink (DO-282B): demodulation, Reed-Solomon correction and message decode (research,
PLAN G14; research/drone_transmitters/16_uat_adsl_fanet_ogn.md section A). Synthetic validation only (no UAT traffic
in India).

PHY: binary CPFSK, 1.041667 Mb/s, +-312.5 kHz (h = 0.6), 1 = shift up. A message is the 36-bit sync 0xEACDDA4E2, then
either a Basic frame (18 data + 12 RS parity bytes) or a Long frame (34 + 14); the receiver tries Long first and keeps
it when the payload type is not 0, else Basic with payload type 0 (as FlightAware dump978, BSD-2). RS over GF(256),
primitive polynomial 0x187, first consecutive root 120 (core.rs, clean-room). Demodulation: the channel at 4 samples /
bit, the phase change over one bit, the slicer threshold centred on the sync's 1 / 0 levels (absorbs the CFO), sync
within 4 bit errors at the best of the 4 timing phases.
Decoded: payload type, address qualifier, 24-bit address, latitude / longitude (23 / 24 bits x 360 / 2^24), altitude
((raw - 41) x 25 ft), NIC, air / ground state, velocities, and from the Mode Status element (payload types 1 and 3):
emitter category (14 = unmanned aerial vehicle) and callsign (base 40). Ground uplinks (FIS-B) are not decoded.
"""
from __future__ import annotations

import numpy as np

from spectrum_intel.core.rs import RS

F_UAT = 978e6
RATE = 1.041667e6
DEV = 312.5e3
SPB = 4
FS = SPB * RATE
SYNC = 0xEACDDA4E2
RS_SHORT, RS_LONG = RS(0x187, 120, 1, 12), RS(0x187, 120, 1, 14)
CATEGORY = {0: "no information", 1: "light", 2: "small", 3: "large", 4: "high-vortex large", 5: "heavy",
            6: "highly manoeuvrable", 7: "rotorcraft", 9: "glider / sailplane", 10: "lighter-than-air",
            11: "parachutist / skydiver", 12: "ultralight / hang-glider / paraglider", 14: "unmanned aerial vehicle",
            15: "space / trans-atmospheric", 17: "surface emergency vehicle", 18: "surface service vehicle",
            19: "point obstacle", 20: "cluster obstacle", 21: "line obstacle"}
QUALIFIER = ("ADS-B ICAO", "ADS-B self-assigned", "TIS-B ICAO", "TIS-B track file", "surface vehicle", "fixed beacon",
             "ADS-R self-assigned", "reserved")
B40 = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ *??"


def _bits(m: bytes, b1, i1, b2, i2):
    """DO-282 bit range (byte, bit), 1-based, bit 1 = MSB, inclusive."""
    v = int.from_bytes(m, "big")
    n = 8 * len(m)
    lo, hi = (b1 - 1) * 8 + (i1 - 1), (b2 - 1) * 8 + (i2 - 1)
    return (v >> (n - 1 - hi)) & ((1 << (hi - lo + 1)) - 1)


def parse(m: bytes) -> dict:
    """Corrected payload (18 or 34 bytes) -> fields."""
    pt, aq = m[0] >> 3, m[0] & 7
    out = dict(payload_type=pt, address_qualifier=QUALIFIER[aq], address=m[1:4].hex().upper(),
               frame="long" if len(m) == 34 else "basic")
    lat_r, lon_r = _bits(m, 5, 1, 7, 7), _bits(m, 7, 8, 10, 7)
    if lat_r or lon_r:
        lat, lon = lat_r * 360 / (1 << 24), lon_r * 360 / (1 << 24)
        out.update(lat=round(lat - 180 if lat > 90 else lat, 5), lon=round(lon - 360 if lon > 180 else lon, 5))
    alt = _bits(m, 11, 1, 12, 4)
    out["alt_ft"] = None if alt == 0 else (alt - 41) * 25
    out["alt_type"] = "geometric" if _bits(m, 10, 8, 10, 8) else "pressure"
    out["nic"] = _bits(m, 12, 5, 12, 8)
    ag = _bits(m, 13, 1, 13, 2)
    out["air_ground"] = ("airborne subsonic", "airborne supersonic", "on ground", "reserved")[ag]
    if ag in (0, 1):
        k = 4 if ag == 1 else 1
        ns, nv = _bits(m, 13, 4, 13, 4), _bits(m, 13, 5, 14, 6)
        es, ev = _bits(m, 14, 7, 14, 7), _bits(m, 14, 8, 16, 1)
        vs, vv = _bits(m, 16, 3, 16, 3), _bits(m, 16, 4, 17, 4)
        if nv:
            out["north_kt"] = (nv - 1) * k * (-1 if ns else 1)
        if ev:
            out["east_kt"] = (ev - 1) * k * (-1 if es else 1)
        if vv:
            out["vert_fpm"] = (vv - 1) * 64 * (-1 if vs else 1)
    if len(m) == 34 and pt in (1, 3):                             # Mode Status
        w = [int.from_bytes(m[17 + 2 * i:19 + 2 * i], "big") for i in range(3)]
        cat = (w[0] // 1600) % 40
        cs = "".join(B40[(x // 1600) % 40] + B40[(x // 40) % 40] + B40[x % 40] for x in w)[1:]
        out.update(emitter_category=CATEGORY.get(cat, f"reserved {cat}"), emitter_category_code=cat, uav=cat == 14,
                   callsign=cs.rstrip(" *"), emergency=_bits(m, 24, 1, 24, 3), mops=_bits(m, 24, 4, 24, 6))
    return out


def frame(m: bytes) -> bytes:
    """Synthesis / tests: payload (18 / 34 bytes) -> data + RS parity."""
    return bytes(m) + (RS_LONG if len(m) == 34 else RS_SHORT).encode(m)


def payload(address: int, lat: float, lon: float, alt_ft: int, category=14, callsign="UAV1", qualifier=1) -> bytes:
    """Synthesis / tests: a Long payload of type 1 (SV + MS + AUXSV)."""
    m = bytearray(34)
    m[0] = 1 << 3 | qualifier
    m[1:4] = address.to_bytes(3, "big")
    v = int.from_bytes(m, "big")

    def put(b1, i1, b2, i2, val):
        nonlocal v
        lo, hi = (b1 - 1) * 8 + (i1 - 1), (b2 - 1) * 8 + (i2 - 1)
        sh = 8 * 34 - 1 - hi
        mask = ((1 << (hi - lo + 1)) - 1) << sh
        v = (v & ~mask) | ((val << sh) & mask)
    put(5, 1, 7, 7, int(round((lat % 180) * (1 << 24) / 360)) & 0x7FFFFF)
    put(7, 8, 10, 7, int(round((lon % 360) * (1 << 24) / 360)) & 0xFFFFFF)
    put(11, 1, 12, 4, alt_ft // 25 + 41)
    cs = (" " + callsign.upper()).ljust(9)[:9]
    idx = [B40.index(c) if c in B40 else 36 for c in cs]
    idx[0] = category
    for i in range(3):
        put(18 + 2 * i, 1, 19 + 2 * i, 8, idx[3 * i] * 1600 + idx[3 * i + 1] * 40 + idx[3 * i + 2])
    return v.to_bytes(34, "big")


def _decide(dphi, thr):
    return (dphi > thr).astype(np.uint8)


def demod(y: np.ndarray, max_err=4):
    """Channel at FS (4 samples / bit) -> [dict(t_s, frame bytes info)] of RS-corrected downlink messages."""
    ph = np.unwrap(np.angle(np.asarray(y, np.complex128)))
    pat = np.array([(SYNC >> (35 - i)) & 1 for i in range(36)], np.int8) * 2 - 1
    out, seen = [], []
    nlong = 8 * 48
    for p in range(SPB):
        idx = np.arange(p, len(ph) - SPB, SPB)
        d = ph[idx + SPB] - ph[idx]                                   # phase change over one bit
        s = np.where(d > 0, 1, -1).astype(np.int8)
        if len(s) < 36 + nlong:
            continue
        c = np.correlate(s, pat, "valid")
        for k in np.flatnonzero(c >= 36 - 2 * max_err):
            if any(abs(k * SPB + p - q) < 30 * SPB for q in seen):
                continue
            sd = d[k:k + 36]
            thr = (sd[pat > 0].mean() + sd[pat < 0].mean()) / 2      # centre on the sync levels: absorbs the CFO
            if np.count_nonzero(_decide(sd, thr) != (pat > 0)) > max_err:
                continue
            b = _decide(d[k + 36:k + 36 + nlong], thr)
            if len(b) < 8 * 30:
                continue
            raw = np.packbits(np.r_[b, np.zeros((-len(b)) % 8, np.uint8)]).tobytes()
            msg, nfix = None, -1
            if len(raw) >= 48:
                cw, n = RS_LONG.decode(raw[:48])
                if cw is not None and cw[0] >> 3 != 0:
                    msg, nfix = cw[:34], n
            if msg is None:
                cw, n = RS_SHORT.decode(raw[:30])
                if cw is not None and cw[0] >> 3 == 0:
                    msg, nfix = cw[:18], n
            if msg is None:
                continue
            seen.append(k * SPB + p)
            d_ = parse(msg)
            d_.update(t_s=(k * SPB + p) / FS, rs_corrected_bytes=nfix,
                      sync_bit_errors=int(np.count_nonzero(_decide(sd, thr) != (pat > 0))), hex=msg.hex())
            out.append(d_)
    return sorted(out, key=lambda m: m["t_s"])


def modulate(m: bytes, fs=FS, lead_bits=40) -> np.ndarray:
    """Synthesis / tests: CPFSK h = 0.6 of sync + frame at fs (rectangular frequency pulses)."""
    by = frame(m)
    bits = [(SYNC >> (35 - i)) & 1 for i in range(36)] + [(v >> (7 - i)) & 1 for v in by for i in range(8)]
    sps = fs / RATE
    n = int(np.ceil((len(bits) + 2 * lead_bits) * sps))
    t = np.arange(n) / sps - lead_bits
    k = np.clip(t.astype(int), -1, len(bits))
    f = np.where((k >= 0) & (k < len(bits)), np.array(bits + [0])[np.clip(k, 0, len(bits))] * 2.0 - 1, 0.0) * DEV
    return np.exp(2j * np.pi * np.cumsum(f) / fs).astype(np.complex64)
