"""AIS codec (ITU-R M.1371-5; research/signal_gaps/06_vhf_below_300.md section 5): HDLC framing (flags 0x7E, bit
stuffing, octets sent LSB first), NRZI, CRC-16/X.25 (good-frame residue 0xF0B8), the 6-bit NMEA armouring and the
fields of message types 1/2/3, 4, 5, 18, 21 and 24. Pure Python."""
from __future__ import annotations

TRAINING = [0, 1] * 12
FLAG = [0, 1, 1, 1, 1, 1, 1, 0]
NAV_STATUS = ["under way using engine", "at anchor", "not under command", "restricted manoeuvrability",
              "constrained by draught", "moored", "aground", "engaged in fishing", "under way sailing",
              "reserved (HSC)", "reserved (WIG)", "power-driven towing astern", "power-driven pushing / towing",
              "reserved", "AIS-SART / MOB / EPIRB active", "not defined"]
MID_FLAG = {419: "India"}


def crc_x25(octets: bytes) -> int:
    """CRC-16/X.25: reflected 0x8408, init 0xFFFF, final xor 0xFFFF."""
    c = 0xFFFF
    for b in octets:
        c ^= b
        for _ in range(8):
            c = (c >> 1) ^ 0x8408 if c & 1 else c >> 1
    return c ^ 0xFFFF


def residue_ok(octets: bytes) -> bool:
    """Data + FCS (low byte first): the register before the final xor ends at 0xF0B8."""
    return len(octets) >= 3 and (crc_x25(octets) ^ 0xFFFF) == 0xF0B8


def armor_to_bits(payload: str, fill: int = 0) -> list[int]:
    bits = []
    for ch in payload:
        v = ord(ch) - 48
        if v > 40:
            v -= 8
        bits += [(v >> (5 - i)) & 1 for i in range(6)]
    return bits[:len(bits) - fill] if fill else bits


def bits_to_octets(bits: list[int]) -> bytes:
    """MSB-first packing (the message's own bit order); a partial last octet is zero-padded."""
    out = bytearray()
    for i in range(0, len(bits), 8):
        v = 0
        for b in bits[i:i + 8]:
            v = (v << 1) | b
        out.append(v << (8 - len(bits[i:i + 8])))
    return bytes(out)


def frame_bits(payload_bits: list[int]) -> list[int]:
    """On-air bits before NRZI: training, flag, stuffed (data + FCS, each octet LSB first), flag."""
    data = bits_to_octets(payload_bits)
    fcs = crc_x25(data)
    octs = data + bytes([fcs & 0xFF, fcs >> 8])
    raw = [(o >> i) & 1 for o in octs for i in range(8)]
    stuffed, run = [], 0
    for b in raw:
        stuffed.append(b)
        run = run + 1 if b else 0
        if run == 5:
            stuffed.append(0)
            run = 0
    return TRAINING + FLAG + stuffed + FLAG


def nrzi(bits: list[int], level: int = 1) -> list[int]:
    """0 = toggle, 1 = hold."""
    out = []
    for b in bits:
        if b == 0:
            level ^= 1
        out.append(level)
    return out


def unnrzi(levels: list[int]) -> list[int]:
    return [1 if levels[i] == levels[i - 1] else 0 for i in range(1, len(levels))]


def frames(bits: list[int]) -> list[tuple[int, bytes]]:
    """(bit index of the opening flag, data octets without the FCS) of every CRC-valid frame in an NRZI-decoded bit
    stream: an opening flag after >= 6 alternating training bits, de-stuffing up to the closing flag."""
    out = []
    n = len(bits)
    i = 6
    while i + 8 <= n:
        if bits[i:i + 8] != FLAG or bits[i - 6:i] != [0, 1, 0, 1, 0, 1] and bits[i - 6:i] != [1, 0, 1, 0, 1, 0]:
            i += 1
            continue
        raw, run, j, ok = [], 0, i + 8, False
        while j < n:
            if bits[j:j + 8] == FLAG:
                ok = True
                break
            b = bits[j]
            if run == 5:
                if b == 1:
                    break                                   # six ones: abort, not data
                run = 0
                j += 1
                continue
            raw.append(b)
            run = run + 1 if b else 0
            j += 1
            if len(raw) > 1100:
                break
        if ok and len(raw) >= 40 and len(raw) % 8 == 0:
            octs = bytes(sum(raw[k + q] << q for q in range(8)) for k in range(0, len(raw), 8))
            if residue_ok(octs):
                out.append((i, octs[:-2]))
                i = j + 8
                continue
        i += 1
    return out


def _u(bits, a, n):
    v = 0
    for b in bits[a:a + n]:
        v = (v << 1) | b
    return v


def _s(bits, a, n):
    v = _u(bits, a, n)
    return v - (1 << n) if v >> (n - 1) & 1 else v


def _text(bits, a, nchars):
    s = ""
    for k in range(nchars):
        if a + 6 * k + 6 > len(bits):
            break
        v = _u(bits, a + 6 * k, 6)
        s += chr(v + 64) if v < 32 else chr(v)
    return s.rstrip("@ ").replace("@", "")


def _pos(bits, a_lon, a_lat):
    lon, lat = _s(bits, a_lon, 28) / 600000.0, _s(bits, a_lat, 27) / 600000.0
    return (None if abs(lon) > 180 else round(lon, 6)), (None if abs(lat) > 90 else round(lat, 6))


def fields(octs: bytes, nbits: int | None = None) -> dict:
    """Message fields from the data octets (MSB-first message bits)."""
    bits = [(o >> (7 - i)) & 1 for o in octs for i in range(8)]
    if nbits:
        bits = bits[:nbits]
    t = _u(bits, 0, 6)
    mmsi = _u(bits, 8, 30)
    out = {"msg_type": t, "repeat": _u(bits, 6, 2), "mmsi": f"{mmsi:09d}"}
    ms = f"{mmsi:09d}"
    # MID: the first three digits of a ship; 00MID... base station, 99MID... aid to navigation, 111MID... SAR aircraft
    mid = int(ms[3:6]) if ms.startswith("111") else int(ms[2:5]) if ms[:2] in ("00", "99") else int(ms[:3])
    if 200 <= mid <= 799:
        out["mid"] = mid
        if mid in MID_FLAG:
            out["flag"] = MID_FLAG[mid]
    if t in (1, 2, 3) and len(bits) >= 168:
        lon, lat = _pos(bits, 61, 89)
        sog, cog, hdg = _u(bits, 50, 10), _u(bits, 116, 12), _u(bits, 128, 9)
        out.update(kind="class A position report", nav_status=NAV_STATUS[_u(bits, 38, 4)],
                   sog_kn=None if sog == 1023 else sog / 10, accuracy=bool(bits[60]), lon=lon, lat=lat,
                   cog_deg=None if cog == 3600 else cog / 10, heading_deg=None if hdg == 511 else hdg,
                   utc_second=_u(bits, 137, 6))
    elif t in (4, 11) and len(bits) >= 168:
        lon, lat = _pos(bits, 79, 107)
        out.update(kind="base station report" if t == 4 else "UTC / date response",
                   utc=f"{_u(bits, 38, 14):04d}-{_u(bits, 52, 4):02d}-{_u(bits, 56, 5):02d}T"
                       f"{_u(bits, 61, 5):02d}:{_u(bits, 66, 6):02d}:{_u(bits, 72, 6):02d}Z", lon=lon, lat=lat)
    elif t == 5 and len(bits) >= 420:
        out.update(kind="static and voyage data", imo=_u(bits, 40, 30), callsign=_text(bits, 70, 7),
                   name=_text(bits, 112, 20), ship_type=_u(bits, 232, 8),
                   length_m=_u(bits, 240, 9) + _u(bits, 249, 9), beam_m=_u(bits, 258, 6) + _u(bits, 264, 6),
                   draught_m=_u(bits, 294, 8) / 10, destination=_text(bits, 302, 20))
    elif t == 18 and len(bits) >= 168:
        lon, lat = _pos(bits, 57, 85)
        sog, cog, hdg = _u(bits, 46, 10), _u(bits, 112, 12), _u(bits, 124, 9)
        out.update(kind="class B position report", sog_kn=None if sog == 1023 else sog / 10, lon=lon, lat=lat,
                   cog_deg=None if cog == 3600 else cog / 10, heading_deg=None if hdg == 511 else hdg)
    elif t == 21 and len(bits) >= 272:
        lon, lat = _pos(bits, 164, 192)
        out.update(kind="aid to navigation", aid_type=_u(bits, 38, 5), name=_text(bits, 43, 20), lon=lon, lat=lat,
                   virtual=bool(bits[269]))
    elif t == 24 and len(bits) >= 160:
        part = _u(bits, 38, 2)
        out["kind"] = "class B static data"
        if part == 0:
            out["name"] = _text(bits, 40, 20)
        elif part == 1 and len(bits) >= 162:
            out.update(ship_type=_u(bits, 40, 8), callsign=_text(bits, 90, 7))
    return out
