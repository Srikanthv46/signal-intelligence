"""HDLC framing and AX.25 / APRS parsing (research, signal gaps S19; research/signal_gaps/06_vhf_below_300.md
section 7): NRZI, flags 0x7E, bit stuffing, octets LSB first, CRC-16/X.25 (good-frame residue 0xF0B8); AX.25 UI
addresses; APRS position (uncompressed, compressed, Mic-E), status, message and comment. Pure Python; mirrored in
C++ core/hdlc."""
from __future__ import annotations

FLAG = [0, 1, 1, 1, 1, 1, 1, 0]


def crc_x25(octets: bytes) -> int:
    """CRC-16/X.25: reflected 0x8408, init 0xFFFF, final xor 0xFFFF."""
    c = 0xFFFF
    for b in octets:
        c ^= b
        for _ in range(8):
            c = (c >> 1) ^ 0x8408 if c & 1 else c >> 1
    return c ^ 0xFFFF


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


def frame_bits(octets: bytes, n_flags: int = 1) -> list[int]:
    """Flags, stuffed (octets + FCS, LSB first), closing flag: the bits before NRZI."""
    fcs = crc_x25(octets)
    raw = [(o >> i) & 1 for o in octets + bytes([fcs & 0xFF, fcs >> 8]) for i in range(8)]
    stuffed, run = [], 0
    for b in raw:
        stuffed.append(b)
        run = run + 1 if b else 0
        if run == 5:
            stuffed.append(0)
            run = 0
    return FLAG * n_flags + stuffed + FLAG


def deframe(bits: list[int], min_octets: int = 18, max_bits: int = 3000) -> list[tuple[int, bytes]]:
    """(index of the opening flag, octets without the FCS) of every CRC-valid frame between flags."""
    out = []
    n = len(bits)
    i = 0
    while i + 8 <= n:
        if bits[i:i + 8] != FLAG:
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
                    break                                   # six ones: abort
                run = 0
                j += 1
                continue
            raw.append(b)
            run = run + 1 if b else 0
            j += 1
            if len(raw) > max_bits:
                break
        if ok and len(raw) >= 8 * min_octets and len(raw) % 8 == 0:
            octs = bytes(sum(raw[k + q] << q for q in range(8)) for k in range(0, len(raw), 8))
            if (crc_x25(octs) ^ 0xFFFF) == 0xF0B8:
                out.append((i, octs[:-2]))
                i = j
                continue
        i = j if ok and not raw else i + 1
    return out


def _call(b: bytes) -> str:
    c = "".join(chr(x >> 1) for x in b[:6]).rstrip()
    ssid = (b[6] >> 1) & 15
    return f"{c}-{ssid}" if ssid else c


def ax25(octs: bytes) -> dict | None:
    """AX.25 UI frame -> {dest, src, path, info}; None when the address field or control / PID is not a UI frame."""
    addrs, k = [], 0
    while k + 7 <= len(octs):
        a = octs[k:k + 7]
        if any(not (0x40 <= x >> 1 < 0x5B or x >> 1 == 0x20 or 0x30 <= x >> 1 <= 0x39) for x in a[:6]):
            return None
        addrs.append(a)
        k += 7
        if a[6] & 1:
            break
    if len(addrs) < 2 or not addrs[-1][6] & 1 or k + 2 > len(octs) or octs[k] != 0x03 or octs[k + 1] != 0xF0:
        return None
    path = [_call(a) + ("*" if a[6] & 0x80 else "") for a in addrs[2:]]
    return {"dest": _call(addrs[0]), "src": _call(addrs[1]), "path": path,
            "info": octs[k + 2:].decode("latin-1")}


def _ddmm(s: str, hemi: str, deg_digits: int) -> float | None:
    s = s.replace(" ", "0")
    try:
        v = int(s[:deg_digits]) + float(s[deg_digits:]) / 60
    except ValueError:
        return None
    return round(-v if hemi in "SW" else v, 5)


def _uncompressed(p: str) -> dict:
    if len(p) < 19 or p[7] not in "NS" or p[17] not in "EW":
        return {}
    lat, lon = _ddmm(p[0:7], p[7], 2), _ddmm(p[9:17], p[17], 3)
    if lat is None or lon is None or abs(lat) > 90 or abs(lon) > 180:
        return {}
    out = {"lat": lat, "lon": lon, "symbol": p[8] + p[18], "format": "uncompressed"}
    rest = p[19:]
    if len(rest) >= 7 and rest[3] == "/" and rest[:3].isdigit() and rest[4:7].isdigit():
        out.update(course_deg=int(rest[:3]), speed_kn=int(rest[4:7]))
        rest = rest[7:]
    out["comment"] = rest
    return out


def _compressed(p: str) -> dict:
    if len(p) < 13:
        return {}
    v = [ord(c) - 33 for c in p[1:9]]
    if any(not 0 <= x <= 90 for x in v):
        return {}
    lat = 90 - (v[0] * 91 ** 3 + v[1] * 91 ** 2 + v[2] * 91 + v[3]) / 380926
    lon = -180 + (v[4] * 91 ** 3 + v[5] * 91 ** 2 + v[6] * 91 + v[7]) / 190463
    out = {"lat": round(lat, 5), "lon": round(lon, 5), "symbol": p[0] + p[9], "format": "compressed"}
    c, s = ord(p[10]) - 33, ord(p[11]) - 33
    if p[10] != " " and 0 <= c <= 89:
        out.update(course_deg=c * 4, speed_kn=round(1.08 ** s - 1, 1))
    out["comment"] = p[13:]
    return out


def _mic_e(dest: str, info: str) -> dict:
    d = dest.split("-")[0]
    if len(d) != 6 or len(info) < 9:
        return {}
    digits = []
    for ch in d:
        if "0" <= ch <= "9":
            digits.append(ord(ch) - 48)
        elif "A" <= ch <= "J":
            digits.append(ord(ch) - 65)
        elif "P" <= ch <= "Y":
            digits.append(ord(ch) - 80)
        elif ch in "KLZ":
            digits.append(0)
        else:
            return {}
    lat = digits[0] * 10 + digits[1] + (digits[2] * 10 + digits[3] + (digits[4] * 10 + digits[5]) / 100) / 60
    north, off100, west = d[3] >= "P", d[4] >= "P", d[5] >= "P"
    deg = ord(info[1]) - 28 + (100 if off100 else 0)
    if 180 <= deg <= 189:
        deg -= 80
    elif 190 <= deg <= 199:
        deg -= 190
    mins = ord(info[2]) - 28
    if mins >= 60:
        mins -= 60
    hund = ord(info[3]) - 28
    lon = deg + (mins + hund / 100) / 60
    sp, dc, se = ord(info[4]) - 28, ord(info[5]) - 28, ord(info[6]) - 28
    speed = sp * 10 + dc // 10
    if speed >= 800:
        speed -= 800
    course = (dc % 10) * 100 + se
    if course >= 400:
        course -= 400
    if not (0 <= lat <= 90 and 0 <= lon <= 180):
        return {}
    return {"lat": round(lat if north else -lat, 5), "lon": round(-lon if west else lon, 5),
            "symbol": info[8] + info[7], "format": "Mic-E", "course_deg": course, "speed_kn": speed,
            "comment": info[9:]}


def _position(body: str) -> dict:
    """An uncompressed position starts with a latitude digit (or an ambiguity space); anything else is the symbol
    table of a compressed one (APRS 1.0.1 ch. 9)."""
    if not body:
        return {}
    return _uncompressed(body) if body[0].isdigit() or body[0] == " " else _compressed(body)


def aprs(frame: dict) -> dict:
    """APRS fields of an AX.25 UI frame's information field (data type identifier first)."""
    info = frame["info"]
    if not info:
        return {}
    t = info[0]
    out: dict = {"data_type": t}
    if t in "!=":
        out.update(_position(info[1:]))
    elif t in "/@" and len(info) >= 8:
        out["timestamp"] = info[1:8]
        out.update(_position(info[8:]))
    elif t in "`'":
        out.update(_mic_e(frame["dest"], info))
    elif t == ">":
        out["status"] = info[1:]
    elif t == ":" and len(info) >= 11 and info[10] == ":":
        out.update(addressee=info[1:10].strip(), message=info[11:])
    elif t == ";" and len(info) >= 18:
        out["object"] = info[1:10].strip()
        out.update(_uncompressed(info[18:]))
    return out
