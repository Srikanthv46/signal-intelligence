"""RDS codec (IEC 62106 / EN 50067; research/signal_gaps/06_vhf_below_300.md section 8): (26,16) block code with
g(x) = x^10 + x^8 + x^7 + x^5 + x^4 + x^3 + 1, offset words A / B / C / C' / D, block sync, group fields (PI, PTY,
TP, PS name, RadioText, clock time). Only error-free blocks are used (no burst correction). Pure Python."""
from __future__ import annotations

G = 0x5B9                          # x^10 + x^8 + x^7 + x^5 + x^4 + x^3 + 1
OFFSETS = {"A": 0x0FC, "B": 0x198, "C": 0x168, "Cp": 0x350, "D": 0x1B4}
BY_SYNDROME = {v: k for k, v in OFFSETS.items()}
NEXT = {"A": ("B",), "B": ("C", "Cp"), "C": ("D",), "Cp": ("D",), "D": ("A",)}
PTY = ["none", "News", "Current affairs", "Information", "Sport", "Education", "Drama", "Culture", "Science",
       "Varied", "Pop music", "Rock music", "Easy listening", "Light classical", "Serious classical", "Other music",
       "Weather", "Finance", "Children's programmes", "Social affairs", "Religion", "Phone-in", "Travel", "Leisure",
       "Jazz music", "Country music", "National music", "Oldies music", "Folk music", "Documentary", "Alarm test",
       "Alarm"]


def _mod(v: int, nbits: int) -> int:
    """v (nbits wide) mod g(x)."""
    for i in range(nbits - 1, 9, -1):
        if v >> i & 1:
            v ^= G << (i - 10)
    return v


def syndrome(word26: int) -> int:
    return _mod(word26, 26)


def encode_block(data16: int, offset: str) -> int:
    """26-bit block: data, then (data x^10 mod g) xor the offset word."""
    return (data16 << 10) | (_mod(data16 << 10, 26) ^ OFFSETS[offset])


def encode_group(blocks: tuple[int, int, int, int], version_b: bool = False) -> list[int]:
    """104 bits (MSB first) of a group."""
    offs = ("A", "B", "Cp" if version_b else "C", "D")
    out = []
    for d, o in zip(blocks, offs):
        w = encode_block(d, o)
        out += [(w >> (25 - i)) & 1 for i in range(26)]
    return out


def blocks(bits: list[int]) -> list[tuple[int, str, int]]:
    """Error-free blocks in sequence: (bit position, offset name, 16-bit data). Sync: a position whose 26-bit window
    has a valid offset and is followed 26 bits later by the next offset in the A, B, C/C', D cycle; from there blocks
    are read every 26 bits while they stay valid and in sequence."""
    n = len(bits)
    words = []
    for i in range(n - 25):
        w = 0
        for b in bits[i:i + 26]:
            w = (w << 1) | b
        words.append(w)
    synd = [syndrome(w) for w in words]
    out, used = [], set()
    i = 0
    while i + 26 < len(words):
        a = BY_SYNDROME.get(synd[i])
        b = BY_SYNDROME.get(synd[i + 26]) if i + 26 < len(words) else None
        if a and b and b in NEXT[a] and i not in used:
            j, cur = i, a
            while j < len(words):
                o = BY_SYNDROME.get(synd[j])
                if o is None or (j != i and o not in NEXT[cur]):
                    break
                out.append((j, o, words[j] >> 10))
                used.add(j)
                cur = o
                j += 26
            i = j
            continue
        i += 1
    return out


def groups(blks: list[tuple[int, str, int]]) -> list[dict]:
    """Groups from blocks in sequence: a group starts at an A block; missing blocks stay None."""
    out = []
    for k, (pos, o, d) in enumerate(blks):
        if o != "A":
            continue
        g = {"pos": pos, "A": d, "B": None, "C": None, "D": None, "version_b": False}
        for pos2, o2, d2 in blks[k + 1:k + 4]:
            step = (pos2 - pos) // 26
            if pos2 - pos != 26 * step or step > 3:
                break
            key = {"B": "B", "C": "C", "Cp": "C", "D": "D"}[o2]
            g[key] = d2
            if o2 == "Cp":
                g["version_b"] = True
        out.append(g)
    return out


def _chars(v: int) -> str:
    """Two characters; 0x0D (RadioText end) is kept as CR, other non-printables become '.'."""
    return "".join(chr(c) if 32 <= c < 127 else ("\r" if c == 13 else ".") for c in ((v >> 8) & 0xFF, v & 0xFF))


def decode(grps: list[dict]) -> dict:
    """PI, PTY, TP, PS name, RadioText and clock time from groups (partial names keep '_' for unseen characters)."""
    pis = [g["A"] for g in grps]
    out: dict = {}
    if not pis:
        return out
    pi = max(set(pis), key=lambda p: (pis.count(p), -pis.index(p)))
    out["pi"] = f"{pi:04X}"
    out["pi_country_nibble"] = pi >> 12
    ps = ["_"] * 8
    rt = ["_"] * 64
    types = []
    for g in grps:
        if g["A"] != pi or g["B"] is None:
            continue
        b = g["B"]
        gt, vb = b >> 12, bool(b >> 11 & 1)
        types.append(f"{gt}{'B' if vb else 'A'}")
        out["pty"] = (b >> 5) & 31
        out["pty_name"] = PTY[out["pty"]]
        out["tp"] = bool(b >> 10 & 1)
        if gt == 0 and g["D"] is not None:
            seg = b & 3
            c = _chars(g["D"])
            ps[2 * seg], ps[2 * seg + 1] = c[0], c[1]
            out["ta"] = bool(b >> 4 & 1)
            out["music"] = bool(b >> 3 & 1)
        elif gt == 2:
            seg = b & 15
            if not vb and g["C"] is not None and g["D"] is not None:
                for k, ch in enumerate(_chars(g["C"]) + _chars(g["D"])):
                    rt[4 * seg + k] = ch
            elif vb and g["D"] is not None:
                for k, ch in enumerate(_chars(g["D"])):
                    rt[2 * seg + k] = ch
        elif gt == 4 and not vb and g["C"] is not None and g["D"] is not None:
            mjd = ((b & 3) << 15) | (g["C"] >> 1)
            hour = ((g["C"] & 1) << 4) | (g["D"] >> 12)
            minute = (g["D"] >> 6) & 63
            off = (g["D"] & 31) * (-1 if g["D"] >> 5 & 1 else 1)
            y1 = int((mjd - 15078.2) / 365.25)
            m1 = int((mjd - 14956.1 - int(y1 * 365.25)) / 30.6001)
            day = mjd - 14956 - int(y1 * 365.25) - int(m1 * 30.6001)
            k = 1 if m1 in (14, 15) else 0
            out["clock_time"] = {"utc": f"{1900 + y1 + k:04d}-{m1 - 1 - 12 * k:02d}-{day:02d}T{hour:02d}:{minute:02d}Z",
                                 "local_offset_h": off / 2}
    if any(c != "_" for c in ps):
        out["ps"] = "".join(ps)
    txt = "".join(rt).split("\r")[0].rstrip("_")
    if txt.strip("_"):
        out["radiotext"] = txt
    out["group_types"] = sorted(set(types))
    return out
