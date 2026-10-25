"""OOK remotes on a pulse list (research, signal gaps S6; research/signal_gaps/04_short_range.md sections 1.1-1.3).
A pulse is (start us, mark us, space us); the last space of a transmission is inf.

    fixed_code(p)  EV1527 / PT2262 (Te = the PT2262 4 alpha): sync 1 Te high + 31 Te low, 24 PWM bits of 4 Te
                   (0 = 1 Te + 3 Te, 1 = 3 Te + 1 Te); no checksum, so repeats are the evidence; PT2262 tri-state when
                   every bit pair is 00 / 11 / 01
    keeloq(p)      Microchip HCS200 / 300 / 301: 12 warm-up pulses (1 Te + 1 Te), a 10 Te header gap, 66 PWM bits of
                   3 Te (short mark 1 Te = 1, long 2 Te = 0, rtl_433 convention), LSB first: bits 32-59 serial,
                   60-63 buttons, 64 battery low, 65 repeat. The 32-bit encrypted hopping part is not reported.
Pure Python; mirrored in C++ receivers/ism_devices/ook.
"""
from __future__ import annotations

TRI = {(0, 0): "0", (1, 1): "1", (0, 1): "F"}


def _median(v: list[float]) -> float:
    s = sorted(v)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def fixed_code(p: list[tuple[float, float, float]]) -> dict | None:
    words = []
    i = 0
    while i < len(p):
        t0, m, s = p[i]
        if not (100.0 <= m <= 1500.0 and 20 * m <= s <= 45 * m):
            i += 1
            continue
        bits, per = [], []
        for k in range(24):
            j = i + 1 + k
            if j >= len(p):
                break
            mk, sp = p[j][1], p[j][2]
            r = mk / m
            if 0.4 <= r <= 1.8:
                b = 0
            elif 2.2 <= r <= 4.2:
                b = 1
            else:
                break
            if k < 23:
                if not 2.6 * m <= mk + sp <= 5.4 * m:
                    break
                per.append(mk + sp)
            bits.append(b)
        if len(bits) == 24:
            code = sum(b << (23 - q) for q, b in enumerate(bits))
            if 0 < code < 0xFFFFFF:
                words.append((t0, code, _median(per) / 4, p[i + 24][0] + p[i + 24][1]))
                i += 25
                continue
        i += 1
    if not words:
        return None
    codes = [w[1] for w in words]
    best = max(set(codes), key=lambda c: (codes.count(c), -codes.index(c)))
    sel = [w for w in words if w[1] == best]
    bits = [(best >> (23 - q)) & 1 for q in range(24)]
    pairs = [(bits[2 * q], bits[2 * q + 1]) for q in range(12)]
    out = {"code": f"{best:06X}", "id": f"{best >> 4:05X}", "data": best & 0xF, "repeats": len(sel),
           "te_us": round(_median([w[2] for w in sel]), 1), "t_first_us": round(sel[0][0], 1),
           "t_last_us": round(sel[-1][0], 1), "t_end_us": round(sel[-1][3], 1), "words": len(words)}
    if all(pq in TRI for pq in pairs):
        out["tristate"] = "".join(TRI[pq] for pq in pairs)
    return out


def keeloq(p: list[tuple[float, float, float]]) -> dict | None:
    pkts = []
    i = 0
    while i + 12 + 66 <= len(p):
        te = p[i][1]
        ok = 200.0 <= te <= 800.0
        for k in range(12):
            if not ok:
                break
            mk, sp = p[i + k][1], p[i + k][2]
            ok = 0.6 * te <= mk <= 1.4 * te and (0.6 * te <= sp <= 1.4 * te if k < 11 else 7 * te <= sp <= 13 * te)
        if not ok:
            i += 1
            continue
        bits = []
        for k in range(66):
            mk, sp = p[i + 12 + k][1], p[i + 12 + k][2]
            r = mk / te
            if 0.5 <= r <= 1.5:
                bits.append(1)
            elif 1.5 < r <= 2.6:
                bits.append(0)
            else:
                break
            if k < 65 and not 2.2 * te <= mk + sp <= 3.8 * te:
                break
        if len(bits) == 66 and not all(bits[:64]):
            serial = sum(bits[32 + q] << q for q in range(28))
            btn = (bits[60] << 3) | (bits[61] << 2) | (bits[62] << 1) | bits[63]
            pkts.append({"t_us": round(p[i][0], 1), "serial": f"{serial:07X}",
                         "button": (btn & 0x08) | ((btn & 0x01) << 2) | (btn & 0x02) | ((btn & 0x04) >> 2),
                         "learn": btn == 0x0F, "battery_low": bool(bits[64]), "repeat": bool(bits[65]),
                         "te_us": round(te, 1)})
            i += 78
            continue
        i += 1
    if not pkts:
        return None
    first = pkts[0]
    same = [q for q in pkts if q["serial"] == first["serial"] and q["button"] == first["button"]]
    return {**{k: first[k] for k in ("serial", "button", "learn", "battery_low", "repeat", "te_us")},
            "packets": len(pkts), "matching_packets": len(same), "t_first_us": first["t_us"]}
