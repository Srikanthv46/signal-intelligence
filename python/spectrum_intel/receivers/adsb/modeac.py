"""SSR Mode A/C replies at 1090 MHz (research, signal-gaps S11; research/signal_gaps/05_aviation_radar_jammers.md).

A reply (ICAO Annex 10 Vol IV 3.1.1.6) is two framing pulses F1 / F2 20.3 us apart with up to 13 information pulse
positions every 1.45 us between them (C1 A1 C2 A2 C4 A4 X B1 D1 B2 D2 B4 D4; X is never transmitted) and an SPI
(ident) pulse 4.35 us after F2. Pulses are 0.45 us. The same format answers Mode A (identity, 4 octal digits ABCD)
and Mode C (pressure altitude in the Gillham code) interrogations, and the interrogation is usually not in view, so
every reply is reported both ways: the code, and the altitude when the bits form a valid Gillham altitude.
Military IFF Mode 3 equals Mode A; Mode 1 / 2 replies use the same pulse format.

No parity protects a reply: a reply is accepted only when its pulses sit on the 1.45 us grid (no pulse between
F1 and F2 off the grid, as overlapping replies or Mode S data would leave), with framing pulses of equal amplitude;
codes seen more than once in the block are the verified ones (an SSR beam gives an aircraft several replies).
Depends on numpy only.
"""
from __future__ import annotations

import numpy as np

FRAME_US = 20.3
STEP_US = 1.45
SPI_US = 24.65
PULSE_US = 0.45
# information pulse positions 1..13 after F1
POSITIONS = ("C1", "A1", "C2", "A2", "C4", "A4", "X", "B1", "D1", "B2", "D2", "B4", "D4")
TOL_US = 0.1                      # leading-edge tolerance on the grid (Annex 10: +-0.1 us), plus one sample
AMP_DB = 3.0                      # F1 / F2 amplitude agreement
WEAK_DB = 9.0                     # pulses this far below the framing pulses belong to noise or a weaker reply
PULSE_DB = 12.0                   # pulse threshold over the envelope median
SPECIAL = {"7500": "unlawful interference (hijack)", "7600": "radio failure", "7700": "emergency"}


def gray_to_bin(g: int) -> int:
    b = 0
    while g:
        b ^= g
        g >>= 1
    return b


def gillham_ft(p: dict) -> int | None:
    """Mode C altitude (ft) from the pulse bits, or None when they are not a valid Gillham code. 500 ft steps:
    Gray code D2 D4 A1 A2 A4 B1 B2 B4 (D1 is unused below 62,700 ft and must be 0); 100 ft steps: C1 C2 C4 in the
    'one-two-four' Gray code (1..5 valid, 7 stands for 5, reflected on odd 500 ft counts); altitude = 500 n500 + 100 n100 - 1300 ft."""
    if p["D1"] or p["X"]:
        return None
    g500 = 0
    for k in ("D2", "D4", "A1", "A2", "A4", "B1", "B2", "B4"):
        g500 = (g500 << 1) | p[k]
    g100 = (p["C1"] << 2) | (p["C2"] << 1) | p["C4"]
    n500, n100 = gray_to_bin(g500), gray_to_bin(g100)
    if n100 in (0, 5, 6):
        return None
    if n100 == 7:
        n100 = 5
    if n500 % 2:
        n100 = 6 - n100
    return n500 * 500 + n100 * 100 - 1300


def squawk(p: dict) -> str:
    """Mode A identity: digits A B C D, each (x4 x2 x1)."""
    return "".join(str(p[d + "4"] * 4 + p[d + "2"] * 2 + p[d + "1"]) for d in "ABCD")


def _pulses(env: np.ndarray, fs: float, thr: float) -> list[tuple[float, float, float]]:
    """(leading edge s, width s, amplitude) of each run above thr; the edge is where the run crosses half its
    peak, interpolated."""
    m = env > np.float32(thr)
    d = np.diff(np.r_[0, m.astype(np.int8), 0])
    st, en = np.flatnonzero(d == 1), np.flatnonzero(d == -1)
    out = []
    for a, b in zip(st, en):
        c = int(a + np.argmax(env[a:b]))
        pk = float(env[c])
        h = 0.5 * pk                  # width and edge at half the peak (filter ringing stays above thr)
        i = c
        while i > 0 and float(env[i - 1]) > h:
            i -= 1
        j = c + 1
        while j < len(env) and float(env[j]) > h:
            j += 1
        lo, hi = (float(env[i - 1]), float(env[i])) if i > 0 else (0.0, 0.0)
        e0 = i - 1 + (h - lo) / max(hi - lo, 1e-12) if i > 0 else float(i)
        out.append((e0 / fs, (j - i) / fs, pk))
    return out


def decode(env: np.ndarray, fs: float, busy: list[tuple[float, float]] = (), noise: float | None = None) -> list[dict]:
    """Mode A/C replies in a magnitude envelope at fs (>= 8 MSps). busy: (t0, t1) spans of decoded Mode S frames,
    whose pulses are not used."""
    noise = float(np.median(env) if noise is None else noise)
    thr = max(noise, 1e-12) * 10 ** (PULSE_DB / 20)
    P = sorted(p for p in _pulses(env, fs, thr)
               if 0.25e-6 <= p[1] <= 0.8e-6 and not any(a - 1e-6 <= p[0] <= b for a, b in busy))
    if len(P) < 2:
        return []
    t = np.array([p[0] for p in P])
    tol = TOL_US * 1e-6 + 1 / fs
    out, used = [], set()
    for i, (t1, _, a1) in enumerate(P):
        if i in used:
            continue
        j = int(np.searchsorted(t, t1 + FRAME_US * 1e-6 - tol))
        if j >= len(P) or abs(t[j] - t1 - FRAME_US * 1e-6) > tol:
            continue
        a2 = P[j][2]
        if abs(20 * np.log10(a1 / a2)) > AMP_DB:
            continue
        inner = range(i + 1, j)
        bits, ok = dict.fromkeys(POSITIONS, 0), True
        weak = min(a1, a2) * 10 ** (-WEAK_DB / 20)
        for k in inner:
            if P[k][2] < weak:
                continue
            q = (t[k] - t1) / (STEP_US * 1e-6)
            n = int(round(q))
            if not 1 <= n <= 13 or abs(t[k] - t1 - n * STEP_US * 1e-6) > tol or bits[POSITIONS[n - 1]]:
                ok = False                    # off-grid pulse: overlapping replies or not a reply at all
                break
            bits[POSITIONS[n - 1]] = 1
        if not ok or bits["X"]:
            continue
        ks = int(np.searchsorted(t, t1 + SPI_US * 1e-6 - tol))
        spi = ks < len(P) and abs(t[ks] - t1 - SPI_US * 1e-6) <= tol
        code = squawk(bits)
        alt = gillham_ft(bits)
        word = "".join(str(bits[k]) for k in POSITIONS)
        out.append({"t_s": round(t1, 7), "code": code, "altitude_ft_if_mode_c": alt, "spi": bool(spi),
                    "pulses": word, "snr_db": round(float(20 * np.log10(0.5 * (a1 + a2) / max(noise, 1e-12))), 1)})
        used.update(range(i, j + 1))
    return out


def summarise(replies: list[dict]) -> dict:
    """Per distinct reply word: count, code, altitude-if-Mode-C, special codes; repeated words are verified."""
    by = {}
    for r in replies:
        g = by.setdefault(r["pulses"], {"code": r["code"], "altitude_ft_if_mode_c": r["altitude_ft_if_mode_c"],
                                        "n": 0, "spi": False, "t_first_s": r["t_s"]})
        g["n"] += 1
        g["spi"] |= r["spi"]
    words = sorted(by.values(), key=lambda g: -g["n"])
    for g in words:
        g["repeated"] = g["n"] >= 2
        if g["code"] in SPECIAL:
            g["special_code"] = SPECIAL[g["code"]]
    return {"n_replies": len(replies), "n_distinct": len(words), "n_repeated": sum(g["repeated"] for g in words),
            "replies": words[:50]}


def synth_reply(code_bits: dict, fs: float, amp: float = 1.0, spi: bool = False) -> np.ndarray:
    """Envelope of one reply (tests): 0.45 us pulses with 0.1 us ramps."""
    n = int((SPI_US + 2) * 1e-6 * fs)
    x = np.zeros(n)
    tt = np.arange(n) / fs

    def pulse(t0):
        r = 0.1e-6
        x[:] = np.maximum(x, amp * np.clip(np.minimum((tt - t0) / r, (t0 + PULSE_US * 1e-6 - tt) / r + 1), 0, 1))

    pulse(0.0)
    pulse(FRAME_US * 1e-6)
    for k, name in enumerate(POSITIONS, 1):
        if code_bits.get(name):
            pulse(k * STEP_US * 1e-6)
    if spi:
        pulse(SPI_US * 1e-6)
    return x


def bits_for_code(code: str) -> dict:
    b = dict.fromkeys(POSITIONS, 0)
    for d, ch in zip("ABCD", code):
        v = int(ch)
        b[d + "4"], b[d + "2"], b[d + "1"] = (v >> 2) & 1, (v >> 1) & 1, v & 1
    return b
