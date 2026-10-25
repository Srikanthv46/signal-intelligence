"""Sub-audible squelch signalling on an FM discriminator output: CTCSS tones and DCS codes.

DCS (Digital-Coded Squelch, 'DPL'): a 23-bit Golay(23,12) codeword sent continuously
as NRZ FSK at 134.4 bit/s, below 300 Hz. Bits 0-8 are the 3-digit octal code (LSB
first), bits 9-11 the fixed marker 100 (bit 11 set), bits 12-22 the parity of the
generator x^11+x^10+x^6+x^5+x^4+x^2+1 (0xC75). Because the code is cyclic, a word
rotated in time is another valid codeword: every DCS code has aliases that radios treat
as the same squelch (023N = 340N = 766N = 047I). This convention reproduces those
published equivalences; the reciprocal polynomial does not (tests/receivers/test_det_voice_fm.py).

Polarity: a 1 bit is taken as a positive frequency shift ('N'); with the opposite
convention N and I swap. Confirmed only against what radios display (README).
"""
from __future__ import annotations

from collections import Counter

import numpy as np
from scipy import signal as ss

from .channels import CTCSS, DCS_CODES

DCS_RATE = 134.4
G = 0xC75
MASK23 = 0x7FFFFF
SUB_FS = 2000.0                    # sub-audio analysis rate


def golay_syndrome(w: int, g: int = G) -> int:
    """Remainder of w(x) / g(x); bit i of w is the coefficient of x^i. 0 = codeword."""
    for i in range(22, 10, -1):
        if w >> i & 1:
            w ^= g << (i - 11)
    return w


def dcs_encode(code_octal: str) -> int:
    """23-bit DCS word (bit 0 sent first) for a code such as "023"."""
    d = int(code_octal, 8) | 0x800
    for p in range(2048):
        if golay_syndrome(d | p << 12) == 0:
            return d | p << 12
    raise ValueError(code_octal)


def _rot(w: int, r: int) -> int:
    return ((w >> r) | (w << (23 - r))) & MASK23


def dcs_aliases(word: int) -> list[tuple[str, str]]:
    """Every (code, polarity) whose transmitted bit stream is the same as ``word``'s."""
    out = []
    for inv in (0, 1):
        v = word ^ (MASK23 if inv else 0)
        for r in range(23):
            x = _rot(v, r)
            if (x >> 9) & 7 == 4 and golay_syndrome(x) == 0:
                out.append((f"{x & 0x1FF:03o}", "I" if inv else "N"))
    return sorted(set(out))


def to_subaudio(disc: np.ndarray, fs: float) -> np.ndarray:
    """Discriminator output (Hz) -> 0-300 Hz band at SUB_FS, mean removed."""
    sos = ss.butter(6, 300, fs=fs, output="sos")
    a = ss.sosfiltfilt(sos, disc - np.mean(disc))
    t = np.arange(0, len(a) / fs, 1 / SUB_FS)
    return np.interp(t, np.arange(len(a)) / fs, a)


def dcs(a: np.ndarray, min_words: int = 2) -> dict | None:
    """Decode DCS from the 0-300 Hz band ``a`` (SUB_FS). Needs >= min_words + 1 codewords
    (~0.35 s). Returns the code, its aliases, and how many 23-bit windows verified."""
    if len(a) / SUB_FS * DCS_RATE < 23 * (min_words + 1):
        return None
    t = np.arange(len(a)) / SUB_FS
    best = None
    for ph in np.linspace(0, 1 / DCS_RATE, 16, endpoint=False):           # bit-clock phase: widest eye
        v = np.interp(np.arange(ph + 0.5 / DCS_RATE, t[-1], 1 / DCS_RATE), t, a)
        m = np.abs(v)
        eye = np.mean(m) / (np.std(m) + 1e-9)
        if best is None or eye > best[0]:
            best = (eye, v)
    v = best[1]
    bits = (v > 0).astype(np.int64)
    votes: Counter = Counter()
    words: dict = {}
    for s in range(len(bits) - 22):
        w = int(np.dot(bits[s:s + 23], 1 << np.arange(23, dtype=np.int64)))
        for inv in (0, 1):
            x = w ^ (MASK23 if inv else 0)
            if (x >> 9) & 7 == 4 and golay_syndrome(x) == 0:
                k = (f"{x & 0x1FF:03o}", "I" if inv else "N")
                votes[k] += 1
                words.setdefault(k, w)                                    # as received (polarity kept)
    if not votes:
        return None
    top, n = votes.most_common(1)[0]
    if n < min_words:
        return None
    family = dcs_aliases(words[top])
    n_family = sum(votes[k] for k in family)
    std = [k for k in family if k[0] in DCS_CODES]
    primary = min(std or family, key=lambda k: (k[1] != "N", k[0]))   # a standard code, normal polarity first
    rep = np.mean(bits[23:] == bits[:-23])
    return {"type": "DCS", "dcs_code": f"D{primary[0]}{primary[1]}",
            "dcs_aliases": [f"D{c}{p}" for c, p in family if (c, p) != primary],
            "dcs_standard": primary[0] in DCS_CODES, "dcs_words_verified": int(n_family),
            "dcs_bits": int(len(bits)), "dcs_bit_repeat_agreement": round(float(rep), 3),
            "dcs_deviation_hz": round(float(np.median(np.abs(v))), 0),
            "dcs_note": "Golay(23,12)-verified; N/I assumes a 1 bit = positive deviation"}


def ctcss(a: np.ndarray, min_s: float = 0.3) -> dict | None:
    """CTCSS tone from the 0-300 Hz band ``a`` (SUB_FS): spectral peak 60-260 Hz, >= 15 dB
    over the band's median, matched to the standard tone list."""
    if len(a) / SUB_FS < min_s:
        return None
    n = max(len(a), int(SUB_FS * 10))                                     # <= 0.1 Hz bins
    A = np.abs(np.fft.rfft(a * np.hanning(len(a)), n)) ** 2
    F = np.fft.rfftfreq(n, 1 / SUB_FS)
    band = (F >= 60) & (F <= 260)
    k = np.flatnonzero(band)[np.argmax(A[band])]
    snr = 10 * np.log10(A[k] / (np.median(A[band]) + 1e-30))
    if snr < 15:
        return None
    y0, y1, y2 = np.log(A[k - 1:k + 2] + 1e-30)                           # parabolic peak interpolation
    f = F[k] + 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2 + 1e-30) * (F[1] - F[0])
    std = min(CTCSS, key=lambda c: abs(c - f))
    res = SUB_FS / len(a)
    ok = abs(std - f) <= max(0.6, 0.3 * res)
    return {"type": "CTCSS", "ctcss_hz": std if ok else None, "ctcss_measured_hz": round(float(f), 2),
            "ctcss_snr_db": round(float(snr), 1), "ctcss_standard": bool(ok),
            "ctcss_deviation_hz": round(float(np.sqrt(2) * np.std(a)), 0)}


def subaudible(disc: np.ndarray, fs: float) -> dict | None:
    """DCS first (its line comb would otherwise read as a 'tone'), then CTCSS."""
    a = to_subaudio(disc, fs)
    return dcs(a) or ctcss(a)
