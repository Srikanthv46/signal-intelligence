"""Bluetooth BR/EDR (Classic) access codes (research, signal gaps S7; research/signal_gaps/04_short_range.md section 2).

The 64-bit sync word is a (64,30) expurgated BCH codeword over the 24-bit LAP (lower address part of BD_ADDR) plus a
6-bit Barker tail, whitened by PN 0x83848D96BBCC54FC; generator 0o260534236651 (Core Spec Vol 2 Part B 6.3.3; checked
against the published codewords of LAP 0x000000 and the GIAC 0x9E8B33). A window of the sliced 1 Mb/s GFSK bits is
an access code when it is within MAX_ERRORS bits of syncword(its own LAP bits): that validates the code (minimum
distance 14) and gives the LAP. The UAP (from the header HEC) needs the master clock and several packets; the NAP is
never sent, so no vendor is claimed.

    syncword(lap) -> int             bit k = the k-th bit sent
    try_decode(y, ofs) -> dict | None
"""
from __future__ import annotations

import numpy as np

PN = 0x83848D96BBCC54FC
G = 0o260534236651                    # 0x585713DA9, degree 34
GIAC, LIAC = 0x9E8B33, 0x9E8B00
MAX_ERRORS = 3                        # chance match of a random window: ~1e-8 (40 checked bits), ~2.5e-5 per burst
FILTER_HZ = 0.6e6
SEARCH_BITS = 160                     # the access code opens the packet: search its first bits only
RATE = 1e6


def _polymod(a: int, g: int) -> int:
    dg = g.bit_length() - 1
    while a.bit_length() - 1 >= dg:
        a ^= g << (a.bit_length() - 1 - dg)
    return a


def syncword(lap: int) -> int:
    barker = 0b101100 if ((lap >> 23) & 1) == 0 else 0b010011
    xt = (lap | (barker << 24)) ^ (PN >> 34)
    return ((xt << 34) | _polymod(xt << 34, G)) ^ PN


TAILS = (0b101100, 0b010011)          # the top 6 bits of every sync word: the Barker tail (the PN cancels there)


def find(bits) -> tuple[int, int, int] | None:
    """(bit position, LAP, bit errors) of the first access code in bits, or None. The Barker tail (top 6 sync bits)
    pre-filters the positions; the few left are checked against the codeword of their own LAP."""
    b = np.asarray(bits, dtype=np.int64)
    n = len(b)
    if n < 64:
        return None
    last = min(n - 64, SEARCH_BITS)
    top = np.zeros(last + 1, np.int64)
    for q in range(6):
        top |= b[58 + q:58 + q + last + 1] << q
    for p in np.flatnonzero((top == TAILS[0]) | (top == TAILS[1])):
        p = int(p)
        w = 0
        for q in range(64):
            w |= int(b[p + q]) << q
        lap = (w >> 34) & 0xFFFFFF
        err = bin(w ^ syncword(lap)).count("1")
        if err <= MAX_ERRORS:
            return p, lap, err
    return None


def try_decode(y: np.ndarray, ofs: float) -> dict | None:
    """Filter to +-0.6 MHz (FFT of the segment), discriminator integrated over one symbol, threshold mid-way between
    its 10th / 90th percentiles, the integer sample phases, either polarity; the first access code found."""
    sps = ofs / RATE
    if len(y) < 72 * sps:
        return None
    z = np.asarray(y).astype(np.complex128)
    f = np.fft.fftfreq(len(z), 1 / ofs)
    z = np.fft.ifft(np.fft.fft(z) * (np.abs(f) <= FILTER_HZ))
    d = np.angle(z[1:] * np.conj(z[:-1]))
    k = max(1, int(round(sps)))
    c = np.concatenate([[0.0], np.cumsum(d)])
    v = (c[k:] - c[:-k]) / k
    srt = np.sort(v)
    thr = 0.5 * (float(srt[int(0.1 * len(srt))]) + float(srt[int(0.9 * len(srt))]))
    n = int((len(v) - sps) / sps)
    for ph in range(int(sps)):
        idx = (ph + np.arange(n) * sps).astype(int)
        b = (v[idx] > thr).astype(np.int64)
        for pol in (0, 1):
            bits = b if pol == 0 else 1 - b
            r = find(bits)
            if r is not None:
                p, lap, err = r
                return {"protocol": "Bluetooth BR/EDR", "lap": f"{lap:06X}", "access_code_errors": err,
                        "inquiry": lap in (GIAC, LIAC), "giac": lap == GIAC, "access_code_ok": True,
                        "inverted": bool(pol), "t_offset_us": round((ph + p * sps) / ofs * 1e6, 2)}
    return None
