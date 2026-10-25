"""
fsk_frame.py -- blind frame structure of any (G)FSK packet (research, "any transmitter": maximum information from an
unknown protocol). Hard decisions at 8 samples / bit (sik.demod_bits) at every timing phase; the phase with the longest
alternating run near the start wins. Reported: preamble length (bits), the 32 bits after the preamble (sync-word
candidate; most packet radios put a 16-32 bit sync word there), the bits that follow (packet length) and the polarity
convention (bit 1 = higher frequency). A sync word repeated across packets identifies the link / protocol even when the
protocol is unknown. A known sync word of a published packet radio (KNOWN_SYNCS, matched exactly within 4 bits of the
preamble end, either polarity, <= 1 bit error for 24 / 32-bit words) is named as a candidate ("-like": a sync word alone does not verify the protocol).
Nothing is decoded or verified.
"""
from __future__ import annotations

import numpy as np

from . import sik as SIK

MIN_PRE = 12                     # alternating bits for a preamble (random data: a 12-run ~ 1 in 2000 per position)
SYNC_BITS = 32
# (name, sync word hex, bits) of published packet-radio sync words; longest first (research files 05, 08, 10)
KNOWN_SYNCS = (("FrSky D16-like (CC2500)", 0xD391D391, 32), ("FlySky AFHDS 2A-like (A7105 ID code)", 0x5475C52A, 32),
               ("TBS Crossfire-like", 0x3E6B3E, 24), ("FrSky D8-like (CC2500)", 0xD391, 16),
               ("SiK-like (Si1000 / Si4432)", 0x2DD4, 16))


def _alt_runs(bits):
    """(start, length) of maximal alternating runs (consecutive bits differ)."""
    d = np.diff(bits.astype(np.int8)) != 0
    runs, i, n = [], 0, len(d)
    while i < n:
        if d[i]:
            j = i
            while j < n and d[j]:
                j += 1
            runs.append((i, j - i + 1))
            i = j
        else:
            i += 1
    return runs


def known_sync(bits, e):
    """The first KNOWN_SYNCS entry found within 4 bits of position e (either polarity; <= 1 bit error for 24 / 32-bit
    words, exact for 16-bit ones), or None."""
    for name, word, n in KNOWN_SYNCS:
        pat = np.array([(word >> (n - 1 - i)) & 1 for i in range(n)], np.uint8)
        max_err = 1 if n >= 24 else 0
        for d in (0, -1, 1, -2, 2, -3, 3, -4, 4):
            p = e + d
            if p < 0 or p + n > len(bits):
                continue
            w = bits[p:p + n]
            for inv in (False, True):
                err = int(np.sum((1 - w if inv else w) != pat))
                if err <= max_err:
                    return dict(name=name, sync_hex=f"{word:0{n // 4}X}", offset_bits=int(d), inverted=inv,
                                bit_errors=err)
    return None


def frame(y, ofs, rate):
    """One burst centred on its FSK carrier at ~8 samples / bit -> dict or None."""
    best = None
    for ph, bits in SIK.demod_bits(y, ofs, rate):
        n = len(bits)
        if n < MIN_PRE + SYNC_BITS:
            continue
        for s, ln in _alt_runs(bits):
            if ln < MIN_PRE or s > max(16, n // 3):
                continue
            if best is None or ln > best[0]:
                best = (ln, s, bits)
            break                                     # the first long run is the preamble
    if best is None:
        return None
    ln, s, bits = best
    e = s + ln
    sync = bits[e:e + SYNC_BITS]
    if len(sync) < SYNC_BITS:
        return None
    out = dict(preamble_bits=int(ln), preamble_last_bit=int(bits[e - 1]),
               sync32_hex=f"{int(''.join(map(str, sync)), 2):08X}", payload_bits=int(len(bits) - e - SYNC_BITS),
               polarity="bit 1 = higher frequency")
    ks = known_sync(bits, e)
    if ks:
        out["known_sync"] = ks
    return out
