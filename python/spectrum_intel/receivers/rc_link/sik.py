"""
sik.py -- SiK telemetry radio packets (3DR / Holybro / mRo / RFD900-compatible; ArduPilot/SiK firmware, BSD-2) from IQ.

Research (research/drone_transmitters/04_telemetry_remoteid.md, 05_sik_mavlink_details.md), synthetic validation only.
On air (Si443x packet handler, SiK ECC=0, MSB first): preamble 1010... (16 nibbles), sync 0x2D 0xD4, then PN9-whitened
[NETID hi, NETID lo, length, data (MAVLink frames + 2-byte TDM trailer), hardware CRC-16].

Not stated by any reachable source, so measured per packet instead of assumed:
  - GFSK polarity (the sync word is searched in both);
  - the PN9 bit convention (x^9 + x^5 + 1, seed 0x1FF, applied MSB-first or LSB-first per byte: both tried);
  - the hardware CRC seed / reflection / byte order (IBM 0x8005 variants tried; the matching one is reported).
A variant is accepted when the MAVLink frames inside pass their own CRC (X.25 with CRC_EXTRA: an integrity check
independent of the radio) or the radio CRC matches; otherwise only the sync / PHY facts are reported.
Golay mode (ECC=1) and AES-encrypted links are not decoded.
"""
from __future__ import annotations

import numpy as np

from spectrum_intel.core import mavlink as MAV

SYNC = 0x2DD4
RATES = (2e3, 4e3, 8e3, 16e3, 19e3, 24e3, 32e3, 48e3, 64e3, 96e3, 128e3, 192e3, 250e3)


def pn9_bits(n: int, seed: int = 0x1FF) -> np.ndarray:
    """PN9 (x^9 + x^5 + 1) output bits, CC1101 / Si443x style (bit 0 of the register is the output)."""
    reg, out = seed, np.empty(n, np.uint8)
    for i in range(n):
        out[i] = reg & 1
        fb = (reg ^ (reg >> 5)) & 1
        reg = (reg >> 1) | (fb << 8)
    return out


_PN9 = pn9_bits(8 * 260)


def _whiten_bytes(nbytes: int, msb_first: bool) -> np.ndarray:
    w = _PN9[:8 * nbytes].reshape(nbytes, 8)
    if not msb_first:                         # bit i of each PN9 byte = i-th generated bit (CC1101 byte table)
        return (w * (1 << np.arange(8))).sum(1).astype(np.uint8)
    return (w * (1 << np.arange(7, -1, -1))).sum(1).astype(np.uint8)


WHITEN = {"pn9_lsb_first": _whiten_bytes(260, False), "pn9_msb_first": _whiten_bytes(260, True)}


def crc16_variants(data: bytes) -> dict:
    """CRC-16/IBM (0x8005) variants: (refin/refout, init) -> value; xorout 0."""
    out = {}
    for init in (0x0000, 0xFFFF):
        c = init                                              # non-reflected (MSB first)
        for b in data:
            c ^= b << 8
            for _ in range(8):
                c = ((c << 1) ^ 0x8005) & 0xFFFF if c & 0x8000 else (c << 1) & 0xFFFF
        out[f"ibm_msb_init{init:04x}"] = c
        c = init                                              # reflected (ARC / MODBUS family)
        for b in data:
            c ^= b
            for _ in range(8):
                c = (c >> 1) ^ 0xA001 if c & 1 else c >> 1
        out[f"ibm_refl_init{init:04x}"] = c
    return out


def demod_bits(y: np.ndarray, ofs: float, rate: float, sps_phases: int = 8):
    """Hard GFSK decisions at every symbol phase: list of (phase, bits). Discriminator smoothed over half a bit."""
    sps = ofs / rate
    ifr = np.angle(y[1:] * np.conj(y[:-1]))
    k = max(1, int(sps / 2))
    ifs = np.convolve(ifr, np.ones(k) / k, "same")
    thr = float(np.median(ifs))
    out = []
    for ph in np.arange(0, sps, max(sps / sps_phases, 1)):
        idx = (ph + np.arange(int((len(ifs) - ph) / sps)) * sps).astype(int)
        out.append((float(ph), (ifs[idx] > thr).astype(np.uint8)))
    return out


def two_level_threshold(v: np.ndarray, iters: int = 6) -> float:
    """Decision level between the two tones of an FSK statistic, from a 1-D 2-means (started from the 5th / 95th
    percentile midpoint): independent of the ratio of ones to zeros."""
    v = np.asarray(v, float)
    thr = 0.5 * float(np.percentile(v, 5) + np.percentile(v, 95))
    for _ in range(iters):
        lo, hi = v[v < thr], v[v >= thr]
        if not len(lo) or not len(hi):
            break
        thr = 0.5 * float(lo.mean() + hi.mean())
    return thr


def demod_bits_mid(y: np.ndarray, ofs: float, rate: float, sps_phases: int = 8):
    """As demod_bits, with the threshold between the two tones (two_level_threshold of the smoothed discriminator)
    instead of the median: for unbalanced payloads (many 0x00 bytes: toy drones, Hubsan) or windows with silence around
    the packet, where the median sits on one tone."""
    sps = ofs / rate
    ifr = np.angle(np.asarray(y[1:]) * np.conj(np.asarray(y[:-1])))
    k = max(1, int(sps / 2))
    ifs = np.convolve(ifr, np.ones(k) / k, "same")
    thr = two_level_threshold(ifs)
    out = []
    for ph in np.arange(0, sps, max(sps / sps_phases, 1)):
        idx = (ph + np.arange(int((len(ifs) - ph) / sps)) * sps).astype(int)
        out.append((float(ph), (ifs[idx] > thr).astype(np.uint8)))
    return out


def _find_sync(bits: np.ndarray, max_err: int = 1):
    """Positions after a 0x2DD4 sync preceded by >= 8 bits of 1010 preamble, in either polarity: [(pos, inverted)]."""
    pat = np.array([(SYNC >> (15 - i)) & 1 for i in range(16)], np.int8) * 2 - 1
    pre = np.array([(i + 1) % 2 for i in range(8)], np.int8) * 2 - 1        # ...1010 ending in 0 before 0x2D (0010 1101)
    s = bits.astype(np.int8) * 2 - 1
    if len(s) < 40:
        return []
    c = np.correlate(s, pat, "valid")
    hits = []
    for inv, cc in ((False, c), (True, -c)):
        for p in np.flatnonzero(cc >= 16 - 2 * max_err):
            if p >= 8:
                q = s[p - 8:p] * (-1 if inv else 1)
                if abs(int(np.sum(q * pre))) >= 6:            # alternating preamble right before the sync
                    hits.append((int(p) + 16, inv))
    return hits


def decode_packet(bits_after_sync: np.ndarray, inverted: bool):
    """-> dict for the best whitening variant, or None. bits: hard decisions after the sync word."""
    b = bits_after_sync ^ (1 if inverted else 0)
    nby = len(b) // 8
    if nby < 6:
        return None
    raw = np.packbits(b[:8 * nby])
    best = None
    for wname, w in WHITEN.items():
        d = raw ^ w[:nby]
        L = int(d[2])
        if 3 + L + 2 > nby or L < 3:
            continue
        hdr_len_data = bytes(d[:3 + L])
        crc_rx_be = (int(d[3 + L]) << 8) | int(d[4 + L])
        crc_rx_le = (int(d[4 + L]) << 8) | int(d[3 + L])
        crc_match = [k + "_" + o for k, v in crc16_variants(hdr_len_data).items()
                     for o, rx in (("hi_first", crc_rx_be), ("lo_first", crc_rx_le)) if v == rx]
        data = bytes(d[3:3 + L])
        frames = MAV.scan(data[:-2] if L > 2 else data)        # last 2 bytes: TDM trailer
        n_ok = sum(1 for f in frames if f.get("crc_ok"))
        score = 2 * n_ok + (3 if crc_match else 0)
        if score == 0:
            continue
        cand = dict(netid=(int(d[0]) << 8) | int(d[1]), length=L, whitening=wname, hw_crc_variant=crc_match or None,
                    hw_crc_ok=bool(crc_match), mavlink_frames=frames, mavlink_crc_ok=n_ok,
                    tdm_trailer_hex=bytes(d[3 + L - 2:3 + L]).hex() if L >= 2 else None, data_hex=data.hex())
        if best is None or score > best[0]:
            best = (score, cand)
    return best[1] if best else None


def sik_try_decode(y: np.ndarray, ofs: float, rate: float):
    """One burst (complex baseband centred on the FSK carrier) -> decoded packet dict, or a PHY-only dict when the
    sync word is found but nothing verifies, or None."""
    phy = None
    for ph, bits in demod_bits(y, ofs, rate):
        for pos, inv in _find_sync(bits):
            r = decode_packet(bits[pos:], inv)
            if r:
                r.update(sync="0x2DD4", polarity_inverted=inv, air_rate_bps=rate)
                return r
            phy = phy or dict(sync="0x2DD4", polarity_inverted=inv, air_rate_bps=rate, decoded=False)
    return phy
