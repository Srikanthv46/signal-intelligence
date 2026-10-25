#!/usr/bin/env python3
"""Bench: decode the TETRA BSCH (SYNC PDU + D-MLE-SYNC: MCC, MNC, colour code, TN/FN/MN) from every
synchronisation burst. Implemented from ETSI EN 300 392-2 (8.2 coding, 9.4 bursts, 21.4.4 SYNC PDU).

    .venv/bin/python bench/live_lo/tetra_bsch.py tetra390 390.0297
"""
import sys
from collections import Counter
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from tetra_probe import load
from tetra_demod import demod, search, Y_SYNC

# --- 8.2.5 scrambling: 32-stage LFSR, taps 32 26 23 22 16 12 11 10 8 7 5 4 2 1; init p(-31)=p(-30)=1, p(-29..0)=e
TAPS = (32, 26, 23, 22, 16, 12, 11, 10, 8, 7, 5, 4, 2, 1)

def scramb_seq(n, e_bits=None):
    p = [1, 1] + list(e_bits if e_bits is not None else [0] * 30)    # p(-31) .. p(0)
    out = []
    for _ in range(n):
        b = 0
        for t in TAPS:
            b ^= p[-t]
        p.append(b); out.append(b)
    return np.array(out, np.int8)

# --- 8.2.4.1 block interleaver: b4(k) = b3(i), k = 1 + (a*i mod K)
def deinterleave(b4, K=120, a=11):
    b3 = np.empty(K, b4.dtype)
    for i in range(1, K + 1):
        b3[i - 1] = b4[(a * i) % K]
    return b3

# --- 8.2.3.1 RCPC: mother rate 1/4, K=5; rate 2/3 puncturing t=3, period 8, P = 1, 2, 5
G = ((1, 1, 0, 0, 1), (1, 0, 1, 1, 1), (1, 1, 1, 0, 1), (1, 1, 0, 1, 1))     # g_i,j for j = 0..4 (D^0..D^4)
P23 = (1, 2, 5)

def depuncture(b3, n_in):
    """-> soft mother-code stream of length 4*n_in (0 = erased)."""
    v = np.zeros(4 * n_in, np.float32)
    t = 3
    for j in range(1, len(b3) + 1):
        k = 8 * ((j - 1) // t) + P23[j - t * ((j - 1) // t) - 1]
        v[k - 1] = b3[j - 1]
    return v

def viterbi(v, n):
    """Soft (+1 = bit 1, -1 = bit 0, 0 = erasure) Viterbi for the 16-state mother code; returns n bits."""
    S = 16
    out_bits = np.zeros((S, 2, 4), np.float32)               # state = last 4 input bits (b(k-1)..b(k-4))
    for s in range(S):
        for u in (0, 1):
            reg = [u] + [(s >> i) & 1 for i in range(4)]      # b(k), b(k-1), ..., b(k-4)
            for i in range(4):
                out_bits[s, u, i] = 2 * (sum(G[i][j] * reg[j] for j in range(5)) % 2) - 1
    nxt = np.array([[((s << 1) | u) & 0xF for u in (0, 1)] for s in range(S)])
    pm = np.full(S, -1e9); pm[0] = 0
    hist = np.zeros((n, S), np.int16); hbit = np.zeros((n, S), np.int8)
    for k in range(n):
        r = v[4 * k:4 * k + 4]; npm = np.full(S, -1e9)
        for s in range(S):
            for u in (0, 1):
                m = pm[s] + float(np.dot(out_bits[s, u], r)); t = nxt[s, u]
                if m > npm[t]:
                    npm[t] = m; hist[k, t] = s; hbit[k, t] = u
        pm = npm
    s = 0 if pm[0] > -1e8 else int(np.argmax(pm))           # tail bits drive the encoder back to 0
    bits = np.zeros(n, np.int8)
    for k in range(n - 1, -1, -1):
        bits[k] = hbit[k, s]; s = hist[k, s]
    return bits

# --- 8.2.2 CRC: CRC-CCITT X^16+X^12+X^5+1, ones-complement
def crc16(bits):
    r = 0xFFFF
    for b in bits:
        fb = ((r >> 15) & 1) ^ int(b)
        r = (r << 1) & 0xFFFF
        if fb:
            r ^= 0x1021
    return r ^ 0xFFFF

def field(bits, pos, n):
    return int("".join(map(str, bits[pos:pos + n])), 2), pos + n

def parse_sync(b):
    f, p = {}, 0
    for name, n in (("system_code", 4), ("colour_code", 6), ("timeslot", 2), ("frame", 5), ("multiframe", 6),
                    ("sharing_mode", 2), ("ts_reserved_frames", 3), ("u_plane_dtx", 1), ("frame18_ext", 1),
                    ("reserved", 1), ("mcc", 10), ("mnc", 14), ("neighbour_cell_bcast", 2),
                    ("cell_service_level", 2), ("late_entry", 1)):
        f[name], p = field(b, p, n)
    f["timeslot"] += 1
    return f

def soft_bits(dd):
    """Soft dibits from the differential phasor: +ve = bit 1. b1 = 1 for negative phase, b2 = 1 for |phase| > pi/2."""
    a = dd / (np.abs(dd) + 1e-12) * np.sqrt(np.abs(dd) / (np.median(np.abs(dd)) + 1e-12))
    s = np.empty(2 * len(dd), np.float32); s[0::2], s[1::2] = -a.imag, -a.real
    return s

def decode_bsch(type5_soft):
    """type5_soft: 120 soft bits (+ve = 1)."""
    sc = scramb_seq(120)
    b4 = np.where(sc == 1, -type5_soft, type5_soft)          # BSCH uses e = 0 (colour code not yet known)
    b3 = deinterleave(b4)
    v = depuncture(b3, 80)
    b2 = viterbi(v, 80)
    info, crc, tail = b2[:60], b2[60:76], b2[76:]
    ok = crc16(info) == int("".join(map(str, crc)), 2)
    return ok, info

def main():
    name, f = sys.argv[1], float(sys.argv[2])
    x, m = load(name)
    bits, dd = demod(x, m, f)[:2]
    soft = soft_bits(dd)
    hits = [h for h, e in search(bits, Y_SYNC, 4) if h >= 120]
    good = Counter()
    for h in hits:
        ok, info = decode_bsch(soft[h - 120:h])
        if ok:
            s = parse_sync(info)
            good[(s["mcc"], s["mnc"], s["colour_code"], s["system_code"])] += 1
            print(f"  bit {h:6d} CRC ok  TN{s['timeslot']} FN{s['frame']:2d} MN{s['multiframe']:2d}  MCC {s['mcc']} MNC {s['mnc']} "
                  f"CC {s['colour_code']}  sys {s['system_code']}  share {s['sharing_mode']}  csl {s['cell_service_level']} "
                  f"nbr {s['neighbour_cell_bcast']}  late {s['late_entry']}")
        else:
            print(f"  bit {h:6d} CRC fail")
    print(f"{f} MHz: {len(hits)} sync bursts, {sum(good.values())} BSCH CRC ok -> {dict(good)}")

if __name__ == "__main__":
    main()
