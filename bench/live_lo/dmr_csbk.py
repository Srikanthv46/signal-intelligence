#!/usr/bin/env python3
"""Bench: decode DMR CSBKs (TS 102 361-1 B.1.1 BPTC(196,96), B.3.12 CRC-CCITT, CSBK mask 0xA5A5) from a dwell.

    .venv/bin/python bench/live_lo/dmr_csbk.py d860 859.8217
"""
import sys
from collections import Counter
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from dmr_probe import fm48k, demod, bits48, SYNC

def bptc_196_96(r):
    m = np.array([r[(13 * i) % 196] for i in range(196)], np.int8)    # deinterleave: matrix bit i sent at 13i mod 196
    M = m[1:].reshape(13, 15)
    # Hamming(15,11) row correction (single error) using syndrome over the 4 parity bits
    H15 = [(0, 1, 2, 3, 5, 7, 8), (1, 2, 3, 4, 6, 8, 9), (2, 3, 4, 5, 7, 9, 10), (0, 1, 2, 4, 6, 7, 10)]
    for row in range(9):
        syn = tuple(int(np.sum(M[row, list(c)]) % 2 != M[row, 11 + k]) for k, c in enumerate(H15))
        if any(syn):
            for col in range(11):
                if tuple(int(col in c) for c in H15) == syn:
                    M[row, col] ^= 1; break
    return np.concatenate([M[0, 3:11], M[1:9, :11].ravel()])          # 8 + 88 = 96 info bits

def crc_ccitt(bits, init=0x0000):
    r = init
    for b in bits:
        fb = ((r >> 15) & 1) ^ int(b); r = (r << 1) & 0xFFFF
        if fb: r ^= 0x1021
    return r

def main():
    name, f = sys.argv[1], float(sys.argv[2])
    d, fo = fm48k(name, f); bits, _, _ = demod(d, fo)
    b = 2 * bits.astype(np.float32) - 1; p = 2 * bits48(SYNC["BS data"]).astype(np.float32) - 1
    hits = np.where((48 - np.correlate(b, p, "valid")) / 2 <= 4)[0]
    ok, seen = 0, Counter()
    for h in hits:
        if h < 108 or h + 156 > len(bits):
            continue
        info = bptc_196_96(np.concatenate([bits[h - 108:h - 10], bits[h + 58:h + 156]]))
        crc = int("".join(map(str, info[80:])), 2)
        calc = crc_ccitt(info[:80]) ^ 0xFFFF ^ 0xA5A5
        if calc != crc:
            continue
        ok += 1
        by = [int("".join(map(str, info[8 * i:8 * i + 8])), 2) for i in range(10)]
        lb, pf, op, fid = by[0] >> 7, (by[0] >> 6) & 1, by[0] & 0x3F, by[1]
        seen[(op, fid, bytes(by[2:10]).hex())] += 1
    print(f"{f} MHz: {len(hits)} bursts, {ok} CSBK CRC ok")
    for (op, fid, data), n in seen.most_common(12):
        print(f"  x{n:2d}  CSBKO {op:2d} (0x{op:02X})  FID 0x{fid:02X}  data {data}")

if __name__ == "__main__":
    main()


# ---- per-burst refined demodulator (RRC matched filter, sync-aided timing/levels)
def rrc(beta, sps, span=6):
    t = np.arange(-span * sps, span * sps + 1) / sps
    h = np.sinc(t) * np.cos(np.pi * beta * t) / (1 - (2 * beta * t) ** 2 + 1e-9)
    return h / h.sum()

SYNC_SYM = None

def sync_symbols(v):
    b = bits48(v); m = {(0, 1): 3, (0, 0): 1, (1, 0): -1, (1, 1): -3}
    return np.array([m[(b[2 * i], b[2 * i + 1])] for i in range(24)], np.float32)

def bursts_refined(name, f, sync="BS data", sr=4800.0):
    d, fo = fm48k(name, f)
    sps = fo / sr
    z = np.convolve(d, rrc(0.2, sps), "same")
    ref = sync_symbols(SYNC[sync])
    # coarse sync positions in samples from the simple demod
    bits, _, _ = demod(d, fo)
    b = 2 * bits.astype(np.float32) - 1; p = 2 * bits48(SYNC[sync]).astype(np.float32) - 1
    hits = np.where((48 - np.correlate(b, p, "valid")) / 2 <= 6)[0]
    out = []
    for h in hits:
        s0 = (h // 2) * sps                         # sync start in samples (approx)
        best = None
        for dt in np.arange(-1.5 * sps, 1.5 * sps, sps / 20):
            idx = (s0 + dt + sps * np.arange(24)).astype(int)
            if idx.min() < 0 or idx.max() >= len(z): continue
            v = z[idx]; c = np.dot(v - v.mean(), ref) / (np.linalg.norm(v - v.mean()) * np.linalg.norm(ref) + 1e-9)
            if best is None or c > best[0]: best = (c, dt, v)
        if best is None: continue
        c, dt, v = best
        A = np.polyfit(ref, v, 1)                    # v = gain * ref + dc
        start = s0 + dt - 54 * sps                   # burst = 54 symbols payload+slot type, 24 sync, 54 after
        idx = (start + sps * np.arange(132)).astype(int)
        if idx.min() < 0 or idx.max() >= len(z): continue
        sym = (z[idx] - A[1]) / A[0]
        q = np.clip(np.round((sym + 3) / 2) * 2 - 3, -3, 3).astype(int)
        mp = {3: (0, 1), 1: (0, 0), -1: (1, 0), -3: (1, 1)}
        bb = np.array([x for s in q for x in mp[s]], np.int8)
        out.append((h, c, bb, float(np.mean((sym - q) ** 2))))
    return out

if __name__ == "__main__" and len(sys.argv) > 3 and sys.argv[3] == "--refined":
    pass
