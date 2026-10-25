#!/usr/bin/env python3
"""Bench: real signal or LO-harmonic image? Compare two captures whose LOs differ by delta.

A mixer driven by a square-ish LO also converts F to baseband through LO harmonics n:
  n=1 real: shown at F                 (display shift 0 per delta)
  n=3 inverted: d = 3LO - F            (display shift +4 delta)
  n=5:          d = F - 5LO            (display shift -4 delta)
  n=7 inverted: d = 7LO - F            (display shift +8 delta)
Each segment with structure is assigned the shift that best aligns the two spectra; the image's true source
frequency follows from that.

    .venv/bin/python bench/live_lo/harmonic_map.py 364 365 [--seg 2]
"""
import argparse, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from overlap import mean_psd

HYP = {0: (1, +1), 4: (3, -1), -4: (5, +1), 8: (7, -1), -8: (9, +1)}

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("a"); ap.add_argument("b"); ap.add_argument("--seg", type=float, default=2.0)
    a = ap.parse_args()
    fa, ca, pa = mean_psd(a.a); fb, cb, pb = mean_psd(a.b); delta = cb - ca
    pa -= np.median(pa); pb -= np.median(pb)
    lo, hi = max(ca + fa[0], cb + fb[0]) + 3e6, min(ca + fa[-1], cb + fb[-1]) - 3e6
    grid = np.arange(lo, hi, 25e3); A = np.interp(grid, ca + fa, pa)
    print(f"LO {ca/1e6:g} vs {cb/1e6:g} (delta {delta/1e6:+g} MHz)\n  segment MHz      peak dB  best   corr   (real-shift corr)  -> source")
    for s0 in np.arange(lo, hi - a.seg * 1e6 + 1, a.seg * 1e6):
        sel = (grid >= s0) & (grid < s0 + a.seg * 1e6)
        if A[sel].max() < 3 or A[sel].std() < 0.5:
            continue
        sc = {}
        for k in HYP:
            B = np.interp(grid[sel] + k * delta, cb + fb, pb)
            sc[k] = np.corrcoef(A[sel], B)[0, 1]
        k = max(sc, key=sc.get); n, sgn = HYP[k]
        mid = s0 + a.seg * 0.5e6; d = mid - ca
        src = mid if n == 1 else (n * ca + d if sgn > 0 else n * ca - d)
        tag = "REAL" if n == 1 else f"H{n} image of {src/1e6:8.1f} MHz"
        print(f"  {s0/1e6:7.1f}-{(s0+a.seg*1e6)/1e6:7.1f}  {A[sel].max():6.1f}  {k:+3d}d  {sc[k]:5.2f}   ({sc[0]:5.2f})   {tag}"
              + ("   [weak evidence]" if sc[k] < 0.5 or sc[k] - sc[0] < 0.15 and n != 1 else ""))

if __name__ == "__main__":
    main()
