#!/usr/bin/env python3
"""Bench: LTE SSS decoder variants (current ZF-like per-subcarrier equaliser vs smoothed-channel MRC).

    .venv/bin/python bench/live_lo/sss_variants.py <sigmf-meta> <f_MHz> [<f_MHz> ...]
"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from lte_probe import Ch
from spectrum_intel.receivers.lte import cellsearch as C
from spectrum_intel.io.sigmf import load

def sss_variant(y, nid2, inst, mode, smooth=7):
    V = C._phy(); sc = V._SC62 % C.N; d = V.lte_pss_freq(C.ROOTS[nid2]); bank = V._SSS_BANK[nid2]
    gaps = {("FDD", "normal"): C.N + C.CP_N, ("TDD", "normal"): 3 * (C.N + C.CP_N) + (C.CP0_N - C.CP_N),
            ("FDD", "extended"): C.N + C.CP_E, ("TDD", "extended"): 3 * (C.N + C.CP_E)}
    res = []
    for (dup, cp), gap in gaps.items():
        acc = np.zeros((2, 168)); used = 0
        for k, p in inst:
            s = p - gap
            if s < 0 or p + C.N > len(y): continue
            Yp = np.fft.fft(y[p:p + C.N])[sc]; Ys = np.fft.fft(y[s:s + C.N])[sc]
            H = Yp / d
            if mode == "zf":
                E = Ys * np.conj(H) / (np.abs(H) ** 2 + np.mean(np.abs(H) ** 2) * 0.05)
            else:
                # sync subcarriers 31 below / 31 above DC: smooth each half separately (DC is not used)
                Hs = np.concatenate([np.convolve(H[:31], np.ones(smooth) / smooth, "same"),
                                     np.convolve(H[31:], np.ones(smooth) / smooth, "same")])
                E = Ys * np.conj(Hs)                                   # MRC: weight by channel SNR
                E = E / (np.sqrt(np.mean(np.abs(E) ** 2)) + 1e-12)
            v = np.real(np.einsum("snk,k->sn", bank, E))
            for par in (0, 1):
                acc[par] += v[0 if (k + par) % 2 == 0 else 1]
            used += 1
        if not used: continue
        flat = np.sort(acc.ravel()); par, nid1 = np.unravel_index(np.argmax(acc), acc.shape)
        med = np.median(flat); noise = 1.4826 * np.median(np.abs(flat - med)) + 1e-12
        res.append(dict(duplex=dup, cp=cp, nid1=int(nid1), z=float((flat[-1] - med) / noise), z2=float((flat[-2] - med) / noise)))
    return max(res, key=lambda r: r["z"] / max(r["z2"], 1e-3))

def main():
    cap = load(sys.argv[1]); ch = Ch(cap.iq, cap.fs, cap.fc)
    for f in map(float, sys.argv[2:]):
        sl = C._Slice(ch, f * 1e6); b = sl.best(np.arange(-30e3, 30.1e3, 5e3))
        fine = sl.best(np.arange(b["df"] - 3e3, b["df"] + 3.01e3, 500.0)); y = sl.shifted(fine["df"])
        inst = C._instances(fine["m"], fine["phase"])
        print(f"{f} MHz: PSS z {fine['z']:.0f} NID2 {fine['nid2']} at {(f*1e6+fine['df'])/1e6:.4f}, {len(inst)} PSS")
        for mode in ("zf", "mrc"):
            r = sss_variant(y, fine["nid2"], inst, mode)
            print(f"   {mode:4s}: {r['duplex']} {r['cp']} NID1 {r['nid1']} PCI {3*r['nid1']+fine['nid2']}  z {r['z']:.2f}  z2 {r['z2']:.2f}"
                  f"  ratio {r['z']/r['z2']:.2f}  decisive {C.sss_decisive(r)}")

if __name__ == "__main__":
    main()
