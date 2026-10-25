#!/usr/bin/env python3
"""Bench: LTE PSS/SSS probe on a stored LO snapshot, using the library's cell-search pieces.

    .venv/bin/python bench/live_lo/lte_probe.py 364 369.5 372.0 [--snap 1] [--step 0.1] [--span 60]
"""
import argparse, json
from pathlib import Path
import numpy as np
from scipy import signal as ss
from spectrum_intel.receivers.lte import cellsearch as C

HERE = Path(__file__).resolve().parent

class Ch:
    def __init__(self, x, fs, fc): self.x, self.fs, self.fc = x, fs, fc
    def slice(self, f, bw, ofs):
        n = np.arange(len(self.x))
        y = self.x * np.exp(-2j * np.pi * (f - self.fc) / self.fs * n).astype(np.complex64)
        up, down = 1, int(round(self.fs / ofs))            # 61.44 / 3.84 = 16
        return ss.resample_poly(y, up, down).astype(np.complex64), ofs

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("freq"); ap.add_argument("lo", type=float); ap.add_argument("hi", type=float)
    ap.add_argument("--snap", type=int, default=1); ap.add_argument("--step", type=float, default=0.1); ap.add_argument("--span", type=float, default=60)
    a = ap.parse_args()
    if a.freq.endswith(".sigmf-meta"):
        from spectrum_intel.io.sigmf import load
        cap = load(a.freq); ch = Ch(cap.iq, cap.fs, cap.fc)
    else:
        c = json.loads((HERE / "out" / a.freq / "capture.json").read_text()); s = np.load(c["snapshots"][a.snap]["file"])
        ch = Ch(s[:, 0].astype(np.float32) + 1j * s[:, 1].astype(np.float32), c["rate"], c["freq_mhz"] * 1e6)
    best = []
    for f in np.arange(a.lo, a.hi + 1e-9, a.step) * 1e6:
        sl = C._Slice(ch, f); b = sl.best(np.arange(-a.span, a.span + 0.1, 5) * 1e3)
        best.append((b["z"], f / 1e6 + b["df"] / 1e6, b["nid2"], b["iso"]))
        if b["z"] >= C.Z_MIN:
            fine = sl.best(np.arange(b["df"] - 3e3, b["df"] + 3.01e3, 500.0))
            inst = C._instances(fine["m"], fine["phase"])
            sss = C.sss_decode(sl.shifted(fine["df"]), fine["nid2"], inst) if len(inst) >= 2 else None
            pci = 3 * sss["nid1"] + fine["nid2"] if sss else None
            print(f"{(f + fine['df'])/1e6:10.4f} MHz  PSS z={fine['z']:.1f} iso={fine['iso']:.1f} NID2={fine['nid2']} "
                  f"SSS={sss and {k: (round(v, 2) if isinstance(v, float) else v) for k, v in sss.items()}} "
                  f"decisive={C.sss_decisive(sss)} PCI={pci}")
    best.sort(reverse=True); print("top z:", [(round(z, 1), round(f, 4), n) for z, f, n, _ in best[:5]])

if __name__ == "__main__":
    main()
