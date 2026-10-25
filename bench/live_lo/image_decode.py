#!/usr/bin/env python3
"""Bench: analyse a stored LO snapshot as an LO-harmonic image (core/spurs: H3 inverted, H5 not inverted).

    .venv/bin/python bench/live_lo/image_decode.py 1156 [--n 3] [--ppm 10]
n = 3: conj(IQ) analysed as centred on 3 x LO; n = 5: IQ as centred on 5 x LO. Prints the decoded cellular signals.
"""
import argparse, json
from pathlib import Path
import numpy as np
from spectrum_intel.api.engine import Engine

HERE = Path(__file__).resolve().parent

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("tag"); ap.add_argument("--n", type=int, default=3)
    ap.add_argument("--ppm", type=float, default=10)
    a = ap.parse_args()
    cap = json.loads((HERE / "out" / a.tag / "capture.json").read_text()); fs = cap["rate"]
    fc = a.n * cap["freq_mhz"] * 1e6; eng = Engine()
    for k, sn in enumerate(cap["snapshots"]):
        z = np.load(sn["file"]); x = (z[:, 0] + 1j * z[:, 1]).astype(np.complex64) / 2048; del z
        if a.n % 4 == 3:
            x = np.conj(x)
        r = eng.analyze_iq(x, fs, fc, 56e6, receiver_ppm=a.ppm); del x
        for s in r["signals"]:
            if s["technology"] in ("lte", "nr", "gsm", "umts") or s["decode"]["status"] != "not_possible":
                F = s["decode"]["fields"]; mib = F.get("mib") if isinstance(F.get("mib"), dict) else {}
                print(k, round(s["frequency"]["center_hz"] / 1e6, 3), s["label"], s["decode"]["status"],
                      s.get("identifiers"), {k2: mib[k2] for k2 in ("systemFrameNumber",) if k2 in mib})

if __name__ == "__main__":
    main()
