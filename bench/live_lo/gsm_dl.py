#!/usr/bin/env python3
"""Bench: GSM downlink scan of a stored LO snapshot on the 200 kHz raster (GSM-900 / DCS-1800 downlink).

    .venv/bin/python bench/live_lo/gsm_dl.py 940 [--band 900|1800] [--snaps 0 1 2] [--min-db 10]
Per ARFCN with power >= --min-db above the scan's 10th-percentile channel: FCCH (library find_fcch), the SCH one TDMA
frame later (library demod_sch + Viterbi + CRC-10) -> BSIC (NCC, BCC) and the frame number; the residual
FCCH tone offset gives the receiver error in ppm.
"""
import argparse, json
from pathlib import Path
import numpy as np

from spectrum_intel.core.dsp import Channelizer
from spectrum_intel.receivers.gsm.sch import OFS, SPS, crc10_ok, demod_sch, find_fcch, sch_fields, viterbi

HERE = Path(__file__).resolve().parent
BANDS = {"900": [(n, 935e6 + 0.2e6 * n) for n in range(1, 125)],
         "1800": [(n, 1805.2e6 + 0.2e6 * (n - 512)) for n in range(512, 886)]}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("tag"); ap.add_argument("--band", default="900")
    ap.add_argument("--snaps", type=int, nargs="+", default=[0, 1, 2]); ap.add_argument("--min-db", type=float, default=10)
    ap.add_argument("--ppm", type=float, default=7.4)
    a = ap.parse_args()
    cap = json.loads((HERE / "out" / a.tag / "capture.json").read_text())
    fs, fc = cap["rate"], cap["freq_mhz"] * 1e6
    chans = [(n, f) for n, f in BANDS[a.band] if abs(f - fc) < 28e6]
    for k in a.snaps:
        z = np.load(cap["snapshots"][k]["file"])
        ch = Channelizer(z[:, 0].astype(np.float32) + 1j * z[:, 1].astype(np.float32), fs, fc)
        del z
        ys = {n: ch.slice(f * (1 + a.ppm * 1e-6), 110e3, OFS)[0] for n, f in chans}
        pw = {n: 10 * np.log10(np.mean(np.abs(y) ** 2) + 1e-30) for n, y in ys.items()}
        floor = float(np.percentile(list(pw.values()), 10))
        print(f"== snap{k}: {len(chans)} ARFCNs, 10th-percentile channel {floor:.1f} dB")
        for n, f in chans:
            if pw[n] - floor < a.min_db:
                continue
            y = ys[n]; line = f"  ARFCN {n:4d} {f / 1e6:8.1f} MHz  {pw[n] - floor:5.1f} dB"
            fcch = find_fcch(y)
            if fcch is None:
                print(line + "  no FCCH"); continue
            s0, cfo, pur = fcch
            line += f"  FCCH {s0 / OFS * 1e3:6.2f} ms purity {pur:.2f} ppm {a.ppm + cfo / f * 1e6:+.2f}"
            y2 = y * np.exp(-2j * np.pi * cfo * np.arange(len(y)) / OFS)
            d = demod_sch(y2, s0 + int(round(1250 * SPS)))
            if d is None:
                print(line + "  SCH: no training sequence"); continue
            bits = viterbi(d[0])
            if crc10_ok(bits[:25], bits[25:35]):
                sf = sch_fields(bits[:25])
                line += f"  SCH OK: BSIC {sf['ncc']}{sf['bcc']} (NCC {sf['ncc']} BCC {sf['bcc']}) FN {sf['frame_number']}"
            else:
                line += f"  SCH CRC failed (TS match {d[1]['ts_match']}/64)"
            print(line)
        del ch, ys


if __name__ == "__main__":
    main()
