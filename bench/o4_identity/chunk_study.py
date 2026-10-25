#!/usr/bin/env python3
"""Bench: how short can a packet chunk be and still give the drone's crystal ppm?

Packets are cut with signal-vault's own DDC (signal_vault.dsp.ddc_) at fs = 61.44 MSps, from the O4 bursts the
library found (out/hopping, out/debug_detection), at two packet rates: the box width (PacketRate.BANDWIDTH) and the
O4-native integer decimation (15.36 / 30.72 / 61.44 MSps). From each packet, non-overlapping chunks of 3 symbols,
1 slot and 2 slots are measured with chunk.chunk_ppm and compared with the full-stretch value (FINDINGS.md).

    cd bench/o4_identity && nice -n 19 ../../.venv/bin/python chunk_study.py
"""
import json
from pathlib import Path

import numpy as np
from signal_vault.dsp import ddc_

import chunk
import o4

FS = 61.44e6
REPO = Path(__file__).resolve().parents[2]
HOP = Path("/Users/srikanth/Documents/modules/signal-vault/signal-vault-data/recordings/hopping")
DBG = Path("/Users/srikanth/Documents/modules/signal-vault/signal-vault-data/recordings/debug_detection")
H, D = "test_signal-20261016_065855", "test_signal_debug-20261016_071042"
NOMINAL = {"2.4G 20MHz c0": 2415.5e6, "2.4G 20MHz c1": 2415.5e6, "2.4G 20MHz c2": 2415.5e6,
           "5.8G 20MHz c3": 5823.5e6, "5.73G 10MHz": 5730.5e6}
SIZES = {"3 symbols": 110e-6, "1 slot": 0.5e-3, "2 slots": 1.0e-3}
MAX_PER_PACKET = 6

# (name, root, prefix, fc MHz, chunk, report dir, full-stretch ppm, template)
LINKS = [("2.4G 20MHz c0", HOP, H, 2428, 0, "hopping", 1.459, "B"),
         ("2.4G 20MHz c1", HOP, H, 2428, 1, "hopping", 1.462, "A"),
         ("2.4G 20MHz c2", HOP, H, 2428, 2, "hopping", 1.423, "A"),
         ("5.8G 20MHz c3", HOP, H, 5840, 3, "hopping", 1.536, "A"),
         ("5.73G 10MHz", DBG, D, 5728, 0, "debug_detection", 1.512, "A")]


def path(root, prefix, fq, ck):
    d = f"{prefix}_{fq}MHz_ch2_chunk{ck:04d}"
    return root / d / f"{d}.sigmf-data", REPO / "out" / ("hopping" if root == HOP else "debug_detection") / f"{d}.json"


def main_box(report, fc):
    """The library's widest in-band O4 video link of the capture: its bursts and band."""
    r = json.loads(report.read_text())
    cands = [s for s in r["signals"] if s["technology"] == "dji_ocusync" and abs(s["frequency"]["offset_hz"]) < FS / 2]
    s = max(cands, key=lambda s: s["frequency"]["bandwidth_hz"])
    return s["frequency"]["low_hz"], s["frequency"]["high_hz"], [(b["start_s"], b["end_s"]) for b in s["time"]["bursts"]]


def main():
    x = o4.load_cf32(path(HOP, H, 2428, 0)[0])
    T = {"A": chunk.build_template(x, 2428e6, 2415.49e6, 18e6, 0, 1_000_000)}
    x = o4.load_cf32(path(HOP, H, 2428, 2)[0])
    T["B"] = chunk.build_template(x, 2428e6, 2415.49e6, 18e6, 100_000, 1_000_000)
    print(f"templates: A {len(T['A']['rows'])} rows (2.4G c0), B {len(T['B']['rows'])} rows (2.4G c2)\n")
    table = {}
    for name, root, prefix, fq, ck, _, truth, tk in LINKS:
        p, rep = path(root, prefix, fq, ck)
        x = o4.load_cf32(p); fc = fq * 1e6
        f_lo, f_hi, bursts = main_box(rep, fc)
        lo_n, hi_n = int(f_lo - fc), int(f_hi - fc)
        mix = fc + (lo_n + hi_n) // 2
        for mode in ("bandwidth", "native"):
            desired = float(hi_n - lo_n) if mode == "bandwidth" else chunk.native_rate(hi_n - lo_n)
            for b0, b1 in bursts:
                s0, s1 = int(b0 * FS), int(b1 * FS)
                if s1 - s0 < int(150e-6 * FS):
                    continue
                pk, _, rate = ddc_(x, (s0, s1, lo_n, hi_n), FS, desired)
                n = len(pk); a, z = int(0.1 * n), int(0.9 * n)
                for sz, dur in SIZES.items():
                    L = int(round(dur * rate))
                    if z - a < L:
                        continue
                    k = min(MAX_PER_PACKET, (z - a) // L)
                    for j in range(k):
                        c0 = a + j * ((z - a - L) // max(k - 1, 1)) if k > 1 else a
                        res = chunk.chunk_ppm(pk[c0:c0 + L], rate, mix, hi_n - lo_n, T[tk])
                        table.setdefault((name, mode, sz), []).append((res, truth, L, rate))
        del x
    print(f"{'link':16s} {'packet rate':12s} {'chunk':10s} {'samples':>8s} {'chunks':>6s} {'n right':>8s} "
          f"{'ppm err median':>14s} {'ppm err sd':>10s} {'95% |err|':>9s}")
    for name, *_ in LINKS:
        for mode in ("bandwidth", "native"):
            for sz in SIZES:
                v = table.get((name, mode, sz))
                if not v:
                    continue
                L = v[0][2]; rate = v[0][3]
                rej = sum(1 for r, *_ in v if not r.get("ok"))
                # tracking mode: n known from the drone's track, every chunk with symbols gives a fine-CFO ppm
                trk = []
                for r, t, *_ in v:
                    if "cfo" not in r and "carrier" not in r:
                        continue
                    nom = NOMINAL[name]; base = r["carrier"] - chunk.SCS * r["n"]          # mix + cfo
                    n_t = round((nom * (1 + t * 1e-6) - base) / chunk.SCS)
                    trk.append((base + chunk.SCS * n_t - nom) / nom * 1e6 - t)
                trk = np.array(trk)
                wrong = sum(1 for r, t, *_ in v if r.get("ok") and abs(r["ppm"] - t) >= 1.0)
                err = np.array([r["ppm"] - t for r, t, *_ in v if r.get("ok") and abs(r["ppm"] - t) < 1.0])
                right = len(err)
                s = (f"{np.median(err):+14.3f} {np.std(err):10.3f} {np.percentile(np.abs(err), 95):9.3f}"
                     if right > 1 else f"{'-':>14s} {'-':>10s} {'-':>9s}")
                print(f"{name:16s} {mode + f' {rate/1e6:.2f}':12s} {sz:10s} {L:8d} {len(v):6d} "
                      f"{right:4d}/{len(v):<3d} {s}   rejected {rej:3d}  WRONG n {wrong:3d}   tracking: {len(trk)} chunks sd {np.std(trk):.3f} 95% {np.percentile(np.abs(trk), 95):.3f}", flush=True)


if __name__ == "__main__":
    main()
