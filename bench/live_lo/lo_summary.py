#!/usr/bin/env python3
"""Bench: one line per signal of an LO, grouped across its snapshots (from lo_step's out/<MHz>/snapNN.json).

    .venv/bin/python bench/live_lo/lo_summary.py 1156 [--min-snr 6]
Columns: centre MHz, width kHz, snapshots seen / taken, max SNR, technology, label, decode status (best), reason.
"""
import argparse, json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RANK = {"success": 2, "partial": 1, "not_possible": 0}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("tag"); ap.add_argument("--min-snr", type=float, default=6)
    a = ap.parse_args()
    files = sorted((HERE / "out" / a.tag).glob("snap*.json"))
    groups = []
    for k, fn in enumerate(files):
        r = json.loads(fn.read_text())
        for s in r.get("signals", []):
            snr = s.get("snr_db") or 0
            if snr < a.min_snr and s["decode"]["status"] == "not_possible":
                continue
            f, bw = s["frequency"]["center_hz"] / 1e6, s["frequency"]["bandwidth_hz"] / 1e3
            g = next((g for g in groups if g["tech"] == s["technology"] and abs(g["f"] - f) * 1e3 < max(bw, g["bw"], 50) / 2), None)
            if g is None:
                g = dict(f=f, bw=bw, tech=s["technology"], label=s["label"], snaps=set(), snr=snr, st=s["decode"]["status"],
                         why=s["decode"]["reason"], ids=s.get("identifiers") or {})
                groups.append(g)
            g["snaps"].add(k); g["snr"] = max(g["snr"], snr)
            if RANK[s["decode"]["status"]] > RANK[g["st"]]:
                g.update(st=s["decode"]["status"], label=s["label"], why=s["decode"]["reason"], ids=s.get("identifiers") or {})
    for g in sorted(groups, key=lambda g: g["f"]):
        ids = f" {json.dumps(g['ids'])}" if g["ids"] else ""
        print(f"{g['f']:10.3f} {g['bw']:8.1f} {len(g['snaps'])}/{len(files)} {g['snr']:5.1f}  {g['tech']:<14} "
              f"{g['label'][:60]:<60} {g['st']:<12}{ids} | {g['why'][:70]}")


if __name__ == "__main__":
    main()
