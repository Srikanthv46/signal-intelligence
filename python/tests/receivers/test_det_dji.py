#!/usr/bin/env python3
"""
test_det_dji.py -- detection / false-alarm / SNR-sweep evaluation of det_dji.

For every DJI capture:
  1. run det_dji.detect on the original -> reference emissions (per label);
     DroneID is checked against the capture's JSON boxes where they exist
  2. add white noise (drone_common.add_awgn) lowering every signal by d dB,
     d in 0, 5, ..., 30, and re-run
  3. per label: detections matching a reference emission (time overlap, same
     band), false alarms (DJI detections matching none), and the in-band SNR of
     a reference burst measured on the noisy capture:
        snr_meas  drone_common.measure_snr (per-bin floor = time percentile; reads
                  low for a link that is on > 50 % of the time, e.g. video)
        snr_therm same box over the receiver thermal floor (det_dji._SNR)
  4. min SNR per label = lowest SNR of the contiguous run of levels (from 0 dB
     added) where >= 50 % of the reference emissions are still found with no
     false alarm of that family
Pure noise and a Herelink capture (LTE-like 15 kHz OFDM, not DJI) must give zero
DJI detections.

  python test_det_dji.py [--levels 0,5,...] [--only mini2,far] [--json out.json]
"""
from __future__ import annotations
import argparse, json, os, sys, time
from collections import Counter, defaultdict
import os
import numpy as np

from spectrum_intel.core import common as C
from spectrum_intel.receivers.dji import detector as D
from _data import DATA as HERE, OUT
from _data import DRONEID_TV as TV
CAPTURES = {
    "mini2": "DJI_Mini_2_20251119-165139/D_20251119-165139.sigmf-meta",
    "mini4pro": "DJI_Mini_Pro_4_20251119-172214/D_20251119-172214.sigmf-meta",
    "phantom4": "DJI_Phantom_4_20251119-182450/D_20251119-182450.sigmf-meta",
    "mini2_2km": "D_20260916-171343/D_20260916-171343.sigmf-meta",
    "tv_2p4_122710": os.path.join(TV, "2P4/D_20260529-122710/D_20260529-122710.sigmf-meta"),
    "tv_5p8_125115": os.path.join(TV, "5P8/D_20260529-125115/D_20260529-125115.sigmf-meta"),
}
NEGATIVE = {"herelink": "Herelink_20251125-110453/D_20251125-110453.sigmf-meta"}
LEVELS = [0, 5, 10, 15, 20, 25, 30]


def _match(e, refs):
    for r in refs:
        ov = min(e.end_us, r.end_us) - max(e.start_us, r.start_us)
        if ov > 0.3 * min(e.end_us - e.start_us, r.end_us - r.start_us) and \
                abs(e.center_mhz - r.center_mhz) < max(1.0, (r.f_hi_mhz - r.f_lo_mhz) / 2):
            return r
    return None


def evaluate(name, path, levels, log=print):
    x, fs, fc, meta = C.load_capture(os.path.join(HERE, path))
    t = time.time()
    ref = D.detect(x, fs, fc, {})
    log(f"\n=== {name}  fc {fc / 1e6:g} MHz: baseline {len(ref)} DJI emissions "
        f"{dict(Counter(e.label for e in ref))} ({time.time() - t:.1f} s)")
    gt = None
    jp = os.path.join(HERE, path).replace(".sigmf-meta", ".json")
    if os.path.exists(jp):
        j = json.load(open(jp))
        gt = [b for v in j["data"].get("drone_id", {}).values() for b in v]
        did = [e for e in ref if e.label == "droneid"]
        hit = sum(1 for b in gt if any(abs(e.start_us / 1e6 - b["start_time"]) < 100e-6 for e in did))
        log(f"    DroneID JSON boxes {len(gt)}, detected {hit}, decoded "
            f"{[e.features.get('serial') for e in did]}")
    by = defaultdict(list)
    for e in ref:
        by[e.label].append(e)
    # reference burst per label: the median-SNR one (a typical burst, not the best)
    refb = {lab: sorted(es, key=lambda e: e.snr_db or -99)[len(es) // 2] for lab, es in by.items()}
    therm = None
    rows = {lab: [] for lab in by}
    fa_any = []
    for d in levels:
        xn = x if d == 0 else C.add_awgn(x, d, np.random.default_rng(d))
        em = D.detect(xn, fs, fc, {"decode_droneid": True})
        therm = D._SNR(xn, fs, fc)
        fa = [e for e in em if _match(e, ref) is None]
        fa_any.append((d, len(fa), dict(Counter(e.label for e in fa))))
        for lab, rs in by.items():
            got = [e for e in em if e.label == lab and _match(e, rs) is not None]
            fa_lab = [e for e in fa if e.family == rs[0].family]
            rb = refb[lab]
            bw = rb.f_hi_mhz - rb.f_lo_mhz
            box = (rb.start_us, rb.end_us, rb.center_mhz - bw / 2.4, rb.center_mhz + bw / 2.4)
            sm = C.measure_snr(xn, fs, fc, *box)
            st = therm(*box)
            dec = sum(1 for e in got if e.decoded)
            rows[lab].append(dict(added_db=d, found=len({id(_match(e, rs)) for e in got}), ref=len(rs),
                                  fa=len(fa_lab), snr_meas=round(sm, 1), snr_therm=st, decoded=dec))
        log(f"  +{d:2d} dB: {dict(Counter(e.label for e in em))}  false alarms {len(fa)} "
            f"{dict(Counter(e.label for e in fa)) if fa else ''}")
    summ = {}
    for lab, rr in rows.items():
        ok_levels = []
        for r in rr:
            if r["found"] >= 0.5 * r["ref"] and r["fa"] == 0:
                ok_levels.append(r)
            else:
                break
        summ[lab] = dict(n_ref=len(by[lab]), rows=rr,
                         min_snr_meas=min((r["snr_meas"] for r in ok_levels), default=None),
                         min_snr_therm=min((r["snr_therm"] for r in ok_levels if r["snr_therm"] is not None), default=None),
                         max_added_db=max((r["added_db"] for r in ok_levels), default=None),
                         fa_total=sum(r["fa"] for r in rr))
    return dict(name=name, labels=summ, false_alarms=fa_any, n_ref=len(ref))


def negative(name, x, fs, fc, log=print):
    t = time.time()
    em = D.detect(x, fs, fc, {})
    log(f"\n=== {name}: {len(em)} DJI detections {dict(Counter(e.label for e in em))} ({time.time() - t:.1f} s)")
    return len(em)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", default=",".join(map(str, LEVELS)))
    ap.add_argument("--only")
    ap.add_argument("--json")
    a = ap.parse_args()
    levels = [int(v) for v in a.levels.split(",")]
    caps = {k: v for k, v in CAPTURES.items() if not a.only or k in a.only.split(",")}
    res = [evaluate(k, v, levels) for k, v in caps.items()]
    # negatives
    neg = {}
    rng = np.random.default_rng(1)
    w = ((rng.standard_normal(2_000_000) + 1j * rng.standard_normal(2_000_000)) / np.sqrt(2)).astype(np.complex64)
    neg["pure_noise_5780"] = negative("pure noise (fc 5780)", w, 61.44e6, 5780e6)
    neg["pure_noise_2440"] = negative("pure noise (fc 2440)", w, 61.44e6, 2440e6)
    for k, v in NEGATIVE.items():
        x, fs, fc, _ = C.load_capture(os.path.join(HERE, v))
        neg[k] = negative(k, x, fs, fc)
        for d in (10, 20):
            neg[f"{k}+{d}dB"] = negative(f"{k} +{d} dB noise", C.add_awgn(x, d, np.random.default_rng(d)), fs, fc)

    print("\n" + "=" * 118)
    print(f"{'capture':15s} {'label':16s} {'n_ref':>5s}  " + " ".join(f"+{d:<2d}dB" for d in levels)
          + f"  {'minSNR meas/therm':>17s} {'FA':>3s}  result")
    print("-" * 118)
    allpass = True
    for r in res:
        if not r["labels"]:
            print(f"{r['name']:15s} {'(no DJI emission found)':16s}")
        for lab, s in r["labels"].items():
            cells = " ".join(f"{x_['found']:>2d}/{x_['ref']:<2d}" + ("*" if x_["fa"] else " ") for x_ in s["rows"])
            ok = s["fa_total"] == 0 and s["rows"][0]["found"] == s["n_ref"]
            allpass &= ok
            ms = f"{s['min_snr_meas']}/{s['min_snr_therm']}"
            print(f"{r['name']:15s} {lab:16s} {s['n_ref']:>5d}  {cells}  {ms:>17s} {s['fa_total']:>3d}  "
                  f"{'PASS' if ok else 'FAIL'}")
    print("-" * 118)
    print("  cells: found/reference at each added-noise level, * = false alarm of that family at that level")
    print("  per-level measured SNR of the reference burst (meas | therm):")
    for r in res:
        for lab, s in r["labels"].items():
            print(f"    {r['name']:15s} {lab:16s} " + "  ".join(
                f"{x_['snr_meas']:6.1f}|{x_['snr_therm'] if x_['snr_therm'] is not None else float('nan'):5.1f}"
                + (f" dec{x_['decoded']}" if lab == "droneid" else "") for x_ in s["rows"]))
    print("-" * 118)
    for k, n in neg.items():
        ok = n == 0
        allpass &= ok
        print(f"  false-alarm test {k:28s} DJI detections {n:3d}   {'PASS' if ok else 'FAIL'}")
    print("=" * 118)
    print("OVERALL", "PASS" if allpass else "FAIL")
    if a.json:
        C.save_emissions(a.json, [], dict(results=res, negatives=neg))


if __name__ == "__main__":
    main()
