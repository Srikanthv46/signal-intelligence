#!/usr/bin/env python3
"""
test_det_rc_link.py -- evaluation of det_rc_link.py (PASS/FAIL table).

  1. ELRS LoRa Pd vs in-band SNR (-20..+10 dB, 2.5 dB steps, >=10 captures per
     point, every packet of a real-timing ELRS link counted) on white noise and
     embedded in real captures (Phantom 4 @2450, Herelink @2390); fingerprint
     precision (cfo_ppm, clock_ppm) from the same runs.
  2. GFSK FHSS Pd vs SNR (-5..+15 dB): ELRS-FLRC-like 650 kBd, ELRS/LR1121-GFSK-
     like 300 kBd, FrSky-D16-like 70 kBd.
  3. False alarms: pure noise and all real captures without injection
     (LoRa / ELRS labels must be zero; other labels are reported).
  4. Grouping: two simultaneous ELRS links (different rates; and same rate with
     different hop sequence / phase / crystal) -> group_hops.
  5. Runtime per capture.

In-band SNR is set against each capture's own per-bin noise floor
(rc_link_synth.embed) and spot-checked with drone_common.measure_snr.

Usage: python test_det_rc_link.py [--quick] [--procs N]
"""
from __future__ import annotations
import os, sys, time, json, argparse
from collections import Counter, defaultdict
from multiprocessing import Pool
import numpy as np

from spectrum_intel.core import common as C
from spectrum_intel.receivers.rc_link import detector as D
from spectrum_intel.receivers.rc_link import synth as S
FS = 61.44e6
CAPS = {
    "p4": ("DJI_Phantom_4_20251119-182450/D_20251119-182450.sigmf-meta"),
    "herelink": ("Herelink_20251125-110453/D_20251125-110453.sigmf-meta"),
    "d5770": ("D_20260916-171343/D_20260916-171343.sigmf-meta"),
    "mini2": ("DJI_Mini_2_20251119-165139/D_20251119-165139.sigmf-meta"),
    "mini4pro": ("DJI_Mini_Pro_4_20251119-172214/D_20251119-172214.sigmf-meta"),
}
from _data import DATA as HERE, OUT
CAPS = {k: os.path.join(HERE, v) for k, v in CAPS.items()}
_BG = {}


def background(name, seed):
    """(x, fc, floor) -- noise ('noise') or a real capture, cached per process."""
    if name == "noise":
        return S.noise_capture(rng=np.random.default_rng(1000 + seed)), 2440e6, None
    if name not in _BG:
        x, fs, fc, _ = C.load_capture(CAPS[name])
        _BG[name] = (x, fc, S.floor_density(x))
    x, fc, fl = _BG[name]
    return x, fc, fl


def _match(e, b, f_tol):
    ft = b["f_hz"] + b["info"].get("cfo_hz", 0.0)
    return abs(e.center_mhz * 1e6 - ft) < f_tol and e.start_us < b["t1_us"] and e.end_us > b["t0_us"]


# ------------------------------------------------------------------ LoRa trial
def lora_trial(args):
    mode, snr, bg, seed = args
    rng = np.random.default_rng(seed * 7919 + int((snr + 50) * 10))
    x, fc, fl = background(bg, seed)
    bw = S.ELRS_LORA_MODES[mode][2]
    b = S.elrs_link(mode, rng=rng, f_window=(fc - 26e6, fc + 26e6), dur_s=len(x) / FS - 1e-3)
    if not b:
        return []
    y, b = S.embed(x, FS, fc, b, snr, floor_bin=fl)
    meas = C.measure_snr(y, FS, fc, b[0]["t0_us"], b[0]["t1_us"],
                         (b[0]["f_hz"] - bw / 2) / 1e6, (b[0]["f_hz"] + bw / 2) / 1e6) if seed < 3 else None
    X = np.fft.fft(y)
    ems = D.detect_lora(X, FS, fc, bws=[(812.5e3, [5, 6, 7, 8])] if bw > 600e3 else None)
    out = []
    for i, bb in enumerate(b):
        m = [e for e in ems if _match(e, bb, bw / 2)]
        full = sorted([e for e in m if e.features.get("cfo_ppm") is not None], key=lambda e: abs(e.start_us - bb["t0_us"]))
        r = dict(det=bool(m), full=bool(full), elrs=any(e.family == "rc_fhss" for e in m),
                 meas_snr=meas if i == 0 else None)
        if full:
            f = full[0].features
            r["cfo_ppm_err"] = f["cfo_ppm"] - bb["info"]["crystal_ppm"]
            if f.get("clock_ppm") is not None:
                r["clock_ppm_err"] = f["clock_ppm"] - bb["info"]["crystal_ppm"]
            r["t_err_us"] = full[0].start_us - bb["t0_us"]
        out.append(r)
    out_fa = [e for e in ems if not any(_match(e, bb, bw) for bb in b)]
    return [dict(o, fa=len(out_fa) if i == 0 else 0) for i, o in enumerate(out)]


# ------------------------------------------------------------------ GFSK trial
GFSK_KINDS = {
    "elrs_flrc_like_650k": dict(bitrate=650e3, dev_hz=162.5e3, n_bytes=30, interval_us=1000, bt=1.0, label="elrs_flrc_like"),
    "elrs_gfsk_like_300k": dict(bitrate=300e3, dev_hz=100e3, n_bytes=12, interval_us=1000, bt=0.5, label="elrs_gfsk_like"),
    "frsky_d16_like_70k": dict(bitrate=70e3, dev_hz=57e3, n_bytes=30, interval_us=9000, bt=0.5, label="frsky_d16_like"),
}


def gfsk_trial(args):
    kind, snr, bg, seed = args
    rng = np.random.default_rng(seed * 104729 + int((snr + 50) * 10))
    x, fc, fl = background(bg, seed)
    kw = dict(GFSK_KINDS[kind]); label = kw.pop("label")
    ch = np.arange(fc - 25e6, fc + 25e6, 1.0e6) + 0.3e6
    b = S.gfsk_fhss_link(channels=ch, rng=rng, f_window=(fc - 26e6, fc + 26e6), dur_s=len(x) / FS - 1e-3, **kw)[:10]
    y, b = S.embed(x, FS, fc, b, snr, floor_bin=fl)
    X = np.fft.fft(y)
    ems = D.detect_fsk(y, X, FS, fc)
    out = []
    for bb in b:
        m = [e for e in ems if _match(e, bb, max(bb["bw_hz"], 300e3))]
        out.append(dict(det=bool(m), cls=any(e.label == label for e in m)))
    return out


def sweep(fn, configs, snrs, ntr, procs):
    jobs = [(c, s, bg, k) for (c, bg) in configs for s in snrs for k in range(ntr)]
    with Pool(procs) as p:
        res = p.map(fn, jobs, chunksize=1)
    table = defaultdict(list)
    for (c, s, bg, k), r in zip(jobs, res):
        table[(c, bg, s)] += r
    return table


def min_snr(rows, snrs, key="det", pd=0.9):
    """Lowest SNR from which Pd stays >= pd at every higher SNR."""
    ok = [np.mean([r[key] for r in rows[s]]) >= pd if rows[s] else False for s in snrs]
    for i in range(len(snrs)):
        if all(ok[i:]):
            return snrs[i]
    return None


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--procs", type=int, default=2)        # 8 GB RAM
    ap.add_argument("--json", default=os.path.join(OUT, "test_det_rc_link.json"))
    a = ap.parse_args()
    ntr = 3 if a.quick else 10
    results, verdicts = {}, []

    def verdict(name, ok, detail):
        verdicts.append((name, "PASS" if ok else "FAIL", detail))
        print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)

    # ---- 1. LoRa Pd vs SNR
    snrs = list(np.arange(-20, 10.01, 2.5))
    cfgs = [("2g4_500hz", "noise"), ("2g4_250hz", "noise"), ("2g4_150hz", "noise"), ("2g4_50hz", "noise"),
            ("2g4_500hz", "p4"), ("2g4_150hz", "p4"), ("2g4_500hz", "herelink")]
    t = time.time()
    tab = sweep(lora_trial, cfgs, snrs, ntr, a.procs)
    print(f"LoRa sweep {time.time() - t:.0f} s")
    print("\nELRS LoRa Pd vs in-band SNR (det = any LoRa label; full = SFD-verified with CFO/clock)")
    print("mode/background        " + " ".join(f"{s:6.1f}" for s in snrs))
    lora_req = {"2g4_500hz": -2.5, "2g4_250hz": -5.0, "2g4_150hz": -7.5, "2g4_50hz": -10.0}
    results["lora"] = {}
    for c, bg in cfgs:
        rows = {s: tab[(c, bg, s)] for s in snrs}
        for key in ("det", "full"):
            print(f"{c:13s} {bg:8s} {key:4s} " + " ".join(f"{np.mean([r[key] for r in rows[s]]) if rows[s] else float('nan'):6.2f}" for s in snrs))
        # measure_snr floors (sig - floor) at 1e-12 -> only meaningful well above 0 dB
        ms = [r["meas_snr"] for s in snrs if s >= 5 for r in rows[s] if r.get("meas_snr") is not None]
        tg = [s for s in snrs if s >= 5 for r in rows[s] if r.get("meas_snr") is not None]
        m_det, m_full = min_snr(rows, snrs, "det"), min_snr(rows, snrs, "full")
        npk = sum(len(rows[s]) for s in snrs) // len(snrs)
        fa = sum(r.get("fa", 0) for s in snrs for r in rows[s])
        results["lora"][f"{c}/{bg}"] = dict(min_snr_det=m_det, min_snr_full=m_full, packets_per_point=npk, fa=fa,
                                            pd={str(s): [float(np.mean([r[k] for r in rows[s]])) if rows[s] else None
                                                         for k in ("det", "full")] for s in snrs},
                                            measure_snr_minus_target=float(np.mean(np.array(ms) - np.array(tg))) if ms else None)
        print(f"   -> Pd>=0.9 from {m_det} dB (any), {m_full} dB (full); {npk} packets/point; "
              f"spurious LoRa dets {fa}; measure_snr - target = {results['lora'][f'{c}/{bg}']['measure_snr_minus_target']}")
        if bg == "noise":
            verdict(f"LoRa {c} Pd>=0.9 by {lora_req[c]} dB (noise)", m_det is not None and m_det <= lora_req[c],
                    f"min SNR {m_det} dB (full-parameter {m_full} dB)")
        else:
            verdict(f"LoRa {c} embedded in {bg}: Pd>=0.9 reached", m_det is not None, f"min SNR {m_det} dB (full {m_full} dB)")
    fa_lora_inj = sum(results["lora"][k]["fa"] for k in results["lora"])
    verdict("No spurious LoRa detections in injection runs", fa_lora_inj == 0, f"{fa_lora_inj}")

    # fingerprint precision vs SNR (noise background, all modes pooled per SNR)
    print("\nFingerprint precision (noise bg): std / bias of cfo_ppm and clock_ppm errors, timing error")
    results["fingerprint"] = {}
    for c in lora_req:
        for s in snrs:
            rr = [r for r in tab[(c, "noise", s)] if "cfo_ppm_err" in r]
            if len(rr) < 5:
                continue
            ce = np.array([r["cfo_ppm_err"] for r in rr]); ke = np.array([r["clock_ppm_err"] for r in rr if "clock_ppm_err" in r])
            te = np.array([r["t_err_us"] for r in rr])
            d = dict(n=len(rr), cfo_ppm_std=float(np.std(ce)), cfo_ppm_bias=float(np.mean(ce)),
                     cfo_ppm_p90=float(np.percentile(np.abs(ce), 90)),
                     clock_ppm_std=float(np.std(ke)) if len(ke) > 2 else None,
                     clock_ppm_bias=float(np.mean(ke)) if len(ke) > 2 else None, t_err_us_rms=float(np.sqrt(np.mean(te ** 2))))
            results["fingerprint"][f"{c}@{s}"] = d
            if s in (-10.0, -5.0, 0.0, 5.0, 10.0) or s == min(x for x in snrs if len(tab[(c, 'noise', x)]) and any('cfo_ppm_err' in r for r in tab[(c, 'noise', x)])):
                print(f"  {c:10s} {s:6.1f} dB n={d['n']:3d} cfo_ppm std {d['cfo_ppm_std']:.3f} p90|err| {d['cfo_ppm_p90']:.3f}  "
                      f"clock_ppm std {d['clock_ppm_std']} bias {d['clock_ppm_bias']}  t_err rms {d['t_err_us_rms']:.2f} us")
    hi = [v for k, v in results["fingerprint"].items() if float(k.split("@")[1]) >= 0]
    verdict("cfo_ppm precision < 0.5 ppm (p90) at SNR >= 0 dB", bool(hi) and max(v["cfo_ppm_p90"] for v in hi) < 0.5,
            f"worst p90 {max(v['cfo_ppm_p90'] for v in hi):.3f} ppm" if hi else "no data")

    # ---- 2. GFSK Pd vs SNR
    gsnrs = list(np.arange(-5, 15.01, 2.5))
    gcfgs = [(k, bg) for k in GFSK_KINDS for bg in ("noise", "p4")]
    t = time.time()
    gtab = sweep(gfsk_trial, gcfgs, gsnrs, ntr, a.procs)
    print(f"\nGFSK sweep {time.time() - t:.0f} s")
    print("GFSK FHSS Pd vs in-band SNR (det = any narrowband/rc label; cls = correct protocol-class label)")
    print("kind/background               " + " ".join(f"{s:6.1f}" for s in gsnrs))
    results["gfsk"] = {}
    for k, bg in gcfgs:
        rows = {s: gtab[(k, bg, s)] for s in gsnrs}
        for key in ("det", "cls"):
            print(f"{k:20s} {bg:5s} {key:3s} " + " ".join(f"{np.mean([r[key] for r in rows[s]]) if rows[s] else float('nan'):6.2f}" for s in gsnrs))
        md, mc = min_snr(rows, gsnrs, "det"), min_snr(rows, gsnrs, "cls", 0.8)
        results["gfsk"][f"{k}/{bg}"] = dict(min_snr_det=md, min_snr_cls80=mc,
                                            pd={str(s): [float(np.mean([r[q] for r in rows[s]])) if rows[s] else None
                                                         for q in ("det", "cls")] for s in gsnrs})
        print(f"   -> Pd>=0.9 from {md} dB; class label >=0.8 from {mc} dB")
        if bg == "noise":
            verdict(f"GFSK {k} Pd>=0.9 by +7.5 dB (noise)", md is not None and md <= 7.5, f"min SNR {md} dB, class {mc} dB")

    # ---- 3. False alarms (no injection)
    print("\nFalse alarms without injection (full detect())")
    results["fa"] = {}
    lora_fa = 0
    for seed in range(4 if a.quick else 8):
        x = S.noise_capture(rng=np.random.default_rng(5000 + seed))
        e = D.detect(x, FS, 2440e6)
        results["fa"][f"noise{seed}"] = dict(Counter(z.label for z in e))
        lora_fa += sum(1 for z in e if z.features.get("modulation") == "lora")
        print(f"  noise{seed}: {dict(Counter(z.label for z in e))}")
    times = {}
    for name, p in CAPS.items():
        x, fs, fc, _ = C.load_capture(p)
        t = time.time(); e = D.detect(x, fs, fc); times[name] = time.time() - t
        cnt = Counter(f"{z.family}/{z.label}" for z in e)
        results["fa"][name] = dict(cnt)
        lora_fa += sum(1 for z in e if z.features.get("modulation") == "lora" or "elrs" in z.label)
        print(f"  {name:9s} fc {fc / 1e6:.0f}: {len(e)} bursts {dict(cnt)}  ({times[name]:.1f} s)")
    verdict("Zero LoRa/ELRS labels on noise + all real captures", lora_fa == 0, f"{lora_fa} LoRa/ELRS false alarms")
    verdict("Runtime < 60 s per capture", max(times.values()) < 60, f"max {max(times.values()):.1f} s")
    results["runtime_s"] = times

    # ---- 4. Grouping
    print("\nGrouping: two simultaneous ELRS links")
    results["grouping"] = {}
    for case, (m1, m2, bg, snr) in {"diff_rate": ("2g4_500hz", "2g4_150hz", "noise", 0.0),
                                     "same_rate": ("2g4_500hz", "2g4_500hz", "p4", 5.0)}.items():
        oks = []
        for seed in range(3):
            rng = np.random.default_rng(900 + seed)
            x, fc, fl = background(bg, seed)
            b1 = S.elrs_link(m1, rng=rng, f_window=(fc - 26e6, fc + 26e6), crystal_ppm=8.0, iq_invert=False)
            b2 = S.elrs_link(m2, rng=rng, f_window=(fc - 26e6, fc + 26e6), crystal_ppm=-11.0, iq_invert=False)
            for bb in b1: bb["link"] = 1
            for bb in b2: bb["link"] = 2
            allb = b1 + b2
            y, allb = S.embed(x, FS, fc, allb, snr, floor_bin=fl)
            ems = [e for e in D.detect(y, FS, fc) if e.features.get("modulation") == "lora"]
            truth = []
            for e in ems:
                mm = [bb for bb in allb if _match(e, bb, 406e3)]
                truth.append(mm[0]["link"] if mm else 0)
            links = D.group_hops(ems)
            good = 0
            desc = []
            for L in links:
                labs = [truth[i] for i in L["members"]]
                maj = Counter(labs).most_common(1)[0]
                purity = maj[1] / len(labs)
                T_true = S.ELRS_LORA_MODES[m1 if maj[0] == 1 else m2][7]
                ok = purity >= 0.9 and abs(L["packet_interval_us"] / T_true - 1) < 1e-3
                good += ok and maj[0] in (1, 2)
                desc.append(f"link{maj[0]} n={len(labs)} purity {purity:.2f} T={L['packet_interval_us']:.1f}us "
                            f"({L['packet_rate_hz']:.1f} Hz) cfo_ppm {L['cfo_ppm']} ch {sorted(set(L['channels_mhz']))[:6]} "
                            f"guess {L['protocol_guess']} timer_ppm {L['timer_ppm']}")
            n_true = len({bb['link'] for bb in allb if any(_match(e, bb, 406e3) for e in ems)})
            ok = good >= n_true == 2
            oks.append(ok)
            print(f"  {case} seed{seed}: {len(ems)} LoRa bursts ({len(b1)}+{len(b2)} sent), {len(links)} links")
            for d in desc: print("     " + d)
        results["grouping"][case] = oks
        verdict(f"group_hops separates two ELRS links ({case})", sum(oks) >= 2, f"{sum(oks)}/3 runs correct")

    print("\n================ SUMMARY ================")
    for n, v, d in verdicts:
        print(f"{v:4s}  {n:60s} {d}")
    results["verdicts"] = verdicts
    try:
        C.save_emissions(a.json, [], extra=dict(results=json.loads(json.dumps(results, default=str))))
    except Exception as ex:   # noqa
        print("could not save json:", ex)


if __name__ == "__main__":
    main()
