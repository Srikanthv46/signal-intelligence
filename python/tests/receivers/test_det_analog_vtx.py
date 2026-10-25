#!/usr/bin/env python3
"""
test_det_analog_vtx.py -- evaluation of det_analog_vtx on synthetic analog VTXs
(analog_vtx_synth) embedded in AWGN and in real captures, plus false-alarm runs on
the real captures with nothing added.

  python test_det_analog_vtx.py [--trials 10] [--quick] [--out DIR]

Prints a PASS/FAIL table and the Pd-vs-SNR curve; writes results JSON, a Pd plot and
a demodulated-picture PNG to --out (default: ./out/avtx).

What counts as a detection: an Emission of family analog_vtx whose sync-tip frequency
or centre lies inside the true VTX's 99 % band (+-1 MHz); how well it is localised
(centre, sync tip, deviation, line rate) is scored separately. Correct classification additionally needs the PAL/NTSC
label to match. Anything else reported is a false alarm.
"""
from __future__ import annotations
import argparse, glob, json, os, sys, time
import os
import numpy as np

from spectrum_intel.core import common as C
from spectrum_intel.receivers.analog_vtx import synth as A
from spectrum_intel.receivers.analog_vtx import detector as D
from _data import DATA as HERE, OUT
CAP_WIFI = os.path.join(HERE, "D_20260916-171343/D_20260916-171343.sigmf-meta")
CAP_DJI4 = glob.glob(os.path.join(HERE, "DJI_Mini_Pro_4_*/*.sigmf-meta"))[0]
# captures with no analog VTX (false-alarm set); Analog_* holds a real one (positive check)
REAL = [CAP_WIFI] + sorted(p for p in glob.glob(os.path.join(HERE, "*/*.sigmf-meta"))
                           if "/Analog_" not in p)
CAP_ANALOG = sorted(glob.glob(os.path.join(HERE, "Analog_*/*.sigmf-meta")))
FS = 61.44e6
N = 2_000_000


def load(path):
    x, fs, fc, _ = C.load_capture(path)
    return x[:N], fs, fc


def make_trial(rng, bgname, bgs, snr_db):
    """One synthetic VTX with random content/standard/deviation/clock/polarity/audio,
    embedded at snr_db (in-band, measure_snr definition) into background bgname."""
    std = rng.choice(["pal", "ntsc"])
    dev = rng.uniform(4e6, 10e6)
    ppm = rng.uniform(-80, 80)
    pol = int(rng.choice([-1, 1]))
    audio = [None, 6.0e6, 6.5e6][rng.integers(3)]
    if bgname == "awgn":
        fc = 5800e6
        bg = A.awgn(N, rng)
        f_off = rng.uniform(-15e6, 15e6)                     # anywhere, not on a channel
        ch_nom = None
    else:
        bg, fs, fc = bgs[bgname]
        ch_nom = 5769e6                                       # RaceBand R4, inside both captures
        f_off = ch_nom - fc + rng.uniform(-0.5e6, 0.5e6)      # VTX carrier offset
    x, info = A.synth_vtx(N, FS, f_off, std, dev, ppm, polarity=pol, audio_hz=audio,
                          preemph=bool(rng.random() < 0.3), rng=rng, return_info=True,
                          image_kind=rng.choice(["mixed", "fpv", "noise"], p=[0.6, 0.3, 0.1]))
    lo, hi = A.band_edges(x, FS)
    y, m, g = A.embed(bg, x, FS, fc, snr_db, (fc + lo) / 1e6, (fc + hi) / 1e6)
    info.update(fc=fc, f_lo=lo, f_hi=hi, measured_snr=m, ch_nom=ch_nom)
    return y, fc, info


def match(ems, fc, info):
    """Split emissions into the one matching the truth (or None) and false alarms."""
    sync_t = (fc + info["sync_hz"]) / 1e6
    cen_t = (fc + info["mean_hz"]) / 1e6
    hit, fa = None, []
    for e in ems:
        # sync tip or centre inside the true VTX band (localisation quality is scored
        # separately); anything else is a false alarm
        lo_t, hi_t = (fc + info["f_lo"]) / 1e6 - 1, (fc + info["f_hi"]) / 1e6 + 1
        ok = lo_t <= e.features["sync_tip_mhz"] <= hi_t or lo_t <= e.center_mhz <= hi_t
        if ok and hit is None:
            hit = e
        elif not ok:
            fa.append(e)
    return hit, fa


def demod_picture(y, fs, fc, em, path, ncols=720):
    """Discriminate the detected VTX and stack its lines into a picture (proof that the
    synth and the detector's timing are right)."""
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    F = em.features
    ch = C.Channelizer(y.astype(np.complex128), fs)
    hb = max(em.f_hi_mhz - em.center_mhz, em.center_mhz - em.f_lo_mhz) * 1e6 + 1e6
    d, ofs = D.fm_discriminate(ch, em.center_mhz * 1e6 - fc, hb, fs)
    T = 1 / F["line_rate_hz"]; t0 = F["sync_phase_us"] * 1e-6
    act0, act1 = (10.5e-6, T - 1.65e-6) if F["standard"] == "pal" else (9.4e-6, T - 1.5e-6)
    nl = int((len(d) / ofs - t0) / T) - 1
    rows = []
    for k in range(nl):
        ts = t0 + k * T + np.linspace(act0, act1, ncols)
        rows.append(np.interp(ts * ofs, np.arange(len(d)), d))
    img = np.array(rows)
    sign = 1 if F["polarity"] == "sync_low" else -1
    img = sign * img
    lo, hi = np.percentile(img, [1, 99])
    plt.figure(figsize=(9, 7))
    plt.imshow(np.clip((img - lo) / (hi - lo), 0, 1), cmap="gray", aspect="auto")
    plt.title(f"demodulated {F['standard'].upper()} lines ({nl} lines, both fields in time order)")
    plt.savefig(path, dpi=80); plt.close()
    return nl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=10)
    ap.add_argument("--quick", action="store_true", help="coarse SNR grid, 4 trials")
    ap.add_argument("--out", default=os.path.join(OUT, "avtx"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    trials = 4 if a.quick else a.trials
    snrs = np.arange(-20, 20.1, 5.0 if a.quick else 2.5)
    rows = []          # (test, PASS/FAIL, detail)
    t_start = time.time()

    bgs = {"wifi": load(CAP_WIFI), "dji4": load(CAP_DJI4)}
    bg_cycle = ["awgn", "wifi", "dji4"]

    # ---------------------------------------------------------------- false alarms
    fa_res = {}
    rng = np.random.default_rng(12345)
    n_noise_fa = 0
    for i in range(10):
        ems = D.detect(A.awgn(N, rng), FS, 5800e6, {})
        n_noise_fa += len(ems)
    fa_res["awgn_x10"] = n_noise_fa
    rows.append(("false alarms: pure noise (10 runs)", n_noise_fa == 0, f"{n_noise_fa} detections"))
    runtimes = []
    for p in REAL:
        x, fs, fc = load(p)
        q, P, t, f = D.tf_power(x, fs)
        score = D.scan(q, t, f)["score"]
        t0 = time.time(); ems = D.detect(x, fs, fc, {}); runtimes.append(time.time() - t0)
        name = os.path.basename(os.path.dirname(p)) if os.path.basename(os.path.dirname(p)).startswith(("DJI_", "Herelink_", "Analog_")) else os.path.basename(p).split(".")[0]
        fa_res[name] = dict(n=len(ems), max_score=round(score, 2), fc_mhz=fc / 1e6, runtime_s=round(runtimes[-1], 1))
        rows.append((f"false alarms: {name} (fc {fc/1e6:g})", len(ems) == 0,
                     f"{len(ems)} detections, max score {score:.1f} (thr {D.THRESHOLD})"))

    for p in CAP_ANALOG:
        x, fs, fc = load(p)
        ems = D.detect(x, fs, fc, {})
        name = os.path.basename(os.path.dirname(p))
        fa_res[name] = dict(n=len(ems), fc_mhz=fc / 1e6,
                            found=[(e.label, round(e.center_mhz, 3), e.features.get("channel")) for e in ems])
        rows.append((f"real analog VTX: {name} (fc {fc/1e6:g})", len(ems) == 1,
                     f"{[(e.label, round(e.center_mhz, 2), e.features.get('channel'), round(e.snr_db or 0, 1)) for e in ems]}"))

    # ---------------------------------------------------------------- Pd vs SNR
    curve = []
    per = []
    rng = np.random.default_rng(2026)
    demo = None
    for snr in snrs:
        hits = cls = fas = 0
        for k in range(trials):
            bgn = bg_cycle[k % 3]
            y, fc, info = make_trial(rng, bgn, bgs, snr)
            t0 = time.time(); ems = D.detect(y, FS, fc, {}); dt = time.time() - t0
            runtimes.append(dt)
            hit, fa = match(ems, fc, info)
            fas += len(fa)
            rec = dict(snr=float(snr), bg=bgn, std=info["std"], dev_mhz=info["dev_hz"] / 1e6, dev_obs_mhz=info["dev_obs_hz"] / 1e6,
                       ppm=info["clock_ppm"], measured_snr=info["measured_snr"], detected=hit is not None,
                       n_fa=len(fa), runtime_s=dt,
                       truth_sync_mhz=(fc + info["sync_hz"]) / 1e6, truth_centre_mhz=(fc + info["mean_hz"]) / 1e6,
                       reports=[(e.label, e.center_mhz, e.features["sync_tip_mhz"], e.features["score"],
                                 e.features["line_rate_ppm"]) for e in ems])
            if hit is not None:
                hits += 1
                F = hit.features
                cls += hit.label == info["std"]
                rec.update(label=hit.label, ppm_est=F["line_rate_ppm"], ppm_sigma=F["line_rate_ppm_sigma"],
                           dev_est=F["deviation_mhz"], dev_src=F.get("deviation_source"),
                           centre_err_mhz=hit.center_mhz - (fc + info["mean_hz"]) / 1e6,
                           sync_err_mhz=F["sync_tip_mhz"] - (fc + info["sync_hz"]) / 1e6,
                           snr_est=hit.snr_db, channel=F.get("channel"), spectrum=F.get("spectrum_source"),
                           audio_true=info["audio_hz"], audio_est=F.get("audio_subcarrier_mhz"),
                           lpf=F.get("lines_per_field"))
                if demo is None and snr >= 15 and bgn == "awgn":
                    demo = (y, FS, fc, hit)
            per.append(rec)
        curve.append(dict(snr=float(snr), pd=hits / trials, pcls=cls / max(hits, 1), fa=fas))
        print(f"  SNR {snr:6.1f} dB  Pd {hits/trials:4.2f}  correct std {cls}/{hits}  extra reports {fas}", flush=True)

    # ---------------------------------------------------------------- two VTXs
    rng = np.random.default_rng(77)
    bg, fs, fc = bgs["wifi"]
    x1, i1 = A.synth_vtx(N, FS, 5752e6 - fc, "pal", 7e6, -20, rng=rng, return_info=True)
    x2, i2 = A.synth_vtx(N, FS, 5790e6 - fc, "ntsc", 7e6, 45, rng=rng, return_info=True)
    two = {}
    for snr in (10.0, 0.0):
        y = bg.copy()
        for xi, ii in ((x1, i1), (x2, i2)):
            lo, hi = A.band_edges(xi, FS)
            y, _, _ = A.embed(y, xi, FS, fc, snr, (fc + lo) / 1e6, (fc + hi) / 1e6)
        ems = D.detect(y, FS, fc, {})
        got = sorted((e.center_mhz, e.label, e.features.get("channel"), e.features["line_rate_ppm"]) for e in ems)
        ok = (len(ems) == 2 and got[0][1] == "pal" and got[1][1] == "ntsc"
              and abs(got[0][0] - 5752) < 3 and abs(got[1][0] - 5790) < 3)
        two[snr] = got
        rows.append((f"two VTXs B2 PAL + B4 NTSC on Wi-Fi capture @ {snr:g} dB", ok, str(got)))

    # ---------------------------------------------------------------- summaries
    det = [r for r in per if r["detected"]]
    pd_at = {c["snr"]: c["pd"] for c in curve}
    tot_fa = sum(c["fa"] for c in curve)
    ok_snrs = [s for s in sorted(pd_at) if all(pd_at[t] >= 0.9 for t in pd_at if t >= s)]
    snr90 = ok_snrs[0] if ok_snrs else None
    rows.append(("Pd >= 0.9 at in-band SNR <= 0 dB", snr90 is not None and snr90 <= 0, f"Pd>=0.9 from {snr90} dB up"))
    rows.append(("Pd >= 0.9 at in-band SNR <= -5 dB (stretch)", snr90 is not None and snr90 <= -5, f"{snr90} dB"))
    rows.append(("no extra (false) reports in any Pd trial", tot_fa == 0, f"{tot_fa}"))
    ncls = sum(r["label"] == r["std"] for r in det)
    rows.append(("PAL/NTSC classification of detections", ncls == len(det), f"{ncls}/{len(det)} correct"))
    hi = [r for r in det if r["snr"] >= 0]
    if hi:
        e_ppm = np.array([r["ppm_est"] - r["ppm"] for r in hi])
        sg = np.array([r["ppm_sigma"] for r in hi])
        within = np.mean(np.abs(e_ppm) <= np.maximum(3 * sg, 1.0))
        rows.append(("line rate err within max(3 x reported sigma, 1 ppm), >=95 %, SNR >= 0", within >= 0.95,
                     f"{within*100:.0f} %, rms {np.sqrt(np.mean(e_ppm**2)):.2f} ppm, max {np.abs(e_ppm).max():.2f}"))
        e_c = np.array([r["centre_err_mhz"] for r in hi])
        rows.append(("centre |err| < 1 MHz (SNR >= 0 dB)", np.all(np.abs(e_c) < 1.0),
                     f"rms {np.sqrt(np.mean(e_c**2)):.2f} MHz, max {np.abs(e_c).max():.2f}"))
        e_d = np.array([r["dev_est"] - r["dev_obs_mhz"] for r in hi if r["dev_est"] is not None])
        rows.append(("deviation (vs picture's actual swing) |err| < 1.5 MHz, SNR >= 0", len(e_d) and np.all(np.abs(e_d) < 1.5),
                     f"bias {e_d.mean():+.2f}, rms {np.sqrt(np.mean(e_d**2)):.2f} MHz, max {np.abs(e_d).max():.2f}"))
        ch_ok = [r for r in hi if r["bg"] != "awgn"]
        nch = sum(r["channel"] == "R4" or (r["channel"] or "").startswith("B3") for r in ch_ok)
        rows.append(("channel R4 (or co-located B3) on real backgrounds, SNR >= 0", nch == len(ch_ok), f"{nch}/{len(ch_ok)}"))
    for bgset, lab in ((("awgn", "dji4"), "AWGN / DJI bg"), (("wifi",), "busy Wi-Fi bg")):
        aud = [r for r in det if r["snr"] >= 10 and r["bg"] in bgset]
        if aud:
            nok = sum((r["audio_true"] is None and r["audio_est"] is None) or
                      (r["audio_true"] is not None and r["audio_est"] is not None and abs(r["audio_est"] * 1e6 - r["audio_true"]) < 1e5)
                      for r in aud)
            rows.append((f"audio subcarrier present/absent + freq, SNR >= 10 dB, {lab}", nok >= 0.9 * len(aud), f"{nok}/{len(aud)}"))
    lp = [r for r in det if r["snr"] >= 0 and r.get("lpf") is not None]
    if lp:
        nlp = sum(abs(r["lpf"] - (312.5 if r["std"] == "pal" else 262.5)) < 0.01 for r in lp)
        rows.append(("lines per field (when 2 field syncs seen), SNR >= 0", nlp >= 0.95 * len(lp), f"{nlp}/{len(lp)}"))
    rows.append(("runtime per capture < 30 s", max(runtimes) < 30, f"max {max(runtimes):.1f} s, mean {np.mean(runtimes):.1f} s"))
    if demo is not None:
        pth = os.path.join(a.out, "demod_picture.png")
        nl = demod_picture(*demo, pth)
        rows.append(("demodulated picture written", True, f"{pth} ({nl} lines)"))

    # ppm precision table
    prec = []
    for snr in snrs:
        rr = [r for r in det if r["snr"] == snr]
        if rr:
            e = np.array([r["ppm_est"] - r["ppm"] for r in rr])
            prec.append((float(snr), len(rr), float(np.sqrt(np.mean(e ** 2))), float(np.median([r["ppm_sigma"] for r in rr]))))

    print("\n=== analog VTX detector: PASS/FAIL ===")
    for name, ok, det_s in rows:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name:60s} {det_s}")
    print("\n=== Pd vs in-band SNR (trials per point: %d; backgrounds cycle AWGN / Wi-Fi capture / DJI Mini 4 Pro capture) ===" % trials)
    print("   SNR dB    Pd   bar")
    for c in curve:
        print(f"   {c['snr']:6.1f}  {c['pd']:4.2f}  {'#' * int(round(c['pd'] * 20))}")
    print("\n=== line-rate error vs SNR ===\n   SNR dB   n   rms err ppm   median reported sigma ppm")
    for s_, n_, e_, sg in prec:
        print(f"   {s_:6.1f}  {n_:2d}   {e_:9.2f}     {sg:9.2f}")
    print(f"\ntotal {time.time() - t_start:.0f} s")

    json.dump(dict(curve=curve, trials=per, false_alarms=fa_res, two_vtx={str(k): v for k, v in two.items()},
                   ppm_precision=prec, rows=[(n, bool(o), d) for n, o, d in rows], threshold=D.THRESHOLD),
              open(os.path.join(a.out, "results.json"), "w"), indent=1, default=float)
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        plt.figure(figsize=(6, 4))
        plt.plot([c["snr"] for c in curve], [c["pd"] for c in curve], "o-")
        plt.axhline(0.9, ls=":", c="gray"); plt.grid(alpha=.3)
        plt.xlabel("in-band SNR (dB)"); plt.ylabel("Pd"); plt.title("analog VTX detection, zero false alarms")
        plt.savefig(os.path.join(a.out, "pd_curve.png"), dpi=90); plt.close()
    except Exception as ex:
        print("plot failed:", ex)
    return 0 if all(o for _, o, _ in rows) else 1


if __name__ == "__main__":
    sys.exit(main())
