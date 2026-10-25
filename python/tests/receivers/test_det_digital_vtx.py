#!/usr/bin/env python3
"""
test_det_digital_vtx.py -- evaluation of det_digital_vtx.detect.

  1. Herelink SNR sweep: C.add_awgn(x, d) for d = 0..50 dB (5 dB steps, the
     brief) plus 55..70 dB to locate the sensitivity limit. In-band SNR is
     measured with C.measure_snr over t 3300-5000 us, 2392-2410 MHz (a video
     burst). Below ~0 dB that single-box estimate is noise-dominated, so the
     nominal value (45.1 dB - d, i.e. the clean measurement minus the added
     noise) is printed beside it.
  2. False alarms: 2M samples of complex white noise, three DJI captures and the
     Wi-Fi-heavy D_20260916 capture -> zero "herelink" emissions; generic OFDM
     hits must carry a sensible numerology label.
  3. Blind CP-OFDM links (blind.py, no protocol receiver): Herelink, DJI Mini 2 (15 kHz),
     Mini 4 Pro (O4, 30 kHz) and the Phantom 4's Lightbridge-like 10.94 kHz link, with their
     numerology, carrier and band; none in noise or in the Wi-Fi + RC-uplink capture.

Usage: python test_det_digital_vtx.py [--quick | --nulls]  (--quick: d in 0,30,50,60,65; --nulls: FA + numerology only)
"""
from __future__ import annotations
import sys, time
from collections import Counter
import os
import numpy as np

from spectrum_intel.core import common as C
from spectrum_intel.receivers.digital_vtx import detector as D
from _data import DATA
HERE = os.path.join(DATA, "Herelink_20251125-110453/D_20251125-110453.sigmf-meta")
NULLS = {
    "DJI_Mini_2 (5780)": os.path.join(DATA, "DJI_Mini_2_20251119-165139/D_20251119-165139.sigmf-meta"),
    "DJI_Mini_4_Pro (5760)": os.path.join(DATA, "DJI_Mini_Pro_4_20251119-172214/D_20251119-172214.sigmf-meta"),
    "DJI_Phantom_4 (2450)": os.path.join(DATA, "DJI_Phantom_4_20251119-182450/D_20251119-182450.sigmf-meta"),
    "D_20260916 Wi-Fi (5770)": os.path.join(DATA, "D_20260916-171343/D_20260916-171343.sigmf-meta"),
}
TRUE_CARRIER_MHZ = 2400.6302      # measured on the clean capture (PSS + CP phase)
TRUE_PCI = 411
N_CONTROL, N_VIDEO = 7, 6         # complete PSS subframes / 2 ms video bursts in the 32.5 ms capture

rows = []


def check(name, ok, detail):
    rows.append((name, "PASS" if ok else "FAIL", detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)


def herelink_sweep(steps):
    x0, fs, fc, _ = C.load_capture(HERE)
    snr0 = C.measure_snr(x0, fs, fc, 3300, 5000, 2392, 2410)
    print(f"\nHerelink SNR sweep (clean in-band SNR {snr0:.1f} dB)")
    print(f"{'d':>3} {'meas':>7} {'nominal':>8} {'ctrl':>5} {'video(meas/pred)':>17} {'ul':>3} {'pci':>5} "
          f"{'carrier err Hz':>14} {'foldZ':>6} {'FA':>3} {'t s':>5}")
    table = []
    for d in steps:
        x = C.add_awgn(x0, d, np.random.default_rng(100 + d)) if d else x0
        snr = C.measure_snr(x, fs, fc, 3300, 5000, 2392, 2410)
        ctx = {}
        t = time.time()
        em = D.detect(x, fs, fc, ctx)
        dt = time.time() - t
        info = ctx["digital_vtx_info"]
        hl = [e for e in em if e.family == "herelink"]
        inband = [e for e in hl if abs(e.center_mhz - TRUE_CARRIER_MHZ) < 1.0]
        fa = len(hl) - len(inband)
        ctrl = [e for e in inband if e.label == "control"]
        vid = [e for e in inband if e.label == "video"]
        vm = sum(1 for e in vid if not e.features.get("predicted"))
        vp = len(vid) - vm
        ul = sum(1 for e in inband if e.label == "uplink")
        pci = ctrl[0].features.get("pci") if ctrl else None
        cerr = (ctrl[0].features["carrier_mhz"] - TRUE_CARRIER_MHZ) * 1e6 if ctrl else None
        fz = info["pss_search"]["best_z"][25]
        nominal = snr0 - d
        print(f"{d:3d} {snr:7.1f} {nominal:8.1f} {len(ctrl):5d} {vm:8d}/{vp:<8d} {ul:3d} {str(pci):>5} "
              f"{'' if cerr is None else f'{cerr:14.0f}':>14} {fz:6.1f} {fa:3d} {dt:5.1f}", flush=True)
        reliable = (len(ctrl) >= N_CONTROL - 1 and len(vid) >= N_VIDEO - 1 and fa == 0
                    and cerr is not None and abs(cerr) < 3e3)
        table.append(dict(d=d, snr=snr, nominal=nominal, ctrl=len(ctrl), vm=vm, vp=vp, ul=ul, pci=pci,
                          cerr=cerr, fz=fz, fa=fa, reliable=reliable, t=dt, em=em))
    return table, snr0


def fingerprint_stability(em):
    """Per-burst spread of the association features on the clean capture."""
    hl = [e for e in em if e.family == "herelink"]
    vid = [e for e in hl if e.label == "video" and not e.features.get("predicted")]
    ctrl = [e for e in hl if e.label == "control"]
    out = {}
    cb = [e.features.get("cfo_burst_hz") for e in vid if e.features.get("cfo_burst_hz") is not None]
    if cb:
        out["cfo_burst_hz"] = (float(np.mean(cb)), float(np.std(cb)), len(cb))
    bu = [e.end_us - e.start_us for e in vid if not e.features.get("truncated")]
    if bu:
        out["video_burst_us"] = (float(np.mean(bu)), float(np.std(bu)), len(bu))
    st = np.array([e.start_us for e in vid if not e.features.get("truncated")])
    if len(st) > 2:
        iv = np.diff(st)
        out["video_start_interval_us"] = (float(np.mean(iv)), float(np.std(iv)), len(iv))
    pt = np.array([e.features["pss_t_us"] for e in ctrl if e.features.get("pss_t_us") is not None])
    if len(pt) > 2:
        k = np.round((pt - pt[0]) / 5000.0)
        T, a = np.polyfit(k, pt, 1)
        out["pss_timing_resid_us"] = (0.0, float(np.std(pt - (a + T * k))), len(pt))
    if ctrl:
        f = ctrl[0].features
        out["link"] = dict((k, f.get(k)) for k in ("carrier_mhz", "cfo_hz", "cfo_ppm", "clock_ppm", "frame_period_ms",
                                                     "frame_phase_us", "pci", "timing_rms_us", "bw_mhz"))
    return out


def null_tests():
    print("\nFalse-alarm tests (zero 'herelink' required)")
    rng = np.random.default_rng(7)
    noise = ((rng.standard_normal(2_000_000) + 1j * rng.standard_normal(2_000_000)) / np.sqrt(2)).astype(np.complex64)
    caps = [("complex noise 2M", (noise, 61.44e6, 2400e6))] + [(k, v) for k, v in NULLS.items()]
    for name, src in caps:
        x, fs, fc = src if isinstance(src, tuple) else C.load_capture(src)[:3]
        ctx = {}
        t = time.time()
        em = D.detect(x, fs, fc, ctx)
        dt = time.time() - t
        info = ctx["digital_vtx_info"]
        hl = [e for e in em if e.family == "herelink"]
        gen = Counter((e.family, e.label, e.features.get("numerology")) for e in em
                      if e.family != "herelink" and e.label.startswith("ofdm_link"))
        hops = [e for e in em if e.label.startswith("ofdm_hop")]
        check(f"{name}: no herelink", len(hl) == 0,
              f"{len(hl)} herelink; PSS fold max Z {max(info['pss_search']['best_z'].values()):.2f} "
              f"(thr {info['pss_search']['threshold']}); rejected cands {len(info.get('rejected', []))}; {dt:.1f}s")
        # blind links: none in noise or in the Wi-Fi + RC-uplink capture (no video link there; Wi-Fi
        # numerology is below the 10 us Tu floor and a 1.1 MHz uplink hop holds too few symbols)
        if name in ("complex noise 2M", "D_20260916 Wi-Fi (5770)"):
            check(f"{name}: no blind OFDM link", not gen, f"{dict(gen)}")
        # hopping links: none in noise or on the video-only drone captures (D_20260916's hops are truth,
        # blind_hops_check)
        if name != "D_20260916 Wi-Fi (5770)":
            check(f"{name}: no OFDM hops", not hops, f"{[(e.start_us, e.center_mhz, e.label) for e in hops]}")
        check(f"{name}: runtime < 60 s", dt < 60, f"{dt:.1f}s")


# blind CP-OFDM links (generic_detect, blind.py): every drone video link of testdata, found without
# its protocol receiver: (capture, numerology, Tu us, CP/Tu, carrier MHz (capture coordinates), bw MHz)
BLIND_TRUTH = [
    ("Herelink (LTE DL PHY, 5 ms TDD)", HERE, "lte_like", 66.667, 0.0714, 2400.63, 18.0),
    ("DJI Mini 2 OcuSync 2, 20 MHz", NULLS["DJI_Mini_2 (5780)"], "lte_like", 66.667, 0.0714, 5770.56, 18.0),
    ("DJI Mini 4 Pro O4, 10 MHz", NULLS["DJI_Mini_4_Pro (5760)"], "dji_o4_like", 33.333, 0.0714, 5770.48, 9.0),
    ("DJI Phantom 4 (Lightbridge-like, label +20 MHz)", NULLS["DJI_Phantom_4 (2450)"], "lightbridge_like", 91.429,
     0.125, 2446.1, 9.6),
]


def blind_links_check():
    print("\nBlind OFDM links (no protocol receiver)")
    for name, path, num, tu, cpf, car, bw in BLIND_TRUTH:
        x, fs, fc, _ = C.load_capture(path)
        ctx = {"herelink": False}
        t = time.time()
        em = D.detect(x, fs, fc, ctx)
        links = ctx["digital_vtx_info"]["generic_links"]
        hit = [l for l in links if abs(l["carrier_mhz"] - car) < 1.0]
        ok = (len(hit) == 1 and hit[0]["numerology"] == num and abs(hit[0]["tu_us"] - tu) < 0.1
              and abs(hit[0]["cp_fraction"] - cpf) < 0.006 and abs(hit[0]["bw_mhz"] - bw) < 1.2)
        check(f"blind: {name}", ok, f"{[(l['numerology'], l['tu_us'], l['cp_fraction'], l['carrier_mhz'], l['bw_mhz'], l['n_bursts'], l['frame_period_ms']) for l in links]} {time.time() - t:.1f}s")
        if num == "dji_o4_like" and hit:
            check("blind: O4 TDD frame 2.5 ms", abs((hit[0]["frame_period_ms"] or 0) - 2.5) < 0.05,
                  f"{hit[0]['frame_period_ms']} ms")


# hopping narrowband CP-OFDM links (RC uplinks, blind_hops): (name, capture, Tu us, hop starts ms, extra hops allowed)
TMP = os.path.join(os.path.dirname(os.path.abspath(DATA)), "tmpdata")
HOP_TRUTH = [
    # DJI Mini 2 RC uplink (testdata README): 5 hops at 5774.48 / 5796.5 / 5796.5 / 5774.48 / 5770.59 MHz;
    # the last sits 0.6 MHz from the LO tone (the dji receiver notches it; not found blind)
    ("D_20260916 OcuSync 2 uplink (15 kHz)", NULLS["D_20260916 Wi-Fi (5770)"], 66.667, [5.589, 9.591, 15.589, 19.589], 0),
    # DJI O4 RC uplink of D_20260609-151042 ch0: 4.4 MHz x 0.41 ms hops, 2.5 ms cadence, re-timed with the
    # video at 16.8 ms (the same hops on ch1-3); 6.38 ms is a further hop at 5792.5 MHz
    ("D_20260609 ch0 O4 uplink (30 kHz)", os.path.join(TMP, "D_20260609-151042/D_20260609-151042_ch0.sigmf-meta"),
     33.333, [1.37, 3.87, 8.88, 17.08, 24.57], 1),
]


def blind_hops_check():
    print("\nBlind OFDM hops (RC uplinks, no protocol receiver)")
    for name, path, tu, starts, extra in HOP_TRUTH:
        if not os.path.exists(path):
            print(f"  [SKIP] {name}: {path} not present")
            continue
        x, fs, fc, _ = C.load_capture(path)
        t = time.time()
        em = D.detect(x, fs, fc, {"herelink": False})
        hops = sorted((e for e in em if e.label.startswith("ofdm_hop") and abs(e.features["tu_us"] - tu) < 0.5),
                      key=lambda e: e.start_us)
        got = [e.start_us / 1e3 for e in hops]
        found = [s for s in starts if any(abs(g - s) < 0.05 for g in got)]
        spare = [g for g in got if not any(abs(g - s) < 0.05 for s in starts)]
        check(f"hops: {name}", len(found) == len(starts) and len(spare) <= extra,
              f"{[(round(e.start_us / 1e3, 2), e.center_mhz, e.features['bw_mhz'], e.features['burst_us']) for e in hops]} "
              f"{time.time() - t:.1f}s")


def main():
    quick = "--quick" in sys.argv
    if "--nulls" in sys.argv:                                     # false-alarm part only
        blind_links_check(); blind_hops_check(); null_tests()
        for n, st, d in rows:
            print(f"{n:<48}  {st}  {d}")
        return sum(1 for r in rows if r[1] == "FAIL")
    steps = [0, 30, 50, 60, 65] if quick else list(range(0, 51, 5)) + [55, 60, 63, 66, 69, 72]
    table, snr0 = herelink_sweep(steps)
    reliable = [r for r in table if r["reliable"]]
    lowest = min(reliable, key=lambda r: r["nominal"]) if reliable else None
    all_ok_above = lowest is not None and all(r["reliable"] for r in table if r["nominal"] >= lowest["nominal"])
    check("Herelink detected at every step d <= 50 dB",
          all(r["reliable"] for r in table if r["d"] <= 50),
          f"{sum(r['reliable'] for r in table if r['d'] <= 50)}/{sum(1 for r in table if r['d'] <= 50)}")
    check("zero out-of-band herelink emissions in the sweep", all(r["fa"] == 0 for r in table),
          f"{sum(r['fa'] for r in table)}")
    if lowest:
        check("lowest reliable SNR", all_ok_above,
              f"nominal {lowest['nominal']:.1f} dB in-band (d={lowest['d']}), every step above also reliable")
    st = fingerprint_stability(table[0]["em"])
    print("\nFingerprint stability on the clean capture (mean, std, n):")
    for k, v in st.items():
        print(f"  {k}: {v}")
    check("PCI decoded on clean capture", table[0]["pci"] == TRUE_PCI, f"{table[0]['pci']}")
    blind_links_check()
    blind_hops_check()
    null_tests()
    print("\n================ SUMMARY ================")
    w = max(len(r[0]) for r in rows)
    for n, s, d in rows:
        print(f"{n:<{w}}  {s}  {d}")
    nf = sum(1 for r in rows if r[1] == "FAIL")
    print(f"\n{len(rows) - nf}/{len(rows)} PASS")
    return nf


if __name__ == "__main__":
    sys.exit(1 if main() else 0)
