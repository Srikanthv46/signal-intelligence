#!/usr/bin/env python3
"""Tests for drone_assoc.py (emitter association: same radio hopping vs. different radios).

  real data   D_20260916-171343 (DJI Mini 2 RC uplink, 5 hops on 3 channels) must be ONE
              emitter on a 10 ms frame (slots +0 / +4 ms) with one crystal ppm
  real-signal injection   the same capture's bursts copied back as a "second drone":
              independent frame phase and/or a different crystal ppm -> must split;
              an exact clone (same phase grid AND same ppm) is physically indistinguishable
  synthetic   two OcuSync-like radios with their own crystals and frame phases -> 2 emitters
  unit        cue logic (type, simultaneity veto, decoded id, energy-edge fallback)

Run:  /Users/srikanth/radioconda/bin/python test_drone_assoc.py      (~2-3 min)
"""
import os, sys, time
import numpy as np

from _data import DATA as HERE, OUT
from spectrum_intel.receivers.emitters import bursts as AB
from spectrum_intel.receivers.emitters import association as DA
from spectrum_intel.core.emission import Emission, load_capture

CAP = os.path.join(HERE, "D_20260916-171343", "D_20260916-171343.sigmf-meta")
PH4 = os.path.join(HERE, "DJI_Phantom_4_20251119-182450", "D_20251119-182450.sigmf-meta")
RESULTS = []


def check(name, cond, info=""):
    RESULTS.append((name, bool(cond)))
    print(f"[{'PASS' if cond else 'FAIL'}] {name}  {info}")


def labels_of(emitters, lab):
    return ["".join(lab[i] for i in e.members) for e in emitters]


def pure(emitters, lab):
    """every emitter holds only one ground-truth source"""
    return all(len(set(s)) == 1 for s in labels_of(emitters, lab))


_cache = {}


def real():
    if "x" not in _cache:
        x, fs, fc, _ = load_capture(CAP)
        em0 = AB.extract_bursts(x, fs, fc)
        _cache.update(x=x, fs=fs, fc=fc, em0=em0)
    return _cache["x"], _cache["fs"], _cache["fc"], _cache["em0"]


# ------------------------------------------------------------------ real capture
def test_real_single_emitter():
    x, fs, fc, em0 = real()
    em = AB.fresh(em0)
    DA.measure_fingerprints(x, fs, fc, em)
    E = DA.associate(em, fc)
    up = [e for e in E if e.label == "ocusync_uplink"]
    check("real: 5 uplink hops found", len(em) == 5, f"{len(em)} bursts")
    check("real: all hops re-timed on learned reference symbols",
          all(e.features.get("timing_source") == "reference_symbols" for e in em))
    check("real: ONE emitter", len(up) == 1 and len(up[0].members) == 5, str([e.members for e in up]))
    if up:
        e = up[0]
        fp = e.fingerprint
        check("real: frame 10 ms, slots at +0/+4 ms", abs(e.frame_period_ms - 10.0) < 1e-3
              and e.hop_pattern.get("slot_offsets_in_frame_ms") == [0.0, 4.0], f"{e.frame_period_ms} {e.hop_pattern.get('slot_offsets_in_frame_ms')}")
        check("real: one crystal (carrier ppm spread < 0.1 ppm)", fp["carrier_ppm_spread"] < 0.1,
              f"{fp['carrier_ppm']:+.3f} ppm, spread {fp['carrier_ppm_spread']:.3f}")
        check("real: lattice residual rms < 0.5 us (incl. cross-channel delay)", fp["timing_residual_rms_ns"] < 500,
              f"{fp['timing_residual_rms_ns']} ns")
        check("real: confidence > 0.99", e.confidence > 0.99, str(e.confidence))
    # same-channel revisits: lattice residual at the ns level
    t = {i: em[i].features["t_start_us"] for i in range(5)}
    c = {i: em[i].features["carrier_hz"] for i in range(5)}
    revisits = [(i, j) for i in range(5) for j in range(i + 1, 5) if abs(c[i] - c[j]) < 50e3]
    res = []
    for i, j in revisits:
        d = t[j] - t[i]
        res.append(abs(d - round(d / 2000.0) * 2000.0) * 1e3)
    check("real: same-channel revisits within 60 ns of the 2 ms grid", len(res) >= 2 and max(res) < 60,
          f"residuals {np.round(res, 1)} ns")
    dc = [abs(c[i] - c[j]) for i, j in revisits]
    check("real: same-channel revisit carrier difference < 600 Hz", max(dc) < 600, f"{np.round(dc, 0)} Hz")


def _inject(dt_us, dppm, gain_db=0.0):
    x, fs, fc, em0 = real()
    y, inj = AB.inject_copy(x, fs, fc, em0, dt_us, dppm, gain_db=gain_db)
    em = AB.fresh(em0) + AB.copy_emissions(em0, inj)
    DA.measure_fingerprints(y, fs, fc, em)
    E = DA.associate(em, fc)
    lab = "o" * len(em0) + "c" * len(inj)
    return E, lab


def test_injection_splits():
    for name, dt, dppm in [("independent phase + 1.5 ppm", 2345.678, 1.5),
                           ("independent phase, SAME ppm (timing alone)", 2345.678, 0.0),
                           ("on the slot grid, 1.5 ppm apart (oscillator alone)", 2500.0, 1.5)]:
        E, lab = _inject(dt, dppm)
        groups = labels_of(E, lab)
        orig = [g for g in groups if "o" in g]
        check(f"inject: {name} -> split", pure(E, lab) and len(orig) == 1 and orig[0] == "ooooo",
              f"emitters {groups}")


def test_injection_clone_is_ambiguous():
    """Same frame-phase grid AND same crystal: physics cannot separate them; association
    must not claim a confident 'different'. (Documented limitation.)"""
    E, lab = _inject(2500.0, 0.0)
    groups = labels_of(E, lab)
    merged = any("o" in g and "c" in g for g in groups)
    check("inject: exact clone (same grid, same ppm) is not confidently split", merged, f"emitters {groups}")


# ------------------------------------------------------------------ synthetic
def test_synthetic_two_radios():
    fc, fs = 5770e6, 61.44e6
    for snr in (12.0, 10.0):
        radios = [dict(ppm=-2.8, phase_ms=5.6, snr_db=snr), dict(ppm=1.3, phase_ms=7.123456, snr_db=snr - 2)]
        x, truth = AB.synth_scene(radios, fc=fc, rng=np.random.default_rng(3))
        em, ch, _ = AB.extract_bursts(x, fs, fc, return_ch=True)
        DA.measure_fingerprints(x, fs, fc, em, ch=ch)
        E = DA.associate(em, fc)
        lab = []
        for e in em:
            tr = min(truth, key=lambda d: abs(d["t_us"] - e.start_us) + abs(d["f_hz"] / 1e6 - e.center_mhz) * 100)
            lab.append(str(tr["radio"]))
        big = [g for g in labels_of(E, lab) if len(g) >= 2]
        check(f"synthetic 2 radios @ {snr:g} dB: pure emitters, 2 multi-burst emitters",
              pure(E, lab) and len(big) == 2, f"{labels_of(E, lab)}")
    radios = [dict(ppm=0.7, phase_ms=3.3, snr_db=10.0)]
    x, truth = AB.synth_scene(radios, fc=fc, rng=np.random.default_rng(5))
    em, ch, _ = AB.extract_bursts(x, fs, fc, return_ch=True)
    DA.measure_fingerprints(x, fs, fc, em, ch=ch)
    E = DA.associate(em, fc)
    check("synthetic 1 radio: one emitter", len(E) == 1, f"{[len(e.members) for e in E]} of {len(em)} bursts")
    if E:
        fp = E[0].fingerprint
        check("synthetic 1 radio: carrier ppm recovered (0.7 +- 0.05; absolute is CP-limited)", abs(fp["carrier_ppm"] - 0.7) < 0.05,
              f"{fp['carrier_ppm']:+.4f}")


def test_phantom4_runs():
    if not os.path.exists(PH4):
        return
    x, fs, fc, _ = load_capture(PH4)
    em = AB.extract_bursts(x, fs, fc)
    DA.measure_fingerprints(x, fs, fc, em)
    E = DA.associate(em, fc)
    allm = sorted(i for e in E for i in e.members)
    check("Phantom 4 capture: every emission in exactly one emitter", allm == list(range(len(em))),
          f"{len(em)} 15 kHz-OFDM bursts -> {len(E)} emitters")


# ------------------------------------------------------------------ unit
def _em(t, f, fam="dji_ocusync", lab="ocusync_uplink", dur=513.0, **feat):
    return Emission(fam, lab, t, t + dur, f - 0.57, f + 0.57, f, 10.0, features=dict(bw_mhz=1.14, scs_khz=15, **feat))


def test_units():
    a = _em(1000.0, 5774.5, t_start_us=1000.0, t_sigma_us=0.01, carrier_hz=5774.5e6 - 16e3, carrier_sigma_hz=20,
            timing_source="reference_symbols")
    b = _em(5000.0, 5796.5, t_start_us=5000.0002, t_sigma_us=0.01, carrier_hz=5796.5e6 - 16.06e3,
            carrier_sigma_hz=20, timing_source="reference_symbols")
    llr, ev = DA.same_emitter_score(a, b)
    check("unit: on-lattice, same ppm -> strongly same", llr > 8, f"{llr:.1f}")
    c = _em(5123.4, 5796.5, t_start_us=5123.4, t_sigma_us=0.01, carrier_hz=5796.5e6 + 6e3, carrier_sigma_hz=20,
            timing_source="reference_symbols")
    llr, ev = DA.same_emitter_score(a, c)
    check("unit: off-lattice, 3.8 ppm apart -> strongly different", llr < -8, f"{llr:.1f}")
    d = _em(1100.0, 5790.5)
    llr, ev = DA.same_emitter_score(a, d)
    check("unit: simultaneous on two channels -> veto", llr <= DA.LLR_VETO, f"{llr:.1f}")
    e = _em(5000.0, 5796.5, fam="dji_droneid", lab="droneid")
    check("unit: different type never same", DA.same_emitter_score(a, e)[0] < -20)
    s1 = _em(1000.0, 5774.5); s1.decoded = {"serial": "1581F4XXX"}
    s2 = _em(7000.0, 5774.5); s2.decoded = {"serial": "1581F4XXX"}
    check("unit: same decoded serial -> same", DA.same_emitter_score(s1, s2)[0] > 10)
    # energy-edge only (no fine features): still informative, never extreme
    p = _em(1000.0, 5774.5); q = _em(5000.0, 5796.5); r = _em(5321.0, 5796.5)
    lq, lr = DA.same_emitter_score(p, q)[0], DA.same_emitter_score(p, r)[0]
    check("unit: energy-edge timing fallback ranks on-lattice above off-lattice", lq > 0 > lr, f"{lq:.1f} {lr:.1f}")
    E = DA.associate([a, b, c, d])
    check("unit: associate groups a+b, separates c and d", sorted(len(x.members) for x in E) == [1, 1, 2],
          str([x.members for x in E]))


if __name__ == "__main__":
    t0 = time.time()
    for fn in (test_units, test_real_single_emitter, test_injection_splits, test_injection_clone_is_ambiguous,
               test_synthetic_two_radios, test_phantom4_runs):
        fn()
    n_ok = sum(ok for _, ok in RESULTS)
    print(f"\n{n_ok}/{len(RESULTS)} passed in {time.time() - t0:.0f} s")
    sys.exit(0 if n_ok == len(RESULTS) else 1)
