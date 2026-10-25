#!/usr/bin/env python3
"""Reference vectors for tests/unit/test_rx_emitters.cpp, written by python receivers.emitters.

    .venv/bin/python c++/tests/vectors/emitters/gen_emitters.py

iq_ci16.bin  synthetic 30.72 MSps capture (int16 I/Q, x = int16 / 2048 as complex64): two OcuSync-uplink-like radios (own ppm, own frame phase)
ref.json     math.hypot pairs; measure_fingerprints features + info on that capture; extract_bursts on it;
             associate() / same_emitter_score() on a hand-built multi-family emission list
"""
import json
import math
from dataclasses import replace
from pathlib import Path

import numpy as np

from spectrum_intel.core.emission import Emission
from spectrum_intel.receivers.emitters import association as A
from spectrum_intel.receivers.emitters.bursts import extract_bursts
from spectrum_intel.receivers.emitters.families import _assoc_models

OUT = Path(__file__).resolve().parent
R, N = 30.72e6, 2048                       # capture rate; the receiver works at 15.36 MSps
CP = (160, 144, 144, 144, 144, 144, 144)
FC = 5770e6


def em_json(e):
    return dict(family=e.family, label=e.label, start_us=e.start_us, end_us=e.end_us, f_lo_mhz=e.f_lo_mhz,
                f_hi_mhz=e.f_hi_mhz, center_mhz=e.center_mhz, snr_db=e.snr_db, confidence=e.confidence,
                decoded=e.decoded, features=e.features)


def slot(rng, n_sc=72, root=25):
    k = np.r_[-n_sc // 2:0, 1:n_sc // 2 + 1]
    n = np.arange(n_sc + 1)
    ref = np.exp(-1j * np.pi * root * n * (n + 1) / (n_sc + 1))[:n_sc]
    out = []
    for i, cp in enumerate(CP):
        X = np.zeros(N, complex)
        X[k % N] = ref if i in (0, 6) else (rng.choice([-1, 1], n_sc) + 1j * rng.choice([-1, 1], n_sc)) / np.sqrt(2)
        s = np.fft.ifft(X) * N / np.sqrt(n_sc)
        out.append(np.r_[s[-cp:], s])
    return np.concatenate(out)


def scene():
    rng = np.random.default_rng(7)
    n = int(0.0092 * R)
    x = (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)
    x += 0.3 * np.exp(2j * np.pi * 1.234e6 / R * np.arange(n))         # a CW spur (find_spurs -> notch)
    radios = [dict(ppm=3.0, phase_ms=0.4, chans=(-4e6, 1e6, 5e6), snr=18.0),
              dict(ppm=-4.0, phase_ms=2.73, chans=(-2e6, 3e6), snr=14.0)]
    ems = []
    for ri, r in enumerate(radios):
        e = r["ppm"] * 1e-6
        g = np.sqrt(10 ** (r["snr"] / 10) * 1.08e6 / R)
        for t_ms in (r["phase_ms"], r["phase_ms"] + 4.0, r["phase_ms"] + (6.0, 1.0)[ri]):
            t0 = t_ms * 1e-3 / (1 + e)
            a = int(round(t0 * R))
            b = slot(rng)
            if a + len(b) >= n:
                continue
            f = (FC + r["chans"][len(ems) % len(r["chans"])]) * (1 + e) - FC
            m = np.arange(a, a + len(b))
            x[a:a + len(b)] += g * b * np.exp(1j * (2 * np.pi * f / R * m + rng.uniform(0, 6.28)))
            cen = round((FC + f) / 1e6, 4)
            jit = rng.uniform(-3, 3)
            ems.append(Emission("dji_ocusync", "uplink", a / R * 1e6 + jit, (a + len(b)) / R * 1e6 + jit,
                                cen - 0.54, cen + 0.54, cen, None, 0.8, None,
                                dict(scs_khz=15.0, bw_mhz=1.08)))
    ems.sort(key=lambda e: (e.start_us, e.center_mhz))
    q = np.round(np.c_[x.real, x.imag] * 2048).clip(-32767, 32767).astype("<i2")
    return q, (q[:, 0] / 2048 + 1j * q[:, 1] / 2048).astype(np.complex64), ems


def assoc_set():
    rng = np.random.default_rng(3)
    out = []
    # herelink control: detector cfo_ppm + frame period, energy-edge timing with t_sigma_us
    for k in range(6):
        t = 101.3 + 5000.0 * k * (1 - 0.2e-6) + rng.normal(0, 0.05)
        c = 2410.597 + (0.5 if k % 2 else 0.0)
        out.append(Emission("herelink", "control", t, t + 1000.0, c - 9, c + 9, c, 40.0 + rng.normal(0, 2), 0.9, None,
                            dict(cfo_ppm=-1.16 + rng.normal(0, 0.02), cfo_ppm_sigma=0.05, frame_period_ms=4.999938,
                                 frame_phase_us=101.3, t_sigma_us=0.08, scs_khz=None if k == 3 else 0,
                                 bw_mhz=18.0, clock_ppm=0.5 + rng.normal(0, 1), clock_ppm_sigma=None)))
    # herelink video: same oscillator -> linked
    for k in range(4):
        t = 3101.8 + 5000.0 * k + rng.normal(0, 0.05)
        out.append(Emission("herelink", "video", t, t + 2013.0, 2406.0, 2424.0, 2415.0, 47.8, 0.9, None,
                            dict(cfo_ppm=-1.15, frame_period_ms=4.999938, frame_phase_us=3101.8, t_sigma_us=0.1,
                                 bw_mhz=18.0)))
    # two RC FHSS radios (clock_ppm ignored for rc_), sync words differ for one
    for k in range(8):
        r = k % 2
        t = 200.0 + 7000.0 * k + 13.0 * r
        c = 2402.0 + 2.0 * ((k * 5) % 20)
        out.append(Emission("rc_fhss", "elrs_like", t, t + 600.0 + 30 * r, c - 0.5, c + 0.5, c, 20.0 + 15 * r, 0.7,
                            {"sync_word": 0x12 if r else 0x34} if k < 6 else None,
                            dict(cfo_hz=1500.0 * (1 + 3 * r) + rng.normal(0, 50), cfo_ref="nominal",
                                 cfo_sigma_hz=100.0, clock_ppm=12.0 * r, clock_ppm_err=3.0)))
    # DroneID: decoded serials (two drones), zc roots
    for k in range(5):
        t = 600.0 + 640000.0 / 10 * k
        out.append(Emission("dji_droneid", "droneid", t, t + 640.0, 2399.5 + 10 * (k % 2), 2408.5 + 10 * (k % 2),
                            2404.0 + 10 * (k % 2), 15.0, 0.9, {"serial_number": "1581F" + ("AB" if k % 3 else "CD")},
                            dict(scs_khz=15.0, bw_mhz=9.0, zc_root=600 if k != 4 else 147,
                                 carrier_hz=(2404.0e6 + 10e6 * (k % 2)) * (1 + 2.5e-6) + rng.normal(0, 30),
                                 carrier_sigma_hz=40.0, timing_source="reference_symbols" if k < 3 else "cp_train",
                                 t_start_us=t + 0.4, t_sigma_us=0.02)))
    # DJI video 20 MHz with carrier and overlapping simultaneous bursts on two channels
    for k in range(4):
        t = 1000.0 + 2500.0 * k
        c = 5760.0 if k < 3 else 5800.0
        out.append(Emission("dji_ocusync", "video_20MHz", t, t + 1400.0 + 300 * k, c - 9, c + 9, c, None, 0.9, None,
                            dict(scs_khz=15.0, bw_mhz=18.0, carrier_hz=c * 1e6 * (1 - 1.1e-6) + 5.0 * k,
                                 carrier_sigma_hz=25.0, frame_period_ms=10.0, clock_ppm=-1.0, clock_ppm_sigma=0.4)))
    out.append(Emission("dji_ocusync", "video_20MHz", 3600.0, 4900.0, 5791.0, 5809.0, 5800.0, None, 0.9, None,
                        dict(scs_khz=15.0, bw_mhz=18.0, carrier_hz=5800e6 * (1 + 4e-6), carrier_sigma_hz=25.0)))
    # DJI uplink with a missing carrier (energy-edge timing) and an analog VTX singleton pair
    for k in range(3):
        t = 500.0 + 4000.0 * k + 2000 * (k == 2)
        out.append(Emission("dji_ocusync", "uplink", t, t + 520.0, 5740.46, 5741.54, 5741.0, 12.0, 0.8, None,
                            dict(scs_khz=15.0, bw_mhz=1.08) if k else dict(scs_khz=15.0)))
    out.append(Emission("analog_vtx", "pal", 0.0, 30000.0, 5795.0, 5805.0, 5800.0, 30.0, 0.9, None, dict(bw_mhz=8.0)))
    out.append(Emission("analog_vtx", "pal", 0.0, 30000.0, 5835.0, 5845.0, 5840.0, 10.0, 0.9, None, dict(bw_mhz=8.0)))
    out.append(Emission("telemetry", "sik", 50.0, 2000.0, 433.0, 433.1, 433.05, None, 0.6, None, {}))
    return out


def main():
    hy = np.random.default_rng(1).standard_normal((300, 2)) * np.repeat(10.0 ** np.arange(-3, 3), 50)[:, None]
    hy[:5, 1] = 0.0
    ref = dict(hypot=[[a, b, math.hypot(a, b)] for a, b in hy.tolist()])

    q, x, ems = scene()
    q.tofile(OUT / "iq_ci16.bin")
    ref["capture"] = dict(fs=R, fc=FC, emissions_in=[em_json(e) for e in ems])
    em = [replace(e, features=dict(e.features)) for e in ems]
    info = A.measure_fingerprints(x, R, FC, em)
    ref["capture"]["info"] = info
    ref["capture"]["emissions_out"] = [em_json(e) for e in em]
    mods = _assoc_models()
    ref["capture"]["emitters"] = [m.to_dict() for m in A.associate(em, FC, models=mods)]
    bursts = extract_bursts(x, R, FC)
    ref["capture"]["bursts"] = [em_json(e) for e in bursts]

    S = assoc_set()
    ref["assoc"] = dict(emissions=[em_json(e) for e in S])
    em, M, EV = A.associate(S, None, models=mods, return_matrix=True)
    ref["assoc"]["emitters"] = [m.to_dict() for m in em]
    ref["assoc"]["M"] = M.tolist()
    ref["assoc"]["pairs"] = [[i, j, A.same_emitter_score(S[i], S[j], mods)[1]]
                             for i in range(len(S)) for j in range(i + 1, len(S))
                             if (S[i].family == S[j].family) and (i + j) % 3 != 1 or j == i + 1]
    ref["assoc"]["emitters_default_models"] = [m.to_dict() for m in A.associate(S)]
    (OUT / "ref.json").write_text(json.dumps(ref, default=str))
    print("bursts", len(bursts), "capture emitters", len(ref["capture"]["emitters"]),
          "assoc emitters", len(ref["assoc"]["emitters"]), "pairs", len(ref["assoc"]["pairs"]))
    for e in ref["capture"]["emissions_out"]:
        f = e["features"]
        print(f"{e['center_mhz']:.4f} {f.get('timing_source')} t {f.get('t_start_us')} c {f.get('carrier_hz')}")


if __name__ == "__main__":
    main()
