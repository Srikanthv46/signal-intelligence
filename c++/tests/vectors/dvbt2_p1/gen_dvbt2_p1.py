#!/usr/bin/env python3
"""Synthetic reference for DVB-T2 P1 (cofdm/dvbt2_p1; research, PLAN D5): P1 symbols at the elementary rate with noise
and carrier offset (find / decode), and cofdm detect on a 15.36 MSps capture (fc 650 MHz) with a DVB-T-like body and
one P1 at 8 dB.

    nice -n 19 .venv/bin/python -I c++/tests/vectors/dvbt2_p1/gen_dvbt2_p1.py
"""
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.signal import resample_poly

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parents[3] / "python"))
from spectrum_intel.receivers.cofdm import detector as D, dvbt2_p1 as P, synth as CS   # noqa: E402

rng = np.random.default_rng(4)
cases, chunks, at = [], [], 0
for s1, s2, snr in ((0, 0, 20), (1, 6, 3), (3, 9, 0), (2, 0, 10)):
    p = P.p1_symbol(s1, s2)
    y = np.r_[np.zeros(2500), p, np.zeros(2500)] * np.exp(2j * np.pi * 0.0071 * np.arange(len(p) + 5000))
    y = (y + (rng.standard_normal(len(y)) + 1j * rng.standard_normal(len(y))) * np.sqrt(10 ** (-snr / 10) / 2 * 384 / 1024))
    y = y.astype(np.complex64)
    f = P.find(y.astype(np.complex128), D.P1_THR)
    cases.append(dict(offset=at, n=len(y), found=[[s, m, d] for s, m, d in f]))
    chunks.append(y); at += len(y)
np.concatenate(chunks).astype(np.complex64).tofile(OUT / "p1.cf32")

fs, fc = 15.36e6, 650e6
body = CS.dvbt_baseband(70, mode="2K", gi="1/8", rng=rng, pilots=False)
yb = np.r_[body[:60000], P.p1_symbol(1, 2), body[60000:]]
r = Fraction(fs / (8e6 * 8 / 7)).limit_denominator(4000)
yb = resample_poly(yb, r.numerator, r.denominator)[:1 << 18]
n = len(yb)
x = yb * np.exp(2j * np.pi * 1.5e6 * np.arange(n) / fs)
x = x / np.sqrt(np.mean(np.abs(x) ** 2)) * np.sqrt(10 ** 0.8 * 8e6 / fs) + (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)
x = x.astype(np.complex64)
x.tofile(OUT / "capture.cf32")
em = D.detect(x, fs, fc)
cap = [dict(label=e.label, f_lo=e.f_lo_mhz, f_hi=e.f_hi_mhz, center=e.center_mhz, conf=e.confidence, features=e.features) for e in em]
print([(c["label"], c["features"].get("p1", {}).get("preamble")) for c in cap], [len(c["found"]) for c in cases])
(OUT / "cases.json").write_text(json.dumps(dict(cases=cases, capture=dict(fs=fs, fc=fc, emissions=cap)), indent=1, default=float))
