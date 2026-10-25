#!/usr/bin/env python3
"""Bench: synthetic NBFM speech (no sub-audible code) detection rate vs SNR, current voice_fm vs a baseline copy.
    PYTHONPATH=<baseline/python> python bench/voice_env_check.py   (baseline)  /  python bench/voice_env_check.py
"""
import numpy as np
from spectrum_intel.receivers.voice_fm import detector as D, synth
fs = 2.4e6
for snr in (8, 10, 12, 15, 20, 30):
    hit = 0
    for seed in range(12):
        rng = np.random.default_rng(seed)
        x = synth.nbfm(fs, 0.034, 150e3, audio="speech", snr_db=snr, seed=seed)
        ems = D.detect(x, fs, 446.0e6)
        hit += any(abs(e.center_mhz * 1e6 - 446.15e6) < 3e3 for e in ems)
    print(f"SNR {snr:2d} dB (12.5 kHz): detected {hit}/12")
