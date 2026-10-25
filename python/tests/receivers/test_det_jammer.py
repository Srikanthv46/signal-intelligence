#!/usr/bin/env python3
"""
test_det_jammer.py -- the swept-carrier jammer detector on synthetic sweeps and look-alikes, and on the bench
capture of a real jammer (testdata/live/LIVE_JAMMER_SWEEP_2460MHz). Also collected by pytest.

  python test_det_jammer.py          prints one line per case with the measured sweep and the run time
"""
from __future__ import annotations

import glob
import os
import sys
import time

import numpy as np
import pytest

from spectrum_intel.core.emission import load_capture
from spectrum_intel.receivers.jammer import detect, measure
from spectrum_intel.receivers.jammer import synth as S
from _data import DATA

LIVE = os.path.join(DATA, "live", "LIVE_JAMMER_SWEEP_2460MHz", "LIVE_JAMMER_SWEEP_2460MHz.sigmf-meta")
LIVE_ALL = [LIVE, os.path.join(DATA, "live", "LIVE_JAMMER_SWEEP_2460MHz_jitter", "LIVE_JAMMER_SWEEP_2460MHz_jitter.sigmf-meta")]
# the whole signal-vault session the testdata blocks come from (10 blocks, 0.32 s apart), when present
SESSION = os.path.expanduser("~/Documents/modules/signal-vault/signal-vault-data/captures/2026-10-19/pluto/"
                             "2460MHz_JAMMING_Height0m_Distance0m_ControllerDist0m_Gain40dB_BW61.44MHz")
FC = 2.45e9

# name, generator, true sweep rate (Hz), rate tolerance (ppm), true lowest / highest frequency (Hz from the LO)
POSITIVES = [
    ("triangle 10 kHz 45 MHz", lambda: S.sweep(10e3, -25e6, 20e6), 10e3, 50, (-25e6, 20e6)),
    ("triangle, jitter 2 us, clipped, beyond the band edge",
     lambda: S.sweep(10e3, -35e6, 15e6, jitter_us=2, clip=0.5, snr_db=25), 10e3, 300, (None, 15e6)),
    ("sawtooth 3 kHz 30 MHz", lambda: S.sweep(3e3, -15e6, 15e6, shape="sawtooth"), 3e3, 50, (-15e6, 15e6)),
    ("sine 25 kHz 50 MHz", lambda: S.sweep(25e3, -25e6, 25e6, shape="sine"), 25e3, 50, (-25e6, 25e6)),
    ("triangle 1.234 kHz 20 MHz, +5 dB/MHz", lambda: S.sweep(1234.0, -10e6, 10e6, snr_db=5), 1234.0, 200,
     (-10e6, 10e6)),
    ("triangle 50 kHz 50 MHz", lambda: S.sweep(50e3, -25e6, 25e6), 50e3, 50, (-25e6, 25e6)),
    ("noise-FM triangle 7 kHz", lambda: S.sweep(7e3, -20e6, 20e6, noise_fm_hz=1e6), 7e3, 50, (-20e6, 20e6)),
    ("triangle 10 kHz, +4 dB/MHz", lambda: S.sweep(10e3, -25e6, 20e6, snr_db=4), 10e3, 50, (-25e6, 20e6)),
    ("triangle 10 kHz over Wi-Fi traffic",
     lambda: S.sweep(10e3, -25e6, 20e6, snr_db=30) + S.wifi_bursts() - S.noise(np.random.default_rng(99)), 10e3,
     50, (-25e6, 20e6)),
]
NEGATIVES = [
    ("noise", lambda: S.noise(np.random.default_rng(9))),
    ("CW tone", S.cw),
    ("pulsed radar", S.pulsed_radar),
    ("pulsed LFM radar (a sweep, but 10 % duty)", S.lfm_radar),
    ("TDD OFDM 1 ms", S.ofdm_tdd),
    ("TDD OFDM 0.1 ms", lambda: S.ofdm_tdd(period=1e-4, duty=0.5)),
    ("periodic hopper", S.hopper),
    ("Bluetooth-like random hopper", S.bt_hopper),
    ("Wi-Fi traffic", S.wifi_bursts),
    ("FM video (PAL lines)", S.fm_video),
    ("narrow sweep 4 MHz", lambda: S.sweep(10e3, -2e6, 2e6)),
]


@pytest.mark.parametrize("name,gen,rate,ppm,band", POSITIVES, ids=[p[0] for p in POSITIVES])
def test_sweep_is_found_and_measured(name, gen, rate, ppm, band):
    em = detect(gen(), S.FS, FC)
    assert len(em) == 1, name
    F = em[0].features
    assert abs(F["sweep_rate_hz"] / rate - 1) < ppm * 1e-6, (name, F["sweep_rate_hz"])
    lo, hi = band
    df = S.FS / 64
    if lo is not None:
        assert abs(F["sweep_low_mhz"] * 1e6 - (FC + lo)) < 1.5 * df, (name, F["sweep_low_mhz"])
    assert abs(F["sweep_high_mhz"] * 1e6 - (FC + hi)) < 1.5 * df, (name, F["sweep_high_mhz"])
    if "beyond" in name:                     # leaves the band: shape may be undetermined, the exit side is known
        assert F["beyond_capture"] == ["below"] and F["shape"] in ("up-down", "undetermined"), (name, F)
    else:
        assert F["beyond_capture"] == [], (name, F["beyond_capture"])
        assert F["shape"] == ("sawtooth" if "sawtooth" in name else "up-down"), (name, F["shape"])


@pytest.mark.parametrize("name,gen", NEGATIVES, ids=[n[0] for n in NEGATIVES])
def test_no_sweep_in_look_alikes(name, gen):
    res, diag = measure(gen(), S.FS, FC)
    assert res is None, (name, diag)


def test_too_weak_sweep_is_not_reported():
    """+2 dB per MHz is below the sensitivity (~ +4 dB/MHz): nothing, rather than a wrong sweep"""
    res, _ = measure(S.sweep(10e3, -25e6, 20e6, snr_db=2), S.FS, FC)
    assert res is None


@pytest.mark.parametrize("path", LIVE_ALL, ids=[os.path.basename(p).split(".")[0] for p in LIVE_ALL])
def test_bench_jammer(path):
    """bench 2026-10-19, confirmed jammer: up-down sweep every 100.00 us, top ~2470-2478 MHz (per block), down to
    below the capture's lower edge (2429.3 MHz); ADC saturated. Long-baseline reference rate: 10000.09-10000.21 Hz.
    The _jitter block's fast transits made the folded peak jump between dwell frequencies (it was refused as a
    hopper before jumps were tested for bridging)."""
    if not os.path.exists(path):
        pytest.skip(f"{path} missing")
    x, fs, fc, _ = load_capture(path)
    em = detect(x, fs, fc)
    assert len(em) == 1
    F = em[0].features
    assert abs(F["sweep_rate_hz"] - 10000.15) < 0.5, F["sweep_rate_hz"]
    assert F["beyond_capture"] == ["below"]
    assert 2468 < F["sweep_high_mhz"] < 2480, F["sweep_high_mhz"]
    assert F["sweep_low_mhz"] < 2437, F["sweep_low_mhz"]
    assert F["shape"] == "up-down" and F["flybacks"] == 0
    assert em[0].f_lo_mhz == pytest.approx((fc - fs / 2) / 1e6)      # continues below the captured band
    assert em[0].f_hi_mhz < F["sweep_high_mhz"] + 4                     # the band stops short of 2479.5 (BLE ch 39)


@pytest.mark.skipif(not os.path.isdir(SESSION), reason="signal-vault 2026-10-19 jamming session not on this machine")
def test_bench_jammer_session():
    """all 10 blocks of the signal-vault jamming session (0.32 s apart): one jammer each, the same sweep rate"""
    rates = []
    for m in sorted(glob.glob(os.path.join(SESSION, "*", "*.sigmf-meta"))):
        x, fs, fc, _ = load_capture(m)
        em = detect(x, fs, fc)
        assert len(em) == 1, m
        rates.append(em[0].features["sweep_rate_hz"])
        assert em[0].features["shape"] == "up-down", m
    assert len(rates) == 10 and max(rates) - min(rates) < 0.3, rates


if __name__ == "__main__":
    ok = True
    for name, gen, rate, ppm, band in POSITIVES:
        x = gen()
        t = time.time()
        em = detect(x, S.FS, FC)
        dt = (time.time() - t) * 1e3
        if em:
            F = em[0].features
            err = (F["sweep_rate_hz"] / rate - 1) * 1e6
            good = abs(err) < ppm
            print(f"{'PASS' if good else 'FAIL'} {name:52s} {dt:5.0f} ms  rate {F['sweep_rate_hz']:.2f} Hz "
                  f"({err:+.0f} ppm)  {F['sweep_low_mhz']:.1f}-{F['sweep_high_mhz']:.1f} MHz  {F['shape']}")
        else:
            good = False
            print(f"FAIL {name:52s} {dt:5.0f} ms  not found")
        ok &= good
    for name, gen in NEGATIVES:
        x = gen()
        t = time.time()
        res, diag = measure(x, S.FS, FC)
        dt = (time.time() - t) * 1e3
        why = "; ".join(d.get("reject", "?") for d in diag["tried"]) or f"no cyclic line ({diag['top_line_db']} dB)"
        print(f"{'PASS' if res is None else 'FAIL'} {name:52s} {dt:5.0f} ms  {why}"[:200])
        ok &= res is None
    for path in LIVE_ALL + sorted(glob.glob(os.path.join(SESSION, "*", "*.sigmf-meta"))):
        if not os.path.exists(path):
            continue
        x, fs, fc, _ = load_capture(path)
        t = time.time()
        em = detect(x, fs, fc)
        dt = (time.time() - t) * 1e3
        F = em[0].features if em else {}
        print(f"{'PASS' if em else 'FAIL'} {os.path.basename(path).split('.')[0]:52s} {dt:5.0f} ms  "
              f"{F.get('sweep_rate_hz')} Hz {F.get('sweep_low_mhz')}-{F.get('sweep_high_mhz')} MHz "
              f"beyond {F.get('beyond_capture')} {F.get('shape')}")
        ok &= bool(em)
    sys.exit(0 if ok else 1)
