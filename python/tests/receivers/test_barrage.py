"""Barrage-noise / noise-FM jammers (research, signal gaps S12) and the services they must not be confused with."""
import numpy as np
import pytest

from spectrum_intel.core.dsp import Channelizer
from spectrum_intel.receivers.jammer import barrage as B
from spectrum_intel.receivers.jammer import synth as S

FC = 1842.5e6


def run(x, f0, f1, fc=FC):
    return B.analyse(Channelizer(x, S.FS, fc), fc + f0, fc + f1)


@pytest.mark.parametrize("f0,f1,snr", [(-10e6, 10e6, 15), (-10e6, 10e6, 8), (-10e6, 20e6, 10)])
def test_barrage(f0, f1, snr):
    r = run(S.barrage(f0, f1, snr, seed=3), f0, f1)
    assert r["class"] == "barrage" and 1.9 < r["envelope_k"] < 2.1 and r["cp_corr"] < 0.025


def test_narrow_plateau_needs_an_allocation():
    r = run(S.barrage(-2.5e6, 2.5e6, 15), -2.5e6, 2.5e6)        # a WCDMA carrier is just as Gaussian
    assert r["plateau_ok"] and r["class"] is None
    assert B.covers_band(1805e6, 1880e6) == ["b3_dl"]


def test_noise_fm():
    r = run(S.noise_fm(5e6), -12e6, 12e6)
    assert r["class"] == "noise_fm" and r["envelope_k"] < 1.2 and 2.2 < r["if_kurtosis"] < 4.5


@pytest.mark.parametrize("name,x,f0,f1", [
    ("LTE 20 MHz 0 dB", lambda: S.ofdm(18e6, snr_db=0), -9e6, 9e6),
    ("LTE 10 MHz 5 dB", lambda: S.ofdm(9e6, snr_db=5, f_off=10e6), 5.5e6, 14.5e6),
    ("Wi-Fi-like 0 dB", lambda: S.ofdm(16.6e6, tu=3.2e-6, cp=0.8e-6, snr_db=0), -8.3e6, 8.3e6),
    ("NR 30 kHz 5 dB", lambda: S.ofdm(38e6, tu=33.33e-6, cp=2.34e-6, snr_db=5), -19e6, 19e6),
    ("TDD bursts", lambda: S.ofdm_tdd(), -10e6, 10e6),
    ("swept + FM noise", lambda: S.sweep(10e3, -15e6, 15e6, noise_fm_hz=2e6), -15e6, 15e6),
    ("analog FM video", lambda: S.fm_video(), -10e6, 10e6),
    ("receiver noise", lambda: S.noise(np.random.default_rng(5)), -10e6, 10e6),
])
def test_negatives(name, x, f0, f1):
    assert run(x(), f0, f1)["class"] is None, name


def test_pipeline_label():
    from spectrum_intel.analyze import analyze_iq
    from spectrum_intel.api import contract as K
    from spectrum_intel.decode.registry import plugins
    a = analyze_iq(S.barrage(-10e6, 10e6, 15), S.FS, FC, plugins=plugins(enable=["jammer"]))
    sig = [s for s in a.signals if s.bw > 10e6]
    assert [s.label for s in sig] == ["Jammer (barrage noise)-like"] and sig[0].decoded["status"] == "partial"
    assert K.technology(sig[0].label, None)[0] == "noise_jammer"
