"""FMCW radio altimeters (research, signal gaps S8) on synthetic sweeps."""
import numpy as np
import pytest

from spectrum_intel.receivers.radalt import detect
from spectrum_intel.receivers.radalt import synth as S

FS = 61.44e6


@pytest.mark.parametrize("fc,bw,rep,amp,period_ms", [
    (4300e6, 104e6, 50, 1.0, 20.0),          # M.2059 type A-like
    (4320e6, 150e6, 143, 0.5, 6.993),
    (4250e6, 133e6, 300, 0.3, 3.333),        # bottom turnaround inside the capture, -10 dB per sample
])
def test_triangle(fc, bw, rep, amp, period_ms):
    em = detect(S.fmcw(FS, 0.034, fc, bw=bw, f_rep=rep, amp=amp), FS, fc)
    assert [e.label for e in em] == ["radio_altimeter"]
    f = em[0].features
    assert f["shape"] == "triangle" and em[0].confidence == 0.8
    assert abs(f["sweep_bw_mhz"] - bw / 1e6) < 0.5
    assert abs(f["period_ms"] - period_ms) < 0.01
    assert abs(f["rate_mhz_per_ms"] - 2 * bw * rep / 1e9) < 0.2


def test_sawtooth_and_single_crossing():
    em = detect(S.fmcw(FS, 0.034, 4300e6, bw=176.8e6, f_rep=1000, shape="sawtooth", amp=0.4), FS, 4300e6)
    f = em[0].features
    assert f["shape"] == "sawtooth up" and f["period_ms"] == 1.0 and em[0].confidence == 0.8
    em = detect(S.fmcw(FS, 0.034, 4300e6, bw=133e6, f_rep=12), FS, 4300e6)      # too slow for a second crossing
    assert em[0].features["shape"] == "single crossing" and em[0].confidence == 0.7


def test_negatives():
    for seed in range(8):
        r = np.random.default_rng(seed)
        x = 0.7 * (r.standard_normal(2 ** 21) + 1j * r.standard_normal(2 ** 21))
        x += 5 * np.exp(2j * np.pi * 7e6 * np.arange(2 ** 21) / FS)            # a CW carrier
        assert detect(x.astype(np.complex64), FS, 4300e6) == []
    x = S.fmcw(FS, 0.034, 4300e6)
    assert detect(x, FS, 3900e6) == [] and detect(x[::4], FS / 4, 4300e6) == []    # out of band, too narrow


def test_pipeline_label():
    from spectrum_intel.analyze import analyze_iq
    from spectrum_intel.decode.registry import plugins
    x = S.fmcw(FS, 0.034, 4300e6, bw=132.8e6, f_rep=150)
    a = analyze_iq(x, FS, 4300e6, plugins=plugins(enable=["radalt"]))
    sig = [s for s in a.signals if s.features.get("decoder") == "radalt"]
    assert [s.label for s in sig] == ["Radio altimeter (FMCW)-like (triangle, 39.84 MHz/ms, 132.8 MHz sweep, "
                                      "period 6.667 ms)"]
    from spectrum_intel.api import contract as K
    assert K.technology(sig[0].label, "radalt")[0] == "radio_altimeter"
    assert K.validation("radio_altimeter", "radalt", {})["level"] == "research"
