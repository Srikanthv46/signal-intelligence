"""ILS glide slope / localizer and VOR (research, signal gaps S9 / S21) on synthetic carriers."""
import numpy as np
import pytest

from spectrum_intel.core.dsp import Channelizer
from spectrum_intel.receivers.navaids import ILS_PAIRS, channels, detect
from spectrum_intel.receivers.navaids import synth as S


def test_channel_tables():
    assert len(ILS_PAIRS) == 40 and all(329.15 <= g <= 335.0 and 108.1 <= l <= 111.95 for l, g in ILS_PAIRS)
    vor = [f for k, f, _ in channels() if k == "vor"]
    assert 108.0 in vor and 108.05 in vor and 108.1 not in vor and 108.15 not in vor and 112.1 in vor
    assert len(vor) == 160


@pytest.mark.parametrize("m90,m150,scale", [(0.4, 0.4, 1.0), (0.35, 0.45, 0.975), (0.48, 0.32, 1.025)])
def test_glide_slope_ddm(m90, m150, scale):
    fs, fc = 2e6, 332.0e6
    x = S.ils(fs, 0.0341, 332.3e6 - fc, m90, m150, scale, snr_db=30, seed=3)    # one 1/30 s window
    em = detect(Channelizer(x, fs, fc))
    assert [e.label for e in em] == ["ils_gs"]
    f = em[0].features
    assert f["channel_mhz"] == 332.3 and f["paired_mhz"] == 111.3 and f["verified_structure"]
    assert abs(f["ddm"] - (m90 - m150)) < 0.005 and abs(f["sdm"] - (m90 + m150)) < 0.01
    assert f["side"] == ("above path" if m90 > m150 else "below path")


def test_localizer_two_frequency_and_ppm():
    fs, fc, ppm = 2e6, 109.5e6, 8.0
    f_true = 109.9e6
    x = S.ils(fs, 0.05, f_true / (1 + ppm * 1e-6) - fc, 0.25, 0.15, two_freq_hz=9e3, seed=4)
    em = detect(Channelizer(x, fs, fc), ppm)
    loc = [e for e in em if e.label == "ils_loc"]
    assert len(loc) == 1
    f = loc[0].features
    assert f["channel_mhz"] == 109.9 and f["paired_mhz"] == 333.8 and abs(f["offset_hz"]) < 50
    assert f["two_frequency_hz"] == pytest.approx(9e3, abs=100) and abs(f["ddm"] - 0.1) < 0.005
    assert f["side"] == "left of course"


@pytest.mark.parametrize("radial", [0.0, 45.0, 123.0, 270.0, 359.0])
def test_vor_radial(radial):
    fs, fc = 1e6, 113.0e6
    x = S.vor(fs, 0.0341, 113.1e6 - fc, radial, snr_db=25, seed=5)
    em = detect(Channelizer(x, fs, fc))
    assert [e.label for e in em] == ["vor"]
    f = em[0].features
    err = (f["radial_deg"] - radial + 180) % 360 - 180
    assert abs(err) < 2.0 and f["verified_structure"] and 400 < f["deviation_hz"] < 560


def test_negatives_carrier_and_fm():
    from spectrum_intel.receivers.voice_fm import synth as V
    fs, fc = 2e6, 332.0e6
    t = np.arange(int(0.05 * fs)) / fs
    rng = np.random.default_rng(6)
    cw = np.exp(2j * np.pi * 300e3 * t) + 0.01 * (rng.standard_normal(len(t)) + 1j * rng.standard_normal(len(t)))
    fm = V.nbfm(fs, 0.05, -400e3, audio="speech", snr_db=30, seed=7)        # 331.6 MHz
    em = detect(Channelizer((cw + fm).astype(np.complex64), fs, fc))
    assert em == []


def test_pipeline_labels():
    from spectrum_intel.analyze import analyze_iq
    from spectrum_intel.decode.registry import plugins
    fs, fc = 2e6, 113.0e6
    x = (S.vor(fs, 0.05, 113.1e6 - fc, 200.0, seed=8) + S.ils(fs, 0.05, 112.6e6 - fc - 0.6e6, 0.2, 0.2, seed=9))
    a = analyze_iq(x.astype(np.complex64), fs, fc, plugins=plugins(enable=["navaids"]))
    labels = sorted(s.label for s in a.signals if s.features.get("decoder") == "navaids")
    assert labels == ["VOR beacon (113.10 MHz)"]
