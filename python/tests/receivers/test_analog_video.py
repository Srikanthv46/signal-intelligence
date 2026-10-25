#!/usr/bin/env python3
"""Line clock and low-SNR behaviour of the analog FM video field demodulator (analog_vtx FieldDemodulator).

  * line clock (H-PLL: every line's sync edge, robust straight-line fit): synthetic NTSC VTX
    (5.3 MHz deviation, an FPV-like picture, clock +2 ppm), the demodulator built 25 ppm off the
    true rate (a detector error / a drifted VTX): clean and 14 dB within 0.5 ppm, 6 dB within
    1 ppm, 3 dB within 3 ppm (the old fold-sharpness refit was 14-39 ppm off at any SNR);
  * noise: the sync-tip noise estimate grows as the SNR falls; the 3x3 median is on at 3 dB and
    off at 14 dB; at 0 dB colour is killed (color_killer "noise"), at 14 dB it is not;
  * the weak 1.2 GHz VTX (testdata/live/LIVE_ANALOG_VTX_1200MHz_weak, sync-tip noise 0.29 of the
    span): the line rate within 4 ppm of the detector's (the old refit: -31.9 for +2.1, a sheared
    picture), median filtered, grey (colour noise killed).

Usage: pytest -q tests/receivers/test_analog_video.py   (or run it as a script)
"""
from __future__ import annotations
import os

import numpy as np

from spectrum_intel.core.emission import load_capture
from spectrum_intel.receivers.analog_vtx import detect, synth
from spectrum_intel.receivers.analog_vtx.video import NOISY, NOISY_COLOR, FieldDemodulator
from _data import DATA

FS = 61.44e6
FC = 1200e6
CLOCK_PPM = 2.0
WEAK = os.path.join(DATA, "live", "LIVE_ANALOG_VTX_1200MHz_weak", "LIVE_ANALOG_VTX_1200MHz_weak.sigmf-meta")


def _vtx(snr_db, seed, rate_err_ppm=0.0):
    rng = np.random.default_rng(seed)
    x, info = synth.synth_vtx(2_000_000, FS, 0.0, "ntsc", 5.3e6, CLOCK_PPM, rng=rng, return_info=True,
                              audio_hz=6.5e6)
    if snr_db is not None:
        half = info["bw99_hz"] / 2
        x = synth.embed(synth.awgn(len(x), rng), x, FS, FC, snr_db, (FC - half) / 1e6, (FC + half) / 1e6)[0]
    dm = FieldDemodulator(FC + info["mean_hz"], info["line_rate_hz"] * (1 + rate_err_ppm * 1e-6), "ntsc",
                          FC + info["sync_hz"], info["bw99_hz"])
    return dm.process(x, FS, FC)


def test_line_clock():
    for snr, tol in ((None, 0.5), (14.0, 0.5), (6.0, 1.0), (3.0, 3.0)):
        for seed in (1, 2, 3):
            fr = _vtx(snr, seed, rate_err_ppm=25.0)
            assert fr is not None, (snr, seed)
            err = fr.meta["line_rate_ppm"] - CLOCK_PPM
            assert abs(err) <= tol, (snr, seed, fr.meta)


def test_noise_and_killers():
    for seed in (1, 2):
        hi, mid, lo = (_vtx(s, seed) for s in (14.0, 3.0, 0.0))
        assert hi.meta["noise"] < mid.meta["noise"] < lo.meta["noise"], (seed, hi.meta, mid.meta, lo.meta)
        assert hi.meta["denoise"] is False and hi.meta["noise"] < NOISY, hi.meta
        assert mid.meta["denoise"] is True, mid.meta
        assert hi.rgb is not None and hi.meta.get("color_killer") is None, hi.meta
        assert lo.meta["noise"] > NOISY_COLOR and lo.rgb is None and lo.meta["color_killer"] == "noise", lo.meta


def test_weak_1200():
    x, fs, fc, _ = load_capture(WEAK)
    em = [e for e in detect(x, fs, fc) if e.family == "analog_vtx"]
    assert len(em) == 1 and em[0].label == "ntsc" and abs(em[0].center_mhz - 1199.84) < 0.5, em
    fr = FieldDemodulator.from_emission(em[0]).process(x, fs, fc)
    assert fr is not None and fr.meta["vsync"] and fr.image.shape == (240, 720), fr and fr.meta
    assert abs(fr.meta["line_rate_ppm"] - em[0].features["line_rate_ppm"]) < 4, fr.meta
    assert fr.meta["denoise"] is True and fr.meta["noise"] > NOISY_COLOR, fr.meta
    assert fr.rgb is None and fr.meta["color_killer"] == "noise", fr.meta


if __name__ == "__main__":
    test_line_clock(); test_noise_and_killers(); test_weak_1200()
    print("ok")
