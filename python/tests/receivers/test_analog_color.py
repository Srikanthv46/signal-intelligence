#!/usr/bin/env python3
"""Colour decoding of analog FM video (analog_vtx FieldDemodulator.rgb).

  * synthetic colour bars (synth.py, standard YUV bars) in NTSC and PAL, clean and at 12 dB
    in-band SNR, clock errors -60 / +80 ppm, audio subcarrier at 6.5 (NTSC) / 5.5 MHz (PAL,
    1.07 MHz above the subcarrier): every bar within 24/255 per RGB channel;
  * no colour burst (monochrome camera), clean and at 12 dB -> rgb None, the colour killer;
  * the real NTSC VTX (testdata/Analog_20260804-161815): colour, and the luma image unchanged;
  * the user's VTX live at 5945 MHz (testdata/live), a 17 dB and a weak 10.7 dB block: colour;
  * color=False -> rgb None.

The demodulator is built from the synthetic VTX's true parameters, so this checks the colour
decoding alone, not the detector's localisation. The noisy case runs seeds 3 and 4: with seed 3
the old fold-sharpness line-rate refit ran to its search edge (118 ppm for 80); the per-line sync
edge fit (the H-PLL) gives 80.0. No burst: seeds 5-8 (the burst gate is tried at three offsets,
so the margin below the killer is checked on more noise draws).

Usage: pytest -q tests/receivers/test_analog_color.py   (or run it as a script)
"""
from __future__ import annotations
import os

import numpy as np

from spectrum_intel.core.emission import load_capture
from spectrum_intel.receivers.analog_vtx import detect, synth
from spectrum_intel.receivers.analog_vtx.video import BURST_MIN, FieldDemodulator
from _data import DATA

FS = 61.44e6
FC = 5800e6
# Y, U, V of 75 % colour bars (white, yellow, cyan, green, magenta, red, blue, black)
BARS = np.array([[1, 0, 0], [0.89, -0.44, 0.1], [0.7, 0.15, -0.62], [0.59, -0.29, -0.52],
                 [0.41, 0.29, 0.52], [0.3, -0.15, 0.62], [0.11, 0.44, -0.1], [0, 0, 0]])
SYNTH_CHROMA = 0.6 / 0.7           # synth.composite: chroma 0.6 V per unit U/V against 0.7 V white
CAP = os.path.join(DATA, "Analog_20260804-161815", "Analog_20260804-161815.sigmf-meta")


def _bars_image(rng, h=576, w=720, kind="mixed"):
    k = np.minimum((np.arange(w) / w * 8).astype(int), 7)
    return tuple(np.repeat(BARS[k, i][None, :], h, 0) for i in range(3))


def _expected_rgb():
    y, u, v = BARS.T
    u, v = u * SYNTH_CHROMA, v * SYNTH_CHROMA
    return np.clip(np.stack([y + 1.140 * v, y - 0.395 * u - 0.581 * v, y + 2.032 * u], 1), 0, 1) * 255


def _vtx(std, ppm, snr_db=None, color=True, seed=3):
    rng = np.random.default_rng(seed)
    make, comp = synth.make_image, synth.composite
    synth.make_image = _bars_image
    if not color:
        synth.composite = lambda *a, **k: comp(*a, **{**k, "color": False})
    try:
        x, info = synth.synth_vtx(2_000_000, FS, 3e6, std, 6e6, ppm, rng=rng, return_info=True,
                                  audio_hz=6.5e6 if std == "ntsc" else 5.5e6)
    finally:
        synth.make_image, synth.composite = make, comp
    if snr_db is not None:
        half = info["bw99_hz"] / 2
        x = synth.embed(synth.awgn(len(x), rng), x, FS, FC, snr_db, (FC + 3e6 - half) / 1e6, (FC + 3e6 + half) / 1e6)[0]
    dm = FieldDemodulator(FC + info["mean_hz"], info["line_rate_hz"], std, FC + info["sync_hz"], info["bw99_hz"])
    return dm.process(x, FS, FC)


def _bar_errors(rgb):
    cols = ((np.arange(8) + 0.5) / 8 * rgb.shape[1]).astype(int)
    got = np.array([np.median(rgb[:, c - 20:c + 20].reshape(-1, 3), 0) for c in cols])
    return np.abs(got - _expected_rgb()).max(1)


def test_color_bars():
    for std in ("ntsc", "pal"):
        for ppm, snr, seed in ((-60.0, None, 3), (80.0, 12.0, 3), (80.0, 12.0, 4)):
            fr = _vtx(std, ppm, snr, seed=seed)
            assert fr is not None and fr.rgb is not None, (std, ppm, snr, fr and fr.meta)
            assert fr.rgb.shape == fr.image.shape + (3,) and fr.rgb.dtype == np.uint8
            assert fr.meta["color"] is True
            err = _bar_errors(fr.rgb)
            assert err.max() <= 24, (std, ppm, snr, err.tolist())


def test_no_burst_is_grey():
    for std in ("pal", "ntsc"):
        for snr, seed in ((None, 5), (12.0, 5), (12.0, 6), (12.0, 7), (8.0, 8)):
            fr = _vtx(std, 30.0, snr, color=False, seed=seed)
            assert fr is not None and fr.rgb is None and fr.meta["color"] is False, (std, snr, fr.meta)
            assert fr.meta["burst_field"] < BURST_MIN / 1.5                # margin below the killer


def test_real_ntsc_capture():
    x, fs, fc, _ = load_capture(CAP)
    e = [e for e in detect(x, fs, fc) if e.family == "analog_vtx"][0]
    grey = FieldDemodulator.from_emission(e)
    grey.color = False
    g = grey.process(x, fs, fc)
    assert g is not None and g.rgb is None
    fr = FieldDemodulator.from_emission(e).process(x, fs, fc)
    assert fr.meta["color"] is True and fr.meta["burst_ratio"] > 10
    assert np.array_equal(fr.image, g.image)                      # luma untouched by the colour path
    rgb = fr.rgb.astype(float)
    sat = rgb.max(2) - rgb.min(2)
    assert np.median(sat) > 15                                     # a coloured picture, not grey
    assert (rgb[..., 0] > rgb[..., 2]).mean() > 0.6                # mostly a hand: skin is red over blue


def test_live_vtx_5945():
    """The user's NTSC VTX off the Pluto: a 17 dB block and a weak 10.7 dB one, both in colour."""
    for name in ("LIVE_ANALOG_VTX_5945MHz", "LIVE_ANALOG_VTX_5945MHz_weak"):
        x, fs, fc, _ = load_capture(os.path.join(DATA, "live", name, name + ".sigmf-meta"))
        em = [e for e in detect(x, fs, fc) if e.family == "analog_vtx"]
        assert len(em) == 1 and em[0].label == "ntsc", name
        fr = FieldDemodulator.from_emission(em[0]).process(x, fs, fc)
        assert fr is not None and fr.rgb is not None, (name, fr and fr.meta)
        assert fr.meta["burst_field"] > 2 * BURST_MIN, (name, fr.meta)


if __name__ == "__main__":
    test_color_bars(); test_no_burst_is_grey(); test_real_ntsc_capture(); test_live_vtx_5945()
    print("ok")
