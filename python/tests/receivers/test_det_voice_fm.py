#!/usr/bin/env python3
"""Narrowband FM voice receiver (voice_fm): detection, channel plan, CTCSS, DCS, audio.

  * DCS convention checked against published code equivalences (023N = 340N = 766N = 047I),
    which the reciprocal Golay polynomial does not reproduce.
  * synthetic NBFM transmitters (synth.py): speech-like audio, CTCSS tones, DCS codes of both
    polarities, offsets, keying, noise; the demodulated audio of a 1 kHz tone.
  * a real PMR446 walkie-talkie off the Pluto (testdata/voice), hard-clipped: one carrier,
    channel 8, DCS 053N family, its clipping products suppressed.
  * no detection on a real LoRa capture, on noise, on an unmodulated carrier, or on 2-FSK data.

Usage: pytest -q tests/receivers/test_det_voice_fm.py   (or run it as a script)
"""
from __future__ import annotations
import json, os

import numpy as np

from spectrum_intel.receivers.voice_fm import VoiceDemodulator, detect, match_channel
from spectrum_intel.receivers.voice_fm import subaudio as SA, synth
from _data import DATA

FS = 250e3
FC = 446.1e6


def _one(x, fs=FS, fc=FC):
    em = detect(x, fs, fc)
    assert len(em) == 1, [(e.label, e.center_mhz) for e in em]
    return em[0]


def test_dcs_convention_matches_published_equivalences():
    fam = set(SA.dcs_aliases(SA.dcs_encode("023")))
    assert {("023", "N"), ("340", "N"), ("766", "N"), ("047", "I")} <= fam
    for code in SA.DCS_CODES:                                   # every standard code is a valid word
        w = SA.dcs_encode(code)
        assert SA.golay_syndrome(w) == 0 and (code, "N") in SA.dcs_aliases(w)


def test_channel_plan():
    assert match_channel(446.09375e6)["channel"] == 8
    c = match_channel(446.0909e6, ppm=7.4)                      # Pluto +7.4 ppm reads 3.3 kHz low
    assert c["service"] == "PMR446" and c["channel"] == 8 and abs(c["offset_hz"] + 2850) < 10
    assert match_channel(433.100e6)["service"] == "LPD433" and match_channel(433.100e6)["channel"] == 2
    assert match_channel(462.5625e6)["service"] == "FRS/GMRS"
    assert match_channel(446.1e6) is None                       # half-way between channels 8 and 9


def test_ctcss_tones_and_channels():
    for k, (tone, ch) in enumerate([(67.0, 1), (88.5, 8), (151.4, 12), (254.1, 16)]):
        f = 446.00625e6 + 12.5e3 * (ch - 1)
        e = _one(synth.nbfm(FS, 0.8, f - FC, ctcss=tone, snr_db=20, seed=k))
        F = e.features
        assert F["ctcss_hz"] == tone and F["channel"] == ch and F["service"] == "PMR446", F
        assert abs(e.center_mhz * 1e6 - f) < 100 and e.confidence >= 0.9


def test_dcs_codes_both_polarities():
    for k, (code, inv) in enumerate([("023", False), ("053", False), ("411", True), ("754", False)]):
        e = _one(synth.nbfm(FS, 0.8, -40e3 + 10e3 * k, dcs=code, dcs_inverted=inv, snr_db=20, seed=k))
        F = e.features
        want = f"D{code}{'I' if inv else 'N'}"
        assert want == F["dcs_code"] or want in F["dcs_aliases"], (want, F["dcs_code"], F["dcs_aliases"])
        assert F["dcs_words_verified"] >= 10


def test_keyed_burst_and_deviation():
    x = synth.nbfm(FS, 0.8, 30e3, audio=1000.0, audio_dev=2000, ctcss=100.0, sub_dev=300, snr_db=25,
                   keyed=(0.2, 0.6))
    e = _one(x)
    F = e.features
    assert abs(F["keyed_s"] - 0.4) < 0.03 and abs(e.start_us / 1e6 - 0.2) < 0.02
    assert 1900 < F["deviation_hz"] < 2700 and F["ctcss_hz"] == 100.0


def test_audio_tone_and_squelch():
    x = synth.nbfm(FS, 0.8, 30e3, audio=1000.0, audio_dev=2000, ctcss=100.0, snr_db=25, keyed=(0.0, 0.4))
    d = VoiceDemodulator(FC + 30e3, deviation_hz=2500)
    a = d.process(x, FS, FC)
    s = a.samples
    on, off = s[int(0.05 * a.fs):int(0.35 * a.fs)], s[int(0.45 * a.fs):]
    P = np.abs(np.fft.rfft(on * np.hanning(len(on)))) ** 2
    f = np.fft.rfftfreq(len(on), 1 / a.fs)
    assert abs(f[np.argmax(P)] - 1000) < 15                   # the tone, not the 100 Hz CTCSS
    assert P[np.abs(f - 100) < 5].max() < 1e-3 * P.max()      # CTCSS removed
    assert 0.2 < np.sqrt(2) * on.std() < 0.6                  # fixed gain: 2 kHz of 2.5 -> ~0.4
    assert np.abs(off).max() < 1e-3                           # squelch closed after the carrier drops
    assert 0.45 < a.meta["squelch_open"] < 0.55


def test_real_pmr446_walkie_talkie():
    meta = json.load(open(os.path.join(DATA, "voice", "real_pmr446_ch8_dcs.json")))
    x = np.fromfile(os.path.join(DATA, "voice", "real_pmr446_ch8_dcs.cf32"), np.complex64)
    e = _one(x, meta["fs_hz"], meta["fc_hz"])                 # clipping products are not extra carriers
    gt, F = meta["ground_truth"], e.features
    assert F["service"] == gt["service"] and F["channel"] == gt["channel"]
    assert F["dcs_code"] in gt["dcs_family"] and F["dcs_words_verified"] >= 10
    assert len(F["suppressed_artifacts"]) >= 3 and e.confidence >= 0.9


def test_no_false_alarms():
    rng = np.random.default_rng(3)
    n = int(0.8 * FS)
    noise = ((rng.standard_normal(n) + 1j * rng.standard_normal(n)) * 0.01).astype(np.complex64)
    assert detect(noise, FS, FC) == []
    cw = noise + np.exp(2j * np.pi * 20e3 * np.arange(n) / FS).astype(np.complex64)
    assert detect(cw, FS, FC) == []                           # an unmodulated carrier is CW, not voice
    bits = rng.integers(0, 2, int(0.8 * 2400) + 1) * 2 - 1    # 2400 Bd 2-FSK, +-3 kHz: data, not voice
    fsk = np.exp(2j * np.pi * np.cumsum(-25e3 + 3e3 * bits[(np.arange(n) / FS * 2400).astype(int)]) / FS)
    assert detect((fsk + noise).astype(np.complex64), FS, FC) == []
    meta = json.load(open(os.path.join(DATA, "lora", "real_434MHz_sf9.json")))
    x = np.fromfile(os.path.join(DATA, "lora", "real_434MHz_sf9.cf32"), np.complex64)
    assert detect(x, meta["fs_hz"], meta["fc_hz"]) == []


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn(); print("PASS", name)
