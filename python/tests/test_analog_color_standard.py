"""analog_vtx colour standard (research: research/drone_transmitters/01_analog_video.md) on synthetic FM VTX signals:
NTSC / PAL / PAL-M / PAL-N with colour, NTSC / PAL monochrome, and below the FM threshold."""
import numpy as np
import pytest

from spectrum_intel.receivers.analog_vtx import detector as D, synth as S

FS, N = 61.44e6, 1 << 21


def run(std, color=True, ppm=0.0, snr_db=None, seed=5):
    rng = np.random.default_rng(seed)
    x = S.synth_vtx(N, FS, 3e6, std=std, dev_hz=4e6, clock_ppm=ppm, rng=rng, color=color)
    if snr_db is not None:                                   # in-band SNR over 20 MHz
        s = np.sqrt(FS / 20e6 * 10 ** (-snr_db / 10) / 2)
        x = x + (rng.standard_normal(N) + 1j * rng.standard_normal(N)).astype(np.complex64) * s
    em = D.detect(x, FS, 5.8e9)
    assert em
    return em[0]


@pytest.mark.parametrize("std,ppm", [("ntsc", 30), ("pal", -40), ("pal_m", 15), ("pal_n", -20)])
def test_color_standard(std, ppm):
    e = run(std, ppm=ppm, snr_db=20)
    F = e.features
    assert F["color"] is True and F["color_standard"] == std
    assert e.label == (std if std in ("pal_m", "pal_n") else std)
    fsc = S.STANDARDS[std]["fsc"] * (1 + ppm * 1e-6)         # the camera crystal scales the subcarrier too
    assert abs(F["chroma_fsc_hz"] - fsc) < 30
    assert (F["burst_swing_deg"] > 25) == (std != "ntsc")


@pytest.mark.parametrize("std", ["ntsc", "pal"])
def test_monochrome(std):
    F = run(std, color=False, snr_db=20).features
    assert F["color"] is False and "color_standard" not in F


def test_low_snr_is_not_black_and_white():
    F = run("pal", snr_db=5).features
    assert F.get("color") in (None, True)
    assert F.get("color_standard") in (None, "pal")


def test_lowband_channel_table():
    assert D.match_channel(5413.0)[0][0] == "LB3"
    assert D.match_channel(5806.0)[0][0] == "R5"


@pytest.mark.parametrize("snr_db", [20, 10])
def test_secam(snr_db):
    """SECAM (no burst; FM chroma Db 4.250 / Dr 4.40625 MHz on alternate lines): line-to-line alternation of the
    back-porch subcarrier frequency (research)."""
    e = run("secam", snr_db=snr_db)
    F = e.features
    assert F["color"] is True and F["color_standard"] == "secam" and e.label == "secam"
    assert abs(F["secam_db_hz"] - 4.25e6) <= 10e3 and abs(F["secam_dr_hz"] - 4.40625e6) <= 10e3


@pytest.mark.parametrize("std", ["pal", "pal_n"])
def test_pal_is_not_secam(std):
    F = run(std, snr_db=10).features
    assert F["color_standard"] == std and (F.get("secam_test") or {}).get("alternation", 0) < 0.25


@pytest.mark.parametrize("tone,snr_db", [(440.0, 20), (1250.0, 12), (180.0, 20)])
def test_audio_features(tone, snr_db):
    """Audio FM subcarrier demodulated (research, PLAN A5): the dominant tone (a propeller's blade-pass on a drone)."""
    rng = np.random.default_rng(3)
    x = S.synth_vtx(N, FS, 3e6, std="pal", dev_hz=4e6, rng=rng, audio_tone_hz=tone)
    s = np.sqrt(FS / 20e6 * 10 ** (-snr_db / 10) / 2)
    x = x + (rng.standard_normal(N) + 1j * rng.standard_normal(N)).astype(np.complex64) * s
    F = D.detect(x, FS, 5.8e9)[0].features
    assert abs(F["audio_peak_hz"] - tone) <= 5 and abs(F["audio_subcarrier_hz"] - 6.5e6) <= 2e3
    assert 5 <= F["audio_dev_rms_khz"] <= 40 and F["audio_flatness"] < 0.6


def test_no_audio_fields_without_subcarrier():
    rng = np.random.default_rng(3)
    x = S.synth_vtx(N, FS, 3e6, std="pal", dev_hz=4e6, rng=rng, audio_hz=None)
    F = D.detect(x, FS, 5.8e9)[0].features
    assert F.get("audio_peak_hz") is None
