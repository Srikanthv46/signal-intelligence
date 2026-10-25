"""AM voice (research, signal gaps S13): aeronautical AM detected and demodulated, FM voice stays FM."""
import numpy as np

from spectrum_intel.receivers.voice_fm import detector as D
from spectrum_intel.receivers.voice_fm import synth as S
from spectrum_intel.receivers.voice_fm.audio import VoiceDemodulator

FS = 1.0e6


def test_am_air_band_detected_with_channel():
    fc = 121.2e6
    x = S.am(FS, 0.3, 300e3, depth=0.8, snr_db=25, seed=1)                # 121.5 MHz emergency
    em = D.detect(x, FS, fc)
    assert len(em) == 1 and em[0].label == "am_voice"
    f = em[0].features
    assert f["service"] == "VHF air band" and f["emergency"].startswith("aeronautical emergency")
    assert f["channel_mhz"] == 121.5 and abs(f["offset_hz"]) < 100 and f["modulation_rms"] > 0.1


def test_fm_voice_not_called_am_and_am_not_fm():
    fc = 446.0e6
    fm = S.nbfm(FS, 0.3, 6.25e3, audio="speech", ctcss=88.5, snr_db=25, seed=2)
    em = D.detect(fm, FS, fc)
    assert em and all(e.label != "am_voice" for e in em)
    am = S.am(FS, 0.3, -300e3, depth=0.7, snr_db=25, seed=3)              # 300.7 MHz: UHF military air band
    em = D.detect(am, FS, 301.0e6)
    assert [e.label for e in em] == ["am_voice"] and em[0].features["service"] == "UHF military air band"


def test_unmodulated_carrier_is_not_am():
    t = np.arange(int(0.3 * FS)) / FS
    rng = np.random.default_rng(4)
    x = np.exp(2j * np.pi * 100e3 * t) + 0.02 * (rng.standard_normal(len(t)) + 1j * rng.standard_normal(len(t)))
    assert D.detect(x.astype(np.complex64), FS, 125e6) == []


def test_am_audio_demodulation_tone():
    x = S.am(FS, 0.3, 50e3, depth=0.5, audio=1000.0, snr_db=30, seed=5)
    d = VoiceDemodulator.from_fields({"center_mhz": 125.05, "audio": "am"})
    ch = d.process(x, FS, 125.0e6)
    a = ch.samples[len(ch.samples) // 4:]
    spec = np.abs(np.fft.rfft(a * np.hanning(len(a))))
    f = np.fft.rfftfreq(len(a), 1 / ch.fs)
    assert abs(f[np.argmax(spec)] - 1000) < 20
    assert 0.15 < np.max(np.abs(a)) < 0.4                                   # 50 % depth -> ~0.25 FS


def test_pipeline_am_signal_label():
    from spectrum_intel.analyze import analyze_iq
    from spectrum_intel.decode.registry import plugins
    x = S.am(2e6, 0.2, 250e3, depth=0.8, snr_db=25, seed=6)                 # 118.25 MHz
    a = analyze_iq(x, 2e6, 118.0e6, plugins=plugins(enable=["voice_fm"]))
    s = [s for s in a.signals if s.label.startswith("Analog AM voice")]
    assert len(s) == 1 and s[0].label == "Analog AM voice (VHF air band, 118.25 MHz)"
    assert s[0].decoded["fields"]["audio"] == "am"
