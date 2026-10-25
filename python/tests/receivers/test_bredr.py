"""Bluetooth BR/EDR access codes (research, signal gaps S7): sync-word construction and LAP recovery."""
import numpy as np

from spectrum_intel.receivers.rc_link import bredr as B
from spectrum_intel.receivers.rc_link import detector as D
from spectrum_intel.receivers.rc_link import synth as S


def test_syncwords():
    assert B.syncword(0x000000) == 0xB0000002C7820E7E          # libbtbb DEFAULT_CODEWORD
    assert B.syncword(B.GIAC) == 0x4E7A2CCE331A3AE2
    for lap in (0, 0x123456, 0x800000, 0xFFFFFF):
        w = B.syncword(lap)
        assert (w >> 34) & 0xFFFFFF == lap and (w >> 58) in B.TAILS


def test_find_with_errors():
    bits = S.bredr_bits(0x5A5A5A, rng=np.random.default_rng(1)).tolist()
    for e in (0, 3):
        b = list(bits)
        for q in range(e):
            b[4 + 7 * q] ^= 1                                     # errors in the parity part of the sync word
        assert B.find(b) == (4, 0x5A5A5A, e)
    b = list(bits)
    for q in range(4):
        b[4 + 7 * q] ^= 1
    assert B.find(b) is None or B.find(b)[1] != 0x5A5A5A
    rng = np.random.default_rng(2)
    assert all(B.find(rng.integers(0, 2, 400).tolist()) is None for _ in range(200))


def test_detect_lap_and_inquiry():
    fs, fc = S.FS, 2440e6
    for sd, lap in ((0, B.GIAC), (1, 0x2C1F07)):
        rng = np.random.default_rng(sd)
        iq, _ = S.bredr_packet(lap, rng=rng)
        y, _ = S.embed(S.noise_capture(600_000, rng), fs, fc, [dict(iq=iq, f_hz=2431e6, t0_us=1000.0, bw_hz=1e6)], 20)
        em = [e for e in D.detect(y, fs, fc) if e.label == "bluetooth_bredr"]
        assert len(em) == 1 and em[0].decoded["lap"] == f"{lap:06X}" and em[0].decoded["inquiry"] == (lap == B.GIAC)


def test_random_gfsk_gives_no_lap():
    fs, fc = S.FS, 2440e6
    for sd in range(10):
        rng = np.random.default_rng(100 + sd)
        iq, _ = S.gfsk_burst(1e6, 160e3, n_bits=400, rng=rng)
        y, _ = S.embed(S.noise_capture(400_000, rng), fs, fc, [dict(iq=iq, f_hz=2431e6, t0_us=500.0, bw_hz=1e6)], 20)
        assert not any((e.decoded or {}).get("lap") for e in D.detect(y, fs, fc))


def test_pipeline_label():
    from spectrum_intel.analyze import analyze_iq
    from spectrum_intel.decode.registry import plugins
    fs, fc = S.FS, 2440e6
    rng = np.random.default_rng(3)
    iq, _ = S.bredr_packet(B.GIAC, rng=rng)
    y, _ = S.embed(S.noise_capture(600_000, rng), fs, fc, [dict(iq=iq, f_hz=2431e6, t0_us=1000.0, bw_hz=1e6)], 20)
    a = analyze_iq(y, fs, fc, plugins=plugins(enable=["rc_link"]))
    sig = [s for s in a.signals if s.features.get("decoder") == "rc_link"]
    assert [s.label for s in sig] == ["Bluetooth BR/EDR (LAP 9E8B33, inquiry)"]
    assert sig[0].decoded["fields"]["bredr"]["laps"] == ["9E8B33"]
