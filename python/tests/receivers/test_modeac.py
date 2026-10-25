"""SSR Mode A/C replies (research, signal gaps S11): Gillham against pyModeS, synthetic replies, pipeline."""
import numpy as np

from spectrum_intel.receivers.adsb import modeac as M

FS = 8e6


def _env(replies, n=40000, noise=0.02, seed=1):
    rng = np.random.default_rng(seed)
    env = np.abs(noise * (rng.standard_normal(n) + 1j * rng.standard_normal(n)))
    for t0, code, amp, spi in replies:
        r = M.synth_reply(M.bits_for_code(code), FS, amp, spi)
        i = int(t0 * FS)
        env[i:i + len(r)] = np.maximum(env[i:i + len(r)], r)
    return env


def test_gillham_matches_pymodes_everywhere():
    from pyModeS import common as c
    keys = ["D2", "D4", "A1", "A2", "A4", "B1", "B2", "B4", "C1", "C2", "C4"]
    for v in range(2 ** 11):
        bits = [(v >> (10 - i)) & 1 for i in range(11)]
        p = dict.fromkeys(M.POSITIONS, 0)
        p.update(zip(keys, bits))
        assert M.gillham_ft(p) == c.gray2alt("".join(map(str, bits)))


def test_replies_codes_spi_and_repetition():
    env = _env([(100e-6, "7700", 1.0, False), (400e-6, "7700", 0.9, False), (900e-6, "1234", 0.5, True),
                (1500e-6, "0000", 0.6, False)])
    rep = M.decode(env, FS)
    assert [r["code"] for r in rep] == ["7700", "7700", "1234", "0000"]
    assert [r["spi"] for r in rep] == [False, False, True, False]
    s = M.summarise(rep)
    top = s["replies"][0]
    assert top["code"] == "7700" and top["n"] == 2 and top["repeated"] and top["special_code"] == "emergency"
    assert s["n_distinct"] == 3 and s["n_repeated"] == 1


def test_garbled_overlap_rejected_and_altitude():
    # two replies 7 us apart overlap: off-grid pulses between the framing pulses -> both refused
    env = _env([(100e-6, "1234", 1.0, False), (107e-6, "4321", 1.0, False)])
    assert M.decode(env, FS) == []
    # Mode C reading of a valid Gillham word
    p = dict.fromkeys(M.POSITIONS, 0)
    p["C2"] = 1
    assert M.gillham_ft(p) == -1000


def test_mode_s_frames_do_not_make_replies():
    rng = np.random.default_rng(3)
    n = 20000
    env = np.abs(0.02 * (rng.standard_normal(n) + 1j * rng.standard_normal(n)))
    # 112-bit PPM frame: 0.5 us pulses on a 0.5 us grid, never on the 1.45 us grid for 20.3 us
    t0 = 200e-6
    for k, b in enumerate(rng.integers(0, 2, 112)):
        a = int((t0 + 8e-6 + k * 1e-6 + (0 if b else 0.5e-6)) * FS)
        env[a:a + 4] = 1.0
    assert M.decode(env, FS, busy=[(t0, t0 + 120e-6)]) == []
    assert M.decode(env, FS) == []


def test_pipeline_mode_ac_only_signal():
    from spectrum_intel.analyze import analyze_iq
    from spectrum_intel.decode.registry import plugins
    fs, fc = 10e6, 1088e6
    n = 2 ** 17
    rng = np.random.default_rng(5)
    x = (0.01 * (rng.standard_normal(n) + 1j * rng.standard_normal(n))).astype(np.complex64)
    t = np.arange(n) / fs
    for t0, code in ((1e-3, "2000"), (4e-3, "2000"), (7e-3, "2000")):
        b = M.bits_for_code(code)
        env = M.synth_reply(b, fs)
        i = int(t0 * fs)
        x[i:i + len(env)] += (env * np.exp(2j * np.pi * 2e6 * t[i:i + len(env)])).astype(np.complex64)
    a = analyze_iq(x, fs, fc, plugins=plugins(enable=["adsb"]))
    sig = [s for s in a.signals if s.label.startswith("SSR Mode A/C")]
    assert len(sig) == 1
    f = sig[0].decoded["fields"]["mode_ac"]
    assert f["replies"][0]["code"] == "2000" and f["replies"][0]["n"] == 3
    assert sig[0].label == "SSR Mode A/C replies" and sig[0].decoded["status"] == "partial"
