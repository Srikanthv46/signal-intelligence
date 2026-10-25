"""GPS L1 C/A acquisition and GNSS indicators (research, signal gaps S3) on synthetic skies."""
import numpy as np
import pytest

from spectrum_intel.core.dsp import Channelizer
from spectrum_intel.receivers.gnss import synth as S
from spectrum_intel.receivers.gnss.acquire import analyse, chi2_sf_even, threshold
from spectrum_intel.receivers.gnss.codes import FIRST10, first10_octal

FS, FC = 4e6, 1575.42e6 - 0.3e6


def test_codes_match_icd_octals():
    assert all(first10_octal(p) == FIRST10[p] for p in range(1, 33))


def test_chi2_threshold():
    from scipy.stats import chi2
    for m in (1, 4, 8):
        g = threshold(m, 1e-7)
        assert chi2.sf(2 * m * g, 2 * m) == pytest.approx(1e-7, rel=1e-6)
        assert chi2_sf_even(10.0, 2 * m) == pytest.approx(chi2.sf(10.0, 2 * m), rel=1e-12)


def _sky(seed, eps=9.5, cn0=(48, 44, 40)):
    rng = np.random.default_rng(seed)
    sats = [{"prn": p, "cn0": c, "doppler": float(rng.uniform(-4000, 4000)), "code_phase": float(rng.uniform(0, 1023))}
            for p, c in zip((3, 7, 12, 19, 25), cn0)]
    return sats, S.sky(FS, 0.034, FC, sats, eps_ppm=eps, seed=seed)


@pytest.mark.parametrize("ppm", [9.5, None])
def test_acquisition_code_doppler_cn0(ppm):
    sats, x = _sky(1)
    r = analyse(Channelizer(x, FS, FC), ppm)
    got = {s["prn"]: s for s in r["satellites"]}
    assert {3, 7} <= set(got) and set(got) <= {3, 7, 12}
    off = 0.0 if ppm is not None else -9.5e-6 * 1575.42e6       # without the ppm the Doppler is apparent
    if ppm is None:
        assert abs(r["receiver_offset_hz_estimate"] - off) < 5e3
    for s in sats[:2]:
        g = got[s["prn"]]
        assert abs((g["code_phase_chips"] - s["code_phase"] + 511.5) % 1023 - 511.5) < 0.5
        assert abs(g["doppler_hz"] - s["doppler"] - off) < 60
        assert abs(g["cn0_dbhz"] - s["cn0"]) < 1.5
    assert not r["indicators"]["spoofing"]["verdict"] and not r["indicators"]["jamming"]["verdict"]


def test_noise_only_finds_nothing():
    x = S.sky(FS, 0.034, FC, [], seed=2)
    r = analyse(Channelizer(x, FS, FC), 9.5)
    assert r["satellites"] == []


def test_spoofed_sky_uniform_and_strong():
    sats = [{"prn": p, "cn0": 52.0, "doppler": d, "code_phase": c}
            for p, d, c in ((2, 300, 10), (5, 900, 200), (9, -1200, 400), (14, 2100, 600), (22, -2500, 800))]
    x = S.sky(FS, 0.034, FC, sats, eps_ppm=0.0, seed=3)
    r = analyse(Channelizer(x, FS, FC), 0.0)
    sp = r["indicators"]["spoofing"]
    assert r["n_satellites"] == 5 and sp["s2_uniform_cn0"] and not sp["s1_cn0_implausible"]
    # a meaconed / spoofed copy of each satellite 6 chips later -> multipeak
    sats2 = sats + [dict(s, code_phase=s["code_phase"] - 6.0, cn0=49.0) for s in sats]
    r = analyse(Channelizer(S.sky(FS, 0.034, FC, sats2, seed=4), FS, FC), 0.0)
    sp = r["indicators"]["spoofing"]
    assert sp["s4_multipeak"] and sp["s2_uniform_cn0"] and sp["verdict"]


def test_cw_jammer_and_pipeline():
    from spectrum_intel.analyze import analyze_iq
    from spectrum_intel.decode.registry import plugins
    n = int(0.034 * FS)
    t = np.arange(n) / FS
    rng = np.random.default_rng(5)
    x = np.sqrt(0.5) * (rng.standard_normal(n) + 1j * rng.standard_normal(n))
    x = (x + 3.0 * np.exp(2j * np.pi * (1575.42e6 + 3e3 - FC) * t)).astype(np.complex64)
    r = analyse(Channelizer(x, FS, FC), 0.0)
    assert r["indicators"]["jamming"]["verdict"] and r["interference"]["lines"]
    sats, sky = _sky(6, eps=0.0)
    a = analyze_iq(sky, FS, FC, plugins=plugins(enable=["gnss"]))
    g = [s for s in a.signals if s.features.get("decoder") == "gnss"]
    assert len(g) == 1 and g[0].label.startswith("GNSS: GPS L1 C/A (")


def test_swept_ppd_jammer():
    from spectrum_intel.receivers.gnss.acquire import swept_ppd
    fs, fc = 61.44e6, 1575.42e6 - 3e6
    for per, span, ctr, amp in ((9e-6, 20e6, 2e6, 3.0), (27e-6, 30e6, -5e6, 1.0), (12e-6, 15e6, 10e6, 1.0)):
        x = S.ppd(fs, 0.017, 1575.42e6 + ctr - fc, span, per, amp=amp)
        r = swept_ppd(Channelizer(x, fs, fc), 1575.42e6)
        assert r and abs(r["period_us"] - per * 1e6) < 0.1 and r["shape"] == "sawtooth up"
        assert abs(r["centre_offset_mhz"] - ctr / 1e6) < 0.5 and r["covers_l1"] == (abs(ctr) < span / 2)
    for sd in range(5):
        assert swept_ppd(Channelizer(S.sky(fs, 0.017, fc, [], seed=sd), fs, fc), 1575.42e6) is None
    x = S.ppd(fs, 0.017, 2e6, 20e6, 9e-6, amp=3.0)
    r = analyse(Channelizer(x, fs, 1575.42e6 - 2e6 + 2e6), None)
    assert r["indicators"]["jamming"]["swept"] and r["indicators"]["jamming"]["verdict"]
