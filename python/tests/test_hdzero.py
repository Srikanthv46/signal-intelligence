"""HDZero 27 MHz-mode detector (receivers/digital_vtx/hdzero.py): slice plan, a synthetic link, and the real capture
testdata/live/LIVE_HDZERO_5658MHz (HDZero Freestyle V2, R1, user-confirmed)."""
from pathlib import Path

import numpy as np
import pytest
from scipy.signal import resample_poly

from spectrum_intel.core import common as C
from spectrum_intel.receivers.digital_vtx import detector as D, hdzero as HZ

LIVE = Path(__file__).resolve().parents[2] / "testdata" / "live" / "LIVE_HDZERO_5658MHz" / "LIVE_HDZERO_5658MHz.sigmf-meta"


def test_slice_plan():
    assert HZ.slice_plan(2097152, 61.44e6) == (921600, 1024)          # 2^21 block: exact 27 MSps
    assert HZ.slice_plan(65536, 61.44e6) == (28800, 1024)
    m, nw = HZ.slice_plan(2000000, 61.44e6)                           # signal-vault block: 878906.25 -> fast length
    assert HZ._smooth(m) and abs(1024 * m / (2000000 * 27e6 / 61.44e6) - nw) <= HZ.NW_TOL


def _synth_capture(snr_db, cfo_hz=20e3, n_sym=90, seed=3):
    g = np.random.default_rng(seed)
    y = resample_poly(HZ.synth(n_sym, g), 256, 225)                   # 27 -> 30.72 MSps
    fs, fc = 30.72e6, 5657e6                                          # R1 inside the usable +-15.36 MHz
    off = 5658e6 - fc + cfo_hz
    y = y * np.exp(2j * np.pi * off * np.arange(len(y)) / fs)
    p = np.mean(np.abs(y) ** 2) * fs / (889 * HZ.SCS)                  # in-band power density -> noise for snr_db
    x = y + (g.standard_normal(len(y)) + 1j * g.standard_normal(len(y))) * np.sqrt(p / 10 ** (snr_db / 10) / 2)
    return x.astype(np.complex64), fs, fc


@pytest.mark.parametrize("snr", [15, 0])
def test_synthetic_link(snr):
    x, fs, fc = _synth_capture(snr)
    ctx = {"report_wifi": False, "generic": False}
    em = D.detect(x, fs, fc, ctx)
    hz = [e for e in em if e.family == "hdzero"]
    assert len(hz) == 1
    F = hz[0].features
    assert F["channel"] == "R1" and abs(F["cfo_hz"] - 20e3) < 300
    assert F["preamble_ms"] and abs(F["preamble_ms"][0] - 30 * HZ.TS / HZ.FS * 1e3) < 0.1


def test_noise_only():
    g = np.random.default_rng(5)
    x = (g.standard_normal(98304) + 1j * g.standard_normal(98304)).astype(np.complex64)
    assert not [e for e in D.detect(x, 30.72e6, 5657e6, {"generic": False}) if e.family == "hdzero"]


@pytest.mark.skipif(not LIVE.exists(), reason="testdata not present")
def test_real_capture():
    x, fs, fc, _ = C.load_capture(str(LIVE))
    ctx = {"generic": False}
    em = [e for e in D.detect(x, fs, fc, ctx) if e.family == "hdzero"]
    assert len(em) == 1
    F = em[0].features
    assert F["channel"] == "R1" and 6.0 < F["cfo_ppm"] < 8.0
    assert F["frame_symbols"] == 426 and abs(F["frame_rate_fps"] - 60.0) < 0.1
    assert F["pilot_coherence"] > 0.8 and F["data_coherence"] < 0.2


def test_wider_search_bands_and_context():
    """PLAN G8: the blind OFDM search covers 1.35-1.444, 1.6-2.53 and 4.4-6.0 GHz; a link outside the ISM / FPV bands
    carries its band context; cellular bands are recognised (the blind search skips candidates there); PLAN G9: the
    analog ranges and channel table include the 3.3 GHz groups."""
    from spectrum_intel.core import bands_india as BI
    from spectrum_intel.receivers.analog_vtx import match_channel
    assert BI.uav_link_band(1437e6).startswith("1.4 GHz UAV data-link band")
    assert BI.uav_link_band(1380e6).startswith("L-band")
    assert BI.uav_link_band(2250e6).startswith("1.6-2.33 GHz")
    assert BI.uav_link_band(4600e6).startswith("4.4-4.9 GHz")
    assert BI.uav_link_band(5800e6) is None
    assert BI.cellular_at(1840e6) and BI.cellular_at(2350e6)
    assert not BI.cellular_at(2250e6) and not BI.cellular_at(1437e6)
    assert any(lo <= 3.3e9 <= hi for lo, hi in BI.ANALOG_VTX_RANGES)
    assert match_channel(3409.0)[0][0] == "3G3A5" and match_channel(3171.0)[0][0] == "3G3B1"
