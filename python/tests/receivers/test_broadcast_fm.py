"""FM broadcast + RDS (research, signal gaps S20) on synthetic stations."""
import numpy as np

from spectrum_intel.core.dsp import Channelizer
from spectrum_intel.receivers.broadcast_fm import detect
from spectrum_intel.receivers.broadcast_fm import rds as R
from spectrum_intel.receivers.broadcast_fm import synth as S

WORKED = ["01010010000000010100010001", "00000001010001001000101110", "11100000000000010111101110",
          "01000001010010010100000010"]                     # research file 06, group 0A PI 5201 "AI"


def test_rds_block_code_worked_example():
    bits = R.encode_group((0x5201, 0x0144, 0xE001, 0x4149))
    assert ["".join(map(str, bits[26 * i:26 * i + 26])) for i in range(4)] == WORKED
    assert [R.syndrome(int(w, 2)) for w in WORKED] == [0x0FC, 0x198, 0x168, 0x1B4]


def test_rds_sync_groups_and_fields():
    rng = np.random.default_rng(1)
    bits = list(rng.integers(0, 2, 17))
    for g in S.ps_groups(0x5201, 10, "AIR FM"):
        bits += R.encode_group(g)
    # 2A RadioText "HELLO" + 4A clock time 2026-10-25 10:30 UTC, +5:30
    bits += R.encode_group((0x5201, (2 << 12) | (10 << 5) | 0, 0x4845, 0x4C4C))
    bits += R.encode_group((0x5201, (2 << 12) | (10 << 5) | 1, 0x4F0D, 0x2020))
    mjd = 61338                                                 # 2026-10-25
    c = ((mjd & 0x7FFF) << 1) | (10 >> 4)
    d = ((10 & 15) << 12) | (30 << 6) | 11
    bits += R.encode_group((0x5201, (4 << 12) | (10 << 5) | (mjd >> 15), c, d))
    bits[300] ^= 1                                              # one damaged block is skipped, not mis-read
    info = R.decode(R.groups(R.blocks(bits)))
    assert info["pi"] == "5201" and info["pi_country_nibble"] == 5 and info["pty_name"] == "Pop music"
    assert info["ps"].replace("_", " ").strip() in ("AIR FM", "AIR")
    assert info["radiotext"] == "HELLO"
    assert info["clock_time"] == {"utc": "2026-10-25T10:30Z", "local_offset_h": 5.5}


def test_stations_stereo_mono_and_rds():
    fs, fc = 2.4e6, 98.0e6
    x = (S.fm_station(fs, 0.5, 300e3, S.ps_groups(0x5201, 10, "AIR FM"), seed=1)
         + S.fm_station(fs, 0.5, -500e3, None, stereo=False, seed=2))
    em = {e.features["channel_mhz"]: e for e in detect(Channelizer(x, fs, fc))}
    assert set(em) == {97.5, 98.3}
    assert not em[97.5].features["stereo"] and "rds" not in em[97.5].features
    f = em[98.3].features
    assert f["stereo"] and f["rds"]["pi"] == "5201" and f["rds"]["ps"] == "AIR FM  " and 50 < f["deviation_khz"] < 80


def test_no_station_in_noise_and_short_block():
    fs, fc = 2.4e6, 98.0e6
    rng = np.random.default_rng(3)
    n = 1 << 17
    x = ((rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)).astype(np.complex64)
    assert detect(Channelizer(x, fs, fc)) == []
    # 34 ms at 61.44 MSps: the station, no RDS (too short)
    y = S.fm_station(61.44e6, 0.034, 5e6, S.ps_groups(0x5201, 10, "AIR FM"), seed=4)
    em = detect(Channelizer(y, 61.44e6, 95.0e6))
    assert [e.features["channel_mhz"] for e in em] == [100.0] and "rds" not in em[0].features


def test_pipeline_label():
    from spectrum_intel.analyze import analyze_iq
    from spectrum_intel.decode.registry import plugins
    fs, fc = 2.4e6, 98.0e6
    x = S.fm_station(fs, 0.5, 300e3, S.ps_groups(0x5201, 10, "AIR FM"), seed=1)
    a = analyze_iq(x, fs, fc, plugins=plugins(enable=["broadcast_fm"]))
    s = [s for s in a.signals if s.features.get("decoder") == "broadcast_fm"]
    assert [x.label for x in s] == ["FM broadcast (98.3 MHz, stereo, RDS 'AIR FM')"]
    assert s[0].decoded["status"] == "success"
