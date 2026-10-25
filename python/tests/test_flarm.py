"""FLARM beacons, detection and clear header only (research, PLAN T5; research/drone_transmitters/
10_cc2500_a7105_xn297_flarm_xbee.md section 5): Manchester, sync, radio ID, address type; rc_link detect."""
import numpy as np

from spectrum_intel.receivers.rc_link import detector as D, flarm as F, synth as S
from spectrum_intel.receivers.rc_link.plugin import PLUGINS


def test_manchester_sync_matches_published_bytes():
    chips = F.manchester([(F.SYNC >> (31 - i)) & 1 for i in range(32)])
    by = bytes(int("".join(map(str, chips[8 * i:8 * i + 8])), 2) for i in range(8))
    assert by.hex() == "5599a5a955666596"                        # the Manchester sync as published


def test_identify():
    rng = np.random.default_rng(1)
    iq, _ = S.flarm_burst(0xDE088B, 0x10, fs=8e5, rng=rng)
    for y in (iq, np.conj(iq)):
        d = F.try_identify(y, 8e5)
        assert d["radio_id"] == "DE088B" and d["address_type"] == "ICAO 24-bit address"
    assert F.try_identify(rng.normal(size=8000) + 1j * rng.normal(size=8000), 8e5) is None


def test_detect():
    rng = np.random.default_rng(2)
    fs, fc = 7.68e6, 868.0e6
    x = S.noise_capture(1 << 18, rng)
    bursts = []
    for k, f in enumerate((868.2e6, 868.4e6)):
        iq, _ = S.flarm_burst(0x123456, cfo_hz=3e3, fs=fs, rng=rng)
        bursts.append(dict(t0_us=2000 + 12000 * k, f_hz=f, iq=iq, bw_hz=200e3))
    y, _ = S.embed(x, fs, fc, bursts, 14)
    em = [e for e in D.detect(y, fs, fc) if e.label == "flarm_like"]
    assert len(em) == 2
    out = PLUGINS[0].summarise(em, None)
    assert out["fields"]["flarm"]["radio_ids"] == ["123456"] and out["fields"]["model"]["family"] == "FLARM device"
    assert PLUGINS[0].pretty("flarm_like", em).startswith("FLARM-like")
