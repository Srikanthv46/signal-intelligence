"""XN297 / Bayang toy drones (research, PLAN R8; research/drone_transmitters/10_cc2500_a7105_xn297_flarm_xbee.md):
scramble table, packet encode / decode, rc_link detect -> address, TX ID, label, fields, model."""
import numpy as np
import pytest

from spectrum_intel.receivers.rc_link import detector as D, synth as S, xn297 as X
from spectrum_intel.receivers.rc_link.plugin import PLUGINS

FS, FC = 61.44e6, 2440e6
TX = bytes([0x11, 0x22, 0x33, 0x44, 0x55])


def test_scramble_table_matches_published_start():
    assert X.SCRAMBLE[:8] == [0xE3, 0xB1, 0x4B, 0xEA, 0x85, 0xBC, 0xE5, 0x66]


@pytest.mark.parametrize("kind", ["bind", "data"])
def test_decode(kind):
    rng = np.random.default_rng(1)
    p = S.bayang_payload(kind, TX, (700, 512, 100, 300))
    iq, _ = S.xn297_burst(TX, p, fs=8e6, rng=rng)
    d = X.try_decode(iq, 8e6)
    assert d["checksum_ok"] and d["address"] == "1122334455" and d["packet_name"] == kind
    if kind == "bind":
        assert d["tx_id"] == "1122334455" and d["hop_channels"] == [0, 0x14, 0x34, 0x54]
    else:
        assert d["sticks"] == dict(aileron=700, elevator=512, throttle=100, rudder=300)
    assert X.try_decode(rng.normal(size=4000) + 1j * rng.normal(size=4000), 8e6) is None


def test_detect():
    rng = np.random.default_rng(5)
    x = S.noise_capture(1 << 21, rng)
    bursts = []
    for k, f in enumerate((2420.0e6, 2436.0e6, 2452.0e6)):        # Bayang hops every 2 ms
        iq, _ = S.xn297_burst(TX, S.bayang_payload("data", TX), cfo_hz=10e3, rng=rng)
        bursts.append(dict(t0_us=1000 + 2000 * k, f_hz=f, iq=iq, bw_hz=1.2e6))
    y, _ = S.embed(x, FS, FC, bursts, 20)
    em = [e for e in D.detect(y, FS, FC) if e.label == "bayang"]
    assert len(em) == 3
    out = PLUGINS[0].summarise(em, None)
    assert out["fields"]["xn297"]["addresses"] == ["1122334455"] and out["fields"]["xn297"]["bayang_checksum_ok"] == 3
    assert out["fields"]["model"]["family"] == "Bayang (XN297) toy drone"
    assert PLUGINS[0].pretty("bayang", em) == "Bayang (XN297)"
