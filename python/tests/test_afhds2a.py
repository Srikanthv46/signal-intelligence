"""FlySky AFHDS 2A (A7105; research, PLAN R7; research/drone_transmitters/10_cc2500_a7105_xn297_flarm_xbee.md):
Hamming FEC, candidate conventions, synthetic packets through try_identify and rc_link detect -> IDs, label, field."""
import numpy as np
import pytest

from spectrum_intel.receivers.rc_link import afhds2a as A, detector as D, synth as S
from spectrum_intel.receivers.rc_link.plugin import PLUGINS

FS, FC = 61.44e6, 2440e6


def _payload(rng, typ=0x58, rx=bytes([5, 6, 7, 8])):
    return bytes([typ, 1, 2, 3, 4]) + rx + bytes(rng.integers(0, 256, 29).tolist())


def test_hamming_corrects_one_bit():
    for lay in A.FEC_LAYOUTS:
        for n in range(16):
            for e in range(7):
                assert A._DEC[lay][A._ENC[lay][n] ^ (1 << e)] == n


@pytest.mark.parametrize("conv", [dict(), dict(layout="classic", high_first=False, crc=(0x8005, 0xFFFF))])
def test_identify(conv):
    rng = np.random.default_rng(1)
    iq, _ = S.afhds2a_burst(_payload(rng), fs=4e6, rng=rng, **conv)
    for y in (iq, np.conj(iq)):
        d = A.try_identify(y, 4e6)
        assert d["crc_ok"] and d["tx_id"] == "01020304" and d["rx_id"] == "05060708" and d["packet_type"] == "0x58"
    assert d["convention"]["fec"] == conv.get("layout", "systematic")
    assert A.try_identify(rng.normal(size=20000) + 1j * rng.normal(size=20000), 4e6) is None


def test_id_code_only():
    """A packet whose FEC / CRC follow none of the candidates still gives the ID-code identification."""
    rng = np.random.default_rng(2)
    bits = np.r_[np.arange(32) % 2, [(A.ID_CODE >> (31 - i)) & 1 for i in range(32)], rng.integers(0, 2, 560),
                 np.zeros(8)].astype(np.uint8)
    iq, _ = S.gfsk_burst(500e3, 186e3, bits=bits, bt=2.0, fs=4e6, rng=rng)
    d = A.try_identify(iq, 4e6)
    assert d["crc_ok"] is None and d["protocol"].startswith("FlySky AFHDS 2A-like")


def test_detect_bind():
    rng = np.random.default_rng(3)
    x = S.noise_capture(1 << 21, rng)
    bursts = []
    for k, f in enumerate((2421.5e6, 2433.0e6, 2452.5e6)):            # 3.85 ms packets on the 500 kHz grid
        iq, _ = S.afhds2a_burst(_payload(rng, 0xBB, b"\xff" * 4), cfo_hz=8e3, rng=rng)
        bursts.append(dict(t0_us=1000 + 3850 * k, f_hz=f, iq=iq, bw_hz=1.2e6))
    y, _ = S.embed(x, FS, FC, bursts, 14)
    em = [e for e in D.detect(y, FS, FC) if e.label.startswith("afhds2a") and e.decoded]
    assert len(em) == 3 and all(e.label == "afhds2a" for e in em)
    f = PLUGINS[0].summarise(em, None)["fields"]["afhds2a"]
    assert f["crc_ok"] == 3 and f["tx_ids"] == ["01020304"] and f["binding"] and f["packet_types"] == ["0xBB"]
    assert PLUGINS[0].pretty("afhds2a", em) == "FlySky AFHDS 2A (A7105)"
