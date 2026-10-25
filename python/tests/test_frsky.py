"""FrSky D8 / D16 (CC2500, research, PLAN R6; research/drone_transmitters/10_cc2500_a7105_xn297_flarm_xbee.md):
CRCs, synthetic packets through try_identify and through rc_link detect -> variant, TX ID, plugin label and fields."""
import numpy as np
import pytest

from spectrum_intel.receivers.rc_link import cc2500 as CC, detector as D, synth as S
from spectrum_intel.receivers.rc_link.plugin import PLUGINS

FS, FC = 61.44e6, 2440e6
LEN = {"D8": 17, "D16 v1": 27, "D16 v2": 27, "D16 LBT": 30}


def test_crc_tables():
    assert CC._T[1] == 0x1189 and CC._T[7] == 0x74BF             # reflected-0x1021 (KERMIT) table
    assert CC.crc16_cc2500(b"") == 0xFFFF
    p = CC.d16_packet(bytes(range(27)), v2=True)
    assert len(p) == 32 and p[0] == 0x1D and CC.crc16_cc2500(p[:30]) == p[30] << 8 | p[31]


@pytest.mark.parametrize("variant", list(LEN))
def test_identify(variant):
    rng = np.random.default_rng(3)
    br = S.FRSKY_PHY[variant][0]
    pl = bytes([0x34, 0x12]) + bytes(rng.integers(0, 256, LEN[variant] - 2).tolist())
    iq, _ = S.frsky_burst(variant, pl, fs=8 * br, rng=rng)
    for y in (iq, np.conj(iq)):                                     # either FSK polarity
        d = CC.try_identify(y, 8 * br, br)
        assert d["crc_ok"] and d["tx_id"] == "0x1234"
        assert d["variant"] == variant.replace("D16 LBT", "D16 LBT v1")
    assert CC.try_identify(rng.normal(size=4000) + 1j * rng.normal(size=4000), 8 * br, br) is None


@pytest.mark.parametrize("variant,snr_db", [("D8", 15), ("D16 v2", 12)])
def test_detect(variant, snr_db):
    rng = np.random.default_rng(7)
    x = S.noise_capture(1 << 21, rng)
    bursts = []
    for k, f in enumerate((2421.5e6, 2431.7e6, 2458.2e6)):         # three hops, 9 ms apart
        pl = bytes([0x34, 0x12]) + bytes(rng.integers(0, 256, LEN[variant] - 2).tolist())
        iq, _ = S.frsky_burst(variant, pl, cfo_hz=5e3, rng=rng)
        bursts.append(dict(t0_us=1000 + 9000 * k, f_hz=f, iq=iq, bw_hz=250e3))
    y, _ = S.embed(x, FS, FC, bursts, snr_db)
    em = [e for e in D.detect(y, FS, FC) if e.label.startswith("frsky_") and e.decoded]   # + undecoded fragments
    assert len(em) == 3 and all(e.decoded["crc_ok"] for e in em)
    f = PLUGINS[0].summarise(em, None)["fields"]["frsky"]
    assert f["crc_ok"] == 3 and f["tx_ids"] == ["0x1234"] and f["variants"] == [variant]
    assert PLUGINS[0].pretty(em[0].label, em) == f"FrSky {variant} (ACCST)"
