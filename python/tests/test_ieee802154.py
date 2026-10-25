"""IEEE 802.15.4 O-QPSK (XBee / Zigbee telemetry; research, PLAN T4; research/drone_transmitters/
10_cc2500_a7105_xn297_flarm_xbee.md section 6): chips, CRC, MAC parse, non-coherent demod vs SNR."""
import numpy as np
import pytest

from spectrum_intel.receivers.rc_link import ieee802154 as Z


def _xbee_frame(rng):
    fcf = 1 | 1 << 5 | 1 << 6 | 2 << 10 | 3 << 14              # data, ack req, intra-PAN, dst short, src extended
    m = bytes([fcf & 0xFF, fcf >> 8, 0x42, 0x34, 0x12, 0xFF, 0xFF]) + bytes.fromhex("0013a20040a1b2c3")[::-1] \
        + bytes(rng.integers(0, 256, 20).tolist())
    c = Z.crc16_kermit(m)
    return m + bytes([c & 0xFF, c >> 8])


def test_tables_and_crc():
    for s in range(1, 8):                                       # symbols 1-7: 4-chip cyclic shifts of symbol 0
        assert np.array_equal(Z.CHIP_ARR[s], np.roll(Z.CHIP_ARR[0], 4 * s))
    odd = np.arange(32) % 2 == 1                                # 8-15: 0-7 with the odd chips inverted
    assert all(np.array_equal(Z.CHIP_ARR[s + 8][odd], 1 - Z.CHIP_ARR[s][odd]) for s in range(8))
    assert Z.crc16_kermit(b"123456789") == 0x2189               # CRC-16/KERMIT check value
    assert Z.channel(2440e6) == 18 and Z.channel(2441e6) is None


@pytest.mark.parametrize("snr_db", [20, 6])
def test_demod(snr_db):
    rng = np.random.default_rng(2)
    m = _xbee_frame(rng)
    iq = Z.modulate(Z.chips_for(m))
    iq = np.r_[np.zeros(100), iq, np.zeros(100)] * np.exp(1j * 1.3)
    w = iq + (rng.normal(size=len(iq)) + 1j * rng.normal(size=len(iq))) * np.sqrt(10 ** (-snr_db / 10) / 2)
    d = Z.demod(w, 8e6)
    assert d["crc_ok"] and d["frame_type"] == "data" and d["dst_pan"] == "1234"
    assert d["src_addr"] == "0013a20040a1b2c3" and d["vendor"] == "Digi International (XBee)"
    assert Z.demod(rng.normal(size=8000) + 1j * rng.normal(size=8000), 8e6) is None


def test_detect_xbee():
    from spectrum_intel.receivers.rc_link import detector as D, synth as S
    from spectrum_intel.receivers.rc_link.plugin import PLUGINS
    FS, FC = 61.44e6, 2440e6
    rng = np.random.default_rng(9)
    x = S.noise_capture(1 << 21, rng)
    bursts = []
    for k in range(3):                                          # channel 17 (2435 MHz), three frames
        iq, _ = S.oqpsk_burst(_xbee_frame(rng), cfo_hz=15e3, rng=rng)
        bursts.append(dict(t0_us=1000 + 6000 * k, f_hz=2435e6, iq=iq, bw_hz=2.6e6))
    y, _ = S.embed(x, FS, FC, bursts, 12)
    em = [e for e in D.detect(y, FS, FC) if e.label == "ieee802154"]
    assert len(em) == 3 and all(e.decoded["channel"] == 17 for e in em)
    out = PLUGINS[0].summarise(em, None)
    f = out["fields"]["ieee802154"]
    assert f["crc_ok"] == 3 and f["pan_ids"] == ["1234"] and "0013a20040a1b2c3" in f["addresses"]
    assert out["fields"]["model"]["family"] == "Digi XBee (802.15.4 / Zigbee)"
    assert PLUGINS[0].pretty("ieee802154", em) == "IEEE 802.15.4 frames (Digi International (XBee))"
