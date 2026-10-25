"""RadioLink (CC2500), FlySky AFHDS 1 and Hubsan (A7105) (research, PLAN G4 / G5;
research/drone_transmitters/14_market_gap_survey.md section 3): hop tables against the published worked examples,
synthetic packets through the decoders and rc_link detect -> labels, fields, model blocks; nothing from noise."""
import numpy as np

from spectrum_intel.receivers.rc_link import afhds2a as A, detector as D, hubsan as H, radiolink as RL, synth as S
from spectrum_intel.receivers.rc_link.plugin import PLUGINS

FS, FC = 61.44e6, 2440e6


def test_radiolink_tables():
    assert RL.pn9(6) == [0xFF, 0xE1, 0x1D, 0x9A, 0xED, 0x85]                    # TI DN509
    hf, rep = RL.hop_table(bytes.fromhex("3A99223A"))
    assert hf == [0x12, 0xBA, 0x66, 0x72, 0x1E, 0x42, 0x06, 0x4E, 0x36, 0xAE, 0x8A, 0xA2, 0x2A, 0x7E, 0x96, 0x5A]
    assert rep == 198
    assert RL.hop_table(bytes.fromhex("FC110D20"))[0][:4] == [0x84, 0x90, 0xB4, 0xA8]


def test_afhds1_hops():
    assert A.afhds1_hops(0x78)[:4] == [73, 153, 63, 143]
    assert A.afhds1_hops(0x2B)[:4] == [98, 18, 128, 48]


def _y(iq, fs_out):
    from fractions import Fraction
    from scipy.signal import resample_poly
    fr = Fraction(fs_out / S.FS).limit_denominator(1000)
    return resample_poly(iq, fr.numerator, fr.denominator), fr.numerator / fr.denominator * S.FS


def test_decoders_and_noise():
    rng = np.random.default_rng(1)
    fifo = RL.data_packet(0x23, bytes.fromhex("fc110d20"), [1024] * 16, pseudo=26567)
    y, o = _y(np.r_[np.zeros(500), S.radiolink_burst(fifo, rng=rng)[0], np.zeros(500)], 8 * RL.RATE)
    d = RL.try_decode(y, o)
    assert d["sub_protocol"] == "air" and d["telemetry_request"] and d["tx_id"] == "fc110d20" and d["hop_slot"] == 7
    pl = A.afhds1_packet(0x12345678, [1500] * 8)
    y, o = _y(np.r_[np.zeros(500), S.a7105_burst(A.ID_CODE, pl, crc=(0x1021, 0x0000), rng=rng)[0], np.zeros(500)],
              4e6)
    d = A.try_identify(y, o)
    assert d["protocol"] == "FlySky AFHDS 1 (A7105)" and d["tx_id"] == "12345678" and d["crc_ok"]
    pk = H.packet(0x20, bytes([0, 0x80, 0, 0x7F, 0, 0x80, 0, 0x80, 0x02, 0x64]) + bytes.fromhex("db042679"))
    y, o = _y(np.r_[np.zeros(500), S.a7105_burst(0x12ABCDEF, pk, rate=100e3, rng=rng)[0], np.zeros(500)], 8e5)
    d = H.try_decode(y, o)
    assert d["id_code"] == "0x12ABCDEF" and d["models"] == ["H107 / H107L / H107C / H107D"]
    n = rng.normal(size=40000) + 1j * rng.normal(size=40000)
    assert RL.try_decode(n, 8 * RL.RATE) is None and H.try_decode(n, 8e5) is None and A.try_identify(n, 4e6) is None


def test_detect_pipeline():
    rng = np.random.default_rng(3)
    x = S.noise_capture(1 << 21, rng)
    tx = bytes.fromhex("fc110d20")
    hf, _ = RL.hop_table(tx)
    bursts, ps = [], 0
    for k in range(3):                                                  # RadioLink air, 19.9 ms apart
        ps = (ps * 0xAA + 3) % 0x7673
        while abs(RL.F0 + hf[ps & 15] * RL.STEP - FC) > 25e6:           # in the usable band
            ps = (ps * 0xAA + 3) % 0x7673
        iq, _ = S.radiolink_burst(RL.data_packet(0x21, tx, [1024] * 16, ps), cfo_hz=5e3, rng=rng)
        bursts.append(dict(t0_us=500 + 9950 * k, f_hz=RL.F0 + hf[ps & 15] * RL.STEP, iq=iq, bw_hz=150e3))
    pl = A.afhds1_packet(0x12345678, [1500, 1500, 1100, 1500, 1000, 2000, 1500, 1500])
    hops = [h for h in A.afhds1_hops(pl[1]) if abs(A.A7105_F0 + h * 0.5e6 - FC) < 25e6]      # in the usable band
    for k in range(3):                                                  # AFHDS 1, 1.51 ms apart
        iq, _ = S.a7105_burst(A.ID_CODE, pl, crc=(0x1021, 0xFFFF), cfo_hz=8e3, rng=rng)
        bursts.append(dict(t0_us=2000 + 1510 * k, f_hz=A.A7105_F0 + hops[k] * 0.5e6, iq=iq, bw_hz=1.0e6))
    pk = H.packet(0x20, bytes([0, 0x80, 0, 0x7F, 0, 0x80, 0, 0x80, 0x02, 0x64]) + bytes.fromhex("db042679"))
    for k in range(2):                                                  # Hubsan H107 on channel 0x46, 10 ms apart
        iq, _ = S.a7105_burst(0x5A1B2C3D, pk, rate=100e3, cfo_hz=5e3, rng=rng)
        bursts.append(dict(t0_us=12000 + 10000 * k, f_hz=H.F0 + 0x46 * 0.5e6, iq=iq, bw_hz=0.6e6))
    y, _ = S.embed(x, FS, FC, bursts, 14)
    em = D.detect(y, FS, FC)
    by = {lab: [e for e in em if e.label.startswith(lab) and e.decoded] for lab in ("radiolink", "afhds1", "hubsan")}
    assert len(by["radiolink"]) == 3 and all(e.decoded["hop_check"] == "match" for e in by["radiolink"])
    assert len(by["afhds1"]) == 3 and all(e.decoded["hop_check"] == "match" for e in by["afhds1"])
    assert sum(e.decoded["crc_ok"] is True for e in by["afhds1"]) >= 2
    assert len(by["hubsan"]) == 2
    P = PLUGINS[0]
    f = P.summarise(by["radiolink"], None)["fields"]
    assert f["radiolink"]["tx_ids"] == ["fc110d20"] and f["model"]["candidates"] == ["T8FB", "T8S"]
    assert P.pretty("radiolink", by["radiolink"]) == "RadioLink (air)"
    f = P.summarise(by["afhds1"], None)["fields"]
    assert f["afhds1"]["tx_ids"] == ["12345678"] and f["model"]["family"] == "FlySky AFHDS 1"
    f = P.summarise(by["hubsan"], None)["fields"]
    assert f["hubsan"]["id_codes"] == ["0x5A1B2C3D"] and f["model"]["model"] == "H107 / H107L / H107C / H107D"


def test_elrs_gfsk_sync():
    import hashlib
    from spectrum_intel.receivers.rc_link import elrs as E
    rng = np.random.default_rng(4)
    uid = hashlib.md5(b'-DMY_BINDING_PHRASE="myquad"').digest()[:6]
    x = S.noise_capture(1 << 21, rng)
    bursts = []
    for k, (wv, f) in enumerate((("pn9_msb", 2431.4e6), ("pn9_msb", 2447.4e6))):
        bits = np.r_[E.gfsk_frame_bits(uid, E.build_sync(5 + k, 17, 9, 3, uid), wv), np.zeros(8, np.uint8)]
        iq, _ = S.gfsk_burst(300e3, 100e3, bits=bits, bt=2.0, pre_bits=0, cfo_hz=6e3, rng=rng)
        bursts.append(dict(t0_us=3000 + 4000 * k, f_hz=f, iq=iq, bw_hz=0.5e6))
    y, _ = S.embed(x, FS, FC, bursts, 14)
    em = [e for e in D.detect(y, FS, FC) if e.label == "elrs_gfsk"]
    assert len(em) == 2
    f = PLUGINS[0].summarise(em, None)["fields"]
    assert (f["elrs_gfsk"]["uid4"], f["elrs_gfsk"]["uid5"]) == (uid[4], uid[5]) and f["elrs_gfsk"]["whitening"] == ["pn9_msb"]
    assert f["model"]["model"] == "LR1121-based module"
    n = rng.normal(size=20000) + 1j * rng.normal(size=20000)
    assert E.gfsk_try_decode(n, 2.4e6) is None


def _mjxq():
    p = bytearray([0x80, 0x00, 0x80, 0x80, 0x40, 0x40, 0x40, 0xA1, 0xB2, 0x00, 0, 0, 0, 0, 0x04])
    return bytes(p) + bytes([sum(p) & 0xFF])


def _h83d():
    p = bytearray([0x13, 0x11, 0x22, 0x33, 0x44, 1, 0x08, 0x03, (0x11 + 0x22 + 0x33 + 0x44) & 0xFF,
                   0x43, 0x7F, 0x7F, 0x7F, 0x20, 0x20, 0x20, 0x20, 0x01, 0x00])
    return bytes(p) + bytes([sum(p[9:19]) & 0xFF])


def _jjrc():
    p = bytearray(16)
    p[1], p[4:8], p[11], p[12], p[14], p[15] = 0x46, bytes([0x80, 0x80, 0x80, 0x80]), 0x46, 0x02, 0xF1, 0x18
    s = 2 + sum(p[:13])
    p[13], p[3] = s & 0xFF, ((((s >> 8) & 0xFF) << 2) + 2) & 0xFF
    return bytes(p)


def test_toy_xn297_and_nrf24():
    from spectrum_intel.receivers.rc_link import toys as T, xn297 as XN
    rng = np.random.default_rng(5)
    cases = [(bytes.fromhex("6d6a777777"), _mjxq(), "MJXQ (WLH08 / E010 class) (XN297)", "a1b200"),
             (bytes.fromhex("c457096521"), _h83d(), "H8_3D / CG023 / H20 (XN297)", "11223344"),
             (bytes.fromhex("cccccccccc"), _jjrc(), "JJRC345 / SkyTumbler (XN297)", "f118")]
    for addr, pay, name, tid in cases:
        iq, _ = S.xn297_burst(addr, pay, fs=8e6, rng=rng)
        d = XN.try_decode(iq, 8e6)
        assert d and d["protocol"] == name and d["tx_id"] == tid, (name, d)
    v2 = bytearray([0x80, 0x00, 0x80, 0x80, 0x40, 0x40, 0x40, 0x12, 0x34, 0x56, 0, 0, 0, 0, 0x00])
    v2 = bytes(v2) + bytes([sum(v2) & 0xFF])
    sy = bytearray([0x80, 0x10, 0x20, 0x30, 0, 0xC0, 0, 0, 0])
    sy = bytes(sy) + bytes([(T._xor(sy) + 0x55) & 0xFF])
    for addr, pay, pcf, name in ((bytes.fromhex("6688686868"), v2, None, "V2x2 (WLtoys V202 / V222 / V262) (nRF24)"),
                                 (bytes.fromhex("a201020304"), sy, 0x0A1, "SymaX (nRF24)")):
        bits = np.r_[T.nrf_packet_bits(addr, pay, pcf), np.zeros(8, np.uint8)]
        iq, _ = S.gfsk_burst(1e6, 160e3, bits=bits, bt=0.5, pre_bits=0, fs=8e6, rng=rng)
        d = T.nrf24_try_decode(np.r_[np.zeros(80), iq, np.zeros(80)], 8e6)
        assert d and d["protocol"] == name and d["esb_control_field"] == (pcf is not None), (name, d)
    n = rng.normal(size=20000) + 1j * rng.normal(size=20000)
    assert T.nrf24_try_decode(n, 8e6) is None


def test_toy_pipeline():
    from spectrum_intel.receivers.rc_link import toys as T
    rng = np.random.default_rng(6)
    x = S.noise_capture(1 << 21, rng)
    bursts = []
    for k in range(2):
        iq, _ = S.xn297_burst(bytes.fromhex("6d6a777777"), _mjxq(), cfo_hz=10e3, rng=rng)
        bursts.append(dict(t0_us=1000 + 4000 * k, f_hz=2422e6 + 10e6 * k, iq=iq, bw_hz=1.2e6))
    v2 = bytearray([0x80, 0x00, 0x80, 0x80, 0x40, 0x40, 0x40, 0x12, 0x34, 0x56, 0, 0, 0, 0, 0x00])
    v2 = bytes(v2) + bytes([sum(v2) & 0xFF])
    bits = np.r_[T.nrf_packet_bits(bytes.fromhex("6688686868"), v2), np.zeros(8, np.uint8)]
    iq, _ = S.gfsk_burst(1e6, 160e3, bits=bits, bt=0.5, pre_bits=0, cfo_hz=-8e3, rng=rng)
    bursts.append(dict(t0_us=12000, f_hz=2446e6, iq=iq, bw_hz=1.2e6))
    y, _ = S.embed(x, FS, FC, bursts, 15)
    em = [e for e in D.detect(y, FS, FC) if e.label == "toy_drone"]
    assert sorted(e.decoded["protocol"] for e in em) == ["MJXQ (WLH08 / E010 class) (XN297)"] * 2 + \
        ["V2x2 (WLtoys V202 / V222 / V262) (nRF24)"]
    f = PLUGINS[0].summarise([e for e in em if "MJXQ" in e.decoded["protocol"]], None)["fields"]
    assert f["toy_drone"]["tx_ids"] == ["a1b200"] and f["model"]["model"] == "MJXQ (WLH08 / E010 class)"
