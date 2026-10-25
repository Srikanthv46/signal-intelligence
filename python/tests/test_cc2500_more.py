"""Graupner HoTT (CC2500 MSK + FEC), Hitec Optima / Minima (MSK, sync = TX ID) and Futaba S-FHSS (2-FSK) (research,
PLAN G12; research/drone_transmitters/17_hott_hitec_r9.md): the CC2500 FEC against the TI DN504 test vector,
synthetic packets through the decoders and rc_link detect -> labels, fields, model blocks; nothing from noise."""
import numpy as np

from spectrum_intel.receivers.rc_link import cc2500 as CC, detector as D, hitec as HI, hott as H, synth as S
from spectrum_intel.receivers.rc_link.plugin import PLUGINS

FS, FC = 61.44e6, 2430e6


def test_fec_dn504_vector():
    d = bytes.fromhex("03010203")
    c = CC.crc16_cc2500(d)
    assert c == 0x303A
    b = H.fec_encode(d + bytes([c >> 8, c & 0xFF]))
    assert CC._bytes(b).hex() == "c83c002084cf3331a2fc404a443047ef"          # DN504 interleaver output
    b[[3, 40, 90]] ^= 1
    assert H.fec_decode(b, 6) == (bytes.fromhex("03010203303a"), 3)


def _y(iq, fs_out):
    from fractions import Fraction
    from scipy.signal import resample_poly
    fr = Fraction(fs_out / S.FS).limit_denominator(1000)
    return resample_poly(iq, fr.numerator, fr.denominator), fr.numerator / fr.denominator * S.FS


def _hott_pkt():
    return H.tx_packet(3, 17, [3000] * 12, bytes.fromhex("7c94000d50"), bytes.fromhex("31bd010090"))


def test_decoders_and_noise():
    rng = np.random.default_rng(1)
    y, o = _y(np.r_[np.zeros(500), S.hott_burst(_hott_pkt(), rng=rng)[0], np.zeros(500)], 8 * H.RATE)
    d = H.try_decode(y, o)
    assert d["tx_id"] == "7c94000d50" and d["rx_id"] == "31bd010090" and d["hop_table"] == 3 and d["hop_index"] == 17
    assert d["channels_us"][0] == 1500.0
    y, o = _y(np.r_[np.zeros(500), S.hitec_burst(HI.hitec_frame(bytes.fromhex("00036a"), [10822] * 9, 1), 0x6A03,
                                                 rng=rng)[0], np.zeros(500)], 8 * HI.RATE)
    d = HI.hitec_try_decode(y, o)
    assert d["tx_id"] == "00036a" and d["sync_matches_tx_id"] and d["repeat_flag"] == "0x40"
    # a CRC-valid frame whose sync is not its TX ID is rejected
    y, o = _y(np.r_[np.zeros(500), S.hitec_burst(HI.hitec_frame(bytes.fromhex("00036a"), [10822] * 9), 0x1234,
                                                 rng=rng)[0], np.zeros(500)], 8 * HI.RATE)
    assert HI.hitec_try_decode(y, o) is None
    pk = HI.sfhss_packet(0x1234, 7, [1520, 1020, 2020, 1520], 13, model=bytes.fromhex("4207"))
    y, o = _y(np.r_[np.zeros(500), S.sfhss_burst(pk, 0x0000, rng=rng)[0], np.zeros(500)], 8 * HI.SF_RATE)
    d = HI.sfhss_try_decode(y, o, f_hz=HI.SF_F0 + 7 * HI.SF_STEP)
    assert d["tx_id"] == "1234" and d["model_hint"] == "T8J" and d["hop_check"] == "match" and d["crc_init"] == "0x0000"
    assert d["channels"] == [1520, 1020, 2020, 1520] and d["hop_code"] == 13
    n = rng.normal(size=40000) + 1j * rng.normal(size=40000)
    assert H.try_decode(n, 8 * H.RATE) is None and HI.hitec_try_decode(n, 8 * H.RATE) is None
    assert HI.sfhss_try_decode(n, 8 * HI.SF_RATE) is None


def test_detect_pipeline():
    rng = np.random.default_rng(3)
    x = S.noise_capture(1 << 21, rng)
    bursts = []
    for k, n in enumerate((20, 25)):                                    # HoTT TX packets, 10 ms apart
        iq, _ = S.hott_burst(_hott_pkt(), cfo_hz=5e3, rng=rng)
        bursts.append(dict(t0_us=500 + 10000 * k, f_hz=H.F0 + n * H.STEP, iq=iq, bw_hz=0.5e6))
    for k in range(2):                                                  # Hitec, 4 repeats 2.736 ms apart
        iq, _ = S.hitec_burst(HI.hitec_frame(bytes.fromhex("00036a"), [10822] * 9, k), 0x6A03, cfo_hz=4e3, rng=rng)
        bursts.append(dict(t0_us=4000 + 2736 * k, f_hz=HI.F0 + 58 * HI.STEP, iq=iq, bw_hz=0.5e6))
    pk = HI.sfhss_packet(0x1234, 9, [1520, 1020, 2020, 1520], 13, model=bytes.fromhex("0f09"))
    for k in range(2):                                                  # S-FHSS DATA1 / DATA2, 1.625 ms apart
        iq, _ = S.sfhss_burst(pk, 0xFFFF, cfo_hz=3e3, rng=rng)
        bursts.append(dict(t0_us=15000 + 1625 * k, f_hz=HI.SF_F0 + 9 * HI.SF_STEP, iq=iq, bw_hz=0.3e6))
    y, _ = S.embed(x, FS, FC, bursts, 14)
    em = D.detect(y, FS, FC)
    by = {lab: [e for e in em if e.label == lab and e.decoded] for lab in ("hott", "hitec", "sfhss")}
    assert len(by["hott"]) == 2 and len(by["hitec"]) == 2 and len(by["sfhss"]) == 2, \
        {k: len(v) for k, v in by.items()}
    P = PLUGINS[0]
    f = P.summarise(by["hott"], None)["fields"]
    assert f["cc2500_rc"]["tx_ids"] == ["7c94000d50"] and f["cc2500_rc"]["hop_table"] == 3
    assert f["model"]["family"] == "Graupner HoTT"
    assert P.pretty("hott", by["hott"]) == "Graupner HoTT"
    f = P.summarise(by["hitec"], None)["fields"]
    assert f["cc2500_rc"]["tx_ids"] == ["00036a"] and f["model"]["family"] == "Hitec"
    f = P.summarise(by["sfhss"], None)["fields"]
    assert f["cc2500_rc"]["hop_checks"] == ["match"] and f["model"]["model"] == "T10J"
