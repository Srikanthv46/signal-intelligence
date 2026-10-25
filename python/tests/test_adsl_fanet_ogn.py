"""ADS-L SRD-860 M-band, OGN tracker (FLARM PHY) and FANET (LoRa payload) (research, PLAN G10;
research/drone_transmitters/16_uat_adsl_fanet_ogn.md sections B-D): the spec's encoding checks, synthetic packets
through the decoders and rc_link detect -> labels, fields, model blocks; nothing from noise."""
import numpy as np

from spectrum_intel.receivers.rc_link import adsl as A, detector as D, lora_phy as P, ogn as O, synth as S
from spectrum_intel.receivers.rc_link.plugin import PLUGINS

FS, FC = 61.44e6, 868.3e6


def test_adsl_encodings():
    assert A._exp_u(0x0140, 12) - 320 == 0 and A._exp_u(0x0528, 12) - 320 == 1000       # spec G examples
    assert A._exp_u(0xC4, 6) / 4 == 120 and A._exp_s(0x048, 6) / 8 == 10 and A._exp_s(0x148, 6) / 8 == -10
    m = bytes.fromhex("8D4840D6202CC371C32CE0576098")                                    # 1090 ES: same CRC-24
    assert A.crc24(m[:11]) == 0x576098
    # XXTEA key 0, 6 rounds and TEA key 0, 8 loops: values from the reference implementation (oracle)
    assert A.xxtea_encrypt(bytes((i * 7 + 3) & 255 for i in range(20))).hex() == \
        "21f4328c6c00c3cbf949ddf6ce895bf5a078f35e"
    assert O.tea_key0(0x12345678, 0x9ABCDEF0) == (0x5363C8DA, 0x45633369)
    pk = A.packet(0xABC123, 7, A.traffic_payload(32.05, 76.71, 1450, 12))
    d = A.parse(pk)
    assert pk[0] == 0x18 and d["uas"] and d["aircraft_category"] == "UAS specific category" and d["alt_m"] == 1450
    assert A.parse(pk[:-1] + bytes([pk[-1] ^ 1])) is None


def test_fanet_parser():
    # tracking: Skytraxx device 01:1234, UAV, 32.05 N 76.71 E, 1450 m, 36 km/h
    lat, lon = int(round(32.05 * 93206)), int(round(76.71 * 46603))
    w = 1450 | 7 << 12 | 1 << 15
    body = lat.to_bytes(3, "little") + lon.to_bytes(3, "little") + w.to_bytes(2, "little") + bytes([72, 5, 64])
    p = bytes([0x01, 0x01, 0x34, 0x12]) + body
    d = P.identify(0xF1, p)
    assert d["protocol"] == "FANET" and d["uav"] and d["device_id"] == "01:1234" and d["alt_m"] == 1450
    assert d["speed_kmh"] == 36.0 and abs(d["lat"] - 32.05) < 1e-4 and d["manufacturer"] == "Skytraxx"
    assert P.identify(0x12, p)["protocol"] == "FANET"                      # deprecated sync, strict check passes
    assert P.identify(0x12, bytes([0x02, 0x01, 0x34, 0x12]) + b"pilot") is None   # strict: tracking only
    assert P.identify(0xF1, bytes([0x02, 0x07, 0x34, 0x12]) + b"Srikanth")["name"] == "Srikanth"


def test_decoders_and_detect():
    rng = np.random.default_rng(2)
    x = S.noise_capture(1 << 21, rng)
    pk = A.packet(0xABC123, 7, A.traffic_payload(32.05, 76.71, 1450, 12))
    og = O.packet(0xDD1234, 3, 32.05, 76.71, 1450, 13)
    bursts = [dict(t0_us=2000, f_hz=868.2e6, iq=S.adsl_burst(pk, cfo_hz=3e3, rng=rng)[0], bw_hz=0.2e6),
              dict(t0_us=12000, f_hz=868.4e6, iq=S.adsl_burst(pk, cfo_hz=3e3, rng=rng)[0], bw_hz=0.2e6),
              dict(t0_us=20000, f_hz=868.2e6, iq=S.ogn_burst(og, cfo_hz=2e3, rng=rng)[0], bw_hz=0.2e6)]
    y, _ = S.embed(x, FS, FC, bursts, 14)
    em = D.detect(y, FS, FC)
    ad = [e for e in em if e.label == "adsl" and e.decoded]
    on = [e for e in em if e.label == "ogn_like" and e.decoded]
    assert len(ad) == 2 and len(on) == 1, [(e.label, e.center_mhz) for e in em]
    P_ = PLUGINS[0]
    f = P_.summarise(ad, None)["fields"]
    assert f["adsl"]["uas"] and f["adsl"]["addresses"] == ["OGN tracker:ABC123"]
    assert P_.pretty("adsl", ad) == "ADS-L iConspicuity (SRD-860), UAS"
    f = P_.summarise(on, None)["fields"]
    assert f["ogn"]["uav"] and f["ogn"]["addresses"] == ["OGN:DD1234"]
    n = rng.normal(size=40000) + 1j * rng.normal(size=40000)
    assert A.try_decode(n, 8e5) is None and O.try_decode(n, 8e5) is None
