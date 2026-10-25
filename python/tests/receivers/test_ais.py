"""AIS (research, signal gaps S15): golden vector (research file 06), synthetic GMSK bursts, pipeline."""
import numpy as np

from spectrum_intel.core.dsp import Channelizer
from spectrum_intel.receivers.marine import ais as A
from spectrum_intel.receivers.marine import detect
from spectrum_intel.receivers.marine import synth as S

P = "177KQJ5000G?tO`K>RA1wUbN0TKH"                     # gpsd AIVDM example
PRE = ("010101010101010101010101" "01111110"
       "001000001000111011011011101000011000010100000010000000001010000011110011100011110101111101101100"
       "001011100001001001000001001111101101011010011110010100000001100010000110110010001000010101" "01111110")
NRZI = ("00110011001100110011001100000001011010100101111000111000011010111010110010101001010101011001010111110111"
        "10100000110000001110001010011110101101101101010010000001110011100100000100110101010001011010111000100101"
        "101011001100000001")


def _sixbit(text, n):
    v = [ord(c) - 64 if 64 <= ord(c) < 96 else ord(c) for c in text.ljust(n, "@")[:n]]
    return [(x >> (5 - i)) & 1 for x in v for i in range(6)]


def _u(v, n):
    v &= (1 << n) - 1
    return [(v >> (n - 1 - i)) & 1 for i in range(n)]


def test_golden_vector():
    bits = A.armor_to_bits(P)
    o = A.bits_to_octets(bits)
    assert o.hex() == "0471db85a1400005cff1fa1b3a2441fe5a9e0246d8" and A.crc_x25(o) == 0xA844
    fb = A.frame_bits(bits)
    assert "".join(map(str, fb)) == PRE and "".join(map(str, A.nrzi(fb))) == NRZI
    m = A.fields(A.frames([0, 1, 0, 1, 0] + A.unnrzi([1] + A.nrzi(fb)))[0][1])
    assert m["mmsi"] == "477553000" and m["nav_status"] == "moored" and m["lon"] == -122.345833
    assert m["lat"] == 47.582833 and m["cog_deg"] == 51.0 and m["heading_deg"] == 181 and m["utc_second"] == 15


def test_message5_and_18_fields():
    b5 = (_u(5, 6) + _u(0, 2) + _u(419012345, 30) + _u(0, 2) + _u(9876543, 30) + _sixbit("VTAB1", 7)
          + _sixbit("SAGAR KANYA", 20) + _u(70, 8) + _u(80, 9) + _u(20, 9) + _u(6, 6) + _u(6, 6) + _u(1, 4)
          + _u(10, 4) + _u(25, 5) + _u(12, 5) + _u(30, 6) + _u(65, 8) + _sixbit("MUMBAI", 20) + [0, 0])
    m = A.fields(A.bits_to_octets(b5), len(b5))
    assert m["flag"] == "India" and m["name"] == "SAGAR KANYA" and m["callsign"] == "VTAB1"
    assert m["destination"] == "MUMBAI" and m["length_m"] == 100 and m["beam_m"] == 12 and m["draught_m"] == 6.5
    b18 = (_u(18, 6) + _u(0, 2) + _u(419000001, 30) + _u(0, 8) + _u(52, 10) + [1] + _u(int(72.85 * 600000), 28)
           + _u(int(18.95 * 600000), 27) + _u(900, 12) + _u(90, 9) + _u(0, 35))
    m = A.fields(A.bits_to_octets(b18), len(b18))
    assert m["kind"] == "class B position report" and m["sog_kn"] == 5.2 and m["heading_deg"] == 90
    assert abs(m["lon"] - 72.85) < 1e-5 and abs(m["lat"] - 18.95) < 1e-5


def test_bursts_both_channels_with_offsets():
    fs, fc = 1e6, 162.0e6
    for cfo in (0.0, 1500.0, -1500.0):
        x = S.capture(fs, 0.1, fc, [(0.005, 161.975e6 + cfo, P), (0.04, 162.025e6 + cfo, P)], snr_db=15, seed=1)
        em = detect(Channelizer(x, fs, fc))
        assert [e.features["channel"] for e in em] == ["AIS 1 (87B)", "AIS 2 (88B)"]
        assert all(e.decoded["mmsi"] == "477553000" for e in em)


def test_noise_and_pipeline():
    fs, fc = 1e6, 162.0e6
    rng = np.random.default_rng(2)
    n = 1 << 17
    x = ((rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)).astype(np.complex64)
    assert detect(Channelizer(x, fs, fc)) == []
    from spectrum_intel.analyze import analyze_iq
    from spectrum_intel.decode.registry import plugins
    y = S.capture(fs, 0.1, fc, [(0.005, 161.975e6, P), (0.04, 162.025e6, P)], snr_db=20, seed=3)
    a = analyze_iq(y, fs, fc, plugins=plugins(enable=["marine"]))
    s = [s for s in a.signals if s.features.get("decoder") == "marine"]
    assert [x.label for x in s] == ["AIS (1 station)"] and s[0].decoded["fields"]["n_messages"] == 2
