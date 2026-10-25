"""UAT 978 MHz ADS-B and China GB 46750-2025 Remote ID (research, PLAN G14;
research/drone_transmitters/16_uat_adsl_fanet_ogn.md sections A and E): the clean-room Reed-Solomon codec against
reference parity, dump978's sample messages, a synthetic UAV message through demodulation and the pipeline plugin, and
the published DJI Mini 5 Pro GB 46750 frame."""
import random

import numpy as np

from spectrum_intel.core import remote_id as R
from spectrum_intel.core.rs import RS
from spectrum_intel.receivers.adsb import uat as U

GB_MINI5 = ("ff2048fffffe3135383146414e4c433235385530323952544e363030303030303030000101d1823b483bf2eb11ed0769833b4822"
            "f0eb1105001c00284700c2083d0902000c050478acc5529e0103")


def test_rs_against_reference_parity():
    # parity from a reference RS implementation (oracle, 8 / 0x187 / fcr 120 / prim 1), random data seed 7
    ref = bytes.fromhex("9197f4299179e8f6d9212a56d92ee5f4a577bfb47d9f1d11a3897b86c0f2")
    assert RS(0x187, 120, 1, 12).encode(ref[:18]) == ref[18:]
    rng = random.Random(3)
    rs = RS(0x187, 120, 1, 14)
    for _ in range(50):
        data = bytes(rng.randrange(256) for _ in range(34))
        cw = bytearray(data + rs.encode(data))
        ne = rng.randrange(8)
        for p in rng.sample(range(48), ne):
            cw[p] ^= rng.randrange(1, 256)
        assert rs.decode(bytes(cw)) == (data + rs.encode(data), ne)
    cw = bytearray(data + rs.encode(data))
    for p in range(9):                                           # beyond capacity: rejected
        cw[p] ^= 0x55
    assert rs.decode(bytes(cw))[0] is None


def test_uat_messages():
    m = U.parse(bytes.fromhex("08a66ef1353e2d525fd4050911882aa038101d06b85d440be2a4c2a0000590000000"))  # dump978 sample
    assert m["address"] == "A66EF1" and m["callsign"] == "N5130E" and m["emitter_category"] == "small"
    assert abs(m["lat"] - 37.43639) < 1e-4 and abs(m["lon"] + 122.08055) < 1e-4 and m["alt_ft"] == 975
    p = U.payload(0xA1B2C3, 37.5, -122.25, 1200, 14, "N123UA")
    rng = np.random.default_rng(1)
    s = U.modulate(p)
    y = np.r_[np.zeros(1000), s * np.exp(2j * np.pi * 8e3 * np.arange(len(s)) / U.FS), np.zeros(1000)]
    y = y + (rng.normal(size=len(y)) + 1j * rng.normal(size=len(y))) * 0.3
    r = U.demod(y)
    assert len(r) == 1 and r[0]["uav"] and r[0]["callsign"] == "N123UA" and r[0]["address"] == "A1B2C3"
    n = (rng.normal(size=200000) + 1j * rng.normal(size=200000))
    assert U.demod(n) == []


def test_uat_pipeline():
    from spectrum_intel.analyze import analyze
    from spectrum_intel.decode import registry
    from spectrum_intel.io.sigmf import Capture
    from scipy.signal import resample_poly
    fs, fc = 15.36e6, 975e6
    rng = np.random.default_rng(2)
    s = U.modulate(U.payload(0xA1B2C3, 37.5, -122.25, 1200, 14, "N123UA"), fs=fs)
    n = 1 << 18
    x = (rng.normal(size=n) + 1j * rng.normal(size=n)) * 0.05
    t = np.arange(len(s)) / fs
    x[50000:50000 + len(s)] += s * np.exp(2j * np.pi * (978e6 - fc) * t)
    a = analyze(Capture.from_iq(x.astype(np.complex64), fs, fc, name="uat"), None, registry.plugins(["adsb"]))
    sig = [s_ for s_ in a.signals if s_.label.startswith("UAT")]
    assert sig, [(s_.label, s_.decoded) for s_ in a.signals]
    f = sig[0].decoded["fields"]
    assert f["aircraft"][0]["uav"] and f["aircraft"][0]["callsign"] == "N123UA"


def test_gb46750_mini5():
    g = R.parse_gb46750(bytes.fromhex(GB_MINI5))
    assert g["uas_id"] == "1581FANLC258U029RTN6" and "registration_mark" not in g
    assert (g["lat"], g["lon"]) == (30.0675106, 121.1859817) and (g["operator_lat"], g["operator_lon"]) == (30.0675643, 121.1859665)
    assert g["track_deg"] == 0.5 and g["speed_mps"] == 2.8 and g["height_m"] == 108.0 and g["geodetic_altitude_m"] == 121.0
    assert g["pressure_altitude_m"] == 182.5 and g["operator_altitude_m"] == 14.5 and g["status"] == "airborne"
    assert g["ua_class"] == "light" and g["manufacturer"]["maker"] == "DJI"
    # carried in the F3411 Wi-Fi beacon element (FA:0B:BC / 0x0D) with a counter byte
    w = R.parse_odid(bytes([7]) + bytes.fromhex(GB_MINI5))
    assert w["standard"] == "China GB 46750-2025" and w["message_counter"] == 7
    bad = bytearray.fromhex(GB_MINI5)
    bad[2] = 0x47                                                # length disagrees with the flags
    assert R.parse_gb46750(bytes(bad)) is None
