"""SiK telemetry radio + MAVLink (research: research/drone_transmitters/05_sik_mavlink_details.md): synthetic SiK packets
carrying HEARTBEAT / GLOBAL_POSITION_INT / VFR_HUD / SYS_STATUS, at 433 MHz, through rc_link detect, from high to low
SNR, for both whitening conventions and both GFSK polarities (the decoder is not told which)."""
import struct
import time

import numpy as np
import pytest

from spectrum_intel.core import mavlink as MAV
from spectrum_intel.receivers.rc_link import detector as D, sik as SIK, synth as S

FS, FC = 61.44e6, 433.0e6
LAT, LON = 12.9716, 77.5946


def frames():
    hb = struct.pack("<IBBBBB", 5, 2, 3, 0x81 | 0x10, 4, 3)                       # LOITER, quadrotor, ArduPilot, armed
    gp = struct.pack("<IiiiihhhH", 123456, int(LAT * 1e7), int(LON * 1e7), 915000, 120000, 150, -30, 0, 9000)
    vf = struct.pack("<ffffhH", 12.5, 11.0, 915.0, 0.4, 90, 47)
    ss = struct.pack("<IIIHHhHHHHHHb", 0, 0, 0, 300, 15800, 1200, 0, 0, 0, 0, 0, 0, 76)
    return (S.mavlink_frame(0, hb, seq=1) + S.mavlink_frame(33, gp, seq=2) + S.mavlink_frame(74, vf, seq=3, v2=True)
            + S.mavlink_frame(1, ss, seq=4))


def test_mavlink_scan():
    fr = MAV.scan(b"\x00\x11" + frames() + b"\x07")
    assert [f["msg"] for f in fr] == ["HEARTBEAT", "GLOBAL_POSITION_INT", "VFR_HUD", "SYS_STATUS"]
    s = MAV.summary(fr)
    assert s["vehicle_type"] == "quadrotor" and s["autopilot"] == "ArduPilot" and s["armed"] and s["flight_mode"] == "LOITER"
    assert abs(s["lat"] - LAT) < 1e-6 and s["alt_rel_m"] == 120.0 and s["groundspeed_mps"] == 11.0 and s["battery_pct"] == 76


@pytest.mark.parametrize("whitening", ["pn9_lsb_first", "pn9_msb_first"])
def test_packet_roundtrip(whitening):
    data = frames() + b"\x12\x34"
    bits = S.sik_packet_bits(0x0019, data, whitening=whitening)
    for inv in (False, True):
        b = bits ^ (1 if inv else 0)
        pos = [p for p, i in SIK._find_sync(b)]
        assert pos, "sync not found"
        r = SIK.decode_packet(b[pos[0]:], inv)
        assert r and r["netid"] == 0x19 and r["whitening"] == whitening and r["mavlink_crc_ok"] == 4
        assert r["hw_crc_ok"] and "ibm_msb_init0000_hi_first" in r["hw_crc_variant"]


def capture(snr_db, seed=4):
    rng = np.random.default_rng(seed)
    n = 1 << 21
    x = ((rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)).astype(np.complex64)
    iq, info = S.sik_burst(0x0019, frames() + b"\x12\x34", rate=64e3, dev_hz=63.75e3, cfo_hz=2.1e3, rng=rng)
    y, _ = S.embed(x, FS, FC, [dict(t0_us=3000.0, f_hz=433.92e6, iq=iq, bw_hz=260e3)], snr_db)
    return y


@pytest.mark.parametrize("snr_db", [20, 12])
def test_sik_in_capture(snr_db):
    y = capture(snr_db)
    t = time.perf_counter()
    em = D.detect(y, FS, FC)
    dt = time.perf_counter() - t
    sik = [e for e in em if e.label == "sik_telemetry"]
    assert sik, f"no SiK decode at {snr_db} dB ({[e.label for e in em]}, {dt:.2f} s)"
    d = sik[0].decoded
    assert d["netid"] == 0x19 and d["mavlink_crc_ok"] == 4
