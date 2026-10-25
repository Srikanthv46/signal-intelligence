"""PLAN W1 (receivers/wifi/field.py): drone-link verdicts from identity and from behaviour; W6 / W7 (receivers/devices):
beacon floods, randomised MACs linked by probe signature + sequence counter, the behaviour finding. Frames are real mac80211.parse output of
hand-built MPDUs; the end-to-end case runs the OFDM plugin on a synthetic capture (ofdm_tx)."""
import struct
from types import SimpleNamespace

import numpy as np

from spectrum_intel.receivers.wifi import field as FD
from spectrum_intel.receivers.wifi import mac80211 as M
from spectrum_intel.receivers.wifi import ofdm_rx as W
from spectrum_intel.receivers.wifi import ofdm_tx as T
from spectrum_intel.receivers.wifi import plugin_ofdm as PO

AP, PHONE = bytes.fromhex("00e04c0a0b0c"), bytes.fromhex("3c2ef9445566")
HOME, S1, S2, S3 = (bytes.fromhex(h) for h in ("ac84c6112233", "a45e60000001", "a45e60000002", "a45e60000003"))


def mgmt(sub, ta, ra, body, bssid=None, seq=0):
    return bytes([sub << 4, 0, 0, 0]) + ra + ta + (bssid or ta) + struct.pack("<H", seq << 4) + body


def beacon(ta, ssid, rsn=True):
    body = bytes(8) + b"\x64\x00" + (b"\x11\x04" if rsn else b"\x01\x04") + bytes([0, len(ssid)]) + ssid \
        + b"\x01\x04\x82\x84\x8b\x96"
    if rsn:
        body += bytes([48, 20, 1, 0]) + b"\x00\x0f\xac\x04" + b"\x01\x00\x00\x0f\xac\x04" + b"\x01\x00\x00\x0f\xac\x02" + b"\x00\x00"
    return mgmt(8, ta, b"\xff" * 6, body)


def data(ta, ra, bssid, to_ds, n, seq=0):
    fc1 = (0x01 if to_ds else 0x02) | 0x40                       # protected (CCMP): payload opaque
    a3 = ra if to_ds else ta
    return bytes([0x88, fc1, 0, 0]) + (bssid if to_ds else ra) + ta + a3 + struct.pack("<H", seq << 4) \
        + b"\x00\x00" + bytes(8) + bytes(n)


def ack(ra):
    return b"\xd4\x00\x00\x00" + ra


def fr(mpdu, t0, dur_us, ta=None, nbytes=None):
    m = M.parse(mpdu)
    ms = [m]
    return {"t0": t0, "t1": t0 + dur_us * 1e-6, "dur_us": dur_us, "bytes": nbytes or len(mpdu) + 4,
            "ta": ta or m.get("ta"), "macs": ms, "ctrl": m.get("type") == "ctrl", "video": False}


def drone_block(n=40):
    """A drone AP streaming to one phone: 1400-byte downlink every 0.8 ms, a 60-byte uplink every 4 ms."""
    out = [fr(beacon(AP, b"Drone-Cam-x", rsn=True), 0.0, 300)]
    for k in range(n):
        out.append(fr(data(AP, PHONE, AP, False, 1400, seq=k), 1e-3 + k * 0.8e-3, 230))
        out.append(fr(ack(AP), 1e-3 + k * 0.8e-3 + 250e-6, 30, ta=PHONE.hex(":")))
        if k % 5 == 0:
            out.append(fr(data(PHONE, AP, AP, True, 60, seq=k), 1e-3 + k * 0.8e-3 + 400e-6, 40))
    for f in out:                                                 # attributed transmitter as links.attribute gives it
        f["ta"] = f["ta"] if isinstance(f["ta"], str) else f["ta"]
    return out


def home_block():
    """A busy home AP: three clients, traffic both ways."""
    out = [fr(beacon(HOME, b"MyHome"), 0.0, 300)]
    for k in range(30):
        c = (S1, S2, S3)[k % 3]
        out.append(fr(data(HOME, c, HOME, False, 900, seq=k), 1e-3 + k * 1e-3, 200))
        out.append(fr(data(c, HOME, HOME, True, 1200, seq=k), 1.4e-3 + k * 1e-3, 250))
    return out


def test_behaviour_verdict_drone_vs_home():
    lk = {r["bssid"]: r for r in FD.links(drone_block() + home_block(), 0.034)}
    d, h = lk[AP.hex(":")], lk[HOME.hex(":")]
    assert d["clients"] == [PHONE.hex(":")] and d["downlink"]["big_ppdus"] == 40
    v = d["drone_link"]
    assert v["basis"] == "behaviour" and v["confidence"] == 0.55 and v["verdict"].endswith("-like")
    assert h["n_clients"] == 3 and "drone_link" not in h


def test_identity_verdict_tello():
    sdk = bytes([0x45, 0, 0, 35, 0, 0, 0, 0, 64, 17, 0, 0, 192, 168, 10, 2, 192, 168, 10, 1]) \
        + struct.pack(">HHHH", 50000, 8889, 15, 0) + b"takeoff"
    cmd = bytes([0x08, 0x01, 0, 0]) + AP + PHONE + AP + b"\x10\x00" + b"\xaa\xaa\x03\x00\x00\x00\x08\x00" + sdk
    frames = drone_block(8) + [fr(cmd, 0.03, 60)]
    v = {r["bssid"]: r for r in FD.links(frames, 0.034)}[AP.hex(":")]["drone_link"]
    assert v["basis"] == "identity" and "Tello" in v["platform"]


def test_beacon_flood_and_sequence_linking():
    """W6 / W7 in receivers/devices: forged beacons from unrelated unregistered prefixes; two randomised MACs with
    one probe signature whose sequence counter continues across the change."""
    from spectrum_intel.receivers.devices import threats as TH
    from spectrum_intel.receivers.devices import tracker as TR
    rng = np.random.default_rng(1)
    flood = [TR.frame(M.parse(beacon(bytes([0x02 | (int(rng.integers(0, 64)) << 2)] + list(rng.integers(0, 256, 5))),
                                     b"free-%d" % k)), t0=k * 1e-3) for k in range(16)]
    home = [TR.frame(M.parse(beacon(bytes([0x02, 0x84, 0xC6, 0x11, 0x22, k]), b"Home-%d" % k)), t0=k * 1e-3)
            for k in range(16)]                                    # 16 virtual BSSIDs of one radio: no flood
    assert [f["kind"] for f in TH.findings(TR.tracks(flood), flood) if f["kind"] == "beacon_flood"] == ["beacon_flood"]
    assert not [f for f in TH.findings(TR.tracks(home), home) if f["kind"] == "beacon_flood"]
    ies = b"\x00\x00" + b"\x01\x08\x02\x04\x0b\x16\x0c\x12\x18\x24" + b"\x7f\x08\x04\x00\x00\x00\x00\x00\x00\x40"
    r1, r2 = bytes([0xDA, 0xA1, 0x19, 0, 0, 1]), bytes([0xDA, 0xA1, 0x19, 0, 0, 2])
    fr_ = [TR.frame(M.parse(mgmt(4, r1, b"\xff" * 6, ies, bssid=b"\xff" * 6, seq=100 + k)), t0=k * 1e-3) for k in range(3)]
    fr_ += [TR.frame(M.parse(mgmt(4, r2, b"\xff" * 6, ies, bssid=b"\xff" * 6, seq=104 + k)), t0=0.01 + k * 1e-3)
            for k in range(3)]
    g = TR.link_randomized(TR.tracks(fr_))
    assert len(g) == 1 and g[0]["sequence_continues"] and g[0]["confidence"] == 0.9
    fr_[3:] = [TR.frame(M.parse(mgmt(4, r2, b"\xff" * 6, ies, bssid=b"\xff" * 6, seq=2000 + k)), t0=0.01 + k * 1e-3)
               for k in range(3)]
    g = TR.link_randomized(TR.tracks(fr_))
    assert len(g) == 1 and not g[0]["sequence_continues"] and g[0]["confidence"] < 0.9


def test_drone_behaviour_finding():
    from spectrum_intel.receivers.devices import threats as TH
    lk = FD.links(drone_block(), 0.034)
    dl = [{"bssid": r["bssid"], "ssid": r.get("ssid"), "clients": r["clients"], **r["drone_link"]}
          for r in lk if r.get("drone_link")]
    f = [x for x in TH.findings([], [], [], dl) if x["kind"] == "drone_link_behaviour"]
    assert len(f) == 1 and f[0]["bssid"] == AP.hex(":") and PHONE.hex(":") in f[0]["devices"]


def test_plugin_end_to_end_synthetic():
    """The OFDM plugin on a synthetic capture: an encrypted drone-style link (24 Mb/s non-HT) gets a behaviour
    verdict in its signal fields and label."""
    fs, fc = 30.72e6, 2437e6
    pp, t = [], 60e-6
    pp.append(dict(samples=T.ppdu(beacon(AP, b"Drone-Cam-x"), fmt="NON_HT", rate_mbps=6, scramble_seed=0x11),
                   f_hz=fc, t_s=t, snr_db=25))
    t += 600e-6
    for k in range(14):
        pp.append(dict(samples=T.ppdu(data(AP, PHONE, AP, False, 1400, seq=k), fmt="NON_HT", rate_mbps=24,
                                      scramble_seed=1 + k), f_hz=fc, t_s=t, snr_db=24))
        t += 560e-6
        pp.append(dict(samples=T.ppdu(ack(AP), fmt="NON_HT", rate_mbps=24, scramble_seed=40 + k), f_hz=fc, t_s=t,
                       snr_db=20))
        t += 80e-6
    n = int((t + 200e-6) * fs)
    x = T.scene(pp, n, fs=fs, fc=fc, seed=4)
    recs = PO._records(W.WifiDetector(fs, fc).run(x), SimpleNamespace(fs=fs))
    res = PO.WifiOfdmDecoder().summarise(recs, len(recs), 0, recs, n / fs)
    f = res["fields"]
    assert f["n_fcs_ok"] >= 14
    dl = f.get("drone_links")
    assert dl and dl[0]["bssid"] == AP.hex(":") and dl[0]["basis"] == "behaviour"
    assert "drone link-like" in res["label"]
    assert np.isclose(f["wifi_links"][0]["downlink"]["median_bytes"], 1400 + 36, atol=8)


def test_wifi_drone_models_by_ssid_and_oui():
    """The Wi-Fi drones of the 2026-10-25 list (research file 18): default network names and registered OUIs."""
    from spectrum_intel.receivers.wifi import droneid as D
    want = {("Spark-1A2B3C", "60:60:1f:00:00:01"): ("DJI Spark", 0.9),
            ("MavicAir-1A2B3C", "60:60:1f:00:00:01"): ("DJI Mavic Air", 0.9),
            ("Mavic-Mini-1a2b3c", None): ("DJI Mavic Mini", 0.7),
            ("PHANTOM3_09e35b", "60:60:1f:00:00:01"): ("DJI Phantom 3 Standard / 4K (Wi-Fi video models)", 0.9),
            ("CGO3_a1b2c3", "e0:b6:f5:80:00:01"): ("Yuneec CGO3 camera video link (Typhoon Q500 4K / Typhoon H)", 0.9),
            ("ANAFI-1234567", "90:3a:e6:00:00:01"): ("Parrot Anafi", 0.9),
            ("RUKO-F11-12345678-5G", None): ("Ruko F11", 0.7),
            ("HolyStoneFPV_12345", None): ("Holy Stone (Wi-Fi FPV: HS720 / HS100 series)", 0.7),
            ("TELLO-ABCDEF", "60:60:1f:00:00:01"): ("DJI Tello", 0.9),
            ("anything", "38:1d:14:00:00:01"): ("Skydio", 0.75),
            ("anything", "98:aa:fc:70:00:01"): ("Hubsan", 0.75)}
    for (ssid, oui), (plat, conf) in want.items():
        s = D.signature(M.vendor(oui) if oui else None, ssid)
        assert s and s["platform"] == plat and s["confidence"] == conf, (ssid, s)
    assert D.signature(None, "MyHome-5G") is None
