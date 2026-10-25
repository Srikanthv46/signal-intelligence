#!/usr/bin/env python3
"""
test_det_devices.py -- device tracking and findings over decoded 802.11 frames
(receivers/devices/tracker.py, threats.py).

  1. Access points and stations are told apart, and each carries its vendor, activity and signal.
  2. Randomised MACs that share a probe-request signature and a transmitter clock are linked.
  3. Findings: deauthentication floods, duplicate SSIDs, weak security, drones, key exchanges.
  4. Rules that would cry wolf do not: virtual BSSIDs on one radio, TKIP offered beside CCMP,
     and a single roaming deauthentication are not raised as attacks.
  5. The plugin reads the generic report contract, so it works for 802.11b and OFDM alike.

pytest collects the test_* functions; `python test_det_devices.py` prints the findings it produces.
"""
from __future__ import annotations

from spectrum_intel.receivers.devices import plugin as PL
from spectrum_intel.receivers.devices import threats as TH
from spectrum_intel.receivers.devices import tracker as T
from spectrum_intel.receivers.wifi import dsss_tx as D
from spectrum_intel.receivers.wifi import mac80211 as M

AP = "24:a3:f0:a5:bd:ab"
STA = "f8:e9:4e:11:22:33"


def _beacon(bssid=AP, ssid="Lab AP", ch=6, **kw):
    return M.parse(D.beacon(bssid, ssid, ch, **kw))


def _mgmt(subtype, ta, ra, bssid, **extra):
    """A small management frame, built as bytes so the dissector does the work."""
    sub = {"deauth": 0xC0, "disassoc": 0xA0, "probe-req": 0x40}[subtype]
    a = (lambda m: bytes(int(x, 16) for x in m.split(":")))
    body = b"\x07\x00" if subtype in ("deauth", "disassoc") else b"\x00\x00"
    if subtype == "probe-req":
        s = extra.get("ssid", "").encode()
        body = bytes([0, len(s)]) + s + bytes([1, 4, 0x82, 0x84, 0x8B, 0x96]) + bytes([50, 1, 0x0C])
    return M.parse(bytes([sub, 0, 0, 0]) + a(ra) + a(ta) + a(bssid) + b"\x00\x00" + body)


def _frames(macs, dt=1e-3, snr=20.0, cfo=1200.0):
    return [T.frame(m, t0=i * dt, dur_s=200e-6, snr_db=snr, cfo_hz=cfo, channel=6,
                    phy="802.11b", source="wifi_dsss")
            for i, m in enumerate(macs)]


def test_access_points_and_stations_are_told_apart():
    fr = _frames([_beacon(), _beacon(seq=2), _mgmt("probe-req", STA, "ff:ff:ff:ff:ff:ff",
                                                   "ff:ff:ff:ff:ff:ff", ssid="HomeNet")])
    devs = T.tracks(fr)
    by = {d["mac"]: d for d in devs}
    assert by[AP]["role"] == "access point" and by[AP]["ssid"] == "Lab AP"
    assert by[AP]["n_frames"] == 2 and by[AP]["security"].startswith("WPA2-Personal")
    assert by[AP]["vendor"].startswith("TeraNXT") and by[AP]["address"] == "global (OUI)"
    assert by[STA]["role"] == "station" and by[STA]["probed_ssids"] == ["HomeNet"]
    assert by[AP]["airtime_ms"] > 0 and by[AP]["snr_db_median"] == 20.0
    assert by[AP]["channels"] == {6: 2}


def test_randomized_macs_are_linked_by_fingerprint_and_clock():
    fr = []
    for i, mac in enumerate(("0e:11:22:33:44:55", "0a:aa:bb:cc:dd:ee", "06:99:88:77:66:55")):
        cfo = 1500.0 if i < 2 else 40000.0            # the third has a different crystal
        fr += [T.frame(_mgmt("probe-req", mac, "ff:ff:ff:ff:ff:ff", "ff:ff:ff:ff:ff:ff", ssid="Cafe"),
                       t0=i * 1e-3, dur_s=1e-4, snr_db=15, cfo_hz=cfo, channel=6, source="wifi_dsss")]
    devs = T.tracks(fr)
    assert all(d["address"].startswith("randomized") for d in devs)
    groups = T.link_randomized(devs)
    assert len(groups) == 1 and groups[0]["n_addresses"] == 2
    assert set(groups[0]["macs"]) == {"0e:11:22:33:44:55", "0a:aa:bb:cc:dd:ee"}
    assert 0 < groups[0]["confidence"] < 1


def test_deauth_flood_is_flagged_and_a_single_one_is_not():
    one = _frames([_mgmt("deauth", AP, STA, AP)])
    assert not [f for f in TH.findings(T.tracks(one), one) if f["kind"] == "deauth_flood"]
    many = _frames([_mgmt("deauth", AP, "ff:ff:ff:ff:ff:ff", AP) for _ in range(30)])
    hit = [f for f in TH.findings(T.tracks(many), many) if f["kind"] == "deauth_flood"]
    assert hit and hit[0]["severity"] == "critical" and hit[0]["n_frames"] == 30
    assert AP in hit[0]["devices"]


def test_duplicate_ssid_versus_virtual_bssids_on_one_radio():
    # Same radio, two virtual BSSs: the locally-administered bit and one byte differ -> not an attack.
    same = _frames([_beacon("f4:f6:47:43:7f:e6", "HomeNet"), _beacon("f6:f6:47:53:7f:e6", "HomeNet")])
    kinds = {f["kind"] for f in TH.findings(T.tracks(same), same)}
    assert "multi_bssid" in kinds and "duplicate_ssid" not in kinds
    # Different hardware, and the copy is open: that is worth raising.
    twin = _frames([_beacon("24:a3:f0:a5:bd:ab", "HomeNet"),
                    _beacon("00:15:6d:01:02:03", "HomeNet", privacy=False)])
    hit = [f for f in TH.findings(T.tracks(twin), twin) if f["kind"] == "duplicate_ssid"]
    assert hit and hit[0]["severity"] == "high" and hit[0]["ssid"] == "HomeNet"


def test_security_findings_are_proportionate():
    openf = _frames([_beacon("24:a3:f0:a5:bd:ab", "FreeWiFi", privacy=False)])
    got = {f["kind"]: f for f in TH.findings(T.tracks(openf), openf)}
    assert got["open_network"]["severity"] == "low"
    devs = T.tracks(openf)
    devs[0]["security"] = "WPA2-Personal (PSK) TKIP/CCMP-128"      # mixed mode: common, not broken
    k = {f["kind"]: f for f in TH.findings(devs, openf)}
    assert k["legacy_cipher_offered"]["severity"] == "low"
    devs[0]["security"] = "WPA1 TKIP"                              # no modern cipher at all
    k = {f["kind"]: f for f in TH.findings(devs, openf)}
    assert k["weak_encryption"]["severity"] == "medium"


def test_drone_finding_carries_the_operator_position():
    import sys, os
    sys.path.insert(0, os.path.dirname(__file__))
    from test_det_droneid import odid_ie_body
    bc = D.beacon("60:60:1f:aa:bb:cc", "Mavic_Air_7721", 6)
    body = odid_ie_body()
    m = M.parse(bc + bytes([221, 4 + len(body)]) + b"\xfa\x0b\xbc\x0d" + body)
    fr = _frames([m])
    hit = [f for f in TH.findings(T.tracks(fr), fr) if f["kind"] == "drone_remote_id"]
    assert hit and hit[0]["severity"] == "high" and hit[0]["confidence"] == 0.99
    assert hit[0]["operator_location"] is not None and hit[0]["uas_id"]
    # without Remote ID the platform is still named, but only as a signature
    fr2 = _frames([M.parse(bc)])
    sig = [f for f in TH.findings(T.tracks(fr2), fr2) if f["kind"] == "drone_signature"]
    assert sig and sig[0]["severity"] == "medium" and sig[0]["confidence"] < 0.99


def test_findings_are_sorted_and_self_describing():
    fr = _frames([_beacon()] + [_mgmt("deauth", AP, "ff:ff:ff:ff:ff:ff", AP) for _ in range(30)])
    found = TH.findings(T.tracks(fr), fr)
    sev = [TH.SEVERITY_ORDER[f["severity"]] for f in found]
    assert sev == sorted(sev)
    for f in found:
        assert f["title"] and f["detail"] and f["evidence"] and 0 < f["confidence"] <= 1
        assert f["severity"] in TH.SEVERITY_ORDER


def test_plugin_reads_the_generic_report_contract():
    """It must work off signal.decoded["fields"]["ppdus"] alone -- it imports no receiver."""
    class Sig:
        def __init__(self, fields):
            self.decoded = {"method": "wifi_dsss receiver", "fields": fields}
            self.features = {}

    dsss = Sig({"ppdus": [{"start_s": 1e-3, "end_s": 1.2e-3, "snr_db": 12.0, "cfo_hz": 900.0,
                           "channel": 6, "rate_mbps": 1.0, "mac": _beacon()}]})
    ofdm = Sig({"ppdus": [{"start_us": 2000.0, "dur_us": 300.0, "snr_db": 20.0, "cfo_khz": 1.2,
                           "channel": "2g6", "fmt": "HT",
                           "macs": [_beacon(seq=5), _mgmt("probe-req", STA, "ff:ff:ff:ff:ff:ff",
                                                          "ff:ff:ff:ff:ff:ff", ssid="HomeNet")]}]})
    out = PL.DeviceTracking().postprocess(None, [dsss, ofdm])
    assert out["device_summary"]["n_frames"] == 3
    assert out["device_summary"]["n_access_points"] == 1
    assert {d["mac"] for d in out["devices"]} == {AP, STA}
    ap = next(d for d in out["devices"] if d["mac"] == AP)
    assert ap["n_frames"] == 2 and sorted(ap["decoded_by"]) == ["wifi_dsss"]
    assert dsss.decoded["fields"]["device_macs"] == [AP]
    assert out["findings"]
    assert PL.DeviceTracking().postprocess(None, []) == {}


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("PASS", name)
    fr = _frames([_beacon()] + [_mgmt("deauth", AP, "ff:ff:ff:ff:ff:ff", AP) for _ in range(30)])
    print("\nfindings on a synthetic deauthentication flood:")
    for f in TH.findings(T.tracks(fr), fr):
        print(f"  [{f['severity']:8s}] {f['title']}\n             {f['evidence'][:100]}")
