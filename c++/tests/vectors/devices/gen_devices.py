"""Reference vectors for receivers/devices (C++ test test_rx_devices.cpp).

Builds synthetic decoded-802.11 signals (the generic report contract: decoded.fields.ppdus) that
exercise every tracker path and every threat rule, runs python receivers.devices on them and
writes cases.json: per case the input signals, the postprocess extras, the per-signal device_macs
and a finalize_report() result. Run: .venv/bin/python c++/tests/vectors/devices/gen_devices.py
"""
import copy
import json
from pathlib import Path
from types import SimpleNamespace

from spectrum_intel.receivers.devices.plugin import DeviceTracking

AP1, AP2, AP3 = "00:11:22:33:44:55", "02:11:22:33:44:56", "a4:2b:b0:10:20:30"
STA1, STA2, RND1, RND2, RND3 = "74:0e:a4:9d:0e:49", "3c:22:fb:01:02:03", "da:a1:19:00:00:01", "de:ad:be:00:00:02", "6a:00:00:00:00:03"
BC = "ff:ff:ff:ff:ff:ff"


def beacon(bssid, ssid, sec, **kw):
    m = {"type": "mgmt", "subtype": "beacon", "bssid": bssid, "ta": bssid, "da": BC, "ssid": ssid,
         "security": sec, "ta_vendor": kw.pop("vendor", None), "seq": 10, "length_bytes": 120}
    m.update(kw)
    return m


def probe(ta, ssid, ids, **kw):
    m = {"type": "mgmt", "subtype": "probe-req", "ta": ta, "da": BC, "ssid": ssid, "element_ids": ids,
         "length_bytes": 80}
    m.update(kw)
    return m


def ofdm(m, t_us, snr=20.0, cfo_khz=1.5, ch=36, rate=6.0, **kw):
    p = {"channel": ch, "fmt": "L-SIG", "rate_mbps": rate, "snr_db": snr, "cfo_khz": cfo_khz,
         "start_us": t_us, "dur_us": 100.0, "mac": m, "fcs_ok": True}
    p.update(kw)
    return p


def dsss(m, t_s, snr=12, cfo=-3000, ch=6, rate=1.0, **kw):
    p = {"channel": ch, "rate_mbps": rate, "snr_db": snr, "cfo_hz": cfo, "start_s": t_s,
         "end_s": t_s + 0.0012, "mac": m}
    p.update(kw)
    return p


def sig(method, ppdus, sid):
    return {"method": method, "fields": {"ppdus": ppdus} if ppdus is not None else {}, "id": sid}


def case_mixed():
    s0 = [ofdm(beacon(AP1, "Home", "WPA2-PSK (CCMP)", vendor="Acme", country={"code": "IN"},
                      ap_uptime=1234, bss_load={"stations": 3}, phy="802.11n"), 100.0),
          ofdm(beacon(AP2, "Home", "WPA2-PSK (CCMP)", vendor="Acme"), 300.0, snr=18.5),
          ofdm(beacon(AP3, "Cafe", "open", vendor="Tenda, Inc."), 500.0, snr=25.25, ch=40, rate=12.0),
          ofdm(beacon("b0:00:00:00:00:01", "Cafe", "WPA2-PSK (CCMP)", vendor="Netgear"), 700.0, snr=9.0),
          ofdm(beacon("b0:00:00:00:00:02", "Old", "WPA (TKIP)", vendor="Linksys",
                      wps={"state": "configured", "selected_registrar": True, "device_name": "WRT54G"}), 900.0),
          ofdm(beacon("b0:00:00:00:00:03", "Mixed", "WPA2 (CCMP+TKIP)", vendor="Linksys"), 1100.0),
          ofdm(beacon("b0:00:00:00:00:04", "", "WPA2-PSK (CCMP)", ssid_hidden=True), 1300.0, snr=None),
          ofdm({"type": "data", "subtype": "qos-data", "ta": STA1, "ra": AP1, "bssid": AP1, "retry": True,
                "pwr_mgmt": True, "encryption": {"pn": 77}, "ta_vendor": "Apple, Inc.", "length_bytes": 1500},
               1500.0, snr=27.6, rate=54.0),
          ofdm({"type": "data", "subtype": "qos-data", "ta": STA1, "ra": AP1, "bssid": AP1,
                "encryption": {"pn": 79}, "length_bytes": 600}, 1700.0, snr=28.3, rate=54),
          ofdm({"type": "ctrl", "subtype": "rts", "ta": STA2, "ra": AP1}, 1900.0, snr=10, rate=24.0),
          ofdm({"type": "data", "subtype": "data", "ta": AP1, "ra": STA2,
                "eapol": {"handshake": "4-way 1/4"}, "sa": AP1, "da": STA2}, 2100.0, rate=None,
               ht_rate_mbps=13.0),
          ofdm({"type": "data", "subtype": "data", "ta": STA2, "ra": AP1,
                "eapol": {"handshake": "4-way 2/4"}}, 2300.0),
          ofdm({"type": "data", "subtype": "data", "ta": STA2, "ra": AP1, "eapol": {"eap": "EAP-Identity"}},
               2500.0, fcs_ok=False),
          ]
    ampdu = ofdm(None, 2700.0, snr=21.0, rate=None, ht_rate_mbps=65.0)
    ampdu.pop("mac")
    ampdu["macs"] = [{"type": "data", "subtype": "qos-data", "ta": STA1, "ra": AP1, "seq": 5, "length_bytes": 100},
                     {"type": "data", "subtype": "qos-data", "ta": STA1, "ra": AP1, "seq": 6, "length_bytes": 101}]
    s0.append(ampdu)
    s1 = [dsss(probe(RND1, "Airport WiFi", [0, 1, 50, 3, 45, 127], ext_element_ids=[35],
                     vendor_elements=["Microsoft", "Apple", "Apple"], ht={"spatial_streams": 2, "40mhz": True},
                     rates_mbps=[1, 2, 5.5, 11]), 0.001),
          dsss(probe(RND1, "<broadcast>", [0, 1, 50, 3, 45, 127], ext_element_ids=[35],
                     vendor_elements=["Microsoft", "Apple", "Apple"], ht={"spatial_streams": 2, "40mhz": True},
                     rates_mbps=[1, 2, 5.5, 11]), 0.004, snr=14),
          dsss(probe(RND2, "Hotel 'Guest'", [0, 1, 50, 3, 45, 127], ext_element_ids=[35],
                     vendor_elements=["Microsoft", "Apple", "Apple"], ht={"spatial_streams": 2, "40mhz": True},
                     rates_mbps=[1, 2, 5.5, 11], p2p={"device_name": "Pixel 7"}), 0.008, cfo=-2000),
          dsss(probe(RND3, None, [0, 1, 50], ht={"spatial_streams": 1}), 0.010, cfo=None, snr=None),
          dsss(probe("6e:00:00:00:00:04", None, [0, 1, 50], ht={"spatial_streams": 1}), 0.011, cfo=None),
          dsss(probe("01:00:5e:00:00:01", "x", [0]), 0.012),
          dsss({"type": "mgmt", "subtype": "probe-resp", "bssid": AP3, "ssid": "Cafe", "da": RND1}, 0.013,
               end_s=None),
          ]
    s1 += [dsss({"type": "mgmt", "subtype": "deauth", "ta": "c0:ff:ee:00:00:01", "da": BC if i == 3 else STA1,
                 "addr2": "x"}, 0.02 + i * 1e-3) for i in range(9)]
    s1 += [dsss({"type": "mgmt", "subtype": "disassoc", "sa": "c0:ff:ee:00:00:02", "ra": f"10:00:00:00:00:{i:02x}"},
                0.04 + i * 1e-3) for i in range(26)]
    s1 += [dsss({"type": "mgmt", "subtype": "deauth", "addr2": "c0:ff:ee:00:00:03", "addr1": STA2}, 0.07 + i * 1e-3)
           for i in range(8)]
    s1 += [dsss({"type": "mgmt", "subtype": "deauth", "ta": "c0:ff:ee:00:00:04"}, 0.08)]
    s2 = [ofdm(beacon("60:60:1f:00:00:01", "TELLO-ABC", "open", vendor="SZ DJI",
                      drone={"platform": "DJI Tello", "confidence": 0.8, "evidence": "OUI + SSID"}), 50.0),
          ofdm({"type": "mgmt", "subtype": "beacon", "bssid": "60:60:1f:00:00:01", "ta": "60:60:1f:00:00:01",
                "drone": {"platform": "DJI Tello", "confidence": 0.6}}, 60.0),
          ofdm(beacon("60:60:1f:00:00:02", "RID-1", "open",
                      drone={"platform": "DJI Mini 4", "confidence": 0.95, "evidence": "Remote ID"},
                      remote_id={"uas_id": "1581F5BK", "lat": 12.9715987, "lon": 77.5945627, "height_m": 42.5,
                                 "operator_lat": 12.97, "operator_lon": 77.59}), 70.0),
          ofdm(beacon("60:60:1f:00:00:03", "RID-2", "open",
                      drone={"platform": "Parrot", "confidence": 0.9},
                      remote_id={"uas_id": None, "lat": 1, "lon": 2}), 80.0),
          ofdm(beacon("60:60:1f:00:00:04", "RID-3", "open",
                      drone={"platform": "X", "confidence": 0.9}, remote_id={"serial": "abc"}), 90.0),
          ]
    return [sig("receivers.wifi ofdm_rx + payload", s0, 0), sig("receivers.wifi_dsss dsss", s1, 1),
            sig("receivers.dji something", None, 2), sig("receivers.wifi ofdm", s2, 3)]


def case_many_networks():
    s = [ofdm(beacon(f"0{i % 4}:aa:bb:cc:dd:{i:02x}", f"net{i:02d}", "WPA3-SAE (GCMP)"), 10.0 * i, snr=5 + i)
         for i in range(11)]
    s.append(ofdm(beacon("0c:aa:bb:cc:dd:ee", "<len 0>", "WPA2-PSK (CCMP)"), 300.0, snr=5))
    return [sig("receivers.wifi a", s, 0), sig("receivers.wifi b", [], 1)]


def case_twins():
    s = [ofdm(beacon("10:00:00:00:00:01", "Corp", "WPA2-PSK (CCMP)", vendor="Cisco, Inc."), 1.0),
         ofdm(beacon("20:00:00:00:00:01", "Corp", "WPA2-PSK (CCMP)", vendor="Cisco Systems"), 2.0),
         ofdm(beacon("30:00:00:00:00:01", "Mesh", "WPA2-PSK (CCMP)", vendor="Eero"), 3.0),
         ofdm(beacon("40:00:00:00:00:01", "Mesh", "WPA2-PSK (CCMP)", vendor="Eero"), 4.0),
         ofdm(beacon("50:00:00:00:00:01", "Corp2", "WPA2-PSK (CCMP)", vendor="A"), 5.0),
         ofdm(beacon("60:00:00:00:00:01", "Corp2", "WEP", vendor="B"), 6.0),
         ofdm(beacon("70:00:00:00:00:01", "Corp2", None, vendor=None), 7.0),
         ofdm(beacon("02:00:00:00:00:01", "Split", "open"), 8.0),
         ofdm(beacon("00:00:00:00:01:02", "Split", "open"), 9.0),
         ofdm(beacon("00:00:00:00:00:02", "Split", "open"), 10.0),
         ]
    return [sig("receivers.wifi x", s, 0)]


def case_w_batch():
    """PLAN W6 / W7 / W1: a beacon flood (unregistered BSSIDs from unrelated prefixes), two randomised MACs with one
    probe signature whose sequence counter continues, a third whose does not, and a behaviour drone link."""
    fl = [ofdm(beacon(f"{(4 * i + 2) & 0xFE | 2:02x}:{i * 37 % 256:02x}:9c:00:00:{i:02x}", f"free{i}", "open",
                      vendor="unknown OUI"), 10.0 * i) for i in range(16)]
    pr = [dsss(probe("da:a1:19:00:00:11", None, [0, 1, 50, 45], rates_mbps=[1, 2], seq=100 + k), 0.001 * k, cfo=-1000)
          for k in range(3)]
    pr += [dsss(probe("da:a1:19:00:00:12", None, [0, 1, 50, 45], rates_mbps=[1, 2], seq=104 + k), 0.01 + 0.001 * k,
                cfo=-1100) for k in range(3)]
    pr += [dsss(probe("da:a1:19:00:00:13", None, [0, 1, 50, 127], rates_mbps=[1], seq=900 + k), 0.02 + 0.001 * k)
           for k in range(2)]
    pr += [dsss(probe("da:a1:19:00:00:14", None, [0, 1, 50, 127], rates_mbps=[1], seq=3000 + k), 0.03 + 0.001 * k)
           for k in range(2)]
    dl = [{"bssid": "00:e0:4c:0a:0b:0c", "ssid": "Drone-Cam-x", "clients": ["3c:2e:f9:44:55:66"],
           "verdict": "drone / FPV video link-like", "platform": None, "basis": "behaviour", "confidence": 0.55,
           "behaviour": {"one_client": True, "downlink_dominant": True, "sustained_downlink": True,
                         "uplink_small": True, "video_payload": False}},
          {"bssid": "60:60:1f:00:00:09", "ssid": None, "clients": [], "verdict": "drone link", "platform": "DJI-like",
           "basis": "identity", "confidence": 0.75, "behaviour": {}}]
    sg = sig("receivers.wifi ofdm", fl, 0)
    sg2 = sig("receivers.wifi_dsss dsss", pr, 1)
    sg3 = sig("receivers.wifi ofdm", [], 2)
    sg3["fields"]["drone_links"] = dl
    return [sg, sg2, sg3]


def case_empty():
    return [sig("receivers.wifi x", [ofdm({"type": "data", "subtype": "data", "ta": BC}, 1.0)], 0)]


def run(signals):
    plug = DeviceTracking()
    sigs = [SimpleNamespace(decoded={"method": s["method"], "fields": copy.deepcopy(s["fields"])}) for s in signals]
    extras = json.loads(json.dumps(plug.postprocess(None, sigs)))
    macs = [s.decoded["fields"].get("device_macs") for s in sigs]
    rep = {"signals": [{"id": s["id"], "decoded": {"fields": ({"device_macs": m} if m else {})}}
                       for s, m in zip(signals, macs)], **copy.deepcopy(extras)}
    plug.finalize_report(rep)
    rep = json.loads(json.dumps(rep))
    return {"signals": signals, "extras": extras, "device_macs": macs, "finalized": rep}


cases = {"mixed": run(case_mixed()), "many_networks": run(case_many_networks()), "twins": run(case_twins()),
         "empty": run(case_empty()), "w_batch": run(case_w_batch())}
out = Path(__file__).with_name("cases.json")
out.write_text(json.dumps(cases, indent=None))
for k, v in cases.items():
    print(k, len(v["extras"].get("devices", [])), "devices",
          [f["kind"] for f in v["extras"].get("findings", [])])
