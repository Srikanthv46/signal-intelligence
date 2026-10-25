#!/usr/bin/env python3
"""Reference for tests/unit/test_wifi_field.cpp (PLAN batch W, Wi-Fi field gaps), written by the Python receiver:

  he_parse          HE-SIG-A bit words (SU / MU / TB, CRC failure) -> ofdm_rx.parse_he_sig_a
  mpdu              W3 / W4: beacons with drone / DIY hotspot SSIDs, Tello SDK command / state, H.264 / H.265 / RTP video
                    datagrams -> mac80211.parse (ipv4, drone, summary)
  links             W1: field.links on frame lists (drone-style link, home AP, Tello identity, empty)
  link_scene.cf32   W1 end to end: an encrypted drone-style link (24 Mb/s) -> records, summarise(block_s): wifi_links,
  + link_scene.json drone_links, label, evidence
  he_scene.cf32     a synthetic capture (30.72 MSps at 2437 MHz) with an HE beacon BSS (HE Operation colour 21) and
  + he_scene.json   HE SU / MU / TB preambles (ofdm_tx.he_preamble) -> WifiDetector packets, records, summary fields

    nice -n 19 .venv/bin/python -I c++/tests/vectors/wifi_field/gen_wifi_field.py
"""
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[3] / "python"))
import struct                                                  # noqa: E402
from spectrum_intel.receivers.wifi import mac80211 as M       # noqa: E402
from spectrum_intel.receivers.wifi import ofdm_rx as W        # noqa: E402
from spectrum_intel.receivers.wifi import ofdm_tx as T        # noqa: E402
from spectrum_intel.receivers.wifi import plugin_ofdm as PO   # noqa: E402


def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items() if not str(k).startswith("_")}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return clean(o.tolist())
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    return o


HE_KEYS = ("he_format", "uplink", "bss_color", "bw", "ldpc", "dcm", "mcs", "nsts", "txop_us")


def he_parse():
    cases = []
    for fmt, kw, mu in (("SU", dict(mcs=11, bss_color=63, uplink=True, nsts=2, txop=1 + 2 * 5, ldpc=1), False),
                        ("SU", dict(mcs=0, bss_color=0, bw=3, gi_ltf=3, txop=40, ldpc=0), False),
                        ("MU", dict(bss_color=9, bw=5, uplink=False, txop=127), True),
                        ("MU", dict(bss_color=33, bw=1, uplink=True, txop=2 * 63), True),
                        ("TB", dict(bss_color=17, bw=2), False)):
        b = T.he_sig_a_bits(fmt, **kw)
        length = 41 if mu else 40
        cases.append({"bits": b, "lsig_length": length, "out": W.parse_he_sig_a(b, length)})
    bad = T.he_sig_a_bits("SU", bss_color=5)
    bad[10] ^= 1                                               # CRC-4 failure
    cases.append({"bits": bad, "lsig_length": 40, "out": None})
    tail = T.he_sig_a_bits("SU", bss_color=5)
    tail[50] = 1                                               # tail failure
    cases.append({"bits": tail, "lsig_length": 40, "out": None})
    return cases


def beacon_he(ta, ssid, color):
    he_op = bytes([255, 7, 36, 0, 0, 0, color, 0xFC, 0xFF])   # HE Operation: params, BSS colour, basic HE-MCS
    body = bytes(8) + b"\x64\x00\x11\x04" + bytes([0, len(ssid)]) + ssid + b"\x01\x04\x8c\x12\x98\x24" + he_op
    return b"\x80\x00\x00\x00" + b"\xff" * 6 + ta + ta + b"\x10\x00" + body


def he_scene():
    fs, fc, n = 30.72e6, 2437e6, 60000
    ap = bytes.fromhex("a0b1c2d3e4f5")
    pp = [dict(samples=T.ppdu(beacon_he(ap, b"si-he", 21), fmt="NON_HT", rate_mbps=6, scramble_seed=0x21),
               f_hz=2437e6, t_s=40e-6, snr_db=25),
          dict(samples=T.he_preamble(T.he_sig_a_bits("SU", mcs=7, bss_color=21, txop=2 * 12), seed=1),
               f_hz=2437e6, t_s=520e-6, snr_db=22, cfo_hz=31e3),
          dict(samples=T.he_preamble(T.he_sig_a_bits("SU", mcs=4, bss_color=21, uplink=True, nsts=1), seed=2),
               f_hz=2437e6, t_s=700e-6, snr_db=18),
          dict(samples=T.he_preamble(T.he_sig_a_bits("MU", bss_color=21, bw=0, txop=1 + 2 * 2), mu=True, seed=3),
               f_hz=2437e6, t_s=880e-6, snr_db=24, cfo_hz=-17e3),
          dict(samples=T.he_preamble(T.he_sig_a_bits("TB", bss_color=9, bw=0), seed=4),
               f_hz=2437e6, t_s=1060e-6, snr_db=20)]
    x = T.scene(pp, n, fs=fs, fc=fc, seed=11)
    det = W.WifiDetector(fs, fc)
    pk = det.run(x)
    recs = PO._records(pk, SimpleNamespace(fs=fs))
    PO._attribute_vht(recs)
    pkts = [{"fmt": p.fmt, "start_sample": p.start_sample, "lsig_length": p.lsig_length, "mcs": p.mcs,
             "sig_crc_ok": p.sig_crc_ok, "bandwidth_mhz": p.bandwidth_mhz, "f_lo_mhz": p.f_lo_mhz,
             "f_hi_mhz": p.f_hi_mhz} for p in pk]
    import copy
    res = PO.WifiOfdmDecoder().summarise(copy.deepcopy(recs), len(recs), 0, recs)
    summ = {k: v for k, v in res["fields"].items() if k.startswith("he_")}
    x.astype(np.complex64).tofile(HERE / "he_scene.cf32")
    return {"fs": fs, "fc": fc, "packets": pkts, "records": recs, "summary": summ}


def ip_udp(payload, sport, dport):
    u = struct.pack(">HHHH", sport, dport, 8 + len(payload), 0) + payload
    return bytes([0x45, 0, 0, 20 + len(u), 0, 0, 0, 0, 64, 17, 0, 0, 192, 168, 10, 2, 192, 168, 10, 1]) + u


def data_frame(ip, tods=True):
    a1, a2 = bytes.fromhex("60601f112233"), bytes.fromhex("5ccf7f010203")
    if not tods:
        a1, a2 = a2, a1
    return (bytes([0x08, 0x01 if tods else 0x02, 0, 0]) + a1 + a2 + bytes.fromhex("60601f112233") + b"\x10\x00"
            + b"\xaa\xaa\x03\x00\x00\x00\x08\x00" + ip)


def beacon(ssid, ta=bytes.fromhex("0a1b2c3d4e5f")):
    body = bytes(8) + b"\x64\x00\x01\x04" + bytes([0, len(ssid)]) + ssid + b"\x01\x04\x82\x84\x8b\x96"
    return b"\x80\x00\x00\x00" + b"\xff" * 6 + ta + ta + b"\x10\x00" + body


STATE = (b"mid:-1;x:0;y:0;z:0;mpry:0,0,0;pitch:1;roll:-2;yaw:45;vgx:0;vgy:0;vgz:0;templ:60;temph:62;tof:52;h:40;"
         b"bat:87;baro:12.35;time:15;agx:-3.00;agy:1.00;agz:-999.00;\r\n")


def mpdus():
    f = {"ssid_" + s.decode().replace(" ", "_"): beacon(s) for s in
         (b"RMTT-9AFF2A", b"HolyStoneFPV-250ce5", b"holystone fpv-1a2b3c", b"WiFi-720P-1B2C3D", b"ExpressLRS TX",
          b"ExpressLRS RX", b"ArduPilot", b"PixRacer", b"CUAVWLINK1234", b"DroneBridge for ESP32", b"ArduPilotX",
          b"MyHome-5G", b"Spark-RC1A2B3C", b"Mavic-Mini-1a2b3c", b"CGO3P_a1b2c3", b"RUKO-F11-12345678-5G",
          b"HolyStoneFPV_12345", b"Mavic-Pro-x")}
    for ssid, ta in ((b"Spark-1A2B3C", "60601f000001"), (b"MavicAir-1A2B3C", "60601f000002"),
                     (b"PHANTOM3_09e35b", "60601f000003"), (b"CGO3_a1b2c3", "e0b6f5800001"),
                     (b"ANAFI-1234567", "903ae6000001"), (b"TELLO-ABCDEF", "60601f000004"),
                     (b"S2-anything", "381d14000001"), (b"H501A-1", "98aafc700001")):
        f["oui_" + ssid.decode()] = beacon(ssid, bytes.fromhex(ta))
    f.update(tello_cmd=data_frame(ip_udp(b"takeoff", 50000, 8889)),
             tello_rc=data_frame(ip_udp(b"rc 0 0 50 -20\r\n", 50000, 8889)),
             tello_wifi=data_frame(ip_udp(b"wifi myssid secret123", 50000, 8889)),
             tello_ok=data_frame(ip_udp(b"ok", 8889, 50000), tods=False),
             tello_batt=data_frame(ip_udp(b"87\r\n", 8889, 50000), tods=False),
             tello_state=data_frame(ip_udp(STATE, 8889, 8890), tods=False),
             tello_state_short=data_frame(ip_udp(b"pitch:1;roll:2;bat:50;", 8889, 8890), tods=False),
             tello_noise=data_frame(ip_udp(bytes(range(60)), 1234, 8889)),
             h264_sps=data_frame(ip_udp(b"\0\0\0\x01\x67\x4d\x40\x1f" + bytes(40), 62512, 11111), tods=False),
             h264_slice3=data_frame(ip_udp(b"\0\0\x01\x41\x9a\x02" + bytes(40), 8800, 40000), tods=False),
             h265_vps=data_frame(ip_udp(b"\0\0\0\x01\x40\x01\x0c" + bytes(40), 5600, 5600), tods=False),
             rtp_fua=data_frame(ip_udp(bytes([0x80, 96, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 0x7C, 0x85]) + bytes(30), 5004, 5004)),
             rtp_pt0=data_frame(ip_udp(bytes([0x80, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 0x7C, 0x85]) + bytes(30), 5004, 5004)))
    out = {}
    for k, b in f.items():
        info = M.parse(b)
        out[k] = {"hex": b.hex(), "out": {key: info.get(key) for key in ("ipv4", "drone", "summary", "ssid")}}
    return out


def link_cases():
    sys.path.insert(0, str(HERE.parents[3] / "python" / "tests" / "receivers"))
    import test_wifi_field as TF                                   # the W1 test's frame builders
    from spectrum_intel.receivers.wifi import field as FD
    sdk = bytes([0x45, 0, 0, 35, 0, 0, 0, 0, 64, 17, 0, 0, 192, 168, 10, 2, 192, 168, 10, 1]) \
        + struct.pack(">HHHH", 50000, 8889, 15, 0) + b"takeoff"
    cmd = bytes([0x08, 0x01, 0, 0]) + TF.AP + TF.PHONE + TF.AP + b"\x10\x00" + b"\xaa\xaa\x03\x00\x00\x00\x08\x00" + sdk
    sets = {"drone_home": TF.drone_block() + TF.home_block(), "tello": TF.drone_block(8) + [TF.fr(cmd, 0.03, 60)],
            "home_only": TF.home_block(), "empty": []}
    out = {}
    for k, fr in sets.items():
        fr = json.loads(json.dumps(clean(fr)))
        out[k] = {"frames": fr, "block_s": 0.034, "out": FD.links(fr, 0.034)}
    return out, TF


def link_scene(TF):
    fs, fc = 30.72e6, 2437e6
    pp, t = [], 60e-6
    pp.append(dict(samples=T.ppdu(TF.beacon(TF.AP, b"Drone-Cam-x"), fmt="NON_HT", rate_mbps=6, scramble_seed=0x11),
                   f_hz=fc, t_s=t, snr_db=25))
    t += 600e-6
    for k in range(14):
        pp.append(dict(samples=T.ppdu(TF.data(TF.AP, TF.PHONE, TF.AP, False, 1400, seq=k), fmt="NON_HT", rate_mbps=24,
                                      scramble_seed=1 + k), f_hz=fc, t_s=t, snr_db=24))
        t += 560e-6
        pp.append(dict(samples=T.ppdu(TF.ack(TF.AP), fmt="NON_HT", rate_mbps=24, scramble_seed=40 + k), f_hz=fc,
                       t_s=t, snr_db=20))
        t += 80e-6
    n = int((t + 200e-6) * fs)
    x = T.scene(pp, n, fs=fs, fc=fc, seed=4)
    recs = PO._records(W.WifiDetector(fs, fc).run(x), SimpleNamespace(fs=fs))
    import copy
    res = PO.WifiOfdmDecoder().summarise(copy.deepcopy(recs), len(recs), 0, recs, n / fs)
    x.astype(np.complex64).tofile(HERE / "link_scene.cf32")
    f = res["fields"]
    return {"fs": fs, "fc": fc, "block_s": n / fs, "label": res["label"], "evidence": res["evidence"],
            "wifi_links": f.get("wifi_links"), "drone_links": f.get("drone_links")}


if __name__ == "__main__":
    lc, TF = link_cases()
    out = {"he_parse": he_parse(), "he_scene": he_scene(), "mpdu": mpdus(), "links": lc, "link_scene": link_scene(TF)}
    json.dump(clean(out), open(HERE / "ref.json", "w"), indent=0)
    print("he scene:", [(p["fmt"], p["sig_crc_ok"]) for p in out["he_scene"]["packets"]], out["he_scene"]["summary"])
