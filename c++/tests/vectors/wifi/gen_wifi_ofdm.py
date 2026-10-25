#!/usr/bin/env python3
"""Reference vectors for tests/unit/test_rx_wifi_ofdm.cpp, written by the Python 802.11 OFDM receiver.

    nice -n 19 .venv/bin/python c++/tests/vectors/wifi/gen_wifi_ofdm.py

  ofdm_blocks.json   viterbi / conv_encode / crc8_sig, Fraction.limit_denominator, channel_plan,
                     payload helpers (deinterleave_indices, depuncture, descramble, llr_exact,
                     smooth_channel, deaggregate), ldpc.params / ldpc.decode on noisy codewords
                     (one that converges, one that runs out of iterations) and ldpc.decode_stream
  ofdm_scene.cf32 + ofdm_scene.json   a synthetic capture (ofdm_tx.scene, 30.72 MSps at 2437 MHz:
                     non-HT 6/24/54 Mb/s, HT MCS7 BCC A-MPDU, HT MCS7 LDPC + STBC, HT MCS3 short GI,
                     one with a CFO) with WifiDetector.run()'s packets and plugin_ofdm._records()
  ofdm_scaled.cf32 + ofdm_scaled.json  802.11 at a 10 MHz clock off the channel grid (ofdm_tx.timescale)
                     with ofdm_scaled.find_bands / numerology / scan and the scaled records
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
sys.path.insert(0, str(REPO / "python"))

from spectrum_intel.receivers.wifi import ldpc as LD          # noqa: E402
from spectrum_intel.receivers.wifi import ofdm_rx as W        # noqa: E402
from spectrum_intel.receivers.wifi import ofdm_scaled as S    # noqa: E402
from spectrum_intel.receivers.wifi import ofdm_tx as T        # noqa: E402
from spectrum_intel.receivers.wifi import payload as P        # noqa: E402
from spectrum_intel.receivers.wifi import plugin_ofdm as PO   # noqa: E402


def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items() if not str(k).startswith("_")}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return clean(o.tolist())
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, complex):
        return [o.real, o.imag]
    return o


def cplx(a):
    a = np.asarray(a, complex)
    return [[float(v.real), float(v.imag)] for v in a]


def blocks():
    rng = np.random.default_rng(7)
    out = {}
    bits = rng.integers(0, 2, 60).tolist()
    soft = (2.0 * np.array(W.conv_encode(np.array(bits))) - 1) + rng.normal(0, 0.8, 120)
    out["viterbi"] = {"soft": soft.tolist(), "free": W.viterbi(soft).tolist(),
                      "terminated": W.viterbi(soft, terminated=True).tolist(),
                      "bits": bits, "encoded": W.conv_encode(np.array(bits)).tolist(),
                      "crc8": W.crc8_sig(bits[:34])}
    from fractions import Fraction
    out["limit_denominator"] = [[x, *(lambda f: [f.numerator, f.denominator])(Fraction(x).limit_denominator(4096))]
                                for x in (20e6 / 61.44e6, 20e6 / 30.72e6, 11e6 / 61.44e6, 10.997e6 / 61.44e6,
                                          20e6 / (61.44e6 * 20 / 10.986e6), 0.5, 1.0 / 3)]
    out["channel_plan"] = [[fc, fs, [c["name"] for c in W.channel_plan(fc, fs, 60e6)]]
                           for fc, fs in ((2437e6, 61.44e6), (5240e6, 61.44e6), (5500e6, 30.72e6))]
    out["deinterleave"] = [[n, b, c, P.deinterleave_indices(n, b, c)[1].tolist()]
                           for n, b, c in ((48, 1, 16), (192, 4, 16), (288, 6, 16), (104, 2, 13), (312, 6, 13))]
    s = rng.normal(size=36)
    out["depuncture"] = [[list(r), P.depuncture(s, r).tolist()] for r in ((3, 4), (2, 3), (5, 6))]
    out["depuncture_in"] = s.tolist()
    sb = np.array(T.scramble([0] * 7 + rng.integers(0, 2, 33).tolist(), 0x35), np.uint8)
    d, seed = P.descramble(sb)
    out["descramble"] = {"in": sb.tolist(), "out": d.tolist(), "seed": list(seed)}
    sym = (rng.normal(size=12) + 1j * rng.normal(size=12)) * 0.4
    var = rng.uniform(0.01, 0.3, 12)
    out["llr"] = {"sym": cplx(sym), "var": var.tolist(),
                  **{str(nb): P.llr_exact(sym, nb, var).reshape(-1).tolist() for nb in (1, 2, 4, 6)}}
    h = np.zeros(64, complex)
    sc = W.HT_USED_SC[2:]
    h[W._bins(sc)] = np.exp(1j * 0.3 * sc) * (1 + 0.1 * rng.normal(size=len(sc)))
    out["smooth"] = {"h": cplx(h), "sc": sc.tolist(), **{str(w): cplx(P.smooth_channel(h, sc, w)) for w in (1, 3)}}
    mp = [bytes(rng.integers(0, 256, n).astype(np.uint8)) for n in (30, 5, 61)]
    agg = bytearray(T.ampdu(mp))
    agg[40] ^= 0xFF                                            # corrupt a byte inside the first MPDU
    out["deaggregate"] = {"in": list(agg), "out": [[list(m), bool(ok)] for m, ok in P.deaggregate(bytes(agg))]}
    out["ldpc_params"] = [[n, nc, list(r), m, LD.params(n, nc, r, m)]
                          for n, nc, r, m in ((27, 312, (5, 6), 2), (100, 52, (1, 2), 1), (1500, 312, (5, 6), 1),
                                              (300, 208, (3, 4), 2), (4000, 156, (2, 3), 1), (8, 104, (1, 2), 1))]
    cws = []
    for L, r, sigma in ((648, (1, 2), 0.75), (1296, (2, 3), 0.55), (1944, (5, 6), 0.9)):
        K = L * r[0] // r[1]
        cw = LD.encode(rng.integers(0, 2, K).astype(np.uint8), L, r)
        llr = (1 - 2.0 * cw) * 2 / sigma ** 2 + rng.normal(0, 2 / sigma, L)
        hard, ok, used = LD.decode(llr, L, r)
        cws.append({"L": L, "rate": list(r), "llr": llr.tolist(), "hard": hard.tolist(), "ok": bool(ok),
                    "used": int(used), "cw": cw.tolist()})
    out["ldpc_decode"] = cws
    LP = LD.params(200, 312, (5, 6), 2)
    bits = rng.integers(0, 2, LP["npld"]).astype(np.uint8)
    tx = LD.encode_stream(bits, LP, (5, 6))
    llr = (1 - 2.0 * tx) * 6 + rng.normal(0, 2.5, len(tx))
    info, ok, used = LD.decode_stream(llr, LP, (5, 6))
    out["ldpc_stream"] = {"length": 200, "ncbps": 312, "rate": [5, 6], "m_stbc": 2, "llr": llr.tolist(),
                          "info": info.tolist(), "ok": [bool(v) for v in ok], "used": [int(v) for v in used],
                          "offsets": [int(v) for v in LD.offsets(LP)]}
    return out


def mac_frame(sub, ta, ra, body=b""):
    fc = {"beacon": b"\x80\x00", "data": b"\x08\x01", "qos": b"\x88\x01", "ack": b"\xd4\x00"}[sub]
    if sub == "ack":
        return fc + b"\x00\x00" + ra
    return fc + b"\x2c\x00" + ra + ta + ta + b"\x10\x00" + (b"\x00\x00" if sub == "qos" else b"") + body


def beacon(ta, ssid):
    body = bytes(8) + b"\x64\x00\x11\x04" + bytes([0, len(ssid)]) + ssid + b"\x01\x04\x82\x84\x8b\x96"
    return mac_frame("beacon", ta, b"\xff" * 6, body)


def scene_main(rng):
    fs, fc, n = 30.72e6, 2437e6, 40000
    ap, sta = bytes.fromhex("00c0ca112233"), bytes.fromhex("3c2ef9445566")
    pp = [
        dict(samples=T.ppdu(beacon(ap, b"si-test"), fmt="NON_HT", rate_mbps=24, scramble_seed=0x2B),
             f_hz=2437e6, t_s=40e-6, snr_db=25),
        dict(samples=T.ppdu(mac_frame("ack", sta, ap), fmt="NON_HT", rate_mbps=6, scramble_seed=0x11),
             f_hz=2437e6, t_s=150e-6, snr_db=22, cfo_hz=37e3),
        dict(samples=T.ppdu(T.ampdu([mac_frame("qos", sta, ap, bytes(range(80))),
                                     mac_frame("qos", sta, ap, bytes(range(40, 140)))]),
                            mcs=7, aggregation=True, scramble_seed=0x47),
             f_hz=2437e6, t_s=260e-6, snr_db=30),
        dict(samples=T.ppdu(mac_frame("data", ap, sta, bytes(range(120))), mcs=7, ldpc=True, stbc=True,
                            scramble_seed=0x5D), f_hz=2437e6, t_s=420e-6, snr_db=30),
        dict(samples=T.ppdu(mac_frame("data", sta, ap, bytes(range(60))), mcs=3, short_gi=True,
                            scramble_seed=0x66), f_hz=2437e6, t_s=640e-6, snr_db=24, cfo_hz=-21e3),
        dict(samples=T.ppdu(mac_frame("data", ap, sta, bytes(range(30))), fmt="NON_HT", rate_mbps=54,
                            scramble_seed=0x0F), f_hz=2437e6, t_s=1000e-6, snr_db=30),
    ]
    x = T.scene(pp, n, fs=fs, fc=fc, seed=3)
    det = W.WifiDetector(fs, fc)
    pk = det.run(x)
    recs = PO._records(pk, SimpleNamespace(fs=fs))
    PO._attribute_vht(recs)
    pkts = [{"channel": p.channel, "center_mhz": p.center_mhz, "start_sample": p.start_sample,
             "end_sample": p.end_sample, "start_us": p.start_us, "end_us": p.end_us, "fmt": p.fmt,
             "lsig_rate_mbps": p.lsig_rate_mbps, "lsig_length": p.lsig_length, "f_lo_mhz": p.f_lo_mhz,
             "f_hi_mhz": p.f_hi_mhz, "snr_db": p.snr_db, "cfo_khz": p.cfo_khz, "stf_metric": p.stf_metric,
             "mcs": p.mcs, "sig_crc_ok": p.sig_crc_ok, "energy_end_us": p.energy_end_us,
             "lsig_evm_db": p.extra["lsig_evm_db"], "lsig_bit_agreement": p.extra["lsig_bit_agreement"]}
            for p in pk]
    x.astype(np.complex64).tofile(HERE / "ofdm_scene.cf32")
    json.dump(clean({"fs": fs, "fc": fc, "n_rejects": len(det.rejects), "packets": pkts, "records": recs}),
              open(HERE / "ofdm_scene.json", "w"), indent=0)
    return len(pk), [r.get("fcs_ok") for r in recs]


def scene_scaled(rng):
    fs, fc, n, clock = 30.72e6, 2100e6, 52000, 10e6
    ap, sta = bytes.fromhex("00c0ca778899"), bytes.fromhex("001122aabbcc")
    pp = []
    t = 30e-6
    for i in range(7):
        if i % 2 == 0:
            s = T.ppdu(mac_frame("data", ap, sta, bytes(range(90 + 10 * i))), fmt="NON_HT", rate_mbps=12,
                       scramble_seed=0x20 + i)
        else:
            s = T.ppdu(mac_frame("data", sta, ap, bytes(range(70))), mcs=2, scramble_seed=0x30 + i)
        y = T.timescale(s, clock, fs=20e6)
        pp.append(dict(samples=y, f_hz=2096.3e6, t_s=t, snr_db=24))
        t += len(y) / 20e6 * 1.0 + 40e-6
    x = T.scene(pp, n, fs=fs, fc=fc, seed=5)
    bands, unions = S.find_bands(x, fs, fc)
    nms = [S.numerology(x, fs, fc, b) for b in (bands[:12] + unions[:4])]
    links = S.scan(x, fs, fc)
    out = {"fs": fs, "fc": fc,
           "bands": [{"f_lo": b["f_lo"], "f_hi": b["f_hi"], "n_frames": b["n_frames"]} for b in bands],
           "unions": [{"f_lo": b["f_lo"], "f_hi": b["f_hi"], "n_frames": b["n_frames"]} for b in unions],
           "numerology": [nm for nm in nms],
           "links": [{"center_hz": ln.center_hz, "clock_hz": ln.clock_hz, "inverted": ln.inverted,
                      "n_rejects": ln.n_rejects, "records": PO._records(ln.packets, SimpleNamespace(fs=fs), ln.clock_hz)}
                     for ln in links]}
    x.astype(np.complex64).tofile(HERE / "ofdm_scaled.cf32")
    json.dump(clean(out), open(HERE / "ofdm_scaled.json", "w"), indent=0)
    return [(ln.clock_hz, len(ln.packets)) for ln in links]


if __name__ == "__main__":
    rng = np.random.default_rng(1)
    json.dump(clean(blocks()), open(HERE / "ofdm_blocks.json", "w"))
    print("scene:", scene_main(rng))
    print("scaled:", scene_scaled(rng))
