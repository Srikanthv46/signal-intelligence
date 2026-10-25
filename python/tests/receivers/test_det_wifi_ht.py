#!/usr/bin/env python3
"""
test_det_wifi_ht.py -- 802.11n (HT) single-stream payload decoding (receivers/wifi/payload.py)
and the reference OFDM transmitter (receivers/wifi/ofdm_tx.py).

  1. Every single-stream rate (MCS 0-7) x both guard intervals: detected, HT-SIG read correctly,
     MPDU byte-exact.
  2. A-MPDU: delimiters, de-aggregation and a separate FCS per subframe, byte-exact.
  3. Sensitivity vs SNR, and tolerance to a carrier offset.
  4. What one antenna cannot do (2+ streams, 40 MHz, extension streams) is declined with a reason; LDPC
     and STBC are decoded (test_det_wifi_scaled.py).
  4b. Real 802.11n frames from the 2026-10-04 sweep decode (HT pilots rotate per symbol).
  5. The Viterbi rewrite decodes identically to the textbook formulation.
  6. Non-HT still decodes through the same transmitter, so ofdm_tx is a check on both chains.

pytest collects the test_* functions; `python test_det_wifi_ht.py` prints a rate/SNR table too.
"""
from __future__ import annotations

import numpy as np

from spectrum_intel.receivers.wifi import dsss_tx as D
from spectrum_intel.receivers.wifi import ofdm_rx as W
from spectrum_intel.receivers.wifi import ofdm_tx as TX
from spectrum_intel.receivers.wifi import payload as P

FS, FC = 20e6, 2437e6
CHAN = dict(name="2g6", ch=6, freq_hz=FC, offset_hz=0.0, usable_half_hz=10e6, partial=False)
MPDU = D.beacon("24:a3:f0:11:22:33", "HT-lab", 6, seq=11)


def _capture(samples, snr_db=30, cfo_hz=0.0, seed=0):
    # Leave the capture mostly idle: the detector takes its noise floor from a power percentile.
    n = int(len(samples) * 1.6) + 8000
    return TX.scene([dict(samples=samples, f_hz=FC, t_s=1000 / FS, snr_db=snr_db, cfo_hz=cfo_hz)],
                    n, FS, FC, seed=seed)


def _detect(samples, **kw):
    det = W.WifiDetector(FS, FC, channels=[CHAN])
    return det.run(_capture(samples, **kw))


def _decode(p):
    rx, e = p.extra["_rx"], p.extra
    return P.decode_ht_payload(rx[0], p.mcs, e["ht_length"], bool(e.get("short_gi")),
                               bool(e.get("aggregation")), n_ess=e.get("n_ess", 0),
                               chan_offset_hz=rx[1], usable_half_hz=rx[2])


def _one_ht(samples, **kw):
    ht = [p for p in _detect(samples, **kw) if p.fmt == "HT"]
    assert len(ht) == 1, f"expected one HT PPDU, got {len(ht)}"
    return ht[0]


def test_every_single_stream_rate_and_guard_interval():
    for mcs in range(8):
        for sgi in (False, True):
            p = _one_ht(TX.ppdu(MPDU, mcs=mcs, short_gi=sgi))
            assert p.mcs == mcs and bool(p.extra["short_gi"]) == sgi, (mcs, sgi)
            assert p.bandwidth_mhz == 20 and p.extra["stbc"] == 0 and p.extra["n_ess"] == 0
            r = _decode(p)
            assert r["n_fcs_ok"] == 1, (mcs, sgi)
            assert bytes.fromhex(r["subframes"][0]["mpdu_hex"]) == MPDU, (mcs, sgi)


def test_ampdu_deaggregation():
    parts = [D.beacon("24:a3:f0:11:22:33", f"agg{i}", 6, seq=100 + i) for i in range(3)]
    parts.append(D.data_frame(n_payload=120, seq=9))
    p = _one_ht(TX.ppdu(TX.ampdu(parts), mcs=5, aggregation=True))
    assert p.extra["aggregation"] == 1
    r = _decode(p)
    assert r["aggregated"] and r["n_subframes"] == 4 and r["n_fcs_ok"] == 4
    assert all(sf["delimiter_crc_ok"] for sf in r["subframes"])
    assert [bytes.fromhex(sf["mpdu_hex"]) for sf in r["subframes"]] == parts
    assert [sf["mac"]["ssid"] for sf in r["subframes"][:3]] == ["agg0", "agg1", "agg2"]
    assert r["subframes"][3]["mac"]["subtype"] == "qos-data"


def test_sensitivity_and_carrier_offset():
    """A 61-byte MCS1 frame (what a real 802.11n beacon looks like) and a long 16-QAM A-MPDU."""
    small = D.data_frame(n_payload=29, seq=3)[:57]
    for snr, need in ((8, 1.0), (10, 1.0)):
        ok = 0
        for t in range(8):
            p = [q for q in _detect(TX.ppdu(small, mcs=1, short_gi=True), snr_db=snr, seed=t)
                 if q.fmt == "HT"]
            ok += bool(p) and _decode(p[0])["n_fcs_ok"] == 1
        assert ok / 8 >= need, (snr, ok)
    big = TX.ampdu([D.data_frame(n_payload=200, seq=i, seed=i) for i in range(8)])
    p = _one_ht(TX.ppdu(big, mcs=3, aggregation=True), snr_db=18, cfo_hz=60e3)
    assert _decode(p)["n_fcs_ok"] == 8


def test_declines_what_one_antenna_cannot_decode():
    for mcs, kw, want in ((9, {}, "spatial"), (14, {}, "spatial"), (3, dict(stbc=2), "STBC=2"),
                          (3, dict(n_ess=1), "extension"), (3, dict(bw_mhz=40), "40 MHz")):
        ok, why = P.ht_decodable(mcs, **kw)
        assert not ok and want in why, (mcs, kw, why)
    # one spatial stream is decodable whatever the coding, and with Alamouti STBC (test_det_wifi_scaled.py)
    assert P.ht_decodable(7)[0] and P.ht_decodable(7, ldpc=1)[0] and P.ht_decodable(7, ldpc=1, stbc=1)[0]
    # the budget guard keeps one PPDU from monopolising the run
    try:
        P.decode_ht_payload({}, 0, P.HT_MAX_LENGTH + 1)
    except ValueError as e:
        assert "budget" in str(e)
    else:
        raise AssertionError("oversized HT payload was not refused")


def test_real_802_11n_frames():
    """Real HT PPDUs (testdata/wifi/real_ht_mcs1_5240): 11 x MCS1 BCC broadcast by one station on
    channel 48. With HT pilots that rotate per symbol (802.11-2020 Eq. 19-58) 8 pass FCS; with the fixed
    [1, 1, 1, -1] pattern the decoder used before -- which our own transmitter shared, so the synthetic
    tests could not see it -- none did."""
    import json
    import os
    from _data import DATA
    meta = json.load(open(os.path.join(DATA, "wifi", "real_ht_mcs1_5240.json")))
    x = np.fromfile(os.path.join(DATA, "wifi", "real_ht_mcs1_5240.cf32"), np.complex64)
    ht = [p for p in W.WifiDetector(meta["fs_hz"], meta["fc_hz"]).run(x) if p.fmt == "HT"]
    gt = meta["ground_truth"]
    assert len(ht) == gt["n_ht"] and all(p.mcs == 1 and p.extra["ldpc"] == 0 for p in ht)
    res = [_decode(p) for p in ht]
    assert sum(r["n_fcs_ok"] for r in res) >= gt["min_fcs_ok"]
    assert {s["mac"]["ta"] for r in res for s in r["subframes"] if s.get("mac")} == {gt["ta"]}
    seeds = [r["scrambler_seed"] for r in res if r["service_ok"]]
    assert np.all(np.diff(seeds) >= 1) and seeds[-1] - seeds[0] <= 2 * len(seeds)   # a sequential counter


def test_constellation_matches_the_standard_not_just_our_transmitter():
    """Pin the QAM bit mapping to IEEE 802.11-2020 Figures 17-10/17-11, written out here by hand.

    This must NOT go through ofdm_tx: the transmitter takes its constellation from the receiver, so
    the two agreeing proves only that they agree. A real 802.11 frame does not share our conventions.
    The bug this guards against (I and Q bits interleaved instead of split) leaves BPSK and QPSK
    perfect and destroys every 16-QAM and 64-QAM frame -- 24/36/48/54 Mb/s and HT MCS 3-7.
    """
    pam2 = {(0, 0): -3, (0, 1): -1, (1, 1): 1, (1, 0): 3}
    pam3 = {(0, 0, 0): -7, (0, 0, 1): -5, (0, 1, 1): -3, (0, 1, 0): -1,
            (1, 1, 0): 1, (1, 1, 1): 3, (1, 0, 1): 5, (1, 0, 0): 7}
    cases = {1: {(0,): -1 + 0j, (1,): 1 + 0j},
             2: {(i, q): (2 * i - 1) + 1j * (2 * q - 1) for i in (0, 1) for q in (0, 1)},
             4: {ib + qb: pam2[ib] + 1j * pam2[qb] for ib in pam2 for qb in pam2},
             6: {ib + qb: pam3[ib] + 1j * pam3[qb] for ib in pam3 for qb in pam3}}
    for nbpsc, table in cases.items():
        pts = np.array(list(table.values())) * P.KMOD[nbpsc]
        llr = P.demap_soft(pts, nbpsc)
        got = (llr > 0).astype(int)
        want = np.array([list(b) for b in table])
        assert np.array_equal(got, want), f"{nbpsc} bits/subcarrier: wrong QAM bit mapping"
        # unit average power, as the standard's K_MOD guarantees
        assert abs(np.mean(np.abs(pts) ** 2) - 1.0) < 1e-9, nbpsc


def test_viterbi_rewrite_matches_the_textbook_form():
    """The fast trellis must decode exactly like the general add-compare-select."""
    def slow(soft, terminated=False):
        soft = np.asarray(soft, float).reshape(-1, 2)
        n = len(soft)
        metric = np.full(64, -1e18)
        metric[0] = 0
        ps_ = np.zeros((n, 64), int)
        pb_ = np.zeros((n, 64), int)
        src = np.repeat(np.arange(64), 2)
        inb = np.tile([0, 1], 64)
        dst = W._NEXT[src, inb]
        for t in range(n):
            cand = metric[src] + W._OUT[src, inb] @ soft[t]
            new = np.full(64, -1e18)
            ps = np.zeros(64, int)
            pb = np.zeros(64, int)
            o = np.argsort(cand)
            new[dst[o]], ps[dst[o]], pb[dst[o]] = cand[o], src[o], inb[o]
            metric, ps_[t], pb_[t] = new, ps, pb
        s = 0 if terminated else int(np.argmax(metric))
        bits = np.zeros(n, int)
        for t in range(n - 1, -1, -1):
            bits[t] = pb_[t, s]
            s = ps_[t, s]
        return bits

    rng = np.random.default_rng(0)
    for trial in range(8):
        bits = rng.integers(0, 2, 160)
        bits[-6:] = 0
        soft = (2 * W.conv_encode(bits) - 1) + rng.normal(scale=0.9, size=2 * len(bits))
        for term in (False, True):
            assert np.array_equal(slow(soft, term), W.viterbi(soft, term)), (trial, term)


def test_non_ht_through_the_same_transmitter():
    for rate in (6, 12, 24, 54):
        pk = _detect(TX.ppdu(MPDU, fmt="NON_HT", rate_mbps=rate))
        nh = [p for p in pk if p.fmt == "NON_HT" and p.lsig_rate_mbps == rate]
        assert nh, rate
        p = nh[0]
        rx = p.extra["_rx"]
        mpdu, ok, _, info, _ = P.decode_nonht_payload(rx[0], rate, p.lsig_length, rx[1], rx[2])
        assert ok, rate
        assert bytes(bytearray(int(v) for v in mpdu))[:-4] == MPDU, rate
        assert info["ssid"] == "HT-lab"


def table():
    big = TX.ampdu([D.data_frame(n_payload=200, seq=i, seed=i) for i in range(8)])
    snrs = (5, 8, 11, 14, 17, 20, 23, 26)
    print(f"{'MCS':>4} {'modulation':>10} {'Mb/s':>6}" + "".join(f"{s:>5}dB" for s in snrs))
    for mcs in range(8):
        row = []
        for snr in snrs:
            ok = 0
            for t in range(4):
                p = [q for q in _detect(TX.ppdu(big, mcs=mcs, aggregation=True), snr_db=snr, seed=t)
                     if q.fmt == "HT" and q.mcs == mcs]
                try:
                    ok += bool(p) and _decode(p[0])["n_fcs_ok"] == 8
                except Exception:
                    pass
            row.append(ok / 4)
        nb, _ = P.HT_MCS[mcs]
        print(f"{mcs:>4} {P.HT_MOD[nb]:>10} {P.ht_rate_mbps(mcs):>6}" + "".join(f"{v:>7.2f}" for v in row))


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("PASS", name)
    print("\nall 8 subframes of a 1952-byte A-MPDU byte-exact, vs in-band SNR (20 MHz):")
    table()
