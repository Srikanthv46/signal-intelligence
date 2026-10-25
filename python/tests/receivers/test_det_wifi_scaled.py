#!/usr/bin/env python3
"""
test_det_wifi_scaled.py -- 802.11 at non-standard clocks (receivers/wifi/ofdm_scaled.py), 802.11n LDPC and
STBC (ldpc.py, payload.decode_ht_payload), and transmitter attribution (links.py).

  1. LDPC: all 12 codes encode/decode; the PPDU layout (shortening, puncturing, repetition, codeword
     count) round-trips for every PSDU length class.
  2. HT payloads byte-exact for BCC/LDPC x with/without STBC (two transmit antennas, standard cyclic
     shifts), MCS 0-7, A-MPDUs; a 40 ppm sampling-clock offset on a 7 kB PPDU.
  3. Synthetic links: a stock 802.11 PPDU stream played at 5, 10 and 11 MHz clocks, at arbitrary
     frequencies, mirrored or not, with a data station and an ACK station (10 dB apart, sequential
     scrambler counters): the numerology, orientation and centre are found, every MPDU decodes, and
     every PPDU is attributed to the right station.
  4. No link from a standard 20 MHz Wi-Fi channel (left to the standard scan), LTE-like OFDM or noise.
  5. Attribution rules on hand-made records: responses, exact counter gaps, two-sided fits, and no
     attribution from a random (non-sequential) seed generator.
  6. Real captures: a Microhard 11 MHz spectrum-inverted 802.11n link (MCS7 LDPC STBC; CRC-valid HT
     MPDUs; two radios told apart) and a DJI 5 MHz 802.11a/g link (testdata/README.md).

pytest collects the test_* functions; `python test_det_wifi_scaled.py` runs them in order.
"""
from __future__ import annotations

import os

import numpy as np

from spectrum_intel.core.emission import load_sigmf
from spectrum_intel.receivers.wifi import dsss_tx as D
from spectrum_intel.receivers.wifi import ldpc as LD
from spectrum_intel.receivers.wifi import links
from spectrum_intel.receivers.wifi import ofdm_rx as W
from spectrum_intel.receivers.wifi import ofdm_scaled as S
from spectrum_intel.receivers.wifi import ofdm_tx as TX
from spectrum_intel.receivers.wifi import payload as P
from _data import DATA

FS20, FC20 = 20e6, 2437e6
CHAN = dict(name="2g6", ch=6, freq_hz=FC20, offset_hz=0.0, usable_half_hz=10e6, partial=False)
STA_A, STA_B = "02:11:22:33:44:55", "02:aa:bb:cc:dd:ee"


# ------------------------------------------------------------------------------------------- helpers
def _ack(ra):
    return bytes([0xD4, 0x00, 0x00, 0x00]) + bytes.fromhex(ra.replace(":", ""))


def _data(seq, n=300, ta=STA_A, ra=STA_B):
    m = bytearray(D.data_frame(bssid=ra, sta=ta, n_payload=n, seq=seq, seed=seq))
    m[4:10], m[10:16], m[16:22] = bytes.fromhex(ra.replace(":", "")), bytes.fromhex(ta.replace(":", "")), \
        bytes.fromhex(ra.replace(":", ""))
    return bytes(m)


def _decode_ht(p):
    rx, e = p.extra["_rx"], p.extra
    return P.decode_ht_payload(rx[0], p.mcs, e["ht_length"], bool(e.get("short_gi")), bool(e.get("aggregation")),
                               n_ess=e.get("n_ess", 0), chan_offset_hz=rx[1], usable_half_hz=rx[2],
                               stbc=e.get("stbc", 0), ldpc=e.get("ldpc", 0), smoothing=e.get("smoothing", 1))


def _one_ht_20(samples, snr_db=32, seed=0, clock_ppm=0.0):
    if clock_ppm:
        samples = TX.timescale(samples, 20e6 * (1 + clock_ppm * 1e-6), 20e6)
    n = int(len(samples) * 1.6) + 8000
    x = TX.scene([dict(samples=samples, f_hz=FC20, t_s=1000 / FS20, snr_db=snr_db)], n, FS20, FC20, seed=seed)
    ht = [p for p in W.WifiDetector(FS20, FC20, channels=[CHAN]).run(x) if p.fmt == "HT"]
    assert len(ht) == 1, len(ht)
    return ht[0]


def link_scene(clock_hz, inverted, offset_hz, fs=61.44e6, n_frames=8, snr_db=26, ht=None, seed=0,
               n_ms=30.0):
    """Station A sends data PPDUs (non-HT 12 Mb/s, or HT with options `ht`), station B answers each with
    an ACK one SIFS later, 10 dB stronger. Both scrambler counters step by one per PPDU. The 20 MSps
    baseband is played at `clock_hz`, optionally mirrored, and placed at `offset_hz` in a capture at fs."""
    rng = np.random.default_rng(seed)
    v = []                                                   # (start sample at 20 MSps, samples, gain)
    t = 400
    for i in range(n_frames):
        payload = _data(100 + i)
        if ht:
            a = TX.ppdu(payload, scramble_seed=10 + i, **ht)
        else:
            a = TX.ppdu(payload, fmt="NON_HT", rate_mbps=12, scramble_seed=10 + i)
        b = TX.ppdu(_ack(STA_A), fmt="NON_HT", rate_mbps=6, scramble_seed=60 + i)
        v.append((t, a, 1.0))
        t += len(a) + 320                                    # SIFS = 16 us at 20 MSps
        v.append((t, b, np.sqrt(10)))
        t += len(b) + int(rng.uniform(2000, 6000))
    base = np.zeros(t + 400, complex)
    for s, x, g in v:
        base[s:s + len(x)] += g * x / np.sqrt(np.mean(np.abs(x) ** 2))
    if inverted:
        base = np.conj(base)
    y = TX.timescale(base, clock_hz, fs)
    n = int(n_ms * 1e-3 * fs)
    assert len(y) < n, (len(y), n)
    noise = (rng.normal(size=n) + 1j * rng.normal(size=n)) / np.sqrt(2)
    # SNR of station A in the occupied band (clock_hz wide)
    amp = np.sqrt(10 ** (snr_db / 10) * clock_hz / fs)
    i0 = int(rng.uniform(0.5e-3, 1.5e-3) * fs)
    ph = np.exp(2j * np.pi * offset_hz / fs * np.arange(len(y)))
    noise[i0:i0 + len(y)] += amp * y * ph
    return noise.astype(np.complex64)


def _recs(link, fs):
    """ofdm_scaled packets -> the records links.attribute reads (what the plugin builds)."""
    from types import SimpleNamespace
    from spectrum_intel.receivers.wifi.plugin_ofdm import _records
    return _records(link.packets, SimpleNamespace(fs=fs), link.clock_hz)


# ------------------------------------------------------------------------------------------- 1. LDPC
def test_ldpc_codes_roundtrip():
    rng = np.random.default_rng(1)
    for L in (648, 1296, 1944):
        for rate in ((1, 2), (2, 3), (3, 4), (5, 6)):
            cd = LD.code(L, rate)
            u = rng.integers(0, 2, cd.K).astype(np.uint8)
            cw = LD.encode(u, L, rate)
            assert cd.syndrome_ok(cw[None])[0], (L, rate)
            llr = (1 - 2.0 * cw) * 3 + rng.normal(0, 1.3, L)
            hard, ok, _ = LD.decode(llr, L, rate)
            assert ok and np.array_equal(hard, cw), (L, rate)


def test_ldpc_ppdu_layout_roundtrip():
    """Every branch of the codeword-count table, with and without STBC, decodes back exactly."""
    rng = np.random.default_rng(2)
    seen = set()
    for length in (10, 40, 90, 150, 220, 300, 1000, 4000):
        for nbpsc, rate in ((1, (1, 2)), (4, (3, 4)), (6, (5, 6))):
            for m_stbc in (1, 2):
                Pm = LD.params(length, 52 * nbpsc, rate, m_stbc)
                assert Pm["navbits"] == Pm["nsym"] * 52 * nbpsc and Pm["nsym"] % m_stbc == 0
                seen.add((Pm["ncw"] > 1, Pm["L"], Pm["npunc"] > 0, Pm["nrep"] > 0))
                bits = rng.integers(0, 2, Pm["npld"]).astype(np.uint8)
                c = LD.encode_stream(bits, Pm, rate)
                assert len(c) == Pm["navbits"]
                info, ok, _ = LD.decode_stream((1 - 2.0 * c) * 4, Pm, rate)
                assert ok.all() and np.array_equal(info[:Pm["npld"]], bits), (length, nbpsc, m_stbc)
    assert {L for _, L, _, _ in seen} == {648, 1296, 1944}
    assert any(p for _, _, p, _ in seen) and any(r for _, _, _, r in seen) and any(m for m, _, _, _ in seen)


# ------------------------------------------------------------------------------------------- 2. HT payloads
def test_ht_bcc_ldpc_stbc_byte_exact():
    mpdu = _data(7, n=200)
    for ldpc in (False, True):
        for stbc in (False, True):
            for mcs in (0, 3, 5, 7):
                p = _one_ht_20(TX.ppdu(mpdu, mcs=mcs, ldpc=ldpc, stbc=stbc))
                assert (p.extra["ldpc"], p.extra["stbc"]) == (int(ldpc), int(stbc)), (mcs, ldpc, stbc)
                assert P.ht_decodable(p.mcs, p.extra["ldpc"], p.extra["stbc"])[0]
                r = _decode_ht(p)
                assert r["n_fcs_ok"] == 1 and bytes.fromhex(r["subframes"][0]["mpdu_hex"]) == mpdu, \
                    (mcs, ldpc, stbc)
                assert r["scrambler_seed"] == 0x5D and r["service_ok"]
    parts = [_data(200 + i, n=150) for i in range(5)]
    p = _one_ht_20(TX.ppdu(TX.ampdu(parts), mcs=7, aggregation=True, ldpc=True, stbc=True))
    r = _decode_ht(p)
    assert r["n_fcs_ok"] == 5 and [bytes.fromhex(s["mpdu_hex"]) for s in r["subframes"]] == parts
    assert r["n_codewords_ok"] == r["n_codewords"] > 1


def test_ht_long_ppdu_with_clock_offset():
    """A 7 kB MCS7 LDPC STBC A-MPDU (0.9 ms) with a 40 ppm sampling-clock offset (0.7 samples of drift
    by the end): the pilot phase-slope tracking keeps the edge subcarriers aligned. The real Microhard
    link drifts 10 ppm against the Pluto."""
    parts = [_data(300 + i, n=1400) for i in range(5)]
    p = _one_ht_20(TX.ppdu(TX.ampdu(parts), mcs=7, aggregation=True, ldpc=True, stbc=True), snr_db=30,
                   clock_ppm=40)
    assert _decode_ht(p)["n_fcs_ok"] == 5


def test_ht_pilots_rotate():
    """802.11-2020 Eq. 19-58: pilot j of HT data symbol n is PSI[(n + j) mod 4] * p_{n+3}."""
    for n in range(8):
        want = np.array([[1, 1, 1, -1][(n + j) % 4] for j in range(4)]) * W.PILOT_POL[(n + 3) % 127]
        assert np.array_equal(P.ht_pilots(n, 1)[0], want)
    assert not np.array_equal(P.ht_pilots(1, 1)[0], P.ht_pilots(0, 1)[0] * W.PILOT_POL[4] / W.PILOT_POL[3])


# ------------------------------------------------------------------------------------------- 3. synthetic links
def _check_link(x, fs, fc, clock_hz, inverted, offset_hz, n_frames=8, ht=None):
    L = S.scan(x, fs, fc)
    assert len(L) == 1, [(ln.center_hz, ln.clock_hz) for ln in L]
    ln = L[0]
    assert ln.clock_hz == clock_hz and ln.inverted == inverted, (ln.clock_hz, ln.inverted)
    assert abs(ln.center_hz - (fc + offset_hz)) < 0.02 * clock_hz, ln.center_hz
    recs = _recs(ln, fs)
    assert len(recs) == 2 * n_frames, len(recs)
    assert sum(1 for r in recs if r.get("fcs_ok")) == 2 * n_frames
    tx = {t["ta"]: t for t in links.attribute(recs, sifs_us=16.0 * 20e6 / clock_hz)}
    assert set(tx) == {STA_A, STA_B}
    assert all(r.get("ta") == (STA_B if r["mac"]["subtype"] == "ack" else STA_A) for r in recs)
    assert tx[STA_A]["scrambler_counter"] == tx[STA_B]["scrambler_counter"] == "sequential"
    assert tx[STA_B]["power_db_median"] - tx[STA_A]["power_db_median"] > 7
    if ht:
        assert all(r["fmt"] == "HT" for r in recs if r["mac"]["subtype"] != "ack")
    return ln, recs


def test_quarter_clock_5mhz_link():
    fs, fc = 61.44e6, 2420e6
    _check_link(link_scene(5e6, False, 12e6, fs), fs, fc, 5e6, False, 12e6)


def test_half_clock_10mhz_mirrored_link():
    fs, fc = 61.44e6, 2060e6
    _check_link(link_scene(10e6, True, -14e6, fs, seed=3), fs, fc, 10e6, True, -14e6)


def test_11mhz_mirrored_ht_ldpc_stbc_link():
    """The Microhard shape: 11 MHz clock, mirrored spectrum, MCS7 LDPC STBC data, legacy ACKs."""
    fs, fc = 61.44e6, 2060e6
    ht = dict(mcs=7, ldpc=True, stbc=True)
    _check_link(link_scene(11e6, True, 9e6, fs, ht=ht, seed=5, snr_db=30), fs, fc, 11e6, True, 9e6, ht=ht)


# ------------------------------------------------------------------------------------------- 4. no false links
def test_no_link_from_standard_wifi_lte_or_noise():
    fs, fc = 61.44e6, 2437e6
    rng = np.random.default_rng(9)
    n = int(30e-3 * fs)
    # standard 20 MHz channel 6: the regular scan owns it
    pp = [dict(samples=TX.ppdu(_data(i), fmt="NON_HT", rate_mbps=24), f_hz=2437e6, t_s=1e-3 + 3e-3 * i, snr_db=25)
          for i in range(8)]
    x = TX.scene(pp, n, fs, fc, seed=1)
    assert S.scan(x, fs, fc) == []
    assert len(W.WifiDetector(fs, fc).run(x)) >= 8
    # LTE-like OFDM bursts: 15 kHz subcarriers, 1024-point FFT at 15.36 MSps, normal CP (72)
    lte = []
    for _ in range(14 * 6):
        X = np.zeros(1024, complex)
        X[np.r_[1:301, 724:1024]] = (rng.choice([-1, 1], 600) + 1j * rng.choice([-1, 1], 600)) / np.sqrt(2)
        s = np.fft.ifft(X) * 32
        lte.append(np.concatenate([s[-72:], s]))
    y = TX.timescale(np.concatenate(lte), 15.36e6, fs)
    x2 = (rng.normal(size=n) + 1j * rng.normal(size=n)).astype(np.complex64)
    for k in range(4):
        i0 = int((2 + 7 * k) * 1e-3 * fs)
        x2[i0:i0 + len(y)] += 4 * y[:n - i0] * np.exp(2j * np.pi * -8e6 / fs * np.arange(min(len(y), n - i0)))
    assert S.scan(x2, fs, fc) == []
    assert S.scan((rng.normal(size=n) + 1j * rng.normal(size=n)).astype(np.complex64), fs, fc) == []


# ------------------------------------------------------------------------------------------- 5. attribution rules
def _rec(t_ms, dur_us, sub=None, ta=None, ra=None, seed=None, fmt="NON_HT", power=40.0):
    mac = {"subtype": sub, **({"ta": ta} if ta else {}), **({"ra": ra} if ra else {})} if sub else None
    r = {"t0": t_ms * 1e-3, "t1": t_ms * 1e-3 + dur_us * 1e-6, "fmt": fmt, "power_db": power, "cfo_khz": 1.0}
    if mac:
        r["mac"] = mac
    if seed is not None:
        r["scrambler_seed"], r["seed_ok"] = seed, True
    return r


def test_attribution_rules():
    A, B = "00:00:00:00:00:0a", "00:00:00:00:00:0b"
    recs = [
        _rec(0.0, 100, "rts", A, B, seed=125),
        _rec(0.120, 30, "cts", None, A, seed=40, power=55),           # 0.12: answers the RTS -> sent by B
        _rec(0.200, 900, None, seed=126, fmt="HT"),                    # ACKed to A below -> sent by A
        _rec(1.120, 30, "ack", None, A, seed=41, power=55),            # its seed fits B's counter 40 .. 42
        _rec(2.000, 900, None, fmt="HT"),                              # 2-stream: no seed, no response
        _rec(3.000, 100, "rts", A, B, seed=1),                         # A: 126 -> 1 = two steps, one PPDU between
        _rec(3.120, 30, "cts", None, A, seed=42, power=55),
        _rec(4.000, 900, None, seed=2, fmt="HT"),                      # fits A's counter on both sides
        _rec(5.000, 100, "rts", A, B, seed=3),
        _rec(5.120, 30, "cts", None, A, seed=43, power=55),
    ]
    rows = {r["ta"]: r for r in links.attribute(recs)}
    got = [(r.get("ta"), r.get("ta_evidence")) for r in sorted(recs, key=lambda r: r["t0"])]
    assert got[1] == (B, "response") and got[2] == (A, "response") and got[6] == (B, "response")
    assert got[3] == (B, "scrambler counter")        # two-sided fit on B's counter
    assert got[4] == (A, "scrambler counter")        # exact gap: A's counter stepped once in between
    assert got[7] == (A, "scrambler counter")
    assert rows[A]["scrambler_counter"] == rows[B]["scrambler_counter"] == "sequential"
    assert rows[B]["power_db_median"] - rows[A]["power_db_median"] == 15
    # a station that picks seeds at random is not treated as a counter
    rng = np.random.default_rng(0)
    rnd = []
    for i in range(12):
        rnd.append(_rec(i * 1.0, 100, "rts", A, B, seed=int(rng.integers(1, 128))))
        rnd.append(_rec(i * 1.0 + 0.5, 300, None, seed=int(rng.integers(1, 128)), fmt="HT"))
    rows = {r["ta"]: r for r in links.attribute(rnd)}
    assert rows[A]["scrambler_counter"] == "not sequential"
    assert not any(r.get("ta_evidence") == "scrambler counter" for r in rnd)


# ------------------------------------------------------------------------------------------- 6. real captures
def _real(name):
    x, fs, fc, _ = load_sigmf(os.path.join(DATA, name, name + ".sigmf-meta"))
    return np.asarray(x, np.complex64), fs, fc


def test_real_microhard_11mhz_inverted_link():
    x, fs, fc = _real("D_20260211-161341")
    L = S.scan(x, fs, fc)
    assert len(L) == 1
    ln = L[0]
    assert ln.clock_hz == 11e6 and ln.inverted and abs(ln.center_hz - 2069.02e6) < 50e3
    recs = _recs(ln, fs)
    assert len(recs) >= 45 and sum(1 for r in recs if r.get("fcs_ok")) >= 25
    ht = [r for r in recs if r["fmt"] == "HT"]
    assert len(ht) >= 20 and sum(1 for r in ht if r.get("fcs_ok")) >= 2      # MCS7 LDPC STBC, CRC-valid
    assert all(m["ta"] == "00:0f:92:fd:58:12" for r in ht for m in r.get("macs", []))
    tx = {t["ta"]: t for t in links.attribute(recs, sifs_us=16.0 * 20 / 11)}
    assert set(tx) == {"00:0f:92:fd:58:12", "00:0f:92:fd:3e:bf"}
    assert all(t["vendor"] == "Microhard Systems Inc." and t["scrambler_counter"] == "sequential"
               for t in tx.values())
    assert sum(1 for r in recs if r.get("ta")) >= 45
    assert all(r["ta"] == "00:0f:92:fd:58:12" for r in ht if r.get("ta"))     # every HT PPDU: the data radio
    assert tx["00:0f:92:fd:3e:bf"]["power_db_median"] - tx["00:0f:92:fd:58:12"]["power_db_median"] > 10


def test_real_dji_5mhz_link():
    x, fs, fc = _real("D_20251119-183758")
    L = S.scan(x, fs, fc)
    assert len(L) == 1
    ln = L[0]
    assert ln.clock_hz == 5e6 and not ln.inverted and abs(ln.center_hz - 2432.03e6) < 50e3
    recs = _recs(ln, fs)
    assert len(recs) >= 11 and sum(1 for r in recs if r.get("fcs_ok")) >= 7
    tx = links.attribute(recs, sifs_us=64.0)
    assert [t["ta"] for t in tx] == ["34:d2:62:53:20:37"] and tx[0]["vendor"] == "Sz Dji Technology Co.,Ltd"
    assert all(r.get("ta") == "34:d2:62:53:20:37" for r in recs)
    assert all(m.get("encryption", {}).get("cipher") == "WEP" for r in recs for m in r.get("macs", [])
               if m.get("subtype") == "data")


if __name__ == "__main__":
    import time
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            t = time.time()
            fn()
            print(f"PASS {name}  ({time.time() - t:.1f} s)")
