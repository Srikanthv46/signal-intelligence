"""HE-SIG-A (802.11ax, PLAN W2): field round-trips on the reference preamble (ofdm_tx.he_preamble), sensitivity,
and the real HE PPDUs of testdata/live/LIVE_WIFI_5755MHz (every one passes CRC-4)."""
import json
from pathlib import Path

import numpy as np
import pytest

from spectrum_intel.receivers.wifi import ofdm_rx as W
from spectrum_intel.receivers.wifi import ofdm_tx as T

REPO = Path(__file__).resolve().parents[3]
FS, FC = 61.44e6, 2437e6


def _rx(sig, snr=20, cfo=0.0, seed=0, n=60000):
    x = T.scene([dict(samples=sig, f_hz=FC + cfo, t_s=1e-4, snr_db=snr)], n=n, fs=FS, fc=FC, seed=seed)
    return [p for p in W.WifiDetector(FS, FC).run(x) if p.fmt == "HE"]


@pytest.mark.parametrize("fmt,kw,want", [
    ("SU", dict(mcs=9, bss_color=42, uplink=True, nsts=2, txop=2 * 10), dict(mcs=9, nsts=2, txop_us=80, bw=20)),
    ("SU", dict(mcs=0, bss_color=1, bw=1, txop=1 + 2 * 2, ldpc=0), dict(mcs=0, nsts=1, txop_us=768, bw=40)),
    ("MU", dict(bss_color=5, bw=2, txop=127), dict(txop_us=None, bw=80)),
    ("TB", dict(bss_color=63, bw=3), dict(txop_us=None, bw=160)),
])
def test_he_sig_a_round_trip(fmt, kw, want):
    pk = _rx(T.he_preamble(T.he_sig_a_bits(fmt, **kw), mu=fmt == "MU"), cfo=23e3)
    assert len(pk) == 1 and pk[0].sig_crc_ok is True
    e = pk[0].extra
    assert e["he_format"] == fmt and e["bss_color"] == kw["bss_color"]
    assert e["uplink"] == (True if fmt == "TB" else kw.get("uplink", False))
    for k, v in want.items():
        assert e.get(k) == v, k
    assert pk[0].bandwidth_mhz == want["bw"]


def test_he_sig_a_crc_rejects():
    b = T.he_sig_a_bits("SU", bss_color=5)
    assert W.parse_he_sig_a(b, 40) is not None
    for i in (0, 13, 30, 44):
        bad = list(b)
        bad[i] ^= 1
        assert W.parse_he_sig_a(bad, 40) is None
    assert W.he_txop_us(127) is None and W.he_txop_us(2 * 7) == 56 and W.he_txop_us(1 + 2 * 1) == 640


def test_he_sig_a_sensitivity():
    """Bench 2026-10-25 (20 trials each): 6 dB 15 %, 8 dB 100 % in-band SNR."""
    ok = 0
    for t in range(10):
        pk = _rx(T.he_preamble(T.he_sig_a_bits("SU", bss_color=t + 1), seed=t), snr=8, cfo=37e3, seed=t)
        ok += any(p.sig_crc_ok and p.extra["bss_color"] == t + 1 for p in pk)
    assert ok >= 9


def test_real_he_ppdus_5755():
    d = REPO / "testdata/live/LIVE_WIFI_5755MHz"
    meta = json.loads((d / "LIVE_WIFI_5755MHz.sigmf-meta").read_text())
    x = np.fromfile(d / "LIVE_WIFI_5755MHz.sigmf-data", "<i2").astype(np.float32).view(np.complex64)
    pk = W.WifiDetector(meta["global"]["core:sample_rate"], meta["captures"][0]["core:frequency"]).run(x)
    he = [p for p in pk if p.fmt == "HE"]
    assert len(he) >= 40 and all(p.sig_crc_ok for p in he)
    down = [p for p in he if not p.extra["uplink"]]
    assert down and all(p.extra["bss_color"] == 14 and p.bandwidth_mhz == 80 for p in down)
    assert {p.extra["he_format"] for p in he} == {"SU"}
