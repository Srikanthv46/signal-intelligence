#!/usr/bin/env python3
"""
test_det_wifi_dsss.py -- 802.11b DSSS/CCK receiver (receivers/wifi/dsss.py) and MAC dissector.

  1. Every rate (1/2/5.5/11 Mb/s) and preamble (long/short): FCS-valid and byte-exact MPDUs.
  2. Sensitivity: 1 Mb/s at 0 dB and 11 Mb/s CCK at 6 dB in-band SNR (22 MHz).
  3. Impairments: CFO +-120 kHz, unlocked chip clock -25 ppm on 1500-byte CCK, indoor multipath.
  4. FCS-guided correction never invents bytes (every corrected frame equals what was sent).
  5. No CRC-valid header on noise; two BSSs on ch1 and ch6 at once, each decoded once.
  6. MAC dissector on a synthetic beacon (SSID, channel, RSN, rates) and data frame.
  7. Real capture (Phantom 4, really 2430 MHz): >= 13 FCS-valid frames with the known SSIDs.

pytest collects the test_* functions; `python test_det_wifi_dsss.py` prints a Pd-vs-SNR table too.
"""
from __future__ import annotations

import os
import sys

import numpy as np

from spectrum_intel.core.bands_india import WIFI_24
from spectrum_intel.core.dsp import IQContext
from spectrum_intel.receivers.wifi import dsss, mac80211
from spectrum_intel.receivers.wifi import dsss_tx as TX

sys.path.insert(0, os.path.dirname(__file__))
from _data import DATA  # noqa: E402

FS, FC = 61.44e6, 2437e6
P4 = os.path.join(DATA, "DJI_Phantom_4_20251119-182450/D_20251119-182450.sigmf-meta")


def _train(rate, pre, snr, n=4, ch=6, size=200, seed=0, **kw):
    pk, t = [], 0.4e-3
    for i in range(n):
        mp = TX.data_frame(n_payload=size, seq=i, seed=seed * 100 + i) if size else TX.beacon(ssid=f"t{i}", channel=ch)
        pk.append(dict(mpdu=mp, rate=rate, preamble=pre, f_hz=WIFI_24[ch], t_s=t, snr_db=snr, **kw))
        t += len(TX.ppdu_chips(mp, rate, pre)) / TX.CHIP_RATE + 250e-6
    return pk, int((t + 1e-3) * FS)


def _run(pk, n, seed=1):
    out = dsss.scan(IQContext(TX.scene(pk, n, FS, FC, seed=seed), FS, FC))
    got = []
    for p in pk:
        m = [r for r in out if r.get("fcs_ok") and abs(r["start_s"] - p["t_s"]) < 20e-6]
        got.append(m[0] if m else None)
    return got, out


def _exact(r, p):
    return r is not None and bytes.fromhex(r["mpdu_hex"]) == bytes(p["mpdu"])[:512]


def test_every_rate_and_preamble():
    for rate, pre in ((1.0, "long"), (2.0, "long"), (2.0, "short"), (5.5, "long"), (5.5, "short"),
                      (11.0, "long"), (11.0, "short")):
        pk, n = _train(rate, pre, 15, n=3, cfo_hz=40e3)
        got, _ = _run(pk, n)
        assert all(_exact(r, p) for r, p in zip(got, pk)), (rate, pre)
        assert all(r["rate_mbps"] == rate and r["preamble"] == pre and r["channel"] == 6 for r in got)


def test_sensitivity():
    for rate, pre, snr, need in ((1.0, "long", 0, 0.9), (2.0, "short", 0, 0.75), (11.0, "short", 6, 0.75)):
        pk, n = _train(rate, pre, snr, n=8, cfo_hz=60e3, seed=3)
        got, _ = _run(pk, n, seed=5)
        assert np.mean([_exact(r, p) for r, p in zip(got, pk)]) >= need, (rate, pre, snr)


def test_impairments():
    cases = [dict(rate=11.0, pre="short", snr=15, cfo_hz=-120e3), dict(rate=1.0, pre="long", snr=6, cfo_hz=120e3),
             dict(rate=11.0, pre="long", snr=15, cfo_hz=10e3, ppm=-25, size=1500),
             dict(rate=11.0, pre="short", snr=20, taps=[(90, 0.7 * np.exp(0.3j)), (200, 0.4 * np.exp(-2j)), (350, 0.2)]),
             dict(rate=5.5, pre="short", snr=20, taps=[(50, 0.5 * np.exp(1j)), (150, 0.25 * np.exp(2.5j))])]
    for c in cases:
        c = dict(c)
        rate, pre, snr = c.pop("rate"), c.pop("pre"), c.pop("snr")
        pk, n = _train(rate, pre, snr, n=3, **c)
        got, _ = _run(pk, n)
        assert all(_exact(r, p) for r, p in zip(got, pk)), (rate, pre, c)


def test_fcs_correction_is_exact():
    fixed = 0
    for rate, pre, snr in ((1.0, "long", -3), (2.0, "short", -2), (5.5, "short", 0), (11.0, "long", 3)):
        pk, n = _train(rate, pre, snr, n=8, seed=7)
        _, out = _run(pk, n, seed=11)
        for r in out:
            if r.get("symbols_corrected"):
                p = min(pk, key=lambda p: abs(p["t_s"] - r["start_s"]))
                assert _exact(r, p)
                fixed += 1
    assert fixed >= 3


def test_noise_and_two_channels():
    x = TX.scene([], int(10e-3 * FS), FS, FC, seed=4)
    assert not [r for r in dsss.scan(IQContext(x, FS, FC)) if r["plcp_crc_ok"]]
    a = TX.beacon("02:00:00:00:00:01", "on-ch1", 1)
    b = TX.beacon("02:00:00:00:00:06", "on-ch6", 6)
    pk = [dict(mpdu=a, f_hz=WIFI_24[1], t_s=1e-3, snr_db=12), dict(mpdu=b, f_hz=WIFI_24[6], t_s=1.2e-3, snr_db=12)]
    out = dsss.scan(IQContext(TX.scene(pk, int(5e-3 * FS), FS, 2425e6), FS, 2425e6))
    ok = [r for r in out if r.get("fcs_ok")]
    assert sorted((r["channel"], r["mac"]["ssid"]) for r in ok) == [(1, "on-ch1"), (6, "on-ch6")]
    assert all(r["channel_source"].startswith("DS parameter") for r in ok)


def test_narrowband_collision_excision():
    """A BLE-like GFSK burst 15 dB over a beacon for 400 us: recovered by STFT excision."""
    ok = 0
    for trial in range(4):
        mp = TX.beacon(ssid=f"victim{trial}", channel=6)
        n = int(4e-3 * FS)
        x = TX.scene([dict(mpdu=mp, f_hz=WIFI_24[6], t_s=0.5e-3, snr_db=8, cfo_hz=20e3)], n, FS, FC,
                     seed=trial).astype(complex)
        rng = np.random.default_rng(trial)
        i0, m = int(1.2e-3 * FS), int(400e-6 * FS)
        fdev = 250e3 * np.repeat(rng.choice([-1, 1], 402), int(FS / 1e6) + 1)[:m]
        ph = 2 * np.pi * np.cumsum(fdev + (WIFI_24[6] + 2e6 - FC + (trial - 1.5) * 1.4e6)) / FS
        x[i0:i0 + m] += np.sqrt(10 ** 0.8 * 22e6 / FS * 10 ** 1.5) * np.exp(1j * ph)
        out = dsss.scan(IQContext(x.astype(np.complex64), FS, FC))
        ok += any(r.get("fcs_ok") and r["mac"]["ssid"] == f"victim{trial}" for r in out)
    assert ok >= 3, ok


def test_mac_dissector():
    m = mac80211.parse(TX.beacon("24:a3:f0:a5:bd:ab", "Lab AP", 6, seq=42, tsf_us=90061_000_000))
    assert (m["subtype"], m["ssid"], m["ds_channel"], m["seq"]) == ("beacon", "Lab AP", 6, 42)
    assert m["security"].startswith("WPA2-Personal (PSK) CCMP-128")
    assert m["rates_mbps"] == [1.0, 2.0, 5.5, 11.0] and m["phy"] == "802.11b"
    assert m["ap_uptime"] == "1d 01:01:01" and m["bssid_vendor"].startswith("TeraNXT")
    d = mac80211.parse(TX.data_frame())
    assert (d["subtype"], d["direction"], d["qos"]["tid"]) == ("qos-data", "AP -> STA", 0)
    assert d["encryption"] == {"cipher": "CCMP/GCMP", "key_id": 0, "pn": 8} and "summary" in d


def test_real_phantom4():
    if not os.path.exists(P4):
        return
    from spectrum_intel.io.sigmf import load
    cap = load(P4)
    out = dsss.scan(IQContext(cap.iq, cap.fs, 2430e6))           # label is +20 MHz wrong (testdata/README)
    ok = [r for r in out if r.get("fcs_ok")]
    ssids = {r["mac"].get("ssid") for r in ok}
    assert len(ok) >= 13, len(ok)
    assert {"Sharma Brothers_2G", "Hindusena", "FTTH", "Futurestride", "SFS consultancy _2G"} <= ssids, ssids


def table():
    print(f"{'rate':>5} {'pre':6}" + "".join(f"{s:>6}dB" for s in (-4, -2, 0, 2, 4, 6, 8)))
    for rate, pre in ((1.0, "long"), (2.0, "short"), (5.5, "short"), (11.0, "short")):
        row = []
        for snr in (-4, -2, 0, 2, 4, 6, 8):
            pk, n = _train(rate, pre, snr, n=10, cfo_hz=50e3, seed=snr + 10)
            got, _ = _run(pk, n, seed=snr + 20)
            row.append(np.mean([_exact(r, p) for r, p in zip(got, pk)]))
        print(f"{rate:5} {pre:6}" + "".join(f"{v:8.2f}" for v in row))


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("PASS", name)
    print("\nFCS-valid, byte-exact fraction vs in-band SNR (22 MHz), 230-byte MPDUs, CFO 50 kHz:")
    table()
