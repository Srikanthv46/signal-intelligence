"""ExpressLRS sub-GHz SYNC packets (research, PLAN R2c): implicit-header LoRa payload + the self-verifying CRC-14
(the packet's own UID4 / UID5 seed it). Implicit-header decoding is checked against gr-lora_sdr
(testdata/lora/tx_implicit_*, tests/receivers/gen_lora_implicit.py)."""
import json
from pathlib import Path

import numpy as np
import pytest

from spectrum_intel.receivers.rc_link import detector as D, elrs as E, lora_phy as P, synth as S

LORA = Path(__file__).resolve().parents[2] / "testdata" / "lora"
UID = bytes([0x9A, 0x01, 0x5C, 0x7E, 0x33, 0xC4])


def test_implicit_header_vs_gr_lora_sdr():
    walks = []
    orig = D.lora_emission
    D.lora_emission = lambda w, fc, fs: (walks.append(w), orig(w, fc, fs))[1]
    try:
        for m in json.load(open(LORA / "tx_implicit.json")):
            x = np.fromfile(LORA / f"tx_implicit_{m['name']}.cf32", np.complex64)
            x = np.r_[np.zeros(20000, np.complex64), x, np.zeros(20000, np.complex64)]
            walks.clear()
            D.detect_lora(np.fft.fft(x), m["fs"], 915e6, bws=[(m["bw"], [m["sf"]])])
            n = len(bytes.fromhex(m["payloads_hex"][0]))
            got = [P.decode_implicit(w["all_vals"], m["sf"], m["cr"], n)["payload"].hex()
                   for w in sorted(walks, key=lambda w: w["start_samp"])]
            assert got == m["payloads_hex"], m["name"]
    finally:
        D.lora_emission = orig


@pytest.mark.parametrize("kw", [dict(ver=4), dict(ver=3), dict(ver=4, binding=True), dict(ver=4, ota8=True)])
def test_sync_crc_roundtrip(kw):
    r = E.parse(E.build_sync(12, 77, 3, 2, UID, **kw))
    assert r["crc_ok"] and r["ota_version"] == kw["ver"] and r["binding"] == kw.get("binding", False)
    assert r["sync"]["uid4"] == UID[4] and r["sync"]["uid5"] == UID[5] and r["sync"]["fhss_index"] == 12


def test_corrupt_sync_not_verified():
    p = bytearray(E.build_sync(12, 77, 3, 2, UID)); p[3] ^= 0x10
    assert E.parse(bytes(p))["crc_ok"] is False


@pytest.mark.parametrize("mode,snr_db", [("900_200hz", 10), ("900_100hz", 10), ("900_50hz", 3), ("900_25hz", 0)])
def test_elrs_900_sync_in_capture(mode, snr_db):
    band, sf, bw, cr_den, pre, pl, hop, ival = D.ELRS_LORA[mode]
    rng = np.random.default_rng(sf)
    fs, fc = 7.68e6, 915e6
    x = S.noise_capture(1 << 18, rng)
    pkt = E.build_sync(40, 9, 2, 3, UID)
    iq, _ = S.lora_packet(sf, bw, n_pre=pre, payload=S.lora_implicit_symbols(pkt, sf, cr_den - 4), sync_word=0x12,
                          cfo_hz=4e3, fs=fs, rng=rng)
    f_ch = 903.5e6 + (926.9e6 - 903.5e6) * 20 / 39                      # FCC915 channel 20 (the sync channel)
    y, _ = S.embed(x, fs, fc, [dict(t0_us=2000, f_hz=f_ch, iq=iq, bw_hz=bw)], snr_db)
    em = [e for e in D.detect(y, fs, fc) if (e.decoded or {}).get("protocol") == "ExpressLRS"]
    assert em, [(e.label, e.decoded) for e in D.detect(y, fs, fc)]
    d = em[0].decoded
    assert d["crc_ok"] and d["packet_type"] == "SYNC" and d["sync"]["uid4"] == UID[4] and d["sync"]["fhss_index"] == 40
    assert em[0].label.startswith("elrs")
