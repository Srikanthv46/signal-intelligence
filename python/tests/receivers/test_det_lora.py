#!/usr/bin/env python3
"""LoRa payload decoding (rc_link detector + lora_phy), checked against references
that do not share our conventions:

  * gr-lora_sdr's TRANSMITTER (testdata/lora/tx_*.cf32, made by gen_lora_vectors.py):
    SF7-SF11, CR 4/5-4/8, CRC on/off, LDRO -- every payload must come back byte-exact.
  * a real 434 MHz packet off the Pluto (testdata/lora/real_434MHz_sf9.cf32) whose
    payload gr-lora_sdr's RECEIVER decoded with a valid CRC.
  * CRC-16 and the whitening sequence against their published check values.

Usage: pytest -q tests/receivers/test_det_lora.py   (or run it as a script)
"""
from __future__ import annotations
import json, os

import numpy as np
from scipy import fft as sfft, signal as ss

from spectrum_intel.receivers.rc_link import detector as D, lora_phy as P
from _data import DATA

LD = os.path.join(DATA, "lora")


def _detect(x, fs, fc):
    return [e for e in D.detect_lora(sfft.fft(x.astype(np.complex64), workers=-1), fs, fc)
            if e.decoded and e.decoded.get("complete")]


def test_crc16_and_whitening_check_values():
    assert P.crc16(b"123456789") == 0x31C3                  # CRC-16/XMODEM check value
    assert list(P.whitening_sequence(8)) == [0xFF, 0xFE, 0xFC, 0xF8, 0xF0, 0xE1, 0xC2, 0x85]


def test_gr_lora_sdr_transmitter_vectors():
    rng = np.random.default_rng(7)
    for m in json.load(open(os.path.join(LD, "tx_vectors.json"))):
        x = np.fromfile(os.path.join(LD, f"tx_{m['name']}.cf32"), np.complex64)
        x = np.r_[np.zeros(20000, np.complex64), x]          # the file starts with a frame
        x = ss.resample_poly(x, 2, 1)                        # 250 -> 500 kSps
        fs = 2 * m["fs"]
        for off in (-90e3, 37e3):
            y = x * np.exp(2j * np.pi * off * np.arange(len(x)) / fs)
            y = y + 0.05 * (rng.standard_normal(len(y)) + 1j * rng.standard_normal(len(y)))
            got = _detect(y, fs, 0.0)
            pays = sorted(e.decoded["payload_hex"] for e in got)
            assert pays == sorted(m["payloads_hex"]), (m["name"], off, pays)
            for e in got:
                d = e.decoded
                assert d["crc_ok"] is (True if m["crc"] else None)
                assert d["cr"] == f"4/{4 + m['cr']}" and d["sync_word"] == "0x12"
                assert abs(e.center_mhz * 1e6 - off) < 2e3


def test_real_pluto_packet_matches_gr_lora_sdr():
    meta = json.load(open(os.path.join(LD, "real_434MHz_sf9.json")))
    x = np.fromfile(os.path.join(LD, "real_434MHz_sf9.cf32"), np.complex64)
    got = _detect(x, meta["fs_hz"], meta["fc_hz"])
    assert len(got) == 1
    gt, d = meta["ground_truth"], got[0].decoded
    assert d["payload_hex"] == gt["payload_hex"] and d["crc_ok"] is True
    assert d["sync_word"] == gt["sync_word"] and d["cr"] == gt["cr"]
    assert got[0].label == "lora_sf9_bw125k" and got[0].confidence >= 0.95
    assert abs(got[0].center_mhz - 434.0023) < 0.002


def test_sf5_sx126x_real_packet():
    """Live SF5 / 500 kHz implicit-header packet at 890.03 MHz (bench 2026-10-19: 12/12 packets CRC-valid,
    and synth.lora_sf5_implicit_symbols re-encodes each payload to the exact on-air symbols)."""
    from spectrum_intel.io import sigmf
    cap = sigmf.load(os.path.join(DATA, "live", "LIVE_LORA_SF5_890MHz"))
    got = _detect(cap.iq, cap.fs, cap.fc)
    assert len(got) == 1
    e, d = got[0], got[0].decoded
    assert e.label == "lora_sf5_bw500k" and e.confidence >= 0.95
    assert d["payload_hex"] == "e883000b3e000032000000aedbd71e14e000e0e80300007001002000000000421f0008000200"
    assert d["crc_ok"] is True and d["implicit_header"] and d["length"] == 38 and d["cr"] == "4/5"
    assert d["cr_ambiguous"] is False
    assert d["sync_word"] == "0x12" and e.features["payload_syms"] == 85
    assert abs(e.center_mhz - 890.0307) < 0.003


def test_sf5_synthetic_roundtrip():
    """Implicit SF5 packets of every CR and assorted lengths, carrier offsets and noise come back byte-exact;
    the (CR, length) search must not settle on a wrong hypothesis."""
    from spectrum_intel.receivers.rc_link import synth
    rng = np.random.default_rng(5)
    fs = 8e6                    # at 2 MSps the synthesised chirp wraps alias in-band (payload peaks -60 %)
    for cr, length, off in ((1, 38, 120e3), (2, 1, -60e3), (3, 17, 0.0), (4, 60, 250e3), (1, 2, -300e3), (2, 101, 33e3)):
        pay = rng.integers(0, 256, length).astype(np.uint8).tobytes()
        iq, _ = synth.lora_packet(5, 500e3, n_pre=8, payload=synth.lora_sf5_implicit_symbols(pay, cr),
                                  cfo_hz=off, fs=fs, rng=rng)
        x = np.r_[np.zeros(12000), iq, np.zeros(12000)].astype(np.complex64)
        x += (0.1 * (rng.standard_normal(len(x)) + 1j * rng.standard_normal(len(x)))).astype(np.complex64)
        got = _detect(x, fs, 0.0)
        assert [e.decoded["payload_hex"] for e in got] == [pay.hex()], (cr, length)
        d = got[0].decoded
        assert d["length"] == length and d["crc_ok"] is True
        assert d["cr"] == f"4/{4 + cr}" or d["cr_ambiguous"]          # a zero-nibble last block fits every CR
        assert abs(got[0].center_mhz * 1e6 - off) < 3e3


def test_lorawan_and_meshtastic_headers():
    up = bytes([0x40, 0x78, 0x56, 0x34, 0x12, 0x80, 0x05, 0x00, 0x0A]) + b"\x00" * 4 + b"\x11\x22\x33\x44"
    h = P.identify(0x34, up)
    assert h["protocol"] == "LoRaWAN" and h["mtype"] == "unconfirmed_data_up"
    assert h["dev_addr"] == "12345678" and h["fcnt"] == 5 and h["adr"] and h["fport"] == 10
    mt = bytes.fromhex("ffffffff" "efbeadde" "01000000") + bytes([0x63, 0x08, 0, 0]) + b"\x00" * 8
    h = P.identify(0x2B, mt)
    assert h["protocol"] == "Meshtastic" and h["to"] == "broadcast" and h["sender"] == "!deadbeef"
    assert h["hop_limit"] == 3 and h["hop_start"] == 3
    assert P.identify(0x12, up) is None                      # private sync word: no protocol claimed


def test_preamble_only_lora_is_not_called_bluetooth():
    from spectrum_intel.core.emission import Emission
    from spectrum_intel.receivers.rc_link.plugin import RcLinkDecoder
    e = Emission(family="lora", label="lora_preamble_sf9_bw125k", start_us=0, end_us=1, f_lo_mhz=433.9,
                 f_hi_mhz=434.1, center_mhz=434.0, snr_db=None, confidence=0.5, features={})
    assert RcLinkDecoder().pretty(e.label, [e]) == "LoRa preamble sf9 bw125k"     # "preamble" contains "ble"


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn(); print("PASS", name)
