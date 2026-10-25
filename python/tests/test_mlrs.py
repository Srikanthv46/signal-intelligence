"""mLRS (research, PLAN R3 / T6; research/drone_transmitters/12_rc_links_more.md section 1): frame CRC, the decodable
SX127x SF6 mode through rc_link detect, mode matching, ExpressLRS kept apart, plugin field, model."""
import numpy as np

from spectrum_intel.receivers.rc_link import detector as D, mlrs as ML, synth as S
from spectrum_intel.receivers.rc_link.plugin import PLUGINS


def test_frame_crc():
    fr = ML.build_frame(0xBEEF, bytes(range(87)))
    d = ML.parse(fr)
    assert d["crc_ok"] and d["frame_sync_word"] == "0xBEEF"
    assert ML.parse(fr[:-1] + bytes([fr[-1] ^ 1])) is None
    assert ML.crc16_mcrf4xx(b"123456789") == 0x6F91                 # CRC-16/MCRF4XX check value


def test_detect_7x_mode():
    rng = np.random.default_rng(3)
    fs, fc = 7.68e6, 915e6
    x = S.noise_capture(1 << 19, rng)
    fr = ML.build_frame(0xBEEF, bytes(rng.integers(0, 256, 87).tolist()))
    iq, _ = S.lora_packet(6, 500e3, n_pre=12, payload=S.lora_implicit_symbols(fr, 6, 1), sync_word=0x12, cfo_hz=3e3,
                          fs=fs, rng=rng)
    y, _ = S.embed(x, fs, fc, [dict(t0_us=2000, f_hz=913.8e6, iq=iq, bw_hz=500e3)], 10)
    em = [e for e in D.detect(y, fs, fc) if e.label.startswith("mlrs")]
    assert len(em) == 1 and em[0].label == "mlrs_lora_sf6_bw500k" and em[0].decoded["frame_sync_word"] == "0xBEEF"
    out = PLUGINS[0].summarise(em, None)
    assert out["fields"]["mlrs"]["frame_sync_words"] == ["0xBEEF"] and out["fields"]["model"]["family"] == "mLRS"
    assert PLUGINS[0].pretty(em[0].label, em) == "mLRS LoRa sf6 bw500k"


def test_mode_match_keeps_elrs_apart():
    exp = D.lora_payload_symbols(6, 91, 5)
    assert [m["mode"] for m in D._mlrs_match(6, 500e3, 12, exp, 915e6)] == ["mlrs_19hz_900", "mlrs_19hz_7x"]
    assert D._mlrs_match(6, 500e3, 8, D.lora_payload_symbols(6, 8, 7), 915e6) == []     # ELRS 200 Hz
