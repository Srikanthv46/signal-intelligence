"""315 / 433 MHz remotes and TPMS (research, signal gaps S6) on synthetic devices."""
import numpy as np
import pytest

from spectrum_intel.core.dsp import Channelizer
from spectrum_intel.receivers.ism_devices import detect
from spectrum_intel.receivers.ism_devices import ook as O
from spectrum_intel.receivers.ism_devices import synth as S
from spectrum_intel.receivers.ism_devices import tpms as T

FS, FC = 2e6, 433.92e6


def test_crc8_and_frames():
    assert T.crc8(b"123456789", 0x07, 0x00) == 0xF4                    # CRC-8/SMBUS check value
    s = S.vdo_symbols(0x21, 0xDEADBEEF, 2, 3, 170, 75, 0x20)
    r = T.decode(s)
    assert r["id"] == "deadbeef" and r["pressure_kpa"] == 233.75 and r["temperature_c"] == 25 and r["repeat"] == 3
    assert T.decode([1 - v for v in s])["id"] == "deadbeef"           # either FSK polarity
    r = T.decode(S.toyota_symbols(0x12345678, 0x80, 150, 65))
    assert r["protocol"].startswith("Toyota") and r["pressure_psi"] == 30.5 and r["temperature_c"] == 25
    bad = S.vdo_symbols(0x21, 0xDEADBEEF, 2, 3, 170, 75, 0x20)
    bad[100] ^= 1
    bad[101] ^= 1                                                     # one flipped data bit: CRC fails
    assert T.decode(bad) is None


def test_ook_pulse_decoders():
    p, t = [], 0.0
    for m, s in S.ev1527(0xA5C3E2, 350, 3):
        p.append((t, m, s))
        t += m + s
    r = O.fixed_code(p)
    assert r["code"] == "A5C3E2" and r["id"] == "A5C3E" and r["data"] == 2 and r["repeats"] == 3
    assert "tristate" not in r
    p, t = [], 0.0
    for m, s in S.keeloq(0x1234567, 0b0010, 0xCAFEBABE, 400, 2):
        p.append((t, m, s))
        t += m + s
    r = O.keeloq(p)
    assert r["serial"] == "1234567" and r["matching_packets"] == 2 and not r["battery_low"]
    assert "CAFEBABE" not in str(r).upper()                          # the hopping code is never reported


@pytest.mark.parametrize("snr", [6, 15])
def test_scene(snr):
    x = S.fsk(FS, 0.25, -0.3e6, [(3000, S.vdo_symbols(0x21, 0xDEADBEEF, 2, 3, 170, 75, 0x20)),
                                 (60000, S.toyota_symbols(0x12345678, 0x80, 150, 65))], snr_db=snr + 6, seed=1)
    x += S.ook(FS, 0.25, 0.1e6, S.ev1527(0xA5C3E2, 350, 4), snr_db=snr + 6, seed=2)
    x += S.ook(FS, 0.25, 0.5e6, S.pt2262("0F1F0011FF10", 100, 5), snr_db=snr + 6, seed=3)
    x += S.ook(FS, 0.25, -0.6e6, S.keeloq(0x0ABCDEF, 0b1000, 0x13572468, 400, 2), snr_db=snr + 6, seed=4)
    em = detect(Channelizer(x, FS, FC))          # four noise floors: +6 dB in the synth keeps snr per device
    got = sorted((e.label, (e.decoded or {}).get("id") or (e.decoded or {}).get("serial")) for e in em)
    assert got == [("keeloq", "0ABCDEF"), ("ook_fixed_code", "1D0F5"), ("ook_fixed_code", "A5C3E"),
                   ("tpms", "12345678"), ("tpms", "deadbeef")]
    pt = [e for e in em if (e.decoded or {}).get("id") == "1D0F5"][0]
    assert pt.decoded["tristate"] == "0F1F0011FF10" and pt.confidence == 0.75


def test_full_rate_block():
    fs = 61.44e6
    x = S.fsk(fs, 2 ** 21 / fs, 433.92e6 - 435e6, [(5000, S.vdo_symbols(0x21, 0x0BADF00D, 0, 1, 160, 70, 16))],
              snr_db=25, seed=4)
    em = detect(Channelizer(x, fs, 435e6))
    assert [(e.label, e.decoded["id"]) for e in em] == [("tpms", "0badf00d")]


def test_negatives():
    from spectrum_intel.receivers.rc_link import synth as RS
    from spectrum_intel.receivers.voice_fm import synth as V
    for sd in range(6):
        r = np.random.default_rng(sd)
        x = 0.7 * (r.standard_normal(400000) + 1j * r.standard_normal(400000))
        assert detect(Channelizer(x.astype(np.complex64), FS, FC)) == []
    iq, _ = RS.lora_packet(9, 125e3, fs=FS, rng=np.random.default_rng(1))
    x = np.zeros(400000, complex)
    x[1000:1000 + len(iq)] = 10 * iq[:399000]
    x = x * np.exp(2j * np.pi * 0.2e6 * np.arange(len(x)) / FS)
    x += np.sqrt(0.5) * (np.random.default_rng(2).standard_normal(len(x)) + 1j * np.random.default_rng(3).standard_normal(len(x)))
    assert detect(Channelizer(x.astype(np.complex64), FS, FC)) == []              # a LoRa chirp is not FSK
    assert detect(Channelizer(V.nbfm(FS, 0.2, 0.3e6, audio="speech", snr_db=25, seed=1), FS, FC)) == []
    x = S.ook(FS, 0.2, 0.1e6, S.ev1527(0xA5C3E2, 350, 3), snr_db=15)
    assert detect(Channelizer(x, FS, 400e6)) == []                                 # out of the bands


def test_pipeline_labels():
    from spectrum_intel.analyze import analyze_iq
    from spectrum_intel.api import contract as K
    from spectrum_intel.decode.registry import plugins
    x = S.fsk(FS, 0.1, -0.3e6, [(3000, S.vdo_symbols(0x21, 0xDEADBEEF, 2, 3, 170, 75, 0x20))], snr_db=20, seed=1)
    x += S.ook(FS, 0.1, 0.3e6, S.ev1527(0xA5C3E2, 300, 3), snr_db=20, seed=2)
    a = analyze_iq(x, FS, FC, plugins=plugins(enable=["ism_devices"]))
    sig = sorted((s.label, s.decoded["status"]) for s in a.signals if s.features.get("decoder") == "ism_devices")
    assert sig == [("OOK fixed-code remote (EV1527 / PT2262-like, ID A5C3E)", "partial"),
                   ("TPMS (Continental / VDO, ID deadbeef)", "success")]
    assert K.technology(sig[1][0], "ism_devices")[0] == "tpms"
    assert K.technology(sig[0][0], "ism_devices")[0] == "ook_remote"
