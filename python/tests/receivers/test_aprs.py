"""APRS / AX.25 packet radio (research, signal gaps S19): HDLC + AX.25 + APRS parsing, Bell 202 AFSK on FM."""
import numpy as np

from spectrum_intel.core import hdlc as H
from spectrum_intel.receivers.voice_fm import afsk
from spectrum_intel.receivers.voice_fm import synth as S

FRAME = bytes.fromhex("82a0a4a64040e0" "acaa64828486" "73" "03f0") + b"!1259.00N/07735.00E>Test"   # research file 06


def ax25(dest, src, path, info):
    def addr(call, last, digi=False):
        c, _, ssid = call.partition("-")
        b = bytes((ord(ch) << 1) for ch in c.ljust(6)[:6])
        return b + bytes([0x60 | (int(ssid or 0) << 1) | (1 if last else 0) | (0x80 if digi else 0)])
    out = addr(dest, False) + addr(src, not path)
    for k, p in enumerate(path):
        out += addr(p.rstrip("*"), k == len(path) - 1, p.endswith("*"))
    return out + bytes([0x03, 0xF0]) + info.encode("latin-1")


def test_worked_example_fcs_and_parse():
    assert H.crc_x25(FRAME) == 0xB7A6
    fr = H.deframe(H.unnrzi([1] + H.nrzi(H.frame_bits(FRAME, 3))))
    assert len(fr) == 1 and fr[0][1] == FRAME
    f = H.ax25(FRAME)
    assert f["src"] == "VU2ABC-9" and f["dest"] == "APRS" and f["path"] == []
    a = H.aprs(f)
    assert a["lat"] == 12.98333 and a["lon"] == 77.58333 and a["symbol"] == "/>" and a["comment"] == "Test"


def test_aprs_formats():
    c = H.aprs({"dest": "APRS", "info": "=/5L!!<*e7>7P[Comment"})                 # APRS 1.0.1 compressed example
    assert c["lat"] == 49.5 and c["lon"] == -72.75 and c["course_deg"] == 88 and c["speed_kn"] == 36.2
    t = H.aprs({"dest": "APRS", "info": "@092345z4903.50N/07201.75W>088/036Test"})
    assert t["timestamp"] == "092345z" and t["course_deg"] == 88 and t["speed_kn"] == 36
    m = H.aprs({"dest": "APRS", "info": ":VU2XYZ-7 :hello{01"})
    assert m["addressee"] == "VU2XYZ-7" and m["message"] == "hello{01"
    assert H.aprs({"dest": "APRS", "info": ">On the air"})["status"] == "On the air"
    e = H.aprs({"dest": "T2SP0W", "info": "`(_fn\"Oj/]=Hi"})
    assert e["format"] == "Mic-E" and e["lat"] == 42.50117 and e["lon"] == -12.129


def test_ax25_path_and_rejects():
    o = ax25("APDW16", "VU3XYZ-5", ["WIDE1-1*", "WIDE2-1"], ">status")
    f = H.ax25(o)
    assert f["path"] == ["WIDE1-1*", "WIDE2-1"] and f["src"] == "VU3XYZ-5"
    assert H.ax25(o[:14] + bytes([0x13, 0xF0]) + o[16:]) is None                 # not a UI frame
    assert H.ax25(bytes([0x01] * 20)) is None


def test_afsk_on_fm_with_noise():
    o = ax25("APRS", "VU2ABC-9", ["WIDE2-2"], "!1259.00N/07735.00E>Bengaluru")
    ok = 0
    for seed in range(4):
        x = S.afsk_packet(48e3, o, 0.0, snr_db=15, seed=seed)
        d = np.angle(x[1:] * np.conj(x[:-1])) * (48e3 / (2 * np.pi))
        r = afsk.decode(d, 48e3)
        ok += bool(r) and r[0]["aprs"]["comment"] == "Bengaluru"
    assert ok >= 3


def test_pipeline_aprs_signal_not_voice():
    from spectrum_intel.analyze import analyze_iq
    from spectrum_intel.decode.registry import plugins
    x = S.afsk_packet(1e6, FRAME, 390e3, snr_db=20)
    a = analyze_iq(x, 1e6, 144.0e6, plugins=plugins(enable=["voice_fm"]))
    s = [s for s in a.signals if s.features.get("decoder") == "voice_fm"]
    assert [x.label for x in s] == ["APRS (VU2ABC-9)"] and s[0].decoded["status"] == "success"
    pk = s[0].decoded["fields"]["packets"][0]
    assert pk["aprs"]["lat"] == 12.98333 and "audio" not in s[0].decoded["fields"]
