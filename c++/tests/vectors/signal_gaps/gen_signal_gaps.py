#!/usr/bin/env python3
"""Reference for tests/unit/test_signal_gaps.cpp (batch S, research/signal_gaps/PLAN.md), written by the Python
receivers. One section per item; each writes <item>.json (+ binary inputs) next to this file.

  modeac   S11 SSR Mode A/C: Gillham table, decode / summarise on a synthetic 8 MSps envelope (repeats, SPI, garble,
           Mode S-like pulses), and a 10 MSps capture through the adsb plugin (Mode A/C-only signal)
  navaids  S9 / S21 ILS glide slope / localizer and VOR: channel table, detect on two captures (VOR + two-frequency
           localizer with a receiver ppm; two glide slopes), and the captures through the navaids plugin (no ppm)
  gnss     S3 GPS L1 C/A: codes, chi2 thresholds, analyse on synthetic skies (ppm / no ppm, spoofed multipeak sky,
           CW jammer), a sky through the gnss plugin, and four swept jammers (PPD, 2 ms at 61.44 MSps) through swept_ppd
  ble_trackers  S7 BLE trackers: parse_ble_ad on 26 advertising payloads (Apple Find My / iBeacon / AirPods /
           Continuity, Swift Pair, SmartTag, FMDN, Eddystone, Fast Pair, Tile, Chipolo, DULT, negatives), and an AirTag
           advert through rc_link detect + pretty + summarise
  broadcast_fm  S20 FM + RDS: block sync / groups / fields on an encoded bitstream (0A, 2A, 4A, 0B, one damaged
           block), detect on a 1.2 MSps capture (stereo station with RDS, mono station), and the capture through the
           broadcast_fm plugin
  ais      S15 AIS: CRC / NRZI / HDLC frames / fields for types 1, 4, 5, 18, 21, 24A / B (golden gpsd sentence
           included), detect on a 1 MSps capture (both channels, carrier offsets), and the capture through the marine
           plugin
  aprs     S19 APRS / AX.25: AX.25 + APRS parsing of 11 frames (uncompressed, compressed, timestamped, Mic-E, status,
           message, object, ambiguity, malformed, non-APRS, not UI), HDLC deframing of two stuffed frames, the AFSK
           demodulator on a discriminator recording, and a 250 kSps capture through voice_fm
  radalt   S8 radio altimeter: detect on three 10 ms captures at 30.72 MSps (triangle with a turnaround inside the
           capture, 1 kHz sawtooth, a single slow crossing), and the first through the radalt plugin
  ism_devices  S6 315 / 433 MHz: TPMS decode on 6 symbol streams (VDO, Toyota, both polarities, a flipped bit,
           a preamble only), EV1527 / PT2262 / KeeLoq on 5 jittered pulse lists, detect on a 2 MSps scene (two TPMS
           frames, EV1527, PT2262, KeeLoq), and the scene through ism_devices + rc_link (arbitrate)
  barrage  S12 noise jammers: analyse on three 2^18-sample captures at 61.44 MSps (barrage, noise-FM, LTE-like
           OFDM), covers_band, and the barrage capture through the jammer plugins
  deinterleave  S10 pulse deinterleaving: 12 pulse lists (constant, interleaved, PW groups, stagger, dual PRF,
           PRI jitter, profiler, Poisson), and a 10 MSps S-band scene (two radars, one staggered) through the pipeline
           without receivers (pulse features + pulsed-radar evidence)
  bredr    S7 Bluetooth BR/EDR: sync words of 6 LAPs, find on 7 bit lists (0-4 errors, random), try_decode on a
           decimated segment, and a 61.44 MSps capture (GIAC inquiry packet) through rc_link
  am       S13 AM voice: air_channel cases, voice_fm detect on a 2 MSps capture (121.5 MHz AM, 121.9 MHz AM tone,
           an FM carrier), and the capture through the voice_fm plugin (labels, fields, AM audio length)

    nice -n 19 .venv/bin/python -I c++/tests/vectors/signal_gaps/gen_signal_gaps.py [item ...]
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[3] / "python"))


def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items() if not str(k).startswith("_")}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return clean(o.tolist())
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    return o


def dump(name, obj):
    with open(HERE / f"{name}.json", "w") as f:
        json.dump(clean(obj), f, indent=1)
        f.write("\n")


def signals_of(x, fs, fc, receivers):
    """The pipeline's signals of these receivers (label, status, fields) for a capture."""
    from spectrum_intel.analyze import analyze_iq
    from spectrum_intel.decode.registry import plugins
    a = analyze_iq(x, fs, fc, plugins=plugins(enable=receivers))
    return [{"label": s.label, "status": s.decoded.get("status"), "confidence": s.confidence,
             "fields": s.decoded.get("fields"), "f0": s.f0, "f1": s.f1}
            for s in a.signals if s.features.get("decoder") in receivers or s.decoded.get("status") != "not_possible"]


# ------------------------------------------------------------------ S11 Mode A/C
def modeac():
    from spectrum_intel.receivers.adsb import modeac as M
    keys = ["D2", "D4", "A1", "A2", "A4", "B1", "B2", "B4", "C1", "C2", "C4"]
    gill = []
    for v in range(2 ** 11):
        p = dict.fromkeys(M.POSITIONS, 0)
        p.update(zip(keys, [(v >> (10 - i)) & 1 for i in range(11)]))
        gill.append([v, M.gillham_ft(p), M.squawk(p)])
    fs = 8e6
    rng = np.random.default_rng(11)
    n = 60000
    env = np.abs(0.02 * (rng.standard_normal(n) + 1j * rng.standard_normal(n)))
    for t0, code, amp, spi in ((100e-6, "7700", 1.0, False), (400e-6, "7700", 0.9, False), (900e-6, "1234", 0.5, True),
                               (1500e-6, "0000", 0.6, False), (2500e-6, "3017", 0.4, False), (2507e-6, "6543", 0.4, False),
                               (4000e-6, "2000", 0.8, False), (4300e-6, "2000", 0.8, True)):
        r = M.synth_reply(M.bits_for_code(code), fs, amp, spi)
        i = int(t0 * fs)
        env[i:i + len(r)] = np.maximum(env[i:i + len(r)], r)
    for k, b in enumerate(rng.integers(0, 2, 112)):             # a Mode S-like PPM frame at 6 ms
        a = int((6e-3 + 8e-6 + k * 1e-6 + (0 if b else 0.5e-6)) * fs)
        env[a:a + 4] = 1.0
    env = env.astype(np.float32)
    env.tofile(HERE / "modeac_env.f32")
    dec = M.decode(env, fs)
    dec_busy = M.decode(env, fs, busy=[(6e-3, 6e-3 + 120e-6)])
    # pipeline: 10 MSps at 1088 MHz, three 2000 replies + one 7700 pair
    sfs, sfc = 10e6, 1088e6
    sn = 2 ** 17
    rng = np.random.default_rng(5)
    x = (0.01 * (rng.standard_normal(sn) + 1j * rng.standard_normal(sn))).astype(np.complex64)
    t = np.arange(sn) / sfs
    for t0, code in ((1e-3, "2000"), (4e-3, "2000"), (7e-3, "2000"), (9e-3, "7700"), (11e-3, "7700")):
        e = M.synth_reply(M.bits_for_code(code), sfs)
        i = int(t0 * sfs)
        x[i:i + len(e)] += (e * np.exp(2j * np.pi * 2e6 * t[i:i + len(e)])).astype(np.complex64)
    x.tofile(HERE / "modeac_scene.cf32")
    dump("modeac", {"gillham": gill, "env_fs": fs, "decode": dec, "summary": M.summarise(dec),
                    "decode_busy": dec_busy,
                    "scene": {"fs": sfs, "fc": sfc, "signals": signals_of(x, sfs, sfc, ["adsb"])}})


# ------------------------------------------------------------------ S13 AM voice
def am():
    from spectrum_intel.receivers.voice_fm import detector as D
    from spectrum_intel.receivers.voice_fm import synth as S
    chans = [[f, p, D.air_channel(f, p)] for f, p in ((121.5e6, 0.0), (121.503e6, 0.0), (118.2625e6, 0.0),
                                                      (243.0e6, 2.0), (300.71e6, 0.0), (399.99e6, 0.0),
                                                      (446.0e6, 0.0), (136.99e6, 10.0))]
    fs, fc = 2e6, 121.7e6
    x = (S.am(fs, 0.25, -200e3, depth=0.8, snr_db=25, seed=1)
         + S.am(fs, 0.25, 200e3, depth=0.6, audio=800.0, snr_db=28, seed=2)
         + S.nbfm(fs, 0.25, 450e3, audio="speech", ctcss=100.0, snr_db=25, seed=3))
    x = x.astype(np.complex64)
    x.tofile(HERE / "am_scene.cf32")
    em = D.detect(x, fs, fc)
    dump("am", {"air_channel": chans, "fs": fs, "fc": fc,
                "emissions": [{"label": e.label, "features": e.features, "confidence": e.confidence,
                               "f_lo_mhz": e.f_lo_mhz, "f_hi_mhz": e.f_hi_mhz} for e in em],
                "scene": {"fs": fs, "fc": fc, "signals": signals_of(x, fs, fc, ["voice_fm"])}})


# ------------------------------------------------------------------ S9 / S21 navaids
def navaids():
    from spectrum_intel.core.dsp import Channelizer
    from spectrum_intel.receivers.navaids import channels, detect
    from spectrum_intel.receivers.navaids import synth as S
    out = {"channels": [[k, f, p] for k, f, p in channels()]}
    scenes = {
        # VHF: a VOR at 109.00 (radial 211), a two-frequency localizer at 109.90, ppm 5
        "nav_vhf": (4e6, 110.0e6, 5.0, lambda fs, fc: (
            S.vor(fs, 0.0341, 109.0e6 / (1 + 5e-6) - fc, 211.0, snr_db=28, seed=1)
            + S.ils(fs, 0.0341, 109.9e6 / (1 + 5e-6) - fc, 0.24, 0.16, two_freq_hz=9e3, snr_db=30, seed=2))),
        # UHF: glide slopes at 330.95 (below path) and 332.30 (odd SDM -> "-like")
        "nav_uhf": (4e6, 331.5e6, 0.0, lambda fs, fc: (
            S.ils(fs, 0.0341, 330.95e6 - fc, 0.33, 0.47, scale=1.012, snr_db=30, seed=3)
            + S.ils(fs, 0.0341, 332.3e6 - fc, 0.2, 0.2, snr_db=30, seed=4))),
    }
    for name, (fs, fc, ppm, mk) in scenes.items():
        x = mk(fs, fc).astype(np.complex64)
        x.tofile(HERE / f"{name}.cf32")
        em = detect(Channelizer(x, fs, fc), ppm)
        out[name] = {"fs": fs, "fc": fc, "ppm": ppm,
                     "emissions": [{"label": e.label, "features": e.features, "confidence": e.confidence,
                                    "f_lo_mhz": e.f_lo_mhz, "f_hi_mhz": e.f_hi_mhz} for e in em],
                     "signals": signals_of(x, fs, fc, ["navaids"])}
    dump("navaids", out)


# ------------------------------------------------------------------ S3 GNSS
def gnss():
    from spectrum_intel.core.dsp import Channelizer
    from spectrum_intel.receivers.gnss import synth as S
    from spectrum_intel.receivers.gnss.acquire import analyse, threshold
    from spectrum_intel.receivers.gnss.codes import ca_code
    fs, fc = 4e6, 1575.42e6 - 0.3e6
    rng = np.random.default_rng(1)
    sky = [{"prn": p, "cn0": c, "doppler": float(rng.uniform(-4000, 4000)), "code_phase": float(rng.uniform(0, 1023))}
           for p, c in ((3, 48), (7, 44), (12, 40), (19, 36))]
    spoof = [{"prn": p, "cn0": 52.0, "doppler": d, "code_phase": c}
             for p, d, c in ((2, 300, 10), (5, 900, 200), (9, -1200, 400), (14, 2100, 600), (22, -2500, 800))]
    spoof += [dict(x, code_phase=x["code_phase"] - 6.0, cn0=49.0) for x in spoof]
    n = int(0.034 * fs)
    t = np.arange(n) / fs
    noise = S.sky(fs, 0.034, fc, [], seed=5)
    cw = (noise + 3.0 * np.exp(2j * np.pi * (1575.42e6 + 3e3 - fc) * t)).astype(np.complex64)
    cases = {"gnss_sky": (S.sky(fs, 0.034, fc, sky, eps_ppm=9.5, seed=1), 9.5),
             "gnss_sky_noppm": (S.sky(fs, 0.034, fc, sky, eps_ppm=9.5, seed=1), None),
             "gnss_spoof": (S.sky(fs, 0.034, fc, spoof, seed=4), 0.0),
             "gnss_cw": (cw, 0.0)}
    out = {"codes": {p: "".join(map(str, ca_code(p))) for p in range(1, 33)},
           "threshold": [[m, pc, threshold(m, pc)] for m in (1, 4, 8) for pc in (1e-3, 1e-7)], "fs": fs, "fc": fc}
    for name, (x, ppm) in cases.items():
        if name != "gnss_sky_noppm":
            x.tofile(HERE / f"{name}.cf32")
        out[name] = {"ppm": ppm, "result": analyse(Channelizer(x, fs, fc), ppm)}
    sc = S.sky(fs, 0.034, fc, sky, eps_ppm=0.0, seed=6)
    sc.tofile(HERE / "gnss_scene.cf32")
    out["scene"] = {"fs": fs, "fc": fc, "signals": signals_of(sc, fs, fc, ["gnss"])}
    # swept jammer (personal privacy device): 2 ms at 61.44 MSps through swept_ppd, 9 us sawtooth over 20 MHz
    from spectrum_intel.receivers.gnss.acquire import swept_ppd
    pfs, pfc = 61.44e6, 1575.42e6 - 3e6
    cases = [(9e-6, 20e6, 2e6, 3.0), (27e-6, 30e6, -5e6, 1.0), (12e-6, 15e6, 10e6, 1.0), (8e-6, 10e6, 0.0, 1.0)]
    ppd = []
    for k, (per, span, ctr, amp) in enumerate(cases):
        px = S.ppd(pfs, 0.002, 1575.42e6 + ctr - pfc, span, per, amp=amp, seed=7 + k)
        px.tofile(HERE / f"gnss_ppd{k}.cf32")
        ppd.append(swept_ppd(Channelizer(px, pfs, pfc), 1575.42e6))
    out["ppd"] = {"fs": pfs, "fc": pfc, "results": ppd}
    dump("gnss", out)


# ------------------------------------------------------------------ S7 BLE trackers
def ble_trackers():
    import struct
    from spectrum_intel.core.remote_id import parse_ble_ad
    from spectrum_intel.receivers.rc_link import detector as D
    from spectrum_intel.receivers.rc_link import synth as S
    from spectrum_intel.receivers.rc_link.plugin import PLUGINS

    def mfg(cid, payload):
        v = struct.pack("<H", cid) + payload
        return bytes([len(v) + 1, 0xFF]) + v

    def svc(uuid, payload):
        v = struct.pack("<H", uuid) + payload
        return bytes([len(v) + 1, 0x16]) + v

    airtag = mfg(0x004C, bytes([0x12, 0x19, 0x10 | 0x04 | 0xC0]) + bytes(22) + bytes([0x00, 0x00]))
    ads = [airtag,
           mfg(0x004C, bytes([0x12, 0x02, 0x24, 0x00])),                                  # near-owner Find My
           mfg(0x004C, bytes([0x02, 0x15]) + bytes(range(16)) + bytes([0, 1, 0, 2, 0xC5])),  # iBeacon
           mfg(0x004C, bytes([0x07, 0x05, 0x01, 0x0E, 0x20, 0x55, 0x66])),                # AirPods Pro
           mfg(0x004C, bytes([0x07, 0x05, 0x01, 0x99, 0x99, 0x55, 0x66])),                # unknown model
           mfg(0x004C, bytes([0x10, 0x05, 0x01, 0x18, 0x44, 0x00, 0x00, 0x0C, 0x01, 0x00])),  # Continuity
           mfg(0x004C, bytes([0x12])),                                                     # truncated
           mfg(0x0006, bytes([0x03, 0x00, 0x80]) + b"Mouse"), mfg(0x0006, bytes([0x03, 0x01, 0x80])),
           mfg(0x0006, bytes([0x03, 0x05, 0x80])),
           svc(0xFD5A, bytes([0x15 | (2 << 5)]) + bytes(11) + bytes([0xE0]) + bytes(7)),
           svc(0xFD59, bytes([0x15 | (4 << 5)]) + bytes(11) + bytes([0x40]) + bytes(7)),
           svc(0xFD5A, bytes.fromhex("156bfa00c84062b28f00e260c3000000ad018b47")),
           svc(0xFEAA, bytes([0x41]) + bytes(20) + bytes([0x20])), svc(0xFEAA, bytes([0x40]) + bytes(32)),
           svc(0xFEAA, bytes([0x10, 0xEB, 0x03]) + b"example"), svc(0xFEAA, bytes([0x20]) + bytes(13)),
           svc(0xFE2C, bytes([0x12, 0x34, 0x56])), svc(0xFE2C, bytes([0x00, 0x40, 0x11, 0x22])),
           svc(0xFEED, b"\x02\x00\x11\x22"), svc(0xFE33, b"\x00\x01"), svc(0xFD69, b"\x01"),
           svc(0xFCB2, bytes([0x05, 0x01, 0xAA])), svc(0xFCB2, bytes([0x05, 0x00])), svc(0xFCB2, bytes([0x05, 0x03])),
           bytes([2, 0x01, 0x06]) + svc(0xFEED, b"\x02\x00") + bytes([5, 0x09]) + b"TILE"]
    cases = [[a.hex(), parse_ble_ad(a)] for a in ads]
    FC, FS = 2420e6, 61.44e6
    rng = np.random.default_rng(3)
    n = 1 << 18
    x = ((rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)).astype(np.complex64)
    addr = bytes([0x11, 0x22, 0x33, 0x44, 0x55, 0xC6])
    iq, _ = S.ble_adv_burst(addr, airtag, ch=37, cfo_hz=8e3, rng=rng)
    y, _ = S.embed(x, FS, FC, [dict(t0_us=1000, f_hz=2402e6, iq=iq, bw_hz=1.0e6)], 20)
    y = y.astype(np.complex64)
    y.tofile(HERE / "ble_tracker_scene.cf32")
    em = [e for e in D.detect(y, FS, FC) if (e.decoded or {}).get("crc_ok")]
    rc = next(p for p in PLUGINS if p.name == "rc_link")
    out = rc.summarise(em, None)
    dump("ble_trackers", {"ad": cases, "fs": FS, "fc": FC, "label": rc.pretty("ble", em),
                          "ble_devices": out["fields"].get("ble_devices"), "evidence": out["evidence"],
                          "decoded": [e.decoded for e in em]})


# ------------------------------------------------------------------ S20 FM broadcast + RDS
def broadcast_fm():
    from spectrum_intel.core.dsp import Channelizer
    from spectrum_intel.receivers.broadcast_fm import detect
    from spectrum_intel.receivers.broadcast_fm import rds as R
    from spectrum_intel.receivers.broadcast_fm import synth as S
    rng = np.random.default_rng(1)
    bits = [int(b) for b in rng.integers(0, 2, 17)]
    for g in S.ps_groups(0x5201, 10, "AIR FM"):
        bits += R.encode_group(g)
    bits += R.encode_group((0x5201, (2 << 12) | (10 << 5) | 0, 0x4845, 0x4C4C))
    bits += R.encode_group((0x5201, (2 << 12) | (10 << 5) | 1, 0x4F0D, 0x2020))
    mjd = 61338
    bits += R.encode_group((0x5201, (4 << 12) | (10 << 5) | (mjd >> 15), ((mjd & 0x7FFF) << 1) | 0,
                            ((10 & 15) << 12) | (30 << 6) | 11))
    bits += R.encode_group((0x5202, 0x0545, 0x1234, 0x4142), version_b=True)
    bits[300] ^= 1
    blks = R.blocks(bits)
    fs, fc = 1.2e6, 98.0e6
    x = (S.fm_station(fs, 0.2, 300e3, S.ps_groups(0x5201, 10, "AIR FM"), seed=1)
         + S.fm_station(fs, 0.2, -400e3, None, stereo=False, seed=2)).astype(np.complex64)
    x.tofile(HERE / "fm_scene.cf32")
    em = detect(Channelizer(x, fs, fc))
    dump("broadcast_fm", {"bits": "".join(map(str, bits)),
                          "blocks": [[p, o, d] for p, o, d in blks], "decode": R.decode(R.groups(blks)),
                          "fs": fs, "fc": fc,
                          "emissions": [{"label": e.label, "features": e.features, "confidence": e.confidence,
                                         "f_lo_mhz": e.f_lo_mhz, "f_hi_mhz": e.f_hi_mhz} for e in em],
                          "scene": {"fs": fs, "fc": fc, "signals": signals_of(x, fs, fc, ["broadcast_fm"])}})


# ------------------------------------------------------------------ S15 AIS
def ais():
    from spectrum_intel.core.dsp import Channelizer
    from spectrum_intel.receivers.marine import ais as A
    from spectrum_intel.receivers.marine import detect
    from spectrum_intel.receivers.marine import synth as S
    P = "177KQJ5000G?tO`K>RA1wUbN0TKH"

    def u(v, n):
        v &= (1 << n) - 1
        return [(v >> (n - 1 - i)) & 1 for i in range(n)]

    def six(t, n):
        v = [ord(c) - 64 if 64 <= ord(c) < 96 else ord(c) for c in t.ljust(n, "@")[:n]]
        return [(x >> (5 - i)) & 1 for x in v for i in range(6)]

    msgs = {
        "t1": A.armor_to_bits(P),
        "t5": (u(5, 6) + u(0, 2) + u(419012345, 30) + u(0, 2) + u(9876543, 30) + six("VTAB1", 7)
               + six("SAGAR KANYA", 20) + u(70, 8) + u(80, 9) + u(20, 9) + u(6, 6) + u(6, 6) + u(1, 4) + u(10, 4)
               + u(25, 5) + u(12, 5) + u(30, 6) + u(65, 8) + six("MUMBAI", 20) + [0, 0, 0, 0]),
        "t18": (u(18, 6) + u(0, 2) + u(419000001, 30) + u(0, 8) + u(52, 10) + [1] + u(int(72.85 * 600000), 28)
                + u(int(18.95 * 600000), 27) + u(900, 12) + u(90, 9) + u(0, 35)),
        "t4": (u(4, 6) + u(0, 2) + u(4190001, 30) + u(2026, 14) + u(10, 4) + u(25, 5) + u(9, 5) + u(15, 6)
               + u(30, 6) + [1] + u(int(72.8 * 600000), 28) + u(int(18.9 * 600000), 27) + u(1, 4) + u(0, 34)),
        "t21": (u(21, 6) + u(0, 2) + u(994190001, 30) + u(14, 5) + six("BUOY NORTH", 20) + [0]
                + u(int(72.9 * 600000), 28) + u(int(19.0 * 600000), 27) + u(0, 53)),
        "t24a": u(24, 6) + u(0, 2) + u(419000002, 30) + u(0, 2) + six("FISHER ONE", 20),
        "t24b": u(24, 6) + u(0, 2) + u(419000002, 30) + u(1, 2) + u(30, 8) + u(0, 42) + six("VTFS1", 7) + u(0, 36),
    }
    codec = {}
    for k, b in msgs.items():
        fb = A.frame_bits(b)
        lv = A.nrzi(fb)
        fr = A.frames([0, 1, 0, 1, 0] + A.unnrzi([1] + lv))
        codec[k] = {"levels": "".join(map(str, [1] + lv)), "frames": [[p, o.hex()] for p, o in fr],
                    "fields": [A.fields(o) for _, o in fr], "crc": A.crc_x25(A.bits_to_octets(b))}
    fs, fc = 1e6, 162.0e6
    x = S.capture(fs, 0.1, fc, [(0.005, 161.975e6 + 1500, P), (0.04, 162.025e6 - 900, P)], snr_db=18, seed=1)
    x.tofile(HERE / "ais_scene.cf32")
    em = detect(Channelizer(x, fs, fc))
    dump("ais", {"codec": codec, "fs": fs, "fc": fc,
                 "emissions": [{"decoded": e.decoded, "features": e.features, "start_us": e.start_us} for e in em],
                 "scene": {"fs": fs, "fc": fc, "signals": signals_of(x, fs, fc, ["marine"])}})


# ------------------------------------------------------------------ S19 APRS / AX.25
def aprs():
    from spectrum_intel.core import hdlc as H
    from spectrum_intel.receivers.voice_fm import afsk
    from spectrum_intel.receivers.voice_fm import synth as S

    def ax(dest, src, path, info):
        def addr(call, last, digi=False):
            c, _, ssid = call.partition("-")
            b = bytes((ord(ch) << 1) for ch in c.ljust(6)[:6])
            return b + bytes([0x60 | (int(ssid or 0) << 1) | (1 if last else 0) | (0x80 if digi else 0)])
        out = addr(dest, False) + addr(src, not path)
        for k, q in enumerate(path):
            out += addr(q.rstrip("*"), k == len(path) - 1, q.endswith("*"))
        return out + bytes([0x03, 0xF0]) + info.encode("latin-1")

    frames = [ax("APRS", "VU2ABC-9", [], "!1259.00N/07735.00E>Test"),
              ax("APDW16", "VU3XYZ-5", ["WIDE1-1*", "WIDE2-1"], ">status text"),
              ax("APRS", "VU2DEF", ["WIDE2-2"], "=/5L!!<*e7>7P[Comment"),
              ax("APRS", "VU2DEF-1", [], "@092345z4903.50N/07201.75W>088/036Test"),
              ax("T2SP0W", "VU2MIC-7", ["WIDE1-1"], "`(_fn\"Oj/]=Hi"),
              ax("APRS", "VU2MSG", [], ":VU2XYZ-7 :hello{01"),
              ax("APRS", "VU2OBJ", [], ";OBJECT   *092345z1259.00N/07735.00E>obj"),
              ax("APRS", "VU2AMB", [], "!12  .  N/077  .  E>ambiguous"),
              ax("APRS", "VU2BAD", [], "!1259.0xN/07735.00E>bad"),
              ax("BEACON", "VU2BCN", [], "plain text \xe9"),
              bytes([0x82, 0xA0, 0xA4, 0xA6, 0x40, 0x40, 0xE0, 0x01]) + bytes(20)]
    cases = []
    for o in frames:
        f = H.ax25(o)
        cases.append([o.hex(), f, H.aprs(f) if f else None])
    stream = H.unnrzi([1] + H.nrzi(H.frame_bits(frames[0], 3) + H.frame_bits(frames[1], 2)))
    d = []
    for seed, o in ((1, frames[2]), (2, frames[4])):
        x = S.afsk_packet(48e3, o, 0.0, snr_db=18, seed=seed)
        d.append((np.angle(x[1:] * np.conj(x[:-1])) * np.float32(48e3 / (2 * np.pi))).astype(np.float32))
    dd = np.concatenate(d)
    dd.tofile(HERE / "aprs_disc.f32")
    fs, fc = 250e3, 144.3e6
    sc = (S.afsk_packet(fs, frames[0], 90e3, snr_db=20, seed=3)).astype(np.complex64)
    sc.tofile(HERE / "aprs_scene.cf32")
    dump("aprs", {"frames": cases, "stream": "".join(map(str, stream)),
                  "deframe": [[p, o.hex()] for p, o in H.deframe(stream)],
                  "afsk": afsk.decode(dd, 48e3), "fs": fs, "fc": fc,
                  "scene": {"fs": fs, "fc": fc, "signals": signals_of(sc, fs, fc, ["voice_fm"])}})


# ------------------------------------------------------------------ S8 radio altimeter
def radalt():
    from spectrum_intel.receivers.radalt import detect
    from spectrum_intel.receivers.radalt import synth as S
    fs = 30.72e6
    cases = []
    for i, (fc, bw, rep, shape, amp) in enumerate(((4250e6, 133e6, 300, "triangle", 0.4),
                                                   (4300e6, 176.8e6, 1000, "sawtooth", 0.6),
                                                   (4300e6, 104e6, 50, "triangle", 1.0))):
        x = S.fmcw(fs, 0.010, fc, bw=bw, f_rep=rep, shape=shape, amp=amp, seed=2)
        x.tofile(HERE / f"radalt_{i}.cf32")
        em = detect(x, fs, fc)
        cases.append({"fs": fs, "fc": fc, "emissions": [
            {"label": e.label, "confidence": e.confidence, "features": e.features, "start_us": e.start_us,
             "end_us": e.end_us, "f_lo_mhz": e.f_lo_mhz, "f_hi_mhz": e.f_hi_mhz} for e in em]})
    x = np.fromfile(HERE / "radalt_0.cf32", np.complex64)
    dump("radalt", {"cases": cases,
                    "scene": {"fs": fs, "fc": 4250e6, "signals": signals_of(x, fs, 4250e6, ["radalt"])}})


# ------------------------------------------------------------------ S6 315 / 433 MHz devices
def ism_devices():
    from spectrum_intel.core.dsp import Channelizer
    from spectrum_intel.receivers.ism_devices import detect
    from spectrum_intel.receivers.ism_devices import ook as O
    from spectrum_intel.receivers.ism_devices import synth as S
    from spectrum_intel.receivers.ism_devices import tpms as T
    vdo = S.vdo_symbols(0x21, 0xDEADBEEF, 2, 3, 170, 75, 0x20)
    bad = list(vdo)
    bad[100] ^= 1
    bad[101] ^= 1
    streams = [vdo, [1 - v for v in vdo], S.toyota_symbols(0x12345678, 0x80, 150, 65),
               [1 - v for v in S.toyota_symbols(0xA1B2C3D4, 0x00, 120, 30)], bad, [0, 1] * 120]
    sym = [{"s": "".join(map(str, st)), "decode": T.decode(st)} for st in streams]

    def plist(train, jitter=0.0, seed=0):
        rng = np.random.default_rng(seed)
        out, t = [], 0.0
        for m, sp in train:
            m2, s2 = m + jitter * rng.standard_normal(), sp + jitter * rng.standard_normal()
            out.append([t, m2, s2])
            t += m2 + s2
        out[-1][2] = float("inf")
        return out
    trains = [plist(S.ev1527(0xA5C3E2, 350, 3), 15, 1), plist(S.pt2262("0F1F0011FF10", 100, 4), 8, 2),
              plist(S.ev1527(0x0F0F0F, 280, 1), 10, 3), plist(S.keeloq(0x1234567, 0b0010, 0xCAFEBABE, 400, 2), 20, 4),
              plist(S.keeloq(0x0ABCDEF, 0b1111, 0x13572468, 300, 1), 10, 5)]
    pulses = [{"p": [[a, b, c if c != float("inf") else -1.0] for a, b, c in tr],
               "fixed_code": O.fixed_code([tuple(q) for q in tr]), "keeloq": O.keeloq([tuple(q) for q in tr])}
              for tr in trains]
    fs, fc = 2e6, 433.92e6
    x = S.fsk(fs, 0.2, -0.3e6, [(3000, S.vdo_symbols(0x21, 0xDEADBEEF, 2, 3, 170, 75, 0x20)),
                                (40000, S.toyota_symbols(0x12345678, 0x80, 150, 65))], snr_db=21, seed=1)
    x += S.ook(fs, 0.2, 0.1e6, S.ev1527(0xA5C3E2, 300, 3), snr_db=21, seed=2)
    x += S.ook(fs, 0.2, 0.5e6, S.pt2262("0F1F0011FF10", 100, 4), snr_db=21, seed=3)
    x += S.ook(fs, 0.2, -0.65e6, S.keeloq(0x0ABCDEF, 0b1000, 0x13572468, 400, 1), snr_db=21, seed=4)
    x = x.astype(np.complex64)
    x.tofile(HERE / "ism_scene.cf32")
    em = detect(Channelizer(x, fs, fc))
    dump("ism_devices", {"symbols": sym, "pulses": pulses, "fs": fs, "fc": fc,
                         "emissions": [{"label": e.label, "confidence": e.confidence, "start_us": e.start_us,
                                        "end_us": e.end_us, "center_mhz": e.center_mhz, "features": e.features,
                                        "decoded": e.decoded} for e in em],
                         "scene": {"fs": fs, "fc": fc, "signals": signals_of(x, fs, fc, ["ism_devices", "rc_link"])}})


# ------------------------------------------------------------------ S12 barrage / noise-FM jammers
def barrage():
    from spectrum_intel.core.dsp import Channelizer
    from spectrum_intel.receivers.jammer import barrage as B
    from spectrum_intel.receivers.jammer import synth as S
    fs, fc, n = S.FS, 1842.5e6, 2 ** 18
    cases = []
    for i, (x, f0, f1) in enumerate(((S.barrage(-10e6, 10e6, 12, seed=3), -10e6, 10e6),
                                     (S.noise_fm(5e6), -12e6, 12e6),
                                     (S.ofdm(18e6, snr_db=3), -9e6, 9e6))):
        x = x[:n].astype(np.complex64)
        x.tofile(HERE / f"barrage_{i}.cf32")
        cases.append({"f0": fc + f0, "f1": fc + f1, "analyse": B.analyse(Channelizer(x, fs, fc), fc + f0, fc + f1)})
    x = np.fromfile(HERE / "barrage_0.cf32", np.complex64)
    dump("barrage", {"fs": fs, "fc": fc, "cases": cases, "covers": [[1805e6, 1880e6, B.covers_band(1805e6, 1880e6)],
                                                                    [2400e6, 2483.5e6, B.covers_band(2400e6, 2483.5e6)]],
                     "scene": {"fs": fs, "fc": fc, "signals": signals_of(x, fs, fc, ["jammer"])}})


# ------------------------------------------------------------------ S10 pulse deinterleaving
def deinterleave():
    from spectrum_intel.analyze import analyze_iq
    from spectrum_intel.core.deinterleave import deinterleave as D
    rng = np.random.default_rng(1)

    def train(pri, n, t0=0.0, pw=1e-6, miss=0.1, stagger=None, prij=0.0):
        t, x, k = [], t0, 0
        while len(t) < n:
            step = stagger[k % len(stagger)] if stagger else pri
            k += 1
            if rng.random() > miss:
                t.append(round(x * 1e6) / 1e6)
            x += step * (1 + prij * (2 * rng.random() - 1))
        return t, [pw] * len(t)
    sets = [[train(1000e-6, 30)], [train(1000e-6, 30), train(1370e-6, 22, t0=123e-6)],
            [train(833e-6, 35), train(500e-6, 60, t0=77e-6, pw=4.7e-6), train(2700e-6, 12, t0=300e-6, pw=20e-6)],
            [train(None, 40, stagger=[900e-6, 1000e-6, 1100e-6])], [train(1 / 900, 30), train(1 / 1200, 40, t0=211e-6)],
            [train(1000e-6, 30, prij=0.05)], [train(100e-6, 300)]]
    cases = []
    for st in sets:
        t = sum((a for a, _ in st), [])
        w = sum((b for _, b in st), [])
        cases.append({"toa": t, "pw": w, "out": D(t, w)})
    for sd in range(5):
        r = np.random.default_rng(100 + sd)
        n = int(r.integers(20, 200))
        t = [float(v) for v in sorted(np.round(r.uniform(0, 0.034, n) * 1e6) / 1e6)]
        cases.append({"toa": t, "pw": [1e-6] * n, "out": D(t, [1e-6] * n)})
    fs, fc, dur = 10e6, 2800e6, 0.034
    n = int(fs * dur)
    x = (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)

    def radar(pri, pw, t0, amp, f_off, miss=0.05, stagger=None):
        t, k = t0, 0
        while t < dur - pw:
            if rng.random() > miss:
                a, b = int(t * fs), int((t + pw) * fs)
                x[a:b] += amp * np.exp(2j * np.pi * f_off * np.arange(a, b) / fs)
            t += stagger[k % len(stagger)] if stagger else pri
            k += 1
    radar(1 / 1000, 1e-6, 100e-6, 30, 0.5e6)
    radar(None, 20e-6, 333e-6, 20, 0.5e6, stagger=[1 / 1100, 1 / 1200, 1 / 1000])
    x = x.astype(np.complex64)
    x.tofile(HERE / "radar_scene.cf32")
    a = analyze_iq(x, fs, fc, plugins=[])
    sig = [{"label": s.label, "confidence": s.confidence, "evidence": s.evidence,
            "features": {k: s.features.get(k) for k in ("pulse_n", "pri_s", "pulse_emitters", "pulse_unassigned")}}
           for s in a.signals if s.features.get("pulse_n")]
    dump("deinterleave", {"cases": cases, "scene": {"fs": fs, "fc": fc, "signals": sig}})


# ------------------------------------------------------------------ S7 Bluetooth BR/EDR access codes
def bredr():
    from spectrum_intel.receivers.rc_link import bredr as B
    from spectrum_intel.receivers.rc_link import synth as S
    laps = [0, B.GIAC, 0x123456, 0x800000, 0xFFFFFF, 0x2C1F07]
    sw = [[lap, f"{B.syncword(lap):016x}"] for lap in laps]
    rng = np.random.default_rng(5)
    finds = []
    for lap, e in ((0x5A5A5A, 0), (0x5A5A5A, 3), (0x5A5A5A, 4), (B.GIAC, 2)):
        b = S.bredr_bits(lap, rng=rng).tolist()
        for q in range(e):
            b[4 + 7 * q] ^= 1
        finds.append({"bits": "".join(map(str, b)), "find": B.find(b)})
    for _ in range(3):
        b = rng.integers(0, 2, 300).tolist()
        finds.append({"bits": "".join(map(str, b)), "find": B.find(b)})
    fs, fc = S.FS, 2440e6
    iq, _ = S.bredr_packet(B.GIAC, rng=np.random.default_rng(3))
    y, _ = S.embed(S.noise_capture(200_000, np.random.default_rng(4)), fs, fc,
                   [dict(iq=iq, f_hz=2431e6, t0_us=1000.0, bw_hz=1e6)], 16)
    y = y.astype(np.complex64)
    y.tofile(HERE / "bredr_scene.cf32")
    from spectrum_intel.receivers.rc_link.detector import _slice
    yb, ob, _ = _slice(np.fft.fft(y), fs, 2431e6 - fc, 0.75e6, 8e6)       # as rc_link re-slices a 1 Mb/s hop
    seg = yb[int(950e-6 * ob):int(1500e-6 * ob)].astype(np.complex64)
    seg.tofile(HERE / "bredr_seg.cf32")
    dump("bredr", {"syncwords": sw, "finds": finds, "seg_ofs": ob, "try_decode": B.try_decode(seg, ob),
                   "scene": {"fs": fs, "fc": fc, "signals": signals_of(y, fs, fc, ["rc_link"])}})


ITEMS = {"modeac": modeac, "am": am, "navaids": navaids, "gnss": gnss, "ble_trackers": ble_trackers,
         "broadcast_fm": broadcast_fm, "ais": ais, "aprs": aprs, "radalt": radalt, "ism_devices": ism_devices, "barrage": barrage, "deinterleave": deinterleave, "bredr": bredr}

if __name__ == "__main__":
    for name in sys.argv[1:] or list(ITEMS):
        ITEMS[name]()
        print("wrote", name)
