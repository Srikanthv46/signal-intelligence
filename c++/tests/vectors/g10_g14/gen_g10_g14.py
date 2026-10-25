#!/usr/bin/env python3
"""Synthetic references for PLAN G10 / G12 / G14 (research): a 2.4 GHz capture (7.68 MSps, 2^18 samples, fc 2430 MHz)
with HoTT, Hitec and S-FHSS packets and an 868 MHz capture (fc 868.3 MHz) with ADS-L and OGN packets, both at 15 dB:
rc_link detect's emissions of those labels (label, decoded), plugin fields and model blocks; FANET payloads through the
LoRa identify(); a UAT message at 4 samples / bit through uat.demod and dump978 sample messages through uat.parse;
Reed-Solomon encode / decode cases; the China GB 46750 DJI Mini 5 Pro frame through parse_odid.

    nice -n 19 .venv/bin/python -I c++/tests/vectors/g10_g14/gen_g10_g14.py
"""
import json
import random
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parents[3] / "python"))
from spectrum_intel.core import remote_id as R                                                  # noqa: E402
from spectrum_intel.core.rs import RS                                                           # noqa: E402
from spectrum_intel.receivers.adsb import uat as U                                              # noqa: E402
from spectrum_intel.receivers.rc_link import adsl as A, detector as D, hitec as HI, hott as H   # noqa: E402
from spectrum_intel.receivers.rc_link import lora_phy as P, ogn as O, synth as S                 # noqa: E402
from spectrum_intel.receivers.rc_link.plugin import PLUGINS                                     # noqa: E402

FS = 7.68e6
KEYS = ("cc2500_rc", "adsl", "ogn", "model")


def run(name, fc, bursts, labels, seed):
    rng = np.random.default_rng(seed)
    x = S.noise_capture(1 << 18, rng)
    y, _ = S.embed(x, FS, fc, bursts, 15)
    y = y.astype(np.complex64)
    y.tofile(OUT / f"{name}.cf32")
    em = [e for e in D.detect(y, FS, fc) if e.label in labels]
    fields = {}
    for lab in labels:
        g = [e for e in em if e.label == lab]
        if g:
            f = PLUGINS[0].summarise(g, None)["fields"]
            fields[lab] = {k: f.get(k) for k in KEYS if f.get(k) is not None}
    print(name, [e.label for e in em])
    return dict(fs=FS, fc=fc, labels=[e.label for e in em], decoded=[e.decoded for e in em], fields=fields)


rng = np.random.default_rng(91)
hp = H.tx_packet(5, 33, [3000, 2800, 3200] + [3000] * 9, bytes.fromhex("ea4d000150"), bytes.fromhex("6eee0100b1"))
b24 = [dict(t0_us=500, f_hz=H.F0 + 25 * H.STEP, iq=S.hott_burst(hp, cfo_hz=4e3, fs=FS, rng=rng)[0], bw_hz=0.5e6)]
for k in range(2):
    iq, _ = S.hitec_burst(HI.hitec_frame(bytes.fromhex("1ca12b"), [10822, 7047, 14597] + [10822] * 6, k), 0x2BA1,
                          cfo_hz=3e3, fs=FS, rng=rng)
    b24.append(dict(t0_us=6000 + 2736 * k, f_hz=HI.F0 + 80 * HI.STEP, iq=iq, bw_hz=0.5e6))
sp = HI.sfhss_packet(0x5A3C, 17, [1520, 1020, 2020, 1600], 9, model=bytes.fromhex("9a06"))
for k in range(2):
    iq, _ = S.sfhss_burst(sp, 0xFFFF, cfo_hz=2e3, fs=FS, rng=rng)
    b24.append(dict(t0_us=14000 + 1625 * k, f_hz=HI.SF_F0 + 17 * HI.SF_STEP, iq=iq, bw_hz=0.3e6))
ref = {"c24": run("c24", 2430e6, b24, ("hott", "hitec", "sfhss"), 92)}

pk = A.packet(0x3D1A22, 6, A.traffic_payload(28.61, 77.21, 250, 11, speed_ms=12.0, track_deg=45.0))
og = O.packet(0x4C2B1A, 3, 30.1, 78.05, 900, 13)
b868 = [dict(t0_us=1500, f_hz=868.2e6, iq=S.adsl_burst(pk, cfo_hz=2e3, fs=FS, rng=rng)[0], bw_hz=0.2e6),
        dict(t0_us=11000, f_hz=868.4e6, iq=S.adsl_burst(pk, cfo_hz=2e3, fs=FS, rng=rng)[0], bw_hz=0.2e6),
        dict(t0_us=20000, f_hz=868.2e6, iq=S.ogn_burst(og, cfo_hz=1e3, fs=FS, rng=rng)[0], bw_hz=0.2e6)]
ref["c868"] = run("c868", 868.3e6, b868, ("adsl", "ogn_like"), 93)

lat, lon = int(round(32.05 * 93206)), int(round(76.71 * 46603))
trk = lat.to_bytes(3, "little") + lon.to_bytes(3, "little") + (1450 | 7 << 12 | 1 << 15).to_bytes(2, "little") + \
    bytes([72, 0x85, 64])
fanet = [[0xF1, "01013412" + trk.hex()], [0x12, "01013412" + trk.hex()], [0x12, "0201341270696c6f74"],
         [0xF1, "02073412" + b"Srikanth".hex()], [0xF1, "c1083412a0" + trk.hex()], [0xF1, "07113412" + trk.hex() + "91"],
         [0xF1, "0400ffff"], [0x2B, "01013412" + trk.hex()]]
ref["fanet"] = [[sw, h, P.identify(sw, bytes.fromhex(h))] for sw, h in fanet]

r = random.Random(5)
rs_cases = []
for nr, k in ((12, 18), (14, 34), (20, 72)):
    c = RS(0x187, 120, 1, nr)
    for ne in (0, 3, nr // 2, nr // 2 + 2):
        data = bytes(r.randrange(256) for _ in range(k))
        cw = bytearray(data + c.encode(data))
        for p in r.sample(range(len(cw)), ne):
            cw[p] ^= r.randrange(1, 256)
        dec, n = c.decode(bytes(cw))
        rs_cases.append([nr, data.hex(), bytes(cw).hex(), dec.hex() if dec is not None else None, n])
ref["rs"] = rs_cases

m = U.payload(0xA1B2C3, 37.5, -122.25, 1200, 14, "N123UA")
s = U.modulate(m)
y = np.r_[np.zeros(700), s * np.exp(2j * np.pi * 8e3 * np.arange(len(s)) / U.FS), np.zeros(900)]
y = (y + (rng.normal(size=len(y)) + 1j * rng.normal(size=len(y))) * 0.3).astype(np.complex64)
y.tofile(OUT / "uat.cf32")
samples = ["00a66ef135445d525a0c0519119021204800",
           "08a66ef1353e2d525fd4050911882aa038101d06b85d440be2a4c2a0000590000000",
           "10a66ef1353ced52614204f911782d00180000000000000000000000000590000000", m.hex()]
ref["uat"] = dict(demod=U.demod(y), parse=[[h, U.parse(bytes.fromhex(h))] for h in samples])
gb = ("ff2048fffffe3135383146414e4c433235385530323952544e363030303030303030000101d1823b483bf2eb11ed0769833b4822"
      "f0eb1105001c00284700c2083d0902000c050478acc5529e0103")
ref["gb46750"] = [["07" + gb, R.parse_odid(bytes.fromhex("07" + gb))], [gb, R.parse_gb46750(bytes.fromhex(gb))],
                  [gb[:4] + "47" + gb[6:], R.parse_gb46750(bytes.fromhex(gb[:4] + "47" + gb[6:]))]]
print(len(ref["uat"]["demod"]), "uat")
(OUT / "cases.json").write_text(json.dumps(ref, indent=1, default=float))
