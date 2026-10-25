#!/usr/bin/env python3
"""Synthetic reference for the blind FSK frame structure (rc_link/fsk_frame; research, PLAN F1): a 7.68 MSps capture
(2^18 samples, fc 2440 MHz) with three 250 kBd packets of one unknown sync word and one 250 kBd packet with the CC2500
D16 sync at 18 dB; rc_link detect's emissions (label, features.frame) and the plugin's fsk_frame field.

    nice -n 19 .venv/bin/python -I c++/tests/vectors/fsk_frame/gen_fsk_frame.py
"""
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parents[3] / "python"))
from spectrum_intel.receivers.rc_link import detector as D, synth as S   # noqa: E402
from spectrum_intel.receivers.rc_link.plugin import PLUGINS            # noqa: E402

FS, FC = 7.68e6, 2440e6
rng = np.random.default_rng(51)


def bits(sync, n_pay, pre=32):
    return np.r_[np.arange(pre) % 2 == 0, [(sync >> (31 - i)) & 1 for i in range(32)], rng.integers(0, 2, n_pay),
                 np.zeros(8)].astype(np.uint8)


x = S.noise_capture(1 << 18, rng)
bursts = []
for k, f in enumerate((2437.6e6, 2439.2e6, 2441.3e6)):
    iq, _ = S.gfsk_burst(250e3, 125e3, bits=bits(0x9A3C5E71, 240), fs=FS, rng=rng)
    bursts.append(dict(t0_us=1500 + 5000 * k, f_hz=f, iq=iq, bw_hz=600e3))
iq, _ = S.gfsk_burst(250e3, 125e3, bits=bits(0xD391D391, 240), fs=FS, rng=rng)    # CC2500 sync at a non-FrSky rate
bursts.append(dict(t0_us=18000, f_hz=2442.5e6, iq=iq, bw_hz=600e3))
y, _ = S.embed(x, FS, FC, bursts, 18)
y = y.astype(np.complex64)
y.tofile(OUT / "capture.cf32")
em = D.detect(y, FS, FC)
fr = [e for e in em if e.features.get("frame")]
ref = dict(fs=FS, fc=FC, labels=[e.label for e in fr], frames=[e.features["frame"] for e in fr],
           field=PLUGINS[0].summarise(fr, None)["fields"].get("fsk_frame"))
print(ref["labels"], ref["field"])
(OUT / "cases.json").write_text(json.dumps(ref, indent=1, default=float))
