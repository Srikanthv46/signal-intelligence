#!/usr/bin/env python3
"""Synthetic reference for XN297 / Bayang (rc_link/xn297; research, PLAN R8): a 7.68 MSps capture (2^18 samples,
fc 2440 MHz) with a Bayang bind packet, two data packets and one XN297 packet of another payload protocol at 22 dB;
rc_link detect's XN297 emissions (label, decoded) and the plugin's xn297 / model fields; plus the scramble table and
CRC cases.

    nice -n 19 .venv/bin/python -I c++/tests/vectors/xn297/gen_xn297.py
"""
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parents[3] / "python"))
from spectrum_intel.receivers.rc_link import detector as D, synth as S, xn297 as X   # noqa: E402
from spectrum_intel.receivers.rc_link.plugin import PLUGINS                        # noqa: E402

FS, FC = 7.68e6, 2440e6
TX = bytes([0x11, 0x22, 0x33, 0x44, 0x55])
rng = np.random.default_rng(71)
x = S.noise_capture(1 << 18, rng)
pls = [S.bayang_payload("bind", TX), S.bayang_payload("data", TX, (700, 512, 100, 300)),
       S.bayang_payload("data", TX, (10, 1000, 900, 512)), bytes(rng.integers(0, 256, 15).tolist())]
bursts = []
for k, (f, p) in enumerate(zip((2437.5e6, 2439.0e6, 2441.5e6, 2442.8e6), pls)):
    iq, _ = S.xn297_burst(TX, p, cfo_hz=8e3, fs=FS, rng=rng)
    bursts.append(dict(t0_us=1000 + 2000 * k, f_hz=f, iq=iq, bw_hz=1.2e6))
y, _ = S.embed(x, FS, FC, bursts, 22)
y = y.astype(np.complex64)
y.tofile(OUT / "capture.cf32")
em = [e for e in D.detect(y, FS, FC) if e.label in ("bayang", "xn297_like")]
fl = PLUGINS[0].summarise([e for e in em if e.decoded], None)["fields"]
data = [bytes(rng.integers(0, 256, n).tolist()) for n in (0, 1, 20)]
# standalone try_decode at low SNR (PLAN L1c soft path): 8 MHz, noise filtered to +-0.75 MHz
sa = []
for k, snr in enumerate((6, 4, 2)):
    iq, _ = S.xn297_burst(TX, S.bayang_payload("data", TX, (100 * k, 512, 300, 700)), fs=8e6, rng=rng)
    iq = np.r_[np.zeros(200), iq, np.zeros(200)]
    w = iq + (rng.normal(size=len(iq)) + 1j * rng.normal(size=len(iq))) * np.sqrt(10 ** (-snr / 10) / 2)
    W = np.fft.fft(w); W[np.abs(np.fft.fftfreq(len(w), 1 / 8e6)) > 0.75e6] = 0
    w = np.fft.ifft(W).astype(np.complex64)
    sa.append(dict(iq=[w.real.tolist(), w.imag.tolist()], out=X.try_decode(w, 8e6)))
ref = dict(fs=FS, fc=FC, standalone=sa, labels=[e.label for e in em], decoded=[e.decoded for e in em], field=fl.get("xn297"),
           model=fl.get("model"), scramble=X.SCRAMBLE, crc=[[d.hex(), X.crc16(d)] for d in data])
print(ref["labels"], ref["field"])
(OUT / "cases.json").write_text(json.dumps(ref, indent=1, default=float))
