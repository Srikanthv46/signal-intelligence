#!/usr/bin/env python3
"""Synthetic reference for the COFDM / DVB-T receiver (research; no capture): two DVB-T signals at 12.288 MSps,
2^18 samples (21 ms: one 2K / 8 MHz TPS frame), 10 dB in-band SNR, and python detect()'s features on each.

    nice -n 19 .venv/bin/python -I c++/tests/vectors/cofdm/gen_cofdm.py
"""
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parents[3] / "python"))
from spectrum_intel.receivers.cofdm import detector as D, synth as S   # noqa: E402

FS, N, FC = 12.288e6, 1 << 18, 1250e6
ref = {"fs": FS, "fc": FC, "cases": []}
for i, (B, gi, qam, cr) in enumerate(((8, "1/8", 4, 0), (7, "1/32", 16, 4))):
    rng = np.random.default_rng(40 + i)
    y = S.dvbt_capture_burst(100, B, FS, gi=gi, qam=qam, rng=rng, cr_hp=cr, constellation={4: 0, 16: 1}[qam],
                             cell_id=0xBEEF if i else None)[:N]
    sig = np.zeros(N, np.complex64)
    sig[:len(y)] = y
    sig *= np.exp(2j * np.pi * 0.7e6 / FS * np.arange(N)).astype(np.complex64)
    x = (sig + (rng.standard_normal(N) + 1j * rng.standard_normal(N)) * np.sqrt(FS / (B * 1e6) * 10 ** (-1.0) / 2)).astype(np.complex64)
    x.tofile(OUT / f"dvbt_{i}.cf32")
    em = D.detect(x, FS, FC)
    ref["cases"].append(dict(file=f"dvbt_{i}.cf32", emissions=[dict(label=e.label, f_lo_mhz=e.f_lo_mhz, f_hi_mhz=e.f_hi_mhz,
                                                                    features=e.features) for e in em]))
    print(i, [(e.label, (e.features.get("tps") or {}).get("constellation")) for e in em])
bits = [0] + [int(b) for b in "1011010" + "0011010111101110" + "010111" + "00" + "00" + "000" + "000" + "000" + "10" + "00" + "0" * 14] * 1
ref["tps_bits"] = bits
ref["tps_out"] = D.tps_decode(np.array(bits, np.uint8))
(OUT / "ref.json").write_text(json.dumps(ref, indent=1, default=float))
