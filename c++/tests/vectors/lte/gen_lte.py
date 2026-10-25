"""Reference vectors for c++/tests/unit/test_rx_lte.cpp, written by the Python receiver
(spectrum_intel.receivers.lte.cellsearch) on a synthetic LTE FDD downlink (PSS/SSS + QPSK fill).

    .venv/bin/python c++/tests/vectors/lte/gen_lte.py   -> lte.json, lte_iq.cf32 next to this script
"""
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from spectrum_intel.core import seq3gpp as V
from spectrum_intel.core.dsp import Channelizer, IQContext
from spectrum_intel.receivers.lte import cellsearch as cs


def synth_lte(rng, nfft, n_slots, nid1, nid2, half_fill, slot0=0):
    """FDD normal-CP downlink at nfft * 15 kHz: PSS (symbol 6) / SSS (symbol 5) in slots 0 and 10."""
    cp0, cp = 20 * nfft // 256, 18 * nfft // 256
    sc = V._SC62 % nfft
    fill = np.r_[np.arange(-half_fill, 0), np.arange(1, half_fill + 1)] % nfft
    out = []
    for s in range(slot0, slot0 + n_slots):
        for l in range(7):
            X = np.zeros(nfft, complex)
            X[fill] = (rng.choice([-1, 1], len(fill)) + 1j * rng.choice([-1, 1], len(fill))) / np.sqrt(2)
            if s % 10 == 0 and l in (5, 6):
                X[np.r_[np.arange(-36, 0), np.arange(1, 37)] % nfft] = 0
                X[sc] = V.lte_pss_freq(cs.ROOTS[nid2]) if l == 6 else V.lte_sss_freq(nid1, nid2, 0 if s % 20 == 0 else 5)
            x = np.fft.ifft(X) * np.sqrt(nfft)
            c = cp0 if l == 0 else cp
            out.append(np.r_[x[-c:], x])
    return np.concatenate(out)


def cells_json(cells):
    out = []
    for c in cells:
        s = c["sss"]
        out.append({"raster_hz": c["raster_hz"], "carrier_hz": c["carrier_hz"], "pss_z": c["pss_z"],
                    "pss_isolation": c["pss_isolation"], "nid2": c["nid2"], "inst": [list(map(int, x)) for x in c["inst"]],
                    "sss_ok": c["sss_ok"], "phase": c["phase"], "err_hz": c["err_hz"],
                    "sss": None if s is None else {k: s[k] for k in ("duplex", "cp", "nid1", "parity", "n", "z", "z2")}})
    return out


rng = np.random.default_rng(11)
fs, fc, f_car = 3.84e6, 944.9e6, 945.0e6 + 300.0
x = synth_lte(rng, 256, 62, 101, 2, 36, slot0=3)            # 31 ms
x = x * np.exp(2j * np.pi * (f_car - fc) * np.arange(len(x)) / fs)
x = x + (rng.normal(size=len(x)) + 1j * rng.normal(size=len(x))) * 0.35
iq = x.astype(np.complex64)
here = Path(__file__).parent
iq.tofile(here / "lte_iq.cf32")

ctx = IQContext(iq, fs, fc)
sig = SimpleNamespace(f0=944.3e6, f1=945.75e6, fc=0.5 * (944.3e6 + 945.75e6))
out = {"fs": fs, "fc": fc, "sig": [sig.f0, sig.f1], "cells": cells_json(cs.search_cells(ctx, sig, 0.0))}
w = cs.wide_ppm_search(ctx.ch, [945.0e6], span=20e3, step=5e3)
out["wide"] = w
tpl = np.conj(np.fft.fft(cs._templates(), 512, axis=1))
out["tpl_512"] = [[float(v.real), float(v.imag)] for v in tpl.ravel()]

# bandwidth_from_energy on a synthetic coarse spectrogram
P = (rng.random((40, 200)) * 0.2 + 1.0).astype(np.float32)
f = np.linspace(940e6, 950e6, 200)
P[:, (f > 942.3e6) & (f < 947.7e6)] *= 30
det = SimpleNamespace(coarse=SimpleNamespace(P=P, f=f), floor_coarse_db=np.full(200, 8.0))
out["bw_P"] = P.ravel().tolist()
out["bw_f"] = f.tolist()
out["bw"] = [{"fcar": fcar, "bw": cs.bandwidth_from_energy(det, fcar)} for fcar in (945.0e6, 944.0e6, 940.5e6)]
(here / "lte.json").write_text(json.dumps(out))
print("wrote lte.json:", [(c["nid2"], c["sss"] and c["sss"]["nid1"], c["sss_ok"]) for c in out["cells"]], out["bw"])
