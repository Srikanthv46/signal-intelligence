#!/usr/bin/env python3
"""Reference vectors for c++/tests/unit/test_rx_analog_vtx.cpp, written by the Python receiver
(python/spectrum_intel/receivers/analog_vtx) on testdata/Analog_20260804-161815 (NTSC VTX) and
testdata/DJI_Mini_2_20251119-165139 (5.78 GHz, no VTX).

    .venv/bin/python c++/tests/vectors/analog_vtx/gen_analog_vtx.py
"""
import json
import time
from pathlib import Path

import numpy as np

from spectrum_intel.core.emission import load_capture
from spectrum_intel.receivers.analog_vtx import detector as D
from spectrum_intel.receivers.analog_vtx.video import FieldDemodulator

OUT = Path(__file__).resolve().parent
REPO = OUT.parents[3]
TD = REPO / "testdata"


def js(o):
    if isinstance(o, dict):
        return {k: js(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [js(v) for v in o]
    if isinstance(o, np.generic):
        return o.item()
    return o


def one(meta, tag):
    x, fs, fc, _ = load_capture(str(meta))
    x = np.asarray(x[:D.MAX_SAMPLES], np.complex64)
    q, P, t, f = D.tf_power(x, fs)
    ref = {"fs": fs, "fc": fc, "n_frames": int(len(t)), "nfft": int(len(f))}
    P[:1000].astype("<f4").tofile(OUT / f"{tag}_P1000.f32")
    q[:1000].astype("<f4").tofile(OUT / f"{tag}_q1000.f32")
    q.std(0).astype("<f4").tofile(OUT / f"{tag}_qsd.f32")
    best = D.scan(q, t, f)
    prof = np.concatenate([best["ppm_profile"][s][1] for s in ("pal", "ntsc")])
    prof.astype("<f8").tofile(OUT / f"{tag}_profile.f64")
    ref["best"] = {"score": best["score"], "std": best["std"], "ppm": best["ppm"], "b0": best["b0"], "k0": best["k0"]}
    F, cnt = D.fold(q, t, best["rate"])
    F.astype("<f8").tofile(OUT / f"{tag}_fold_best.f64")
    t0 = time.time()
    em = D.detect(x, fs, fc, {})
    ref["detect_s"] = time.time() - t0
    ref["emissions"] = [js(e.to_dict()) for e in em]
    if em:
        img = D.video_frame(x, fs, fc, em[0])
        ref["video_shape"] = list(img.shape)
        img.tofile(OUT / f"{tag}_video.u8")
        demod = FieldDemodulator.from_emission(em[0])
        frames = []
        for k, (a, b) in enumerate(((0, len(x)), (500000, len(x)))):
            fr = demod.process(x[a:b], fs, fc)
            frames.append(None if fr is None else {"shape": list(fr.image.shape), "meta": js(fr.meta)})
            if fr is not None:
                fr.image.tofile(OUT / f"{tag}_field{k}.u8")
                if fr.rgb is not None:
                    fr.rgb.tofile(OUT / f"{tag}_field{k}_rgb.u8")
                frames[-1]["rgb"] = fr.rgb is not None
        ref["fields"] = frames
    (OUT / f"{tag}.json").write_text(json.dumps(ref, indent=1))


if __name__ == "__main__":
    one(TD / "Analog_20260804-161815" / "Analog_20260804-161815.sigmf-meta", "analog")
    one(next(TD.rglob("D_20251119-165139.sigmf-meta")), "dji")
