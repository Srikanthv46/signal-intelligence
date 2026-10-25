#!/usr/bin/env python3
"""Synthetic reference for the colour-standard classifier (analog_vtx detector.chroma_standard; research:
research/drone_transmitters/01_analog_video.md). No capture is used: composite video of each standard from
receivers/analog_vtx/synth.py (5 ms at 20.48 MSps, camera clock offset, white noise) stands in for the FM
discriminator output; float32 on disk, and the reference is computed on the float32-rounded samples.

    nice -n 19 .venv/bin/python -I c++/tests/vectors/analog_vtx/gen_chroma.py
"""
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parents[3] / "python"))
from spectrum_intel.receivers.analog_vtx import detector as D, synth as S   # noqa: E402

FS = 20.48e6
N = int(FS * 5e-3)
ref = {"fs": FS, "cases": []}
for i, (std, ppm, color) in enumerate((("ntsc", 30, True), ("pal_m", 15, True), ("pal", -40, True), ("pal_n", -20, True),
                                       ("pal", 0, False), ("secam", 0, True), ("secam", 25, True))):
    rng = np.random.default_rng(100 + i)
    v = S.composite(N, FS, std, ppm, rng, t0_lines=3.0, color=color)
    v = (v + 0.02 * rng.standard_normal(N)).astype("<f4")
    v.tofile(OUT / f"chroma_{i}.f32")
    rate = S.STANDARDS[std]["fl"] * (1 + ppm * 1e-6)
    t0 = (np.ceil(3.0) - 3.0) / rate                    # first sync leading edge in the block
    line_std = "ntsc" if S.STANDARDS[std]["lines"] == 525 else "pal"
    r = D.chroma_standard(v.astype(np.float64), FS, rate, t0, line_std)
    ref["cases"].append(dict(file=f"chroma_{i}.f32", std=std, line_std=line_std, rate=rate, t0=t0,
                             out=json.loads(json.dumps(r, default=float))))
    print(std, color, r.get("color_standard"), r.get("chroma_fsc_hz"), r.get("secam_test"))
(OUT / "chroma.json").write_text(json.dumps(ref, indent=1))

# audio subcarrier (research, PLAN A5): composite + a 6.5 MHz FM subcarrier carrying a 440 Hz tone + noise, as the
# FM discriminator would output it (20.48 MSps, 34 ms); reference = detector.audio_features on the float32 samples
rng = np.random.default_rng(7)
na = int(FS * 34e-3)
v = S.composite(na, FS, "pal", 0.0, rng, t0_lines=1.0)
t = np.arange(na) / FS
tone = 0.8 * np.sin(2 * np.pi * 440.0 * t) + 0.2 * rng.standard_normal(na)
v = v + 0.05 * np.cos(2 * np.pi * 6.5e6 * t + 2 * np.pi * 40e3 * np.cumsum(tone) / FS)
v = (v + 0.01 * rng.standard_normal(na)).astype("<f4")
v.tofile(OUT / "audio.f32")
au = D.audio_features(v.astype(np.float64), FS, 6.5e6)
print("audio", au)
(OUT / "audio.json").write_text(json.dumps({"fs": FS, "fa": 6.5e6, "out": au}, indent=1, default=float))
