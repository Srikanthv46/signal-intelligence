#!/usr/bin/env python3
"""Reference vectors for c++/tests/unit/test_rx_voice_fm.cpp, written by the Python receiver
(python/spectrum_intel/receivers/voice_fm): the real PMR446 capture (testdata/voice, DCS D053N)
and a synthetic speech + CTCSS 88.5 Hz transmitter (synth.nbfm, saved as synth.cf32).

    .venv/bin/python c++/tests/vectors/voice_fm/gen_voice_fm.py
"""
import json
import time
from pathlib import Path

import numpy as np

from spectrum_intel.core.dsp import Channelizer
from spectrum_intel.receivers.voice_fm import detector as D
from spectrum_intel.receivers.voice_fm import subaudio as S
from spectrum_intel.receivers.voice_fm.audio import VoiceDemodulator, demodulate
from spectrum_intel.receivers.voice_fm.synth import nbfm

OUT = Path(__file__).resolve().parent
TD = OUT.parents[3] / "testdata"


def js(o):
    if isinstance(o, dict):
        return {k: js(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [js(v) for v in o]
    if isinstance(o, np.generic):
        return o.item()
    return o


def one(x, fs, fc, tag):
    ch = Channelizer(x, fs, fc)
    ref = {"fs": fs, "fc": fc}
    cands = D.candidates(ch)
    ref["candidates"] = js(cands)
    ref["measure"] = []
    for c in cands:
        m = D.measure(ch, c)
        ref["measure"].append(None if m is None else js({k: v for k, v in m.items() if k != "segs_s"}))
    t0 = time.time()
    em = D.detect(x, fs, fc, {"ch": ch, "ppm": 0.0})
    ref["detect_s"] = time.time() - t0
    ref["emissions"] = [js(e.to_dict()) for e in em]
    if em:
        f = {"center_mhz": em[0].features["center_mhz"], "deviation_hz": em[0].features["deviation_hz"],
             "channel_bw_khz": em[0].features["channel_bw_khz"]}
        ref["fields"] = f
        a = demodulate(x, fs, fc, f)
        a.samples.astype("<f4").tofile(OUT / f"{tag}_audio.f32")
        ref["audio"] = {"fs": a.fs, "meta": js(a.meta), "n": len(a.samples)}
        demod = VoiceDemodulator.from_fields(f)
        blocks, metas = [], []
        for k in range(0, len(x), 50000):
            c = demod.process(x[k:k + 50000], fs, fc)
            metas.append(None if c is None else js(c.meta))
            if c is not None:
                blocks.append(c.samples)
        np.concatenate(blocks).astype("<f4").tofile(OUT / f"{tag}_stream.f32")
        ref["stream"] = metas
    (OUT / f"{tag}.json").write_text(json.dumps(ref, indent=1))


if __name__ == "__main__":
    x = np.fromfile(TD / "voice" / "real_pmr446_ch8_dcs.cf32", np.complex64)
    one(x, 250000.0, 446.1e6, "real")
    fs, fc = 250000.0, 446.1e6
    y = nbfm(fs, 0.8, 446.04375e6 - fc, audio="speech", ctcss=88.5, snr_db=25.0, keyed=(0.1, 0.7), seed=3)
    y.tofile(OUT / "synth.cf32")
    one(y, fs, fc, "synth")
    codes = S.DCS_CODES
    (OUT / "dcs.json").write_text(json.dumps({
        "codes": codes, "words": [S.dcs_encode(c) for c in codes],
        "aliases_023": S.dcs_aliases(S.dcs_encode("023")),
        "syndromes": [S.golay_syndrome(w) for w in (0, 1, 0x7FFFFF, 0x123456, 0x400C75)]}, indent=1))
