#!/usr/bin/env python3
"""Bench: GSM FCCH search over a CFO grid (no receiver ppm needed) on a SigMF capture, every ARFCN in band.

    .venv/bin/python bench/live_lo/gsm_fcch_cfo.py testdata/live/LIVE_940MHz [--max-ppm 30] [--min-db 3]
Stage 1: STFT of the ARFCN channel (window = 0.4 ms, zero-padded) -> per frame, the fraction of power in the
strongest bin within F_FCCH +- max CFO ("purity"). Stage 2: the library's find_fcch refined at that tone. Then SCH
one frame later (library demod_sch: training match, Viterbi, CRC-10). Prints ppm = (tone - 67.708 kHz) / f.
"""
import argparse, json, sys
from pathlib import Path
import numpy as np

from spectrum_intel.core.dsp import Channelizer
from spectrum_intel.receivers.gsm.sch import OFS, SPS, F_FCCH, crc10_ok, demod_sch, find_fcch, sch_fields, viterbi
from spectrum_intel.receivers.gsm.sch import _arfcn_freq

WIN = int(0.4e-3 * OFS)


def coarse_fcch(y, max_cfo, nfft=2048, hop=WIN // 4):
    """Best (purity, cfo, frame start) of a pure tone near +F_FCCH over STFT frames of length WIN."""
    m = (len(y) - WIN) // hop
    if m <= 0:
        return None
    idx = np.arange(WIN)[None, :] + hop * np.arange(m)[:, None]
    fr = y[idx]
    X = np.abs(np.fft.fft(fr, nfft, axis=1)) ** 2
    f = np.fft.fftfreq(nfft, 1 / OFS)
    sel = np.abs(f - F_FCCH) <= max_cfo
    e = np.sum(np.abs(fr) ** 2, axis=1) * WIN + 1e-30
    pur = X[:, sel] / e[:, None]
    k = np.unravel_index(np.argmax(pur), pur.shape)
    j = int(np.flatnonzero(sel)[k[1]])
    a, b, c = X[k[0], j - 1], X[k[0], j], X[k[0], (j + 1) % nfft]          # parabolic peak (log power)
    la, lb, lc = np.log(a + 1e-30), np.log(b + 1e-30), np.log(c + 1e-30)
    d = 0.5 * (la - lc) / (la - 2 * lb + lc) if (la - 2 * lb + lc) < 0 else 0.0
    return float(pur[k]), float(f[j] + d * OFS / nfft - F_FCCH), int(k[0] * hop)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("path"); ap.add_argument("--max-ppm", type=float, default=30)
    ap.add_argument("--min-db", type=float, default=3)
    a = ap.parse_args()
    d = Path(a.path); meta = json.loads(next(d.glob("*.sigmf-meta")).read_text())
    fs, fc = meta["global"]["core:sample_rate"], meta["captures"][0]["core:frequency"]
    raw = next(d.glob("*.sigmf-data"))
    x = (np.fromfile(raw, dtype="<f4") if meta["global"]["core:datatype"] == "cf32_le"
         else np.fromfile(raw, dtype="<i2").astype(np.float32)).view(np.complex64)
    ch = Channelizer(x, fs, fc); del x
    rows = []
    for f in np.arange(np.ceil((fc - 27.5e6) / 200e3) * 200e3, fc + 27.5e6, 200e3):
        arfcn, fr = _arfcn_freq(f)
        if arfcn is None:
            continue
        y, _ = ch.slice(fr, 110e3, OFS)
        rows.append((arfcn, fr, y, 10 * np.log10(np.mean(np.abs(y) ** 2) + 1e-30)))
    floor = np.percentile([r[3] for r in rows], 10)
    for arfcn, fr, y, p in rows:
        if p - floor < a.min_db:
            continue
        c = coarse_fcch(y, a.max_ppm * 1e-6 * fr + 2e3)
        if c is None or c[0] < 0.3:
            continue
        pur0, cfo0, k0 = c
        y2 = y * np.exp(-2j * np.pi * cfo0 * np.arange(len(y)) / OFS)
        a0 = max(0, k0 - WIN)                                    # refine around the coarse frame only
        fine = find_fcch(y2[a0:k0 + 2 * WIN + 148 * SPS])
        if fine is not None:
            fine = (fine[0] + a0, fine[1], fine[2])
        line = f"ARFCN {arfcn:4d} {fr / 1e6:8.1f} {p - floor:5.1f} dB  coarse pur {pur0:.2f} ppm {cfo0 / fr * 1e6:+6.2f}"
        if fine is None:
            print(line + "  fine: none"); continue
        s0, _, pur = fine                     # find_fcch's lag-1 CFO is biased by in-channel residuals: keep the FFT one
        cfo = cfo0
        y3 = y * np.exp(-2j * np.pi * cfo * np.arange(len(y)) / OFS)
        dd = demod_sch(y3, s0 + int(round(1250 * SPS)))
        sch = "SCH: none"
        if dd is not None:
            bits = viterbi(dd[0]); ok = crc10_ok(bits[:25], bits[25:35])
            sch = f"SCH ts {dd[1]['ts_match']}/64 " + (f"CRC ok BSIC {sch_fields(bits[:25])['ncc']}{sch_fields(bits[:25])['bcc']}" if ok else "CRC fail")
        print(line + f"  fine pur {pur:.2f} t {s0 / OFS * 1e3:6.2f} ms ppm {cfo / fr * 1e6:+6.2f}  {sch}")


if __name__ == "__main__":
    main()
