"""Radio altimeters, 4200-4400 MHz (research, signal gaps S8; research/signal_gaps/05_aviation_radar_jammers.md
section 1). FMCW altimeters sweep 104-177 MHz (ITU-R M.2059), wider than any capture, so a capture sees crossings.

    detect(x, fs, fc) -> list[Emission]      family "radalt", label "radio_altimeter" (at most one per block)

1. Ridge: 128-point Hann FFT frames every 64 samples; per frame the peak bin (parabolic on log power) where it stands
   >= 10 dB over the frame's median.
2. Crossings: ridge frames joined into tracks (<= 4 bins per elapsed frame, gaps of up to 8 frames bridged, so fades
   and noise peaks do not cut a sweep); a track with >= 5 frames, on in >= half of the
   frames it spans, and a linear fit with R^2 >= 0.99 is a crossing (a
   track failing the fit is split once at its frequency extreme: a turnaround inside the capture)
   (time at the capture centre, rate, direction); same-direction pieces within 20 us of each other merge; the
   crossing must cover >= 15 MHz.
3. Sweep: rate = median |slope|; same-sign crossing spacing (in whole periods) -> period; alternating signs -> triangle, turnarounds
   f_top = f_c + r (t_down - t_up) / 2 and f_bot = f_c - r (t_up - t_down) / 2 -> sweep bandwidth and the
   triangle period 2 B / r; one sign only -> sawtooth.
4. Validity: rate 2-1500 MHz/ms. Confidence 0.7, 0.8 when the bandwidth is 95-190 MHz or the period matches an
   M.2059 type (20 ms, 6.67 ms, 7 ms, 1 ms, +-5 %).
Depends on numpy and spectrum_intel.core only.
"""
from __future__ import annotations

import numpy as np

from spectrum_intel.core.emission import Emission

BAND = (4200e6, 4400e6)
NFFT, HOP = 128, 64
RIDGE_DB = 10.0
MAX_STEP_BINS = 4
MIN_RUN = 5
MAX_GAP = 8
MERGE_S = 20e-6
MIN_FILL = 0.5
MIN_SPAN_HZ = 15e6
MIN_R2 = 0.99
RATE_RANGE = (2e9, 1500e9)            # Hz/s: 2-1500 MHz/ms
TABLE_PERIODS_S = (20e-3, 6.667e-3, 7e-3, 1e-3)


def ridge(x: np.ndarray, fs: float):
    """(frame centre times s, ridge frequency Hz relative to fc, ridge on mask)."""
    nf = (len(x) - NFFT) // HOP + 1
    if nf < MIN_RUN:
        return None
    idx = np.arange(nf)[:, None] * HOP + np.arange(NFFT)[None, :]
    F = np.fft.fft(np.asarray(x)[idx].astype(np.complex128) * np.hanning(NFFT)[None, :], axis=1)
    P = np.abs(F) ** 2
    k = np.argmax(P, axis=1)
    rows = np.arange(nf)
    pk = P[rows, k]
    on = 10 * np.log10(pk / (np.median(P, axis=1) + 1e-30) + 1e-30) >= RIDGE_DB
    a = np.log(P[rows, (k - 1) % NFFT] + 1e-30)
    b = np.log(pk + 1e-30)
    c = np.log(P[rows, (k + 1) % NFFT] + 1e-30)
    den = a - 2 * b + c
    frac = np.where(den < 0, 0.5 * (a - c) / np.where(den < 0, den, -1.0), 0.0)
    kk = np.where(k < NFFT // 2, k, k - NFFT).astype(np.float64)
    t = (rows * HOP + NFFT / 2) / fs
    return t, (kk + frac) * (fs / NFFT), on


def _fit(t: np.ndarray, f: np.ndarray, r: list[int]) -> dict | None:
    tt, ff = t[r], f[r]
    tm, fm = float(np.mean(tt)), float(np.mean(ff))
    sxx = float(np.sum((tt - tm) ** 2))
    if sxx <= 0:
        return None
    slope = float(np.sum((tt - tm) * (ff - fm))) / sxx
    if slope == 0:
        return None
    res = ff - (fm + slope * (tt - tm))
    syy = float(np.sum((ff - fm) ** 2))
    r2 = 1 - float(np.sum(res ** 2)) / syy if syy > 0 else 0.0
    return {"t_s": tm - fm / slope, "rate_hz_s": slope, "span_hz": abs(slope) * float(tt[-1] - tt[0]), "r2": r2,
            "f_lo_hz": float(ff.min()), "f_hi_hz": float(ff.max()), "frames": list(r)}


def crossings(t: np.ndarray, f: np.ndarray, on: np.ndarray, fs: float) -> list[dict]:
    """Ridge frames joined into tracks (a frame joins the first open track it continues within MAX_STEP_BINS per
    elapsed frame, up to MAX_GAP frames back; otherwise it opens one), fitted; same-direction tracks crossing f_c
    within MERGE_S are one crossing split by a fade and are refitted together."""
    step = MAX_STEP_BINS * fs / NFFT
    open_, done = [], []
    for i in np.flatnonzero(on):
        i = int(i)
        keep = []
        for r in open_:
            (keep if i - r[-1] <= MAX_GAP else done).append(r)
        open_ = keep
        for r in open_:
            if abs(f[i] - f[r[-1]]) <= step * (i - r[-1]):
                r.append(i)
                break
        else:
            open_.append([i])
    done += open_
    pieces = []
    for r in done:
        if len(r) < MIN_RUN or len(r) < MIN_FILL * (r[-1] - r[0] + 1):
            continue                              # a sweep is on in most frames it spans; chance noise chains are sparse
        c = _fit(t, f, r)
        if c is not None and c["r2"] >= MIN_R2:
            pieces.append(r)
            continue
        fr = f[r]                                 # a turnaround inside the capture: split the V at its extreme
        for k in (int(np.argmax(fr)), int(np.argmin(fr))):
            if MIN_RUN <= k + 1 and len(r) - k >= MIN_RUN:
                pieces += [r[:k + 1], r[k:]]
                break
    out = []
    for r in sorted(pieces, key=lambda r: r[0]):
        c = _fit(t, f, r)
        if c is None or c["r2"] < MIN_R2:
            continue
        if out and (out[-1]["rate_hz_s"] > 0) == (c["rate_hz_s"] > 0) and abs(out[-1]["t_s"] - c["t_s"]) < MERGE_S:
            m = _fit(t, f, sorted(out[-1]["frames"] + r))
            if m is not None and m["r2"] >= MIN_R2:
                out[-1] = m
                continue
        out.append(c)
    return [c for c in out if c["span_hz"] >= MIN_SPAN_HZ]


def detect(x: np.ndarray, fs: float, fc: float) -> list[Emission]:
    if not (fc - fs / 2 < BAND[1] and fc + fs / 2 > BAND[0]) or fs < 20e6:
        return []
    rd = ridge(x, fs)
    if rd is None:
        return []
    cr = crossings(*rd, fs)
    if not cr:
        return []
    rate = float(np.median([abs(c["rate_hz_s"]) for c in cr]))
    if not RATE_RANGE[0] <= rate <= RATE_RANGE[1]:
        return []
    signs = [1 if c["rate_hz_s"] > 0 else -1 for c in cr]
    feats = {"n_crossings": len(cr), "rate_mhz_per_ms": round(rate / 1e9, 3),
             "r2_min": round(min(c["r2"] for c in cr), 4)}
    period = None
    for sg in (1, -1):
        ts = [c["t_s"] for c, s in zip(cr, signs) if s == sg]
        if len(ts) >= 2:                      # spacings are whole periods when a crossing is missed
            d = np.diff(ts)
            period = float(np.median(d / np.maximum(np.round(d / d.min()), 1.0)))
            break
    tops, bots = [], []
    for a, b in zip(cr[:-1], cr[1:]):
        if a["rate_hz_s"] > 0 > b["rate_hz_s"]:
            tops.append(fc + rate * (b["t_s"] - a["t_s"]) / 2)
        elif a["rate_hz_s"] < 0 < b["rate_hz_s"]:
            bots.append(fc - rate * (b["t_s"] - a["t_s"]) / 2)
    if len(set(signs)) == 2:
        feats["shape"] = "triangle"
    elif len(cr) >= 2:
        feats["shape"] = "sawtooth up" if signs[0] > 0 else "sawtooth down"
    else:
        feats["shape"] = "single crossing"
    bw = None
    if tops:
        feats["f_top_mhz"] = round(float(np.median(tops)) / 1e6, 2)
    if bots:
        feats["f_bottom_mhz"] = round(float(np.median(bots)) / 1e6, 2)
    if tops and bots:
        bw = float(np.median(tops)) - float(np.median(bots))
        feats["sweep_bw_mhz"] = round(bw / 1e6, 1)
        if feats["shape"] == "triangle":
            period = 2 * bw / rate
    if period is not None:
        feats["period_ms"] = round(period * 1e3, 3)
    conf = 0.7
    if (bw is not None and 95e6 <= bw <= 190e6) or (
            period is not None and any(abs(period - p) <= 0.05 * p for p in TABLE_PERIODS_S)):
        conf = 0.8
    seen_lo = min(c["f_lo_hz"] for c in cr) + fc
    seen_hi = max(c["f_hi_hz"] for c in cr) + fc
    feats["seen_mhz"] = [round(seen_lo / 1e6, 2), round(seen_hi / 1e6, 2)]
    t0 = float(rd[0][0])
    t1 = float(rd[0][-1])
    return [Emission(family="radalt", label="radio_altimeter", start_us=t0 * 1e6, end_us=t1 * 1e6,
                     f_lo_mhz=seen_lo / 1e6, f_hi_mhz=seen_hi / 1e6, center_mhz=(seen_lo + seen_hi) / 2e6,
                     confidence=conf, features=feats)]
