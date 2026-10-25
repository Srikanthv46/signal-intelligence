"""Aeronautical navigation aids: ILS glide slope, ILS localizer, VOR (research, signal gaps S9 / S21;
research/signal_gaps/05_aviation_radar_jammers.md section 2, 06_vhf_below_300.md section 9).

    detect(ch, ppm=0.0) -> list[Emission]      family "navaids", label ils_gs | ils_loc | vor

1. Channels: the 40 ILS glide-path frequencies (329.15-335.00 MHz) and their paired localizers (108.10-111.95 MHz,
   odd tenths), the VOR channels (108-117.95 MHz, 50 kHz; even tenths below 112 MHz). Each channel inside the capture
   is checked on the block FFT: a carrier line >= 15 dB over the median of +-35 kHz around it (receiver ppm applied);
   a second line 4-32 kHz away (>= -15 dB) marks a two-frequency (capture-effect) system.
2. ILS: the strongest carrier's envelope (+-1.5 kHz at 6 kSps) over whole 1/30 s windows (90 and 150 Hz are both
   harmonics of 30 Hz), least-squares fit of [1, 90 s Hz, 150 s Hz] with a common tone scale s in 0.97-1.03 (the
   ratio is locked to 5/3): m90, m150 (depths), DDM = m90 - m150, SDM = m90 + m150.
3. VOR: the envelope (+-12.5 kHz at 32 kSps) over whole 30 Hz periods; 30 Hz AM by least squares; the 9960 Hz
   subcarrier cut out of the envelope's FFT (+-700 Hz), its instantaneous frequency fitted with 30 Hz -> deviation and
   phase; bearing (radial from the station) = AM phase - FM phase.
4. Validity: ILS needs both tones >= 5 sigma, s inside the search range and the SDM of its kind (glide path
   0.6-1.0, localizer 0.25-0.65); VOR needs the subcarrier >= 10 dB over its surroundings, AM depth 0.15-0.45 and
   deviation 300-700 Hz. Short of that the emission is "-like" (confidence 0.55).
Depends on numpy and spectrum_intel.core only.
"""
from __future__ import annotations

import numpy as np

from spectrum_intel.core.dsp import Channelizer
from spectrum_intel.core.emission import Emission

# (localizer MHz, glide path MHz): ICAO Annex 10 Vol I pairing (research file 05, table 2.2)
ILS_PAIRS = [
    (108.10, 334.70), (108.15, 334.55), (108.30, 334.10), (108.35, 333.95), (108.50, 329.90), (108.55, 329.75),
    (108.70, 330.50), (108.75, 330.35), (108.90, 329.30), (108.95, 329.15), (109.10, 331.40), (109.15, 331.25),
    (109.30, 332.00), (109.35, 331.85), (109.50, 332.60), (109.55, 332.45), (109.70, 333.20), (109.75, 333.05),
    (109.90, 333.80), (109.95, 333.65), (110.10, 334.40), (110.15, 334.25), (110.30, 335.00), (110.35, 334.85),
    (110.50, 329.60), (110.55, 329.45), (110.70, 330.20), (110.75, 330.05), (110.90, 330.80), (110.95, 330.65),
    (111.10, 331.70), (111.15, 331.55), (111.30, 332.30), (111.35, 332.15), (111.50, 332.90), (111.55, 332.75),
    (111.70, 333.50), (111.75, 333.35), (111.90, 331.10), (111.95, 330.95)]
LINE_DB = 15.0
SEARCH_HZ = 35e3
TWO_FREQ_HZ = (4e3, 32e3)
SCALES = np.linspace(0.97, 1.03, 61)
ILS_FS, VOR_FS = 6e3, 32e3
SUB_HZ, SUB_BW = 9960.0, 700.0


def vor_channels() -> list[float]:
    """VOR frequencies (MHz): 108.00-117.95 on 50 kHz, but below 112 MHz only even tenths (+0.05)."""
    out = []
    for k in range(200):
        f = round(108.0 + 0.05 * k, 2)
        if f < 112.0 and (int(round(f * 100)) // 10) % 2:
            continue
        out.append(f)
    return out


def channels() -> list[tuple[str, float, float | None]]:
    """(kind, MHz, paired MHz or None) of every channel."""
    return ([("ils_gs", g, l) for l, g in ILS_PAIRS] + [("ils_loc", l, g) for l, g in ILS_PAIRS]
            + [("vor", f, None) for f in vor_channels()])


def _lines(ch: Channelizer, f_meas: float) -> list[tuple[float, float]] | None:
    """Carrier lines (Hz, dB over the window median) within +-SEARCH_HZ of f_meas, strongest first."""
    k = ch._bins(f_meas - ch.fc, SEARCH_HZ)
    if len(k) < 64:
        return None
    p = np.abs(ch.X[k]).astype(np.float64) ** 2
    f = ch.f[k] + ch.fc
    o = np.argsort(f)
    p, f = p[o], f[o]
    floor = float(np.median(p)) + 1e-30
    db = 10 * np.log10(p / floor + 1e-30)
    out = []
    for i in np.argsort(-db, kind="stable"):
        if db[i] < LINE_DB or len(out) >= 2:
            break
        if all(abs(f[i] - g) >= 2e3 for g, _ in out):
            out.append((float(f[i]), float(db[i])))
    return out


def _fit(a: np.ndarray, fs: float, tones: tuple[float, ...]) -> tuple[np.ndarray, float]:
    """Least squares of a on [1, cos, sin per tone]: coefficients and residual power."""
    t = np.arange(len(a)) / fs
    cols = [np.ones(len(a))]
    for f in tones:
        cols += [np.cos(2 * np.pi * f * t), np.sin(2 * np.pi * f * t)]
    H = np.stack(cols, 1)
    c = np.linalg.solve(H.T @ H, H.T @ a)
    r = a - H @ c
    return c, float(np.mean(r ** 2))


def _whole(n: int, fs: float, rate: float = 30.0) -> int:
    """Samples in the largest whole number of 1/rate periods that fit in n samples."""
    per = fs / rate
    k = int(n / per + 1e-6)
    return int(round(k * per)) if k >= 1 else 0


def measure_ils(ch: Channelizer, f_line: float) -> dict | None:
    y, fs = ch.slice(f_line, 1.5e3, ILS_FS)
    n = _whole(len(y), fs)
    if n < 64:
        return None
    a = np.abs(y[:n]).astype(np.float64)
    res = []
    for s in SCALES:
        res.append(_fit(a, fs, (90.0 * s, 150.0 * s))[1])
    res = np.array(res)
    i = int(np.argmin(res))
    s = float(SCALES[i])
    step = float(SCALES[1] - SCALES[0])
    if 0 < i < len(SCALES) - 1:                       # parabola through the three best
        r0, r1, r2 = res[i - 1], res[i], res[i + 1]
        den = r0 - 2 * r1 + r2
        if den > 0:
            s += float(0.5 * (r0 - r2) / den * step)
    c, rp = _fit(a, fs, (90.0 * s, 150.0 * s))
    a0 = float(c[0])
    m90, m150 = float(np.hypot(c[1], c[2]) / a0), float(np.hypot(c[3], c[4]) / a0)
    snr_env = a0 ** 2 / (rp + 1e-30)
    sigma = float(np.sqrt(2.0 / (n * snr_env)))
    return {"m90": m90, "m150": m150, "ddm": m90 - m150, "sdm": m90 + m150, "tone_scale": s, "sigma": sigma,
            "edge": i in (0, len(SCALES) - 1), "window_s": n / fs}


def measure_vor(ch: Channelizer, f_line: float) -> dict | None:
    y, fs = ch.slice(f_line, 12.5e3, VOR_FS)
    n = _whole(len(y), fs)
    if n < 256:
        return None
    a = np.abs(y[:n]).astype(np.float64)
    a = a / float(np.mean(a))
    c, _ = _fit(a, fs, (30.0,))
    am_depth = float(np.hypot(c[1], c[2]))
    am_phase = float(np.arctan2(c[2], c[1]))
    A = np.fft.fft(a)
    fr = np.fft.fftfreq(n, 1 / fs)
    sub = np.abs(fr - SUB_HZ) <= SUB_BW
    side = (np.abs(fr - SUB_HZ) > 1.5 * SUB_BW) & (np.abs(fr - SUB_HZ) <= 3 * SUB_BW)
    if sub.sum() < 8 or side.sum() < 8:
        return None
    sub_db = float(10 * np.log10(np.mean(np.abs(A[sub]) ** 2) / (np.mean(np.abs(A[side]) ** 2) + 1e-30) + 1e-30))
    Z = np.zeros(n, complex)
    k0 = int(round(SUB_HZ * n / fs))
    idx = np.flatnonzero(sub)
    Z[(idx - k0) % n] = A[idx]
    z = np.fft.ifft(Z)
    fi = np.angle(z[1:] * np.conj(z[:-1])) * fs / (2 * np.pi)
    t = (np.arange(len(fi)) + 0.5) / fs
    H = np.stack([np.ones(len(fi)), np.cos(2 * np.pi * 30 * t), np.sin(2 * np.pi * 30 * t)], 1)
    cf = np.linalg.solve(H.T @ H, H.T @ fi)
    dev = float(np.hypot(cf[1], cf[2]))
    fm_phase = float(np.arctan2(cf[2], cf[1]))
    bearing = float(np.degrees(am_phase - fm_phase) % 360.0)
    return {"am_depth": am_depth, "deviation_hz": dev, "subcarrier_db": sub_db, "bearing_deg": bearing,
            "subcarrier_offset_hz": float(cf[0]), "window_s": n / fs}


def detect(ch: Channelizer, ppm: float = 0.0) -> list[Emission]:
    lo, hi = ch.fc - 0.45 * ch.fs, ch.fc + 0.45 * ch.fs
    out = []
    dur_us = ch.L / ch.fs * 1e6
    for kind, f_mhz, pair in channels():
        f_nom = f_mhz * 1e6
        f_meas = f_nom / (1 + ppm * 1e-6)
        if not (lo + SEARCH_HZ < f_meas < hi - SEARCH_HZ):
            continue
        ln = _lines(ch, f_meas)
        if not ln:
            continue
        f0, db0 = ln[0]
        two = None
        if len(ln) > 1 and TWO_FREQ_HZ[0] <= abs(ln[1][0] - f0) <= TWO_FREQ_HZ[1] and ln[1][1] >= db0 - 15:
            two = round(abs(ln[1][0] - f0))
        feats = {"kind": kind, "channel_mhz": f_mhz, "carrier_mhz": round(f0 / 1e6, 6),
                 "offset_hz": round(f0 * (1 + ppm * 1e-6) - f_nom), "line_db": round(db0, 1),
                 "receiver_ppm_applied": ppm}
        if two:
            feats["two_frequency_hz"] = two
        if kind == "vor":
            m = measure_vor(ch, f0)
            if m is None:
                continue
            ok = m["subcarrier_db"] >= 10 and 0.15 <= m["am_depth"] <= 0.45 and 300 <= m["deviation_hz"] <= 700
            if m["subcarrier_db"] < 6:
                continue                                  # no 9960 Hz subcarrier: not a VOR
            feats.update({"radial_deg": round(m["bearing_deg"], 1), "am_depth": round(m["am_depth"], 3),
                          "deviation_hz": round(m["deviation_hz"]), "subcarrier_db": round(m["subcarrier_db"], 1),
                          "window_s": round(m["window_s"], 4)})
            half = 12.5e3
        else:
            m = measure_ils(ch, f0)
            if m is None:
                continue
            tones = min(m["m90"], m["m150"]) >= 5 * m["sigma"]
            if not tones and max(m["m90"], m["m150"]) < 5 * m["sigma"]:
                continue                                  # a plain carrier on the grid: not ILS
            sdm_ok = (0.6 <= m["sdm"] <= 1.0) if kind == "ils_gs" else (0.25 <= m["sdm"] <= 0.65)
            ok = tones and sdm_ok and not m["edge"]
            side = (("above path" if m["ddm"] > 0 else "below path") if kind == "ils_gs"
                    else ("left of course" if m["ddm"] > 0 else "right of course"))
            feats.update({"paired_mhz": pair, "m90": round(m["m90"], 4), "m150": round(m["m150"], 4),
                          "ddm": round(m["ddm"], 4), "sdm": round(m["sdm"], 4), "side": side,
                          "tone_scale": round(m["tone_scale"], 4), "depth_sigma": round(m["sigma"], 5),
                          "window_s": round(m["window_s"], 4)})
            half = (two or 0) / 2 + 2e3
        feats["verified_structure"] = bool(ok)
        out.append(Emission(family="navaids", label=kind, start_us=0.0, end_us=dur_us,
                            f_lo_mhz=(f0 - half) / 1e6 if not two else (min(f0, ln[1][0]) - 2e3) / 1e6,
                            f_hi_mhz=(f0 + half) / 1e6 if not two else (max(f0, ln[1][0]) + 2e3) / 1e6,
                            center_mhz=f0 / 1e6, snr_db=round(db0, 1), confidence=0.8 if ok else 0.55,
                            features=feats))
    return out
