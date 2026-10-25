"""315 / 433 MHz short-range devices (research, signal gaps S6; research/signal_gaps/04_short_range.md section 1):
OOK fixed-code remotes (EV1527 / PT2262), KeeLoq code-hopping remotes (clear fields only), FSK TPMS (Continental /
VDO, Toyota), and unrecognised OOK / FSK devices.

    detect(ch) -> list[Emission]     family "ism_devices", label tpms | ook_fixed_code | keeloq | ook_device | fsk_device

1. Candidates: per band (314-316, 433.05-434.79 MHz) inside the capture, the band at >= 2 MSps cut into 256-point
   Hann frames; per frame the bin with the largest power over that bin's median across frames, kept when >= 12 dB;
   the kept frequencies sorted and split at gaps > 100 kHz (an FSK pair stays one device); a group of >= 3 frames
   is a candidate at its median frequency; the 4 with the most frames.
2. Per candidate: +-100 kHz at 1 MSps; envelope = 8-sample moving average of |y|; on when above the geometric mean
   of its 10th and 99.5th percentiles (which must be >= 12 dB apart); runs shorter than 30 us flipped (debounce);
   runs -> pulses (start, mark, space) in us.
3. OOK: EV1527 / PT2262 and KeeLoq on the pulse list (ook.py).
4. FSK: a mark >= 2 ms whose discriminator (8-sample average, centred mid-way between its 10th / 90th percentiles)
   is two-level: the median |f| (the deviation) >= 10 kHz and >= 60 % of the samples within 0.6-1.4 of it (a chirp
   or voice is not); sign runs rounded to 52 us symbols; TPMS frames on the symbols, both polarities (tpms.py).
5. Nothing decoded: >= 16 pulses with marks 50-2000 us -> an OOK device; an FSK mark -> an FSK device (0.55, "-like").
Depends on numpy and spectrum_intel.core only.
"""
from __future__ import annotations

import numpy as np

from spectrum_intel.core.emission import Emission

from . import ook as O
from . import tpms as T

BANDS = ((314.0e6, 316.0e6), (433.05e6, 434.79e6))
SCAN_FS = 2e6
NSTFT = 256
ON_DB = 12.0
CLUSTER_HZ = 100e3
MIN_FRAMES = 3
MAX_CAND = 4
SLICE_HALF = 100e3
SLICE_FS = 1e6
SMOOTH = 8
MIN_RUN_US = 30.0
ENV_DB = 12.0
FSK_MIN_US = 2000.0
FSK_DEV_HZ = 10e3
FSK_FRAC = 0.6
SYMBOL_US = 52.0


def candidates(ch) -> list[tuple[float, int, float]]:
    """(centre Hz, frames, peak dB over the bin median) of bursty carriers in the bands, most frames first."""
    out = []
    for lo, hi in BANDS:
        a, b = max(lo, ch.fc - 0.45 * ch.fs), min(hi, ch.fc + 0.45 * ch.fs)
        if b - a < 50e3:
            continue
        mid, half = (a + b) / 2, (b - a) / 2
        y, fs = ch.slice(mid, half, max(SCAN_FS, b - a))
        nf = len(y) // NSTFT
        if nf < 8:
            continue
        Y = np.fft.fft(y[:nf * NSTFT].astype(np.complex128).reshape(nf, NSTFT) * np.hanning(NSTFT)[None, :], axis=1)
        P = np.abs(Y) ** 2
        fr = np.fft.fftfreq(NSTFT, 1 / fs)
        R = P / (np.median(P, axis=0)[None, :] + 1e-30)
        R[:, np.abs(fr) > half] = 0.0
        k = np.argmax(R, axis=1)
        r = R[np.arange(nf), k]
        on = 10 * np.log10(r + 1e-30) >= ON_DB
        if not on.any():
            continue
        f_on = fr[k[on]] + mid
        d_on = 10 * np.log10(r[on])
        order = np.argsort(f_on, kind="stable")
        f_on, d_on = f_on[order], d_on[order]
        cuts = [0] + [int(i) + 1 for i in np.flatnonzero(np.diff(f_on) > CLUSTER_HZ)] + [len(f_on)]
        for c0, c1 in zip(cuts[:-1], cuts[1:]):
            if c1 - c0 >= MIN_FRAMES:
                out.append((float(np.median(f_on[c0:c1])), c1 - c0, round(float(d_on[c0:c1].max()), 1)))
    out.sort(key=lambda c: -c[1])
    return out[:MAX_CAND]


def _mavg(a: np.ndarray, n: int) -> np.ndarray:
    c = np.concatenate([[0.0], np.cumsum(a)])
    return (c[n:] - c[:-n]) / n


def _runs(on: np.ndarray) -> list[tuple[int, int, bool]]:
    """(start, length, value) of the runs of a boolean array."""
    if len(on) == 0:
        return []
    e = np.flatnonzero(on[1:] != on[:-1]) + 1
    b = [0] + [int(i) for i in e] + [len(on)]
    return [(b[i], b[i + 1] - b[i], bool(on[b[i]])) for i in range(len(b) - 1)]


def pulses(y: np.ndarray, fs: float):
    """(pulse list [(start us, mark us, space us)], mark sample ranges, envelope dB) or None."""
    env = _mavg(np.abs(y).astype(np.float64), SMOOTH)
    noise, peak = float(np.percentile(env, 10)), float(np.percentile(env, 99.5))
    if noise <= 0 or 20 * np.log10(peak / noise) < ENV_DB:
        return None
    on = env > np.sqrt(noise * peak)
    short = int(round(MIN_RUN_US * 1e-6 * fs))
    for s0, n, v in _runs(on):                    # debounce: noise spikes / dropouts join their neighbours
        if n < short:
            on[s0:s0 + n] = not v
    us = 1e6 / fs
    runs = _runs(on)
    out, marks = [], []
    for i, (s, n, v) in enumerate(runs):
        if not v:
            continue
        sp = runs[i + 1][1] * us if i + 1 < len(runs) else float("inf")
        out.append((s * us, n * us, sp))
        marks.append((s + SMOOTH // 2, s + n + SMOOTH // 2))
    return out, marks, round(20 * float(np.log10(peak / noise)), 1)


def symbols(seg: np.ndarray, fs: float) -> list[int] | None:
    """PCM symbols of an FSK mark, or None when it is not FSK."""
    if len(seg) < 2 * SMOOTH:
        return None
    d = np.angle(seg[1:] * np.conj(seg[:-1])).astype(np.float64) * (fs / (2 * np.pi))
    d = _mavg(d, SMOOTH)
    d = d - 0.5 * (float(np.percentile(d, 10)) + float(np.percentile(d, 90)))   # mid-way between the tones
    a = np.abs(d)
    dev = float(np.median(a))                    # two-level: |d| sits near one deviation (a chirp spreads it)
    if dev < FSK_DEV_HZ or float(np.mean((a >= 0.6 * dev) & (a <= 1.4 * dev))) < FSK_FRAC:
        return None
    spb = SYMBOL_US * 1e-6 * fs
    out = []
    for _, n, v in _runs(d > 0):
        out += [int(v)] * int(round(n / spb))
    return out


def detect(ch) -> list[Emission]:
    out = []
    for f, frames, peak_db in candidates(ch):
        y, fs = ch.slice(f, SLICE_HALF, SLICE_FS)
        pr = pulses(y, fs)
        if pr is None:
            continue
        p, marks, env_db = pr
        base = {"carrier_mhz": round(f / 1e6, 4), "frames": frames, "scan_db": peak_db, "envelope_db": env_db,
                "n_pulses": len(p)}
        f_lo, f_hi = (f - SLICE_HALF) / 1e6, (f + SLICE_HALF) / 1e6
        found = False
        fsk_marks = 0
        for (t0, mk, _), (s0, s1) in zip(p, marks):
            if mk < FSK_MIN_US:
                continue
            sy = symbols(y[s0:s1], fs)
            if sy is None:
                continue
            fsk_marks += 1
            r = T.decode(sy)
            if r is not None:
                found = True
                out.append(Emission(family="ism_devices", label="tpms", start_us=t0, end_us=t0 + mk, f_lo_mhz=f_lo,
                                    f_hi_mhz=f_hi, center_mhz=f / 1e6, snr_db=env_db, confidence=0.9, decoded=r,
                                    features={**base, "symbols": len(sy), "mark_us": round(mk, 1)}))
        fc_ = O.fixed_code(p)
        if fc_ is not None:
            found = True
            out.append(Emission(family="ism_devices", label="ook_fixed_code", start_us=fc_["t_first_us"],
                                end_us=fc_["t_end_us"], f_lo_mhz=f_lo, f_hi_mhz=f_hi, center_mhz=f / 1e6,
                                snr_db=env_db, confidence=0.75 if fc_["repeats"] >= 2 else 0.45, decoded=fc_,
                                features=base))
        kq = O.keeloq(p)
        if kq is not None:
            found = True
            out.append(Emission(family="ism_devices", label="keeloq", start_us=kq["t_first_us"],
                                end_us=p[-1][0] + p[-1][1], f_lo_mhz=f_lo, f_hi_mhz=f_hi, center_mhz=f / 1e6,
                                snr_db=env_db, confidence=0.75 if kq["matching_packets"] >= 2 else 0.6, decoded=kq,
                                features=base))
        if found:
            continue
        mk = [q[1] for q in p]
        if fsk_marks:
            out.append(Emission(family="ism_devices", label="fsk_device", start_us=p[0][0],
                                end_us=p[-1][0] + p[-1][1], f_lo_mhz=f_lo, f_hi_mhz=f_hi, center_mhz=f / 1e6,
                                snr_db=env_db, confidence=0.55, features={**base, "fsk_bursts": fsk_marks}))
        elif len(p) >= 16 and all(50.0 <= m <= 2000.0 for m in mk):
            out.append(Emission(family="ism_devices", label="ook_device", start_us=p[0][0],
                                end_us=p[-1][0] + p[-1][1], f_lo_mhz=f_lo, f_hi_mhz=f_hi, center_mhz=f / 1e6,
                                snr_db=env_db, confidence=0.55,
                                features={**base, "mark_short_us": round(float(np.percentile(mk, 10)), 1),
                                          "mark_long_us": round(float(np.percentile(mk, 90)), 1)}))
    return out
