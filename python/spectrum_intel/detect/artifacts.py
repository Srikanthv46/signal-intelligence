"""Receiver artifacts that look like signals: LO leakage at DC, IQ-imbalance
images, ADC clipping."""
from __future__ import annotations

import numpy as np

from ..core import spurs
from ..model import Signal

ADC_FULL_SCALE = 2048.0          # AD936x 12-bit, as written by signal-vault (raw counts)


def clipping(x: np.ndarray, full_scale: float = ADC_FULL_SCALE) -> dict:
    a = np.maximum(np.abs(x.real), np.abs(x.imag))
    peak = float(a.max())
    # raw counts if the peak is in the ADC range; else assume normalised [-1, 1]
    fs_ = full_scale if peak > 4 else 1.0
    # raw counts reach full_scale - 1 (2047); normalised samples reach ~1.0
    frac = float(np.mean(a >= 0.995 * (fs_ - 1 if fs_ > 1 else fs_)))
    rms = float(np.sqrt(np.mean(np.abs(x) ** 2)))
    return {"peak_abs": round(peak, 1), "full_scale": fs_, "clipped_fraction": frac,
            "clipping": frac > 1e-6,
            "rms_dbfs": round(20 * np.log10(rms / fs_ + 1e-12), 1)}


def flag_lo_leakage(sigs: list[Signal], fc: float, df_fine: float):
    for s in sigs:
        if s.bw <= 6 * df_fine + 1 and abs(s.fc - fc) <= 2 * df_fine:
            s.flags.append("lo_leakage")


def flag_receiver_spurs(sigs: list[Signal], lines: list[dict], fs: float, ppm: float | None,
                        max_bw: float = 25e3):
    """Mark detections of the receiver's own lines (core.spurs.receiver_lines): a narrow detection on a line is
    flagged receiver_spur; a wider one that contains a line only records it."""
    for s in sigs:
        if s.flags:
            continue
        on = [ln for ln in lines if s.f0 - 1e3 <= ln["freq_hz"] <= s.f1 + 1e3]
        if not on:
            continue
        ln = max(on, key=lambda x: x["peak_db"])
        if s.bw <= max_bw:
            s.flags.append("receiver_spur")
            s.features.update(spur_kind=ln["kind"], spur_freq_hz=ln["freq_hz"], spur_sample_rate_hz=fs)
            if ln["harmonic"]:
                s.features["spur_harmonic"] = ln["harmonic"]
            if ln["kind"] == "reference_locked":
                s.features["spur_receiver_ppm"] = round(float(ppm), 2)
        else:
            s.features["contains_receiver_spur_hz"] = [x["freq_hz"] for x in on]


def flag_iq_images(sigs: list[Signal], fc: float, min_rej_db: float = 15.0):
    """An image is the mirror (2fc - f) of a stronger signal with the same burst
    timing, at least min_rej_db weaker (AD936x image rejection is ~40-50 dB, so
    only very strong emitters leave visible images)."""
    for s in sigs:
        if "lo_leakage" in s.flags:
            continue
        m0, m1 = 2 * fc - s.f1, 2 * fc - s.f0
        for o in sigs:
            if o is s or o.snr_db - s.snr_db < min_rej_db:
                continue
            fov = max(0.0, min(m1, o.f1) - max(m0, o.f0))
            if fov < 0.6 * max(s.bw, 1.0):
                continue
            tov = sum(max(0.0, min(a.t1, b.t1) - max(a.t0, b.t0)) for a in s.bursts for b in o.bursts)
            if tov >= 0.7 * sum(b.t1 - b.t0 for b in s.bursts):
                s.flags.append("iq_image")
                s.features["image_of_f_center_hz"] = round(o.fc)
                break


def group_images(sigs: list[Signal]) -> list[Signal]:
    """One artifact entry per mirrored source instead of one per fragment."""
    keep, groups = [], {}
    for s in sigs:
        if "iq_image" in s.flags:
            groups.setdefault(s.features.get("image_of_f_center_hz"), []).append(s)
        else:
            keep.append(s)
    for src, g in groups.items():
        if len(g) == 1:
            keep.append(g[0])
            continue
        m = Signal(bursts=[b for x in g for b in x.bursts], f0=min(x.f0 for x in g), f1=max(x.f1 for x in g),
                   snr_db=max(x.snr_db for x in g))
        m.flags = ["iq_image"]
        m.features["image_of_f_center_hz"] = src
        m.features["n_fragments"] = len(g)
        m.detection = dict(g[0].detection, n_bursts=len(m.bursts),
                           description=f"{len(g)} energy fragments mirrored about the LO from the signal at "
                                       f"{src / 1e6:.3f} MHz, merged into one artifact entry")
        keep.append(m)
    return keep
