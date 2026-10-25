"""Detection stage: IQ -> list[Signal] with provenance."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..io.sigmf import Capture
from ..model import Signal
from . import artifacts, cfar, merge, noise
from .stft import Spectrogram, stft

COARSE_NFFT = 1024          # 60 kHz x 16.7 us at 61.44 MSps
FINE_NFFT = 16384           # 3.75 kHz x 267 us
FINE_HOP = 16384           # no overlap: CFAR statistics assume independent cells


@dataclass
class DetectResult:
    signals: list[Signal]
    coarse: Spectrogram
    floor_coarse_db: np.ndarray
    fine: Spectrogram
    floor_fine_db: np.ndarray
    usable_hz: tuple[float, float]
    info: dict = field(default_factory=dict)


def detect(cap: Capture, model: noise.NoiseModel | None = None, usable_frac: float = 0.95) -> DetectResult:
    x, fs, fc = cap.iq, cap.fs, cap.fc
    C = stft(x, fs, fc, COARSE_NFFT)
    prof = noise.profile_db(C.P)
    floor_c, ninfo = noise.estimate_floor(prof, fc, C.df, model)
    F = stft(x, fs, fc, FINE_NFFT, FINE_HOP)
    floor_f = 10 * np.log10(np.interp(F.f, C.f, 10 ** (floor_c / 10)))

    prof_f = noise.profile_db(F.P)
    # strong switched emitters lift the whole floor while on: per-frame gain g(t)
    g_c = noise.time_gain(C.P, floor_c, prof)
    g_f = np.interp(F.t + F.frame_len / 2, C.t + C.frame_len / 2, g_c).astype(np.float32)
    C.P /= g_c[:, None]
    F.P /= g_f[:, None]
    # always-on emitters first, then bursts above max(floor, always-on level)
    cc = cfar.detect_continuous(C, prof, floor_c, "coarse")
    cf = cfar.detect_continuous(F, prof_f, floor_f, "fine", min_bins=1, flank_db=3.0)
    bc = cfar.detect_bursts(C, floor_c, np.maximum(floor_c, prof), "coarse")
    bf = cfar.detect_bursts(F, floor_f, np.maximum(floor_f, prof_f), "fine",
                            passes=cfar.BURST_PASSES[:1], flank_db=3.0)
    boxes = merge.fuse_resolutions(cc, cf) + merge.fuse_resolutions(bc, bf)
    boxes = merge.merge_fragments(boxes)
    split = []
    for b in boxes:
        split += merge.split_valleys(b, C, floor_c) if b.res == "coarse" else [b]
    sigs = merge.absorb_into_continuous(merge.group_signals(split))

    half = usable_frac * fs / 2
    if cap.bandwidth_hz:                    # the analog filter, when narrower than the sampled band
        half = min(half, cap.bandwidth_hz / 2)
    lo, hi = fc - half, fc + half
    for s in sigs:
        if s.f0 < lo or s.f1 > hi:
            s.flags.append("truncated_at_capture_edge")
        res = sorted({b.res for b in s.bursts})
        s.detection = {
            "method": "STFT energy detection, hysteresis CFAR",
            "spectrogram": ("coarse 1024-pt Hann, 60 kHz x 16.7 us" if "coarse" in res else "") +
                           (" / " if len(res) == 2 else "") +
                           ("fine 16384-pt Hann, 3.75 kHz x 267 us" if "fine" in res else ""),
            "noise_floor": ninfo["method"] + "; per-frame gain g(t) from noise-only bins",
            "thresholds": (f"continuous: 20th-pct-over-time profile > {cfar.CONT_SEED_DB} dB (grow "
                           f"{cfar.CONT_GROW_DB} dB) above floor; bursts: hysteresis CFAR on "
                           f"P/max(floor, continuous level), passes (t x f cells, seed/grow dB) "
                           f"{[p for p in cfar.BURST_PASSES]}; fine-res bursts need 3 dB flank contrast"),
            "kinds": sorted({b.kind for b in s.bursts}),
            "n_bursts": len(s.bursts),
            "description": (f"{len(s.bursts)} time-frequency region(s) rose above the estimated noise "
                            f"floor (peak {max(b.peak_snr_db for b in s.bursts)} dB) and were "
                            f"connected/merged into one emitter by shared frequency span."),
        }
    artifacts.flag_lo_leakage(sigs, fc, F.df)
    artifacts.flag_iq_images(sigs, fc)
    sigs = artifacts.group_images(sigs)
    sigs.sort(key=lambda s: (s.f0, s.t0))
    ninfo["time_gain_db_p95"] = round(float(10 * np.log10(np.percentile(g_c, 95))), 2)
    info = {"noise": ninfo, "n_continuous": (len(cc), len(cf)), "n_bursts": (len(bc), len(bf)),
            "noise_floor_dbfs_per_bin_median": round(float(np.median(floor_c)), 1)}
    return DetectResult(sigs, C, floor_c, F, floor_f, (lo, hi), info)
