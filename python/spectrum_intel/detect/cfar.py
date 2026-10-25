"""Energy detection on the noise-normalised spectrogram -> time-frequency boxes.

Two complementary detectors:

continuous  per-bin 20th-percentile-over-time profile vs. the noise floor.
            Over ~2000 frames that percentile is stable to ~0.1 dB, so a
            carrier >= 3 dB above the floor in 80% of frames is certain. Each
            contiguous run of such bins is one always-on emitter (t = whole capture).

bursts      hysteresis CFAR on P / baseline, baseline = max(floor, continuous
            level) per bin: a burst must rise above whatever is always there.
            Averaging k cells of exponential noise gives Gamma(k, 1/k):
              3x3 cells, seed 6 dB / grow 3 dB  -> ~1e-7 false seeds per cell
              15x5 cells, seed 3 dB / grow 1.5 dB (weak wideband bursts)
            Both masks are OR-ed; components need at least one seed.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from ..model import Box
from .stft import Spectrogram

# (time cells, freq cells, seed dB, grow dB)
BURST_PASSES = ((3, 3, 6.0, 3.0), (15, 5, 3.0, 1.5))
CONT_SEED_DB, CONT_GROW_DB = 3.0, 1.5
SEED_DB, GROW_DB = 6.0, 3.0


def _runs(mask: np.ndarray):
    d = np.diff(np.concatenate([[0], mask.astype(np.int8), [0]]))
    return list(zip(np.where(d == 1)[0], np.where(d == -1)[0]))


def detect_continuous(S: Spectrogram, prof_db: np.ndarray, floor_db: np.ndarray, res: str,
                      min_bins: int = 2, flank_db: float | None = None) -> list[Box]:
    excess = ndimage.uniform_filter1d(prof_db - floor_db, 3)
    seed, grow = excess > CONT_SEED_DB, excess > CONT_GROW_DB
    T = S.t[-1] + S.frame_len
    out = []
    for a, b in _runs(grow):
        if b - a < min_bins or not seed[a:b].any():
            continue
        if flank_db is not None:
            # narrow tone test: the run must stand above the profile beside it
            w = max(b - a, 4)
            fl = [excess[max(0, a - 1 - w):a - 1].mean() if a - 1 - w >= 0 else -np.inf,
                  excess[b + 1:b + 1 + w].mean() if b + 1 + w <= len(excess) else -np.inf]
            if excess[a:b].mean() - max(fl) < flank_db:
                continue
        lin = 10 ** (excess[a:b] / 10)
        out.append(Box(0.0, float(T), float(S.f[a] - S.df / 2), float(S.f[b - 1] + S.df / 2),
                       round(float(excess[a:b].max()), 1),
                       round(float(10 * np.log10(max(np.mean(lin) - 1, 1e-3))), 1), res, int(b - a),
                       kind="continuous"))
    return out


def detect_bursts(S: Spectrogram, floor_db: np.ndarray, baseline_db: np.ndarray, res: str,
                  passes=BURST_PASSES, min_cells: int = 6, strong_db: float = 15.0,
                  flank_db: float | None = None) -> list[Box]:
    rel = S.P / (10 ** (baseline_db / 10))[None, :]
    seed = np.zeros(rel.shape, bool)
    grow = np.zeros(rel.shape, bool)
    sm0 = None
    for nt, nf, sd, gr in passes:
        sm = ndimage.uniform_filter(rel, size=(nt, nf), mode="nearest")
        sm0 = sm if sm0 is None else sm0
        seed |= sm > 10 ** (sd / 10)
        grow |= sm > 10 ** (gr / 10)
    lab, n = ndimage.label(grow, structure=np.ones((3, 3)))
    if n == 0:
        return []
    has_seed = np.zeros(n + 1, bool)
    has_seed[np.unique(lab[seed])] = True
    has_seed[0] = False
    snr_floor = ndimage.uniform_filter(S.P / (10 ** (floor_db / 10))[None, :], size=(3, 3), mode="nearest")
    boxes = []
    for k, sl in enumerate(ndimage.find_objects(lab), start=1):
        if not has_seed[k]:
            continue
        m = lab[sl] == k
        cells = int(m.sum())
        pk_rel = float(10 * np.log10(sm0[sl][m].max()))
        if cells < min_cells and pk_rel < strong_db:
            continue
        ts, fsl = sl
        if flank_db is not None and not _flank_ok(sm0, ts, fsl, flank_db):
            continue
        lin = snr_floor[sl][m]
        boxes.append(Box(
            t0=float(S.t[ts.start]), t1=float(S.t[ts.stop - 1] + S.frame_len),
            f0=float(S.f[fsl.start] - S.df / 2), f1=float(S.f[fsl.stop - 1] + S.df / 2),
            peak_snr_db=round(float(10 * np.log10(lin.max())), 1),
            mean_snr_db=round(float(10 * np.log10(max(np.mean(lin) - 1, 1e-3))), 1),
            res=res, n_cells=cells))
    return boxes


def _flank_ok(sm, ts, fsl, min_db):
    """Narrowband test: the box must stand min_db above the spectrum on both
    sides over the same time span (rejects fragments of a weak wideband signal)."""
    w = max(fsl.stop - fsl.start, 4)
    F = sm.shape[1]
    inner = sm[ts, fsl].mean()
    flanks = []
    for a, b in ((fsl.start - 1 - w, fsl.start - 1), (fsl.stop + 1, fsl.stop + 1 + w)):
        a, b = max(0, a), min(F, b)
        if b - a >= 2:
            flanks.append(sm[ts, a:b].mean())
    return not flanks or inner >= 10 ** (min_db / 10) * max(flanks)
