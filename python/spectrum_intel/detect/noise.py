"""Noise-floor estimation.

Per bin, the q-th percentile over time of the coarse spectrogram is the floor
wherever a signal is off at least (100-q)% of the time; for exponential
(complex-Gaussian) noise that percentile sits -ln(1-q) below the mean, which we
correct. Always-on wideband carriers (LTE/NR downlink, TV) never switch off, so
on its own the percentile would call them noise. Two ways round that:

  * sweep template (preferred): the receiver's noise *shape* across baseband
    (AD936x filters, LO-dependent edge roll-up) is shared by neighbouring sweep
    steps, while a real emitter sits at a different baseband offset in each of
    them. A low percentile across neighbours of their level-normalised profiles
    gives the shape with the emitters removed; one offset per snapshot fits the
    level.
  * morphological opening (fallback, single file): removes anything narrower
    than the window from the profile; misses always-on signals wider than it.
    Those are bridged: where the opening sits >= 10 dB above the band's quiet
    level (10th percentile of the opening over the central 2/3 of the band) it is
    standing on a wide always-on carrier (live: a 20 MHz LTE B3 carrier 40 dB up),
    so across each such run (grown down to 3 dB above the quiet level) the floor
    is interpolated from the quiet bins on either side.

The final floor is the lower of the two (each over-estimates in a different
situation).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.ndimage import grey_opening, median_filter, uniform_filter1d

Q = 20.0                                    # percentile over time
Q_BIAS_DB = -10 * np.log10(-np.log(1 - Q / 100))   # +6.52 dB: percentile -> mean for exp. noise
NFFT = 1024


def profile_db(P: np.ndarray) -> np.ndarray:
    """Per-bin noise-mean estimate (dB) from a (T, F) power spectrogram."""
    return 10 * np.log10(np.percentile(P, Q, axis=0) + 1e-30) + Q_BIAS_DB


def _patch_dc(v: np.ndarray, half: int = 2) -> np.ndarray:
    """LO leakage sits on the DC bin of every snapshot; keep it out of the floor."""
    v = v.copy()
    c = len(v) // 2
    lo, hi = c - half - 1, c + half + 1
    v[lo + 1:hi] = np.linspace(v[lo], v[hi], hi - lo + 1)[1:-1]
    return v


def opening_floor(prof_db: np.ndarray, df: float, win_hz: float = 3e6) -> np.ndarray:
    k = max(3, int(win_hz / df) | 1)
    fl = grey_opening(_patch_dc(prof_db), size=k, mode="nearest")
    return uniform_filter1d(fl, k, mode="nearest")


PLATEAU_SEED_DB = 10.0                      # opening this far above the quiet level: a wide always-on carrier
PLATEAU_GROW_DB = 3.0                       # ... and its run extends while above this


def bridge_plateaus(op: np.ndarray) -> tuple[np.ndarray, int]:
    """Opening floor with wide always-on carriers bridged (see module doc); returns (floor, bins bridged)."""
    n = len(op)
    quiet = float(np.percentile(op[n // 6:n - n // 6], 10))
    seed, grow = op - quiet > PLATEAU_SEED_DB, op - quiet > PLATEAU_GROW_DB
    if not seed.any():
        return op, 0
    out = op.copy()
    edges = np.flatnonzero(np.diff(np.r_[0, grow.astype(np.int8), 0]))
    nb = 0
    for i0, i1 in zip(edges[::2], edges[1::2]):     # runs [i0, i1) of grow
        if not seed[i0:i1].any():
            continue
        if i0 == 0 and i1 == n:
            return op, 0                            # no quiet bin to anchor on
        if i0 == 0:
            out[i0:i1] = op[i1]
        elif i1 == n:
            out[i0:i1] = op[i0 - 1]
        else:
            out[i0:i1] = np.interp(np.arange(i0, i1), [i0 - 1, i1], [op[i0 - 1], op[i1]])
        nb += i1 - i0
    return out, nb


def _mode(v: np.ndarray, step: float = 0.25) -> float:
    v = v[np.isfinite(v)]
    h, e = np.histogram(v, bins=np.arange(v.min(), v.max() + step, step))
    h = uniform_filter1d(h.astype(float), 5)
    return float(e[np.argmax(h)] + step / 2)


@dataclass
class NoiseModel:
    """Sweep-wide store of per-snapshot noise profiles (coarse, NFFT bins)."""
    fcs: np.ndarray                     # (N,) centre frequencies
    profiles: np.ndarray                # (N, NFFT) profile_db of each snapshot
    names: list[str] = field(default_factory=list)
    neighbour_hz: float = 100e6

    def save(self, path: Path):
        np.savez_compressed(path, fcs=self.fcs, profiles=self.profiles, names=np.array(self.names))

    @classmethod
    def load(cls, path: Path) -> "NoiseModel":
        d = np.load(path)
        return cls(d["fcs"], d["profiles"], list(d["names"]))

    def template(self, fc: float):
        # other tuning steps only: repeats at the same fc share the target's emitters
        d = np.abs(self.fcs - fc)
        idx = np.where((d <= self.neighbour_hz) & (d > 1e3))[0]
        if len(idx) < 6 or len(set(np.round(self.fcs[idx]))) < 3:
            return None, 0
        P = self.profiles[idx]
        F = P.shape[1]
        central = slice(F // 6, F - F // 6)
        rel = P - np.percentile(P[:, central], 10, axis=1, keepdims=True)
        tmpl = np.percentile(rel, 20, axis=0)
        tmpl = _patch_dc(median_filter(tmpl, 9, mode="nearest"))
        return tmpl, len(idx)


def estimate_floor(prof_db: np.ndarray, fc: float, df: float, model: NoiseModel | None = None) -> tuple[np.ndarray, dict]:
    """Noise floor (dB per bin, same scale as profile_db) + provenance dict."""
    op = opening_floor(prof_db, df)
    info = {"percentile_over_time": Q, "percentile_bias_db": round(Q_BIAS_DB, 2),
            "opening_window_hz": 3e6}
    tmpl, n = (None, 0) if model is None else model.template(fc)
    if tmpl is None:
        op, nb = bridge_plateaus(op)
        info["method"] = ("per-snapshot: per-bin 20th percentile over time, morphological opening (3 MHz), "
                          "wide always-on carriers (>= 10 dB above the quiet level) bridged")
        info["bridged_hz"] = round(nb * df)
        return op, info
    r = _patch_dc(prof_db) - tmpl
    off = _mode(r[r <= np.median(r)])
    tf = tmpl + off
    floor = np.minimum(tf, op)
    info.update(method="sweep template (neighbour snapshots within ±100 MHz, 20th pct of level-"
                       "normalised profiles) + per-snapshot offset (histogram mode); "
                       "min() with per-snapshot opening floor",
                n_neighbour_snapshots=n, template_offset_db=round(off, 2),
                frac_bins_template_below_opening=round(float(np.mean(tf < op)), 3))
    return floor, info


def time_gain(P: np.ndarray, floor_db: np.ndarray, prof_db: np.ndarray, q: float = 30.0,
              smooth: int = 5, max_excess_db: float = 1.0, min_bins: int = 32) -> np.ndarray:
    """Per-frame noise-level multiplier g(t) >= 1.

    A strong (often out-of-band) emitter that switches on and off compresses
    the front end / ADC and lifts the whole floor with it, so the percentile
    floor (which sees the quiet frames) under-reads during its on-time. g(t)
    is the q-th percentile, per frame, of P/floor over noise-only bins
    (profile within max_excess_db of the floor), bias-corrected for
    exponential noise and smoothed over `smooth` frames.
    """
    noise_bins = np.where(prof_db - floor_db < max_excess_db)[0]
    if len(noise_bins) < min_bins:
        return np.ones(P.shape[0], np.float32)
    r = P[:, noise_bins] / (10 ** (floor_db[noise_bins] / 10))[None, :]
    g = np.percentile(r, q, axis=1) / (-np.log(1 - q / 100))
    g = uniform_filter1d(g, smooth, mode="nearest")
    return np.maximum(g, 1.0).astype(np.float32)
