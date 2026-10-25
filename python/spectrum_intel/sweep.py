"""Sweep-level context shared by all snapshots of a folder.

1. Noise profiles of every snapshot (one cheap pass) -> noise-floor templates
   from neighbouring tunings.
2. Receiver frequency-error calibration (ppm) from the strongest LTE PSS.

The capture metadata is trusted: each snapshot's centre frequency is the LO it
was captured at (acquisition must flush the receiver buffers after retuning).
"""
from __future__ import annotations

import json
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .detect import noise
from .detect.stft import stft
from .io.sigmf import find_captures, load

CACHE = "sweep_context.npz"


def _profile(p):
    c = load(p)
    S = stft(c.iq, c.fs, c.fc, noise.NFFT)
    return c.fc, c.t0.timestamp() if c.t0 else 0.0, noise.profile_db(S.P), c.name


@dataclass
class SweepContext:
    names: list[str]
    fc: np.ndarray                                    # centre frequency from each snapshot's metadata
    t: np.ndarray
    profiles: np.ndarray
    ppm: dict | None = None                           # receiver frequency error calibration

    @property
    def noise_model(self) -> noise.NoiseModel:
        return noise.NoiseModel(self.fc, self.profiles, self.names)

    def save(self, path: Path):
        np.savez_compressed(path, names=np.array(self.names), fc=self.fc, t=self.t, profiles=self.profiles,
                            ppm=json.dumps(self.ppm))

    @classmethod
    def load(cls, path: Path) -> "SweepContext":
        d = np.load(path)
        fc = d["fc"] if "fc" in d else d["fc_label"]  # caches written before 2026-10-05 store the label
        return cls(list(d["names"]), fc, d["t"], d["profiles"],
                   json.loads(str(d["ppm"])) if "ppm" in d else None)


def build_context(folder: str | Path, cache_dir: Path | None = None, workers: int = 8) -> SweepContext:
    cache = (Path(cache_dir) / CACHE) if cache_dir else None
    if cache and cache.exists():
        return SweepContext.load(cache)
    paths = find_captures(folder)
    with ProcessPoolExecutor(workers) as ex:
        r = list(ex.map(_profile, paths, chunksize=4))
    fc = np.array([a[0] for a in r])
    t = np.array([a[1] for a in r])
    prof = np.array([a[2] for a in r])
    names = [a[3] for a in r]
    ctx = SweepContext(names, fc, t, prof)
    ctx.ppm = calibrate_ppm(ctx, paths)
    if cache:
        cache.parent.mkdir(parents=True, exist_ok=True)
        ctx.save(cache)
    return ctx


def calibrate_ppm(ctx: SweepContext, paths: list[Path]) -> dict | None:
    """Receiver frequency error from the strongest LTE PSS in the sweep.

    One TCXO serves every tuning, so the error in ppm is common to the sweep.
    The tuning whose LTE-downlink bins rise most above the floor is searched on a
    1 MHz grid with a +-500 kHz offset search; the PSS carrier's distance from
    the 100 kHz EARFCN raster gives the error."""
    from .core.bands_india import band_keys
    from .core.dsp import Channelizer
    try:
        from .receivers.lte import cellsearch as lte
    except ImportError:                                     # LTE receiver removed: no sweep ppm
        return None
    by_name = {p.name[:-len(".sigmf-meta")]: p for p in paths}
    best = None
    for i, n in enumerate(ctx.names):
        fc = ctx.fc[i]
        f = fc + (np.arange(ctx.profiles.shape[1]) - ctx.profiles.shape[1] // 2) * 61.44e6 / ctx.profiles.shape[1]
        sel = np.array([bool(set(band_keys(v, v)) & set(lte.LTE_BANDS)) for v in f])
        if sel.sum() < 100:
            continue
        score = float(np.percentile(ctx.profiles[i][sel], 60) - np.percentile(ctx.profiles[i], 5))
        if best is None or score > best[0]:
            best = (score, n, fc, f[sel])
    if best is None or best[1] not in by_name:
        return None
    cap = load(by_name[best[1]])
    ch = Channelizer(cap.iq, cap.fs, cap.fc)
    fs_ = best[3]
    centres = np.arange(np.ceil(fs_.min() / 1e6) * 1e6, fs_.max(), 1e6)
    cal = lte.wide_ppm_search(ch, centres)
    if cal:
        cal.update(source=f"sweep: strongest LTE PSS in {best[1]} (tuned {best[2] / 1e6:.2f} MHz), "
                          f"+-500 kHz search on a 1 MHz grid", snapshot=best[1])
    return cal
