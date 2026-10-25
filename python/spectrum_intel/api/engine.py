"""Production entry point: one IQ block in, one response dict out.

    from spectrum_intel.api import Engine, IQBlock
    engine = Engine()                                   # load receivers once, reuse for every block
    resp = engine.analyze(IQBlock(iq, sample_rate_hz=61.44e6, center_frequency_hz=2.437e9,
                                  bandwidth_hz=56e6, start_time="2026-10-09T10:00:00Z"))

The response follows schema/response.schema.json. The engine never raises for bad
input or a failing receiver: it returns status "insufficient_input" / "error" with the
reason in errors[], or "ok" with receiver failures listed in errors[].
"""
from __future__ import annotations

import time
import traceback
from dataclasses import dataclass, field
from datetime import datetime, timezone

import numpy as np

from .. import __version__
from ..analyze import analyze
from ..core import parallel
from ..decode import registry
from ..io.sigmf import Capture
from ..report import writer
from . import contract as K
from . import response as R

IQ_FORMATS = ("cf32", "cf64", "ci16", "ci8", "cu8")


def to_complex64(iq, iq_format: str | None = None) -> np.ndarray:
    """Complex64 samples from what an SDR pipeline delivers: a complex array, or
    interleaved I/Q (1-D real array, or shape (N, 2)) as int16/int8/uint8/float.
    Values are not rescaled: raw ADC counts keep the clipping check meaningful."""
    a = np.asarray(iq)
    if np.iscomplexobj(a):
        return np.ascontiguousarray(a.ravel(), np.complex64)
    if iq_format == "cu8" or a.dtype == np.uint8:
        a = a.astype(np.float32) - 127.5
    a = a.astype(np.float32, copy=False)
    if a.ndim == 2 and a.shape[1] == 2:
        a = a.ravel()
    if a.ndim != 1 or len(a) % 2:
        raise ValueError(f"real IQ must be interleaved I/Q of even length or shape (N, 2); got shape {np.shape(iq)}")
    return (a[0::2] + 1j * a[1::2]).astype(np.complex64)


def _parse_time(t) -> datetime | None:
    if t is None or isinstance(t, datetime):
        return t
    if isinstance(t, (int, float)):
        return datetime.fromtimestamp(float(t), tz=timezone.utc)
    return datetime.fromisoformat(str(t).replace("Z", "+00:00"))


@dataclass
class IQBlock:
    """One block from the acquisition pipeline.

    iq                   complex samples, or interleaved I/Q (see to_complex64)
    sample_rate_hz       complex sample rate
    center_frequency_hz  RF frequency of baseband 0 Hz (the LO); trusted as given
    bandwidth_hz         analog front-end bandwidth (None: the whole sampled band is usable)
    start_time           UTC time of the first sample: datetime, ISO 8601 string or Unix seconds
    receiver_ppm         the SDR's frequency error in ppm (measured LO/TCXO offset), if known: LTE, NR,
                         GSM and voice receivers search around it; unknown, they estimate it per block
                         where they can, and weak cells may be missed
    gain_db, source_id, request_id  passed through to the response
    """
    iq: object
    sample_rate_hz: float
    center_frequency_hz: float
    bandwidth_hz: float | None = None
    start_time: object = None
    gain_db: float | None = None
    source_id: str | None = None
    request_id: str | None = None
    iq_format: str | None = None
    receiver_ppm: float | None = None
    _samples: np.ndarray | None = field(default=None, repr=False)

    def __post_init__(self):
        self.start_time = _parse_time(self.start_time)

    @property
    def samples(self) -> np.ndarray:
        if self._samples is None:
            self._samples = to_complex64(self.iq, self.iq_format)
        return self._samples

    @property
    def n_samples(self) -> int:
        try:
            return len(self.samples)
        except Exception:                                   # noqa: BLE001 — malformed input reported by validate()
            return 0


def validate(block: IQBlock) -> str | None:
    """Why this block cannot be analysed, or None."""
    for name in ("sample_rate_hz", "center_frequency_hz"):
        v = getattr(block, name)
        if not isinstance(v, (int, float, np.integer, np.floating)) or not np.isfinite(v):
            return f"{name} must be a finite number, got {v!r}"
    if block.sample_rate_hz <= 0:
        return f"sample_rate_hz must be > 0, got {block.sample_rate_hz}"
    if block.bandwidth_hz is not None and not (np.isfinite(block.bandwidth_hz) and block.bandwidth_hz > 0):
        return f"bandwidth_hz must be > 0 or None, got {block.bandwidth_hz}"
    if block.receiver_ppm is not None and not (np.isfinite(block.receiver_ppm) and abs(block.receiver_ppm) < 1000):
        return f"receiver_ppm must be a finite number of ppm (|x| < 1000), got {block.receiver_ppm}"
    try:
        x = block.samples
    except Exception as e:                                  # noqa: BLE001
        return f"iq: {e}"
    need = K.requirements()["pipeline"]["min_samples"]
    if len(x) < need:
        return (f"{len(x)} samples is below the minimum of {need} for analysis "
                f"(recommended {K.requirements()['pipeline']['recommended_samples']}); "
                f"accumulate blocks with spectrum_intel.api.StreamAnalyzer")
    if not np.all(np.isfinite(x[:: max(1, len(x) // 4096)])) or not np.isfinite(np.sum(x)):
        return "iq contains NaN or Inf"
    return None


class Engine:
    """Holds the receivers (loaded once) and analyses blocks. Thread-compatible: one
    engine may serve one block at a time per thread; receivers keep no per-block state."""

    def __init__(self, enable: list[str] | None = None, disable: list[str] | None = None,
                 detail: str = "full"):
        if detail not in ("full", "compact"):
            raise ValueError("detail must be 'full' or 'compact'")
        self.plugins = registry.plugins(enable, disable)
        self.detail = detail
        self.packages: dict[str, list[str]] = {}
        for p in self.plugins:
            m = type(p).__module__.split(".")
            pkg = m[m.index("receivers") + 1] if "receivers" in m else p.name
            self.packages.setdefault(pkg, []).append(p.name)

    @property
    def info(self) -> dict:
        return {"name": "spectrum_intel", "implementation": "python", "version": __version__,
                "receivers_enabled": list(self.packages), "receivers_unavailable": registry.unavailable(),
                "threads": parallel.threads()}

    def analyze(self, block: IQBlock) -> dict:
        t = time.time()
        why = validate(block)
        if why:
            return R.failure(block, self.info, "insufficient_input", "INVALID_INPUT", why, time.time() - t)
        try:
            cap = Capture.from_iq(block.samples, float(block.sample_rate_hz), float(block.center_frequency_hz),
                                  t0=block.start_time, name=block.source_id or "block", gain_db=block.gain_db,
                                  bandwidth_hz=block.bandwidth_hz)
            ppm = None if block.receiver_ppm is None else {
                "ppm": float(block.receiver_ppm), "source": "caller (IQBlock.receiver_ppm)"}
            a = analyze(cap, None, self.plugins, ppm=ppm)
            for s in a.signals:                              # binary products are for files, not the response
                s.features.pop("_video_frame", None)
                s.features.pop("_audio", None)
            rep = writer.build(a)
            return R.build(rep, analysis=a, block=block, engine=self.info, packages=self.packages,
                           total_s=time.time() - t, detail=self.detail)
        except Exception as e:                              # noqa: BLE001 — the caller gets a response, always
            return R.failure(block, self.info, "error", "ANALYSIS_FAILED",
                             f"{type(e).__name__}: {e} | {traceback.format_exc()[-600:]}", time.time() - t)

    def analyze_iq(self, iq, sample_rate_hz: float, center_frequency_hz: float, bandwidth_hz: float | None = None,
                   **kw) -> dict:
        return self.analyze(IQBlock(iq, sample_rate_hz, center_frequency_hz, bandwidth_hz, **kw))


_DEFAULT: Engine | None = None


def analyze_block(iq, sample_rate_hz: float, center_frequency_hz: float, bandwidth_hz: float | None = None,
                  **kw) -> dict:
    """One-call form of Engine().analyze_iq(...) with a process-wide default engine."""
    global _DEFAULT
    if _DEFAULT is None:
        _DEFAULT = Engine()
    return _DEFAULT.analyze_iq(iq, sample_rate_hz, center_frequency_hz, bandwidth_hz, **kw)
