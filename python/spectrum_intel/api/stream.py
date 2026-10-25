"""Accumulate small SDR blocks (2048, 4096, ... samples) into analysable ones.

    stream = StreamAnalyzer(engine, target_samples=2**21)
    for buf in sdr:                                      # contiguous blocks of one tuning
        for resp in stream.push(buf, fs, fc, bw, start_time=t):
            publish(resp)
    resp = stream.flush()                                # what is left (None if too short)

A change of centre frequency, sample rate or bandwidth, or a gap in time (when
start times are given), closes the current block: samples of different tunings or
discontinuous samples are never analysed together.
"""
from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np

from . import contract as K
from .engine import Engine, IQBlock, _parse_time, to_complex64


class StreamAnalyzer:
    def __init__(self, engine: Engine | None = None, target_samples: int | None = None,
                 max_gap_samples: int = 1):
        self.engine = engine or Engine()
        req = K.requirements()["pipeline"]
        self.target = int(target_samples or req["recommended_samples"])
        if self.target < req["min_samples"]:
            raise ValueError(f"target_samples must be >= {req['min_samples']}")
        self.max_gap = max_gap_samples
        self._key = None                    # (fs, fc, bw) of the samples being collected
        self._t0: datetime | None = None    # time of the first buffered sample (or of the next one, when empty)
        self._buf: list[np.ndarray] = []
        self._n = 0
        self._meta: dict = {}

    @property
    def buffered(self) -> int:
        return self._n

    def push(self, iq, sample_rate_hz: float, center_frequency_hz: float, bandwidth_hz: float | None = None,
             start_time=None, iq_format: str | None = None, **meta) -> list[dict]:
        """Add one block; returns the responses of every block it completed (usually 0 or 1)."""
        x = to_complex64(iq, iq_format)
        t = _parse_time(start_time)
        key = (float(sample_rate_hz), float(center_frequency_hz), bandwidth_hz)
        out = []
        if self._key is not None and (key != self._key or self._gap(t)):
            r = self.flush()
            if r is not None:
                out.append(r)
        if self._key is None or self._n == 0:
            continues = self._key == key            # same tuning: the next sample time is already known
            self._key, self._meta = key, meta
            self._t0 = t if t is not None else (self._t0 if continues else None)
        pos = 0
        while pos < len(x):
            take = min(len(x) - pos, self.target - self._n)
            self._buf.append(x[pos:pos + take])
            self._n += take
            pos += take
            if self._n >= self.target:
                out.append(self._analyse(keep_going=True))
        return out

    def flush(self) -> dict | None:
        """Analyse what is buffered if it is long enough (else drop it) and forget the tuning."""
        r = self._analyse(keep_going=False) if self._n >= K.requirements()["pipeline"]["min_samples"] else None
        self._key, self._t0, self._buf, self._n, self._meta = None, None, [], 0, {}
        return r

    def _gap(self, t: datetime | None) -> bool:
        if t is None or self._t0 is None or self._n == 0:
            return False
        expected = self._t0 + timedelta(seconds=self._n / self._key[0])
        return abs((t - expected).total_seconds()) * self._key[0] > self.max_gap

    def _analyse(self, keep_going: bool) -> dict:
        fs, fc, bw = self._key
        x = np.concatenate(self._buf) if len(self._buf) > 1 else self._buf[0]
        resp = self.engine.analyze(IQBlock(x, fs, fc, bw, start_time=self._t0, **self._meta))
        if keep_going:                      # the next sample continues the same tuning
            self._t0 = None if self._t0 is None else self._t0 + timedelta(seconds=self._n / fs)
            self._buf, self._n = [], 0
        return resp
