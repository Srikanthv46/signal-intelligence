"""Shared DSP helpers."""
from __future__ import annotations

import numpy as np
import scipy.fft as sfft


class Channelizer:
    """One FFT of the whole capture; any (centre, width, rate) slice is then a
    small IFFT (brick-wall filter — fine for detection and feature extraction)."""

    def __init__(self, x: np.ndarray, fs: float, fc: float):
        self.fs, self.fc, self.L = fs, fc, len(x)
        self.X = sfft.fft(x.astype(np.complex64), workers=-1)
        self.f = sfft.fftfreq(self.L, 1 / fs)

    def _bins(self, f_off: float, half_bw: float) -> np.ndarray:
        """Indices (ascending) of the bins with |f - f_off| <= half_bw: the same test as a scan of all L bins,
        evaluated only on the index range around the window (2 bins of margin)."""
        df = self.fs / self.L
        lo, hi = int(np.floor((f_off - half_bw) / df)) - 2, int(np.ceil((f_off + half_bw) / df)) + 2
        lo, hi = max(lo, -((self.L) // 2) - 1), min(hi, (self.L - 1) // 2 + 1)
        k = np.arange(lo, hi + 1)
        idx = np.unique(k % self.L)
        return idx[np.abs(self.f[idx] - f_off) <= half_bw]

    def slice(self, f_center_abs: float, half_bw: float, out_fs: float) -> tuple[np.ndarray, float]:
        """Complex baseband of [f_center ± half_bw] at ~out_fs. Returns (y, actual_fs)."""
        f_off = f_center_abs - self.fc
        M = max(16, int(round(self.L * out_fs / self.fs)))
        sel = self._bins(f_off, half_bw)
        k0 = int(round(f_off / self.fs * self.L))
        k = np.where(sel > self.L // 2, sel - self.L, sel)
        Y = np.zeros(M, np.complex64)
        Y[(k - k0) % M] = self.X[sel]
        return sfft.ifft(Y, workers=-1) * (M / self.L), M * self.fs / self.L


def burst_mask(n: int, fs: float, bursts, pad_s: float = 0.0) -> np.ndarray:
    m = np.zeros(n, bool)
    for b in bursts:
        m[max(0, int((b.t0 - pad_s) * fs)):min(n, int((b.t1 + pad_s) * fs))] = True
    return m


class IQContext:
    """Minimal stand-in for the pipeline's DecodeContext, for using a receiver outside
    the pipeline: ``dsss.scan(IQContext(iq, fs, fc))``. Provides .cap (iq, fs, fc, n,
    duration), .ch (Channelizer), .cache and .det (None, or an object with .signals)."""

    def __init__(self, iq: np.ndarray, fs: float, fc: float, det=None):
        from types import SimpleNamespace
        iq = np.asarray(iq, np.complex64)
        self.cap = SimpleNamespace(iq=iq, fs=float(fs), fc=float(fc), n=len(iq), duration=len(iq) / fs)
        self.ch = Channelizer(iq, fs, fc)
        self.det = det
        self.cache: dict = {}
