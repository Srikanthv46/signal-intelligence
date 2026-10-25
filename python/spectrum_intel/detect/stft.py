"""Short-time power spectra on absolute frequency / time axes."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.fft as sfft
from numpy.lib.stride_tricks import sliding_window_view


@dataclass
class Spectrogram:
    P: np.ndarray       # (T, F) float32 power per bin; white noise of variance s2 -> E[P] = s2
    f: np.ndarray       # (F,) absolute bin-centre frequency, Hz
    t: np.ndarray       # (T,) frame start time, s
    nfft: int
    hop: int
    fs: float
    fc: float

    @property
    def df(self) -> float:
        return self.fs / self.nfft

    @property
    def dt(self) -> float:
        return self.hop / self.fs

    @property
    def frame_len(self) -> float:
        return self.nfft / self.fs

    def f_index(self, f_hz):
        return np.clip(np.round((np.asarray(f_hz) - self.f[0]) / self.df).astype(int), 0, len(self.f) - 1)

    def t_index(self, t_s):
        return np.clip(np.round(np.asarray(t_s) / self.dt).astype(int), 0, len(self.t) - 1)


def stft(x: np.ndarray, fs: float, fc: float, nfft: int, hop: int | None = None) -> Spectrogram:
    hop = hop or nfft
    frames = sliding_window_view(x, nfft)[::hop]
    w = np.hanning(nfft).astype(np.float32)
    S = sfft.fft(frames * w, axis=1, workers=-1)
    P = (np.abs(sfft.fftshift(S, axes=1)) ** 2 / np.sum(w ** 2)).astype(np.float32)
    f = fc + (np.arange(nfft) - nfft // 2) * fs / nfft
    t = np.arange(P.shape[0]) * hop / fs
    return Spectrogram(P, f, t, nfft, hop, fs, fc)
