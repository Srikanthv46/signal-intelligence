"""Narrowband FM voice -> audio, one block at a time (live or file playback).

    demod = VoiceDemodulator.from_fields(signal_fields)     # report decoded.fields
    chunk = demod.process(iq, fs, fc)   # -> AudioChunk(samples float32 [-1, 1] at demod.audio_fs, meta) or None

Per block:
  1. slice the channel (+-channel_bw/2) out of one FFT of the block, at 3 x audio_fs;
  2. carrier squelch per 10 ms frame: channel power >= 8 dB over the noise density
     measured 15-100 kHz either side (median, so neighbouring carriers do not bias it);
  3. FM discriminator (Hz); 300-3400 Hz band-pass removes the CTCSS/DCS signalling;
     750 us de-emphasis; fixed gain from the deviation (no AGC, so silence stays silent);
  4. squelch gate with 5 ms ramps, decimate by 3.

Filter states carry across blocks (live blocks are near-contiguous); a block that does
not cover the channel returns None. Depends on numpy/scipy and spectrum_intel.core only.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import signal as ss

from spectrum_intel.core.dsp import Channelizer

DEEMPH_S = 750e-6
SQUELCH_DB = 8.0


@dataclass
class AudioChunk:
    samples: np.ndarray               # float32 mono, [-1, 1]
    fs: float
    meta: dict = field(default_factory=dict)


class VoiceDemodulator:
    def __init__(self, center_hz: float, deviation_hz: float = 2500.0, channel_bw_hz: float = 12.5e3,
                 audio_fs: int = 16000, squelch: bool = True, mode: str = "fm"):
        self.mode = mode                                          # "fm" | "am" (envelope, research)
        self.center_hz = float(center_hz)
        self.channel_bw_hz = float(channel_bw_hz)
        # level reference: the channel's rated deviation (2.5 kHz at 12.5 kHz spacing, 5 kHz at 25)
        self.ref_dev = max(float(deviation_hz), 2500.0 if channel_bw_hz <= 12.5e3 else 5000.0)
        self.audio_fs = int(audio_fs)
        self.squelch = squelch
        self._mid = 3 * self.audio_fs
        self._bp = ss.butter(4, [300, 3400], "bandpass", fs=self._mid, output="sos")
        self._zbp = np.zeros((self._bp.shape[0], 2))
        a = np.exp(-1 / (self._mid * DEEMPH_S))                   # one-pole de-emphasis
        self._de = ([1 - a], [1, -a])
        self._zde = np.zeros(1)
        self._dec = ss.firwin(63, 0.9 * self.audio_fs / 2, fs=self._mid)
        self._zdec = np.zeros(len(self._dec) - 1)
        g1k = 1 / np.sqrt(1 + (2 * np.pi * 1000 * DEEMPH_S) ** 2)
        self.gain = 0.5 / (self.ref_dev * g1k)                    # 1 kHz tone at full deviation -> 0.5 FS
        self._open = False

    @classmethod
    def from_fields(cls, f: dict, audio_fs: int = 16000) -> "VoiceDemodulator":
        """From a report signal's decoded.fields (voice_fm plugin)."""
        am = f.get("audio") == "am"
        return cls(center_hz=float(f["center_mhz"]) * 1e6, deviation_hz=float(f.get("deviation_hz") or 2500),
                   channel_bw_hz=(8e3 if am else float(f.get("channel_bw_khz") or 12.5) * 1e3), audio_fs=audio_fs,
                   mode="am" if am else "fm")

    def covers(self, fs: float, fc: float) -> bool:
        return abs(self.center_hz - fc) + self.channel_bw_hz / 2 < 0.48 * fs

    def process(self, iq: np.ndarray, fs: float, fc: float) -> AudioChunk | None:
        if not self.covers(fs, fc) or len(iq) < fs * 0.005:
            return None
        ch = Channelizer(np.asarray(iq, np.complex64), fs, fc)
        y, fs_y = ch.slice(self.center_hz, self.channel_bw_hz / 2, self._mid)
        # noise density beside the channel (median: robust to neighbouring carriers)
        off = np.abs(ch.f - (self.center_hz - fc))
        side = (off > 15e3) & (off < 100e3) & (np.abs(ch.f) < 0.48 * fs)
        nd = float(np.median(np.abs(ch.X[side]) ** 2)) / np.log(2) / (fs * len(iq)) if side.any() else 0.0
        n_frame = max(1, int(0.01 * fs_y))
        nf = len(y) // n_frame
        p = np.mean(np.abs(y[:nf * n_frame].reshape(nf, n_frame)) ** 2, 1) if nf else np.zeros(0)
        snr = 10 * np.log10(p / (nd * self.channel_bw_hz + 1e-30) + 1e-30)
        gate = np.ones(len(y)) if not self.squelch else self._gate(snr, n_frame, len(y))
        if self.mode == "am":
            # envelope; the band-pass removes the carrier (DC); 100 % modulation -> 0.5 FS
            e = np.abs(y).astype(np.float64)
            v, self._zbp = ss.sosfilt(self._bp, e, zi=self._zbp)
            v = v * (0.5 / (float(np.mean(e)) + 1e-30)) * gate
        else:
            d = np.angle(y[1:] * np.conj(y[:-1])) * (fs_y / (2 * np.pi))
            d = np.r_[d[:1], d] if len(d) else np.zeros(len(y))
            v, self._zbp = ss.sosfilt(self._bp, d, zi=self._zbp)
            v, self._zde = ss.lfilter(*self._de, v, zi=self._zde)
            v = v * self.gain * gate
        v, self._zdec = ss.lfilter(self._dec, 1, v, zi=self._zdec)
        out = np.clip(v[::3], -1, 1).astype(np.float32)
        return AudioChunk(out, fs_y / 3, {"squelch_open": float(np.mean(gate > 0.5)),
                                          "snr_db": round(float(np.median(snr)), 1) if len(snr) else None,
                                          "center_mhz": self.center_hz / 1e6})

    def _gate(self, snr_db: np.ndarray, n_frame: int, n: int) -> np.ndarray:
        g = np.zeros(len(snr_db))
        st = self._open
        for i, s in enumerate(snr_db):                            # 3 dB hysteresis
            st = s >= (SQUELCH_DB - 3 if st else SQUELCH_DB)
            g[i] = st
        self._open = bool(st)
        g = np.repeat(g, n_frame)
        g = np.r_[g, np.full(n - len(g), g[-1] if len(g) else 0.0)]
        r = max(1, int(0.005 * self._mid))
        return np.convolve(g, np.ones(r) / r, mode="same")       # 5 ms ramps: no clicks


def demodulate(iq: np.ndarray, fs: float, fc: float, fields: dict, audio_fs: int = 16000) -> AudioChunk | None:
    """Whole capture -> audio (one-shot)."""
    return VoiceDemodulator.from_fields(fields, audio_fs).process(iq, fs, fc)
