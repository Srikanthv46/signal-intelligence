"""Synthetic narrowband FM voice transmitters for tests: audio (tones or speech-like
syllables), CTCSS or DCS, deviation, carrier offset, keying and noise."""
from __future__ import annotations

import numpy as np
from scipy import signal as ss

from .subaudio import DCS_RATE, dcs_encode


def speechlike(n: int, fs: float, rng) -> np.ndarray:
    """Band-limited noise (300-3000 Hz) gated at a 4 Hz syllable rate: peaky, like speech."""
    v = ss.sosfilt(ss.butter(4, [300, 3000], "bandpass", fs=fs, output="sos"), rng.standard_normal(n))
    t = np.arange(n) / fs
    env = np.clip(np.sin(2 * np.pi * 4 * t + rng.uniform(0, 6)), 0, None) ** 2
    v *= env
    return v / (np.percentile(np.abs(v), 99.5) + 1e-12)


def dcs_wave(code: str, inverted: bool, n: int, fs: float) -> np.ndarray:
    w = dcs_encode(code)
    bits = np.array([(w >> i) & 1 for i in range(23)])
    if inverted:
        bits = 1 - bits
    k = (np.arange(n) / fs * DCS_RATE).astype(int) % 23
    x = 2.0 * bits[k] - 1
    return ss.sosfiltfilt(ss.butter(4, 250, fs=fs, output="sos"), x)


def nbfm(fs: float, dur: float, f_off: float, audio: str | float | None = "speech", audio_dev: float = 2000.0,
         ctcss: float | None = None, dcs: str | None = None, dcs_inverted: bool = False, sub_dev: float = 500.0,
         snr_db: float | None = 25.0, keyed: tuple[float, float] | None = None, seed: int = 0) -> np.ndarray:
    """Complex baseband at ``fs``; SNR is in a 12.5 kHz channel. ``audio``: "speech",
    a tone frequency (Hz), or None. ``keyed``: (t0, t1) seconds the carrier is on."""
    rng = np.random.default_rng(seed)
    fa = 48000.0
    na = int(dur * fa)
    m = np.zeros(na)
    if audio == "speech":
        m += audio_dev * speechlike(na, fa, rng)
    elif audio is not None:
        m += audio_dev * np.sin(2 * np.pi * float(audio) * np.arange(na) / fa)
    if ctcss:
        m += sub_dev * np.sin(2 * np.pi * ctcss * np.arange(na) / fa)
    if dcs:
        m += sub_dev * dcs_wave(dcs, dcs_inverted, na, fa)
    n = int(dur * fs)
    mi = np.interp(np.arange(n) / fs, np.arange(na) / fa, m)
    x = np.exp(1j * (2 * np.pi * np.cumsum(f_off + mi) / fs + rng.uniform(0, 6)))
    if keyed:
        t = np.arange(n) / fs
        x *= (t >= keyed[0]) & (t < keyed[1])
    if snr_db is not None:
        npow = 10 ** (-snr_db / 10) * fs / 12.5e3
        x = x + np.sqrt(npow / 2) * (rng.standard_normal(n) + 1j * rng.standard_normal(n))
    return x.astype(np.complex64)


def am(fs: float, dur: float, f_off: float, depth: float = 0.8, audio: str | float = "speech",
       snr_db: float | None = 25.0, keyed: tuple[float, float] | None = None, seed: int = 0) -> np.ndarray:
    """Full-carrier AM voice (aeronautical style) at ``fs``; SNR of the carrier in a 12.5 kHz channel."""
    rng = np.random.default_rng(seed)
    fa = 48000.0
    na = int(dur * fa)
    m = speechlike(na, fa, rng) if audio == "speech" else np.sin(2 * np.pi * float(audio) * np.arange(na) / fa)
    n = int(dur * fs)
    mi = np.interp(np.arange(n) / fs, np.arange(na) / fa, m)
    t = np.arange(n) / fs
    x = (1 + depth * mi) * np.exp(1j * (2 * np.pi * f_off * t + rng.uniform(0, 6)))
    if keyed:
        x *= (t >= keyed[0]) & (t < keyed[1])
    if snr_db is not None:
        npow = 10 ** (-snr_db / 10) * fs / 12.5e3
        x = x + np.sqrt(npow / 2) * (rng.standard_normal(n) + 1j * rng.standard_normal(n))
    return x.astype(np.complex64)


def afsk_packet(fs: float, octets: bytes, f_off: float, dev: float = 3000.0, n_flags: int = 40, lead_s: float = 0.02,
                snr_db: float | None = 25.0, seed: int = 0) -> np.ndarray:
    """An AX.25 frame as Bell 202 AFSK (1200 bit/s, mark 1200 / space 2200 Hz, continuous phase) on narrowband FM,
    with n_flags of TXDELAY and lead_s of unmodulated carrier before and after."""
    from spectrum_intel.core import hdlc as H
    rng = np.random.default_rng(seed)
    lv = H.nrzi(H.frame_bits(octets, n_flags))
    fa = 48000.0
    spb = fa / 1200.0
    na = int(len(lv) * spb)
    tone = np.where(np.array(lv)[(np.arange(na) / spb).astype(int)] == 1, 1200.0, 2200.0)
    audio = np.sin(2 * np.pi * np.cumsum(tone) / fa)
    lead = np.zeros(int(lead_s * fa))
    m = dev * np.concatenate([lead, audio, lead])
    n = int(len(m) / fa * fs)
    mi = np.interp(np.arange(n) / fs, np.arange(len(m)) / fa, m)
    x = np.exp(1j * (2 * np.pi * np.cumsum(f_off + mi) / fs + rng.uniform(0, 6)))
    if snr_db is not None:
        npow = 10 ** (-snr_db / 10) * fs / 12.5e3
        x = x + np.sqrt(npow / 2) * (rng.standard_normal(n) + 1j * rng.standard_normal(n))
    return x.astype(np.complex64)
