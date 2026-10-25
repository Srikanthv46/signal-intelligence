"""Synthetic 315 / 433 MHz devices (tests): OOK pulse trains and FSK TPMS frames at an offset, unit noise."""
from __future__ import annotations

import numpy as np

from . import tpms as T


def ev1527(code24: int, te: float = 350.0, repeats: int = 4) -> list[tuple[float, float]]:
    """(mark us, space us) of `repeats` frames: sync then 24 PWM bits."""
    frame = [(te, 31 * te)] + [((3 * te, te) if (code24 >> (23 - q)) & 1 else (te, 3 * te)) for q in range(24)]
    return frame * repeats


def pt2262(tristate: str, alpha: float = 100.0, repeats: int = 4) -> list[tuple[float, float]]:
    half = {"0": (0, 0), "1": (1, 1), "F": (0, 1)}
    word = []
    for ch in tristate:
        for b in half[ch]:
            word.append((12 * alpha, 4 * alpha) if b else (4 * alpha, 12 * alpha))
    return ([(4 * alpha, 124 * alpha)] + word) * repeats


def keeloq(serial: int, button: int, hop: int, te: float = 400.0, packets: int = 2, vlow: int = 0, rpt: int = 0):
    bits = [(hop >> q) & 1 for q in range(32)] + [(serial >> q) & 1 for q in range(28)]
    bits += [(button >> (3 - q)) & 1 for q in range(4)] + [vlow, rpt]
    out = []
    for _ in range(packets):
        out += [(te, te)] * 11 + [(te, 10 * te)]
        out += [((te, 2 * te) if b else (2 * te, te)) for b in bits[:-1]]
        out.append(((te, 39 * te) if bits[-1] else (2 * te, 39 * te)))
    return out


def ook(fs: float, dur: float, f_off: float, train: list[tuple[float, float]], t0_us: float = 500.0,
        snr_db: float = 20.0, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n = int(dur * fs)
    env = np.zeros(n)
    t = t0_us
    for m, s in train:
        a, b = int(t * 1e-6 * fs), int((t + m) * 1e-6 * fs)
        env[a:min(b, n)] = 1.0
        t += m + s
    k = max(1, int(10e-6 * fs))                  # 10 us rise / fall, as a real transmitter's keying
    env = np.convolve(env, np.ones(k) / k, mode="same")
    x = 10 ** (snr_db / 20) * env * np.exp(2j * np.pi * f_off * np.arange(n) / fs)
    return (x + np.sqrt(0.5) * (rng.standard_normal(n) + 1j * rng.standard_normal(n))).astype(np.complex64)


def vdo_symbols(state: int, sid: int, flags: int, rep: int, p_raw: int, t_raw: int, batt: int) -> list[int]:
    b = bytes([state, sid >> 24 & 255, sid >> 16 & 255, sid >> 8 & 255, sid & 255, (flags << 4) | rep, p_raw, t_raw,
               batt])
    b += bytes([T.crc8(b, 0x07, 0xAA)])
    inv = list(T.VDO_PRE)
    for byte in b:
        for q in range(8):
            v = (byte >> (7 - q)) & 1
            inv += [1 - v, v]
    return [1 - v for v in inv] + [0, 1] * 4


def toyota_symbols(sid: int, status: int, p1: int, temp: int) -> list[int]:
    b = bytes([sid >> 24 & 255, sid >> 16 & 255, sid >> 8 & 255, sid & 255, (status & 0x80) | (p1 >> 1),
               ((p1 & 1) << 7) | (temp >> 1), ((temp & 1) << 7) | (status & 0x7F), p1 ^ 0xFF])
    b += bytes([T.crc8(b, 0x07, 0x80)])
    s = [1, 0] * 4 + T.TOYOTA_PRE[:-1]
    prev = s[-1]
    for byte in b:
        for q in range(8):
            v = (byte >> (7 - q)) & 1
            h1 = 1 - prev
            h2 = h1 if v else 1 - h1
            s += [h1, h2]
            prev = h2
    s += [1 - prev, 1 - prev, prev]
    return s


def fsk(fs: float, dur: float, f_off: float, bursts: list[tuple[float, list[int]]], dev: float = 40e3,
        sym_us: float = 52.0, snr_db: float = 20.0, seed: int = 0) -> np.ndarray:
    """bursts: (start us, symbols); a 1 symbol is +dev."""
    rng = np.random.default_rng(seed)
    n = int(dur * fs)
    fi = np.zeros(n)
    on = np.zeros(n)
    for t0, sy in bursts:
        for q, v in enumerate(sy):
            a, b = int((t0 + q * sym_us) * 1e-6 * fs), int((t0 + (q + 1) * sym_us) * 1e-6 * fs)
            fi[a:min(b, n)] = dev if v else -dev
            on[a:min(b, n)] = 1.0
    x = 10 ** (snr_db / 20) * on * np.exp(1j * (2 * np.pi * np.cumsum(f_off + fi) / fs))
    return (x + np.sqrt(0.5) * (rng.standard_normal(n) + 1j * rng.standard_normal(n))).astype(np.complex64)
