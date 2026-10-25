"""Synthetic DVB-T (EN 300 744) transmitter for testing the COFDM receiver: 2K / 8K, any guard interval, any channel
bandwidth B (elementary period T = 7 / (8 B) us, research/drone_transmitters/06_dvbt_cofdm.md), continual / scattered
pilots (4/3 boost, w_k), TPS (DBPSK, sync word per frame, BCH(67,53)), random QPSK / 16-QAM / 64-QAM data cells.
"""
from __future__ import annotations

from fractions import Fraction

import numpy as np
from scipy.signal import resample_poly

from . import tables as T


def _qam(m, n, rng):
    lv = {4: [-1, 1], 16: [-3, -1, 1, 3], 64: [-7, -5, -3, -1, 1, 3, 5, 7]}[m]
    lv = np.array(lv, float)
    return (rng.choice(lv, n) + 1j * rng.choice(lv, n)) / np.sqrt(2 * np.mean(lv ** 2))


def tps_block(frame: int, constellation=0, hierarchy=0, cr_hp=0, cr_lp=0, gi="1/8", mode="2K", cell_id=None) -> list:
    """s0..s67 of frame 1..4."""
    sync = list(T.SYNC) if frame in (1, 3) else [1 - b for b in T.SYNC]
    length = 0b011111 if cell_id is not None else 0b010111
    gi_code = {v: k for k, v in T.GI_TPS.items()}[gi]
    mode_code = {v: k for k, v in T.MODE_TPS.items()}[mode]
    cid = 0 if cell_id is None else ((cell_id >> 8) & 0xFF if frame in (1, 3) else cell_id & 0xFF)
    info = (sync + T.bits(length, 6) + T.bits(frame - 1, 2) + T.bits(constellation, 2) + T.bits(hierarchy, 3)
            + T.bits(cr_hp, 3) + T.bits(cr_lp, 3) + T.bits(gi_code, 2) + T.bits(mode_code, 2) + T.bits(cid, 8) + [0] * 6)
    return [0] + info + T.bch_parity(info)


def dvbt_baseband(n_sym, mode="2K", gi="1/8", qam=4, start_symbol=0, start_frame=1, rng=None, pilots=True, **tps):
    """Complex baseband at the elementary rate 1/T (N samples per useful symbol), carriers k' = k - (K-1)/2
    around DC. start_symbol: symbol index l (0..67) of the first symbol in its frame."""
    rng = rng or np.random.default_rng()
    M = T.MODES[mode]
    N, K = M["N"], M["K"]
    Ng = N // T.GUARD[gi]
    w = T.prbs(K)
    cpset, tpsset = np.array(M["cp"]), np.array(M["tps"])
    kc = (K - 1) // 2
    out = []
    frame, l = start_frame, start_symbol
    tps_bits = tps_block(frame, gi=gi, mode=mode, **tps)
    # DBPSK state of the TPS carriers at the first generated symbol of its frame (re-initialised at l = 0)
    ref = 2 * (0.5 - w[tpsset]).astype(float)
    for b in tps_bits[1:l]:
        if b:
            ref = -ref
    for _ in range(n_sym):
        if l == 68:
            l, frame = 0, frame % 4 + 1
            tps_bits = tps_block(frame, gi=gi, mode=mode, **tps)
        if l == 0:
            ref = 2 * (0.5 - w[tpsset]).astype(float)
        elif tps_bits[l]:
            ref = -ref
        c = _qam(qam, K, rng)
        if not pilots:                                   # DVB-T numerology without its pilots / TPS (negative test)
            X = np.zeros(N, complex)
            X[(np.arange(K) - kc) % N] = c
            s = np.fft.ifft(X) * (N / np.sqrt(K))
            out.append(np.r_[s[-Ng:], s])
            l += 1
            continue
        sp = np.arange(3 * (l % 4), K, 12)
        c[sp] = 4 / 3 * 2 * (0.5 - w[sp])
        c[cpset] = 4 / 3 * 2 * (0.5 - w[cpset])
        c[tpsset] = ref
        X = np.zeros(N, complex)
        X[(np.arange(K) - kc) % N] = c
        s = np.fft.ifft(X) * (N / np.sqrt(K))                      # unit mean power
        out.append(np.r_[s[-Ng:], s])
        l += 1
    return np.concatenate(out)


def dvbt_capture_burst(n_sym, bw_mhz, fs, mode="2K", gi="1/8", **kw):
    """DVB-T signal at capture rate fs: the elementary rate 8 B / 7 MHz resampled to fs (rational)."""
    y = dvbt_baseband(n_sym, mode=mode, gi=gi, **kw)
    r = Fraction(fs / (8e6 * bw_mhz / 7)).limit_denominator(4000)
    return resample_poly(y, r.numerator, r.denominator).astype(np.complex64)
