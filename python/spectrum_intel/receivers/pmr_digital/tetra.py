"""TETRA downlink (ETSI EN 300 392-2): pi/4-DQPSK demodulation, synchronisation bursts, BSCH decoding.

    demod(y, fs)            -> (hard bits, soft bits, phase error) of a carrier sliced to 8 samples/symbol
    find_training(bits)     -> sync (38-bit y) and normal-burst (22-bit n/p) training-sequence hits
    decode_bsch(soft120)    -> (crc_ok, SYNC PDU fields) of one synchronisation burst's block 1

BSCH chain (8.2): type-5 bits -> descramble (e = 0: the colour code is not known yet) -> block de-interleave
(K = 120, a = 11) -> de-puncture the rate-2/3 RCPC code -> 16-state rate-1/4 mother code (K = 5) soft Viterbi
-> 60 information bits + CRC-CCITT (16) + 4 tail bits. The SYNC PDU (21.4.4.2) + D-MLE-SYNC carry the system
code, colour code, TN/FN/MN, MCC and MNC.
Depends on numpy only.
"""
from __future__ import annotations

import numpy as np

SR = 18000.0                        # symbol rate
SPS = 8
SLOT_BITS = 510                     # 255 symbols per timeslot
Y_SYNC = np.array([1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1,
                   0, 0, 1, 1, 1], np.int8)                                 # 9.4.4.3.4 synchronisation training
N_TS = np.array([1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0], np.int8)   # normal training n
P_TS = np.array([0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0], np.int8)   # normal training p
BSCH_BITS = 120                     # SB block 1 (type-5 bits), immediately before the sync training sequence
SYNC_MAX_ERR = 4
NTS_MAX_ERR = 2


# ------------------------------------------------------------------------------------------- demodulation
def rrc(beta: float, sps: int, span: int = 8) -> np.ndarray:
    t = np.arange(-span * sps, span * sps + 1) / sps
    den = 1 - (4 * beta * t) ** 2
    h = np.where(np.abs(den) < 1e-8, np.pi / 4 * np.sinc(1 / (4 * beta)),
                 np.sinc(t) * np.cos(np.pi * beta * t) / np.where(np.abs(den) < 1e-8, 1.0, den))
    return h / np.sqrt(np.sum(h ** 2))


_RRC = rrc(0.35, SPS)


def symbol_line(y: np.ndarray, fs: float) -> float:
    """18 ksym/s line of the matched-filtered envelope power over the median of its +-1..6 kHz neighbourhood."""
    z = np.convolve(y, _RRC, "same")
    e = np.abs(z) ** 2
    E = np.abs(np.fft.rfft(e - e.mean()))
    fe = np.fft.rfftfreq(len(e), 1 / fs)
    on, bg = np.abs(fe - SR) <= 3 * fs / len(e) + 5, (np.abs(fe - SR) > 1e3) & (np.abs(fe - SR) < 6e3)
    return float(E[on].max() / (np.median(E[bg]) + 1e-30)) if on.any() and bg.any() else 0.0


def demod(y: np.ndarray, quarter: int = 0, block_syms: int = 2000) -> tuple[np.ndarray, np.ndarray, float, float]:
    """y: carrier at SPS samples/symbol. Matched filter, carrier rotation removed with the 4th power of the 1-symbol
    differential (per block of block_syms symbols, which also re-picks the timing phase, so a sample-clock error is
    tracked), differential detection. The 4th power leaves the rotation ambiguous by quarter turns per symbol
    (4.5 kHz): `quarter` adds that many; the training sequences tell which is right.
    Returns (hard bits, soft bits (+ = 1), rms phase error in rad, carrier offset of y in Hz)."""
    z = np.convolve(y, _RRC, "same")
    d = z[SPS:] * np.conj(z[:-SPS])
    rot = (np.angle(np.mean(d ** 4)) - np.pi) / 4 + quarter * np.pi / 2
    z = z * np.exp(-1j * rot / SPS * np.arange(len(z)))
    out = []
    blk = block_syms * SPS
    for b0 in range(0, max(1, len(z) - blk // 4), blk):
        seg = z[b0:b0 + blk + SPS]
        ph = max(range(SPS), key=lambda p: float(np.mean(np.abs(seg[p::SPS]) ** 2)))
        v = seg[ph::SPS][:block_syms + 1]
        if len(v) < 2:
            continue
        dd = v[1:] * np.conj(v[:-1])
        r = (np.angle(np.mean(dd ** 4)) - np.pi) / 4
        out.append(dd * np.exp(-1j * r))
    dd = np.concatenate(out) if out else np.zeros(0, complex)
    ph = np.angle(dd)
    # dibits (5.5.2): 00 -> +pi/4, 01 -> +3pi/4, 10 -> -pi/4, 11 -> -3pi/4
    bits = np.empty(2 * len(ph), np.int8)
    bits[0::2], bits[1::2] = ph < 0, np.abs(ph) > np.pi / 2
    a = dd / (np.abs(dd) + 1e-12) * np.sqrt(np.abs(dd) / (np.median(np.abs(dd)) + 1e-12))
    soft = np.empty(2 * len(dd), np.float64)
    soft[0::2], soft[1::2] = -a.imag, -a.real
    err = float(np.sqrt(np.mean(np.angle(np.exp(4j * (ph - np.pi / 4))) ** 2)) / 4) if len(ph) else np.pi
    return bits, soft, err, float(rot * SR / (2 * np.pi))


def _hits(bits: np.ndarray, seq: np.ndarray, max_err: int) -> list[tuple[int, int]]:
    if len(bits) < len(seq):
        return []
    c = np.correlate(2 * bits.astype(np.float32) - 1, 2 * seq.astype(np.float32) - 1, "valid")
    err = (len(seq) - c) / 2
    return [(int(i), int(round(err[i]))) for i in np.flatnonzero(err <= max_err + 1e-6)]


def find_training(bits: np.ndarray) -> dict:
    """Sync (y) and normal (n, p) training-sequence hits; a TETRA downlink puts them at fixed slot phases."""
    return {"sync": _hits(bits, Y_SYNC, SYNC_MAX_ERR),
            "normal": _hits(bits, N_TS, NTS_MAX_ERR) + _hits(bits, P_TS, NTS_MAX_ERR)}


def slot_consistent(pos: list[int]) -> int:
    """Largest number of hits sharing one phase modulo the 510-bit slot."""
    if not pos:
        return 0
    return int(np.bincount(np.asarray(pos) % SLOT_BITS).max())


# ------------------------------------------------------------------------------------------- BSCH codec
_TAPS = (32, 26, 23, 22, 16, 12, 11, 10, 8, 7, 5, 4, 2, 1)       # 8.2.5 scrambling polynomial


def scrambling_sequence(n: int, e_bits=None) -> np.ndarray:
    """8.2.5: p(-31) = p(-30) = 1, p(-29..0) = extended colour code e(1..30); BSCH uses e = 0."""
    p = [1, 1] + list(e_bits if e_bits is not None else [0] * 30)
    out = np.empty(n, np.int8)
    for i in range(n):
        b = 0
        for t in _TAPS:
            b ^= p[-t]
        p.append(b)
        out[i] = b
    return out


_SCR_BSCH = scrambling_sequence(BSCH_BITS)
_DEINT = np.array([(11 * i) % 120 for i in range(1, 121)])        # 8.2.4.1: b4(k) = b3(i), k = 1 + (a i mod K)
_G = np.array([(1, 1, 0, 0, 1), (1, 0, 1, 1, 1), (1, 1, 1, 0, 1), (1, 1, 0, 1, 1)], np.int8)  # 8.2.3.1.1 g_i,j
_P23 = (1, 2, 5)                                                  # rate 2/3 puncturing, period 8


def _depuncture_index(n_type3: int) -> np.ndarray:
    j = np.arange(1, n_type3 + 1)
    q = (j - 1) // 3
    return 8 * q + np.array(_P23)[j - 3 * q - 1] - 1


_DEPUNCT = _depuncture_index(BSCH_BITS)


def _trellis():
    out = np.zeros((16, 2, 4), np.float64)
    nxt = np.zeros((16, 2), np.int16)
    for s in range(16):
        for u in (0, 1):
            reg = np.array([u] + [(s >> i) & 1 for i in range(4)], np.int8)
            out[s, u] = 2 * ((_G @ reg) % 2) - 1
            nxt[s, u] = ((s << 1) | u) & 0xF
    return out, nxt


_TR_OUT, _TR_NXT = _trellis()


def viterbi(v: np.ndarray, n: int) -> np.ndarray:
    """Soft Viterbi of the 16-state rate-1/4 mother code; v: 4n soft values (+ = 1, 0 = erased). Terminated in 0."""
    pm = np.full(16, -1e9)
    pm[0] = 0.0
    prev = np.zeros((n, 16), np.int16)
    pbit = np.zeros((n, 16), np.int8)
    for k in range(n):
        m = pm[:, None] + _TR_OUT @ v[4 * k:4 * k + 4]                  # [state, input]
        new = np.full(16, -1e9)
        for u in (0, 1):
            t = _TR_NXT[:, u]
            for s in range(16):
                if m[s, u] > new[t[s]]:
                    new[t[s]] = m[s, u]
                    prev[k, t[s]] = s
                    pbit[k, t[s]] = u
        pm = new
    s = 0
    bits = np.zeros(n, np.int8)
    for k in range(n - 1, -1, -1):
        bits[k] = pbit[k, s]
        s = prev[k, s]
    return bits


def crc16(bits) -> int:
    """8.2.2: CRC-CCITT (x^16 + x^12 + x^5 + 1), preset to ones, ones-complemented."""
    r = 0xFFFF
    for b in bits:
        fb = ((r >> 15) & 1) ^ int(b)
        r = (r << 1) & 0xFFFF
        if fb:
            r ^= 0x1021
    return r ^ 0xFFFF


_SYNC_PDU = (("system_code", 4), ("colour_code", 6), ("timeslot", 2), ("frame", 5), ("multiframe", 6),
             ("sharing_mode", 2), ("ts_reserved_frames", 3), ("u_plane_dtx", 1), ("frame18_extension", 1),
             ("reserved", 1), ("mcc", 10), ("mnc", 14), ("neighbour_cell_broadcast", 2), ("cell_service_level", 2),
             ("late_entry", 1))


def parse_sync_pdu(b: np.ndarray) -> dict:
    f, p = {}, 0
    for name, n in _SYNC_PDU:
        f[name] = int("".join(map(str, b[p:p + n])), 2)
        p += n
    f["timeslot"] += 1                       # TN coded 0..3 for timeslots 1..4
    f.pop("reserved")
    return f


def decode_bsch(soft120: np.ndarray) -> tuple[bool, dict | None]:
    """soft120: the 120 type-5 soft bits (+ = 1) of a synchronisation burst's block 1."""
    b4 = np.where(_SCR_BSCH == 1, -soft120, soft120)
    b3 = b4[_DEINT]
    v = np.zeros(4 * 80, np.float64)
    v[_DEPUNCT] = b3
    b2 = viterbi(v, 80)
    info, crc = b2[:60], int("".join(map(str, b2[60:76])), 2)
    if crc16(info) != crc:
        return False, None
    return True, parse_sync_pdu(info)
