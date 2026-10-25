"""3GPP LTE synchronisation sequences (PSS Zadoff-Chu, SSS m-sequences) shared by
the LTE cell search and the Herelink receiver."""
from __future__ import annotations

import numpy as np


# ============================================================== LTE primitives
def lte_pss_freq(u):
    """3GPP 36.211 6.11.1.1 PSS: 62-length ZC, root u in {25, 29, 34} <- N_ID2 0/1/2."""
    n = np.arange(62)
    return np.where(n <= 30, np.exp(-1j * np.pi * u * n * (n + 1) / 63),
                    np.exp(-1j * np.pi * u * (n + 1) * (n + 2) / 63))


def _mseq(fb):
    x = [0, 0, 0, 0, 1]
    for i in range(26):
        x.append(fb(x, i) % 2)
    return 1 - 2 * np.array(x)


_S = _mseq(lambda x, i: x[i + 2] + x[i])
_C = _mseq(lambda x, i: x[i + 3] + x[i])
_Z = _mseq(lambda x, i: x[i + 4] + x[i + 2] + x[i + 1] + x[i])


def lte_sss_freq(nid1, nid2, sub):
    """36.211 6.11.2.1 SSS for subframe 0 (sub=0) or 5 (sub=5)."""
    qp = nid1 // 30; q = (nid1 + qp * (qp + 1) // 2) // 30; mp = nid1 + q * (q + 1) // 2
    m0 = mp % 31; m1 = (m0 + mp // 31 + 1) % 31
    n = np.arange(31)
    s0, s1 = _S[(n + m0) % 31], _S[(n + m1) % 31]
    c0, c1 = _C[(n + nid2) % 31], _C[(n + nid2 + 3) % 31]
    z0, z1 = _Z[(n + m0 % 8) % 31], _Z[(n + m1 % 8) % 31]
    d = np.zeros(62)
    if sub == 0:
        d[0::2] = s0 * c0; d[1::2] = s1 * c1 * z0
    else:
        d[0::2] = s1 * c0; d[1::2] = s0 * c1 * z1
    return d


_SC62 = np.r_[np.arange(-31, 0), np.arange(1, 32)]          # PSS/SSS subcarriers (DC skipped)
_SSS_BANK = {nid2: np.array([[lte_sss_freq(n1, nid2, s) for n1 in range(168)] for s in (0, 5)])
             for nid2 in range(3)}                             # [sub, nid1, 62]


# ============================================================== 3GPP Gold sequence (36.211 7.2 / 38.211 5.2.1)
def _lfsr31(init: np.ndarray, taps: tuple, n: int) -> np.ndarray:
    """x(k+31) = sum_t x(k+t) mod 2. Every new bit depends on bits >= 28 places back
    (the largest tap is 3), so blocks of 28 are computed at once."""
    x = np.zeros(n + 31, np.int8)
    x[:31] = init
    k = 0
    while k < n:
        b = min(28, n - k)
        acc = x[k:k + b].copy()
        for t in taps[1:]:
            acc ^= x[k + t:k + t + b]
        x[k + 31:k + 31 + b] = acc
        k += b
    return x


def gold(c_init: int, n: int, nc: int = 1600) -> np.ndarray:
    """Pseudo-random sequence c(0..n-1) for c_init: x1 = 1,0,...,0 with x1(k+31) = x1(k+3)+x1(k),
    x2 = bits of c_init with x2(k+31) = x2(k+3)+x2(k+2)+x2(k+1)+x2(k), c(k) = x1(k+Nc)+x2(k+Nc).
    Bit-identical to py3gpp's _calc_gold / nrPRBS (int8 0/1)."""
    x1 = _lfsr31(np.r_[1, np.zeros(30, np.int8)], (0, 3), n + nc)
    x2 = _lfsr31(np.array([(c_init >> i) & 1 for i in range(31)], np.int8), (0, 1, 2, 3), n + nc)
    return (x1[nc:nc + n] ^ x2[nc:nc + n]).astype(np.int8)


def qpsk(bits: np.ndarray) -> np.ndarray:
    """38.211 5.1.3: ((1-2b(2i)) + j(1-2b(2i+1))) / sqrt(2)."""
    b = np.asarray(bits).reshape(-1, 2)
    return ((1 - 2 * b[:, 0]) + 1j * (1 - 2 * b[:, 1])) / np.sqrt(2)
