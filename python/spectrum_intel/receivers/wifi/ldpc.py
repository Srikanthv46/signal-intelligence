"""IEEE 802.11n/ac LDPC code (802.11-2020 19.3.11.7): parity-check matrices, PPDU encoding
parameters, a batched normalised min-sum decoder and a systematic encoder.

    params(length_bytes, ncbps, rate, m_stbc)  -> dict  N_CW, L_LDPC, shortening, puncturing, repetition, N_SYM
    decode_stream(llr, params, rate)            -> (info bits, per-codeword ok, per-codeword iterations)
    encode_stream(bits, params, rate)           -> coded bits in transmit order (reference transmitter)
    decode(llr, L, rate) / encode(u, L, rate)   -> one codeword

LLRs are log P(0)/P(1). The twelve prototype matrices (block lengths 648/1296/1944, rates 1/2, 2/3,
3/4, 5/6; Table F-1..F-3 of the standard) are as published in Saurabh Tavildar's open-source LDPC
library (github.com/tavildar/LDPC, LdpcC/WiFiLDPC.h, MIT licence, Copyright (c) Saurabh Tavildar).
Entry e >= 0 is the Z x Z identity cyclically shifted right by e (row i has its one in column
(i + e) mod Z); -1 is the zero block. That convention, the codeword layout (information bits, then
parity), shortening, puncturing and repetition were all confirmed on real frames: 802.11n MCS7 STBC
PPDUs of a Microhard link decode to CRC-32-valid MPDUs (testdata/D_20260211-161341).
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy import sparse

_BASE_TEXT = {
    (648, 1, 2): """
          0  -1  -1  -1   0   0  -1  -1   0  -1  -1   0   1   0  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1
         22   0  -1  -1  17  -1   0   0  12  -1  -1  -1  -1   0   0  -1  -1  -1  -1  -1  -1  -1  -1  -1
          6  -1   0  -1  10  -1  -1  -1  24  -1   0  -1  -1  -1   0   0  -1  -1  -1  -1  -1  -1  -1  -1
          2  -1  -1   0  20  -1  -1  -1  25   0  -1  -1  -1  -1  -1   0   0  -1  -1  -1  -1  -1  -1  -1
         23  -1  -1  -1   3  -1  -1  -1   0  -1   9  11  -1  -1  -1  -1   0   0  -1  -1  -1  -1  -1  -1
         24  -1  23   1  17  -1   3  -1  10  -1  -1  -1  -1  -1  -1  -1  -1   0   0  -1  -1  -1  -1  -1
         25  -1  -1  -1   8  -1  -1  -1   7  18  -1  -1   0  -1  -1  -1  -1  -1   0   0  -1  -1  -1  -1
         13  24  -1  -1   0  -1   8  -1   6  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   0   0  -1  -1  -1
          7  20  -1  16  22  10  -1  -1  23  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   0   0  -1  -1
         11  -1  -1  -1  19  -1  -1  -1  13  -1   3  17  -1  -1  -1  -1  -1  -1  -1  -1  -1   0   0  -1
         25  -1   8  -1  23  18  -1  14   9  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   0   0
          3  -1  -1  -1  16  -1  -1   2  25   5  -1  -1   1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   0
    """,
    (648, 2, 3): """
         25  26  14  -1  20  -1   2  -1   4  -1  -1   8  -1  16  -1  18   1   0  -1  -1  -1  -1  -1  -1
         10   9  15  11  -1   0  -1   1  -1  -1  18  -1   8  -1  10  -1  -1   0   0  -1  -1  -1  -1  -1
         16   2  20  26  21  -1   6  -1   1  26  -1   7  -1  -1  -1  -1  -1  -1   0   0  -1  -1  -1  -1
         10  13   5   0  -1   3  -1   7  -1  -1  26  -1  -1  13  -1  16  -1  -1  -1   0   0  -1  -1  -1
         23  14  24  -1  12  -1  19  -1  17  -1  -1  -1  20  -1  21  -1   0  -1  -1  -1   0   0  -1  -1
          6  22   9  20  -1  25  -1  17  -1   8  -1  14  -1  18  -1  -1  -1  -1  -1  -1  -1   0   0  -1
         14  23  21  11  20  -1  24  -1  18  -1  19  -1  -1  -1  -1  22  -1  -1  -1  -1  -1  -1   0   0
         17  11  11  20  -1  21  -1  26  -1   3  -1  -1  18  -1  26  -1   1  -1  -1  -1  -1  -1  -1   0
    """,
    (648, 3, 4): """
         16  17  22  24   9   3  14  -1   4   2   7  -1  26  -1   2  -1  21  -1   1   0  -1  -1  -1  -1
         25  12  12   3   3  26   6  21  -1  15  22  -1  15  -1   4  -1  -1  16  -1   0   0  -1  -1  -1
         25  18  26  16  22  23   9  -1   0  -1   4  -1   4  -1   8  23  11  -1  -1  -1   0   0  -1  -1
          9   7   0   1  17  -1  -1   7   3  -1   3  23  -1  16  -1  -1  21  -1   0  -1  -1   0   0  -1
         24   5  26   7   1  -1  -1  15  24  15  -1   8  -1  13  -1  13  -1  11  -1  -1  -1  -1   0   0
          2   2  19  14  24   1  15  19  -1  21  -1   2  -1  24  -1   3  -1   2   1  -1  -1  -1  -1   0
    """,
    (648, 5, 6): """
         17  13   8  21   9   3  18  12  10   0   4  15  19   2   5  10  26  19  13  13   1   0  -1  -1
          3  12  11  14  11  25   5  18   0   9   2  26  26  10  24   7  14  20   4   2  -1   0   0  -1
         22  16   4   3  10  21  12   5  21  14  19   5  -1   8   5  18  11   5   5  15   0  -1   0   0
          7   7  14  14   4  16  16  24  24  10   1   7  15   6  10  26   8  18  21  14   1  -1  -1   0
    """,
    (1296, 1, 2): """
         40  -1  -1  -1  22  -1  49  23  43  -1  -1  -1   1   0  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1
         50   1  -1  -1  48  35  -1  -1  13  -1  30  -1  -1   0   0  -1  -1  -1  -1  -1  -1  -1  -1  -1
         39  50  -1  -1   4  -1   2  -1  -1  -1  -1  49  -1  -1   0   0  -1  -1  -1  -1  -1  -1  -1  -1
         33  -1  -1  38  37  -1  -1   4   1  -1  -1  -1  -1  -1  -1   0   0  -1  -1  -1  -1  -1  -1  -1
         45  -1  -1  -1   0  22  -1  -1  20  42  -1  -1  -1  -1  -1  -1   0   0  -1  -1  -1  -1  -1  -1
         51  -1  -1  48  35  -1  -1  -1  44  -1  18  -1  -1  -1  -1  -1  -1   0   0  -1  -1  -1  -1  -1
         47  11  -1  -1  -1  17  -1  -1  51  -1  -1  -1   0  -1  -1  -1  -1  -1   0   0  -1  -1  -1  -1
          5  -1  25  -1   6  -1  45  -1  13  40  -1  -1  -1  -1  -1  -1  -1  -1  -1   0   0  -1  -1  -1
         33  -1  -1  34  24  -1  -1  -1  23  -1  -1  46  -1  -1  -1  -1  -1  -1  -1  -1   0   0  -1  -1
          1  -1  27  -1   1  -1  -1  -1  38  -1  44  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   0   0  -1
         -1  18  -1  -1  23  -1  -1   8   0  35  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   0   0
         49  -1  17  -1  30  -1  -1  -1  34  -1  -1  19   1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   0
    """,
    (1296, 2, 3): """
         39  31  22  43  -1  40   4  -1  11  -1  -1  50  -1  -1  -1   6   1   0  -1  -1  -1  -1  -1  -1
         25  52  41   2   6  -1  14  -1  34  -1  -1  -1  24  -1  37  -1  -1   0   0  -1  -1  -1  -1  -1
         43  31  29   0  21  -1  28  -1  -1   2  -1  -1   7  -1  17  -1  -1  -1   0   0  -1  -1  -1  -1
         20  33  48  -1   4  13  -1  26  -1  -1  22  -1  -1  46  42  -1  -1  -1  -1   0   0  -1  -1  -1
         45   7  18  51  12  25  -1  -1  -1  50  -1  -1   5  -1  -1  -1   0  -1  -1  -1   0   0  -1  -1
         35  40  32  16   5  -1  -1  18  -1  -1  43  51  -1  32  -1  -1  -1  -1  -1  -1  -1   0   0  -1
          9  24  13  22  28  -1  -1  37  -1  -1  25  -1  -1  52  -1  13  -1  -1  -1  -1  -1  -1   0   0
         32  22   4  21  16  -1  -1  -1  27  28  -1  38  -1  -1  -1   8   1  -1  -1  -1  -1  -1  -1   0
    """,
    (1296, 3, 4): """
         39  40  51  41   3  29   8  36  -1  14  -1   6  -1  33  -1  11  -1   4   1   0  -1  -1  -1  -1
         48  21  47   9  48  35  51  -1  38  -1  28  -1  34  -1  50  -1  50  -1  -1   0   0  -1  -1  -1
         30  39  28  42  50  39   5  17  -1   6  -1  18  -1  20  -1  15  -1  40  -1  -1   0   0  -1  -1
         29   0   1  43  36  30  47  -1  49  -1  47  -1   3  -1  35  -1  34  -1   0  -1  -1   0   0  -1
          1  32  11  23  10  44  12   7  -1  48  -1   4  -1   9  -1  17  -1  16  -1  -1  -1  -1   0   0
         13   7  15  47  23  16  47  -1  43  -1  29  -1  52  -1   2  -1  53  -1   1  -1  -1  -1  -1   0
    """,
    (1296, 5, 6): """
         48  29  37  52   2  16   6  14  53  31  34   5  18  42  53  31  45  -1  46  52   1   0  -1  -1
         17   4  30   7  43  11  24   6  14  21   6  39  17  40  47   7  15  41  19  -1  -1   0   0  -1
          7   2  51  31  46  23  16  11  53  40  10   7  46  53  33  35  -1  25  35  38   0  -1   0   0
         19  48  41   1  10   7  36  47   5  29  52  52  31  10  26   6   3   2  -1  51   1  -1  -1   0
    """,
    (1944, 1, 2): """
         57  -1  -1  -1  50  -1  11  -1  50  -1  79  -1   1   0  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1
          3  -1  28  -1   0  -1  -1  -1  55   7  -1  -1  -1   0   0  -1  -1  -1  -1  -1  -1  -1  -1  -1
         30  -1  -1  -1  24  37  -1  -1  56  14  -1  -1  -1  -1   0   0  -1  -1  -1  -1  -1  -1  -1  -1
         62  53  -1  -1  53  -1  -1   3  35  -1  -1  -1  -1  -1  -1   0   0  -1  -1  -1  -1  -1  -1  -1
         40  -1  -1  20  66  -1  -1  22  28  -1  -1  -1  -1  -1  -1  -1   0   0  -1  -1  -1  -1  -1  -1
          0  -1  -1  -1   8  -1  42  -1  50  -1  -1   8  -1  -1  -1  -1  -1   0   0  -1  -1  -1  -1  -1
         69  79  79  -1  -1  -1  56  -1  52  -1  -1  -1   0  -1  -1  -1  -1  -1   0   0  -1  -1  -1  -1
         65  -1  -1  -1  38  57  -1  -1  72  -1  27  -1  -1  -1  -1  -1  -1  -1  -1   0   0  -1  -1  -1
         64  -1  -1  -1  14  52  -1  -1  30  -1  -1  32  -1  -1  -1  -1  -1  -1  -1  -1   0   0  -1  -1
         -1  45  -1  70   0  -1  -1  -1  77   9  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   0   0  -1
          2  56  -1  57  35  -1  -1  -1  -1  -1  12  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   0   0
         24  -1  61  -1  60  -1  -1  27  51  -1  -1  16   1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1   0
    """,
    (1944, 2, 3): """
         61  75   4  63  56  -1  -1  -1  -1  -1  -1   8  -1   2  17  25   1   0  -1  -1  -1  -1  -1  -1
         56  74  77  20  -1  -1  -1  64  24   4  67  -1   7  -1  -1  -1  -1   0   0  -1  -1  -1  -1  -1
         28  21  68  10   7  14  65  -1  -1  -1  23  -1  -1  -1  75  -1  -1  -1   0   0  -1  -1  -1  -1
         48  38  43  78  76  -1  -1  -1  -1   5  36  -1  15  72  -1  -1  -1  -1  -1   0   0  -1  -1  -1
         40   2  53  25  -1  52  62  -1  20  -1  -1  44  -1  -1  -1  -1   0  -1  -1  -1   0   0  -1  -1
         69  23  64  10  22  -1  21  -1  -1  -1  -1  -1  68  23  29  -1  -1  -1  -1  -1  -1   0   0  -1
         12   0  68  20  55  61  -1  40  -1  -1  -1  52  -1  -1  -1  44  -1  -1  -1  -1  -1  -1   0   0
         58   8  34  64  78  -1  -1  11  78  24  -1  -1  -1  -1  -1  58   1  -1  -1  -1  -1  -1  -1   0
    """,
    (1944, 3, 4): """
         48  29  28  39   9  61  -1  -1  -1  63  45  80  -1  -1  -1  37  32  22   1   0  -1  -1  -1  -1
          4  49  42  48  11  30  -1  -1  -1  49  17  41  37  15  -1  54  -1  -1  -1   0   0  -1  -1  -1
         35  76  78  51  37  35  21  -1  17  64  -1  -1  -1  59   7  -1  -1  32  -1  -1   0   0  -1  -1
          9  65  44   9  54  56  73  34  42  -1  -1  -1  35  -1  -1  -1  46  39   0  -1  -1   0   0  -1
          3  62   7  80  68  26  -1  80  55  -1  36  -1  26  -1   9  -1  72  -1  -1  -1  -1  -1   0   0
         26  75  33  21  69  59   3  38  -1  -1  -1  35  -1  62  36  26  -1  -1   1  -1  -1  -1  -1   0
    """,
    (1944, 5, 6): """
         13  48  80  66   4  74   7  30  76  52  37  60  -1  49  73  31  74  73  23  -1   1   0  -1  -1
         69  63  74  56  64  77  57  65   6  16  51  -1  64  -1  68   9  48  62  54  27  -1   0   0  -1
         51  15   0  80  24  25  42  54  44  71  71   9  67  35  -1  58  -1  29  -1  53   0  -1   0   0
         16  29  36  41  44  56  59  37  50  24  -1  65   4  65  52  -1   4  -1  73  52   1  -1  -1   0
    """,
}


def base_matrix(L: int, rate: tuple[int, int]) -> np.ndarray:
    p, q = rate
    rows = [[int(v) for v in ln.split()] for ln in _BASE_TEXT[(L, p, q)].strip().splitlines()]
    return np.array(rows, int)


class _Code:
    """One lifted code: edge lists for the decoder, dense parity generator for the encoder (lazy)."""

    def __init__(self, L: int, rate: tuple[int, int]):
        B = base_matrix(L, rate)
        Z = L // 24
        self.L, self.Z, self.M = L, Z, B.shape[0] * Z
        self.K = L - self.M
        rows, cols = [], []
        for i, j in zip(*np.nonzero(B >= 0)):
            t = np.arange(Z)
            rows.append(i * Z + t)
            cols.append(j * Z + (t + B[i, j]) % Z)
        r, c = np.concatenate(rows), np.concatenate(cols)
        order = np.lexsort((c, r))
        r, c = r[order], c[order]
        deg = np.bincount(r, minlength=self.M)
        self.dmax = int(deg.max())
        # check-major dense layout, padded with a dummy variable (index L) that is always a sure 0
        self.idx = np.full((self.M, self.dmax), L, int)
        pos = np.arange(len(r)) - np.repeat(np.cumsum(deg) - deg, deg)
        self.idx[r, pos] = c
        self.pad = self.idx == L
        self.E = sparse.csr_matrix((np.ones(self.M * self.dmax), (np.arange(self.M * self.dmax),
                                                                     self.idx.reshape(-1))),
                                   shape=(self.M * self.dmax, L + 1))
        self._r, self._c = r, c
        self._gen = None

    def syndrome_ok(self, hard: np.ndarray) -> np.ndarray:
        h = np.concatenate([hard, np.zeros(hard.shape[:-1] + (1,), hard.dtype)], -1)
        return ~(np.bitwise_xor.reduce(h[..., self.idx], axis=-1).any(-1))

    def generator(self) -> np.ndarray:
        """P such that parity = P @ info (mod 2): Gauss-Jordan on [H_p | H_i] over GF(2)."""
        if self._gen is None:
            H = np.zeros((self.M, self.L), np.uint8)
            H[self._r, self._c] = 1
            A = np.concatenate([H[:, self.K:], H[:, :self.K]], 1)
            for col in range(self.M):
                piv = col + int(np.flatnonzero(A[col:, col])[0])
                if piv != col:
                    A[[col, piv]] = A[[piv, col]]
                hit = np.flatnonzero(A[:, col])
                hit = hit[hit != col]
                A[hit] ^= A[col]
            self._gen = A[:, self.M:]
        return self._gen


@lru_cache(maxsize=None)
def code(L: int, rate: tuple[int, int]) -> _Code:
    return _Code(L, tuple(rate))


def encode(u, L: int, rate) -> np.ndarray:
    cd = code(L, tuple(rate))
    u = np.asarray(u, np.uint8)
    return np.concatenate([u, (cd.generator().astype(np.int64) @ u % 2).astype(np.uint8)])


def decode(llr, L: int, rate, iters: int = 60, alpha: float = 0.8):
    """Normalised min-sum, flooding schedule, batched: llr is (L,) or (B, L).
    Returns (hard bits, syndrome satisfied, iterations used), stopping each codeword at its first
    satisfied syndrome."""
    cd = code(L, tuple(rate))
    one = np.ndim(llr) == 1
    lam = np.atleast_2d(np.asarray(llr, float))
    nb = lam.shape[0]
    lam = np.concatenate([np.clip(lam, -1e4, 1e4), np.full((nb, 1), 1e6)], 1)   # dummy var: sure 0
    R = np.zeros((nb, cd.M, cd.dmax))
    done = np.zeros(nb, bool)
    out = (lam[:, :L] < 0).astype(np.uint8)
    used = np.full(nb, iters)
    for it in range(iters + 1):
        tot = lam + (R.reshape(nb, -1) @ cd.E)
        hard = (tot[:, :L] < 0).astype(np.uint8)
        ok = cd.syndrome_ok(hard)
        new = ok & ~done
        out[new], used[new] = hard[new], it
        done |= ok
        if done.all() or it == iters:
            break
        Q = tot[:, cd.idx] - R
        Q[:, cd.pad] = 1e6
        a, s = np.abs(Q), np.where(Q < 0, -1.0, 1.0)
        sp = np.prod(s, axis=2, keepdims=True)
        am = np.argmin(a, axis=2)
        m1 = np.take_along_axis(a, am[..., None], 2)
        np.put_along_axis(a, am[..., None], np.inf, 2)
        m2 = a.min(axis=2, keepdims=True)
        first = np.zeros(a.shape, bool)
        np.put_along_axis(first, am[..., None], True, 2)
        R = alpha * sp * s * np.where(first, m2, m1)
        R[:, cd.pad] = 0.0
    out[~done] = (tot[~done, :L] < 0)
    return (out[0], bool(done[0]), int(used[0])) if one else (out, done, used)


# ------------------------------------------------------------------ PPDU encoding (19.3.11.7.5)
def params(length_bytes: int, ncbps: int, rate, m_stbc: int = 1) -> dict:
    """Codeword count and size, shortening, puncturing and repetition for an HT/VHT LDPC PSDU."""
    p, q = rate
    R = p / q
    npld = 8 * length_bytes + 16
    navbits = ncbps * m_stbc * int(np.ceil(npld / (ncbps * R * m_stbc)))
    if navbits <= 648:
        ncw, L = 1, (1296 if navbits >= npld + 912 * (1 - R) else 648)
    elif navbits <= 1296:
        ncw, L = 1, (1944 if navbits >= npld + 1464 * (1 - R) else 1296)
    elif navbits <= 1944:
        ncw, L = 1, 1944
    elif navbits <= 2592:
        ncw, L = 2, (1944 if navbits >= npld + 2916 * (1 - R) else 1296)
    else:
        ncw, L = int(np.ceil(npld / (1944 * R))), 1944
    K = L * p // q
    nshrt = max(0, ncw * K - npld)
    npunc = max(0, ncw * L - navbits - nshrt)
    if ((npunc > 0.1 * ncw * L * (1 - R) and nshrt < 1.2 * npunc * R / (1 - R))
            or npunc > 0.3 * ncw * L * (1 - R)):
        navbits += ncbps * m_stbc
        npunc = max(0, ncw * L - navbits - nshrt)
    nrep = max(0, int(round(navbits - ncw * L * (1 - R) - npld)))
    return dict(npld=npld, navbits=navbits, nsym=navbits // ncbps, ncw=ncw, L=L, K=K,
                nshrt=nshrt, npunc=npunc, nrep=nrep)


def _split(n: int, ncw: int) -> list[int]:
    b = n // ncw
    return [b + (1 if i < n % ncw else 0) for i in range(ncw)]


def _layout(P):
    """Per codeword: (information bits sent, parity bits sent, repeated bits)."""
    ncw, L, K = P["ncw"], P["L"], P["K"]
    return [(K - s, L - K - u, r) for s, u, r in zip(_split(P["nshrt"], ncw), _split(P["npunc"], ncw),
                                                    _split(P["nrep"], ncw))]


def encode_stream(bits, P: dict, rate) -> np.ndarray:
    """Scrambled SERVICE+PSDU bits (length npld) -> coded bits in transmit order (length navbits)."""
    bits = np.asarray(bits, np.uint8)
    out, pos = [], 0
    for k_i, par, rep in _layout(P):
        u = np.zeros(P["K"], np.uint8)
        u[:k_i] = bits[pos:pos + k_i]
        pos += k_i
        cw = encode(u, P["L"], rate)
        tx = np.concatenate([cw[:k_i], cw[P["K"]:P["K"] + par]])
        out.append(np.concatenate([tx, np.resize(tx, rep)]) if rep else tx)
    return np.concatenate(out)


def transmitted(cw, P: dict, i: int) -> np.ndarray:
    """The bits codeword i (full, length L) puts on the air, in transmit order."""
    k_i, par, rep = _layout(P)[i]
    tx = np.concatenate([cw[:k_i], cw[P["K"]:P["K"] + par]])
    return np.concatenate([tx, np.resize(tx, rep)]) if rep else tx


def offsets(P: dict) -> list[int]:
    """Start of each codeword's transmitted bits in the coded stream."""
    return list(np.cumsum([0] + [k + p + r for k, p, r in _layout(P)])[:-1])


def decode_stream(llr, P: dict, rate, iters: int = 60, full: bool = False):
    """LLRs of the coded bits in transmit order -> (information bits, ok per codeword, iterations)
    [+ the full decoded codewords (ncw, L) when full=True].
    Shortened bits are known zeros, punctured bits are erasures, repeated bits are combined."""
    llr = np.asarray(llr, float)
    L, K = P["L"], P["K"]
    lay = _layout(P)
    batch, pos = np.zeros((len(lay), L)), 0
    for i, (k_i, par, rep) in enumerate(lay):
        seg = llr[pos:pos + k_i + par + rep]
        pos += k_i + par + rep
        seg = np.pad(seg, (0, k_i + par + rep - len(seg)))
        cw = batch[i]
        cw[:k_i] = seg[:k_i]
        cw[k_i:K] = 1e3
        cw[K:K + par] = seg[k_i:k_i + par]
        sent = np.concatenate([np.arange(k_i), np.arange(K, K + par)])
        if rep:
            np.add.at(cw, np.resize(sent, rep), seg[k_i + par:])
    hard, ok, used = decode(batch, L, rate, iters=iters)
    info = np.concatenate([hard[i, :k_i] for i, (k_i, _, _) in enumerate(lay)])
    return (info, ok, used, hard) if full else (info, ok, used)
