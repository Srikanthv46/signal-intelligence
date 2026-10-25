"""
rs.py -- Reed-Solomon codes over GF(2^8), systematic and shortened (clean-room implementation: field from the primitive
polynomial, generator g(x) = prod (x - alpha^(prim (fcr + i))), syndromes, Berlekamp-Massey, Chien search, Forney).
Parameters as in the common (symsize, gfpoly, fcr, prim, nroots) convention; data first, parity after it, the first
byte is the highest-degree coefficient. Used by UAT 978 (gfpoly 0x187, fcr 120, prim 1, nroots 12 / 14 / 20).
"""
from __future__ import annotations


class RS:
    def __init__(self, gfpoly: int, fcr: int, prim: int, nroots: int):
        self.nroots, self.fcr, self.prim = nroots, fcr, prim
        self.exp = [0] * 512
        self.log = [0] * 256
        x = 1
        for i in range(255):
            self.exp[i] = x
            self.log[x] = i
            x <<= 1
            if x & 0x100:
                x ^= gfpoly
        for i in range(255, 512):
            self.exp[i] = self.exp[i - 255]
        # prim-th root of alpha (prim = 1: alpha itself); iprim = 1 / prim mod 255
        self.iprim = next(k for k in range(1, 256) if (k * prim) % 255 == 1)
        g = [1]
        for i in range(nroots):
            r = self.exp[(prim * (fcr + i)) % 255]
            ng = [0] * (len(g) + 1)
            for j, c in enumerate(g):
                ng[j] ^= c
                ng[j + 1] ^= self.mul(c, r)
            g = ng
        self.gen = g                                           # highest degree first, g[0] = 1

    def mul(self, a, b):
        return 0 if a == 0 or b == 0 else self.exp[self.log[a] + self.log[b]]

    def encode(self, data: bytes) -> bytes:
        """Parity bytes of data (systematic: codeword = data + parity)."""
        par = [0] * self.nroots
        for d in data:
            fb = d ^ par[0]
            par = par[1:] + [0]
            if fb:
                for j in range(self.nroots):
                    par[j] ^= self.mul(fb, self.gen[j + 1])
        return bytes(par)

    def decode(self, block: bytes):
        """Codeword (data + parity, shortened) -> (corrected block, number of corrected bytes), or (None, -1) when the
        errors exceed the code's capacity."""
        n = len(block)
        r = list(block)
        nr = self.nroots
        # syndromes S_i = r(alpha^(prim (fcr + i)))
        S = []
        for i in range(nr):
            a = self.exp[(self.prim * (self.fcr + i)) % 255]
            s = 0
            for v in r:
                s = self.mul(s, a) ^ v
            S.append(s)
        if not any(S):
            return bytes(r), 0
        # Berlekamp-Massey: error locator Lambda(x), lowest degree first
        lam, b = [1] + [0] * nr, [1] + [0] * nr
        L, m = 0, 1
        bb = 1
        for k in range(nr):
            d = S[k]
            for i in range(1, L + 1):
                d ^= self.mul(lam[i], S[k - i])
            if d == 0:
                m += 1
                continue
            coef = self.mul(d, self.exp[255 - self.log[bb]])
            t = lam[:]
            for i in range(m, nr + 1):
                lam[i] ^= self.mul(coef, b[i - m])
            if 2 * L <= k:
                L, b, bb, m = k + 1 - L, t, d, 1
            else:
                m += 1
        deg = max(i for i, c in enumerate(lam) if c)
        if deg != L or L > nr // 2:
            return None, -1
        # Chien search over the codeword positions: position p (0 = first byte) has locator X = alpha^(prim (n-1-p))
        pos, Xs = [], []
        for p in range(n):
            e = (self.prim * (n - 1 - p)) % 255
            xinv = self.exp[(255 - e) % 255]
            v, xp = 0, 1
            for c in lam[:L + 1]:
                v ^= self.mul(c, xp)
                xp = self.mul(xp, xinv)
            if v == 0:
                pos.append(p)
                Xs.append(self.exp[e])
        if len(pos) != L:
            return None, -1
        # Forney: Omega(x) = S(x) Lambda(x) mod x^nr; e = X^(1 - fcr) Omega(X^-1) / Lambda'(X^-1)
        om = [0] * nr
        for i in range(nr):
            for j in range(min(i, L) + 1):
                om[i] ^= self.mul(S[i - j], lam[j])
        for p, X in zip(pos, Xs):
            xi = self.exp[255 - self.log[X]]
            num, xp = 0, 1
            for c in om:
                num ^= self.mul(c, xp)
                xp = self.mul(xp, xi)
            den, xp = 0, 1                                    # Lambda'(x): odd terms
            for i in range(1, L + 1, 2):
                den ^= self.mul(lam[i], self.exp[(self.log[xi] * (i - 1)) % 255] if i > 1 else 1)
            if den == 0:
                return None, -1
            # X^(1 - fcr) (prim = 1 convention; for prim != 1 the locator powers above already carry it)
            sc = self.exp[(self.log[X] * (1 - self.fcr)) % 255]
            r[p] ^= self.mul(self.mul(num, sc), self.exp[255 - self.log[den]])
        return bytes(r), len(pos)
