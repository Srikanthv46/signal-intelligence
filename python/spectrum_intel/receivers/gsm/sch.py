"""GSM BCCH carrier receiver: FCCH tone search, SCH differential-GMSK demod, Viterbi, CRC-10,
BSIC / frame number fields."""
from __future__ import annotations

import numpy as np

from spectrum_intel.core.bands_india import gsm_arfcn

RSYM = 1625e3 / 6                    # 270.833 ksym/s
SPS = 4
OFS = RSYM * SPS
F_FCCH = RSYM / 4                    # +67.708 kHz
TS_SCH = np.array([1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
                   0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1],
                  np.uint8)
GSM_BANDS = {"b8_dl", "b3_dl"}


def _arfcn_freq(f):
    a = gsm_arfcn(f)
    if a is None:
        return None, None
    fr = 935e6 + 0.2e6 * a if a < 512 else 1805.2e6 + 0.2e6 * (a - 512)
    return a, fr


def find_fcch(y, win=int(0.4e-3 * OFS)):
    """Sliding fraction of power in the +67.7 kHz tone; returns (start, residual
    CFO = tone - 67.708 kHz, peak purity) or None."""
    n = np.arange(len(y))
    z = y * np.exp(-2j * np.pi * F_FCCH * n / OFS)
    c = np.convolve(z, np.ones(win), "valid")
    e = np.convolve(np.abs(y) ** 2, np.ones(win), "valid") * win
    r = np.abs(c) ** 2 / (e + 1e-30)
    k = int(np.argmax(r))
    if r[k] < 0.5:
        return None
    # the window lies fully inside the 148-symbol burst for k in [S, S + 148*SPS - win]:
    # S = start of the >= 95% purity plateau (exact burst start)
    above = np.where(r > 0.95 * r[k])[0]
    seg = above[(above >= k - 148 * SPS) & (above <= k + 148 * SPS)]
    s0 = int(seg[0])
    zz = z[k:k + win]
    cfo = float(np.angle(np.sum(zz[1:] * np.conj(zz[:-1]))) * OFS / (2 * np.pi))
    return s0, cfo, float(r[k])


FCCH_WIN = int(0.4e-3 * OFS)
FCCH_NFFT = 2048                     # 529 Hz bins at 4 sps; parabolic interpolation refines the peak


def _peak_hz(p, j, nfft=FCCH_NFFT):
    """Frequency of the power-spectrum peak at bin j, parabolic interpolation on log power."""
    la, lb, lc = (float(np.log(p[(j + d) % nfft] + 1e-30)) for d in (-1, 0, 1))
    den = la - 2 * lb + lc
    d = 0.5 * (la - lc) / den if den < 0 else 0.0
    f = j if j < nfft // 2 else j - nfft
    return (f + d) * OFS / nfft


def find_fcch_cfo(y, max_cfo):
    """FCCH with an unknown carrier offset (|cfo| <= max_cfo): (start, cfo, purity) or None.

    The tone is searched on 0.4 ms STFT frames (hop 0.1 ms, 2048-point FFT) as the largest fraction of frame power
    in one bin within F_FCCH +- max_cfo; find_fcch then gives the burst start and purity at that tone, and the CFO
    is the interpolated FFT peak of the 0.4 ms window from the burst start (the lag-1 phase of find_fcch is biased
    by residual adjacent-carrier power in the 220 kHz channel)."""
    win, hop = FCCH_WIN, FCCH_WIN // 4
    m = (len(y) - win) // hop + 1
    if m < 1:
        return None
    fr = y[np.arange(win)[None, :] + hop * np.arange(m)[:, None]]
    X = np.abs(np.fft.fft(fr, FCCH_NFFT, axis=1)) ** 2
    f = np.fft.fftfreq(FCCH_NFFT, 1 / OFS)
    sel = np.flatnonzero(np.abs(f - F_FCCH) <= max_cfo)
    e = np.sum(np.abs(fr) ** 2, axis=1) * win + 1e-30
    pur = X[:, sel] / e[:, None]
    i, j = np.unravel_index(int(np.argmax(pur)), pur.shape)
    if pur[i, j] < 0.3:                               # screen: the frame purity is near the burst's (>= 0.5 below)
        return None
    n = np.arange(len(y))
    cfo0 = _peak_hz(X[i], int(sel[j])) - F_FCCH
    y2 = y * np.exp(-2j * np.pi * cfo0 * n / OFS)
    a0 = max(0, i * hop - win)
    fc = find_fcch(y2[a0:i * hop + 2 * win + 148 * SPS])
    if fc is None:
        return None
    s0 = fc[0] + a0
    P = np.abs(np.fft.fft(y2[s0:s0 + win], FCCH_NFFT)) ** 2
    near = np.flatnonzero(np.abs(f - F_FCCH) <= 2 * OFS / FCCH_NFFT)
    k = int(near[np.argmax(P[near])])
    return s0, cfo0 + _peak_hz(P, k) - F_FCCH, fc[2]


def _conv_encode_tables():
    nxt = np.zeros((16, 2), int)
    out = np.zeros((16, 2, 2), int)
    for s in range(16):
        for b in (0, 1):
            reg = (b << 4) | s                       # reg bit4 = D^0 (new), bit0 = D^4
            d = [(reg >> (4 - i)) & 1 for i in range(5)]  # d[i] = D^i
            g0 = d[0] ^ d[3] ^ d[4]
            g1 = d[0] ^ d[1] ^ d[3] ^ d[4]
            nxt[s, b] = reg >> 1
            out[s, b] = (g0, g1)
    return nxt, out


_NXT, _OUT = _conv_encode_tables()


def viterbi(soft):
    """soft: 78 values, positive = bit 0. Returns 39 bits."""
    n = len(soft) // 2
    INF = 1e18
    pm = np.full(16, INF)
    pm[0] = 0
    hist = np.zeros((n, 16), int)
    hbit = np.zeros((n, 16), int)
    for i in range(n):
        s0, s1 = soft[2 * i], soft[2 * i + 1]
        npm = np.full(16, INF)
        for s in range(16):
            if pm[s] >= INF:
                continue
            for b in (0, 1):
                g0, g1 = _OUT[s, b]
                cost = (s0 if g0 else -s0) + (s1 if g1 else -s1)      # penalise disagreement
                t = _NXT[s, b]
                if pm[s] + cost < npm[t]:
                    npm[t] = pm[s] + cost
                    hist[i, t] = s
                    hbit[i, t] = b
        pm = npm
    s = 0                                                  # 4 zero tail bits terminate in state 0
    bits = []
    for i in range(n - 1, -1, -1):
        bits.append(hbit[i, s])
        s = hist[i, s]
    return np.array(bits[::-1], np.uint8)


def crc10_ok(bits25, par10):
    poly = 0b10101110101                                   # D10+D8+D6+D5+D4+D2+1
    reg = list(bits25) + [0] * 10
    for i in range(25):
        if reg[i]:
            for j in range(11):
                reg[i + j] ^= (poly >> (10 - j)) & 1
    rem = np.array(reg[25:], np.uint8) ^ 1                 # remainder is inverted
    return bool(np.array_equal(rem, np.asarray(par10, np.uint8)))


def sch_fields(b):
    """SCH information from the 25 decoded bits (LSB-first octet packing)."""
    by = [sum(int(b[8 * k + j]) << j for j in range(8) if 8 * k + j < 25) for k in range(4)]
    bsic = (by[0] >> 2) & 0x3F
    t1 = ((by[0] & 3) << 9) | (by[1] << 1) | (by[2] >> 7)
    t2 = (by[2] >> 2) & 0x1F
    t3p = ((by[2] & 3) << 1) | (by[3] & 1)
    t3 = 10 * t3p + 1
    fn = 51 * ((t3 - t2) % 26) + t3 + 51 * 26 * t1
    return {"bsic": bsic, "ncc": bsic >> 3, "bcc": bsic & 7, "t1": t1, "t2": t2, "t3_prime": t3p,
            "frame_number": fn, "fields_plausible": t2 <= 25 and t3p <= 4}


def demod_sch(y, start):
    """Differential GMSK detection of the SCH burst around `start` (sample index of
    the burst's first tail bit). Returns (bits39, info) or None."""
    best = None
    for off in range(-20 * SPS, 20 * SPS + 1):
        s = start + off
        if s < 0 or s + 150 * SPS >= len(y):
            continue
        r = y[s:s + 149 * SPS:SPS]
        dph = np.angle(r[1:] * np.conj(r[:-1]))            # ~ +-pi/2 per symbol
        for pol in (1, -1):
            a = pol * np.sign(dph)                         # +1 -> d=0
            d = (a < 0).astype(np.uint8)
            b = np.zeros(len(d) + 1, np.uint8)
            for i in range(len(d)):
                b[i + 1] = d[i] ^ b[i]
            tsb = b[42:106]                                # 3 tail + 39 data -> TS at 42
            m = int(np.sum(tsb == TS_SCH))
            # training sequence is matched up to the unknown reference bit (complement)
            m = max(m, 64 - m)
            if best is None or m > best[0]:
                best = (m, s, pol, dph, b)
    if best is None or best[0] < 56:
        return None
    m, s, pol, dph, b = best
    if np.sum(b[42:106] == TS_SCH) < 32:
        b = b ^ 1
    # soft values for the 78 coded bits (positive = 0): reliability from |phase step|
    rel = np.abs(np.abs(dph) - 0) / (np.pi / 2)
    rel = np.r_[1.0, np.clip(rel, 0.05, 1.5)]
    idx = np.r_[np.arange(3, 42), np.arange(106, 145)]
    soft = (1 - 2 * b[idx].astype(float)) * rel[idx]
    return soft, {"ts_match": m, "sample": s}
