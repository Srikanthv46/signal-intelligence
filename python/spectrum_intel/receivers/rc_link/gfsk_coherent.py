"""
gfsk_coherent.py -- coherent detection of h ~ 0.5 GFSK (BLE 1M) for low SNR, plus CRC-24 Chase correction.

Bench: bench/ble_lowsnr (synthetic only; PLAN L1). BER at 8 dB in-band SNR (1 MHz): 2e-4 coherent vs 2.3e-2 one-bit
differential vs 5e-2 for the limiter-discriminator slicer; h 0.45 .. 0.55 and drift up to 100 Hz/us cost < 1 dB.

Laurent: h 0.5 GFSK ~ sum_k b_k C0(t - kT), b_k = b_{k-1} j alpha_k, so b_k^2 = (-1)^(k+1) b_0^2 whatever the data.
1. first pass: one-bit differential decisions (LPF, z[n+T] z*[n]) -> alpha_k; h = the value whose modulation-removed
   products are most coherent over 8 segments (ISI-aware per-bit phase table TAB);
2. matched filter C0; the h walk pi (h - 0.5) cumsum(alpha) and the matched-filter ISI (b_k (1 + j P1 (alpha_{k+1} -
   alpha_k))) removed; s_k = r_k^2 (-1)^(k+1) is a tone at 2 x CFO: its frequency track from segment FFT peaks
   (8, 4, then 1 segments; quadratic fit -> CFO + drift); timing from the strongest tone over +-3 samples;
3. phase from a +-W moving sum of s_k, decisions on the axis of each b_k (raw matched filter: the ISI is on the
   other axis), alpha_k = Im(b_k b*_(k-1)); a second iteration redoes 1-2 with these decisions.
A wrong b_k flips alpha_k and alpha_{k+1}: Chase correction flips those pairs at the least reliable b_k.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np
from scipy import signal

SPS = 8                                   # samples per bit (the BLE slice is 8 MSps)
# Laurent C0 of h 0.5 BT 0.5 GFSK (Gaussian pulse truncated to 3 bits), 8 samples/bit; regenerate with
# python/tests/test_gfsk_coherent.py::c0_pulse
C0 = np.array([
    1.5493647186e-08, 1.2515359415e-05, 1.1692026245e-04, 7.2607693471e-04, 3.3765011577e-03, 1.2201251491e-02,
    3.5143701476e-02, 8.2708419697e-02, 1.6331527480e-01, 2.7794917074e-01, 4.1790894564e-01, 5.6710940193e-01,
    7.0748264018e-01, 8.2411817347e-01, 9.0807272557e-01, 9.5695397637e-01, 9.7267765168e-01, 9.5676223260e-01,
    9.0767404231e-01, 8.2351009211e-01, 7.0669872824e-01, 5.6622299350e-01, 4.1702312740e-01, 2.7717237872e-01,
    1.6272616787e-01, 8.2330007259e-02, 3.4943224682e-02, 1.2116025896e-02, 3.3481098154e-03, 7.1881877861e-04,
    1.1552896953e-04, 1.2323470268e-05])
MF_DELAY = 9                              # matched-filter sample of bit k = bit start + 8 k - MF_DELAY
_AC = np.correlate(C0, C0, "full")
P1 = float(_AC[len(C0) - 1 + SPS] / _AC[len(C0) - 1])
H_LPF = signal.firwin(33, 0.6e6, fs=8e6)
# per-bit phase change (rad) of clean h 0.5 GFSK as the first pass measures it, for alpha_k = +1, by
# [previous bit == this bit][next bit == this bit]
TAB = np.array([[0.6064282430011425, 1.3034506592017188], [0.8962552450020563, 1.5725184171543038]])
H_GRID = np.round(np.arange(0.42, 0.581, 0.02), 2)
W = 16                                    # phase window +-W bits
AXIS = np.array([1j, -1, -1j, 1])               # j^(k+1) for k % 4
SIGN = np.array([-1.0, 1.0])                    # (-1)^(k+1) for k % 2


def _km(v, it=20):
    a, b = np.percentile(v, 25), np.percentile(v, 75)
    for _ in range(it):
        m = np.abs(v - a) < np.abs(v - b)
        if m.all() or (~m).all():
            break
        a, b = v[m].mean(), v[~m].mean()
    return min(a, b), max(a, b)


def diff_products(y):
    """LPF + one-bit differential products p[i] = z[i + T] z*[i] (phase change over the bit starting at i)."""
    z = np.convolve(np.r_[y, np.zeros(SPS, complex)], H_LPF, "same")
    return z[SPS:] * np.conj(z[:-SPS])


def find_access_address(p, aa_bits, top=3, min_ncorr=0.55):
    """Sample positions where the 32 access-address bits start: soft correlation of the differential phase with
    +-1 at every sample, normalised to [-1, 1]; up to `top` peaks >= min_ncorr, 32 bits apart."""
    d = np.angle(p)
    t = 2.0 * np.asarray(aa_bits, float) - 1
    L = SPS * (len(t) - 1) + 1
    if len(d) < L:
        return []
    m = len(d) - L + 1
    num, den = np.zeros(m), np.zeros(m)
    ad = np.abs(d)
    for r in range(SPS):                            # polyphase: the template is non-zero every SPS samples
        dr, ar = d[r::SPS], ad[r::SPS]
        nr = len(range(r, m, SPS))
        if nr <= 0:
            continue
        num[r::SPS] = np.correlate(dr, t, "valid")[:nr]
        den[r::SPS] = np.correlate(ar, np.ones(len(t)), "valid")[:nr]
    nc = num / (den + 1e-12)
    out = []
    for i in np.argsort(-nc, kind="stable"):
        if nc[i] < min_ncorr or len(out) >= top:
            break
        if all(abs(i - j) >= 32 * SPS for j, _ in out):
            out.append((int(i), float(nc[i])))
    return out


def _first_pass(p, s0, n, a_in=None):
    d = np.angle(p[s0 + np.arange(n) * SPS])
    a = np.where(d > np.mean(_km(d)), 1, -1)
    if a_in is not None:
        a = a_in
    an = np.r_[a[0], a, a[-1]]
    t = TAB[(an[:-2] == a).astype(int), (an[2:] == a).astype(int)] * a
    q0 = p[s0 + np.arange(n) * SPS]
    e = np.linspace(0, n, 9).astype(int)
    sc = np.array([sum(abs(q[e[i]:e[i + 1]].sum()) for i in range(8))
                   for q in (q0 * np.exp(-1j * (h / 0.5) * t) for h in H_GRID)])
    j = int(np.argmax(sc))
    h = float(H_GRID[j])
    if 0 < j < len(H_GRID) - 1:                    # parabolic interpolation between grid points
        den = sc[j - 1] - 2 * sc[j] + sc[j + 1]
        if den < 0:
            h += float(np.clip(0.5 * (sc[j - 1] - sc[j + 1]) / den, -0.5, 0.5)) * (H_GRID[1] - H_GRID[0])
    return a, h


def _tone_freq(s, pad=8):
    nf = 1 << int(np.ceil(np.log2(max(len(s), 16) * pad)))
    F = np.abs(np.fft.fft(s, nf))
    j = int(np.argmax(F))
    a, b, c = F[j - 1], F[j], F[(j + 1) % nf]
    dj = 0.5 * (a - c) / (a - 2 * b + c + 1e-30)
    return 2 * np.pi * (((j + dj) / nf + 0.5) % 1 - 0.5), float(b)


def _track(s, stages):
    n = len(s)
    ph = np.zeros(n)
    k = np.arange(n)
    for ns in stages:
        res = s * np.exp(-1j * ph)
        e = np.linspace(0, n, ns + 1).astype(int)
        fr, wt, kc = [], [], []
        for i in range(ns):
            if e[i + 1] - e[i] < 8:
                continue
            w, m = _tone_freq(res[e[i]:e[i + 1]])
            fr.append(w); wt.append(m); kc.append((e[i] + e[i + 1] - 1) / 2)
        if not fr:
            break
        deg = min(2, len(fr) - 1)
        fk = np.polyval(np.polyfit(kc, fr, deg, w=np.sqrt(wt)), k) if deg > 0 else np.full(n, fr[0])
        ph = ph + np.cumsum(fk) - fk[0]
    return ph


def _coherent(r, p, s0, n, a_in):
    a, h = _first_pass(p, s0, n, a_in)
    walk = np.pi * (h - 0.5) * np.cumsum(a)
    isi = 1 + 1j * P1 * (np.r_[a[1:], a[-1]] - a)
    best = None
    for dt in range(-3, 4):
        idx = s0 - MF_DELAY + dt + np.arange(n) * SPS
        rk = r[np.clip(idx, 0, len(r) - 1)] * np.exp(-1j * walk)
        s = (rk / isi) ** 2 * SIGN[np.arange(n) % 2]
        ph = _track(s, (8, 4))
        m = _tone_freq(s * np.exp(-1j * ph))[1]
        if best is None or m > best[0]:
            best = (m, dt, rk, s, ph)
    _, dt, rk, s, ph = best
    ph = ph + _track(s * np.exp(-1j * ph), (1,))
    cs = np.r_[0, np.cumsum(s * np.exp(-1j * ph))]
    k = np.arange(n)
    ph2 = np.unwrap(np.angle(cs[np.clip(k + W + 1, 0, n)] - cs[np.clip(k - W, 0, n)])) + ph
    v = rk * np.exp(-0.5j * ph2)
    ax = AXIS[k % 4]                                # j^(k+1), exact
    c = np.real(v * np.conj(ax))
    bh = np.where(c >= 0, 1, -1) * ax
    alpha = np.r_[a[0], np.where(np.imag(bh[1:] * np.conj(bh[:-1])) > 0, 1, -1)]
    return alpha, np.abs(c), h, dt, ph


def demodulate(y, s0, n, iters=2):
    """Bits 0..n-1 of the h ~ 0.5 GFSK burst whose bit 0 starts at sample s0 of y (SPS samples / bit).
    Returns (bits uint8, reliability of each b_k, h, cfo_hz_at_mid (8 MSps), timing offset in samples of the
    last iteration)."""
    n = int(min(n, (len(y) - s0) // SPS - 1))
    if n < 48 or s0 < 0:
        return None
    pad = 2 * SPS
    r = np.convolve(np.r_[np.zeros(pad, complex), y, np.zeros(pad + len(C0), complex)], C0[::-1], "full")
    r = r[len(C0) - 1 + pad:]                    # r[m] = sum_j y[m + j] C0[j]
    p = diff_products(y)
    # bit timing: the one-bit differential swing |Im p| peaks with the window on the bit (40 / 40 exact at 6 dB;
    # the access-address search is +-3 samples and the first-pass phase table is timing-sensitive)
    k = np.arange(n) * SPS
    offs = [o for o in range(-SPS // 2, SPS // 2 + 1) if s0 + o >= 0 and s0 + o + k[-1] < len(p)]
    s0 += max(offs, key=lambda o: np.mean(np.abs(np.imag(p[s0 + o + k]))))
    a_in = None
    for _ in range(iters):
        alpha, rel, h, dt, ph = _coherent(r, p, s0, n, a_in)
        a_in = alpha
        if dt and s0 + dt >= 0 and (len(y) - s0 - dt) // SPS - 1 >= n:
            s0 += dt                                # the next first pass on the refined bit timing
            dt = 0
    mid = n // 2
    w2 = (ph[min(mid + 1, n - 1)] - ph[mid])           # rad / bit of s = 2 x carrier
    return (alpha > 0).astype(np.uint8), rel, h, w2 / 2 / (2 * np.pi) * 1e6, dt


# ---- CRC-24 (BLE): register as an int, bit p = reg[p] of detector._ble_crc_ok
CRC_MASK = sum(1 << q for q in (1, 3, 4, 6, 9, 10))
_NMAX = 16 + 8 * 255 + 24


def _clock0(S):
    fb = (S >> 23) & 1
    S = ((S << 1) & 0xFFFFFF) | fb
    return S ^ CRC_MASK if fb else S


def _tables():
    M = np.zeros(_NMAX, np.int64)                       # contribution of a lone 1 that has j zero bits after it
    S = 1 ^ CRC_MASK
    for j in range(_NMAX):
        M[j] = S; S = _clock0(S)
    I = np.zeros(_NMAX + 1, np.int64)                   # init 0x555555 clocked over N zero bits
    S = 0
    init = 0x555555
    for p in range(24):
        if (init >> p) & 1:
            S |= 1 << p
    for N in range(_NMAX + 1):
        I[N] = S; S = _clock0(S)
    return M, I


_M, _I = _tables()


def crc_syndrome_parts(dw, n_pdu):
    """(syndrome, per-position contributions) of dewhitened bits dw: syndrome 0 <=> CRC ok. A flip at position i
    (< n_pdu + 24) changes the syndrome by contrib[i]."""
    i = np.arange(n_pdu)
    contrib = np.r_[_M[n_pdu - 1 - i], [1 << (23 - t) for t in range(24)]].astype(np.int64)
    bits = dw[:n_pdu + 24].astype(bool)
    syn = np.bitwise_xor.reduce(contrib[bits]) if bits.any() else np.int64(0)
    return int(syn ^ _I[n_pdu]), contrib


def chase(dw, n_pdu, rel, off, k_least=16, max_units=3, protect=()):
    """Flip <= max_units b-symbol units (each flips bits k and k+1) among the k_least least reliable b_k of the
    PDU + CRC (bit index = off + position) to zero the CRC syndrome. Units touching `protect` positions (the
    length field) are excluded. Returns (corrected dewhitened bits, n units) or None."""
    syn, contrib = crc_syndrome_parts(dw, n_pdu)
    if syn == 0:
        return dw, 0
    npos = n_pdu + 24
    ks = np.arange(off, off + npos)                 # b index k: flips positions k-off and k-off+1 (bits k, k+1)
    r = rel[ks].copy()
    r[-1] = -1.0                                    # the last b has no successor (ISI model, TX ramp-down): always tried
    order = ks[np.argsort(r, kind="stable")]
    units = []
    for k in order:
        pos = [q for q in (k - off, k - off + 1) if 0 <= q < npos]
        if any(q in protect for q in pos):
            continue
        units.append((pos, int(np.bitwise_xor.reduce(contrib[pos]))))
        if len(units) >= k_least:
            break
    for m in range(1, max_units + 1):
        for comb in combinations(range(len(units)), m):
            v = 0
            for u in comb:
                v ^= units[u][1]
            if v == syn:
                out = dw.copy()
                for u in comb:
                    for q in units[u][0]:
                        out[q] ^= 1
                return out, m
    return None
