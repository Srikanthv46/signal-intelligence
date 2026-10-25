"""Coherent GFSK (h 0.5) detection via the Laurent C0 pulse: matched filter, data-aided CFO / phase from known bits,
decision-directed PLL; alpha_k = Im(b_k b*_{k-1}). Bench only (PLAN L1)."""
import numpy as np
from scipy.special import erf

SPS = 8


def c0_pulse(bt=0.5, sps=SPS, L=3):
    """Laurent's main pulse C0 for h = 0.5 GFSK (Gaussian frequency pulse truncated to L bits)."""
    T = 1.0
    t = np.arange(0, (L + 1) * sps) / sps          # 0 .. (L+1)T
    s = np.sqrt(np.log(2)) / (2 * np.pi * bt)
    # frequency pulse g(t) centred at L/2, phase pulse q(t) = int g, q(LT) = 1/2
    tt = np.linspace(0, L, 4000)
    g = 0.5 * (erf((tt - L / 2 + 0.5) / (np.sqrt(2) * s)) - erf((tt - L / 2 - 0.5) / (np.sqrt(2) * s)))
    q = np.cumsum(g) * (tt[1] - tt[0]); q *= 0.5 / q[-1]
    qf = lambda x: np.interp(x, tt, q, left=0.0, right=0.5)
    h = 0.5
    def S(x):
        x = np.abs(x)
        out = np.zeros_like(x)
        m1 = x < L * T
        out[m1] = np.sin(2 * np.pi * h * qf(x[m1])) / np.sin(np.pi * h)
        m2 = (x >= L * T) & (x < 2 * L * T)
        out[m2] = np.sin(np.pi * h - 2 * np.pi * h * qf(x[m2] - L * T)) / np.sin(np.pi * h)
        return out
    c = np.ones_like(t)
    for i in range(L):
        c *= S(t + i * T)
    return c


C0 = c0_pulse()
C0_PEAK = int(np.argmax(C0))


def laurent_symbols(bits, b0=1 + 0j):
    """b_k = b_{k-1} * j * alpha_k, alpha = 2 bit - 1."""
    a = 2 * np.asarray(bits, float) - 1
    return b0 * np.cumprod(1j * a)


def mf(y):
    return np.convolve(y, C0[::-1], "full")[len(C0) - 1:]     # r[n] = sum y[n+m] C0[m]


def detect_coherent(y, known_bits, start, n_bits, loop_bw=0.01, cfo_known=None):
    """y at SPS samples/bit, start = sample index of bit 0 of known_bits (+ the pulse delay handled here).
    Returns (bits_hat, reliab) for all n_bits from bit 0."""
    pad = 2 * SPS
    r = mf(np.r_[np.zeros(pad, complex), y])
    idx = pad + start + np.arange(n_bits) * SPS     # start: MF sample of bit 0 (its frequency pulse begins 1 bit early)
    idx = idx[idx < len(r)]
    rk = r[idx]
    nk = len(known_bits)
    bk = laurent_symbols(known_bits)
    z = rk[1:nk] * np.conj(bk[1:nk])                # data-removed
    w = np.angle(np.sum(z[1:] * np.conj(z[:-1])))   # rad / bit
    if cfo_known is not None:
        w = cfo_known
    rk = rk * np.exp(-1j * w * np.arange(len(rk)))
    th = np.angle(np.sum(rk[1:nk] * np.conj(bk[1:nk])))
    # DD loop with ISI removal: v_k ~ A (b_k P0 + (b_{k-1} + b_{k+1}) P1) e^{j theta}; b_{k+1} tentative
    ac = np.correlate(C0, C0, "full"); P0 = ac[len(C0) - 1]; P1 = ac[len(C0) - 1 + SPS] / P0
    A = np.abs(np.sum(rk[1:nk] * np.conj(bk[1:nk]))) / (nk - 1)
    n = len(rk)
    bh = np.zeros(n, complex); rel = np.zeros(n)
    k1 = 2 * np.pi * loop_bw; k2 = k1 ** 2 / 4; nu = 0.0

    def dec(v, k):
        ax = 1j ** (k + 1)                           # the axis of b_k (b_0 = j alpha_0)
        c = np.real(v * np.conj(ax))
        return (1.0 if c >= 0 else -1.0) * ax, abs(c)

    for k in range(n):
        v = rk[k] * np.exp(-1j * th)
        d, rel[k] = dec(v, k)
        if k < nk:
            d = bk[k]
        bh[k] = d
        nx = dec(rk[k + 1] * np.exp(-1j * (th + nu)), k + 1)[0] if k + 1 < n else 0
        if k + 1 < nk:
            nx = bk[k + 1]
        pv = bh[k - 1] if k else 0
        ex = A * (d + P1 * (pv + nx))
        e = np.imag(v * np.conj(ex)) / (abs(ex) ** 2 + 1e-12)
        nu += k2 * e; th += k1 * e + nu
    alpha = np.imag(bh[1:] * np.conj(bh[:-1]))
    bits = np.r_[known_bits[0], (alpha > 0).astype(np.uint8)]
    relb = np.r_[np.inf, np.minimum(rel[1:], rel[:-1])]
    return bits, relb, rel


def detect_nda(y, start, n_bits, W=16, fine=True):
    """Coherent h 0.5 GFSK detection without known bits or a loop: b_k^2 = (-1)^(k+1), so s_k = r_k^2 (-1)^(k+1)
    is a tone at twice the CFO with phase 2 theta. CFO from its FFT over the whole burst, phase from a +-W
    moving sum (unwrapped), decisions on the axis of each b_k, alpha_k = Im(b_k b*_(k-1))."""
    pad = 2 * SPS
    r = mf(np.r_[np.zeros(pad, complex), y, np.zeros(pad, complex)])
    best = None
    for dt in (range(-3, 4) if fine else (0,)):
        idx = pad + start + dt + np.arange(n_bits) * SPS
        idx = idx[idx < len(r)]
        rk = r[idx]
        s = rk ** 2 * (-1.0) ** (np.arange(len(rk)) + 1)
        nf = 1 << int(np.ceil(np.log2(len(s) * 8)))
        F = np.abs(np.fft.fft(s, nf))
        j = int(np.argmax(F))
        a, b, c = F[j - 1], F[j], F[(j + 1) % nf]
        dj = 0.5 * (a - c) / (a - 2 * b + c + 1e-30)
        w2 = 2 * np.pi * (((j + dj) / nf + 0.5) % 1 - 0.5)     # rad / bit of s = 2 w
        if best is None or b > best[0]:
            best = (b, dt, rk, w2)
    _, dt, rk, w2 = best
    s = rk ** 2 * (-1.0) ** (np.arange(len(rk)) + 1) * np.exp(-1j * w2 * np.arange(len(rk)))
    cs = np.r_[0, np.cumsum(s)]
    n = len(s)
    lo = np.clip(np.arange(n) - W, 0, n); hi = np.clip(np.arange(n) + W + 1, 0, n)
    ph2 = np.unwrap(np.angle(cs[hi] - cs[lo])) + w2 * np.arange(n)
    v = rk * np.exp(-0.5j * ph2)
    ax = 1j ** (np.arange(n) + 1)
    c = np.real(v * np.conj(ax))
    bh = np.sign(c + 1e-30) * ax
    alpha = np.imag(bh[1:] * np.conj(bh[:-1]))
    return (alpha > 0).astype(np.uint8), np.minimum(np.abs(c[1:]), np.abs(c[:-1])), np.abs(c), dt


def _lpf_taps(ofs=8e6, fc=0.6e6):
    from scipy import signal
    return signal.firwin(33, fc, fs=ofs)


H_LPF = _lpf_taps()


def _pattern_table():
    """Per-bit phase change (rad) of clean h 0.5 BT 0.5 GFSK for alpha_k = +1, by neighbours (alpha_{k-1}, alpha_{k+1}),
    measured exactly as the first pass measures it (LPF, one-bit differential at the bit centre)."""
    from spectrum_intel.receivers.rc_link import synth as S
    rng = np.random.default_rng(0)
    bits = rng.integers(0, 2, 4000).astype(np.uint8)
    iq, _ = S.gfsk_burst(1e6, 250e3, bits=bits, bt=0.5, fs=8e6, rng=rng)
    z = np.convolve(iq, H_LPF, "same")
    d = np.angle(z[SPS:] * np.conj(z[:-SPS]))[np.arange(4000 - 1) * SPS + SPS // 2 - SPS // 2]
    a = 2 * bits.astype(int) - 1
    tab = np.zeros((2, 2))
    for i in (0, 1):
        for j in (0, 1):
            m = np.zeros(len(a), bool); m[1:-2] = True
            k = np.flatnonzero(m)
            k = k[((a[k - 1] == a[k]) == bool(i)) & ((a[k + 1] == a[k]) == bool(j))]   # i, j: neighbour = this bit
            tab[i, j] = np.mean(d[k] * a[k])
    return tab


TAB = _pattern_table()


def first_pass(y, start, n_bits, nseg=8, hs=np.arange(0.42, 0.581, 0.02), a_in=None):
    """Differential decisions; then for each h: q_k = p_k exp(-j (h/0.5) alpha_k TAB_k) (modulation removed), complex
    sums over nseg segments, weighted quadratic fit of their unwrapped phases -> frequency track (rad / bit);
    h = the most coherent. Complex sums: no FM clicks (an angle fit broke at 8 dB)."""
    z = np.convolve(np.r_[y, np.zeros(SPS, complex)], H_LPF, "same")
    p = z[SPS:] * np.conj(z[:-SPS])
    idx = start + 9 + np.arange(n_bits) * SPS
    idx = idx[(idx >= 0) & (idx < len(p))]
    pk = p[idx]
    d = np.angle(pk)
    n = len(d)
    a = np.where(d > np.mean(_km(d)), 1, -1)
    if a_in is not None:
        a = np.r_[a_in, a[len(a_in):]][:n]
    an = np.r_[a[0], a, a[-1]]
    t = TAB[(an[:-2] == a).astype(int), (an[2:] == a).astype(int)] * a
    edges = np.linspace(0, n, nseg + 1).astype(int)
    kc = (edges[:-1] + edges[1:]) / 2 - n / 2
    best = None
    for h in hs:
        q = pk * np.exp(-1j * (h / 0.5) * t)
        S_ = np.array([q[edges[i]:edges[i + 1]].sum() for i in range(nseg)])
        sc = np.abs(S_).sum()
        if best is None or sc > best[0]:
            best = (sc, h, S_)
    _, h, S_ = best
    ph = np.unwrap(np.angle(S_))
    wgt = np.abs(S_)
    coef = np.polyfit(kc, ph, 2, w=wgt) if nseg >= 4 else np.polyfit(kc, ph, 0)
    ftrack = np.polyval(coef, np.arange(n) - n / 2)
    return a, h, ftrack


def _km(v, it=20):
    a, b = np.percentile(v, 25), np.percentile(v, 75)
    for _ in range(it):
        m = np.abs(v - a) < np.abs(v - b)
        if m.all() or (~m).all():
            break
        a, b = v[m].mean(), v[~m].mean()
    return min(a, b), max(a, b)


def detect_coh2(y, start, n_bits, W=16, use_ft=True, use_walk=True):
    a, h, ft = first_pass(y, start, n_bits)
    if not use_ft:
        ft = ft * 0
    if not use_walk:
        h = 0.5
    # remove the frequency track (sample level) and the h walk (symbol level)
    ph_bits = np.r_[0, np.cumsum(ft)]               # phase at bit boundaries
    tt = (np.arange(len(y)) - (start + 9)) / SPS
    ph = np.interp(tt, np.arange(len(ph_bits)), ph_bits, left=None, right=None)
    ph = np.where(tt < 0, ft[0] * tt, ph)
    ph = np.where(tt > len(ft), ph_bits[-1] + ft[-1] * (tt - len(ft)), ph)
    y2 = y * np.exp(-1j * ph)
    pad = 2 * SPS
    r = mf(np.r_[np.zeros(pad, complex), y2, np.zeros(pad, complex)])
    best = None
    walk = np.pi * (h - 0.5) * np.r_[np.cumsum(a)]
    for dt in range(-3, 4):
        idx = pad + start + dt + np.arange(n_bits) * SPS
        idx = idx[idx < len(r)]
        rk = r[idx] * np.exp(-1j * walk[:len(idx)])
        s = rk ** 2 * (-1.0) ** (np.arange(len(rk)) + 1)
        nf = 1 << int(np.ceil(np.log2(len(s) * 8)))
        F = np.abs(np.fft.fft(s, nf))
        j = int(np.argmax(F))
        aa, bb, cc = F[j - 1], F[j], F[(j + 1) % nf]
        dj = 0.5 * (aa - cc) / (aa - 2 * bb + cc + 1e-30)
        w2 = 2 * np.pi * (((j + dj) / nf + 0.5) % 1 - 0.5)
        if best is None or bb > best[0]:
            best = (bb, dt, rk, w2)
    _, dt, rk, w2 = best
    n = len(rk)
    s = rk ** 2 * (-1.0) ** (np.arange(n) + 1) * np.exp(-1j * w2 * np.arange(n))
    cs = np.r_[0, np.cumsum(s)]
    lo = np.clip(np.arange(n) - W, 0, n); hi = np.clip(np.arange(n) + W + 1, 0, n)
    ph2 = np.unwrap(np.angle(cs[hi] - cs[lo])) + w2 * np.arange(n)
    v = rk * np.exp(-0.5j * ph2)
    ax = 1j ** (np.arange(n) + 1)
    c = np.real(v * np.conj(ax))
    bh = np.sign(c + 1e-30) * ax
    alpha = np.imag(bh[1:] * np.conj(bh[:-1]))
    return (alpha > 0).astype(np.uint8), np.minimum(np.abs(c[1:]), np.abs(c[:-1])), h, dt, a


def _tone_freq(s, pad=8):
    nf = 1 << int(np.ceil(np.log2(max(len(s), 16) * pad)))
    F = np.abs(np.fft.fft(s, nf))
    j = int(np.argmax(F))
    a, b, c = F[j - 1], F[j], F[(j + 1) % nf]
    dj = 0.5 * (a - c) / (a - 2 * b + c + 1e-30)
    return 2 * np.pi * (((j + dj) / nf + 0.5) % 1 - 0.5), b


def track_phase(s, stages=(8, 4, 1)):
    """Phase track of a tone with drift: per stage, FFT-peak frequency of each segment (on the residual), polynomial
    fit of frequency vs time (deg <= 2), integrated into the phase."""
    n = len(s)
    ph = np.zeros(n)
    k = np.arange(n)
    for ns in stages:
        res = s * np.exp(-1j * ph)
        e = np.linspace(0, n, ns + 1).astype(int)
        fr = []; wt = []; kc = []
        for i in range(ns):
            seg = res[e[i]:e[i + 1]]
            if len(seg) < 8:
                continue
            w, m = _tone_freq(seg)
            fr.append(w); wt.append(m); kc.append((e[i] + e[i + 1] - 1) / 2)
        if not fr:
            break
        deg = min(2, len(fr) - 1)
        cf = np.polyfit(kc, fr, deg, w=np.sqrt(wt)) if deg > 0 else np.array([fr[0]])
        fk = np.polyval(cf, k)
        ph = ph + np.cumsum(fk) - fk[0]
        ph -= 0
    # constant phase handled by the windowed estimate
    return ph


def detect_coh3(y, start, n_bits, W=16, iters=2):
    a_in = None
    for it in range(iters):
        out = _coh3(y, start, n_bits, W, a_in)
        a_in = np.r_[1 if it == 0 else a_in[0], 2 * out[0].astype(int) - 1]   # alpha_0 is not decided coherently
        a_in[0] = out[4][0]
    return out[:4]


def _coh3(y, start, n_bits, W, a_in):
    a, h, ft = first_pass(y, start, n_bits, a_in=a_in)
    pad = 2 * SPS
    r = mf(np.r_[np.zeros(pad, complex), y, np.zeros(pad, complex)])
    walk = np.pi * (h - 0.5) * np.cumsum(a)
    ac = np.correlate(C0, C0, "full"); P1 = ac[len(C0) - 1 + SPS] / ac[len(C0) - 1]
    an = np.r_[a[1:], a[-1]]
    isi = 1 + 1j * P1 * (an - a)                    # r_k = b_k (P0 + j P1 (alpha_{k+1} - alpha_k)): divide it out
    best = None
    for dt in range(-3, 4):
        idx = pad + start + dt + np.arange(n_bits) * SPS
        idx = idx[idx < len(r)]
        rk = r[idx] * np.exp(-1j * walk[:len(idx)])   # decisions: raw (the ISI is on the other axis)
        s = (rk / isi[:len(idx)]) ** 2 * (-1.0) ** (np.arange(len(rk)) + 1)   # phase: ISI removed
        ph = track_phase(s, (8, 4))
        _, m = _tone_freq(s * np.exp(-1j * ph))
        if best is None or m > best[0]:
            best = (m, dt, rk, s, ph)
    _, dt, rk, s, ph = best
    ph = ph + track_phase(s * np.exp(-1j * ph), (1,))
    n = len(rk)
    s2 = s * np.exp(-1j * ph)
    cs = np.r_[0, np.cumsum(s2)]
    lo = np.clip(np.arange(n) - W, 0, n); hi = np.clip(np.arange(n) + W + 1, 0, n)
    ph2 = np.unwrap(np.angle(cs[hi] - cs[lo])) + ph
    v = rk * np.exp(-0.5j * ph2)
    ax = 1j ** (np.arange(n) + 1)
    c = np.real(v * np.conj(ax))
    bh = np.sign(c + 1e-30) * ax
    alpha = np.imag(bh[1:] * np.conj(bh[:-1]))
    return (alpha > 0).astype(np.uint8), np.minimum(np.abs(c[1:]), np.abs(c[:-1])), h, dt, a
