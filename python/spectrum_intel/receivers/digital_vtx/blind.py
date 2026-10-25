"""
blind.py -- vendor-independent CP-OFDM link detector (digital_vtx "generic").

A drone video / data link of any make (DJI OcuSync / O4 / Lightbridge, Herelink, Autel, Walksnail,
an unknown one) is, at the physical layer, a CP-OFDM signal: every symbol ends with a copy of its
first Tu seconds. That repetition is found without knowing the vendor, the numerology or the band,
and well below the noise floor (the whole capture is summed coherently at lag Tu):

1. scan: each slice (10 / 20 MHz, half-overlapping) is time-AGC'd (a 20 us envelope, so strong Wi-Fi
   bursts carry no more weight than the quiet stretches between them) and its autocorrelation over
   lags 8-120 us is computed by FFT; peaks above the detrended lag profile are candidate Tu.
2. confirm: a candidate is an OFDM link when its CP metric has a symbol-rate line (cyclic spectrum,
   4 ms blocks) or, after the period refinement, a sharp per-window fold in >= 3 windows of 1 ms
   (high-SNR LTE-like links have many slot-rate lines that hide the symbol line). Tu is refined to
   a fraction of a sample (the lag must match within ~0.1 sample: a one-sample error decorrelates
   the band edges of an 18 MHz link) and every later slice is taken at a rate where Tu is an integer
   number of samples.
3. band: the CP cross-spectrum -- the FFT of each symbol's CP and of the same window Tu later,
   multiplied and summed coherently over symbols -- holds only what repeats at lag Tu: overlapping
   Wi-Fi, an RC hop at another timing or a mirrored IQ image 17 dB down does not widen the band.
4. bursts: the symbol grid is found per 2 ms window (1 ms hop) after refining the symbol period from
   the window-to-window phase drift, so a link that re-times mid-capture (D_20260609-151042: the
   grid jumps by 10.3 us at 16.8 ms) keeps all its bursts. Per-symbol CP evidence, runs, then the
   TDD frame (activity ACF + fold contrast) completes weak bursts frame by frame.
5. gates: Tu >= 10 us (SCS <= 100 kHz: 802.11-family numerologies are the Wi-Fi receivers' and an
   FM video line rate folds onto a short "Tu"), >= 40 symbols in bursts, one burst >= 10 symbols;
   duplicates (same Tu, bands overlapping > 50 %) and IQ images (same Tu, mirrored about fc, >= 10 dB
   weaker) are dropped.

6. hops (blind_hops): narrowband CP-OFDM bursts that hop -- a drone's RC uplink (DJI O4: 30 kHz,
   ~4.3 MHz, 0.41 ms, once per 2.5 ms frame, TDD with the video; OcuSync 2 / Autel: 15 kHz,
   1.1 MHz). Too short and too spread in frequency for the capture-level scan, they are searched with
   a windowed CP fold (core.ofdm15k_scan) at every numerology found in the capture (a drone's uplink
   shares its video numerology) and at the drone families' (15 / 30 kHz). Each hit gets its edges from
   the in-band envelope, its band from the on-vs-off spectrum (rectangle fit), and is kept only when
   the CP fold at that numerology over the refined burst is significant; hops the scan missed are
   re-tested at the frame cadence.

blind_links(sl, fs, fc, usable_half, noise_pps) -> list of link dicts (offsets relative to fc).
blind_hops(sl, fs, fc, usable_half, links) -> list of hop-link dicts.
sl(f_off, half_bw, out_fs) -> (complex slice, its exact rate) is the caller's channelizer.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage
from scipy.signal import find_peaks

from spectrum_intel.core import bursts as BU, ofdm15k as O15

WIDTHS = (10e6, 20e6)                # scan slice widths (step W/2)
AGC_S = 20e-6                        # time-AGC envelope
TU_SCAN = (8e-6, 120e-6)             # lag range of the scan
TU_MIN = 10e-6                       # shortest reported Tu (SCS <= 100 kHz)
K_ACF = 7.0                          # ACF peak threshold (robust z over the detrended lag profile)
LINE_Z = 10.0                        # symbol-rate line threshold (noise / Wi-Fi peaks <= 5)
WIN_Z, WIN_N = 8.0, 3                # or: >= WIN_N windows of 1 ms whose fold z > WIN_Z
MAX_CLUSTERS = 24
CP_RATIOS = (1 / 32, 1 / 16, 1 / 14, 1 / 12, 1 / 8, 1 / 6, 1 / 4)
SEG_WIN_S = 2e-3                     # grid-phase window (hop half of it)
SEG_ZW = 4.0                         # a window holds the link when its fold z > this
MIN_SYMS, MIN_LONGEST = 40, 10


def _movsum(v, k):
    c = np.concatenate(([0], np.cumsum(v)))
    return c[k:] - c[:-k]


def _cp_metric(y, N, w):
    p = y[:-N] * np.conj(y[N:])
    e = (np.abs(y[:-N]) ** 2 + np.abs(y[N:]) ** 2) / 2
    cs = _movsum(p, w)
    return np.abs(cs) / (_movsum(e, w) + 1e-20), cs


def _agc(y, ofs):
    env = ndimage.uniform_filter1d(np.abs(y) ** 2, max(4, int(AGC_S * ofs)))
    return y / np.sqrt(env + 1e-30)


def _acf(y, Lmax):
    """Energy-normalised autocorrelation sum_n y[n+L] conj(y[n]) for L = 0..Lmax."""
    n = len(y)
    nf = 1 << int(np.ceil(np.log2(n + Lmax + 1)))
    r = np.fft.ifft(np.abs(np.fft.fft(y, nf)) ** 2)[:Lmax + 1]
    cs = np.concatenate(([0], np.cumsum(np.abs(y) ** 2)))
    L = np.arange(Lmax + 1)
    en = 0.5 * ((cs[n - L] - cs[0]) + (cs[n] - cs[L]))
    return np.conj(r) / (en + 1e-30)


def rate_for(tu, B, cp=0.0):
    """Output rate >= 1.25 B (>= 2.5 MHz; >= 16 samples per CP when cp is given: the per-symbol CP
    metric of a narrow link with a short CP is too noisy on a handful of samples) at which Tu is an
    integer number of samples."""
    k = int(np.ceil(tu * max(2.5e6, 1.25 * B, 16 / cp if cp > 0 else 0.0)))
    return k / tu, k


def scan(sl, usable, skip=None):
    """Lag-Tu autocorrelation peaks of every W-wide slice. skip(off): a slice whose whole band (sampled every 1 MHz)
    is skipped is not scanned (PLAN G8: a capture inside a cellular band costs no scan)."""
    hits = []
    for W in [w for w in WIDTHS if w <= 2 * usable] or [2 * usable]:     # a narrow capture: one slice
        for off in np.arange(-usable + W / 2, usable - W / 2 + 1, W / 2):
            if skip is not None and all(skip(float(f)) for f in np.arange(off - W / 2, off + W / 2 + 1, 1e6)):
                continue
            y, ofs = sl(float(off), W / 2, 1.25 * W)
            y = _agc(y, ofs)
            Lmin, Lmax = int(TU_SCAN[0] * ofs), int(TU_SCAN[1] * ofs)
            if Lmax >= len(y) // 2:
                continue
            a = np.abs(_acf(y, Lmax))[Lmin:]
            d = a - ndimage.median_filter(a, max(7, int(1.5e-6 * ofs)) | 1)
            sd = 1.4826 * np.median(np.abs(d - np.median(d))) + 1e-12
            pk, _ = find_peaks(d, height=K_ACF * sd, distance=max(3, int(0.5e-6 * ofs)))
            for p in pk:
                hits.append(dict(W=W, off=float(off), N=int(p + Lmin), ofs=ofs, tu=(p + Lmin) / ofs,
                                 z=float(d[p] / sd)))
    return hits


def _symrate(m, N, ofs, seg_s=4e-3, lo=1.02, hi=1.30):
    """Symbol period (samples) from the strongest cyclic line of the CP metric in
    [1/(hi N), 1/(lo N)], power-summed over seg_s blocks; z against the median over a wider range."""
    L = min(int(seg_s * ofs), len(m))
    nf = 1 << int(np.ceil(np.log2(L)))
    w = np.hanning(L)
    acc = np.zeros(nf // 2 + 1)
    for i in range(0, len(m) - L + 1, L):
        s = m[i:i + L] - np.mean(m[i:i + L])
        acc += np.abs(np.fft.rfft(s * w, nf)) ** 2
    a = np.arange(len(acc)) / nf
    win = (a >= 1 / (hi * N)) & (a <= 1 / (lo * N))
    wide = (a >= 1 / (2.5 * N)) & (a <= 1 / (0.8 * N))
    i = int(np.argmax(np.where(win, acc, -1.0)))
    z = acc[i] / (np.median(acc[wide]) + 1e-30)
    d = 0.0
    if 0 < i < len(acc) - 1:
        l, c, r = acc[i - 1:i + 2]
        den = l - 2 * c + r
        d = 0.5 * (l - r) / den if den != 0 else 0.0
    return nf / (i + d), float(z)


def _tu_refine(y, ofs, N):
    """Sub-sample Tu: parabola through |ACF| around its maximum in N-1..N+1 (the lag with the best
    symbol-rate line need not be the ACF peak)."""
    r = np.array([abs(np.vdot(y[:-L], y[L:])) for L in range(N - 2, N + 3)])
    j = 1 + int(np.argmax(r[1:4]))
    den = r[j - 1] - 2 * r[j] + r[j + 1]
    d = 0.5 * (r[j - 1] - r[j + 1]) / den if den < 0 else 0.0
    return (N - 2 + j + float(np.clip(d, -0.5, 0.5))) / ofs


def window_phases(m, P, Wn, H):
    """Per window (length Wn, hop H): the CP metric folded modulo P (3-bin smoothed) -> phase of
    its peak and the peak's robust z."""
    nW = max(1, (len(m) - Wn) // H + 1)
    nb = int(np.ceil(P))
    ph = np.zeros(nW, np.int64)
    zw = np.zeros(nW)
    for w in range(nW):
        seg = m[w * H:w * H + Wn]
        idx = np.floor(np.mod(np.arange(w * H, w * H + len(seg)), P)).astype(np.int64)
        f = np.bincount(idx, seg, nb) / np.maximum(np.bincount(idx, minlength=nb), 1)
        f = f + np.roll(f, 1) + np.roll(f, -1)
        med = np.median(f)
        k = int(np.argmax(f))
        ph[w] = k
        zw[w] = (f[k] - med) / (1.4826 * np.median(np.abs(f - med)) + 1e-12)
    return ph, zw


def refine_period(m, P, ofs, win_s=1e-3, iters=4):
    """Symbol period from the drift of the per-window grid phase: the median step between
    consecutive strong windows (a timing jump is an outlier)."""
    Wn = int(win_s * ofs)
    for _ in range(iters):
        ph, zw = window_phases(m, P, Wn, Wn)
        d = [float(np.mod(ph[i + 1] - ph[i] + P / 2, P) - P / 2) for i in range(len(ph) - 1)
             if zw[i] > 6 and zw[i + 1] > 6]
        d = [v for v in d if abs(v) < P / 4]
        if len(d) < 3:
            break
        st = float(np.median(d))
        if abs(st) < 0.5:
            break
        P = P + st / (Wn / P)
    return P


def confirm(sl, off, W, N0, ofs0):
    """-> (tu_s, ts_s, line_z, strong_windows) or None."""
    yr, ofs0 = sl(off, W / 2, ofs0)
    y = _agc(yr, ofs0)
    best = None
    for NN in (N0 - 1, N0, N0 + 1):
        m, _ = _cp_metric(y, NN, max(3, NN // 16))
        P, z = _symrate(m, NN, ofs0)
        if best is None or z > best[1]:
            best = (P, z, NN)
    P0, zl, NN = best
    tu = _tu_refine(yr, ofs0, NN)
    ofs, N = rate_for(tu, W)
    yr, ofs = sl(off, W / 2, ofs)
    m, _ = _cp_metric(yr, N, max(3, N // 16))
    Wn = int(1e-3 * ofs)
    nwin, Pw = 0, None
    for Pc in [P0 / ofs0 * ofs] + [N * (1 + r) for r in CP_RATIOS]:
        if not 1.02 * N <= Pc <= 1.30 * N:
            continue
        Pr = refine_period(m, Pc, ofs)
        k = int(np.sum(window_phases(m, Pr, Wn, Wn)[1] > WIN_Z))
        if k > nwin:
            nwin, Pw = k, Pr
    if zl >= LINE_Z:
        return tu, P0 / ofs0, zl, nwin
    if nwin >= WIN_N:
        return tu, Pw / ofs, zl, nwin
    return None


def segment(sl, off, B, tu, ts):
    """Bursts on a local symbol grid. Returns (runs [(t0_s, t1_s, n_sym)], symbol times (s),
    per-symbol z, refined symbol period (s), number of grids, CP-phase carrier fraction (Hz))."""
    ofs, N = rate_for(tu, B, ts - tu)
    y, ofs = sl(off, B / 2, ofs)
    P = ts * ofs
    cp = max(3, int(round(P - N)))
    m, cs = _cp_metric(y, N, max(3, int(0.8 * cp)))
    P = refine_period(m, P, ofs)
    Wn = int(SEG_WIN_S * ofs)
    H = Wn // 2
    ph, zw = window_phases(m, P, Wn, H)
    act = np.where(zw > SEG_ZW)[0]
    if act.size == 0:
        return [], np.zeros(0), np.zeros(0), P / ofs, 0, 0.0
    # grids: consecutive strong windows whose phases differ by more than a CP start a new one
    sp = ph[zw > 6].astype(float)
    ngrid = 1 + int(np.sum(np.abs(np.mod(np.diff(sp) + P / 2, P) - P / 2) > cp)) if sp.size > 1 else 1
    smp = lambda a: m[max(int(round(a)) - 2, 0):int(round(a)) + 3].max()
    tk, on, off_, q = [], [], [], []
    for bk in range(len(m) // H):
        cand = [w for w in (bk - 1, bk) if 0 <= w < len(zw) and zw[w] > SEG_ZW]
        if not cand:
            j = int(act[np.argmin(np.abs(act - bk))])
            if abs(j - bk) > 2:
                continue
            cand = [j]
        w = max(cand, key=lambda i: (zw[i], -i))
        p0 = ph[w]
        k = int(np.ceil((bk * H - p0) / P))
        while p0 + k * P < (bk + 1) * H and p0 + k * P + P + 3 < len(m):
            a = p0 + k * P
            tk.append(a); on.append(smp(a)); off_.append(smp(a + P / 2)); q.append(smp(a + P / 4))
            k += 1
    tk, on, off_, q = map(np.asarray, (tk, on, off_, q))
    nref = off_ - q
    zr = (on - off_) / (1.4826 * np.median(np.abs(nref - np.median(nref))) + 1e-9)
    mask = ndimage.binary_closing(ndimage.uniform_filter1d(zr, 4) * 2 > 4.0, np.ones(3))
    lab, _ = ndimage.label(mask)
    runs = []
    for sl_ in ndimage.find_objects(lab):
        a, b = sl_[0].start, sl_[0].stop
        while a < b and zr[a] < 1.0: a += 1
        while b > a and zr[b - 1] < 1.0: b -= 1
        if b - a >= 3 and zr[a:b].mean() * np.sqrt(b - a) >= 6:
            runs.append((tk[a] / ofs, (tk[b - 1] + P) / ofs, b - a))
    mr = []
    for r in runs:
        if mr and r[0] - mr[-1][1] <= 3.5 * P / ofs:
            mr[-1] = (mr[-1][0], r[1], mr[-1][2] + r[2])
        else:
            mr.append(r)
    good = [int(t) for t, z in zip(tk, zr) if z > 3 and int(t) < len(cs)]
    frac = float(-np.angle(np.sum(cs[good])) / (2 * np.pi) * (ofs / N)) if good else 0.0
    return mr, tk / ofs, zr, P / ofs, ngrid, frac


def _fold_contrast(tk, e, T, ts):
    nb = max(8, int(round(T / ts)))
    bi = np.minimum((np.mod(tk, T) / T * nb).astype(np.int64), nb - 1)
    M = np.bincount(bi, e, nb) / np.maximum(np.bincount(bi, minlength=nb), 1)
    zt = M * np.sqrt(max(tk[-1] / T, 1.0))
    return float(np.mean(np.sort(zt)[-max(nb // 4, 1):])), zt


def frame_period(tk, zr, ts, lo=0.9e-3, hi=12e-3):
    """TDD frame period: ACF of the clipped per-symbol evidence (smoothed: bursts are long, so its
    peak is broad), smallest lag within 0.7 of the best, then a fine search (+-6 %) on the
    frame-fold contrast. -> (T_s, contrast) or (None, 0)."""
    if len(tk) < 50:
        return None, 0.0
    e = np.clip(zr, -3, 6)
    n = int(np.ceil(tk[-1] / ts)) + 1
    a = np.zeros(n)
    a[np.round(tk / ts).astype(np.int64)] = e
    a = a - a.mean()
    A = np.fft.irfft(np.abs(np.fft.rfft(a, 2 * n)) ** 2, 2 * n)[:n]
    A = ndimage.uniform_filter1d(A / (A[0] + 1e-30), max(3, int(0.2e-3 / ts)))
    L = np.arange(n)
    ok = (L * ts >= lo) & (L * ts <= hi) & (L < n // 3)
    pk = [int(p) for p in find_peaks(np.where(ok, A, -1.0))[0] if ok[p] and A[p] > 0.1]
    if not pk:
        return None, 0.0
    best = max(A[p] for p in pk)
    T0 = min(p for p in pk if A[p] >= 0.7 * best) * ts
    Ts = np.arange(0.94 * T0, 1.06 * T0, ts / 4)
    # fine search on the variance the fold explains (R^2, fixed bin count): the fold contrast scales with
    # sqrt(frames in the block), so it always grew towards the search's short edge (Walksnail-like 5805 MHz:
    # 6.673 ms reported, 7.163 ms from the data-run starts; PLAN D7)
    nb = max(8, int(round(T0 / ts)))
    ec = e - e.mean()
    tot = float(np.sum(ec ** 2)) + 1e-30
    r2 = []
    for T in Ts:
        bi = np.minimum((np.mod(tk, T) / T * nb).astype(np.int64), nb - 1)
        cnt = np.bincount(bi, minlength=nb)
        sm = np.bincount(bi, ec, nb)
        r2.append(float(np.sum(sm ** 2 / np.maximum(cnt, 1))) / tot)
    k = int(np.argmax(r2))
    return float(Ts[k]), float(_fold_contrast(tk, e, Ts[k], ts)[0])


def frame_complete(runs, tk, zr, ts, T, min_con=6.0):
    """Weak bursts from the TDD frame (as dji._frame_complete, on local grids): the folded
    per-symbol evidence gives the burst regions of the frame; in every frame whose region holds
    evidence (mean z sqrt(n) >= 4) its strongest contiguous stretch (gaps <= 4 filled) is added,
    edges trimmed to symbols on by their own evidence; dips <= 6 symbols inside a region joined."""
    e = np.clip(zr, -3, 6)
    con, zt = _fold_contrast(tk, e, T, ts)
    if con < min_con:
        return runs
    nb = len(zt)
    sh = int(np.argmin(ndimage.uniform_filter1d(zt, 5, mode="wrap")))
    lab, _ = ndimage.label(np.roll(zt, -sh) > 3)
    regs = [((s[0].start + sh) * T / nb, (s[0].stop + sh) * T / nb)
            for s in ndimage.find_objects(lab) if s[0].stop - s[0].start >= 3]
    add, present = [], []
    for f in range(-1, int(tk[-1] / T) + 1):
        for a, b in regs:
            sel = np.where((tk >= f * T + a) & (tk < f * T + b))[0]
            if sel.size < 3:
                continue
            seg = e[sel]
            if seg.mean() * np.sqrt(seg.size) < 4:
                continue
            present.append((tk[sel[0]], tk[sel[-1]] + ts))
            on = ndimage.uniform_filter1d(seg, 4) >= 1
            idx = np.where(on)[0]
            for i0, i1 in zip(idx[:-1], idx[1:]):
                if 1 < i1 - i0 <= 5:
                    on[i0:i1] = True
            lab2, _ = ndimage.label(on)
            best = None
            for s in ndimage.find_objects(lab2):
                a2, b2 = s[0].start, s[0].stop
                if b2 - a2 >= 3:
                    key = (float(seg[a2:b2].sum()), a2, b2)
                    if best is None or key > best:
                        best = key
            if best is None:
                continue
            _, a2, b2 = best
            while a2 < b2 and seg[a2] < 1: a2 += 1
            while b2 > a2 and seg[b2 - 1] < 1: b2 -= 1
            if b2 - a2 >= 3:
                add.append((tk[sel[a2]], tk[sel[b2 - 1]] + ts, b2 - a2))
    out = []
    for r in sorted(runs + add):
        if out and (r[0] <= out[-1][1] + 1.5 * ts or (r[0] - out[-1][1] <= 6.5 * ts and
                    any(p0 <= out[-1][1] and r[0] <= p1 for p0, p1 in present))):
            out[-1] = (out[-1][0], max(out[-1][1], r[1]), out[-1][2] + r[2])
        else:
            out.append(r)
    return out


def band_cp(sl, off, W, tu, ts, tsym):
    """CP cross-spectrum over [off-W, off+W]: per symbol, FFT of a CP-long window at the CP and of
    the same window Tu later; Y2 conj(Y1) summed coherently (the CFO phase is common to all
    symbols). Returns (frequency offsets, real part rotated onto the mean phase, per symbol)."""
    ofs, N = rate_for(tu, 2 * W)
    y, ofs = sl(off, W, ofs)
    cp = int(round((ts - tu) * ofs))
    M = max(8, int(0.9 * cp))
    nf = max(64, 1 << int(np.ceil(np.log2(4 * M))))
    win = np.hanning(M)
    acc = np.zeros(nf, complex)
    n = 0
    for t in tsym:
        a = int(round(t * ofs)) + max(0, (cp - M) // 2)
        if a < 0 or a + N + M > len(y):
            continue
        acc += np.fft.fft(y[a + N:a + N + M] * win, nf) * np.conj(np.fft.fft(y[a:a + M] * win, nf))
        n += 1
    acc = np.fft.fftshift(acc)
    c = np.real(acc * np.exp(-1j * np.angle(np.sum(acc))))
    return off + (np.arange(nf) - nf / 2) * ofs / nf, c / max(n, 1)


def measure_band(sl, off, W, tu, ts, tsym, usable):
    """Occupied band: the CP cross-spectrum within 6 dB of its in-band level, gaps < 1 MHz closed,
    the component holding the most repeating power; retried on a 2x / 4x wider slice, centred on
    the estimate, while it touches the slice edge. -> (lo, hi) offsets or None."""
    r = None
    for hw in (W, 2 * W, 4 * W):
        hw = min(hw, usable)
        if r is not None:
            off = float(np.clip((r[0] + r[1]) / 2, -usable + hw, usable - hw)) if hw < usable else 0.0
        fr, c = band_cp(sl, off, hw, tu, ts, tsym)
        sm = ndimage.uniform_filter1d(c, 3)
        top = sm[sm > 0.5 * np.percentile(sm, 95)]
        if top.size == 0 or np.median(top) <= 0:
            return None
        df = fr[1] - fr[0]
        occ = ndimage.binary_closing(sm > 0.25 * np.median(top), np.ones(max(3, int(1e6 / df))))
        lab, nl = ndimage.label(occ)
        if nl == 0:
            return None
        bi = max(range(1, nl + 1), key=lambda i: (float(sm[lab == i].sum()), -i))
        idx = np.where(lab == bi)[0]
        r = (float(fr[idx[0]]), float(fr[idx[-1]]))
        if (r[0] - (off - hw) > 0.04 * 2 * hw and (off + hw) - r[1] > 0.04 * 2 * hw) or hw >= usable:
            break
    return r


def blind_links(sl, fs, fc, usable, noise_pps, log=None, tf=None, keep_slice_snr=None, skip=None):
    """All CP-OFDM links (Tu >= 10 us) in the capture; see the module docstring. tf: the capture's
    core.bursts.TFPower (band and burst extent from in-band power, seeded by the CP runs). keep_slice_snr: links
    at or above this SNR keep their slice (centre, band, rate_for(Tu, B)) as lk['slice'] = (y, rate). skip(off): True
    leaves a candidate at that offset unconfirmed (cheap band gate)."""
    dur = None
    hits = sorted(scan(sl, usable, skip), key=lambda h: (-h["z"], h["W"], h["off"], h["N"]))
    cl = []
    for h in hits:
        for c in cl:
            if abs(h["tu"] - c["tu"]) < 0.015 * c["tu"] + 0.1e-6 and abs(h["off"] - c["off"]) < min(h["W"], c["W"]) / 2:
                break
        else:
            cl.append(h)
    out = []
    for c in cl[:MAX_CLUSTERS]:
        if c["tu"] < TU_MIN or (skip is not None and skip(c["off"])):
            continue
        if any(abs(c["tu"] - o["tu"]) < 0.015 * o["tu"] and o["lo"] <= c["off"] <= o["hi"] for o in out):
            continue
        W = c["W"]
        cf = confirm(sl, c["off"], W, c["N"], c["ofs"])
        if cf is None:
            if log:
                log(f"  Tu {c['tu'] * 1e6:.2f} us (ACF z {c['z']:.1f}) at {(fc + c['off']) / 1e6:.1f} MHz: not confirmed")
            continue
        tu, ts, zl, nwin = cf
        runs, tk, zr, ts, ng, frac = segment(sl, c["off"], W, tu, ts)
        if not runs:
            continue
        bd = measure_band(sl, c["off"], W, tu, ts, tk[zr > 2], usable)
        if bd is None:
            continue
        lo, hi = bd
        B, cen = hi - lo, (lo + hi) / 2
        runs, tk, zr, ts, ng, frac = segment(sl, cen, 1.02 * B, tu, ts)
        T, con = frame_period(tk, zr, ts) if len(tk) else (None, 0.0)
        if T:
            runs = frame_complete(runs, tk, zr, ts, T)
        nsym = [int(round((r[1] - r[0]) / ts)) for r in runs]
        if not runs or sum(nsym) < MIN_SYMS or max(nsym) < MIN_LONGEST:
            if log:
                log(f"  Tu {tu * 1e6:.2f} us at {(fc + cen) / 1e6:.2f} MHz: too little OFDM "
                    f"({sum(nsym)} symbols, longest {max(nsym) if nsym else 0})")
            continue
        einfo = {}
        if tf is not None:
            # band and burst extent from in-band power, seeded by the CP runs (core.bursts): the CP
            # cross-spectrum band of D_20260609-151042 ch0 reached 2.6 MHz into Wi-Fi, and the per-symbol
            # CP evidence splits a weak burst
            seeds = [(a, b) for a, b, _ in runs]
            for it in range(2):
                fb = BU.fit_band(tf, seeds, lo, hi, max(lo - 0.5 * B, -usable), min(hi + 0.5 * B, usable))
                if fb is not None:
                    lo, hi = fb["lo"], fb["hi"]
                    B, cen = hi - lo, (lo + hi) / 2
                if it == 0:
                    seeds, einfo = BU.energy_bursts(tf, seeds, lo, hi, usable, min_sym_s=ts)
                    if not einfo.get("used"):
                        break
            if einfo.get("used"):
                runs = [(a, b, int(round((b - a) / ts))) for a, b in seeds]
        runs = [(float(a), float(b), n) for a, b, n in runs]
        if any(abs(o["tu"] - tu) < 0.015 * tu and min(o["hi"], hi) - max(o["lo"], lo) > 0.5 * min(o["hi"] - o["lo"], B)
               for o in out):
            continue
        ofs2, _ = rate_for(tu, B)
        yb, ofs2 = sl(cen, B / 2, ofs2)
        inb = np.zeros(len(yb), bool)
        for a, b, _ in runs:
            inb[int(a * ofs2):int(b * ofs2)] = True
        p_on = float(np.mean(np.abs(yb[inb]) ** 2)) if inb.any() else 0.0
        snr = 10 * np.log10(max(p_on / (noise_pps * B / fs + 1e-30) - 1, 1e-3))
        dur = len(yb) / ofs2
        lk = dict(tu=tu, ts=ts, lo=lo, hi=hi, cen=cen, bw=B, line_z=zl, win_z=nwin, acf_z=c["z"],
                  runs=runs, n_grids=ng, frame_s=T, frame_contrast=con, snr_db=float(snr), p_on=p_on,
                  duty=sum(r[1] - r[0] for r in runs) / dur, cfo_frac_hz=frac)
        if keep_slice_snr is not None and snr >= keep_slice_snr:
            lk["slice"] = (yb, ofs2)
        if log:
            log(f"  link Tu {tu * 1e6:.3f} us Ts {ts * 1e6:.3f} us at {(fc + cen) / 1e6:.3f} MHz, {B / 1e6:.2f} MHz, "
                f"line {zl:.1f}, windows {nwin}, {len(runs)} bursts, duty {lk['duty']:.2f}, frame "
                f"{None if T is None else round(T * 1e3, 3)} ms, {ng} grid(s), SNR {snr:.1f} dB, energy {einfo}")
        out.append(lk)
    # IQ image: a mirrored copy about fc (same Tu), >= 10 dB weaker
    for o in out:
        o["image_of"] = next(((fc + q["cen"]) / 1e6 for q in out if q is not o and abs(q["tu"] - o["tu"]) < 0.015 * o["tu"]
                              and abs(q["cen"] + o["cen"]) < 0.5e6 and q["p_on"] > 10 * o["p_on"]), None)
    return [o for o in out if o["image_of"] is None]


# ============================================================== hopping narrowband links (RC uplinks)
HOP_FAMILIES = ((66.667e-6, 144 / 2048), (33.333e-6, 72 / 1024))   # LTE-like 15 kHz, DJI O4 30 kHz (normal CP)
HOP_WIN_S = 0.25e-3                  # scan window (half-overlapping)
HOP_K = 5.0                          # scan threshold (robust sigma of the window scores)
HOP_TRACK_K = 3.5                    # threshold at the times the frame cadence predicts
HOP_Z = 6.0                          # CP fold z of a verified hop (fixed numerology; noise <= ~4.5)
HOP_CON_DB = 3.0                     # on-vs-off band contrast of a hop
HOP_MIN_SYMS = 4
HOP_MAX_S = 1.0e-3                   # longer: not a hop (OFDM RC uplink hops: O4 0.41 ms, OcuSync 2 0.51 ms;
                                     # video bursts >= 1.2 ms)
HOP_BW = (0.5e6, 8e6)                # hop width (video modes: >= 9 MHz)
HOP_PAD_S = 1e-3                     # context either side of a hit for its edges and off-spectrum
HOP_DUR_TOL, HOP_BW_TOL = 0.2, 0.3   # a link's hops: dwell / width within these fractions of the median


def hop_width(tu):
    """Scan slice width for a numerology: ~75 subcarriers, >= 1 MHz (OcuSync / O4 uplinks: 72+ used)."""
    return max(1.0e6, 75.0 / tu)


def _fold_z(y, N, cp, P):
    """Robust z of the sharpest bin of the CP metric folded modulo the (fractional) period P."""
    m, _ = _cp_metric(y, N, cp)
    K = int(len(m) // P)
    if K < 3:
        return 0.0
    idx = (np.arange(K)[:, None] * P + np.arange(int(P))[None, :]).astype(np.int64)
    f = m[idx].mean(0)
    med = np.median(f)
    mad = 1.4826 * np.median(np.abs(f - med)) + 1e-12
    return float((f.max() - med) / mad)


def refine_hop(sl, f_off, t0_s, t1_s, tu, cpf, W):
    """Edges, band and CP evidence of one hop candidate. Returns dict(t0, t1, lo, hi, z, contrast)
    or None. Edges: the in-band envelope (10 us blocks) above the geometric mean of its 90th / 10th
    percentiles within +-1 ms; band: on-burst vs off-burst median spectrum (pads either side) over
    3 W, rectangle fit; evidence: CP fold z at (tu, cpf) over the burst on a slice of that band."""
    ts = tu * (1 + cpf)
    yn, o2 = sl(f_off, W / 2, 1.25 * W)
    a2 = max(0, int((t0_s - HOP_PAD_S) * o2)); b2 = min(len(yn), int((t1_s + HOP_PAD_S) * o2))
    k = max(4, int(10e-6 * o2))
    nb = (b2 - a2) // k
    if nb < 5:
        return None
    e = (np.abs(yn[a2:a2 + nb * k]) ** 2).reshape(nb, k).mean(1)
    on = e > np.sqrt(np.percentile(e, 90) * np.percentile(e, 10))
    on = ndimage.binary_closing(on, np.ones(3, bool))
    pk = int(np.clip(((t0_s + t1_s) / 2 * o2 - a2) // k, 0, nb - 1))
    if not on[pk]:
        w = np.where(on)[0]
        if not w.size:
            return None
        pk = int(w[np.argmin(np.abs(w - pk))])
    s0 = s1 = pk
    while s0 > 0 and on[s0 - 1]: s0 -= 1
    while s1 < nb - 1 and on[s1 + 1]: s1 += 1
    t0, t1 = (a2 + s0 * k) / o2, (a2 + (s1 + 1) * k) / o2
    if t1 - t0 < HOP_MIN_SYMS * ts or t1 - t0 > HOP_MAX_S:
        return None
    # band: on-burst vs off-burst spectrum over a window of +-1.5 W; a band that reaches the window edge
    # is re-centred on itself (window widened to its width) once, and rejected if it still does: the
    # signal is wider than a hop (a slice of a video burst)
    c, hw = f_off, 1.5 * W
    for _ in range(2):
        y, ofs = sl(c, hw, 2.5 * hw)
        nf = 1 << int(np.ceil(np.log2(ofs * 32e-6)))
        win = np.hanning(nf)

        def frames(u0, u1):
            i0, i1 = max(int(u0 * ofs), 0), min(int(u1 * ofs), len(y))
            n = (i1 - i0) // nf
            if n < 1:
                return None
            return np.abs(np.fft.fftshift(np.fft.fft(y[i0:i0 + n * nf].reshape(n, nf) * win, axis=1), axes=1)) ** 2
        fon = frames(t0, t1)
        foff = [f for f in (frames(t0 - HOP_PAD_S, t0 - nf / ofs), frames(t1 + nf / ofs, t1 + HOP_PAD_S))
                if f is not None]
        if fon is None or not foff:
            return None
        fb = c + (np.arange(nf) - nf // 2) * ofs / nf
        inw = np.where(np.abs(fb - c) <= hw)[0]                # the slice is brick-wall: bins inside it
        S_on = ndimage.uniform_filter1d(np.median(fon, 0)[inw], 3)
        S_off = ndimage.uniform_filter1d(np.median(np.concatenate(foff), 0)[inw], 3)
        d = 10 * np.log10(np.maximum(S_on, 1e-30) / np.maximum(S_off, 1e-30))
        a, b, inn, out = BU.rect_fit(d)
        lo, hi = float(fb[inw[a]] - ofs / nf / 2), float(fb[inw[b - 1]] + ofs / nf / 2)
        edge = a <= 1 or b >= len(inw) - 1
        if not edge:
            break
        c, hw = (lo + hi) / 2, max(hw, hi - lo)
    if edge:
        return None
    bw = hi - lo
    if not (HOP_BW[0] <= bw <= HOP_BW[1]):
        return None
    o3, N = rate_for(tu, bw, tu * cpf)
    cp = max(1, int(round(N * cpf)))
    y3, o3 = sl((lo + hi) / 2, bw / 2, o3)
    z = _fold_z(y3[int(t0 * o3):int(t1 * o3)], N, cp, N * (1 + cpf))
    return dict(t0=t0, t1=t1, lo=lo, hi=hi, z=z, contrast=inn - out)


def hop_period(starts, lo=1e-3, hi=20e-3, tol=30e-6):
    """Hop cadence: the longest period (5 us grid) that the most start-time differences are whole
    multiples of (within tol); None below 3 such pairs (two pairs fit every common divisor). A link
    that re-times mid-capture keeps its period (pairs across the jump simply do not count)."""
    s = np.sort(np.asarray(starts, float))
    d = np.array([b - a for i, a in enumerate(s) for b in s[i + 1:]])
    if d.size == 0:
        return None
    best, Tb = 0, None
    for T in np.arange(lo, hi + 1e-12, 5e-6):
        q = np.round(d / T)
        n = int(np.sum((q >= 1) & (np.abs(d - q * T) < tol)))
        if n >= best:
            best, Tb = n, float(T)
    return Tb if best >= 3 else None


class _Ch:
    """core.ofdm15k's channelizer interface (slice -> samples) over a slicer."""

    def __init__(self, sl):
        self.sl = sl

    def slice(self, f, h, r):
        return self.sl(f, h, r)[0]


def _overlap(a0, a1, b0, b1):
    return min(a1, b1) - max(a0, b0)


def blind_hops(sl, fs, fc, usable, links, log=None, skip=None):
    """Hopping narrowband CP-OFDM links (see the module docstring, 6). links: blind_links' output (their
    numerologies are searched first, and a hit inside one of their bursts is that link's symbols). skip(off): True
    leaves a scan hit at that offset unrefined (cheap band gate, as blind_links); when the whole usable band is
    skipped there is nothing to scan (the scan's threshold is a whole-grid statistic, so a partly skipped band is
    still scanned in full)."""
    if skip is not None and all(skip(float(f)) for f in np.arange(-usable, usable + 1, 1e6)):
        return []
    nums = []
    for tu, cpf, T in [(lk["tu"], (lk["ts"] - lk["tu"]) / lk["tu"], lk["frame_s"]) for lk in links] + \
                      [(t, c, None) for t, c in HOP_FAMILIES]:
        if tu >= TU_MIN and not any(abs(tu - n[0]) < 0.015 * n[0] for n in nums):
            nums.append((float(tu), float(cpf), T))
    ch = _Ch(sl)
    dur = None
    out, found = [], []
    for tu, cpf, T in nums:
        W = hop_width(tu)
        ofs, N = rate_for(tu, W, tu * cpf)
        cp = max(1, int(round(N * cpf)))
        ts = tu * (1 + cpf)
        hits, st = O15.ofdm15k_scan(ch, fc, [("hop", ofs, N, cp, W, W / 2)], win_s=HOP_WIN_S, k_sigma=HOP_K,
                                    usable_half=usable)
        own = [lk for lk in links if abs(lk["tu"] - tu) < 0.015 * tu]
        if dur is None:
            dur = len(ch.slice(0.0, W / 2, 1.25 * W)) / (1.25 * W)

        def inside_link(r):
            c = (r["lo"] + r["hi"]) / 2
            return any(lk["lo"] <= c <= lk["hi"] and any(_overlap(r["t0"], r["t1"], a, b) > 0.5 * (r["t1"] - r["t0"])
                                                         for a, b, _ in lk["runs"]) for lk in own)

        def dup(r, hops):
            return any(_overlap(r["t0"], r["t1"], h["t0"], h["t1"]) > 0.5 * min(r["t1"] - r["t0"], h["t1"] - h["t0"])
                       and _overlap(r["lo"], r["hi"], h["lo"], h["hi"]) > 0 for h in hops)

        hops = []
        for h in sorted(hits, key=lambda h: (-h["score"], h["t0_s"], h["f_off"])):
            if any(h["t0_s"] < g["t1"] and h["t1_s"] > g["t0"] and g["lo"] <= h["f_off"] <= g["hi"] for g in hops):
                continue
            if skip is not None and skip(h["f_off"]):
                continue
            r = refine_hop(sl, h["f_off"], h["t0_s"], h["t1_s"], tu, cpf, W)
            if r is None or r["z"] < HOP_Z or r["contrast"] < HOP_CON_DB or inside_link(r) or dup(r, hops):
                continue
            r["tracked"] = False
            hops.append(r)
        # hops the scan missed: re-test the times the frame cadence predicts, at a lower threshold
        if len(hops) >= 2:
            Tc = hop_period([h["t0"] for h in hops]) or T
            if Tc:
                dmed = float(np.median([h["t1"] - h["t0"] for h in hops]))
                bmed = float(np.median([h["hi"] - h["lo"] for h in hops]))
                cand = sorted({round(h["t0"] + k * Tc, 6) for h in hops for k in range(-int(dur / Tc) - 1, int(dur / Tc) + 2)
                               if 0 <= h["t0"] + k * Tc <= dur - dmed})
                times = []
                for t in cand:
                    if any(_overlap(t, t + dmed, g["t0"], g["t1"]) > 0 for g in hops):
                        continue
                    if not times or t - times[-1] > 0.1e-3:
                        times.append(t)
                thr = st["hop"]["median"] + HOP_TRACK_K * st["hop"]["sigma"]
                Wn = int(dmed * ofs)
                offs = np.arange(-usable + W / 2, usable - W / 2 + 1, W / 2)
                slices = [(float(o), ch.slice(float(o), W / 2 * 1.05, ofs)) for o in offs]
                for t in times:
                    sc = []
                    for o, y in slices:
                        best = (-1.0, 0.0)
                        for dt in np.arange(-0.15e-3, 0.151e-3, 0.05e-3):
                            a = int((t + dt) * ofs)
                            if a < 0 or a + Wn + N + cp >= len(y):
                                continue
                            m, _ = _cp_metric(y[a:a + Wn + N + cp], N, cp)
                            f = O15._folded_scores(m, N + cp, Wn, Wn)
                            if f.size and f[0] > best[0]:
                                best = (float(f[0]), t + dt)
                        if best[0] > thr:
                            sc.append((best[0], o, best[1]))
                    for _, o, t0 in sorted(sc, reverse=True)[:4]:
                        r = refine_hop(sl, o, t0, t0 + dmed, tu, cpf, W)
                        if (r is None or r["z"] < HOP_Z or r["contrast"] < HOP_CON_DB or inside_link(r) or dup(r, hops)
                                or abs((r["t1"] - r["t0"]) - dmed) > 0.25 * dmed
                                or abs((r["hi"] - r["lo"]) - bmed) > 0.35 * bmed):
                            continue
                        r["tracked"] = True
                        hops.append(r)
                        break
        if log:
            log(f"  hops Tu {tu * 1e6:.2f} us (W {W / 1e6:.2f} MHz): {len(hits)} scan hits, {len(hops)} verified "
                + str([(round(h['t0'] * 1e3, 2), round((fc + (h['lo'] + h['hi']) / 2) / 1e6, 2), round(h['z'], 1))
                       for h in sorted(hops, key=lambda h: h['t0'])]))
        for h in hops:
            h["num"] = len(found)
        found.append(dict(tu=tu, cpf=cpf, T=T, ts=ts, own=own, hops=hops))
    # one hop, one numerology: where hops of two numerologies coincide (O4's comb-4 pilots also fold at
    # twice its Tu), the stronger CP fold wins
    kept = []
    for h in sorted((h for f in found for h in f["hops"]), key=lambda h: (-h["z"], h["t0"], h["lo"])):
        if not any(_overlap(h["t0"], h["t1"], g["t0"], g["t1"]) > 0.5 * min(h["t1"] - h["t0"], g["t1"] - g["t0"])
                   and _overlap(h["lo"], h["hi"], g["lo"], g["hi"]) > 0 for g in kept):
            kept.append(h)
    for i, f in enumerate(found):
        tu, cpf, T, ts, own = f["tu"], f["cpf"], f["T"], f["ts"], f["own"]
        hops = sorted((h for h in kept if h["num"] == i), key=lambda h: h["t0"])
        # the hops of one link share their dwell and width (O4: 0.41 ms x 4.4 MHz, OcuSync 2: 0.51 ms x
        # 1.1 MHz): keep those within HOP_DUR_TOL / HOP_BW_TOL of the median (twice, so that one outlier
        # does not set the reference)
        for _ in range(2):
            if len(hops) < 2:
                break
            dm = float(np.median([h["t1"] - h["t0"] for h in hops]))
            bm = float(np.median([h["hi"] - h["lo"] for h in hops]))
            hops = [h for h in hops if abs(h["t1"] - h["t0"] - dm) <= HOP_DUR_TOL * dm
                    and abs(h["hi"] - h["lo"] - bm) <= HOP_BW_TOL * bm]
        if len(hops) < 2:
            continue                          # one verified burst is not a hopping link
        # hop channels: hops whose bands overlap by more than half
        chans = []
        for h in sorted(hops, key=lambda h: (h["lo"] + h["hi"]) / 2):
            for c in chans:
                if _overlap(h["lo"], h["hi"], c["lo"], c["hi"]) > 0.5 * min(h["hi"] - h["lo"], c["hi"] - c["lo"]):
                    c["hops"].append(h)
                    break
            else:
                chans.append(dict(lo=h["lo"], hi=h["hi"], hops=[h]))
        for c in chans:
            cen = float(np.median([(h["lo"] + h["hi"]) / 2 for h in c["hops"]]))
            for h in c["hops"]:
                h["channel"] = cen
        Tc = hop_period([h["t0"] for h in hops]) or T
        out.append(dict(tu=tu, ts=ts, cpf=cpf, hops=hops, n_channels=len(chans), frame_s=Tc,
                        paired=[lk["cen"] for lk in own],
                        dur_s=float(np.median([h["t1"] - h["t0"] for h in hops])),
                        bw=float(np.median([h["hi"] - h["lo"] for h in hops]))))
    return out
