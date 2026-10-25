"""Bench: DJI O4 (OcuSync 4) link measurements for emitter identity.

O4 at 61.44 MSps (measured 2026-10-16, FINDINGS.md):
- 30 kHz SCS, N = 2048; a 0.5 ms slot (30720 samples) has 13 symbols with a 144-sample CP plus one symbol with a
  176-sample CP (NR-like), so the symbol grid is found by folding the CP metric on the slot, not on one symbol;
- reference symbols at the odd slot positions, pilots on every second subcarrier (comb-2), the same symbol set
  recurring on the 2.5 ms frame (5 slots); 10 MHz mode carries a 300-subcarrier slice of the 20 MHz sequence.

Measurements per link (one contiguous stretch of samples, no recorder discontinuity inside):
- slot grid and fine CFO (CP phase, +-15 kHz);
- comb-2 differential pilot spectra X(k+2) X*(k): channel-, timing- and power-independent, used to align a link's
  subcarrier grid on another link's (integer shift) and to test that two links carry the same sequence;
- symbol-clock offset from the timing drift between recurring reference symbols 5 slots apart;
- carrier = f0 + CFO + 30 kHz x n, n the integer bin offset of the transmitter's centre.
"""
from __future__ import annotations

from collections import Counter

import numpy as np
from scipy.signal import find_peaks
from scipy.stats import binom

FS = 61.44e6
N = 2048
CP = 144
SLOT = 30720
SCS = FS / N                       # 30 kHz
FRAME_SLOTS = 5                    # 2.5 ms


def load_cf32(path, count=-1, offset=0):
    return np.fromfile(path, np.complex64, count=count, offset=8 * offset)


def baseband(x, fc, f0, bw):
    """Shift f0 to DC and keep |f| <= bw/2 (FFT brick wall)."""
    X = np.fft.fft(x); f = np.fft.fftfreq(len(x), 1 / FS)
    X[np.abs(f - (f0 - fc)) > bw / 2] = 0
    return np.fft.ifft(X) * np.exp(-2j * np.pi * (f0 - fc) * np.arange(len(x)) / FS)


def cp_metric(y):
    p = y[:-N] * np.conj(y[N:]); c = np.cumsum(np.concatenate([[0], p])); P = c[CP:] - c[:-CP]
    e = np.abs(y[:-N]) ** 2 + np.abs(y[N:]) ** 2; ce = np.cumsum(np.concatenate([[0], e])); E = (ce[CP:] - ce[:-CP]) / 2
    return P, E


def slot_grid(y, a, b, on_thr=0.35):
    """Symbol starts in y[a:b] from the CP metric folded on the slot; on-flags; fine CFO (Hz, +-SCS/2)."""
    P, E = cp_metric(y[a:b + N + CP])
    ph = (a + np.arange(len(P))) % SLOT
    m = np.bincount(ph, weights=np.abs(P), minlength=SLOT) / (np.bincount(ph, weights=E, minlength=SLOT) + 1e-30)
    pk, _ = find_peaks(np.concatenate([m, m[:3000]]), distance=1800)
    pk = np.unique(pk % SLOT)
    pk = np.sort(pk[np.argsort(-m[pk])][:14])
    starts = np.array([sl + q for sl in range(a - a % SLOT, b, SLOT) for q in pk if a <= sl + q < b - N - 200])
    rho = np.array([abs(P[s - a]) / (E[s - a] + 1e-30) for s in starts])
    on = rho > on_thr
    z = np.sum(P[starts[on] - a]) if on.any() else 0j
    return dict(starts=starts, on=on, rho=rho, cfo=float(-np.angle(z) * FS / (2 * np.pi * N)),
                quality=float(m[pk].min() / np.median(m)))


def symbol_spectra(y, starts, cfo, back=8):
    """FFT of every symbol (window starts `back` samples inside the CP), CFO removed, centred (bin N/2 = f0)."""
    yc = y * np.exp(-2j * np.pi * cfo * np.arange(len(y)) / FS)
    out = np.zeros((len(starts), N), complex); ok = np.zeros(len(starts), bool)
    for i, s in enumerate(starts):
        j = int(s) + CP - back
        if j + N <= len(y):
            out[i] = np.fft.fftshift(np.fft.fft(yc[j:j + N])); ok[i] = True
    return out, ok


def diff_rows(X, sel, K, SH):
    """Unit-modulus comb-2 differential spectra over bins -K-SH..K+SH of the selected symbols."""
    v = X[sel][:, N // 2 - K - SH - 2: N // 2 + K + SH + 2]
    d = v[:, 2:] * np.conj(v[:, :-2])
    return d / (np.abs(d) + 1e-30)


def match(A, KA, B, KB, SH, nfft=4096):
    """For every row of A: best normalised |corr| over the rows of B and integer shifts -SH..SH (B bin = A bin + s)."""
    K = min(KA, KB)
    a = A[:, SH + KA - K: SH + KA + K]; b = B[:, KB - K: KB + K + 2 * SH]
    FA = np.fft.fft(a, nfft, axis=1); FB = np.fft.fft(b, nfft, axis=1)
    best = np.zeros(len(a)); shift = np.zeros(len(a), int)
    for i in range(len(a)):
        c = np.abs(np.fft.ifft(FB * np.conj(FA[i]), axis=1)[:, :2 * SH + 1]) / (2 * K)
        j, s = np.unravel_index(np.argmax(c), c.shape)
        best[i] = c[j, s]; shift[i] = s - SH
    return best, shift


def match_stats(best, shift, SH):
    """Most common shift among the strongest quarter of A symbols, and its binomial p under 'no common sequence'."""
    top = np.argsort(-best)[:max(10, len(best) // 4)]
    s, c = Counter(shift[top].tolist()).most_common(1)[0]
    return dict(shift=int(s), hits=int(c), n=len(top), p=float(binom.sf(c - 1, len(top), 1 / (2 * SH + 1))),
                best_median=float(np.median(best)), best_max=float(best.max()))


def sco_ppm(X, starts, sel, K, min_corr=0.12, max_tau=4.0, max_slots=40):
    """Symbol-clock offset of the transmitter relative to our sample clock, in ppm (+ = transmitter faster).

    Reference symbols recur at whole-slot lags (2 slots, the 5-slot frame, ...): X_j(k) = X_i(k) exp(-j 2 pi k tau
    / N), tau (samples) the delay of the second against the nominal grid. tau from the phase slope across the
    pilot bins; every recurring pair (|corr| >= min_corr, any lag up to max_slots) gives tau / dt = -(eps_tx -
    eps_rx); weighted LS through the origin over all pairs."""
    idx = {int(s): i for i, s in enumerate(starts)}
    taus, wts, dts = [], [], []
    for i, s in enumerate(starts):
      for m in range(1, max_slots + 1):
        j = idx.get(int(s) + m * SLOT)
        if j is None or not (sel[i] and sel[j]):
            continue
        a = X[i][N // 2 - K: N // 2 + K]; b = X[j][N // 2 - K: N // 2 + K]
        cor = abs(np.vdot(a, b)) / np.sqrt(np.vdot(a, a).real * np.vdot(b, b).real)
        if cor < min_corr:
            continue
        D = a * np.conj(b)
        k = np.arange(-K, K)
        # pilot parity: the comb-2 bins whose step-2 products are coherent
        zq = [np.sum(D[q + 2::2][:len(D[q::2]) - 1] * np.conj(D[q::2][:len(D[q::2]) - 1])) for q in (0, 1)]
        q = int(np.argmax(np.abs(zq)))
        tau_c = np.angle(zq[q]) * N / (4 * np.pi)                    # coarse: phase step over 2 bins
        # fine: weighted LS of the residual phase across the whole band (slope = 2 pi tau / N per bin)
        Dp = D[q::2] * np.exp(-2j * np.pi * k[q::2] * tau_c / N); kp = k[q::2]
        Dp = Dp * np.exp(-1j * np.angle(np.sum(Dp)))
        ph = np.angle(Dp); w = np.abs(Dp)
        sl = np.sum(w * kp * ph) / np.sum(w * kp ** 2)
        tau = tau_c + sl * N / (2 * np.pi)
        if abs(tau) < max_tau:
            taus.append(tau); wts.append(cor ** 2); dts.append(m * SLOT)
    if len(taus) < 3:
        return dict(ppm=None, sigma=None, n=len(taus))
    taus = np.array(taus); w = np.array(wts); dt = np.array(dts, float)
    keep = np.ones(len(taus), bool)
    for _ in range(3):                                   # LS through the origin, 5-MAD outlier rejection
        r_ = taus[keep] / dt[keep]
        d = np.sum(w[keep] * dt[keep] * taus[keep]) / np.sum(w[keep] * dt[keep] ** 2)
        res = taus - d * dt
        mad = np.median(np.abs(res[keep])) + 1e-4
        keep = np.abs(res) < 5 * mad
    d = np.sum(w[keep] * dt[keep] * taus[keep]) / np.sum(w[keep] * dt[keep] ** 2)
    res = taus[keep] - d * dt[keep]
    # conservative sigma: residual scatter, pairs sharing symbols counted as n_symbols / 2 independent ones
    n_eff = max(len(set(np.flatnonzero(keep).tolist())) // 4, 2)
    sd = np.sqrt(np.sum(w[keep] * res ** 2) / np.sum(w[keep])) / np.sqrt(np.mean(dt[keep] ** 2)) / np.sqrt(n_eff)
    return dict(ppm=float(-d * 1e6), sigma=float(sd * 1e6), n=int(keep.sum()))


def analyse_link(x, fc, f0, bw, a, b, SH=200):
    """Everything measured on one link in samples [a, b) of capture x (centre fc)."""
    y = baseband(x, fc, f0, bw + 2e6)
    g = slot_grid(y, a, b)
    X, ok = symbol_spectra(y, g["starts"], g["cfo"])
    sel = g["on"] & ok
    K = int(bw / SCS / 2) - 8
    return dict(f0=f0, bw=bw, K=K, grid_quality=g["quality"], n_on=int(sel.sum()), cfo=g["cfo"],
                rows=diff_rows(X, sel, K, SH), sco=sco_ppm(X, g["starts"], sel, K))


def carrier_ppm(f0, cfo, n, f_nom):
    car = f0 + cfo + SCS * n
    return car, (car - f_nom) / f_nom * 1e6


def resolve_carrier(f0, cfo, sco, raster=0.5e6, nmax=4):
    """Integer bin n and nominal centre such that the carrier ppm equals the symbol-clock ppm (one crystal):
    for each n, f_nom = carrier / (1 + sco); keep the n whose f_nom lies closest to the raster."""
    best = None
    for n in range(-nmax, nmax + 1):
        car = f0 + cfo + SCS * n
        fn = car / (1 + sco * 1e-6)
        r = fn - raster * round(fn / raster)
        if best is None or abs(r) < abs(best[2]):
            best = (n, raster * round(fn / raster), r)
    return dict(n=best[0], f_nom=best[1], raster_residual_hz=float(best[2]))
