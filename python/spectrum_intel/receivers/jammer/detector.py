"""Swept-carrier jammer detector.

A barrage jammer of the common kind is one carrier (a VCO) swept back and forth across a band by a periodic
tuning voltage. In a short STFT every frequency bin it crosses lights up once or twice per sweep period, so
the per-bin occupancy has a strong cyclic line at the sweep rate, and the STFT folded on that period shows
one carrier whose frequency moves continuously with phase.

    detect(x, fs, fc) -> list[Emission]      family "jammer", label "swept_carrier" (at most one per block)
    measure(x, fs, fc) -> (result | None, diagnostics)

1. Occupancy: Hann STFT power, fs/64 bins, hop 32, summed over frame pairs (0.96 MHz x 1.04 us at 61.44 MSps);
   the outer 1/16 of the band each side (anti-alias roll-off) is dropped.
2. Trigger: cyclic spectrum of the per-bin log occupancy, whitened per bin and summed over bins; up to N_CAND
   lines >= TRIG_DB over their local median, between max(300 Hz, 8 periods per block) and row_rate / 8, each
   refined on its harmonics.
3. Period: a multiple d <= MAX_DIV of the line's period. A fold on a fraction of the true period averages
   shifted copies of the pattern, a fold on a multiple repeats it, so the period is the shortest whose phase
   structure (variance of fold / per-bin mean - 1) is >= CONTRAST_FRAC of the best; refused when 2-5 times
   that period shows clearly more structure (a pulsed signal). The rate is then corrected from the drift of
   2-4 block segments' folds against the whole-block fold.
4. Confirmation on the fold (<= 128 phases): per phase the peak bin over that bin's floor (10th percentile
   over phases). The carrier is in the captured band for >= MIN_ACTIVE of the period, spans >= MIN_SPAN_HZ,
   moves continuously (<= MAX_FLYBACK jumps and <= MAX_EXITS band exits per period, after a 5-phase median)
   and dominates the span it sweeps (75th percentile of the peak over the span median >= MIN_DOM_DB). Analog
   video line rates (and harmonics) are refused.

Validated on synthetic triangle / sawtooth / sinusoidal / noise-FM sweeps (1.2-50 kHz, 20-50 MHz, down to
+4 dB per MHz) against noise, CW, pulsed and pulsed-LFM radar, TDD OFDM, periodic and random hoppers and FM
video (tests/receivers/test_det_jammer.py), and on the bench capture testdata/live/LIVE_JAMMER_SWEEP_2460MHz.
"""
from __future__ import annotations

import numpy as np
from scipy import fft as sfft
from scipy import sparse

from spectrum_intel.core.emission import Emission

NFFT = 64            # STFT size: fs/64 bins, hop NFFT/2
TRIG_DB = 10.0       # cyclic line over its local median
HARM_DB = 10.0       # a harmonic counts for the rate refinement when this far over the local median
N_CAND = 3           # cyclic-line candidates tried on the fold
MAX_DIV = 8          # the true period is up to this many periods of a cyclic line
CONTRAST_FRAC = 0.8  # the period: shortest fold with this fraction of the best phase structure
MAX_PHASES = 128
N_SEG = 4            # block segments for the rate (drift) fit, at most ...
SEG_PERIODS = 20     # ... and at least this many periods each (>= 2 segments)
ACTIVE_DB = 3.0      # folded phase holds the carrier when its peak bin is this far over the bin floor ...
ACTIVE_REL_DB = 10.0  # ... and >= (median peak of the strongest half of the phases) - this
MIN_ACTIVE = 0.3     # carrier inside the capture for >= 30 % of the period
MIN_SPAN_HZ = 10e6   # narrower periodic FM (analog video lines, FSK) is not a sweep
MAX_FLYBACK = 1      # frequency jumps per period between consecutive phases (a sawtooth's flyback)
MAX_EXITS = 2        # band exits per period (beyond the captured band; clipping images can split the gap)
MIN_DOM_DB = 4.0     # folded peak over the median of the swept span (75th pct of phases; clipping spreads power)
SMOOTH = 5           # circular median over phases (trajectory, jump count); none below 40 phases per period
MAX_OVERSHOOT = 3    # occupied band: at most this many bins beyond the trajectory's range
BRIDGE_DB = 10.0     # a jump is bridged (a fast transit) when the bins between stay within this of the peaks
MIN_RATE_HZ = 300.0
TV_LINES = (15625.0, 15734.264)
TRAJ_POINTS = 32     # trajectory points reported


def occupancy(x, nfft=NFFT):
    """Hann STFT power (hop nfft/2) summed over frame pairs: (n // nfft - 1) rows, one per nfft/fs seconds,
    x nfft bins in FFT order (bin 0 = DC)."""
    x = np.asarray(x, np.complex64)
    h = nfft // 2
    m = len(x) // nfft - 1
    F = np.empty((2 * m, nfft), np.complex64)
    F[0::2] = x[:m * nfft].reshape(m, nfft)
    F[1::2] = x[h:h + m * nfft].reshape(m, nfft)
    F *= np.hanning(nfft).astype(np.float32)
    X = sfft.fft(F, axis=1, overwrite_x=True)
    P = X.real ** 2
    P += X.imag ** 2
    return P.reshape(m, 2, nfft).sum(axis=1)


def _interp(y, k):
    """parabolic interpolation (on log y, in float64) of a peak at index k -> fractional index"""
    if 0 < k < len(y) - 1:
        a, b, c = (float(np.log(float(y[i]) + 1e-30)) for i in (k - 1, k, k + 1))
        d = a - 2 * b + c
        if d < 0:
            return float(k + 0.5 * (a - c) / d)
    return float(k)


def _local_median(S, w=64):
    """median of S over a sliding +-w window (on a w-spaced grid, interpolated)"""
    idx = np.arange(0, len(S), w)
    m = np.array([np.median(S[max(0, i - w):i + w + 1]) for i in idx])
    return np.interp(np.arange(len(S)), idx, m) + 1e-30


def candidates(P, fr, T):
    """P: occupancy rows at rate fr -> (refined line rates, strongest line contrast in dB)"""
    m = P.shape[0]
    med = np.median(P[::8], axis=0) + 1e-20
    Z = np.log(np.maximum(P / med, 1e-3))                # log occupancy: strong random bursts do not dominate
    Z -= Z.mean(axis=0)
    n = sfft.next_fast_len(m, real=True)
    R = sfft.rfft(Z, n=n, axis=0)
    amin, amax = max(MIN_RATE_HZ, 8.0 / T), fr / 8
    lo, hi = int(np.ceil(amin * n / fr)), min(R.shape[0] - 3, int(amax * n / fr))
    if hi - lo < 16:
        return [], 0.0
    A = R.real ** 2 + R.imag ** 2
    S = (A / (A[lo:hi].mean(axis=0) + 1e-30)).sum(axis=1).astype(np.float64)  # per-bin whitening: bursty
                                                                              # bins do not dominate
    C = S / _local_median(S)
    peaks = []
    for k in lo + np.argsort(C[lo:hi])[::-1]:
        if C[k] < 10 ** (TRIG_DB / 10) or len(peaks) >= N_CAND:
            break
        if all(abs(k - j) > 3 and abs(k / j - round(k / j)) * j > 3 and abs(j / k - round(j / k)) * k > 3
               for j in peaks):                           # not a harmonic of a line already taken
            peaks.append(int(k))
    rates = []
    for k in peaks:
        best = _interp(S, k) * fr / n
        for h in range(2, 40):                           # refine on the highest strong harmonic
            if h * best > amax:
                break
            j = int(round(h * best * n / fr))
            if j + 3 > len(S):
                break
            jj = j - 2 + int(np.argmax(S[j - 2:j + 3]))
            if C[jj] >= 10 ** (HARM_DB / 10):
                best = _interp(S, jj) * fr / n / h
        rates.append(best)
    return rates, float(10 * np.log10(C[lo:hi].max()))


def fold(P, fr, period, nph, t0=0.0):
    """mean occupancy per (phase, bin) over all periods (rows at rate fr, the first at time t0)"""
    m = P.shape[0]
    ph = ((t0 + (np.arange(m) + 0.5) / fr) / period % 1.0 * nph).astype(np.int32)
    np.minimum(ph, nph - 1, out=ph)
    M = sparse.csr_matrix((np.ones(m, np.float32), (ph, np.arange(m, dtype=np.int32))), shape=(nph, m))
    return np.asarray(M @ P, np.float64) / np.maximum(np.bincount(ph, minlength=nph), 1)[:, None]


def _norm(F):
    return F / (F.mean(axis=0) + 1e-20) - 1.0


def _period(P, fr, T, period, info):
    """true period from a cyclic line's period (see module doc, step 3) -> (period, phases) or None"""
    amin = max(MIN_RATE_HZ, 8.0 / T)
    cache = {}

    def contrast(p):
        key = round(p * 1e9)
        if key not in cache:
            n = int(min(MAX_PHASES, max(8, p * fr)))
            F = _norm(fold(P, fr, p, n))               # linear: a pattern diluted k times scores 1/k
            cache[key] = (float(np.mean(F * F)), n)
        return cache[key][0]
    ds = [d for d in range(1, MAX_DIV + 1) if 1.0 / (d * period) >= amin]
    if not ds:
        info["reject"] = "sweep period too long for the block"
        return None
    v = {d: contrast(d * period) for d in ds}
    vmax = max(v.values())
    d = next(d for d in ds if v[d] >= CONTRAST_FRAC * vmax)
    period *= d
    nph = cache[round(period * 1e9)][1]
    longer = [contrast(k * period) for k in range(2, 6) if 1.0 / (k * period) >= amin]
    info["contrast"] = round(v[d], 4)
    if longer and max(longer) * CONTRAST_FRAC > v[d]:
        info["reject"] = "the pattern repeats only on a longer period (pulsed, not a continuous sweep)"
        return None
    # rate correction: lag of each segment's fold against the whole-block fold, fitted linearly over time
    ref = _norm(fold(P, fr, period, nph))
    nseg = int(min(N_SEG, max(2, T / period / SEG_PERIODS)))
    q = P.shape[0] // nseg
    lags, mids = [], []
    for i in range(nseg):
        B = _norm(fold(P[i * q:(i + 1) * q], fr, period, nph, t0=i * q / fr))
        cc = np.array([np.sum(ref * np.roll(B, s, axis=0)) for s in range(nph)])
        s0 = int(np.argmax(cc))
        c3 = cc[[(s0 - 1) % nph, s0, (s0 + 1) % nph]] - cc.min() + 1e-9
        sf = s0 + _interp(c3, 1) - 1
        lags.append(sf - nph if sf > nph / 2 else sf)
        mids.append((i + 0.5) * q / fr)
    slope = float(np.polyfit(mids, lags, 1)[0])           # phases per second the pattern runs early
    if abs(slope) * P.shape[0] / fr <= nph / 4:          # a larger drift means the fold did not lock
        period *= 1 - slope * period / nph
    return period, nph


def _confirm(P, fr, T, line_hz, df, info):
    """P: occupancy rows at rate fr, bins low -> high frequency (see module doc, steps 3-4)"""
    if any(abs(line_hz / (h * tv) - 1) < 2e-3 for tv in TV_LINES for h in range(1, 17)):
        info["reject"] = "analog video line rate"
        return None
    got = _period(P, fr, T, 1.0 / line_hz, info)
    if got is None:
        return None
    period, nph = got
    info["rate_hz"] = 1.0 / period
    Pf = fold(P, fr, period, nph)
    E = Pf / (np.percentile(Pf, 10, axis=0) + 1e-20)      # excess over each bin's floor
    k = np.argmax(E, axis=1)
    pk_db = 10 * np.log10(E[np.arange(nph), k])
    strong = np.sort(pk_db)[nph // 2:]
    thr_db = max(ACTIVE_DB, float(np.median(strong)) - ACTIVE_REL_DB)
    active = pk_db >= thr_db
    info["active_frac"] = round(float(active.mean()), 3)
    if active.mean() < MIN_ACTIVE:
        info["reject"] = f"carrier in band for {active.mean():.0%} of the period (< {MIN_ACTIVE:.0%})"
        return None
    act = np.flatnonzero(active)
    # the trajectory, smoothed by a circular median over SMOOTH phases: within a fast transit (smeared by the
    # sweep-to-sweep jitter) a single phase's peak can land on a clipping product or another emitter
    h = min(SMOOTH // 2, nph // 40)                    # fast sweeps (few phases) keep their sharp turns
    fk0 = np.array([_interp(E[i], k[i]) for i in range(nph)])
    fk = np.array([np.median(np.take(fk0, range(i - h, i + h + 1), mode="wrap")) for i in range(nph)])
    ks = np.array([np.median(np.take(k, range(i - h, i + h + 1), mode="wrap")) for i in range(nph)])
    span = (fk[act].max() - fk[act].min()) * df
    lo_b, hi_b = int(ks[act].min()), int(ks[act].max()) + 1
    dom = pk_db - 10 * np.log10(np.median(E[:, lo_b:hi_b], axis=1) + 1e-20)
    gap = ((np.roll(act, -1) - act) % nph) > 1           # the carrier leaves the captured band after act[i]
    # a jump between consecutive phases is a flyback (or a hop) only when some bin between the two peaks is dark in
    # both phases: below the weaker peak - BRIDGE_DB. A fast transit smeared by the sweep jitter keeps the bins
    # between them within 1-6 dB of the peaks (bench jammer); a strong hopper's switching splatter stays ~30 dB down
    flyback = 0
    jump = max(4, 0.3 * span / df)
    for j, a in enumerate(act):
        b = act[(j + 1) % len(act)]
        if gap[j] or abs(ks[b] - ks[a]) <= jump:
            continue
        lo_k, hi_k = int(min(ks[a], ks[b])), int(max(ks[a], ks[b]))
        floor_db = max(ACTIVE_DB, min(pk_db[a], pk_db[b]) - BRIDGE_DB)
        if (np.maximum(E[a, lo_k:hi_k + 1], E[b, lo_k:hi_k + 1]) < 10 ** (floor_db / 10)).any():
            flyback += 1
    exits = int(np.sum(gap))
    dom75 = float(np.percentile(dom[act], 75))
    info.update(span_hz=round(span), flyback=flyback, exits=exits, dom_db=round(dom75, 1))
    if span < MIN_SPAN_HZ:
        info["reject"] = f"sweep span {span / 1e6:.1f} MHz (< {MIN_SPAN_HZ / 1e6:.0f})"
        return None
    if flyback > MAX_FLYBACK or exits > MAX_EXITS:
        info["reject"] = f"{flyback} frequency jumps, {exits} band exits per period (hopping, not a sweep)"
        return None
    if dom75 < MIN_DOM_DB:
        info["reject"] = f"no single dominant carrier per phase ({dom75:.1f} dB < {MIN_DOM_DB})"
        return None
    # occupied extent: from the trajectory's range outward over bins the fold shows lit at some phase (the
    # per-sweep extremes wander ~2 MHz around the folded mean trajectory), at most MAX_OVERSHOOT bins: clipping
    # spreads a saturating jammer's power far wider, and the band must not take in a neighbouring emitter
    nb = P.shape[1]
    lit = E.max(axis=0) >= 10 ** (thr_db / 10)
    lo_lim, hi_lim = max(0, lo_b - MAX_OVERSHOOT), min(nb, hi_b + MAX_OVERSHOOT)
    while lo_b > lo_lim and lit[lo_b - 1]:
        lo_b -= 1
    while hi_b < hi_lim and lit[hi_b]:
        hi_b += 1
    # where the carrier leaves the band: the trajectory on both sides of each gap is near one edge
    beyond = set()
    for i in np.flatnonzero(gap):
        a, b = act[i], act[(i + 1) % len(act)]
        beyond.add("below" if ks[a] + ks[b] < nb - 1 else "above")
    return dict(period=period, nph=nph, fk=fk, active=active, pk_db=pk_db, span=span, flyback=flyback,
                exits=exits, dom75=dom75, beyond=sorted(beyond), occ=(lo_b, hi_b))


def measure(x, fs, fc):
    """-> (sweep dict with absolute frequencies, or None; diagnostics)"""
    P = occupancy(x)
    fr = fs / NFFT                                       # occupancy row rate
    T = len(x) / fs
    df = fs / NFFT
    edge = max(1, NFFT // 16)                            # anti-alias roll-off at the capture edges
    P = np.ascontiguousarray(np.fft.fftshift(P, axes=1)[:, edge:NFFT - edge])  # bins low -> high frequency
    rates, top_db = candidates(P, fr, T)
    diag = {"top_line_db": round(top_db, 1), "tried": []}
    for r in rates:
        info = {"line_hz": round(float(r), 1)}
        diag["tried"].append(info)
        res = _confirm(P, fr, T, r, df, info)
        if res is not None:
            res["freqs"] = fc + (res["fk"] + edge - NFFT // 2) * df
            lo_b, hi_b = res["occ"]                      # occupied bins [lo_b, hi_b) -> band edges
            res["occupied_hz"] = (fc + (lo_b + edge - NFFT // 2 - 0.5) * df, fc + (hi_b + edge - NFFT // 2 - 0.5) * df)
            return res, diag
    return None, diag


def detect(x, fs, fc, ctx=None):
    """-> [Emission] (family jammer, label swept_carrier) or []"""
    res, _ = measure(x, fs, fc)
    if res is None:
        return []
    T = len(x) / fs
    period = res["period"]
    act = np.flatnonzero(res["active"])
    f = res["freqs"]
    f_min, f_max = float(f[act].min()), float(f[act].max())
    lo = fc - fs / 2 if "below" in res["beyond"] else res["occupied_hz"][0]
    hi = fc + fs / 2 if "above" in res["beyond"] else res["occupied_hz"][1]
    nph = res["nph"]
    step = max(1, len(act) // TRAJ_POINTS)
    traj = [[round(float(i) / nph * period * 1e6, 2), round(float(f[i]) / 1e6, 3)] for i in act[::step]]
    feats = {
        "sweep_rate_hz": round(1.0 / period, 2),
        "sweep_period_us": round(period * 1e6, 3),
        "sweep_low_mhz": round(f_min / 1e6, 3),          # lowest / highest carrier frequency seen in the capture
        "sweep_high_mhz": round(f_max / 1e6, 3),
        "sweep_span_mhz": round((f_max - f_min) / 1e6, 3),
        "beyond_capture": res["beyond"],                  # the sweep continues below / above the captured band
        # a sawtooth has one flyback; up-down (triangle / sine) none; a sweep partly outside the band (or broken
        # up by clipping images) may show either, so the shape is only named when the whole sweep is seen
        "shape": "up-down" if not res["flyback"] else "sawtooth" if not res["exits"] else "undetermined",
        "in_band_fraction": round(float(res["active"].mean()), 3),
        "flybacks": res["flyback"],
        "band_exits": res["exits"],
        "dominance_db": round(res["dom75"], 1),
        "fold_peak_db": round(float(np.median(res["pk_db"][act])), 1),
        "trajectory_us_mhz": traj,
        "modulation": "swept CW (FM by a periodic tuning voltage)",
    }
    conf = 0.9 if res["dom75"] >= 6.0 and res["span"] >= 20e6 else 0.75
    return [Emission(family="jammer", label="swept_carrier", start_us=0.0, end_us=T * 1e6,
                     f_lo_mhz=lo / 1e6, f_hi_mhz=hi / 1e6, center_mhz=(lo + hi) / 2e6,
                     snr_db=feats["fold_peak_db"], confidence=conf, features=feats)]
