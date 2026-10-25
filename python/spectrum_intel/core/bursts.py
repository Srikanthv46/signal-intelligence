"""Occupied band and burst extent of a KNOWN link from its in-band power, whatever the waveform.

A receiver finds a link by its structure (CP, sync sequence) and knows roughly where it is; the
structure test is weak per symbol at low SNR, so its bursts come out split, short or late, and an
edge walk on a faded spectrum cuts the band. In-band power over many subcarriers is far more
sensitive (an 18 MHz OFDM link at 0 dB SNR stands ~17 sigma above the floor in every 36 us symbol),
but it is blind to identity. Here power gives the extent, the receiver's structure evidence gives
the membership:

- `TFPower`: one short-time power map of the capture (Hann frames of ~16.7 us).
- `fit_band`: the per-bin median power while the link is on over the median while it is off
  (bursty interferers are rejected by the medians); the band is the contiguous run that maximises
  the summed excess over a level midway between in-band and out-of-band (in-band multipath dips
  and pilot/null patterns do not end it).
- `energy_bursts`: per frame the median over sub-bands of the power over each sub-band's own
  floor (an interferer over less than half the band does not move it). A frame is ON above half
  the link's level, OFF below it, and UNKNOWN when something wider or much stronger covers the band
  (guard bands either side lit, or > 4x the link level). Bursts are ON stretches joined across
  UNKNOWN frames and < 1-symbol dips; a burst is the link's only if the receiver's own evidence
  (seed runs) touches it. Result: the union of those bursts and the seeds.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

FRAME_S = 16.67e-6                   # power-map frame (1024 samples at 61.44 MSps)
MAX_SB = 8                           # sub-bands of the link band (median over them)
SB_HZ = 2e6                          # nominal sub-band width
GUARD_GAP_HZ = 0.3e6                 # guard bands start this far outside the band
GUARD_HZ = 1.0e6                     # and are this wide (>= 10 % of the band)
MAX_BRIDGE_S = 0.6e-3                # longest UNKNOWN stretch joined inside one burst
SMOOTH = 3                           # frames of time smoothing for the ON/OFF decision


class TFPower:
    """|STFT|^2 of the capture: Hann frames of L = 2^round(log2(fs * FRAME_S)) samples, no overlap,
    bins fftshifted (f: offset from fc, Hz). Frame k covers samples [k L, (k + 1) L)."""

    def __init__(self, x, fs, frame_s=FRAME_S):
        x = np.asarray(x)
        self.fs = float(fs)
        self.L = L = max(16, int(2 ** round(np.log2(fs * frame_s))))
        n = len(x) // L
        w = np.hanning(L)
        F = np.fft.fft(x[:n * L].astype(np.complex128).reshape(n, L) * w, axis=1)
        self.P = np.abs(np.fft.fftshift(F, axes=1)) ** 2
        self.dt = L / self.fs
        self.f = (np.arange(L) - L // 2) * (self.fs / L)
        self.n = n

    def cols(self, lo, hi):
        return np.where((self.f >= lo) & (self.f <= hi))[0]

    def frames(self, runs):
        """Boolean per frame: frame entirely inside one of runs [(t0, t1)] seconds."""
        on = np.zeros(self.n, bool)
        for t0, t1 in runs:
            a = int(np.ceil(t0 / self.dt - 1e-9)); b = int(np.floor(t1 / self.dt + 1e-9))
            on[max(a, 0):min(b, self.n)] = True
        return on

    def touched(self, runs, pad=0):
        """Boolean per frame: frame overlaps one of runs, widened by pad frames."""
        m = np.zeros(self.n, bool)
        for t0, t1 in runs:
            a = int(np.floor(t0 / self.dt)) - pad; b = int(np.ceil(t1 / self.dt)) + pad
            m[max(a, 0):min(b, self.n)] = True
        return m


def rect_fit(d):
    """Contiguous [lo, hi) maximising sum(d - tau), tau midway between the mean inside and the median
    outside (iterated from half the smoothed maximum). Returns (lo, hi, inside, outside)."""
    d = np.asarray(d, float)
    n = len(d)
    tau = float(np.max(ndimage.uniform_filter1d(d, 5))) / 2
    lo, hi, inn, out = 0, n, 0.0, 0.0
    for _ in range(4):
        best, s, a = -np.inf, 0.0, 0
        for i in range(n):
            if s <= 0:
                s, a = 0.0, i
            s += d[i] - tau
            if s > best:
                best, lo, hi = s, a, i + 1
        inn = float(np.mean(d[lo:hi]))
        rest = np.concatenate((d[:lo], d[hi:]))
        out = float(np.median(rest)) if rest.size else 0.0
        tau = (inn + out) / 2
    return lo, hi, inn, out


def fit_band(tf, runs, lo, hi, search_lo, search_hi, smooth_hz=0.25e6):
    """Occupied band of a link on during runs [(t0, t1)] s, searched in [search_lo, search_hi] Hz
    (offsets from fc). The fit must contain the centre of [lo, hi]. Returns dict(lo, hi, contrast_db,
    z) or None (no idle or no active frames, or no clear step)."""
    on = tf.frames(runs)
    idle = ~tf.touched(runs, pad=2)
    if on.sum() < 3 or idle.sum() < 3:
        return None
    c = tf.cols(search_lo, search_hi)
    if c.size < 8:
        return None
    A = np.median(tf.P[on][:, c], axis=0)
    I = np.median(tf.P[idle][:, c], axis=0)
    k = max(1, int(round(smooth_hz / (tf.fs / tf.L))))
    d = 10 * np.log10(np.maximum(ndimage.uniform_filter1d(A, k), 1e-30)
                      / np.maximum(ndimage.uniform_filter1d(I, k), 1e-30))
    a, b, inn, out = rect_fit(d)
    f = tf.f[c]
    cen = (lo + hi) / 2
    if not (f[a] <= cen <= f[b - 1]):
        return None
    rest = np.concatenate((d[:a], d[b:]))
    sd = 1.4826 * float(np.median(np.abs(rest - np.median(rest)))) + 1e-9 if rest.size >= 8 else 1e-9
    z = (inn - out) / sd * np.sqrt((b - a) / k)        # step vs the error of the in-band mean
    if inn - out < 1.0 or z < 6.0:
        return None
    df = tf.fs / tf.L
    return dict(lo=float(f[a] - df / 2), hi=float(f[b - 1] + df / 2), contrast_db=round(inn - out, 2), z=float(z))


def _subband_ratio(tf, lo, hi):
    """Per frame (time-smoothed over SMOOTH frames): power of each sub-band of [lo, hi]."""
    c = tf.cols(lo, hi)
    nsb = int(np.clip(round((hi - lo) / SB_HZ), 2, MAX_SB))
    groups = np.array_split(c, nsb)
    p = np.stack([tf.P[:, g].mean(1) for g in groups], 1)
    return p, ndimage.uniform_filter1d(p, SMOOTH, axis=0)


def _guard(tf, a, b, edge):
    c = tf.cols(a, b)
    if c.size < 3 or abs(b) > edge or abs(a) > edge:
        return None
    return tf.P[:, c].mean(1)


def _lower_sigma(v):
    """Noise scale from the lower half (median - 15.87th percentile): the frames away from the seeds
    still hold missed burst parts and interferers, which only ever push values up."""
    return float(np.median(v) - np.percentile(v, 15.87)) + 1e-12


def energy_bursts(tf, seeds, lo, hi, usable_half, min_sym_s=0.0):
    """Bursts of the link in band [lo, hi] (Hz offsets from fc), seeded by the receiver's own runs
    seeds [(t0, t1)] s. Returns (runs [(t0, t1)] s, info). When the link is not visible in power
    (too weak, or never off) the seeds come back unchanged (info["used"] False)."""
    seeds = sorted((float(a), float(b)) for a, b in seeds)
    info = dict(used=False)
    if not seeds or hi - lo <= 0:
        return seeds, info
    raw, sm = _subband_ratio(tf, lo, hi)
    seed_on = tf.frames(seeds)
    away = ~tf.touched(seeds, pad=2)
    if away.sum() < 0.2 * tf.n or seed_on.sum() < 3:
        return seeds, info
    floor = np.maximum(np.percentile(sm[away], 25, axis=0), 1e-30)
    q, qr = sm / floor - 1.0, raw / floor - 1.0
    nsb = q.shape[1]
    # a lit guard band: something just outside the band, which may reach into it; the quarter of the
    # sub-bands next to that guard is left out of the statistic for that frame
    lit = []
    gw = max(GUARD_HZ, 0.1 * (hi - lo))
    for a, b in ((lo - GUARD_GAP_HZ - gw, lo - GUARD_GAP_HZ), (hi + GUARD_GAP_HZ, hi + GUARD_GAP_HZ + gw)):
        g = _guard(tf, a, b, usable_half)
        if g is None:
            lit.append(np.zeros(tf.n, bool))
            continue
        gs = ndimage.uniform_filter1d(g, SMOOTH)
        gr = gs / max(float(np.percentile(gs[away], 25)), 1e-30) - 1.0
        lit.append(gr > 6 * _lower_sigma(gr[away]))
    e = -(-nsb // 4)
    keep_lo, keep_hi = np.where(lit[0], e, 0), np.where(lit[1], nsb - e, nsb)
    both = keep_hi - keep_lo < 1                     # a narrow band with both guards lit: keep all
    keep_lo, keep_hi = np.where(both, 0, keep_lo), np.where(both, nsb, keep_hi)
    # per frame: the lower 25th percentile (a rank, no interpolation) over the kept sub-bands: a frame
    # is on only when three quarters of them are (an interferer must cover most of the band to make a
    # frame look on, however few sub-bands are left)
    r, rr = np.empty(tf.n), np.empty(tf.n)
    for kl in np.unique(keep_lo):
        for kh in np.unique(keep_hi):
            m = (keep_lo == kl) & (keep_hi == kh)
            if m.any():
                r[m] = np.percentile(q[m, kl:kh], 25, axis=1, method="lower")
                rr[m] = np.percentile(qr[m, kl:kh], 25, axis=1, method="lower")
    level = float(np.median(r[seed_on]))
    sig = _lower_sigma(r[away])
    info.update(level_db=round(10 * np.log10(1 + max(level, 0.0)), 2), sigma=round(sig, 4))
    if level < 6 * sig:
        return seeds, info
    thr = max(level / 2, 4 * sig)
    # UNKNOWN: something much stronger over most of the band (the link cannot be judged under it)
    unk = np.median(q, axis=1) > 4 * level + 4 * sig
    # widened by the smoothing length: at an interferer's edges the smoothed median drops below the
    # threshold a frame early while the sub-bands are still lifted
    unk = ndimage.binary_dilation(unk, np.ones(SMOOTH, bool))
    on = (r > thr) & ~unk
    # runs of ON frames joined across gaps of <= MAX_BRIDGE_S made of UNKNOWN frames and at most 2
    # symbols of dip (a fade next to a strong interferer; gaps between TDD bursts are >= 0.4 ms)
    max_unk = int(MAX_BRIDGE_S / tf.dt)
    max_dip = max(1, int(2 * min_sym_s / tf.dt))
    idx = np.where(on)[0]
    runs = []
    for k in idx:
        if runs:
            a, b = runs[-1]
            gap = range(b, k)
            if len(gap) <= max_unk and sum(1 for g in gap if not unk[g]) <= max_dip:
                runs[-1] = (a, k + 1)
                continue
        runs.append((k, k + 1))
    # edges at single-frame resolution: trim / extend to frames on by their own (unsmoothed) power
    edged = []
    for a, b in runs:
        while b - a > 1 and rr[a] < thr and not unk[a]: a += 1
        while b - a > 1 and rr[b - 1] < thr and not unk[b - 1]: b -= 1
        while a > 0 and rr[a - 1] > thr and not unk[a - 1] and not (edged and a - 1 < edged[-1][1]): a -= 1
        while b < tf.n and rr[b] > thr and not unk[b]: b += 1
        edged.append((a, b))
    # membership: the receiver's evidence touches the burst
    touch = tf.touched(seeds)
    mine = [(float(a * tf.dt), float(b * tf.dt)) for a, b in edged if touch[a:b].any()]
    # a seed where the power says clearly OFF (below half the threshold, nothing over the band) is
    # not this link (a short structure false alarm next to an interferer); a seed's edge frames that are
    # clearly OFF are trimmed (frame completion / the symbol grid can run past the burst)
    off = (r < thr / 2) & ~unk
    keep = []
    for t0, t1 in seeds:
        a, b = int(np.floor(t0 / tf.dt)), min(int(np.ceil(t1 / tf.dt)), tf.n)
        if b > a and not unk[a:b].any() and not on[a:b].any() and float(np.mean(rr[a:b])) < thr / 2:
            continue
        while b - a > 1 and off[a]:
            a += 1; t0 = float(a * tf.dt)
        while b - a > 1 and off[b - 1]:
            b -= 1; t1 = float(b * tf.dt)
        keep.append((t0, t1))
    out = []
    for a, b in sorted(mine + keep):
        if out and a <= out[-1][1] + max(tf.dt, min_sym_s):
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    info.update(used=True, thr=round(thr, 4), n_energy=len(mine))
    return out, info
