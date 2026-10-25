"""GPS L1 C/A acquisition from one block, C/N0, jamming and spoofing indicators (research, signal gaps S3;
research/signal_gaps/02_gnss.md sections 7-9).

    analyse(ch, ppm=None, exclude_hz=()) -> dict | None     None when L1 is not in the block

1. Gate: 1575.42 MHz +-1.1 MHz inside the usable band, fs >= 2.1 MSps, >= 16 ms of signal.
2. Decimate: the block FFT cropped to +-1.024 MHz around L1 / (1 + ppm) -> 2.048 MSps; the gate runs on every second
   sample (1.024 MSps, 1 sample / chip).
3. Gate search, all 32 PRNs: Doppler bins of 125 Hz over +-(5 + 3.2) kHz with the receiver ppm (+-25 kHz without);
   per bin the 16 ms are wiped off and folded 4 ms coherently (4 x 1 ms code periods summed) into 4 segments; one
   FFT per segment shared by all PRNs, times each PRN's conjugate replica FFT (code sampled at fs (1 + ppm)), IFFT,
   |.|^2 summed over the segments (non-coherent).
4. Threshold: the noise power from the PRN's whole grid (peak region excluded); 2 M S ~ chi2(2M) for noise, P_fa 1e-3
   per PRN per block split over the cells.
5. Refine each hit at 2.048 MSps over the whole block (4 ms coherent, Doppler +-187.5 Hz in 62.5 Hz steps): code
   phase, Doppler, C/N0 = 10 log10((S - 1) / T_coh) with the parabolic code-straddle correction, plus the
   implementation loss measured on synthetic skies (L_IMPL_DB).
6. Indicators: in-band J/N (L1 +-1.023 MHz power density over the median density of the capture outside 1559-1610
   MHz), narrow lines in band (not receiver lines), spoofing S1 (C/N0 > 55 dB-Hz), S2 (>= 4 PRNs within 1.5 dB at
   > 45 dB-Hz), S3 (J/N >= 1.5 dB with >= 4 PRNs above 45 dB-Hz), S4 (a second correlation peak >= -15 dB in
   >= 2 PRNs).
Depends on numpy and spectrum_intel.core only.
"""
from __future__ import annotations

import math

import numpy as np

from spectrum_intel.core.dsp import Channelizer

from .codes import CHIPS, F_L1, G2_TAPS, RATE, ca_code

FS_DEC = 2.048e6
HALF_BW = 1.024e6
T_COH_MS = 4
M_GATE = 4
BIN_HZ = 125.0
PFA_PRN = 1e-3
SECOND_DB = 15.0
CN0_HIGH, CN0_IMPLAUSIBLE = 45.0, 55.0
# implementation loss added back to C/N0 (crop to +-1.024 MHz ~0.5 dB, nav-bit edges in 4 ms folds, code drift over the
# block, residual straddle): mean -1.8 dB on synthetic skies at 40-48 dB-Hz (tests/receivers/test_gnss.py)
L_IMPL_DB = 1.8
XCORR_NEAR_DB = 15.0
XCORR_DB = 18.0                  # a hit this much weaker than a PRN at its Doppler (mod 1 kHz) is a cross-correlation


def chi2_sf_even(x: float, dof: int) -> float:
    """Survival function of chi2 with an even number of degrees of freedom: e^-x/2 sum_k (x/2)^k / k!."""
    h = x / 2
    term, s = 1.0, 1.0
    for k in range(1, dof // 2):
        term *= h / k
        s += term
    return math.exp(-h) * s


def threshold(m: int, p_cell: float) -> float:
    """S (mean over m of |R|^2 / sigma^2) with P(noise > S) = p_cell: bisection on chi2(2m)."""
    lo, hi = 0.0, 1000.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if chi2_sf_even(2 * m * mid, 2 * m) > p_cell:
            lo = mid
        else:
            hi = mid
    return hi


def _replica_fft(prn: int, n: int, fs_eff: float) -> np.ndarray:
    c = 1.0 - 2.0 * ca_code(prn).astype(np.float64)
    idx = (np.floor(np.arange(n) * (RATE / fs_eff)).astype(np.int64)) % CHIPS
    return np.conj(np.fft.fft(c[idx].astype(np.complex128)))


def _search(y: np.ndarray, fs: float, fs_eff: float, prns: list[int], dopplers: np.ndarray, m: int) -> np.ndarray:
    """Non-coherent grid S[prn, doppler, code sample] (mean over m of |R|^2), 4 ms coherent folds."""
    n = int(round(fs * 1e-3))
    k = T_COH_MS
    x = y[:m * k * n].astype(np.complex128)
    t = np.arange(len(x)) / fs
    R = np.stack([_replica_fft(p, n, fs_eff) for p in prns])            # (P, n)
    out = np.empty((len(prns), len(dopplers), n))
    for j, fd in enumerate(dopplers):
        seg = (x * np.exp(-2j * np.pi * fd * t)).reshape(m, k, n).sum(1)    # (m, n): 4 ms coherent
        F = np.fft.fft(seg, axis=1)
        c = np.fft.ifft(F[None, :, :] * R[:, None, :], axis=2)              # (P, m, n)
        out[:, j, :] = np.mean(np.abs(c) ** 2, axis=1)
    return out


def _peak(g: np.ndarray, chip_samples: float, bin_excl: int):
    """Peak, noise mean (peak region excluded) and the strongest cell outside the peak region."""
    j, i = np.unravel_index(int(np.argmax(g)), g.shape)
    n = g.shape[1]
    dcode = np.abs((np.arange(n) - i + n // 2) % n - n // 2)
    near_code = dcode <= 1.5 * chip_samples
    near_bin = np.abs(np.arange(g.shape[0]) - j) <= bin_excl
    region = near_bin[:, None] & near_code[None, :]
    noise = float(np.mean(g[~region]))
    second = float(np.max(np.where(region, -np.inf, g)))
    return int(j), int(i), float(g[j, i]), noise, second


def _band_power(ch: Channelizer, f0: float, half: float) -> float | None:
    """Mean |X|^2 per bin over [f0 +- half] (absolute Hz), None when outside the capture."""
    sel = np.abs(ch.f + ch.fc - f0) <= half
    sel &= np.abs(ch.f) < 0.45 * ch.fs
    if sel.sum() < 16:
        return None
    return float(np.mean(np.abs(ch.X[sel]).astype(np.float64) ** 2))


PPD_MIN_FS = 30e6              # 32-point frames must stay short (~0.5 us) against a fast sweep
PPD_WINDOW_S = 400e-6           # ~44 sweeps of 9 us: long enough to fold, short enough for a slowly rippling period
PPD_NFFT, PPD_HOP = 32, 16
PPD_RIDGE_DB = 6.0
PPD_COHERENCE = 0.6             # share of the ridge variance the fold explains (noise folded on ~44 periods: ~0.02)


def swept_ppd(ch: Channelizer, f_l1: float) -> dict | None:
    """A fast swept carrier (GNSS 'personal privacy device', research file 05 section 5.2) around L1. The first
    400 us of the band around L1 in 32-point Hann FFT frames every 16 samples: the ridge (peak bin per frame,
    parabolic on log power) is the sweep's frequency over time. Its autocorrelation peak over 0.5-50 us lags
    (parabolic refinement) gives the sweep period, the fold on that period the trajectory. Valid when the ridge
    stands >= 6 dB over each frame's median, the fold explains >= 60 % of the ridge variance and spans >= 1 MHz.
    Shape from the fold's steps over 10 % of the period."""
    if ch.fs < PPD_MIN_FS:
        return None
    half = min(25e6, 0.45 * ch.fs - abs(f_l1 - ch.fc))
    if half < 2e6:
        return None
    y, fs = ch.slice(f_l1, half, ch.fs)
    n = min(len(y), int(PPD_WINDOW_S * fs))
    nf = (n - PPD_NFFT) // PPD_HOP + 1
    if nf < 64:
        return None
    idx = np.arange(nf)[:, None] * PPD_HOP + np.arange(PPD_NFFT)[None, :]
    F = np.fft.fft(y[idx].astype(np.complex128) * np.hanning(PPD_NFFT)[None, :], axis=1)
    P = np.abs(F) ** 2
    k = np.argmax(P, axis=1)
    pk = P[np.arange(nf), k]
    if float(np.median(10 * np.log10(pk / (np.median(P, axis=1) + 1e-30) + 1e-30))) < PPD_RIDGE_DB:
        return None
    a = np.log(P[np.arange(nf), (k - 1) % PPD_NFFT] + 1e-30)
    b = np.log(pk + 1e-30)
    c = np.log(P[np.arange(nf), (k + 1) % PPD_NFFT] + 1e-30)
    den = a - 2 * b + c
    frac = np.where(den < 0, 0.5 * (a - c) / np.where(den < 0, den, -1.0), 0.0)
    kk = np.where(k < PPD_NFFT // 2, k, k - PPD_NFFT).astype(np.float64)
    ridge = (kk + frac) * (fs / PPD_NFFT)
    fr = fs / PPD_HOP
    x = ridge - float(np.mean(ridge))
    var = float(np.mean(x * x))
    m = len(x)
    lo, hi = max(2, int(round(0.5e-6 * fr))), min(int(50e-6 * fr), m // 4)
    if var <= 0 or hi <= lo + 2:
        return None
    r = np.fft.irfft(np.abs(np.fft.rfft(x, 2 * m)) ** 2, 2 * m)[:m]
    top = float(np.max(r[lo + 1:hi - 1])) if hi - lo > 2 else 0.0
    # the shortest lag with a local maximum within 85 % of the strongest: a period that is not a whole number of
    # frames can peak higher at twice the period
    peaks = [q for q in range(lo + 1, hi - 1) if r[q] > r[q - 1] and r[q] >= r[q + 1] and r[q] >= 0.85 * top]
    if not peaks:
        return None                     # no interior maximum: the zero-lag lobe still decays (no repetition)
    j = peaks[0]
    if r[j] / r[0] < 0.5:
        return None
    T = float(j)
    ra, rb, rc = float(r[j - 1]), float(r[j]), float(r[j + 1])        # j is an interior maximum
    d2 = ra - 2 * rb + rc
    if d2 < 0:
        T += 0.5 * (ra - rc) / d2
    nb = max(4, int(round(T)))
    ph = (np.floor((np.arange(m) % T) / T * nb).astype(np.int64)) % nb
    cnt = np.bincount(ph, minlength=nb).astype(np.float64)
    fold = np.bincount(ph, weights=x, minlength=nb) / np.maximum(cnt, 1)
    coh = float(np.mean(fold[ph] ** 2)) / var
    span = float(fold.max() - fold.min())
    if coh < PPD_COHERENCE or span < 1e6:
        return None
    centre = float(np.mean(ridge)) + 0.5 * float(fold.max() + fold.min())
    lag = max(1, nb // 10)
    up = float(np.mean((np.roll(fold, -lag) - fold) > 0))
    shape = "sawtooth up" if up >= 0.85 else "sawtooth down" if up <= 0.15 else "triangle"
    return {"period_us": round(T / fr * 1e6, 3), "span_mhz": round(span / 1e6, 2),
            "centre_offset_mhz": round(centre / 1e6, 2),
            "covers_l1": bool(centre - span / 2 <= 0 <= centre + span / 2), "shape": shape,
            "coherence": round(coh, 3)}


def interference(ch: Channelizer, f_l1: float, exclude_hz=()) -> dict:
    """In-band J/N and narrow lines near L1 (receiver lines excluded)."""
    inband = _band_power(ch, f_l1, 1.023e6)
    f = ch.f + ch.fc
    use = (np.abs(ch.f) < 0.45 * ch.fs) & ((f < 1559e6) | (f > 1610e6))
    ref = float(np.median(np.abs(ch.X[use]).astype(np.float64) ** 2)) / math.log(2) if use.sum() >= 4096 else None
    out = {"jn_db": None if not (inband and ref) else round(10 * math.log10(inband / ref), 2)}
    sel = np.flatnonzero((np.abs(f - f_l1) <= 1.023e6) & (np.abs(ch.f) < 0.45 * ch.fs))
    if len(sel):
        p = np.abs(ch.X[sel]).astype(np.float64) ** 2
        med = float(np.median(p)) + 1e-30
        lines = []
        for k in sel[np.argsort(-p, kind="stable")[:8]]:
            pk = float(np.abs(ch.X[k])) ** 2
            if 10 * math.log10(pk / med) < 25:
                break
            fk = float(f[k])
            if any(abs(fk - e) < 2e3 for e in exclude_hz) or any(abs(fk - g) < 5e3 for g, _ in lines):
                continue
            lines.append((fk, round(10 * math.log10(pk / med), 1)))
        out["lines"] = [{"freq_mhz": round(a / 1e6, 6), "db_over_band": b} for a, b in lines]
    sw = swept_ppd(ch, f_l1)
    if sw:
        out["swept"] = sw
    return out


def analyse(ch: Channelizer, ppm: float | None = None, exclude_hz=()) -> dict | None:
    eps = (ppm or 0.0) * 1e-6
    f_meas = F_L1 / (1 + eps)
    usable = 0.45 * min(ch.fs, ch.fs)
    if ch.fs < 2.1e6 or abs(f_meas - ch.fc) + 1.1e6 > usable or ch.L / ch.fs < 0.016:
        return None
    y, fs = ch.slice(f_meas, HALF_BW, FS_DEC)
    fs_eff = fs * (1 + eps)
    prns = sorted(G2_TAPS)
    span = 5e3 + (3.2e3 if ppm is not None else 20e3)
    dop = np.arange(-span, span + BIN_HZ / 2, BIN_HZ)
    # gate at 1 sample / chip
    yg, fg = y[::2], fs / 2
    if len(yg) < M_GATE * T_COH_MS * int(round(fg * 1e-3)):
        return None
    G = _search(yg, fg, fs_eff / 2, prns, dop, M_GATE)
    n_g = G.shape[2]
    gamma = threshold(M_GATE, PFA_PRN / (n_g * len(dop)))
    hits = []
    for p, g in zip(prns, G):
        j, i, pk, noise, second = _peak(g, fg / RATE, 2)
        s = pk / noise
        if s >= gamma:
            hits.append((p, float(dop[j]), i * 2, s))
    # refine at 2 samples / chip over the whole block
    m_full = max(1, len(y) // (T_COH_MS * int(round(fs * 1e-3))))
    sats = []
    for p, fd, i0, _s in hits:
        fine = fd + np.arange(-3, 4) * (BIN_HZ / 2)
        g = _search(y, fs, fs_eff, [p], fine, m_full)[0]
        j, i, pk, noise, second = _peak(g, fs / RATE, 4)            # +-250 Hz: the 4 ms main lobe
        gam = threshold(m_full, PFA_PRN / (g.shape[1] * len(fine)))
        s = pk / noise
        if s < gam:
            continue
        n = g.shape[1]
        a0, a1, a2 = (math.sqrt(max(g[j, (i - 1) % n], 0)), math.sqrt(pk), math.sqrt(max(g[j, (i + 1) % n], 0)))
        den = a0 - 2 * a1 + a2
        dlt = 0.5 * (a0 - a2) / den if den < 0 else 0.0
        a_pk = a1 - 0.25 * (a0 - a2) * dlt
        s_corr = (a_pk ** 2) / noise
        cn0 = 10 * math.log10(max(s_corr - 1, 1e-9) / (T_COH_MS * 1e-3)) + L_IMPL_DB
        sats.append({"prn": p, "cn0_dbhz": round(cn0, 1), "doppler_hz": round(float(fine[j]) * (1 + eps)),
                     "code_phase_chips": round((-(i + dlt) * RATE / fs_eff) % CHIPS, 2),
                     "peak_over_noise": round(s, 1), "threshold": round(gam, 2),
                     "second_peak_db": round(10 * math.log10(max(second, 1e-30) / pk), 1),
                     "multipeak": bool(second >= gam * noise and 10 * math.log10(second / pk) >= -SECOND_DB)})
    sats.sort(key=lambda d: -d["cn0_dbhz"])
    sats = _xcorr_guard(sats)
    offset = None
    if ppm is None and sats:                          # the receiver's offset: satellites spread +-5 kHz around it
        offset = int(np.median([s["doppler_hz"] for s in sats]))
    itf = interference(ch, f_meas, exclude_hz)
    return {"system": "GPS L1 C/A", "satellites": sats, "n_satellites": len(sats),
            "receiver_offset_hz_estimate": offset,
            "search": {"prns": len(prns), "doppler_span_hz": span, "doppler_bin_hz": BIN_HZ, "coherent_ms": T_COH_MS,
                       "gate_segments": M_GATE, "refine_segments": m_full, "receiver_ppm_applied": ppm},
            "interference": itf, "indicators": indicators(sats, itf)}


def _xcorr_guard(sats: list[dict]) -> list[dict]:
    """Drop a hit that is a cross-correlation of a much stronger PRN: Gold-code cross-correlation peaks are at most
    -24 dB and sit at the strong signal's Doppler (modulo the 1 kHz code-period lines). Several strong signals
    together also leave scattered peaks just over the threshold: a hit >= 15 dB under a PRN and < 1.5 x its
    threshold is dropped too (a real satellite that weak beside strong ones is not claimed)."""
    out = []
    for s in sats:                                    # strongest first
        ghost = any((b["cn0_dbhz"] - s["cn0_dbhz"] >= XCORR_DB and
                     abs((s["doppler_hz"] - b["doppler_hz"] + 500) % 1000 - 500) <= 100) or
                    (b["cn0_dbhz"] - s["cn0_dbhz"] >= XCORR_NEAR_DB and s["peak_over_noise"] < 1.5 * s["threshold"])
                    for b in out)
        if not ghost:
            out.append(s)
    return out


def indicators(sats: list[dict], itf: dict) -> dict:
    cn0 = [s["cn0_dbhz"] for s in sats]
    strong = [c for c in cn0 if c > CN0_HIGH]
    jn = itf.get("jn_db")
    lines = itf.get("lines") or []
    spoof = {"s1_cn0_implausible": bool(cn0 and max(cn0) > CN0_IMPLAUSIBLE),
             "s2_uniform_cn0": bool(len(strong) >= 4 and float(np.std(strong)) < 1.5),
             "s3_power_rise_with_tracking": bool(jn is not None and jn >= 1.5 and len(strong) >= 4),
             "s4_multipeak": sum(s["multipeak"] for s in sats) >= 2}
    n_spoof = sum(spoof.values())
    spoofing = spoof["s3_power_rise_with_tracking"] or (len(sats) >= 4 and
                                                         sum(v for k, v in spoof.items() if k != "s3_power_rise_with_tracking") >= 2)
    jam_floor = jn is not None and jn >= 6.0
    swept = bool(itf.get("swept"))
    jamming = swept or ((jam_floor or bool(lines)) and not spoofing and len(strong) < 4)
    return {"spoofing": {**spoof, "n": n_spoof, "verdict": bool(spoofing)},
            "jamming": {"floor_rise": bool(jam_floor), "narrowband_lines": len(lines), "swept": swept,
                        "verdict": bool(jamming)}}
