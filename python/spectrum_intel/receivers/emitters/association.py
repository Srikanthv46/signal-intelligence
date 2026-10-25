#!/usr/bin/env python3
"""
drone_assoc.py -- emitter association: when bursts hop between frequencies, is it
ONE radio hopping or SEVERAL radios?

Public API
    measure_fingerprints(x, fs, fc, emissions, ch=None, notch_hz=()) -> info
        refines 15 kHz-SCS OFDM emissions in place (features t_start_us, t_sigma_us,
        carrier_hz, carrier_sigma_hz, timing_source, ref_corr, snr_meas_db, ...)
    same_emitter_score(a, b, models=None, prior_same=0.5) -> (llr, evidence)
    associate(emissions, fc=None, models=None, prior_same=0.5, return_matrix=False)
        -> list[Emitter]
    report(emitters) -> str;   analyse_capture(meta_path) -> (emissions, emitters)
    python drone_assoc.py capture.sigmf-meta

Physics (each cue is a log-likelihood ratio  log p(obs|same radio)/p(obs|different
radio of the same type); cues are treated as conditionally independent and summed)

 1 Timing lattice. A TDD radio starts every burst on its own frame clock:
   t = t0 + n*T + slot offset.  A second radio has an independent t0, uniform over
   the lattice spacing.  Burst starts are measured on LEARNED REFERENCE SYMBOLS:
   symbols that repeat burst-to-burst (on the DJI uplink, symbols 0 and 6, Zadoff-
   Chu-like) are found by cross-burst correlation, averaged into a template, and
   every burst is re-timed by the phase slope across subcarriers (~10 ns at
   8-12 dB in-band SNR; CP-train fallback ~0.1-0.5 us; energy edge 5 us).
   Model: Bayesian average of the strict 10 ms frame with slots at +0/+4 ms and a
   0.5 ms slot grid (robust to other sub-patterns); sigma = measurement (+) 30 ns
   jitter (+) 0.02 us/MHz cross-channel group delay (+) |dt| x clock-ppm uncertainty
   (the carrier ppm predicts the frame-clock ppm when one crystal drives both).
 2 Oscillator. Carrier and sample clock come from one crystal, so the carrier
   offset in ppm vs. the nominal raster is the same on every hop at every carrier.
   Carrier = CP-phase fine CFO + DC-null (or matched-filter) integer subcarrier,
   refined by the reference-symbol phase drift (~20 Hz = 0.004 ppm).  Same radio:
   0 +- 0.06 ppm (measured hop-to-hop wander); independent crystals: 0 +- 2.8 ppm.
   Outliers modelled as +-1/2 subcarrier slips plus a small broad term, applied
   to BOTH hypotheses so a gross error cannot vote either way.
 3 Clock vs carrier (emitter level): frame period fitted over the bursts gives
   the frame-clock ppm; consistency with the carrier ppm is reported.
 4 Simultaneity veto: a single transmit chain cannot be on two channels at once.
 5 Power (weak), shape mismatch (negative only), decoded identity / ZC root.

Grouping: emissions of different (family, label) never merge; within a type,
average-linkage clustering on the pairwise LLR matrix (+ prior log-odds), never
across a vetoed pair.  Emitter.confidence = sigmoid(weakest member's mean LLR to
the rest + prior); for singletons, 1 - max P(same) to any other emitter.
Emitter.separation gives, per other emitter of the same type, P(same) and the
per-cue mean LLRs: that is the "same drone or different drone" answer.
"""
from __future__ import annotations
import math
from dataclasses import dataclass, field, asdict
import numpy as np
from scipy.signal import find_peaks

# =================================================================== measurement
# 15 kHz-SCS OFDM (OcuSync / DroneID / LTE-like) at a 15.36 MSps working rate
_R15 = 15.36e6
_N15 = 1024
_CP15 = 72
_SYM15 = _N15 + _CP15
_TS15 = _SYM15 / _R15


def _is_ofdm15(e):
    f = e.features or {}
    return abs(float(f.get("scs_khz", 0)) - 15.0) < 0.5 or e.family in ("dji_ocusync", "dji_droneid")


def _cp_metric(y, N, cp):
    c = y[:-N] * np.conj(y[N:])
    en = (np.abs(y[:-N]) ** 2 + np.abs(y[N:]) ** 2) / 2
    C = np.convolve(c, np.ones(cp), "valid")
    return np.abs(C) / (np.convolve(en, np.ones(cp), "valid") + 1e-30), C


def _parabolic(v, i):
    if i <= 0 or i >= len(v) - 1:
        return 0.0
    a, b, c = v[i - 1], v[i], v[i + 1]
    d = a - 2 * b + c
    return 0.0 if d == 0 else 0.5 * (a - c) / d


def _integer_centre(P, kk, bw_bins, lev_in, lev_out, k_live):
    """Carrier bin of an OFDM symbol power profile P (bins on the subcarrier grid
    after fine-CFO removal). 1) matched filter: the half-integer centre c whose
    width-W window holds the most power relative to its shoulders (W = measured
    -6 dB width); 2) if the band has an unused DC subcarrier (LTE-DL / OcuSync
    style), the single deep notch within +-4 bins of c IS the carrier: exact."""
    o = np.argsort(kk); k, p = kk[o], P[o]
    live = np.abs(k) <= k_live                                 # inside the channelizer slice
    k, p = k[live], p[live]
    lp = 10 * np.log10(p + 1e-30)
    # 1) matched filter with the nominal width, centre on the half-subcarrier grid
    W = int(round(bw_bins))
    best, c_best = -np.inf, None
    lim = int(0.35 * bw_bins)
    for c in np.arange(-lim, lim + 0.01, 0.5):
        inside = np.abs(k - c) <= W / 2
        sh = (np.abs(k - c) > W / 2) & (np.abs(k - c) <= W / 2 + 6)
        if sh.sum() < 4:
            continue
        sc = np.mean(lp[inside]) - np.mean(lp[sh])
        if sc > best:
            best, c_best = sc, c
    if c_best is None:
        return None, 0, False
    # 2) unused DC subcarrier (LTE-DL / OcuSync style): a single-bin notch near the
    #    centre, clearly deeper than any other bin there, pins the carrier exactly
    depth = np.full(len(k), -np.inf)
    for i in range(4, len(k) - 4):
        if abs(k[i]) <= 0.3 * bw_bins and lp[i] <= lp[i - 1] and lp[i] <= lp[i + 1]:
            depth[i] = np.median(np.r_[lp[i - 4:i - 1], lp[i + 2:i + 5]]) - lp[i]
    o2 = np.argsort(depth)[::-1]
    i0 = o2[0]
    edge_db = 10 * np.log10(lev_in / max(lev_out, 1e-30))
    if depth[i0] >= min(7.0, 0.5 * edge_db) and depth[i0] - max(depth[o2[1]], 0) >= 2.5:
        return float(k[i0]), W, True
    return float(c_best), W, False


def _cp_stage(ch, fc, e, half_bw):
    """Symbol grid from the CP train + fine CFO from the CP phase + integer
    subcarrier from the occupied-bin edges. Returns dict or None."""
    R, N, cp, sym = _R15, _N15, _CP15, _SYM15
    f0 = e.center_mhz * 1e6 - fc
    y = ch.slice(f0, half_bw, R)
    a = max(0, int((e.start_us - 150) * 1e-6 * R)); b = min(len(y), int((e.end_us + 150) * 1e-6 * R))
    seg = y[a:b]
    if len(seg) < 3 * sym:
        return None
    m, C = _cp_metric(seg, N, cp)
    hi = max(0.25, 0.5 * float(np.max(m)))
    pk, _ = find_peaks(m, height=hi, distance=int(0.8 * sym))
    # only peaks inside the burst's energy box
    s_rel = (e.start_us * 1e-6 * R) - a
    e_rel = (e.end_us * 1e-6 * R) - a
    pk = pk[(pk >= s_rel - 0.6 * sym) & (pk <= e_rel - 0.5 * sym)]
    if len(pk) < 2:
        return None
    fr = np.array([p + _parabolic(m, p) for p in pk])
    # uniform grid at the nominal symbol period: robust phase
    ph = (fr - fr[0]) % sym
    ph = np.where(ph > sym / 2, ph - sym, ph)
    base = fr[0] + np.median(ph)
    k = np.round((fr - base) / sym)
    res = fr - (base + k * sym)
    good = np.abs(res) < 0.08 * sym
    if good.sum() < 2:
        return None
    base = base + np.mean(res[good])
    res = fr[good] - (base + k[good] * sym)
    kmin = int(np.round((s_rel - base) / sym))                 # first symbol: grid point at the energy start
    t_cp = (a + base + kmin * sym) / R                         # seconds, start of first symbol's CP
    nsym = int(np.round((e_rel - (base + kmin * sym)) / sym))
    sig_t = max(float(np.std(res)) if good.sum() > 2 else 0.25 * cp, 0.05 * cp) / np.sqrt(good.sum()) / R
    # fine CFO from the CP correlation phase (mod 15 kHz)
    cc = C[pk[good]]
    phase = np.angle(np.sum(cc))
    nu = -phase / (2 * np.pi * N / R)
    circ = np.abs(np.sum(cc / (np.abs(cc) + 1e-30))) / len(cc)
    sig_phi = math.sqrt(max(-2 * math.log(max(circ, 1e-6)), 1e-4)) / math.sqrt(len(cc))
    sig_nu = sig_phi / (2 * np.pi * N / R)
    # integer subcarrier: occupied bins after removing the fine CFO
    y2 = ch.slice(f0 + nu, half_bw * 1.6, R)
    nsym = max(nsym, 1)
    P = np.zeros(N)
    for j in range(nsym):
        s = int(round(t_cp * R)) + j * sym + cp // 2
        if 0 <= s and s + N <= len(y2):                        # grid can start before the capture
            P += np.abs(np.fft.fft(y2[s:s + N])) ** 2
    kk = np.fft.fftfreq(N) * N
    bw_bins = (e.f_hi_mhz - e.f_lo_mhz) * 1e6 / (R / N)
    inb = np.abs(kk) < 0.3 * bw_bins
    outb = (np.abs(kk) > 0.75 * bw_bins) & (np.abs(kk) < 1.5 * half_bw / (R / N))
    lev_in = np.median(P[inb]); lev_out = np.median(P[outb]) if outb.any() else lev_in * 1e-3
    snr_db = 10 * np.log10(max(lev_in / max(lev_out, 1e-30) - 1, 1e-3))
    k_c, n_occ, dc_null = _integer_centre(P, kk, bw_bins, lev_in, lev_out, 1.6 * half_bw / (R / N) - 3)
    if k_c is None:
        return None
    carrier = fc + f0 + nu + k_c * R / N
    return dict(y=None, f_slice=f0 + nu + k_c * R / N, t_cp=t_cp, sig_t=sig_t, nsym=nsym, carrier=carrier,
                sig_nu=sig_nu, snr_db=snr_db, n_occ=n_occ, dc_null=dc_null, cp_peaks=int(good.sum()),
                edge_db=float(10 * np.log10(lev_in / max(lev_out, 1e-30))))


def _symbols(ch, f_slice, t_cp, nsym, half_bw, kmax):
    """FFT of each symbol (window starting mid-CP), occupied bins |k| <= kmax only."""
    R, N, cp, sym = _R15, _N15, _CP15, _SYM15
    y = ch.slice(f_slice, half_bw, R)
    s0 = t_cp * R
    kk = np.fft.fftfreq(N) * N
    sel = np.abs(kk) <= kmax
    S = []
    for j in range(nsym):
        s = s0 + j * sym + cp / 2
        a = int(np.floor(s)); frac = s - a
        if a + N + 1 > len(y) or a < 0:
            break
        X = np.fft.fft(y[a:a + N])
        X *= np.exp(2j * np.pi * kk * frac / N)               # fractional-sample window start
        S.append(X[sel])
    return np.array(S), kk[sel]


def _delay_xcorr(A, B, k, up=16):
    """delay (samples @R, A later than B => positive) maximising |sum A B* e^{+j2pi k d/N}|, + normalised peak"""
    N = _N15
    Z = np.zeros(N * up, complex)
    idx = (k.astype(int)) % (N * up)
    np.add.at(Z, idx, A * np.conj(B))
    r = np.fft.ifft(Z)
    lim = int(6e-6 * _R15 * up)                                # search +-6 us
    cand = np.r_[np.arange(0, lim), np.arange(N * up - lim, N * up)]
    i = cand[np.argmax(np.abs(r[cand]))]
    d = (i if i < N * up / 2 else i - N * up) / up
    return d, float(np.abs(r[i]) * N * up / (np.linalg.norm(A) * np.linalg.norm(B) + 1e-30))


def _fit_delay_phase(prods, k, tj):
    """prods[j] = X_b,j * conj(T_j) per bin. Fit common delay d (samples) and
    phase phi_j = phi0 + 2 pi dnu t_j. Returns d, dnu, sig_d, sig_dnu, rho."""
    N = _N15
    P = np.sum(prods, 0)
    # coarse delay from the summed product
    up = 16
    Z = np.zeros(N * up, complex); np.add.at(Z, k.astype(int) % (N * up), P)
    r = np.fft.ifft(Z)
    lim = int(6e-6 * _R15 * up)
    cand = np.r_[np.arange(0, lim), np.arange(N * up - lim, N * up)]
    i = cand[np.argmax(np.abs(r[cand]))]
    d = (i if i < N * up / 2 else i - N * up) / up             # product ~ exp(-j 2 pi k d / N)
    # fine delay: weighted phase slope across bins of the de-rotated product
    for _ in range(2):
        Q = P * np.exp(2j * np.pi * k * d / N)
        q = Q * np.exp(-1j * np.angle(np.sum(Q)))
        w = np.abs(q) + 1e-30
        ang = np.angle(q)
        A = np.vstack([k, np.ones_like(k)]).T
        W = w / w.sum()
        sl, ic = np.linalg.lstsq(A * np.sqrt(W)[:, None], ang * np.sqrt(W), rcond=None)[0]
        d -= sl * N / (2 * np.pi)
    Q = P * np.exp(2j * np.pi * k * d / N)
    q = Q * np.exp(-1j * np.angle(np.sum(Q)))
    resid = np.angle(q)
    sig_sc = float(np.sqrt(np.sum(np.abs(q) * resid ** 2) / np.sum(np.abs(q))))
    sig_d = sig_sc / (2 * np.pi / N * np.sqrt(np.sum((k - k.mean()) ** 2)))
    # per-symbol phases at the common delay -> dnu
    ramp = np.exp(2j * np.pi * k * d / N)
    ph = np.array([np.angle(np.sum(p * ramp)) for p in prods])
    rho = float(np.abs(np.sum(Q)) / (np.sum(np.abs(prods)) + 1e-30))
    if len(prods) >= 2:
        phu = np.unwrap(ph)
        dnu, _ = np.polyfit(tj, phu, 1) / (2 * np.pi)
        sig_ph = sig_sc / np.sqrt(len(k)) * np.sqrt(len(prods))
        sig_dnu = sig_ph / (2 * np.pi * np.sqrt(np.sum((tj - tj.mean()) ** 2)))
    else:
        dnu, sig_dnu = 0.0, np.inf
    return float(d), float(dnu), float(sig_d), float(sig_dnu), rho


def _template_stage(ch, items, half_bw, min_ref_corr=0.5, min_burst_corr=0.35):
    """Learn the symbols that repeat burst-to-burst (reference / pilot symbols) and
    re-time every burst against them: preamble-quality timing (ns) and CFO (Hz).
    items: list of dicts from _cp_stage (same numerology). Updates in place."""
    R, N, sym = _R15, _N15, _SYM15
    if len(items) < 2:
        return dict(ref_symbols=[], note="single burst: no learned reference")
    kmax = int(np.median([it["n_occ"] for it in items]) // 2)
    nsym = int(min(np.median([it["nsym"] for it in items]), min(it["nsym"] for it in items)))
    if nsym < 2:
        return dict(ref_symbols=[], note="bursts too short")
    S = []
    for it in items:
        s, k = _symbols(ch, it["f_slice"], it["t_cp"], nsym, half_bw, kmax)
        S.append(s)
    nb = len(items)
    if any(len(s) < nsym for s in S):
        return dict(ref_symbols=[], note="burst runs off the capture")
    # which symbol indices repeat across bursts?
    order = np.argsort([-it["snr_db"] for it in items])
    pairs = [(order[0], j) for j in order[1:]] + [(order[1], j) for j in order[2:6]]
    rho = np.zeros(nsym)
    for jj in range(nsym):
        rho[jj] = np.median([_delay_xcorr(S[a][jj], S[b][jj], k)[1] for a, b in pairs])
    J = [int(j) for j in np.where(rho > min_ref_corr)[0]]
    if not J:
        return dict(ref_symbols=[], symbol_repeat_corr=np.round(rho, 2).tolist(), note="no repeating symbol")
    tj = np.array(J) * sym / R
    T = np.array([S[order[0]][j] for j in J])
    # integer-subcarrier repair: a burst whose carrier bin was mis-resolved does not
    # correlate with the reference; re-cut it +-1..3 subcarriers off and keep the best
    n_fixed = 0
    for b in range(nb):
        r0 = _fit_delay_phase([S[b][j] * np.conj(T[i]) for i, j in enumerate(J)], k, tj)[4]
        if r0 >= 0.5:
            continue
        best = (r0, 0, None)
        for m in (-1, 1, -2, 2, -3, 3):
            sm, _ = _symbols(ch, items[b]["f_slice"] + m * R / N, items[b]["t_cp"], nsym, half_bw, kmax)
            if len(sm) < nsym:
                continue
            r = _fit_delay_phase([sm[j] * np.conj(T[i]) for i, j in enumerate(J)], k, tj)[4]
            if r > best[0]:
                best = (r, m, sm)
        if best[1] and best[0] > 0.5 and best[0] > 1.5 * r0:
            m = best[1]
            items[b]["f_slice"] += m * R / N
            items[b]["carrier"] += m * R / N
            items[b]["int_fix"] = m
            S[b] = best[2]
            n_fixed += 1
    fits = None
    for it_ in range(3):
        fits = []
        acc = np.zeros_like(T)
        for b in range(nb):
            prods = [S[b][j] * np.conj(T[i]) for i, j in enumerate(J)]
            d, dnu, sd, sdn, r = _fit_delay_phase(prods, k, tj)
            fits.append((d, dnu, sd, sdn, r))
            # align burst b onto the template and accumulate
            ramp = np.exp(2j * np.pi * k * d / N)
            ph0 = np.angle(np.sum([p * ramp for p in prods]))
            for i, j in enumerate(J):
                acc[i] += S[b][j] * ramp * np.exp(-1j * (ph0 + 2 * np.pi * dnu * (tj[i] - tj.mean())))
        T = acc / nb
    d = np.array([f[0] for f in fits]); dnu = np.array([f[1] for f in fits])
    # the template inherits the reference bursts' CP errors; CP errors are zero-mean,
    # so referencing to the group mean of the CP estimates removes that common bias
    t_cp = np.array([it["t_cp"] for it in items])
    sig_cp = np.array([it["sig_t"] for it in items])
    t_new = t_cp + d / R
    w = 1 / sig_cp ** 2
    bias_t = np.sum(w * (t_new - t_cp)) / np.sum(w)
    car_cp = np.array([it["carrier"] for it in items])
    sig_nu_cp = np.array([it["sig_nu"] for it in items])
    car_new = car_cp + dnu
    w = 1 / sig_nu_cp ** 2
    bias_f = np.sum(w * (car_new - car_cp)) / np.sum(w)
    for b, it in enumerate(items):
        dt = t_new[b] - bias_t - t_cp[b]
        df = car_new[b] - bias_f - car_cp[b]
        ok = (fits[b][4] > min_burst_corr and abs(dt) < max(5 * sig_cp[b], 2e-6)
              and abs(df) < max(5 * sig_nu_cp[b], 1200.0))
        it["ref_corr"] = round(fits[b][4], 3)
        if ok:
            it["t_fine"] = t_new[b] - bias_t
            it["sig_t_fine"] = math.hypot(fits[b][2] / R, 3e-9)
            if fits[b][4] >= 0.6:      # weak match: timing holds up, the phase-drift CFO slips -> keep CP CFO
                it["carrier_fine"] = car_new[b] - bias_f
                it["sig_nu_fine"] = math.hypot(fits[b][3], 2.0)
            else:
                it["carrier_fine"], it["sig_nu_fine"] = it["carrier"], it["sig_nu"]
        it["ref_ok"] = bool(ok)
    return dict(ref_symbols=J, symbol_repeat_corr=np.round(rho, 2).tolist(), integer_fixes=n_fixed,
                cp_vs_ref_timing_rms_us=round(float(np.std(t_new - t_cp) * 1e6), 3))


def measure_fingerprints(x, fs, fc, emissions, ch=None, notch_hz=()):
    """Refine the association features of 15 kHz-SCS OFDM emissions in place:
        t_start_us / t_sigma_us     first-symbol CP start: CP-train fit (~0.1-0.5 us),
                                    then re-timed on learned reference symbols (ns)
        carrier_hz / carrier_sigma_hz   absolute carrier (CP phase + occupied-bin
                                    edges, refined by reference-symbol phase drift)
        snr_db, n_occ_sc, ref_corr, timing_source
    Bursts are grouped by (family, label) for the learned reference."""
    from spectrum_intel.core import spectral as SS
    if ch is None:                     # notch CW spurs / LO leakage: a tone correlates at every lag
        notch = list(notch_hz) or [(sp["offset_hz"], 30e3) for sp in SS.find_spurs(x, fs, fc)]
        ch = SS.Channelizer(np.asarray(x, np.complex128), fs, notch)
    groups = {}
    for i, e in enumerate(emissions):
        if _is_ofdm15(e):
            groups.setdefault((e.family, e.label, round((e.f_hi_mhz - e.f_lo_mhz) / 2)), []).append(i)
    info = {}
    for key, idx in groups.items():
        bw = np.median([emissions[i].f_hi_mhz - emissions[i].f_lo_mhz for i in idx]) * 1e6
        half_bw = min(bw / 2 * 1.08, _R15 / 2 * 0.95)
        items, ids = [], []
        for i in idx:
            it = _cp_stage(ch, fc, emissions[i], half_bw)
            if it is not None:
                items.append(it); ids.append(i)
        tinfo = _template_stage(ch, items, half_bw) if len(items) >= 2 else {}
        info[str(key)] = tinfo
        for i, it in zip(ids, items):
            f = emissions[i].features
            fine = it.get("ref_ok", False)
            f["t_start_us"] = (it["t_fine"] if fine else it["t_cp"]) * 1e6
            f["t_sigma_us"] = (it["sig_t_fine"] if fine else it["sig_t"]) * 1e6
            f["carrier_hz"] = it["carrier_fine"] if fine else it["carrier"]
            f["carrier_sigma_hz"] = it["sig_nu_fine"] if fine else it["sig_nu"]
            f["t_cp_us"] = it["t_cp"] * 1e6
            f["carrier_cp_hz"] = it["carrier"]
            f["timing_source"] = "reference_symbols" if fine else "cp_train"
            f["ref_corr"] = it.get("ref_corr")
            f["n_occ_sc"] = it["n_occ"]
            f["cp_peaks"] = it["cp_peaks"]
            if emissions[i].snr_db is None:
                emissions[i].snr_db = round(float(it["snr_db"]), 1)
            f["snr_meas_db"] = round(float(it["snr_db"]), 1)
    return info


# =================================================================== physics models
# Per (family, label) association model. Every number is a physical prior with
# its origin noted; override per call through associate(..., models={...}).
#
#   lattices   [(frame_us, offsets_us, weight), ...]  burst starts of ONE radio lie on
#              t0 + n*frame + offset (its own frame clock). Several lattices = Bayesian
#              model average (a strict sub-pattern plus a looser slot grid that holds
#              even if the sub-pattern differs on another firmware/model).
#   jitter_us  start-time scatter of one radio around its lattice, same channel
#   chan_delay_us_per_mhz  extra start offset between DIFFERENT channels (receiver +
#              transmitter group delay varies across the band)
#   raster     (offset_hz, step_hz) nominal channel raster -> per-burst carrier ppm
#   ppm_h1     carrier-ppm spread of ONE radio between hops (AFC wander, Doppler)
#   ppm_pop    crystal error spread across DIFFERENT radios (per radio, 1 sigma)
#   single_tx  one transmit chain: cannot occupy two channels at the same time
DJI_UPLINK = dict(
    # Real DJI Mini 2 RC uplink (D_20260916-171343): 0.5 ms 15 kHz-OFDM slots every
    # 4 / 6 ms (10 ms frame, slots at +0 and +4 ms). Slot grid 0.5 ms = LTE-like slot.
    lattices=[(10000.0, (0.0, 4000.0), 0.6), (500.0, (0.0,), 0.4)],
    # measured on the same capture with reference-symbol timing: same-channel revisits
    # sit 13-30 ns off the lattice; channels 22 MHz apart differ by ~0.37 us
    jitter_us=0.03, chan_delay_us_per_mhz=0.02,
    # carriers of all 5 hops (4 channels) = raster x.5 MHz - 16 kHz; visit-to-visit
    # wander 300-400 Hz at 5.8 GHz = 0.05-0.07 ppm
    raster=(0.0, 0.5e6), ppm_h1=0.06, ppm_drift_per_s=0.05, ppm_pop=2.0,
    snr_h1_db=4.0, single_tx=True)
DJI_DRONEID = dict(lattices=[], jitter_us=1.0, chan_delay_us_per_mhz=0.02, raster=(0.0, 0.5e6),
                   ppm_h1=0.06, ppm_drift_per_s=0.05, ppm_pop=2.0, snr_h1_db=4.0, single_tx=True)
DEFAULT_MODEL = dict(variable_length=False, lattices=[], jitter_us=2.0, chan_delay_us_per_mhz=0.05, raster=None,
                     ppm_h1=0.1, ppm_drift_per_s=0.1, ppm_pop=5.0, snr_h1_db=5.0, single_tx=True)
MODELS = {("dji_ocusync", "ocusync_uplink"): DJI_UPLINK,
          ("dji_droneid", "droneid"): DJI_DRONEID,
          ("dji_ocusync", "ocusync_video"): dict(DJI_DRONEID, variable_length=True)}
# det_dji video labels: variable-length bursts, link-level carrier from the detector
for _lab in ("video_10MHz", "video_20MHz", "video_o4_10MHz", "video_o4_20MHz"):
    MODELS[("dji_ocusync", _lab)] = dict(DJI_DRONEID, variable_length=True)
MODELS[("dji_ocusync", "uplink")] = DJI_UPLINK                 # det_dji's label for the RC uplink

EPS_OUTLIER = 0.03          # per-cue probability that a same-radio measurement is an outlier
                            # (mis-resolved subcarrier, collision, energy-edge start): caps negative LLRs
EPS_TIMING_FINE = 0.001     # ... for two reference-symbol timed bursts (frame-clock slip in one capture)
EPS_BROAD = 0.003           # broad gross-error component of the oscillator cue
EPS_BROAD_FINE = 3e-4       # ... when both carriers come from the reference-symbol fit
LLR_VETO = -12.0            # a pair this unlikely is never merged


def _model(e, models=None):
    m = dict(DEFAULT_MODEL)
    m.update((models or MODELS).get((e.family, e.label), MODELS.get((e.family, e.label), {})))
    f = e.features or {}
    if f.get("frame_period_ms") and not m["lattices"]:          # detector-measured frame clock
        m["lattices"] = [(float(f["frame_period_ms"]) * 1e3, (0.0,), 1.0)]
    return m


# ------------------------------------------------------------------ feature access
def _t(e):
    f = e.features or {}
    if "t_start_us" in f:
        return float(f["t_start_us"]), float(f.get("t_sigma_us", 0.5))
    if "frame_phase_us" in f and "t_sigma_us" in f:
        return float(e.start_us), float(f["t_sigma_us"])
    return float(e.start_us), 5.0                                  # energy edge (+-several us)


def _carrier(e):
    f = e.features or {}
    if "carrier_hz" in f:
        return float(f["carrier_hz"]), float(f.get("carrier_sigma_hz", 200.0))
    return None, None


def _ppm(e, m):
    """Carrier offset in ppm vs the nominal channel (same crystal => same ppm on every hop)."""
    f = e.features or {}
    c, sc = _carrier(e)
    # a carrier refined on learned reference symbols has its integer subcarrier locked;
    # detector CFOs can be whole subcarriers off (e.g. a hop beside the LO spur), so it wins
    refined = c is not None and m.get("raster") and f.get("timing_source") == "reference_symbols"
    if f.get("cfo_ppm") is not None and not refined:
        return float(f["cfo_ppm"]), float(f.get("cfo_ppm_sigma", 0.05)), "detector"
    if c is not None and m.get("raster"):
        off, step = m["raster"]
        nom = off + round((c - off) / step) * step
        return (c - nom) / nom * 1e6, sc / nom * 1e6, f"raster {step/1e6:g} MHz"
    if f.get("cfo_hz") is not None and f.get("cfo_ref") == "nominal":
        cen = e.center_mhz * 1e6
        return float(f["cfo_hz"]) / cen * 1e6, float(f.get("cfo_sigma_hz", 200.0)) / cen * 1e6, "detector cfo_hz"
    return None, None, None


def _ids(e):
    out = {}
    for src in (e.decoded or {}), (e.features or {}):
        for k in ("serial", "serial_number", "drone_serial", "mac", "uuid"):
            if src.get(k):
                out["id"] = str(src[k]).lower()
        for k in ("zc_root", "pilot_id", "sync_word", "bind_id"):
            if src.get(k) is not None:
                out[k] = src[k]
    return out


def _bw(e):
    return float((e.features or {}).get("bw_mhz", e.f_hi_mhz - e.f_lo_mhz))


# ------------------------------------------------------------------ cue likelihoods
_SQ2PI = math.sqrt(2 * math.pi)


def _npdf(r, s):
    return math.exp(-0.5 * (r / s) ** 2) / (s * _SQ2PI)


def _lnpdf(r, s):
    return -0.5 * (r / s) ** 2 - math.log(s * _SQ2PI)


def _lmix(d, terms):
    """log of a Gaussian mixture sum(w * N(d - mu; s)), computed in the log domain: for a
    large d every term underflows to 0 and a likelihood ratio becomes 0/0."""
    v = [math.log(w) + _lnpdf(d - mu, sg) for w, mu, sg in terms if w > 0]
    m = max(v)
    return m + math.log(sum(math.exp(x - m) for x in v))


def _lattice_ratio(dt, sig, frame, offsets):
    """p(dt | same radio, lattice) / p(dt | independent radios). Under H0 the phase of
    dt modulo the frame is uniform; under H1 it is a wrapped normal around the allowed
    differences D = {o_j - o_i}. Returns (ratio, residual_us, best_d)."""
    D = sorted({(b - a) % frame for a in offsets for b in offsets})
    ph = dt % frame
    tot, best = 0.0, (np.inf, 0.0)
    for d in D:
        for kk in range(-3, 4):
            r = ph - d + kk * frame
            tot += _npdf(r, sig)
            if abs(r) < abs(best[0]):
                best = (r, d)
    return tot * frame / len(D), best[0], best[1]


def _cue_timing(a, b, ma, ppm_pair):
    if not ma["lattices"]:
        return None
    ta, sa = _t(a); tb, sb = _t(b)
    dt = tb - ta
    if ppm_pair is not None:                       # one crystal: the frame clock runs ppm off
        e, se = ppm_pair, 1.0                      # like the carrier (SoC-dependent: allow 1 ppm)
    else:
        e, se = 0.0, 25.0
    dfm = abs(a.center_mhz - b.center_mhz)
    same_ch = dfm < 0.5 * max(_bw(a), 0.2)
    sig = math.sqrt(sa ** 2 + sb ** 2 + ma["jitter_us"] ** 2
                    + (0 if same_ch else (ma["chan_delay_us_per_mhz"] * dfm) ** 2)
                    + (abs(dt) * se * 1e-6) ** 2)
    # a radio slipping its frame clock between two bursts of one capture is rare; a bad
    # start estimate is the realistic outlier, so the floor depends on timing quality
    fine = all((e.features or {}).get("timing_source") == "reference_symbols" for e in (a, b))
    eps = EPS_TIMING_FINE if fine else EPS_OUTLIER
    ratio, best = 0.0, None
    for frame, offs, w in ma["lattices"]:
        fr = frame * (1 - e * 1e-6)
        rr, res, d = _lattice_ratio(dt, sig, fr, [o * (1 - e * 1e-6) for o in offs])
        ratio += w * ((1 - eps) * rr + eps)
        if best is None or frame < best[3]:          # report vs. the finest (slot) grid
            best = (rr, res, d, frame)
    llr = math.log(max(ratio, 1e-300))
    rr, res, d, frame = best
    n = round((dt - d) / frame) if frame else 0
    W = min(fr_ / len(o_) for fr_, o_, _ in ma["lattices"])
    return dict(cue="timing_lattice", llr=llr,
                detail=(f"dt {dt/1e3:+.4f} ms = {n:d} x {frame/1e3:g} ms + {d/1e3:g} ms "
                        f"{res*1e3:+.0f} ns (sigma {sig*1e3:.0f} ns{'' if same_ch else ', cross-channel'}); "
                        f"an independent radio's phase is uniform over ~{W/1e3:g} ms"),
                residual_us=res, sigma_us=sig)


def _cue_oscillator(a, b, ma):
    pa, sa, srca = _ppm(a, ma); pb, sb, srcb = _ppm(b, ma)
    ca, sca = _carrier(a); cb, scb = _carrier(b)
    dtau = abs(_t(b)[0] - _t(a)[0]) * 1e-6
    if pa is None and ca is not None and cb is not None and abs(ca - cb) < 100e3:
        # no raster known, but a revisit of the same channel compares carriers directly
        d = (cb - ca) / ca * 1e6; s = math.hypot(sca, scb) / ca * 1e6
        how = "same-channel carrier difference"
    elif pa is not None and pb is not None:
        d = pb - pa; s = math.hypot(sa, sb)
        how = f"carrier ppm vs {srca}"
    else:
        return None
    s1 = math.sqrt(s ** 2 + ma["ppm_h1"] ** 2 + (ma["ppm_drift_per_s"] * dtau) ** 2)
    s0 = math.sqrt(s ** 2 + 2 * ma["ppm_pop"] ** 2)
    # a gross measurement outlier (mis-resolved subcarrier, collision) is equally likely
    # under both hypotheses: same broad component on both sides, so it cannot vote
    # outliers: a mis-resolved subcarrier (an exact +-1/+-2 x SCS carrier error) and a
    # small broad component; both equally likely under either hypothesis
    scs = float((a.features or {}).get("scs_khz") or 0) * 1e3
    fine = all((e.features or {}).get("timing_source") == "reference_symbols" for e in (a, b))
    out = [(EPS_BROAD_FINE if fine else EPS_BROAD, 0.0, 10.0 * ma["ppm_pop"])]
    if scs:
        step = scs / (a.center_mhz * 1e6) * 1e6
        out += [(EPS_OUTLIER / 4, k * step, s1) for k in (-2, -1, 1, 2)]
    llr = _lmix(d, [(1 - EPS_OUTLIER, 0.0, s1)] + out) - _lmix(d, [(1 - EPS_OUTLIER, 0.0, s0)] + out)
    det = (f"{how}: {pa:+.3f} vs {pb:+.3f} ppm, " if pa is not None else f"{how}: ") + \
          (f"diff {d:+.3f} ppm ({d * a.center_mhz:+.0f} Hz); one crystal => 0 +- {s1:.3f}, "
           f"independent crystals => 0 +- {s0:.2f} ppm")
    return dict(cue="oscillator_ppm", llr=llr, detail=det, dppm=d)


def _cue_clock(a, b, ma):
    fa, fb = a.features or {}, b.features or {}
    if fa.get("clock_ppm") is None or fb.get("clock_ppm") is None:
        return None
    if a.family.startswith("rc_") or b.family.startswith("rc_"):
        return None                     # det_rc_link's chip-clock estimate is biased (10-50 ppm scatter)
    d = float(fb["clock_ppm"]) - float(fa["clock_ppm"])
    s = math.hypot(float(fa.get("clock_ppm_sigma", fa.get("clock_ppm_err", 2.0)) or 2.0),
                   float(fb.get("clock_ppm_sigma", fb.get("clock_ppm_err", 2.0)) or 2.0))
    s1 = math.hypot(s, 0.2); s0 = math.sqrt(s ** 2 + 2 * ma["ppm_pop"] ** 2)
    out = [(EPS_OUTLIER, 0.0, 10 * ma["ppm_pop"])]
    return dict(cue="sample_clock_ppm", llr=_lmix(d, [(1 - EPS_OUTLIER, 0.0, s1)] + out)
                - _lmix(d, [(1 - EPS_OUTLIER, 0.0, s0)] + out),
                detail=f"symbol-clock ppm diff {d:+.2f} (sigma {s:.2f})")


def _cue_overlap(a, b, ma):
    ov = min(a.end_us, b.end_us) - max(a.start_us, b.start_us)
    short = max(min(a.end_us - a.start_us, b.end_us - b.start_us), 1.0)
    sep = abs(a.center_mhz - b.center_mhz) > 0.5 * (_bw(a) + _bw(b)) / 2
    if ma["single_tx"] and ov > 0.2 * short and sep:
        return dict(cue="simultaneous", llr=-20.0,
                    detail=f"on air at the same time ({ov:.0f} us overlap) on channels "
                           f"{abs(a.center_mhz - b.center_mhz):.2f} MHz apart: one transmit chain cannot")
    return None


def _cue_power(a, b, ma):
    if a.snr_db is None or b.snr_db is None:
        return None
    d = float(b.snr_db) - float(a.snr_db)
    s1, s0 = ma["snr_h1_db"], 10.0
    out = [(EPS_OUTLIER, 0.0, 25.0)]
    return dict(cue="power", llr=_lmix(d, [(1 - EPS_OUTLIER, 0.0, s1)] + out) - _lmix(d, [(1 - EPS_OUTLIER, 0.0, s0)] + out),
                detail=f"SNR {a.snr_db:.1f} vs {b.snr_db:.1f} dB (hop-to-hop fading ~{s1:g} dB; "
                       f"different ranges ~{s0:g} dB)")


def _cue_shape(a, b, ma=None):
    """Only mismatches are informative (two radios of one model look identical).
    Burst duration is compared only for fixed-length bursts (uplink slots, DroneID);
    video bursts carry a variable amount of data, so their length says nothing."""
    ba, bb = _bw(a), _bw(b)
    z = abs(ba - bb) / (0.08 * max(ba, bb) + 0.05)
    da, db = a.end_us - a.start_us, b.end_us - b.start_us
    if not (ma or {}).get("variable_length"):
        z = max(z, abs(da - db) / (0.1 * max(da, db) + 20.0))
    if z <= 2:
        return None
    return dict(cue="shape", llr=-min(0.5 * (z ** 2 - 4), 8.0),
                detail=f"bandwidth {ba:.2f} vs {bb:.2f} MHz, duration {da:.0f} vs {db:.0f} us")


def _cue_identity(a, b):
    ia, ib = _ids(a), _ids(b)
    out = []
    if "id" in ia and "id" in ib:
        same = ia["id"] == ib["id"]
        out.append(dict(cue="decoded_id", llr=15.0 if same else -15.0,
                        detail=f"decoded identity {ia['id']} vs {ib['id']}"))
    for k in ("zc_root", "pilot_id", "sync_word", "bind_id"):
        if k in ia and k in ib and ia[k] != ib[k]:
            out.append(dict(cue=k, llr=-8.0, detail=f"{k} {ia[k]} vs {ib[k]} (radio configuration differs)"))
    return out


def same_emitter_score(a, b, models=None, prior_same=0.5, detail=True):
    """Log-likelihood ratio log p(obs | same radio) / p(obs | different radios of the
    same type) for two emissions, summed over (conditionally independent) physical
    cues, plus the evidence list. Returns (llr, evidence). Posterior P(same) =
    sigmoid(llr + logit(prior_same)) is included in the evidence."""
    ev = []
    if (a.family, a.label) != (b.family, b.label):
        ev.append(dict(cue="type", llr=-30.0,
                       detail=f"{a.family}/{a.label} vs {b.family}/{b.label}: different transmitter types"))
        llr = -30.0
    else:
        ma = _model(a, models)
        ev += _cue_identity(a, b)
        c = _cue_overlap(a, b, ma)
        if c: ev.append(c)
        osc = _cue_oscillator(a, b, ma)
        if osc: ev.append(osc)
        pa, _, _ = _ppm(a, ma); pb, _, _ = _ppm(b, ma)
        ppm_pair = None if pa is None or pb is None else 0.5 * (pa + pb)
        for c in (_cue_timing(a, b, ma, ppm_pair), _cue_clock(a, b, ma), _cue_power(a, b, ma), _cue_shape(a, b, ma)):
            if c: ev.append(c)
        llr = float(sum(c["llr"] for c in ev))
    lp = llr + math.log(prior_same / (1 - prior_same))
    post = 1 / (1 + math.exp(-max(min(lp, 700), -700)))
    ev.append(dict(cue="total", llr=llr, detail=f"sum LLR {llr:+.1f} nats -> P(same | prior {prior_same:g}) = {post:.4g}",
                   p_same=post))
    for c in ev:
        c["llr"] = round(float(c["llr"]), 2)
    return llr, ev


# =================================================================== grouping
@dataclass
class Emitter:
    emitter_id: str
    family: str
    label: str
    members: list                       # indices into the input emission list (time order)
    frequencies: list                   # distinct channel centres (MHz)
    hop_pattern: dict
    frame_period_ms: float | None
    fingerprint: dict
    confidence: float                   # P(members are one radio) (or P(distinct) for a singleton)
    evidence: list = field(default_factory=list)
    separation: list = field(default_factory=list)   # vs. other emitters of the same type
    linked: list = field(default_factory=list)       # other-type emitters of the same drone system

    def to_dict(self):
        return asdict(self)


def _sigmoid(v):
    return 1 / (1 + math.exp(-max(min(v, 700), -700)))


def _cluster(idx, M, prior_logit):
    """Average-linkage agglomeration on the pairwise LLR matrix: merge the pair of
    clusters whose mean cross LLR (+ prior) is highest and positive, never across a
    vetoed pair (e.g. simultaneous transmission on two channels)."""
    cl = [[i] for i in idx]
    while len(cl) > 1:
        best = (0.0, None)
        for i in range(len(cl)):
            for j in range(i + 1, len(cl)):
                sub = M[np.ix_(cl[i], cl[j])]
                if sub.min() <= LLR_VETO:
                    continue
                s = sub.mean() + prior_logit
                if s > best[0]:
                    best = (s, (i, j))
        if best[1] is None:
            break
        i, j = best[1]
        cl[i] = sorted(cl[i] + cl[j]); del cl[j]
    return cl


def _frame_fit(es, ma, ppm):
    """Joint lattice fit of one emitter's bursts: t_i = t0 + u_i * (1 - eps), u_i the
    nominal lattice position (slot grid). eps = frame-clock offset in ppm (same sign
    convention as carrier ppm). Returns dict or None."""
    if not ma["lattices"] or len(es) < 2:
        return None
    frame, offs, _ = ma["lattices"][0]
    grid = min(fr / len(o) if len(o) == 1 else fr for fr, o, _ in ma["lattices"])
    grid = min([fr for fr, o, _ in ma["lattices"] if len(o) == 1] or [frame])
    t = np.array([_t(e)[0] for e in es]); s = np.array([_t(e)[1] for e in es])
    f = np.array([e.center_mhz for e in es])
    e0 = (ppm or 0.0) * 1e-6
    u = np.round((t - t[0]) / (grid * (1 - e0))) * grid
    var = s ** 2 + ma["jitter_us"] ** 2 + (ma["chan_delay_us_per_mhz"] * (f - np.median(f))) ** 2
    w = 1 / var
    A = np.vstack([np.ones_like(u), u]).T
    if len(es) >= 3 and np.ptp(u) > 0:
        cov = np.linalg.inv(A.T @ (A * w[:, None]))
        p = cov @ (A.T @ (w * t))
        res = t - A @ p
        eps = -(p[1] - 1) * 1e6
        sig_eps = math.sqrt(cov[1, 1]) * 1e6
    else:
        p = np.array([t[0] - u[0] * (1 - e0), 1 - e0])
        res = t - A @ p
        eps, sig_eps = None, None
    # sub-pattern: positions within the strict frame
    sub = sorted({round(float(v % frame) / 1e3, 3) for v in u})
    t0 = float(p[0]) % frame
    return dict(t0_us=float(p[0]), frame_phase_us=round(t0, 3),
                frame_period_ms=frame * (1 - (eps if eps is not None else (ppm or 0.0)) * 1e-6) / 1e3,
                clock_ppm=None if eps is None else round(eps, 3),
                clock_ppm_sigma=None if sig_eps is None else round(sig_eps, 3),
                residual_rms_ns=round(float(np.sqrt(np.mean(res ** 2))) * 1e3, 1),
                residual_ns=[round(float(v) * 1e3, 1) for v in res],
                slot_offsets_in_frame_ms=sub, lattice_grid_ms=grid / 1e3, frame_ms=frame / 1e3)


def _summarise(eid, es, idx, M_sub, ma, prior_logit):
    es_sorted = sorted(zip(idx, es), key=lambda p: _t(p[1])[0])
    idx = [i for i, _ in es_sorted]; es = [e for _, e in es_sorted]
    t = np.array([_t(e)[0] for e in es])
    freqs = sorted({round((_carrier(e)[0] or e.center_mhz * 1e6) / 1e6, 3) for e in es})
    # merge near-identical centres (same channel, CFO-level differences)
    chans = []
    for f_ in freqs:
        if not chans or f_ - chans[-1] > 0.2:
            chans.append(f_)
    pp = [_ppm(e, ma) for e in es]
    pv = np.array([p[0] for p in pp if p[0] is not None])
    ps = np.array([max(p[1], 1e-4) for p in pp if p[0] is not None])
    fp = {}
    ppm = None
    if len(pv):
        w = 1 / (ps ** 2 + ma["ppm_h1"] ** 2)
        ppm = float(np.sum(w * pv) / np.sum(w))
        fp.update(carrier_ppm=round(ppm, 4), carrier_ppm_sigma=round(float(1 / np.sqrt(np.sum(w))), 4),
                  carrier_ppm_spread=round(float(np.sqrt(np.sum(w * (pv - ppm) ** 2) / np.sum(w))), 4),  # weighted
                  carrier_ppm_source=pp[0][2],
                  carrier_offset_hz_per_hop=[None if p[0] is None else round(p[0] * e.center_mhz, 0)
                                             for p, e in zip(pp, es)])
    ff = _frame_fit(es, ma, ppm)
    if ff:
        fp.update(timing_residual_rms_ns=ff["residual_rms_ns"], timing_residual_ns=ff["residual_ns"])
        if ff["clock_ppm"] is not None:
            fp.update(clock_ppm=ff["clock_ppm"], clock_ppm_sigma=ff["clock_ppm_sigma"])
            if ppm is not None:
                z = (ff["clock_ppm"] - ppm) / math.hypot(ff["clock_ppm_sigma"], 1.0)
                fp.update(clock_minus_carrier_ppm=round(ff["clock_ppm"] - ppm, 3), clock_carrier_z=round(z, 2),
                          clock_carrier_consistent=bool(abs(z) < 3))
    fp.update(timing_source=dict((s, sum(1 for e in es if (e.features or {}).get("timing_source", "energy_edge") == s))
                                 for s in {(e.features or {}).get("timing_source", "energy_edge") for e in es}),
              snr_db=[None if e.snr_db is None else round(float(e.snr_db), 1) for e in es],
              bw_mhz=round(float(np.median([_bw(e) for e in es])), 3),
              burst_us=round(float(np.median([e.end_us - e.start_us for e in es])), 1))
    iv = np.diff(t) / 1e3
    hop = dict(n_bursts=len(es), n_channels=len(chans),
               sequence=[(round(float(ti) / 1e3, 4), round(e.center_mhz, 3)) for ti, e in zip(t, es)],
               intervals_ms=[round(float(v), 4) for v in iv])
    if ff:
        hop.update(frame_ms=ff["frame_ms"], slot_offsets_in_frame_ms=ff["slot_offsets_in_frame_ms"],
                   frame_phase_us=ff["frame_phase_us"])
    # confidence + evidence
    ev = []
    if len(es) >= 2:
        coh = [float(np.mean(np.delete(M_sub[k], k))) for k in range(len(es))]
        conf = _sigmoid(min(coh) + prior_logit)
        ev.append(f"cohesion: weakest member's mean LLR to the rest {min(coh):+.1f} nats "
                  f"(mean pairwise {float(np.mean(M_sub[np.triu_indices(len(es), 1)])):+.1f})")
        if ff:
            ev.append(f"timing: {len(es)} burst starts on ONE frame clock "
                      f"({ff['frame_ms']:g} ms frame, slots at {ff['slot_offsets_in_frame_ms']} ms, "
                      f"phase {ff['frame_phase_us']/1e3:.4f} ms): residual rms {ff['residual_rms_ns']:.0f} ns "
                      f"over {np.ptp(t)/1e3:.1f} ms; an independent radio's start is uniform over "
                      f"~{min(fr / len(o) for fr, o, _ in ma['lattices'])/1e3:g} ms")
        if len(pv) >= 2:
            ev.append(f"oscillator: carrier offset {ppm:+.3f} ppm on all {len(pv)} hops "
                      f"(weighted spread {fp['carrier_ppm_spread']:.3f} ppm) across {chans[0]:.1f}-{chans[-1]:.1f} MHz: "
                      f"one crystal; two radios differ by ~{math.sqrt(2) * ma['ppm_pop']:.1f} ppm rms")
        if "clock_ppm" in fp and ppm is not None:
            ev.append(f"clock/carrier: frame clock {fp['clock_ppm']:+.2f} +- {fp['clock_ppm_sigma']:.2f} ppm vs "
                      f"carrier {ppm:+.3f} ppm ({'consistent' if fp['clock_carrier_consistent'] else 'INCONSISTENT'}"
                      f" with one reference crystal)")
    else:
        conf = None
        ev.append("single emission: identity rests on the separation evidence")
    return Emitter(eid, es[0].family, es[0].label, idx, chans, hop,
                   None if not ff else round(ff["frame_period_ms"], 7), fp,
                   None if conf is None else round(conf, 4), ev)


def associate(emissions, fc=None, models=None, prior_same=0.5, return_matrix=False):
    """Group emissions into physical emitters (one radio each, possibly hopping).
    Emissions of different (family, label) never merge; within a type, pairwise
    same_emitter_score LLRs drive average-linkage clustering. Each Emitter carries
    its hop pattern, frame clock, oscillator fingerprint, a confidence and the
    evidence; .separation says, for every other emitter of the same type, how
    strongly the physics says DIFFERENT radio (P(same) and the deciding cues)."""
    n = len(emissions)
    prior_logit = math.log(prior_same / (1 - prior_same))
    M = np.full((n, n), -30.0)
    EV = {}
    for i in range(n):
        M[i, i] = 0.0
        for j in range(i + 1, n):
            a, b = emissions[i], emissions[j]
            if (a.family, a.label) != (b.family, b.label):
                continue
            llr, ev = same_emitter_score(a, b, models, prior_same)
            M[i, j] = M[j, i] = llr
            EV[(i, j)] = ev
    classes = {}
    for i, e in enumerate(emissions):
        classes.setdefault((e.family, e.label), []).append(i)
    out = []
    for key, idx in classes.items():
        ma = _model(emissions[idx[0]], models)
        for cl in sorted(_cluster(idx, M, prior_logit), key=lambda c: min(_t(emissions[i])[0] for i in c)):
            em = _summarise(f"E{len(out) + 1}", [emissions[i] for i in cl], cl, M[np.ix_(cl, cl)], ma, prior_logit)
            out.append(em)
    # separation between emitters of the same type
    for a in out:
        for b in out:
            if a is b or (a.family, a.label) != (b.family, b.label):
                continue
            sub = M[np.ix_(a.members, b.members)]
            mean = float(sub.mean())
            cues = {}
            for i in a.members:
                for j in b.members:
                    for c in EV.get((min(i, j), max(i, j)), []):
                        if c["cue"] != "total":
                            cues.setdefault(c["cue"], []).append(c["llr"])
            a.separation.append(dict(other=b.emitter_id, mean_llr=round(mean, 2),
                                     p_same=round(_sigmoid(mean + prior_logit), 6),
                                     verdict="DIFFERENT radio" if mean + prior_logit < 0 else "ambiguous",
                                     cue_mean_llr={k: round(float(np.mean(v)), 2) for k, v in cues.items()}))
        if a.confidence is None:                    # singleton: confidence it is distinct
            ps = [s["p_same"] for s in a.separation]
            a.confidence = round(1 - max(ps), 4) if ps else None
    # same drone system: different-type emitters sharing the oscillator (an RC locks its
    # carrier to the drone's downlink, so uplink, video and DroneID show one ppm)
    for a in out:
        for b in out:
            if a is b or a.family.split("_")[0] != b.family.split("_")[0] or a.label == b.label:
                continue
            pa, pb = a.fingerprint.get("carrier_ppm"), b.fingerprint.get("carrier_ppm")
            if pa is None or pb is None:
                continue
            sa = math.hypot(a.fingerprint["carrier_ppm_sigma"], 0.1)
            sb = math.hypot(b.fingerprint["carrier_ppm_sigma"], 0.1)
            z = (pa - pb) / math.hypot(sa, sb)
            if abs(z) < 3:
                a.linked.append(dict(other=b.emitter_id, label=b.label, dppm=round(pa - pb, 4),
                                     note="same reference oscillator within error: likely the same "
                                          "drone system (RC uplink locks to the drone)"))
    if return_matrix:
        return out, M, EV
    return out


def report(emitters, emissions=None):
    lines = []
    for e in emitters:
        hp = e.hop_pattern
        lines.append(f"{e.emitter_id} {e.family}/{e.label}: {hp['n_bursts']} bursts on {hp['n_channels']} "
                     f"channel(s) {e.frequencies} MHz, confidence {e.confidence}")
        if e.frame_period_ms:
            lines.append(f"   frame {e.frame_period_ms:.6f} ms, slots {hp.get('slot_offsets_in_frame_ms')} ms, "
                         f"intervals {hp['intervals_ms']}")
        fp = e.fingerprint
        if "carrier_ppm" in fp:
            lines.append(f"   carrier {fp['carrier_ppm']:+.3f} +- {fp['carrier_ppm_sigma']:.3f} ppm "
                         f"(spread {fp['carrier_ppm_spread']:.3f})"
                         + (f", timing residual rms {fp['timing_residual_rms_ns']:.0f} ns" if "timing_residual_rms_ns" in fp else ""))
        for s in e.evidence:
            lines.append("   + " + s)
        for s in e.separation:
            lines.append(f"   vs {s['other']}: {s['verdict']} (mean LLR {s['mean_llr']:+.1f}, P(same) {s['p_same']:.2g}) "
                         f"cues {s['cue_mean_llr']}")
        for s in e.linked:
            lines.append(f"   linked to {s['other']} ({s['label']}): dppm {s['dppm']:+.3f} -- {s['note']}")
    return "\n".join(lines)


def analyse_capture(meta_path, log=print):
    """Capture -> 15 kHz-OFDM bursts (assoc_bursts) -> fingerprints -> emitters."""
    from . import bursts as AB
    from spectrum_intel.core.emission import load_capture
    x, fs, fc, meta = load_capture(meta_path)
    em, ch, spurs = AB.extract_bursts(x, fs, fc, return_ch=True)
    measure_fingerprints(x, fs, fc, em, ch=ch)
    emitters = associate(em, fc)
    return em, emitters


if __name__ == "__main__":
    import sys
    em, emitters = analyse_capture(sys.argv[1])
    for i, e in enumerate(em):
        f = e.features
        print(f"[{i}] {e.label:15s} t {f.get('t_start_us', e.start_us):12.4f} us +-{f.get('t_sigma_us', 0):.3f}  "
              f"carrier {f.get('carrier_hz', 0)/1e6:.6f} MHz +-{f.get('carrier_sigma_hz', 0):.0f} Hz  "
              f"({f.get('timing_source')})")
    print(report(emitters, em))
