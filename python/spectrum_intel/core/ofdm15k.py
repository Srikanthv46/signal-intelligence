"""15 kHz-SCS OFDM burst scan and hop tracker (LTE-like links: DJI RC uplink and
similar), numerology labels. From the research spectrum_survey module, code unchanged."""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from .spectral import _cp_metric, _folded_scores


# ============================================================== 15 kHz OFDM
# (name, output rate, FFT size N, normal CP at that rate, occupied BW, scan step)
PROFILES = [
    ("1.4MHz", 1.92e6, 128, 9, 1.08e6, 0.25e6),
    ("10MHz", 15.36e6, 1024, 72, 9.0e6, 1.0e6),
    ("20MHz", 30.72e6, 2048, 144, 18.0e6, 2.0e6),
]

def ofdm15k_scan(ch, fc, profiles=PROFILES, win_s=0.5e-3, k_sigma=7.0, usable_half=28e6):
    """Returns raw hits: dict(profile, f_off, t0_s, t1_s, score) per connected
    region of the (centre, time-window) score grid above threshold, plus the
    per-profile score statistics (median, robust sigma) for later passes."""
    hits, stats = [], {}
    for name, ofs, N, cp, occ, step in profiles:
        sym = N + cp
        W = int(win_s * ofs); hop = W // 2
        offs = np.arange(-usable_half + occ / 2, usable_half - occ / 2 + 1, step)
        grid = []
        for off in offs:
            y = ch.slice(off, occ / 2 * 1.05, ofs)
            m, _ = _cp_metric(y, N, cp)
            grid.append(_folded_scores(m, sym, W, hop))
        G = np.array(grid)
        med = np.median(G); mad = np.median(np.abs(G - med)) * 1.4826
        thr = med + k_sigma * mad
        lab, n = ndimage.label(G > thr, np.ones((3, 3)))
        for i, sl in enumerate(ndimage.find_objects(lab)):
            sub = np.where(lab[sl] == i + 1, G[sl], -np.inf)
            ci, wi = np.unravel_index(np.argmax(sub), sub.shape)
            ci += sl[0].start; wi += sl[1].start
            hits.append(dict(profile=name, ofs=ofs, N=N, cp=cp, occ=occ, f_off=float(offs[ci]),
                             t0_s=wi * hop / ofs, t1_s=(wi * hop + W) / ofs,
                             score=float(G[ci, wi]), thr=float(thr)))
        stats[name] = dict(median=float(med), sigma=float(mad))
    return hits, stats

def refine_ofdm15k(ch, hit, pad_s=1.0e-3):
    """Burst edges from the envelope, centre from the CP phase (fine CFO),
    -10 dB occupied bandwidth from a min-hold PSD (rejects overlapping bursty Wi-Fi)."""
    ofs, N, cp, occ = hit["ofs"], hit["N"], hit["cp"], hit["occ"]
    y = ch.slice(hit["f_off"], occ / 2 * 1.05, ofs)
    a = max(0, int((hit["t0_s"] - pad_s) * ofs)); b = min(len(y), int((hit["t1_s"] + pad_s) * ofs))
    seg = y[a:b]
    k = max(int(20e-6 * ofs), 4)
    env = np.convolve(np.abs(seg) ** 2, np.ones(k) / k, "same")
    pk = int(np.clip(((hit["t0_s"] + hit["t1_s"]) / 2) * ofs - a, 0, len(env) - 1))
    hi, lo = np.percentile(env, 90), np.percentile(env, 10)
    on = env > np.sqrt(hi * lo)
    s = pk
    while s > 0 and on[s - 1]: s -= 1
    e = pk
    while e < len(on) - 1 and on[e + 1]: e += 1
    s0, e0 = a + s, a + e
    if e0 - s0 < 2 * (N + cp):                            # shorter than two symbols: not an OFDM burst
        return None
    # occupied band (-10 dB from the in-band level) at 3x the profile width: gives the
    # bandwidth AND the coarse centre (the scan grid is up to step/2 off; the CP phase
    # below only resolves +-SCS/2 = +-7.5 kHz)
    wide = ch.slice(hit["f_off"], occ * 1.5, ofs * 2)
    ws, we = s0 * 2, e0 * 2
    nf = 4 * N
    frames = [np.abs(np.fft.fftshift(np.fft.fft(wide[j:j + nf] * np.hanning(nf)))) ** 2
              for j in range(ws, max(ws + 1, we - nf), nf // 2)]
    bw, shift = occ, 0.0
    if frames:
        ps = 10 * np.log10(np.maximum(ndimage.uniform_filter1d(np.median(np.array(frames), 0), 9), 1e-20))
        ff = (np.arange(nf) - nf // 2) / nf * ofs * 2
        inb = np.abs(ff) < occ * 0.3
        occ_mask = ps > np.median(ps[inb]) - 10
        L_ = R_ = int(np.argmax(np.where(inb, ps, -np.inf)))       # strongest in-band bin

        while L_ > 0 and occ_mask[L_ - 1]: L_ -= 1
        while R_ < nf - 1 and occ_mask[R_ + 1]: R_ += 1
        bw, shift = float(ff[R_] - ff[L_]), float((ff[R_] + ff[L_]) / 2)
    f_off = hit["f_off"] + shift
    # fine CFO from the CP correlation phase over the burst, at the corrected centre
    y = ch.slice(f_off, occ / 2 * 1.05, ofs)
    m, c = _cp_metric(y[s0:e0], N, cp)
    strong = m > np.percentile(m, 90)
    ph = np.angle(np.sum(np.convolve(c, np.ones(cp), "valid")[strong]))
    cfo = -ph / (2 * np.pi * N / ofs)
    f_off += cfo
    if bw > 1.5 * occ:                                    # band edges not found: not a burst of this width
        return None
    # verification: a genuine 15 kHz OFDM burst shows a TRAIN of CP peaks one symbol
    # period apart. The folded score alone can be pushed up by edges of strong bursts.
    from scipy.signal import find_peaks
    # Fit the peaks to a uniform grid: the fitted period must be the symbol period. Per-peak
    # scatter is tolerated (a residual tone, e.g. LO leakage, flattens the CP plateau and the
    # peak wanders inside it); random-spaced artefact peaks do not fit one period.
    sym = N + cp
    pk, _ = find_peaks(m, height=0.5, distance=int(0.8 * sym))
    if len(pk) < 3:
        return None
    k = np.round((pk - pk[0]) / sym)
    if len(np.unique(k)) < 3:
        return None
    period, a0 = np.polyfit(k, pk, 1)
    scatter = float(np.sqrt(np.mean((pk - (a0 + period * k)) ** 2)))
    if abs(period - sym) > 0.03 * sym or scatter > 0.12 * sym:
        return None
    nsym = int(round((e0 - s0) / (N + cp)))
    return dict(start_s=s0 / ofs, end_s=e0 / ofs, f_off=f_off, cfo_hz=cfo, bw_hz=bw, nsym=nsym,
                score=hit["score"], profile=hit["profile"], occ=occ, cp_peaks=int(len(pk)),
                period_us=round(float(period / ofs * 1e6), 2), cp_scatter_us=round(scatter / ofs * 1e6, 2))

# ============================================================== blind numerology
# useful-symbol duration -> what it is. Wi-Fi: 64-pt FFT at 20 MSps (+ 0.8/0.4 us GI);
# 802.11ax data: 256-pt at 20 MSps; LTE / DJI OcuSync: 15 kHz SCS.
NUMEROLOGIES = [(3.2e-6, "wifi_like_ofdm", "312.5 kHz SCS (802.11a/g/n/ac symbol)"),
                (12.8e-6, "wifi_he_like_ofdm", "78.125 kHz SCS (802.11ax data symbol)"),
                (66.67e-6, "lte_like_ofdm", "15 kHz SCS (LTE / OcuSync)")]

def ofdm_numerology(ch, f_off, half_bw, t0_s, t1_s, max_tu=80e-6):
    """Blind useful-symbol length of whatever occupies (f, t): the lag where the
    normalised autocorrelation of the segment peaks (the CP repeats the symbol
    tail Tu later). Returns (Tu_s, peak, baseline) or None if no clear peak."""
    ofs = max(4 * half_bw, 2.5e6)
    y = ch.slice(f_off, half_bw, ofs)
    a, b = int(t0_s * ofs), int(t1_s * ofs)
    z = y[a:b]
    E = np.sum(np.abs(z) ** 2) + 1e-20
    lags = np.arange(max(2, int(2.5e-6 * ofs)), min(int(max_tu * ofs), len(z) // 2))   # below ~2.5 us the
    # band-limited autocorrelation decay dominates; shortest symbol of interest is Wi-Fi 3.2 us
    if lags.size < 8:
        return None
    r = np.array([np.abs(np.vdot(z[:-L], z[L:])) / E for L in lags])
    # an OFDM peak is narrow (CP-length); remove the slowly varying bandlimited-noise trend
    base = ndimage.median_filter(r, max(5, int(2e-6 * ofs)) | 1)
    d = r - base
    i = int(np.argmax(d))
    noise = 1.4826 * np.median(np.abs(d - np.median(d))) + 1e-9
    if d[i] < 6 * noise or d[i] < 0.04:
        return None
    return float(lags[i] / ofs), float(r[i]), float(base[i])

def classify_numerology(tu):
    for ref, lab, why in NUMEROLOGIES:
        if abs(tu - ref) < 0.04 * ref + 0.15e-6:
            return lab, why
    return "ofdm_unknown", f"OFDM, Tu {tu*1e6:.2f} us -> SCS {1/tu/1e3:.1f} kHz (no match in table)"

def label_ofdm15k(b):
    dur_us = (b["end_s"] - b["start_s"]) * 1e6
    bw = b["bw_hz"]
    if bw < 2.5e6 and 300 <= dur_us <= 1200:
        return "ocusync_uplink_like", "LTE-like 15 kHz OFDM, ~1.4 MHz, slot-length burst (RC->drone hop)"
    if 7e6 <= bw <= 11e6 and 500 <= dur_us <= 800:
        return "droneid_candidate", "LTE-like 15 kHz OFDM, ~10 MHz, ~9-symbol burst (DJI DroneID shape)"
    if bw >= 7e6:
        return "ocusync_video_like", "LTE-like 15 kHz OFDM, wideband (drone->RC video downlink shape)"
    return "ofdm15k_other", "15 kHz-SCS OFDM, unrecognised shape"

def merge_ofdm15k(bursts):
    """A narrow burst also lights wider profiles whose window contains it, and
    neighbouring scan centres see the same burst: keep one record per burst,
    preferring the narrowest profile."""
    order = {p[0]: i for i, p in enumerate(PROFILES)}
    keep = []
    for b in sorted(bursts, key=lambda b: (order[b["profile"]], -b["score"])):
        dup = False
        for k in keep:
            ov = min(b["end_s"], k["end_s"]) - max(b["start_s"], k["start_s"])
            if ov <= 0.5 * min(b["end_s"] - b["start_s"], k["end_s"] - k["start_s"]):
                continue
            if abs(b["f_off"] - k["f_off"]) <= (b["occ"] + k["occ"]) / 2:
                dup = True; break
        if not dup:
            keep.append(b)
    return keep

def _cycle(starts_s):
    """Repeat period of a burst train: constant spacing, or an alternating
    pair (e.g. 4 + 6 ms) whose sum is constant. None if neither."""
    iv = np.diff(np.sort(starts_s))
    if len(iv) >= 3:
        pair = iv[:-1] + iv[1:]
        if np.std(pair) < 0.15 * np.mean(pair) and np.std(iv) > 0.1 * np.mean(iv):
            return float(np.mean(pair))
    if len(iv) >= 2 and np.std(iv) < 0.1 * np.mean(iv):
        return float(np.mean(iv))
    return None

def track_hops(ch, fc, bursts, stats, dur_s, k_sigma=4.0, usable_half=28e6):
    """Hop-train confirm pass. Once >= 4 narrow 15 kHz bursts establish a cycle,
    re-test only the times the cycle predicts, across every centre, at a lower
    threshold. Recovers hops the blind pass missed (e.g. next to the LO spur)
    without lowering the blind threshold everywhere."""
    name, ofs, N, cp, occ, step = PROFILES[0]
    narrow = [b for b in bursts if b["profile"] == name]
    if len(narrow) < 4:
        return []
    starts = np.array([b["start_s"] for b in narrow])
    T = _cycle(starts)
    if T is None:
        return []
    cand = sorted({round(t + k * T, 5) for t in starts for k in range(-4, 5) if k})
    cand = [t for t in cand if 0 <= t <= dur_s - 0.5e-3 and np.min(np.abs(starts - t)) > 1e-3]
    merged = []
    for t in cand:                                          # predictions from different anchors coincide
        if not merged or t - merged[-1] > 0.3e-3:
            merged.append(t)
    thr = stats[name]["median"] + k_sigma * stats[name]["sigma"]
    sym, W = N + cp, int(0.5e-3 * ofs)
    offs = np.arange(-usable_half + occ / 2, usable_half - occ / 2 + 1, step / 2)   # few times: afford a finer grid
    # a lower threshold over ~10^3 (centre, offset) trials per prediction needs a second,
    # independent test: the candidate must look like the train (narrow, same burst length)
    dur0 = float(np.median([b["end_s"] - b["start_s"] for b in narrow]))
    slices = [(off, ch.slice(off, occ / 2 * 1.05, ofs)) for off in offs]
    found = []
    for t in merged:
        cands = []
        for off, y in slices:
            best = (-1.0, None)
            for dt in np.arange(-0.25e-3, 0.11e-3, 0.05e-3):       # hop timing jitters by ~0.1-0.2 ms
                a = int((t + dt) * ofs)
                if a < 0 or a + W + N + cp >= len(y):
                    continue
                m, _ = _cp_metric(y[a:a + W + N + cp], N, cp)
                sc = _folded_scores(m, sym, W, W)
                if sc.size and sc[0] > best[0]:
                    best = (float(sc[0]), t + dt)
            if best[0] > thr:
                cands.append((best[0], off, best[1]))
        for score, off, t0 in sorted(cands, reverse=True)[:8]:
            hit = dict(profile=name, ofs=ofs, N=N, cp=cp, occ=occ, f_off=off,
                       t0_s=t0, t1_s=t0 + 0.5e-3, score=score, thr=thr)
            b = refine_ofdm15k(ch, hit)
            if b is None or b["bw_hz"] > 2.5e6 or abs((b["end_s"] - b["start_s"]) - dur0) > 0.25 * dur0:
                continue
            b["tracked"] = True
            found.append(b)
            break
    return found

def hop_summary(bursts, fc):
    """Frequency/time pattern of a hopping link."""
    if not bursts:
        return {}
    bs = sorted(bursts, key=lambda b: b.start_us)
    t = np.array([b.start_us for b in bs]) / 1e3
    iv = np.diff(t)
    out = dict(n_bursts=len(bs), centres_mhz=sorted({round(b.center_mhz, 2) for b in bs}),
               sequence=[(round(b.start_us / 1e3, 2), round(b.center_mhz, 2)) for b in bs],
               intervals_ms=[round(v, 2) for v in iv],
               burst_us=round(float(np.median([b.end_us - b.start_us for b in bs])), 0),
               bw_mhz=round(float(np.median([b.f_hi_mhz - b.f_lo_mhz for b in bs])), 2))
    T = _cycle(t)
    if T is not None:
        out["cycle_ms"] = round(T, 2)
        out["bursts_per_cycle_in_band"] = int(round(T / np.mean(iv)))
    out["tracked_bursts"] = sum(1 for b in bs if b.info.get("tracked"))
    out["note"] = ("only hops landing inside this capture's band are seen; the full hop set "
                   "is wider than one 56 MHz capture")
    return out
