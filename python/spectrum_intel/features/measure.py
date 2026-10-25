"""Physical-layer features of one detected signal, measured on its own IQ.

Each feature records how it was measured, so the classifier evidence and the
report can cite it.
"""
from __future__ import annotations

import numpy as np

from ..core import bands_india as B
from ..core.deinterleave import deinterleave
from ..core.dsp import Channelizer, burst_mask
from ..model import Signal

MAX_SAMPLES = 1 << 19


def timing(sig: Signal, capture_s: float) -> dict:
    b = sorted(sig.bursts, key=lambda x: x.t0)
    on = sum(x.t1 - x.t0 for x in b)
    out = {"n_bursts": len(b), "duty_cycle": round(min(1.0, on / capture_s), 3),
           "burst_dur_s_median": float(np.median([x.t1 - x.t0 for x in b])),
           "burst_period_s": None, "continuous": any(x.kind == "continuous" for x in b)}
    if len(b) >= 3:
        d = np.diff([x.t0 for x in b])
        med = float(np.median(d))
        if med > 0 and np.std(d) / med < 0.15:
            out["burst_period_s"] = round(med, 6)
    return out


PULSE_THR_DB = 10.0        # 1 us envelope above the band's median level
PULSE_GAP_S = 20e-6        # gaps shorter than this join one pulse group (DME pair, radar short + long pulse)
PRI_MIN_S = 100e-6         # group-to-group intervals below this are fragments of one group
PRI_TOL = 0.15             # intervals within +-15 % of the shortest cluster are one PRI (stagger included)
PRI_MIN_FRAC = 0.5         # ... and that cluster holds at least half of the intervals
LFM_MIN_S = 10e-6          # chirp test on pulses at least this long
LFM_R2 = 0.8               # linear fit of the instantaneous frequency (1 us blocks)
OFDM_MIN_LAGS = 4          # a CP-OFDM autocorrelation lag (one useful symbol) fits >= 4 times into its bursts
DEINT_MIN = 6              # pulses before the deinterleaver runs (core.deinterleave)
DEINT_MAX = 2000           # ... and the most it takes (the first ones)
PULSE_MIN_BW = 100e3       # a pulse of width tau occupies >~ 1/tau: below 100 kHz no pulse < 10 us fits the band


def pulse_train(sig: Signal, ch: Channelizer) -> dict:
    """Pulse groups of a bursty signal over the whole capture, on its own band at >= 1 MSps with a 1 us
    envelope: the spectrogram's burst boxes (16.7 us frames, CFAR growth) are far longer than radar/DME pulses.
    Returns pulse_n, pulse_width_s_median, pulse_duty, and when the groups repeat pri_s / pri_spread; for the
    widest pulse a linear instantaneous-frequency fit (LFM chirp: lfm_sweep_hz over lfm_pulse_s)."""
    if sig.bw < PULSE_MIN_BW:
        return {}
    y, ofs = ch.slice(sig.fc, 0.5 * sig.bw, min(ch.fs, max(1.25 * sig.bw, 1e6)))
    m = max(1, int(round(ofs * 1e-6)))
    k = len(y) // m
    if k < 64:
        return {}
    e = (np.abs(y[:k * m]) ** 2).reshape(k, m).mean(1)
    idx = np.flatnonzero(e > np.median(e) * 10 ** (PULSE_THR_DB / 10))
    if len(idx) == 0:
        return {}
    dt = m / ofs
    br = np.flatnonzero(np.diff(idx) * dt > PULSE_GAP_S)
    st, en = idx[np.r_[0, br + 1]], idx[np.r_[br, len(idx) - 1]]
    w = (en - st + 1) * dt
    out = {"pulse_n": int(len(st)), "pulse_width_s_median": round(float(np.median(w)), 7),
           "pulse_duty": round(float(w.sum() / (k * dt)), 4),
           "_pulse_t_s": [float(v) for v in st * dt]}          # pulse starts (detect.merge_pulse_trains)
    if len(st) >= DEINT_MIN:                                     # emitters by PRI (research, signal gaps S10)
        de = deinterleave([float(v) for v in st[:DEINT_MAX] * dt], [float(v) for v in w[:DEINT_MAX]])
        out["pulse_emitters"] = de["emitters"]
        out["pulse_unassigned"] = de["unassigned"]
    d = np.diff(st) * dt
    d = d[d >= PRI_MIN_S]
    if len(d) >= 2:
        # weak pulses below the threshold leave 2 x PRI gaps: the PRI is the shortest interval cluster
        p0 = float(np.percentile(d, 25))
        c = d[np.abs(d - p0) <= PRI_TOL * p0]
        if len(c) >= max(2, PRI_MIN_FRAC * len(d)):            # random intervals also leave a few in the window
            pri = float(np.median(c))
            out["pri_s"] = round(pri, 7)
            out["pri_spread"] = round(float((c.max() - c.min()) / pri), 3)
    j = int(np.argmax(w))
    if w[j] >= LFM_MIN_S:
        a, b = st[j] * m, (en[j] + 1) * m
        a, b = a + (b - a) // 10, b - (b - a) // 10                 # central 80 %: edges ride on the rise/fall
        # instantaneous frequency per 1 us block (phase step summed over the block: per-sample phase noise
        # at ~15 dB SNR leaves R^2 ~0.75 on a clean chirp, 1 us blocks ~0.9)
        q = y[a + 1:b] * np.conj(y[a:b - 1])
        nb = len(q) // m
        if nb >= 8:
            fi = np.angle(q[:nb * m].reshape(nb, m).sum(1)) * ofs / (2 * np.pi)
            t = (np.arange(nb) + 0.5) * m / ofs
            tc, fc = t - t.mean(), fi - fi.mean()                   # least-squares line, closed form
            stt, sff, stf = float(np.sum(tc * tc)), float(np.sum(fc * fc)) + 1e-30, float(np.sum(tc * fc))
            r2 = stf * stf / (stt * sff)
            sweep = stf / stt * (b - a) / ofs
            if r2 >= LFM_R2 and abs(sweep) > 2 * ofs / (b - a):
                out["lfm_sweep_hz"] = round(sweep)
                out["lfm_pulse_s"] = round(float(w[j]), 7)
                out["lfm_r2"] = round(r2, 3)
    return out


def _peaks(spec, freqs, fmin, fmax, k=3):
    s = (freqs >= fmin) & (freqs <= fmax)
    if s.sum() < 8:
        return []
    v, f = spec[s], freqs[s]
    med = np.median(v) + 1e-30
    idx = np.argsort(v)[::-1][:k]
    return [(float(f[i]), float(v[i] / med)) for i in idx]


def iq_features(sig: Signal, ch: Channelizer, n: int) -> tuple[dict, list[str]]:
    """Envelope, cyclostationary lines, OFDM CP lag, FSK deviation."""
    notes = []
    bw = max(sig.bw, 5e3)
    ofs_req = min(ch.fs, max(1.25 * bw, 20e3))
    y, ofs = ch.slice(sig.fc, 0.55 * bw, ofs_req)
    m = burst_mask(n, ch.fs, sig.bursts)
    idx = (np.arange(len(y)) * ch.fs / ofs).astype(int)
    y = y[m[np.clip(idx, 0, n - 1)]]
    if len(y) > MAX_SAMPLES:
        y = y[:MAX_SAMPLES]
    out = {"analysis_rate_hz": round(ofs), "analysis_samples": int(len(y))}
    if len(y) < 64:
        return out, ["too few on-burst samples for IQ features"]
    a = np.abs(y)
    out["envelope_cv"] = round(float(a.std() / (a.mean() + 1e-12)), 3)
    out["envelope_kurtosis"] = round(float(np.mean(a ** 4) / (np.mean(a ** 2) ** 2 + 1e-30)), 2)
    notes.append("envelope CV/kurtosis on on-burst samples (Gaussian noise/OFDM: CV~0.52, kurt~2; "
                 "constant envelope: CV->0, kurt->1)")

    # cyclic lines of |y|^2 -> symbol / chip / frame clocks of linear modulations
    e = a ** 2 - np.mean(a ** 2)
    E = np.abs(np.fft.rfft(e))
    fr = np.fft.rfftfreq(len(e), 1 / ofs)
    pk = _peaks(E, fr, max(2e3, 4 * ofs / len(e)), 0.98 * ofs / 2)
    if pk and pk[0][1] > 25:
        out["cyclic_line_hz"] = round(pk[0][0], 1)
        out["cyclic_line_ratio"] = round(pk[0][1], 1)
    # carrier lines of y^2 / y^4 (BPSK / QPSK signature)
    for p in (2, 4):
        Z = np.abs(np.fft.fft(y ** p))
        k = int(np.argmax(Z))
        r = float(Z[k] / (np.median(Z) + 1e-30))
        if r > 60:
            out[f"x{p}_line_ratio"] = round(r, 1)
            out[f"x{p}_line_hz"] = round(float(np.fft.fftfreq(len(y), 1 / ofs)[k]), 1)
    notes.append("cyclic lines: FFT of |x|^2 (symbol clock), x^2 and x^4 (BPSK/QPSK carrier lines); "
                 "peak/median > 25 (|x|^2) or > 60 (x^p) reported")

    # OFDM: cyclic-prefix autocorrelation at lag Tu, searched beyond 3/BW
    N = len(y)
    R = np.fft.ifft(np.abs(np.fft.fft(y, 2 * N)) ** 2)[:N]
    R = np.abs(R) / (np.abs(R[0]) + 1e-30)
    lags = np.arange(N) / ofs
    s = (lags > max(3 / bw, 2e-6)) & (lags < min(4e-3, N / ofs / 4))
    if s.sum() > 16:
        k = int(np.argmax(R[s]))
        r, lag = float(R[s][k]), float(lags[s][k])
        noise = float(np.median(R[s]) + 1e-12)
        out["acf_peak_lag_s"] = round(lag, 9)
        out["acf_peak"] = round(r, 4)
        out["acf_peak_over_median"] = round(r / noise, 1)
        # a real CP-OFDM signal: strong, isolated peak and >= 12 subcarriers inside the band
        if r / noise > 15 and r > 0.06 and bw * lag >= 12 and bw >= 150e3:
            out["ofdm_scs_hz"] = round(1 / lag, 1)
    notes.append("OFDM test: normalised autocorrelation peak at lag > 3/BW; a CP-OFDM signal "
                 "correlates at its useful symbol time Tu (SCS = 1/Tu)")

    # FSK: instantaneous-frequency clusters (only meaningful for constant envelope)
    if out["envelope_cv"] < 0.35:
        ifr = np.angle(y[1:] * np.conj(y[:-1])) * ofs / (2 * np.pi)
        lo, hi = np.percentile(ifr, [10, 90])
        out["inst_freq_p10_p90_hz"] = [round(float(lo)), round(float(hi))]
        h, _ = np.histogram(ifr, bins=64, range=(lo - 0.2 * (hi - lo), hi + 0.2 * (hi - lo)))
        out["inst_freq_bimodal"] = bool(h[:32].max() > 1.5 * h[28:36].min() and
                                        h[32:].max() > 1.5 * h[28:36].min())
        notes.append("instantaneous-frequency 10/90th percentiles and bimodality (2-FSK test)")
    return out, notes


def infer_modulation(sig: Signal, f: dict, df_fine: float) -> tuple[str, list[str]]:
    why = []
    if sig.bw <= 3 * df_fine + 1:
        return "CW (unmodulated carrier)", [f"occupied bandwidth {sig.bw:.0f} Hz <= 3 fine bins"]
    if sig.snr_db < 3:
        return "unknown (SNR too low)", [f"in-band SNR {sig.snr_db} dB: envelope/cyclic/CP statistics not reliable"]
    if f.get("lfm_sweep_hz"):
        why.append(f"pulses with a linear frequency sweep ({f['lfm_sweep_hz'] / 1e3:+.0f} kHz over "
                   f"{f['lfm_pulse_s'] * 1e6:.0f} us, fit R^2 {f['lfm_r2']}) -> LFM chirp (pulse compression)")
        return "LFM chirp pulses", why
    pw = f.get("pulse_width_s_median") if f.get("pulse_n", 0) >= 2 else None
    if f.get("ofdm_scs_hz") and pw and pw < OFDM_MIN_LAGS * f["acf_peak_lag_s"]:
        why.append(f"autocorrelation peak at {f['acf_peak_lag_s'] * 1e6:.2f} us but pulses only {pw * 1e6:.1f} us "
                   "long -> not OFDM (pulse train)")
    elif f.get("ofdm_scs_hz"):
        why.append(f"autocorrelation peak {f['acf_peak']} ({f['acf_peak_over_median']}x median) "
                   f"at {f['acf_peak_lag_s'] * 1e6:.2f} us -> OFDM, SCS {f['ofdm_scs_hz'] / 1e3:.2f} kHz")
        return "OFDM", why
    cl = f.get("cyclic_line_hz")
    if cl and abs(cl - 1e6) < 2e3 and sig.bw > 8e6:
        why.append(f"|x|^2 cyclic line at {cl / 1e6:.4f} MHz (x{f['cyclic_line_ratio']}) -> 1 Msym/s "
                   "Barker-spread DSSS clock")
        return "DSSS (DBPSK/DQPSK, Barker-11)", why
    cv = f.get("envelope_cv")
    if cv is not None and cv < 0.25 and sig.snr_db > 10:
        if f.get("inst_freq_bimodal"):
            why.append(f"constant envelope (CV {cv}) and bimodal instantaneous frequency -> 2-FSK/GFSK")
            return "GFSK/2-FSK", why
        why.append(f"constant envelope (CV {cv}) -> FM / FSK / GMSK family")
        return "constant-envelope (FM/FSK/GMSK)", why
    if f.get("x2_line_ratio") and sig.snr_db > 10:
        why.append(f"x^2 carrier line (x{f['x2_line_ratio']}) -> BPSK-like")
        return "BPSK-like", why
    if f.get("x4_line_ratio") and sig.snr_db > 10:
        why.append(f"x^4 carrier line (x{f['x4_line_ratio']}) -> QPSK-like")
        return "QPSK-like", why
    if cl:
        why.append(f"|x|^2 cyclic line at {cl / 1e3:.3f} kHz -> single-carrier linear modulation, "
                   f"symbol rate ~{cl / 1e3:.1f} ksym/s (or a multiple)")
        return "single-carrier linear (PSK/QAM)", why
    return "unknown", ["no decisive envelope, cyclostationary or OFDM signature"]


def extract(sig: Signal, ch: Channelizer, n: int, capture_s: float, df_fine: float):
    f = timing(sig, capture_s)
    iq, notes = iq_features(sig, ch, n)
    f.update(iq)
    if not f["continuous"] and B.pulsed_band(sig.f0, sig.f1, sig.fc):
        f.update(pulse_train(sig, ch))
        notes.append("pulse train: 1 us envelope of the signal's band over the whole capture, pulses > median + "
                     f"{PULSE_THR_DB:g} dB, gaps < {PULSE_GAP_S * 1e6:g} us joined; PRI = shortest interval "
                     "cluster holding >= half of the intervals; LFM = linear instantaneous-frequency fit on the "
                     "widest pulse; pulse_emitters: core.deinterleave (PW groups, CDIF + sequence search, stagger, "
                     "jitter mode; research)")
    if f.get("cyclic_line_hz"):
        f["symbol_rate"] = f["cyclic_line_hz"]
    mod, why = infer_modulation(sig, f, df_fine)
    sig.features.update(f)
    sig.modulation = mod
    sig.features.setdefault("symbol_rate", None)
    sig.detection["feature_methods"] = notes
    sig.classify_trace += why
