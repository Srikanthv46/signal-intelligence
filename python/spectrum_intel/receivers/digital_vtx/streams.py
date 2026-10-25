"""
streams.py -- how many transmit streams share each subcarrier of a blind CP-OFDM link (one receive channel).

Evidence, not identification: the result is a "-like" verdict with its measurements (bench/tx_antennas/FINDINGS.md).

On every run of >= MIN_SEG contiguous full-band symbols (symbol grid folded from the CP metric, fractional CFO
from the CP phase; training / narrowband symbols dropped by their occupancy), the raw FFT bins of the band give
  kurtosis   per-subcarrier E|X|^4 / (E|X|^2)^2 over symbols, after per-symbol power normalisation (a static
             one-tap channel cancels): one stream per subcarrier (1 TX, CDD, antenna switching) <= 1.40
             (QPSK 1.00 .. 256-QAM 1.40); two streams (STBC, spatial multiplexing) ~1.45-1.65; Gaussian-like ~1.9;
  pair delta the Alamouti fingerprint. In a pair |y1|^2 + |y2|^2 = (|h1|^2 + |h2|^2)(|a|^2 + |b|^2), so the
             pair's powers are anti-correlated. The statistic is the power correlation of symbol pairs (2i, 2i+1)
             minus that of (2i+1, 2i+2), counted from the run start. 802.11n STBC (Microhard, real):
             -0.33; synthetic: -0.16 (2nd branch -15 dB) .. -0.44; one stream / spatial multiplexing: 0.
             The same over subcarrier pairs detects SFBC.
Verdicts: 2tx_stbc_like / 2tx_sfbc_like (pair delta beyond PAIR_MIN on every run that the capture does not cut
off, with the same parity), 1_stream_like (kurtosis <= K1_MAX and no pairing), unknown otherwise (with the
reason). Gated on link SNR >= SNR_MIN and >= MIN_SC subcarriers.
"""
from __future__ import annotations

import numpy as np

SNR_MIN = 20.0          # dB: below, noise lifts the kurtosis of a single stream toward 2
MIN_SEG = 20            # contiguous full-band symbols per run
MIN_SC = 200            # subcarriers: fewer leave the pair statistic too noisy (null spread ~0.05)
OCC_MIN = 0.75          # fraction of the band's subcarriers above 0.1 x their mean power (data ~0.9;
                        # comb-2 training 0.5, a narrowband burst on the grid < 0.5)
CP_MIN = 0.5            # per-symbol normalised CP correlation
BAND_FRAC = 0.9         # use the central 90 % of the measured band
K1_MAX = 1.43           # one stream per subcarrier: 256-QAM 1.40 at 20 dB, 1.42 at 15 dB
PAIR_MIN = 0.10         # |pair delta|: noise <= 0.03 at 1700 subcarriers x 50 symbols, ~0.05 at 200 x 20
PAIR_NULL = 0.07        # below: no pairing


def _movsum(v, k):
    c = np.concatenate(([0], np.cumsum(v)))
    return c[k:] - c[:-k]


def _kurtosis(Z):
    p = np.abs(Z) ** 2
    p = p / np.mean(p, 1, keepdims=True)
    k = np.mean(p * p, 0) / np.mean(p, 0) ** 2
    return float(np.median(k))


def _rho(a, b, axis):
    a = a - a.mean(axis, keepdims=True)
    b = b - b.mean(axis, keepdims=True)
    den = np.sqrt(np.sum(a * a, axis) * np.sum(b * b, axis)) + 1e-30
    return float(np.median(np.sum(a * b, axis) / den))


def _pairs(Z):
    """(time delta, freq delta). Time: power correlation of symbols (2i, 2i+1) minus (2i+1, 2i+2), counted from the
    run start (signed: the paired parity). Freq: the same over subcarrier pairs, per block of subcarriers between
    empty ones (a null DC restarts the pairs), as -|delta| of the weakest block (pairing must hold in every block)."""
    p = np.abs(Z) ** 2
    p = p / np.mean(p, 1, keepdims=True)
    q = p / np.mean(p, 0, keepdims=True)            # time pairs: the channel's |H(k)|^2 removed first
    m = p.shape[0] // 2 * 2
    dt = _rho(q[0:m:2], q[1:m:2], 1) - _rho(q[1:m - 1:2], q[2:m:2], 1)
    mp = np.mean(p, 0)
    live = mp > 0.1 * np.median(mp)
    edges = np.flatnonzero(np.diff(np.concatenate(([0], live.astype(int), [0]))))
    blocks = [(a, b) for a, b in zip(edges[0::2], edges[1::2]) if b - a >= 64]
    df = []
    for a, b in blocks:
        n = (b - a) // 2 * 2
        pb = p[:, a:a + n]
        df.append(abs(_rho(pb[:, 0:n:2], pb[:, 1:n:2], 0) - _rho(pb[:, 1:n - 1:2], pb[:, 2:n:2], 0)))
    return dt, (-min(df) if df else 0.0)


def segments(y, N, ncp, Lf, a, b):
    """Symbol FFTs of the contiguous full-band runs inside y[a:b] (y: the link at N / Tu samples per second,
    centred). Returns [(first sample, Z (symbols x N, fftshifted))], one per contiguous run of kept symbols."""
    a, b = max(0, a), min(len(y), b)
    seg = y[a:b]
    if len(seg) < (MIN_SEG + 1) * Lf + N:
        return []
    c = seg[:-N] * np.conj(seg[N:])
    e = (np.abs(seg[:-N]) ** 2 + np.abs(seg[N:]) ** 2) / 2
    cs, es = _movsum(c, ncp), _movsum(e, ncp)
    m = np.abs(cs) / (es + 1e-30)
    n_sym = int((len(m) - 1) // Lf)
    i = np.arange(n_sym)
    G = np.minimum(np.round(np.arange(int(np.ceil(Lf)))[:, None] + i[None, :] * Lf).astype(int), len(m) - 1)
    phase = int(np.argmax(m[G].sum(1)))                   # symbol grid: fold the CP metric at the symbol period
    st = G[phase]
    off = ncp // 2                                         # mid-CP window: only |X|^2 is used
    st = st[st + off + N <= len(seg)]
    ok = m[st] >= CP_MIN
    if ok.sum() < MIN_SEG:
        return []
    eps = np.angle(np.sum(cs[st[ok]])) / (2 * np.pi * N)    # cycles / sample: y[n] ~ exp(-j 2 pi eps n)
    W = st[:, None] + off + np.arange(N)[None, :]
    X = np.fft.fftshift(np.fft.fft(seg[W] * np.exp(2j * np.pi * eps * W), axis=1), axes=1)
    return [(a + int(st[0]), X, ok)]


def tx_streams(y, ofs, tu, ts, bw, runs, snr_db, dur_s):
    """y: the link's slice (centred on the band, `ofs` samples/s with Tu an integer number of samples);
    runs: [(t0_s, t1_s, n)] bursts; returns the decoded-field dict."""
    out = dict(verdict="unknown", tx_antennas=None)
    if snr_db < SNR_MIN:
        out["note"] = f"link SNR {snr_db:.1f} dB < {SNR_MIN:.0f} dB: noise masks the constellation statistics"
        return out
    N = int(round(tu * ofs))
    ncp = int(round((ts - tu) * ofs))
    Lf = ts * ofs
    if ncp < 4:
        out["note"] = "cyclic prefix too short at this rate"
        return out
    k = np.arange(-(N // 2), N - N // 2)
    band = np.abs(k) <= BAND_FRAC * bw / 2 * tu
    if band.sum() < MIN_SC:
        out["note"] = f"{int(band.sum())} subcarriers < {MIN_SC}: too few for the statistics"
        return out
    res = []
    for t0, t1, _ in runs:
        a, b = int(t0 * ofs) - int(Lf), int(t1 * ofs) + int(Lf)
        edge = t0 < 2 * ts or t1 > dur_s - 2 * ts               # cut off by the capture: parity unknown
        for s0, X, ok in segments(y, N, ncp, Lf, a, b):
            Z = X[:, band]
            P = np.abs(Z) ** 2
            ref = np.mean(P[ok], 0) + 1e-30
            occ = np.mean(P > 0.1 * ref, 1)
            keep = ok & (occ >= OCC_MIN)
            # contiguous runs of kept symbols
            idx = np.flatnonzero(keep)
            if len(idx) == 0:
                continue
            brk = np.flatnonzero(np.diff(idx) != 1)
            for g in np.split(idx, brk + 1):
                if len(g) < MIN_SEG:
                    continue
                Zg = Z[g]
                dt, df = _pairs(Zg)
                res.append(dict(n=len(g), kappa=_kurtosis(Zg), dt=dt, df=df,
                                edge=bool(edge or (g[0] == 0 and t0 < 2 * ts))))
    if not res:
        out["note"] = f"no run of >= {MIN_SEG} contiguous full-band data symbols"
        return out
    kap = float(np.median([r["kappa"] for r in res]))
    inner = [r for r in res if not r["edge"]] or res
    dts = np.array([r["dt"] for r in inner])
    dfs = np.array([r["df"] for r in inner])
    dt, df = float(np.median(dts)), float(np.median(dfs))
    out.update(kurtosis=round(kap, 3), time_pair_delta=round(dt, 3), freq_pair_delta=round(df, 3),
               segments=len(res), symbols=int(sum(r["n"] for r in res)), subcarriers=int(band.sum()))

    def paired(d):
        return abs(np.median(d)) >= PAIR_MIN and (np.all(d <= -PAIR_MIN / 2) or np.all(d >= PAIR_MIN / 2))
    if paired(dts) and abs(dt) >= abs(df):
        out.update(verdict="2tx_stbc_like", tx_antennas=2,
                   note=f"Alamouti over OFDM symbol pairs ({'even' if dt < 0 else 'odd'} pairs from the run start): "
                        f"pair power correlation delta {dt:+.2f} on {len(inner)} run(s), kurtosis {kap:.2f}")
    elif np.all(dfs <= -PAIR_MIN) and abs(df) >= PAIR_MIN:
        out.update(verdict="2tx_sfbc_like", tx_antennas=2,
                   note=f"Alamouti over subcarrier pairs: pair power correlation delta {df:+.2f}, kurtosis {kap:.2f}")
    elif kap <= K1_MAX and abs(dt) < PAIR_NULL and abs(df) < PAIR_NULL:
        out.update(verdict="1_stream_like",
                   note=f"kurtosis {kap:.2f} <= {K1_MAX}: one stream per subcarrier (one transmit antenna, or "
                        "several with cyclic delay / switching), no Alamouti pairing")
    else:
        out["note"] = (f"kurtosis {kap:.2f} and pair deltas {dt:+.2f} / {df:+.2f} match no verdict: above a single "
                       "stream without the Alamouti pairing (spatial multiplexing, a varying allocation or "
                       "non-QAM data)")
    return out
