"""How many streams share each subcarrier of a CP-OFDM link: the channel-free kurtosis test (bench tool).

With one transmit stream per subcarrier, X_i(k) = H(k) s_i(k) and the static one-tap channel cancels from
    kappa(k) = mean_i |X_i(k)|^4 / (mean_i |X_i(k)|^2)^2
leaving the constellation's own kurtosis (QPSK 1.00, 16-QAM 1.32, 64-QAM 1.38, 256-QAM 1.40). Two antennas
sending different symbols on the same subcarrier (Alamouti STBC/SFBC, spatial multiplexing) give
X = h1 s1 + h2 s2, whose kurtosis is a^4 k1 + b^4 k2 + 4 a^2 b^2 (a^2 = |h1|^2 / (|h1|^2 + |h2|^2)): about 1.5-1.7
for comparable branch powers. Cyclic delay diversity and antenna switching keep one stream per subcarrier and
stay at the QAM value (CDD shows as two taps in the channel instead: see `channel_taps`).

Each symbol is first normalised to unit mean power (TX power steps between symbols would otherwise inflate
kappa), and the finite-sample bias of the estimator is calibrated by Monte Carlo at the same symbol count and
SNR (`reference`).
"""
import numpy as np

QAM = {'QPSK': 4, '16QAM': 16, '64QAM': 64, '256QAM': 256}


def qam(M, n, rng):
    m = int(np.sqrt(M))
    lv = np.arange(-(m - 1), m, 2, dtype=float)
    s = rng.choice(lv, n) + 1j * rng.choice(lv, n)
    return s / np.sqrt(np.mean(np.abs(lv) ** 2) * 2)


def kurtosis(Z, per_symbol_norm=True):
    """Z: symbols x subcarriers (occupied data subcarriers only). Median over subcarriers of the per-subcarrier
    kurtosis; per-symbol power normalisation first."""
    Z = np.asarray(Z)
    if per_symbol_norm:
        Z = Z / np.sqrt(np.mean(np.abs(Z) ** 2, 1, keepdims=True))
    p = np.abs(Z) ** 2
    k = np.mean(p ** 2, 0) / np.mean(p, 0) ** 2
    return float(np.median(k)), k


def synth(kind, n_sym, n_sc, snr_db, rng, M=16, p2=0.5, taps=6):
    """Frequency-domain symbols (n_sym x n_sc) through a static random multipath channel per TX antenna.
    kind: '1tx', 'cdd' (same symbols, 2nd antenna delayed), 'sfbc' (Alamouti over subcarrier pairs),
    'stbc' (Alamouti over symbol pairs), 'sm' (two independent streams), 'gauss' (DFT-spread-like)."""
    def chan():
        h = (rng.standard_normal(taps) + 1j * rng.standard_normal(taps)) * np.exp(-np.arange(taps) / 2)
        return np.fft.fft(h, 2048)[:n_sc] / np.sqrt(np.sum(np.abs(h) ** 2))
    H1, H2 = chan() * np.sqrt(1 - p2), chan() * np.sqrt(p2)
    S1 = qam(M, n_sym * n_sc, rng).reshape(n_sym, n_sc)
    S2 = qam(M, n_sym * n_sc, rng).reshape(n_sym, n_sc)
    if kind == '1tx':
        X = H1 * S1
    elif kind == 'cdd':
        d = 8
        X = (H1 + H2 * np.exp(-2j * np.pi * np.arange(n_sc) * d / 2048)) * S1
    elif kind == 'sm':
        X = H1 * S1 + H2 * S2
    elif kind == 'sfbc':
        X = np.empty_like(S1)
        a, b = S1[:, 0::2], S1[:, 1::2]
        X[:, 0::2] = H1[0::2] * a + H2[0::2] * b
        X[:, 1::2] = -H1[1::2] * np.conj(b) + H2[1::2] * np.conj(a)
    elif kind == 'stbc':
        X = np.empty_like(S1)
        a, b = S1[0::2], S1[1::2]
        X[0::2] = H1 * a + H2 * b
        X[1::2] = -H1 * np.conj(b) + H2 * np.conj(a)
    elif kind == 'gauss':
        X = H1 * (rng.standard_normal(S1.shape) + 1j * rng.standard_normal(S1.shape)) / np.sqrt(2)
    else:
        raise ValueError(kind)
    pw = np.mean(np.abs(X) ** 2)
    nz = (rng.standard_normal(X.shape) + 1j * rng.standard_normal(X.shape)) * np.sqrt(pw / 2 * 10 ** (-snr_db / 10))
    return X + nz


def reference(n_sym, n_sc, snr_db, trials=4, seed=1):
    """Median kurtosis of each hypothesis at this symbol count / SNR (mean over trials)."""
    rng = np.random.default_rng(seed)
    out = {}
    hyp = [('1tx ' + q, '1tx', M, 0.5) for q, M in QAM.items()] + \
          [('cdd 64QAM', 'cdd', 64, 0.5)] + \
          [(f'{k} {q} p2={p:.2f}', k, M, p) for k in ('sfbc', 'sm') for q, M in (('QPSK', 4), ('16QAM', 16), ('64QAM', 64))
           for p in (0.5, 0.2)] + [('gauss', 'gauss', 16, 0.5)]
    for name, kind, M, p in hyp:
        out[name] = float(np.mean([kurtosis(synth(kind, n_sym, n_sc, snr_db, rng, M=M, p2=p))[0] for _ in range(trials)]))
    return out


def channel_taps(H_comb, comb=2, N=2048):
    """Power-delay profile from a channel sampled on every `comb`-th subcarrier (contiguous array).
    Returns |h(n)|^2 for n in samples at N*SCS, normalised to the strongest tap; aliasing period N/comb."""
    M = N // comb
    h = np.fft.ifft(np.pad(H_comb, (0, M - len(H_comb))))
    p = np.abs(h) ** 2
    return p / p.max()


def pair_correlation(Z):
    """Alamouti fingerprint. Within an Alamouti pair |y1|^2 + |y2|^2 = (|h1|^2 + |h2|^2)(|a|^2 + |b|^2), so the
    pair's powers are negatively correlated (rho ~ 2|h1 h2|^2 (kappa_s - 2) / var ~ -0.5 for equal branches);
    unpaired neighbours and single-antenna / spatially multiplexed symbols give rho ~ 0.
    Z: symbols x contiguous occupied subcarriers. Returns the mean power correlation of neighbour pairs:
    {'freq_even': (k, k+1) for even k, 'freq_odd': odd k, 'time_even': (i, i+1) for even i, 'time_odd': odd i}."""
    Z = Z / np.sqrt(np.mean(np.abs(Z) ** 2, 1, keepdims=True))
    p = np.abs(Z) ** 2

    def rho(a, b, axis):
        a = a - a.mean(axis, keepdims=True)
        b = b - b.mean(axis, keepdims=True)
        return float(np.median(np.sum(a * b, axis) / np.sqrt(np.sum(a * a, axis) * np.sum(b * b, axis))))
    n = p.shape[1] // 2 * 2
    m = p.shape[0] // 2 * 2
    q = p / p.mean(0, keepdims=True)          # time axis: remove the channel's |H(k)|^2 first
    return {'freq_even': rho(p[:, 0:n:2], p[:, 1:n:2], 0),
            'freq_odd': rho(p[:, 1:n - 1:2], p[:, 2:n:2], 0),
            'time_even': rho(q[0:m:2], q[1:m:2], 1),
            'time_odd': rho(q[1:m - 1:2], q[2:m:2], 1)}
