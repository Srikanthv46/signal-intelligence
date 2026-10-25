"""receivers.digital_vtx.streams on synthetic CP-OFDM bursts (time domain, 10.94 kHz SCS, CP 1/8, 2x2 multipath)."""
import numpy as np
import pytest

from spectrum_intel.receivers.digital_vtx import streams as ST

N, CP, OFS = 1024, 128, 11.2e6          # Tu 91.43 us (802.16e numerology)
USED = np.r_[-400:0, 1:401]             # 800 subcarriers, DC empty
NSYM, NBURST, GAP = 50, 4, 3000


def qam(M, n, rng):
    m = int(np.sqrt(M))
    lv = np.arange(-(m - 1), m, 2, dtype=float)
    return (rng.choice(lv, n) + 1j * rng.choice(lv, n)) / np.sqrt(np.mean(lv ** 2) * 2)


def chan(rng, taps=6):
    h = (rng.standard_normal(taps) + 1j * rng.standard_normal(taps)) * np.exp(-np.arange(taps) / 2)
    return h / np.sqrt(np.sum(np.abs(h) ** 2))


def make(kind, rng, M=16, snr_db=30.0, p2=0.5):
    h1, h2 = chan(rng) * np.sqrt(1 - p2), chan(rng) * np.sqrt(p2)
    H1, H2 = np.fft.fft(h1, N)[USED % N], np.fft.fft(h2, N)[USED % N]
    out, runs, t = [np.zeros(GAP, complex)], [], GAP
    for _ in range(NBURST):
        S1 = qam(M, NSYM * len(USED), rng).reshape(NSYM, -1)
        S2 = qam(M, NSYM * len(USED), rng).reshape(NSYM, -1)
        if kind == "1tx":
            F = H1 * S1
        elif kind == "cdd":
            F = (H1 + H2 * np.exp(-2j * np.pi * USED * 8 / N)) * S1
        elif kind == "sm":
            F = H1 * S1 + H2 * S2
        elif kind == "stbc":
            F = np.empty_like(S1)
            F[0::2] = H1 * S1[0::2] + H2 * S1[1::2]
            F[1::2] = -H1 * np.conj(S1[1::2]) + H2 * np.conj(S1[0::2])
        elif kind == "sfbc":
            F = np.empty_like(S1)
            a, b = S1[:, 0::2], S1[:, 1::2]
            F[:, 0::2] = H1[0::2] * a + H2[0::2] * b
            F[:, 1::2] = -H1[1::2] * np.conj(b) + H2[1::2] * np.conj(a)
        G = np.zeros((NSYM, N), complex)
        G[:, USED % N] = F
        x = np.fft.ifft(G, axis=1)
        x = np.hstack([x[:, -CP:], x]).ravel()
        runs.append((t / OFS, (t + len(x)) / OFS, NSYM))
        out += [x, np.zeros(GAP, complex)]
        t += len(x) + GAP
    y = np.concatenate(out)
    pw = np.mean(np.abs(y[np.abs(y) > 0]) ** 2)
    y = y + (rng.standard_normal(len(y)) + 1j * rng.standard_normal(len(y))) * np.sqrt(pw / 2 * 10 ** (-snr_db / 10))
    y = y * np.exp(2j * np.pi * 3.1e3 / OFS * np.arange(len(y)))          # fractional CFO
    return y, runs


def verdict(kind, seed, **kw):
    rng = np.random.default_rng(seed)
    y, runs = make(kind, rng, **kw)
    bw = (USED.max() - USED.min()) * OFS / N
    return ST.tx_streams(y, OFS, N / OFS, (N + CP) / OFS, bw, runs, 30.0, len(y) / OFS)


@pytest.mark.parametrize("kind,M,want", [("1tx", 64, "1_stream_like"), ("1tx", 4, "1_stream_like"),
                                         ("cdd", 64, "1_stream_like"), ("stbc", 16, "2tx_stbc_like"),
                                         ("stbc", 64, "2tx_stbc_like"), ("sfbc", 16, "2tx_sfbc_like"),
                                         ("sm", 16, "unknown")])
def test_verdict(kind, M, want):
    r = verdict(kind, 7, M=M)
    assert r["verdict"] == want, r
    assert (r["tx_antennas"] == 2) == want.startswith("2tx")


def test_weak_second_branch_still_stbc():
    assert verdict("stbc", 11, M=64, p2=0.08)["verdict"] == "2tx_stbc_like"


def test_low_snr_is_unknown():
    rng = np.random.default_rng(1)
    y, runs = make("stbc", rng)
    r = ST.tx_streams(y, OFS, N / OFS, (N + CP) / OFS, 8.75e6, runs, 12.0, len(y) / OFS)
    assert r["verdict"] == "unknown" and "SNR" in r["note"]
