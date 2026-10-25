"""Transmit-antenna test for an unknown CP-OFDM link from one receive channel (bench tool).

    .venv/bin/python -I bench/tx_antennas/link_streams.py <capture> --offset-hz 5.07e6 --nfft 2048 --cp 256 \
        --fs2 22.4e6 [--band-sc 864]

Finds every run of contiguous symbols (CP correlation, fixed grid per run), then per run of >= 20 symbols:
  kappa      per-subcarrier kurtosis (streams.kurtosis): <= 1.40 one stream per subcarrier (1 TX, CDD, switching),
             1.45-1.65 two streams (2 TX: STBC or spatial multiplexing), ~1.9 Gaussian-like data;
  delta      time-pair power correlation, paired parity minus the other (streams.pair_correlation): ~0 for one
             stream and for spatial multiplexing, clearly negative (-0.15..-0.45) for Alamouti STBC over symbol
             pairs; the paired parity must follow the frame (pairs start at the first data symbol);
  freq_delta the same over subcarrier pairs (SFBC).
Real-capture reference: wifi_reference.py (802.11 PPDUs with SIG-stated modes).
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from ofdm_sync import load, baseband, symbols          # noqa: E402
from streams import kurtosis, pair_correlation          # noqa: E402


def runs_of(st, L, tol=4, min_len=20):
    out, cur = [], [st[0]]
    for a, b in zip(st[:-1], st[1:]):
        if abs(b - a - L) <= tol:
            cur.append(b)
        else:
            out.append(cur)
            cur = [b]
    out.append(cur)
    return [np.array(r) for r in out if len(r) >= min_len]


def analyse(path, offset_hz, N, cp, fs2, band_sc):
    x, fs, fc = load(path)
    y = baseband(x, fs, offset_hz, fs2)
    st, _, eps, _ = symbols(y, N, cp)
    y = y * np.exp(-2j * np.pi * eps * np.arange(len(y)))
    L = N + cp
    k = np.arange(-N // 2, N // 2)
    res = []
    for r in runs_of(st, L):
        i = np.arange(len(r))
        A = np.polyfit(i, r, 1)
        s = np.round(A[1] + A[0] * i).astype(int)
        D = np.array([np.fft.fftshift(np.fft.fft(y[q + cp:q + L])) for q in s])
        P = np.mean(np.abs(D) ** 2, 0)
        occ = (P > 0.1 * np.median(P[np.abs(k) <= band_sc])) & (np.abs(k) <= band_sc + 8)
        Z = D[:, occ]
        pc = pair_correlation(Z)
        res.append(dict(t_ms=s[0] / fs2 * 1e3, n_sym=len(s), n_sc=int(occ.sum()), kappa=kurtosis(Z)[0],
                        delta_even=pc['time_even'] - pc['time_odd'],
                        freq_delta=pc['freq_even'] - pc['freq_odd'], **pc))
    return res, eps * fs2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('path')
    ap.add_argument('--offset-hz', type=float, required=True)
    ap.add_argument('--nfft', type=int, default=2048)
    ap.add_argument('--cp', type=int, default=256)
    ap.add_argument('--fs2', type=float, default=22.4e6)
    ap.add_argument('--band-sc', type=int, default=864)
    a = ap.parse_args()
    res, cfo = analyse(a.path, a.offset_hz, a.nfft, a.cp, a.fs2, a.band_sc)
    print(f'fractional CFO {cfo:.0f} Hz; {len(res)} runs of >= 20 symbols')
    for r in res:
        print(f"{r['t_ms']:8.3f} ms {r['n_sym']:3d} sym {r['n_sc']:4d} sc  kappa {r['kappa']:.3f}  "
              f"time pairs even {r['time_even']:+.3f} odd {r['time_odd']:+.3f} (delta {r['delta_even']:+.3f})  "
              f"freq delta {r['freq_delta']:+.3f}")


if __name__ == '__main__':
    main()
