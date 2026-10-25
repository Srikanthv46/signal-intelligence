"""Real-capture check of the stream tests (streams.kurtosis / pair_correlation) on 802.11 PPDUs whose transmit
mode the SIG fields state: non-HT and HT single stream (1 stream per subcarrier), HT STBC (Alamouti over
symbol pairs, 2 TX), HT MCS 8-15 (2 spatial streams, 2 TX).

    .venv/bin/python -I bench/tx_antennas/wifi_reference.py <capture.sigmf-meta | .cf32> [--scaled] [--min-sym 16]

Uses the library's detectors (receivers.wifi.ofdm_rx, ofdm_scaled) only to find the PPDUs and their state; the
data symbols are the raw FFT of every data OFDM symbol on the data subcarriers (pilots excluded), no
equalisation, so the tests see exactly what they see on an unknown link.
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', 'python'))
from streams import kurtosis, pair_correlation          # noqa: E402
from spectrum_intel.receivers.wifi import ofdm_rx as W, payload as WP, ldpc as LD   # noqa: E402


def load(path):
    if path.endswith('.cf32'):
        m = json.load(open(path[:-5] + '.json'))
        return np.fromfile(path, np.complex64), m['fs_hz'], m['fc_hz']
    from spectrum_intel.io import sigmf
    cap = sigmf.load(path)
    return cap.iq, cap.fs, cap.fc


def data_symbols(p):
    """(mode, n_sym, Z) for a decoded-preamble PPDU: Z = n_sym x data subcarriers, raw FFT."""
    rx = p.extra.get('_rx')
    if not rx or rx[0] is None:
        return None
    st = rx[0]
    r, start, cfo, back = st['r'], st['start'], st['cfo'], st['back']
    e = p.extra
    if p.fmt == 'NON_HT':
        if p.lsig_rate_mbps not in WP.NONHT_MCS:
            return None
        nbpsc, cr = WP.NONHT_MCS[p.lsig_rate_mbps]
        ndbps = 48 * nbpsc * cr[0] // cr[1]
        n_sym = int(np.ceil((16 + 8 * p.lsig_length + 6) / ndbps))
        data0, sym, gi, bins = 400, 80, 16, W._bins(W.DATA_SC)
        try:                                   # ground truth only when the payload is what L-SIG says
            ok = WP.decode_nonht_payload(st, p.lsig_rate_mbps, p.lsig_length, rx[1], rx[2])[1]
        except Exception:                      # noqa: BLE001
            ok = False
        mode = f'NON_HT {p.lsig_rate_mbps:g} Mb/s (1 stream)' + (' FCS ok' if ok else ' FCS FAIL')
    elif p.fmt == 'HT' and p.mcs is not None and 'ht_length' in e:
        nss = p.mcs // 8 + 1
        stbc = int(e.get('stbc', 0))
        nbpsc, cr = WP.HT_MCS[p.mcs % 8]
        ncbps = 52 * nbpsc * nss
        ndbps = ncbps * cr[0] // cr[1]
        m_stbc = 2 if stbc else 1
        L = e['ht_length']
        if e.get('ldpc'):
            n_sym = LD.params(L, ncbps, cr, m_stbc)['nsym']
        else:
            n_sym = m_stbc * int(np.ceil((16 + 8 * L + 6) / (m_stbc * ndbps)))
        nsts = nss + stbc
        n_ltf = {1: 1, 2: 2, 3: 4, 4: 4}[nsts]
        data0 = 640 + 80 * (n_ltf + e.get('n_ess', 0))
        sym, gi = (72, 8) if e.get('short_gi') else (80, 16)
        bins = W._bins(W.HT_DATA_SC)
        mode = (f'HT MCS{p.mcs} ' + ('STBC (2 TX, Alamouti)' if stbc else
                                     f'{nss} spatial streams (2 TX)' if nss > 1 else '1 stream'))
    else:
        return None
    need = start + data0 + sym * n_sym
    if need > len(r):
        return None
    seg = r[start:need] * np.exp(-2j * np.pi * cfo / W.FS_WIFI * (start + np.arange(need - start)))
    Z = np.array([np.fft.fft(seg[data0 + sym * k + gi - back:data0 + sym * k + gi - back + 64])[bins]
                  for k in range(n_sym)])
    return mode, n_sym, Z


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('path')
    ap.add_argument('--scaled', action='store_true', help='802.11 at a vendor clock (receivers.wifi.ofdm_scaled)')
    ap.add_argument('--min-sym', type=int, default=16)
    a = ap.parse_args()
    x, fs, fc = load(a.path)
    if a.scaled:
        from spectrum_intel.receivers.wifi import ofdm_scaled as S
        pk = [p for ln in S.scan(x, fs, fc) for p in ln.packets]
    else:
        pk = W.WifiDetector(fs, fc).run(x)
    rows = []
    for p in pk:
        d = data_symbols(p)
        if d is None or d[1] < a.min_sym:
            continue
        mode, n, Z = d
        kap = kurtosis(Z)[0]
        pc = pair_correlation(Z)
        rows.append((mode, n, kap, pc))
        print(f'{p.start_us / 1e3:8.3f} ms  {mode:34s} {n:4d} sym  kappa {kap:.3f}  ' +
              '  '.join(f'{k} {v:+.3f}' for k, v in pc.items()))
    by = {}
    for mode, n, kap, pc in rows:
        by.setdefault(mode, []).append((kap, pc['time_even'], pc['time_odd'], pc['freq_even']))
    print('\nper mode (median): kappa, time_even, time_odd, freq_even')
    for mode, v in by.items():
        v = np.array(v)
        print(f'  {mode:34s} n={len(v):3d}  ' + '  '.join(f'{x:+.3f}' for x in np.median(v, 0)))


if __name__ == '__main__':
    main()
