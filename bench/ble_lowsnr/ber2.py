import sys
import numpy as np
from ber import chan_known, demod_disc_int, SPS, _kmeans2
from coherent import detect_coherent

def run(snr, trials=60, seed=1, delta=0):
    rng = np.random.default_rng(seed)
    ec = ed = n = 0
    for _ in range(trials):
        y, bits = chan_known(snr, rng)
        b, rel, _ = detect_coherent(y, bits[:40], -9, len(bits))
        m = min(len(b), len(bits)) - 4
        ec += np.sum(b[40:m] != bits[40:m])
        v = demod_disc_int(y); idx = np.arange(len(bits)) * SPS + 4; idx = idx[idx < len(v)]
        vv = v[idx]; hb = vv > np.mean(_kmeans2(vv))
        ed += np.sum(hb[40:m] != bits[40:m].astype(bool)); n += m - 40
    return ec / n, ed / n

for snr in (4, 6, 8, 10):
    print(snr, "coherent %.5f diff %.5f" % run(snr))
