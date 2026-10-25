"""How strong is the STBC pair signature when the second antenna's branch is weaker (p2 = its power share)?"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from streams import synth, pair_correlation, kurtosis
rng = np.random.default_rng(5)
for M in (16, 64):
    for p2 in (0.5, 0.35, 0.25, 0.15, 0.08, 0.03):
        v = []
        for _ in range(3):
            Z = synth('stbc', 50, 1728, 28, rng, M=M, p2=p2)
            pc = pair_correlation(Z)
            v.append((kurtosis(Z)[0], pc['time_even'] - pc['time_odd']))
        k, d = np.mean(v, 0)
        print(f'STBC {M:2d}QAM p2={p2:.2f} ({10*np.log10(p2/(1-p2)):+5.1f} dB): kappa {k:.3f}  delta {d:+.3f}')
