import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from streams import synth, pair_correlation, kurtosis
rng = np.random.default_rng(3)
for kind, M in (('1tx', 64), ('cdd', 64), ('sm', 16), ('sfbc', 16), ('sfbc', 64), ('stbc', 16), ('gauss', 16)):
    Z = synth(kind, 50, 1728, 28, rng, M=M)
    pc = pair_correlation(Z)
    print(f'{kind:5s} {M:3d}QAM: kappa {kurtosis(Z)[0]:.3f}  ' + '  '.join(f'{k} {v:+.3f}' for k, v in pc.items()))
