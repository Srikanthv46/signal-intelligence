"""Synthetic validation: does the kurtosis test separate 1-stream from 2-stream transmissions at this link's
symbol count (50 per run) and SNR (25-31 dB)?"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from streams import reference
for snr in (31, 25, 20, 15):
    r = reference(50, 1728, snr)
    print(f'SNR {snr} dB, 50 symbols x 1728 subcarriers:')
    for k, v in r.items():
        print(f'   {k:22s} {v:.3f}')
