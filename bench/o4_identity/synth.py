"""Bench: synthetic O4-like transmitter (structure of FINDINGS.md) with a known crystal offset, for testing o4.py.

One crystal: the carrier is offset by f_nom x eps and every symbol starts at its nominal time / (1 + eps) (exact
fractional delay applied in the frequency domain). The pilot sequence is defined on the 20 MHz grid (bins
-300..300, odd bins, odd slot positions, depending on the slot within the 2.5 ms frame); 10 MHz mode carries bins
-300..0 of it, re-centred (slice offset 150 bins).
"""
from __future__ import annotations

import numpy as np

from o4 import FS, N, SLOT, FRAME_SLOTS

CPS = [176] + [144] * 13
OFFS = np.cumsum([0] + [N + c for c in CPS[:-1]])          # symbol starts within the slot
REF_POS = (1, 3, 5, 7, 9, 11, 13)


def pilot_table(seed=7):
    """Pilot values on the 20 MHz grid: [slot-in-frame][position] -> dict(bin -> QPSK)."""
    r = np.random.default_rng(seed)
    bins = np.arange(-299, 300, 2)                             # odd bins
    return {(sl, p): dict(zip(bins, np.exp(1j * np.pi / 2 * r.integers(0, 4, len(bins)) + 1j * np.pi / 4)))
            for sl in range(FRAME_SLOTS) for p in REF_POS}


def o4_tx(n_samples, bw_mhz=20, eps_ppm=0.0, f_off_hz=0.0, snr_db=20.0, t0=12345.6, on_slots=(0, 1, 2, 3),
          seed=1, pilots=None, amp_ramp_db=0.0, noise=True):
    """Complex baseband at 61.44 MSps: one O4-like link centred at f_off_hz (+ its crystal offset applied by the
    caller through f_off_hz), symbol clock (1 + eps)."""
    r = np.random.default_rng(seed)
    pilots = pilots or pilot_table()
    eps = eps_ppm * 1e-6
    half = 300 if bw_mhz == 20 else 150
    used = np.array([k for k in range(-half, half + 1) if k != 0])
    y = np.zeros(n_samples + 2 * SLOT, complex)
    nslot = int(n_samples / SLOT) + 2
    for m in range(-1, nslot):
        sl_in_frame = m % FRAME_SLOTS
        if sl_in_frame not in on_slots:
            continue
        for p in range(14):
            t_nom = t0 + m * SLOT + OFFS[p]
            t = t_nom / (1 + eps)                                # transmitter clock (1 + eps) faster
            i0 = int(np.floor(t)); frac = t - i0
            if i0 < 0 or i0 + N + CPS[p] >= len(y):
                continue
            X = np.zeros(N, complex)
            X[used % N] = np.exp(1j * np.pi / 2 * r.integers(0, 4, len(used)) + 1j * np.pi / 4)
            if p in REF_POS:
                tab = pilots[(sl_in_frame, p)]
                for k in used:
                    kk = k if bw_mhz == 20 else k - 150          # 10 MHz: bins -300..0 of the 20 MHz sequence
                    if kk in tab:
                        X[k % N] = tab[kk]
            k = np.fft.fftfreq(N, 1 / N)
            s = np.fft.ifft(X * np.exp(-2j * np.pi * k * frac / N)) * (N / np.sqrt(len(used)))   # unit power
            y[i0:i0 + CPS[p] + N] += np.concatenate([s[-CPS[p]:], s])
    y = y[:n_samples]
    n = np.arange(n_samples)
    y = y * np.exp(2j * np.pi * f_off_hz * n / FS)
    if amp_ramp_db:
        y = y * 10 ** (np.linspace(0, amp_ramp_db, n_samples) / 20)
    if noise:
        bw_frac = len(used) * FS / N / FS
        sig = 1 / np.sqrt(bw_frac)                                # in-band SNR -> full-band noise level
        y = y + sig * 10 ** (-snr_db / 20) / np.sqrt(2) * (r.standard_normal(n_samples) + 1j * r.standard_normal(n_samples))
    return y.astype(np.complex64)
