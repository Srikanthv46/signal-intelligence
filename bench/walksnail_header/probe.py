"""Bench (PLAN D2): first OFDM symbols of each Walksnail-like frame: BPSK test and frame-to-frame sign flips."""
import sys
import numpy as np
sys.path.insert(0, "python")
from spectrum_intel.core import common as C
from spectrum_intel.receivers.digital_vtx import detector as D

x, fs, fc, _ = C.load_capture("testdata/live/LIVE_OFDM_VTX_5805MHz/LIVE_OFDM_VTX_5805MHz.sigmf-meta")
x = x[:1 << 21]
X = np.fft.fft(x.astype(np.complex128))
N, CP = 2048, 256
ofs = N / 91.4286e-6                                    # 22.4 MSps
y, ofs = D._slice(X, fs, 5805.07e6 - fc, 10.0e6, ofs)
L = N + CP
# CP metric per sample
m = np.abs(y[:-N] * np.conj(y[N:]))
cpm = np.convolve(y[:-N] * np.conj(y[N:]), np.ones(CP), "valid")
# frames 1 and 2 agree with the 7.1629 ms period: anchor every frame on that grid (sample of frame 1's first symbol)
P = 387516 - 227067
k = np.arange(-N // 2, N // 2)
used = (np.abs(k) <= 864) & (k != 0)


def sym(first, j):
    s0 = first + CP + j * L
    return np.fft.fftshift(np.fft.fft(y[s0:s0 + N]))[used]


def realness(A, B):
    r = B * np.conj(A)
    g = np.sum(r[1:] * np.conj(r[:-1]))                 # linear phase ramp (timing offset between frames)
    r = r * np.exp(-1j * np.angle(g) * np.arange(len(r)))
    r *= np.exp(-1j * np.angle(np.sum(r ** 2)) / 2)     # common phase, sign-blind (BPSK flips allowed)
    ang = np.angle(r)
    return np.mean(np.minimum(np.abs(ang), np.pi - np.abs(ang)) < 0.4), np.mean(np.abs(ang) > np.pi / 2)


for j in range(-6, 9):
    vals = [realness(sym(227067 + kf * P, j), sym(227067 + (kf + 1) * P, j)) for kf in (-1, 0, 1)]
    print(f"symbol {j:+d}: " + "  ".join(f"real {r:.2f} flip {f:.2f}" for r, f in vals))
