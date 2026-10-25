"""Synthetic swept carriers and look-alikes for the jammer detector tests (complex64, unit-variance noise)."""
from __future__ import annotations

import numpy as np

FS = 61.44e6
N = 2_000_000


def noise(rng, n=N):
    return ((rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)).astype(np.complex64)


def _fm(f, amp):
    return amp * np.exp(1j * 2 * np.pi * np.cumsum(f) / FS)


def sweep(rate, f_lo, f_hi, shape="triangle", snr_db=10.0, jitter_us=0.0, noise_fm_hz=0.0, clip=None, seed=0):
    """carrier swept between f_lo and f_hi (Hz from the LO) `rate` times a second; snr_db = carrier power over
    the noise in 1 MHz. jitter_us: slow random wander of the sweep timing; clip: I/Q hard clip at clip x rms."""
    rng = np.random.default_rng(seed)
    t = np.arange(N) / FS
    if jitter_us:
        w = np.interp(np.arange(N), np.arange(0, N, 6144), np.cumsum(rng.standard_normal(len(range(0, N, 6144)))))
        t = t + w / (np.std(w) + 1e-9) * jitter_us * 1e-6
    ph = (t * rate) % 1.0
    u = {"triangle": lambda: 1 - np.abs(2 * ph - 1), "sawtooth": lambda: ph,
         "sine": lambda: 0.5 - 0.5 * np.cos(2 * np.pi * ph)}[shape]()
    f = f_lo + (f_hi - f_lo) * u
    if noise_fm_hz:
        f = f + noise_fm_hz * rng.standard_normal(N)
    # beyond +-0.45 fs the receiver's anti-alias filter removes the carrier (it would alias otherwise)
    x = _fm(f, 10 ** (snr_db / 20) * np.sqrt(1e6 / FS)) * (np.abs(f) < 0.45 * FS) + noise(rng)
    if clip:
        c = clip * np.std(x)
        x = np.clip(x.real, -c, c) + 1j * np.clip(x.imag, -c, c)
    return x.astype(np.complex64)


def pulsed_radar(prf=1000.0, pw=1e-6, f=5e6, seed=1):
    rng = np.random.default_rng(seed)
    t = np.arange(N) / FS
    on = ((t * prf) % 1.0) < pw * prf
    return (80 * on * np.exp(2j * np.pi * f * t) + noise(rng)).astype(np.complex64)


def lfm_radar(prf=2000.0, pw=50e-6, bw=20e6, seed=2):
    """pulsed chirp radar: a sweep, but only during short pulses"""
    rng = np.random.default_rng(seed)
    tp = (np.arange(N) / FS * prf) % 1.0 / prf
    return (80 * (tp < pw) * _fm(-bw / 2 + bw * tp / pw, 1.0) + noise(rng)).astype(np.complex64)


def ofdm_tdd(period=1e-3, duty=0.6, bw=20e6, snr_db=15, seed=3):
    """wideband noise-like bursts (OFDM-like) on a fixed TDD period"""
    from scipy import signal
    rng = np.random.default_rng(seed)
    s = signal.lfilter(signal.firwin(129, bw / 2, fs=FS), 1, noise(rng)) * 10 ** (snr_db / 20) * np.sqrt(FS / bw)
    on = ((np.arange(N) / FS / period) % 1.0) < duty
    return (s * on + noise(rng)).astype(np.complex64)


def hopper(hop_s=100e-6, freqs=(-20e6, 5e6, -8e6, 17e6, 0.0, -14e6, 11e6), seed=4):
    """fixed periodic hop pattern of narrow 2-FSK carriers"""
    rng = np.random.default_rng(seed)
    k = (np.arange(N) / FS / hop_s).astype(int) % len(freqs)
    f = np.array(freqs)[k] + 250e3 * np.sign(rng.standard_normal(N // 64 + 1)).repeat(64)[:N]
    return (_fm(f, 80) + noise(rng)).astype(np.complex64)


def bt_hopper(seed=8):
    """Bluetooth-like: 625 us slots on pseudo-random 1 MHz channels, half of them used"""
    rng = np.random.default_rng(seed)
    slot = int(625e-6 * FS)
    f, on = np.zeros(N), np.zeros(N, bool)
    for s in range(0, N, slot):
        f[s:s + slot] = rng.integers(-28, 28) * 1e6 + 250e3 * np.sign(rng.standard_normal())
        on[s:s + int(366e-6 * FS)] = rng.random() < 0.5
    return (80 * on * _fm(f, 1.0) + noise(rng)).astype(np.complex64)


def fm_video(line=15625.0, dev=4e6, seed=5):
    """analog-video-like FM: a sync tip each line, a random picture level per line"""
    rng = np.random.default_rng(seed)
    t = np.arange(N) / FS
    ln = (t * line).astype(int)
    pos = (t * line) % 1.0
    lvl = rng.uniform(-0.3, 1.0, ln.max() + 1)[ln] + 0.2 * np.sin(2 * np.pi * pos * 7)
    return (_fm(dev * np.where(pos < 0.073, -1.0, lvl), 80) + noise(rng)).astype(np.complex64)


def wifi_bursts(seed=6, snr_db=20):
    """two 20 MHz OFDM-like channels with independent random bursts (CSMA traffic)"""
    from scipy import signal
    rng = np.random.default_rng(seed)
    out = noise(rng).astype(np.complex128)
    t = np.arange(N) / FS
    for fc in (-15e6, 10e6):
        s = signal.lfilter(signal.firwin(129, 8.3e6, fs=FS), 1, noise(rng)) * 10 ** (snr_db / 20) \
            * np.sqrt(FS / 16.6e6) * np.exp(2j * np.pi * fc * t)
        on = np.zeros(N, bool)
        p = 0
        while p < N:
            p += int(rng.exponential(300e-6) * FS)
            d = int(rng.uniform(40e-6, 2e-3) * FS)
            on[p:p + d] = True
            p += d
        out += s * on
    return out.astype(np.complex64)


def cw(seed=7):
    rng = np.random.default_rng(seed)
    return (30 * np.exp(2j * np.pi * 3e6 * np.arange(N) / FS) + noise(rng)).astype(np.complex64)


def _band(x, f0, f1):
    """x band-limited to [f0, f1] (Hz from the LO) by a brick wall in frequency."""
    X = np.fft.fft(x)
    f = np.fft.fftfreq(len(x), 1 / FS)
    return np.fft.ifft(X * ((f >= f0) & (f <= f1)))


def barrage(f0, f1, snr_db=15.0, seed=10):
    """band-limited complex Gaussian noise on [f0, f1] (Hz from the LO), snr_db over the noise per Hz."""
    rng = np.random.default_rng(seed)
    s = _band(noise(rng), f0, f1) * 10 ** (snr_db / 20)
    return (s + noise(rng)).astype(np.complex64)


def noise_fm(dev_hz, centre=0.0, cutoff=2e6, snr_db=15.0, seed=11):
    """VCO driven by low-pass Gaussian noise: instantaneous frequency centre + dev_hz * N(0, 1) (filtered)."""
    from scipy import signal
    rng = np.random.default_rng(seed)
    m = signal.lfilter(signal.firwin(255, cutoff, fs=FS), 1, rng.standard_normal(N))
    f = centre + dev_hz * m / np.std(m)
    x = _fm(f, 10 ** (snr_db / 20)) * (np.abs(f) < 0.45 * FS)
    return (x + noise(rng)).astype(np.complex64)


def ofdm(bw, tu=66.67e-6, cp=4.69e-6, snr_db=10.0, f_off=0.0, seed=12):
    """continuous CP-OFDM (random QPSK on every subcarrier inside bw)."""
    rng = np.random.default_rng(seed)
    nfft = int(round(tu * FS))
    ncp = int(round(cp * FS))
    k = np.fft.fftfreq(nfft, 1 / FS)
    use = (np.abs(k - f_off) <= bw / 2) & (k != f_off)
    sym = []
    while sum(len(s) for s in sym) < N:
        X = np.zeros(nfft, complex)
        X[use] = (rng.choice([-1, 1], use.sum()) + 1j * rng.choice([-1, 1], use.sum())) / np.sqrt(2)
        s = np.fft.ifft(X) * (nfft / np.sqrt(use.sum()))           # unit power
        sym.append(np.concatenate([s[-ncp:], s]))
    x = np.concatenate(sym)[:N] * 10 ** (snr_db / 20) * np.sqrt(bw / FS)     # snr_db in band
    return (x + noise(rng)).astype(np.complex64)
