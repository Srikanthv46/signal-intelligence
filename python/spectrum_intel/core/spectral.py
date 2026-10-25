"""Shared DSP for the receivers: FFT channelizer with spur notches, CW spur finder,
CP metrics, PSD noise floor, spectrogram, SNR measurement (from the research
spectrum_survey / drone_common modules, code unchanged)."""
from __future__ import annotations

import numpy as np
from scipy import ndimage


# ============================================================== channelizer
class Channelizer:
    """FFT-domain channelizer: one FFT of the whole capture, then any
    (centre, width, output rate) slice is a small IFFT. Brick-wall response,
    which is fine for detection and numerology measurements."""

    def __init__(self, x, fs, notch_hz=()):
        self.fs, self.L = fs, len(x)
        self.X = np.fft.fft(x)
        self.f = np.fft.fftfreq(self.L, 1 / fs)
        for f0, half in notch_hz:                         # zero spur bins (CW correlates at every lag)
            self.X[np.abs(self.f - f0) <= half] = 0

    def slice(self, f_off, half_bw, out_fs):
        M = int(round(self.L * out_fs / self.fs))
        sel = np.where(np.abs(self.f - f_off) <= half_bw)[0]
        k0 = int(round(f_off / self.fs * self.L))
        k = np.where(sel > self.L // 2, sel - self.L, sel)       # signed bin: raw index scrambles negative
        Y = np.zeros(M, complex)                                 # bands whenever L is not a multiple of M
        Y[(k - k0) % M] = self.X[sel]
        return np.fft.ifft(Y) * (M / self.L)

# ============================================================== spurs
def find_spurs(x, fs, fc, nfft=8192, min_prom_db=15.0):
    """Persistent narrow tones: median-over-time PSD bins standing well above
    their neighbourhood. Median, so bursty signals do not register."""
    n = len(x) // nfft
    S = np.abs(np.fft.fftshift(np.fft.fft(x[:n * nfft].reshape(n, nfft) * np.hanning(nfft), axis=1), axes=1)) ** 2
    psd = 10 * np.log10(np.median(S, 0) + 1e-20)
    base = ndimage.median_filter(psd, 63)
    f = (np.arange(nfft) - nfft // 2) / nfft * fs
    out = []
    prom = psd - base
    for k in np.where(prom > min_prom_db)[0]:
        if prom[k] < prom[max(k - 1, 0)] or prom[k] < prom[min(k + 1, nfft - 1)]:
            continue                                              # local maxima only
        out.append(dict(offset_hz=float(f[k]), freq_mhz=round((fc + f[k]) / 1e6, 4),
                        prominence_db=round(float(prom[k]), 1),
                        note="LO leakage (DC)" if abs(f[k]) < 2 * fs / nfft else "CW tone"))
    return out

def _cp_metric(y, N, cp):
    c = y[:-N] * np.conj(y[N:])
    e = (np.abs(y[:-N]) ** 2 + np.abs(y[N:]) ** 2) / 2
    k = np.ones(cp)
    return np.abs(np.convolve(c, k, "valid")) / (np.convolve(e, k, "valid") + 1e-20), c

def _folded_scores(m, sym, W, hop):
    """Fold the CP metric modulo the symbol period inside each window: an OFDM
    burst stacks its CP peaks into one sharp bin; noise/Wi-Fi stay flat."""
    out = []
    for a in range(0, len(m) - W, hop):
        k = W // sym
        prof = m[a:a + k * sym].reshape(k, sym).mean(0)
        out.append(prof.max() - np.median(prof))
    return np.array(out)


def add_awgn(x, extra_db, rng=None):
    """Lower the SNR of every signal in x by extra_db by adding white noise.
    The capture's own noise floor N0 is estimated (10th pct of per-bin median PSD)
    and noise of power N0*(10^(extra_db/10)-1) is added, so after the call the
    floor is extra_db higher and every signal's SNR is extra_db lower.
    Caveat: on captures whose floor is not thermal (very low gain, strong-signal skirts)
    the nominal drop is approximate -- evaluations should measure each target's actual
    in-band SNR (measure_snr) rather than trust extra_db."""
    rng = rng or np.random.default_rng(0)
    n0 = noise_floor_psd(x)["floor_power_per_sample"]
    add = n0 * (10 ** (extra_db / 10) - 1)
    w = (rng.standard_normal(len(x)) + 1j * rng.standard_normal(len(x))) * np.sqrt(add / 2)
    return (x + w).astype(np.complex64)

def measure_snr(x, fs, fc, t0_us, t1_us, f_lo_mhz, f_hi_mhz, nfft=1024):
    """In-band SNR of a known (time, freq) box: mean power in the box over the
    per-bin noise floor in the same bins. Use this to report true SNR in evals."""
    P, t, f = spectrogram(x, fs, fc, nfft)
    P /= np.sum(np.hanning(nfft) ** 2)
    fl = noise_floor_psd(x, nfft)["floor_bin"]
    ti = (t >= t0_us) & (t <= t1_us); fi = (f >= f_lo_mhz) & (f <= f_hi_mhz)
    sig = P[np.ix_(ti, fi)].mean(0) - fl[fi]
    return float(10 * np.log10(max(sig.mean(), 1e-12) / fl[fi].mean()))

def noise_floor_psd(x, nfft=1024):
    """Per-bin noise floor (10th percentile over time of the periodogram) and the
    implied per-sample complex noise power across the full band."""
    n = len(x) // nfft
    S = np.abs(np.fft.fftshift(np.fft.fft(x[:n * nfft].reshape(n, nfft) * np.hanning(nfft), axis=1), axes=1)) ** 2
    S /= np.sum(np.hanning(nfft) ** 2)
    floor_bin = np.percentile(S, 25, axis=0) / 0.2877   # 25th pct of exp(1) = -ln(0.75) x mean
    floor_bin = ndimage.median_filter(floor_bin, 31)
    # receiver thermal floor: low percentile across bins -- some bins are never quiet
    # (always-on emitters, filter skirts), so the median over-estimates it
    return dict(floor_bin=floor_bin, floor_power_per_sample=float(np.percentile(floor_bin, 20)), nfft=nfft)

def spectrogram(x, fs, fc, nfft=1024, hop=None):
    """Power spectrogram (time x freq, linear), time axis us, freq axis MHz."""
    hop = hop or nfft
    n = (len(x) - nfft) // hop + 1
    idx = np.arange(nfft)[None, :] + hop * np.arange(n)[:, None]
    P = np.abs(np.fft.fftshift(np.fft.fft(x[idx] * np.hanning(nfft), axis=1), axes=1)) ** 2
    t_us = (np.arange(n) * hop + nfft / 2) / fs * 1e6
    f_mhz = ((np.arange(nfft) - nfft // 2) / nfft * fs + fc) / 1e6
    return P, t_us, f_mhz

def robust_threshold(v, k_sigma):
    med = np.median(v)
    return med + k_sigma * 1.4826 * np.median(np.abs(v - med)), med
