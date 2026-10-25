"""BER of BLE 1M GFSK demodulators vs in-band SNR (1 MHz), synthetic only (PLAN L1)."""
import sys
import numpy as np
from scipy import ndimage, signal
from spectrum_intel.receivers.rc_link import synth as S
from spectrum_intel.receivers.rc_link.detector import _kmeans2

OFS, SPS = 8e6, 8


def chan(snr_db, rng, nbits=400):
    bits = rng.integers(0, 2, nbits).astype(np.uint8)
    iq, _ = S.gfsk_burst(1e6, 250e3, bits=bits, bt=0.5, cfo_hz=rng.uniform(-40e3, 40e3), fs=OFS, rng=rng)
    n0 = 1 / (1e6 * 10 ** (snr_db / 10))
    w = (rng.standard_normal(len(iq)) + 1j * rng.standard_normal(len(iq))) * np.sqrt(n0 * OFS / 2)
    y = iq + w
    Y = np.fft.fft(y); f = np.fft.fftfreq(len(y), 1 / OFS)
    Y[np.abs(f) > 0.75e6] = 0                       # the pipeline's brick-wall slice
    return np.fft.ifft(Y), bits


def demod_current(y):
    ifr = np.angle(y[1:] * np.conj(y[:-1]))
    v = ndimage.uniform_filter1d(ifr, SPS // 2)
    return v                                          # sample at k*SPS + SPS/2 - 1 (diff index)


def lpf(y, fc):
    h = signal.firwin(33, fc, fs=OFS)
    return np.convolve(y, h, "same")


def demod_diff(y, fc=0.6e6):
    z = lpf(y, fc)
    d = np.zeros(len(z)); h = SPS // 2
    p = z[SPS:] * np.conj(z[:-SPS])                 # phase change over one bit
    d[h:h + len(p)] = np.imag(p) / (np.abs(p) + 1e-12)
    return d


def demod_disc_int(y, fc=0.6e6):
    z = lpf(y, fc)
    d = np.zeros(len(z)); h = SPS // 2
    p = z[SPS:] * np.conj(z[:-SPS])
    d[h:h + len(p)] = np.angle(p)
    return d


def ber(fn, snr, trials=60, seed=1, off=0):
    rng = np.random.default_rng(seed)
    e = n = 0
    for _ in range(trials):
        y, bits = chan(snr, rng)
        v = fn(y)
        best = None
        for o in range(SPS):                          # best timing phase (genie) and CFO-removed threshold
            idx = (np.arange(len(bits)) * SPS + o)
            idx = idx[idx < len(v)]
            vv = v[idx]; thr = np.mean(_kmeans2(vv))
            err = np.sum(((vv > thr) != bits[:len(vv)].astype(bool))[4:-4])     # ramp / filter edges excluded
            best = err if best is None else min(best, err)
        e += best; n += len(bits) - 8
    return e / n


if __name__ == "__main__":
    for snr in (6, 8, 10, 12):
        print(snr, "current %.4f" % ber(demod_current, snr), "diff_im %.4f" % ber(demod_diff, snr),
              "diff_ang %.4f" % ber(demod_disc_int, snr),
              "diff_im0.5 %.4f" % ber(lambda y: demod_diff(y, 0.5e6), snr),
              "diff_im0.75 %.4f" % ber(lambda y: demod_diff(y, 0.75e6), snr))


def chan_known(snr_db, rng, nbits=400):
    """Like chan, but bits[0:40] are the BLE preamble + AA (known), and the burst starts at sample 0."""
    aa = [(0x8E89BED6 >> i) & 1 for i in range(32)]
    pre = [aa[0] ^ (i % 2 == 0) for i in range(8)]
    bits = np.r_[pre, aa, rng.integers(0, 2, nbits - 40)].astype(np.uint8)
    iq, _ = S.gfsk_burst(1e6, 250e3, bits=bits, bt=0.5, cfo_hz=rng.uniform(-40e3, 40e3), fs=OFS, rng=rng)
    n0 = 1 / (1e6 * 10 ** (snr_db / 10))
    w = (rng.standard_normal(len(iq)) + 1j * rng.standard_normal(len(iq))) * np.sqrt(n0 * OFS / 2)
    Y = np.fft.fft(iq + w); f = np.fft.fftfreq(len(iq), 1 / OFS)
    Y[np.abs(f) > 0.75e6] = 0
    return np.fft.ifft(Y), bits
