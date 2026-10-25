"""core.bursts: band and burst extent of a known link from in-band power (synthetic)."""
import numpy as np

from spectrum_intel.core import bursts as BU

FS = 30.72e6
N = 1 << 20                                   # 34.1 ms
LO, HI = -1.5e6, 7.5e6                        # the link: 9 MHz at +3 MHz
BURSTS = [(1.0e-3 + k * 2.5e-3, 1.0e-3 + k * 2.5e-3 + 1.5e-3) for k in range(13)]


def _capture(link_db=0.0, wifi=True, seed=1):
    rng = np.random.default_rng(seed)
    noise = (rng.standard_normal(N) + 1j * rng.standard_normal(N)) / np.sqrt(2)
    f = np.fft.fftfreq(N, 1 / FS)
    # the link: band-limited noise with a multipath ripple (+-3 dB), on in BURSTS
    S = (rng.standard_normal(N) + 1j * rng.standard_normal(N)) * ((f >= LO) & (f <= HI))
    S *= 10 ** (1.5 * np.cos(2 * np.pi * f / 1.7e6) / 20)
    s = np.fft.ifft(S)
    s *= np.sqrt(10 ** (link_db / 10) * (HI - LO) / FS / np.mean(np.abs(s) ** 2))   # in-band SNR link_db
    on = np.zeros(N, bool)
    for a, b in BURSTS:
        on[int(a * FS):int(b * FS)] = True
    x = noise + s * on
    if wifi:
        # a strong 16 MHz burst over the lower half of the band, right after burst 4 (it must not extend it)
        W = (rng.standard_normal(N) + 1j * rng.standard_normal(N)) * ((f >= -10e6) & (f <= 6e6))
        w = np.fft.ifft(W)
        w *= np.sqrt(10 * 16e6 / FS / np.mean(np.abs(w) ** 2))                   # 10 dB over the noise
        a, b = int((BURSTS[4][1] + 20e-6) * FS), int((BURSTS[4][1] + 0.3e-3) * FS)
        x[a:b] += w[a:b]
    return x.astype(np.complex64)


def _seeds():
    """What a weak structure test gives: every burst split in two, starting late."""
    out = []
    for a, b in BURSTS:
        m = (a + b) / 2
        out += [(a + 0.2e-3, m - 0.15e-3), (m + 0.1e-3, b - 0.1e-3)]
    return out


def test_rect_fit_ignores_dips():
    d = np.r_[np.zeros(40), 3 + np.tile([0, 0, -2.5, 0, 0], 12), np.zeros(40)]
    lo, hi, inn, out = BU.rect_fit(d)
    assert (lo, hi) == (40, 100)


def test_band_and_bursts_of_a_weak_link():
    tf = BU.TFPower(_capture(0.0), FS)
    seeds = _seeds()
    fb = BU.fit_band(tf, seeds, 0.0, 5e6, -5e6, 12e6)
    assert fb is not None
    assert abs(fb["lo"] - LO) < 0.15e6 and abs(fb["hi"] - HI) < 0.15e6
    runs, info = BU.energy_bursts(tf, seeds, fb["lo"], fb["hi"], FS / 2 * 0.9, min_sym_s=35e-6)
    assert info["used"]
    assert len(runs) == len(BURSTS)                           # no split burst
    for (a, b), (ta, tb) in zip(runs, BURSTS):
        assert abs(a - ta) < 40e-6 and abs(b - tb) < 40e-6    # edges to a frame or two (16.7 us)


def test_noise_only_returns_the_seeds():
    tf = BU.TFPower(_capture(-60.0, wifi=False), FS)
    seeds = _seeds()
    runs, info = BU.energy_bursts(tf, seeds, LO, HI, FS / 2 * 0.9)
    assert not info["used"] and runs == sorted(seeds)
