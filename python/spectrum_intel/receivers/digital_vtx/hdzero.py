"""
hdzero.py -- HDZero (Divimath DM6300) digital FPV video link, 27 MHz mode.

    hdzero_links(sl, fs, fc, usable, noise_pps, n_samples) -> list[dict]

The PHY as measured on a real HDZero Freestyle V2 (testdata/live/LIVE_HDZERO_5658MHz, Raceband R1, 2026-10-22; the
numerology is not published):
  - CP-OFDM, 1024-point FFT at 27 MSps: Tu 37.926 us, SCS 26.367 kHz; CP 32 samples (1/32), symbol 39.11 us;
  - 889 used subcarriers (-444..444, the centre one is a pilot) = 23.44 MHz;
  - 26 continual pilots (CONTINUAL), scattered pilots every 32 carriers shifting 8 per symbol (4-symbol cycle);
  - continuous (no TDD gaps), one transmit antenna;
  - a block of ~22 near-identical preamble symbols once per video frame (426 symbols = 16.66 ms at 60 fps);
  - the data carriers hold the video samples directly (JSCC): no bit format to decode.

Search: on every HDZero channel centre inside the capture, the band is cut out at 27 MSps (`sl`, from the block FFT
the receiver already has), the fractional CFO is taken from the CP phase at the CP-fold peak, and the integer carrier
offset (Pluto / VTX ppm) from the continual-pilot coherence. A link is accepted when the 26 continual pilots are
coherent from symbol to symbol well above the data carriers next to them (PILOT_MIN, PILOT_CONTRAST): random data,
noise, Wi-Fi or another OFDM numerology does not give that. The 17 MHz narrow mode (RunCam Nano 90) is not measured
yet and is not searched.
"""
from __future__ import annotations

import numpy as np

FS = 27e6
N, CP = 1024, 32
TS = N + CP
HALF = 444                                   # used carriers -HALF..HALF
SCS = FS / N
CONTINUAL = np.array([-444, -400, -360, -312, -280, -240, -224, -211, -200, -168, -136, -96, -56, 0, 56, 104, 128,
                      160, 184, 216, 264, 280, 312, 364, 404, 444])
# data carriers beside the pilots (same channel, same fading): the coherence reference
_REF = np.array(sorted(set((CONTINUAL[:, None] + np.array([-2, 2])[None, :]).ravel()) - set(CONTINUAL)))
_REF = _REF[np.abs(_REF) <= HALF]
_REF = _REF[(_REF - 4) % 8 != 0]             # not on a scattered-pilot position (every 8th carrier over 4 symbols)
CHANNELS = {**{f"R{i + 1}": f for i, f in enumerate((5658, 5695, 5732, 5769, 5806, 5843, 5880, 5917))},
            "E1": 5705, "F1": 5740, "F2": 5760, "F4": 5800,
            **{f"L{i + 1}": 5362 + 37 * i for i in range(8)}}
MAX_PPM = 25.0                               # VTX + receiver clock error searched
PILOT_MIN = 0.2                              # mean continual-pilot coherence (real link 0.95 at 14 dB, 0.28 at -6 dB)
PILOT_CONTRAST = 4.0                         # ... over the coherence of the data carriers beside them
MIN_SYMBOLS = 24                             # ~1 ms
PREAMBLE_RHO = 0.3                           # preamble symbols repeat (|rho| 0.96 at 14 dB); data symbols ~0.035
FRAME_SYMS = (300, 800)                      # 30-90 fps
NW_RANGE, NW_TOL = (1000, 1048), 0.02         # symbol FFT sizes tried; allowed non-integer part


def _smooth(n):
    for p in (2, 3, 5, 7, 11, 13):
        while n % p == 0:
            n //= p
    return n == 1


def slice_plan(n_samples, fs):
    """(slice length, per-symbol FFT size nw): the band is cut out at nw / 1024 x 27 MSps, with a 13-smooth slice
    length (a fast IFFT) for which nw = 1024 x length / (n_samples x 27 MSps / fs) is an integer to within NW_TOL.
    nw = 1024 nearest. 2^n blocks at 61.44 / 30.72 MSps are exact (921600, 1024). Fallback: the exact rate."""
    m0 = n_samples * FS / fs
    best = None
    for nw in sorted(range(NW_RANGE[0], NW_RANGE[1] + 1), key=lambda v: abs(v - N)):
        if best is not None and abs(nw - N) > abs(best[1] - N):
            break
        m = m0 * nw / N
        for c in range(int(np.ceil(m * (1 - NW_TOL / N))), int(np.floor(m * (1 + NW_TOL / N))) + 1):
            if _smooth(c):
                best = (c, nw)
                break
    return best if best is not None else (int(round(m0)), N)


def _fold_phase(y, lag, tsf, n_sym):
    """CP metric (lag-Tu product summed over a CP) folded on the (fractional) symbol grid -> phase, sum, contrast."""
    c = y[lag:] * np.conj(y[:-lag])
    cs = np.concatenate(([0], np.cumsum(c)))
    m = cs[CP:] - cs[:-CP]                    # CP-length moving sum
    base = np.round(np.arange(n_sym) * tsf).astype(int)
    f = m[base[:, None] + np.arange(int(tsf))[None, :]].sum(0)
    ph = int(np.argmax(np.abs(f)))
    return ph, f[ph], float(np.abs(f[ph]) / (np.median(np.abs(f)) + 1e-30))


def _coherence(S, cols):
    q = S[1:, cols] * np.conj(S[:-1, cols])
    return np.abs(q.mean(0)) / (np.abs(q).mean(0) + 1e-30)


def _preambles(S, cols):
    """Start symbols of the per-frame preamble blocks: symbols that repeat one of the next four. cols(k) -> bins."""
    used = np.arange(-HALF, HALF + 1)
    keep = ~np.isin(used, CONTINUAL) & ((used - 4) % 8 != 0)
    X = S[:, cols(used[keep])]
    X = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-30)
    M = len(X)
    rep = np.zeros(M, bool)
    for d in range(1, 5):
        r = np.abs(np.sum(X[d:] * np.conj(X[:-d]), 1)) > PREAMBLE_RHO
        rep[:-d] |= r
        rep[d:] |= r
    idx = np.flatnonzero(rep)
    if not len(idx):
        return []
    blocks, s, p = [], idx[0], idx[0]
    for i in idx[1:]:
        if i - p > 6:
            blocks.append((s, p))
            s = i
        p = i
    blocks.append((s, p))
    return [(int(a), int(b)) for a, b in blocks if b - a >= 3]


def hdzero_links(sl, fs, fc, usable, n_samples, log=None):
    """sl(f_off, half_bw, out_fs) -> (y, rate). One dict per accepted HDZero link (p_band: mean in-band power of
    the slice, for the caller's SNR)."""
    out = []
    span = (HALF + 0.5) * SCS
    # latency: the slice length at exactly 27 MSps can have a large prime factor (2e6 samples at 61.44 MSps ->
    # 878906 = 2 x 439453: a 0.21 s IFFT); slice_plan picks a fast length whose symbol FFT size is an integer
    # (894348 -> 1042 samples, 15 ms). The symbol grid is tracked fractionally.
    m1, nw = slice_plan(n_samples, fs)
    out_fs = m1 * fs / n_samples
    for ch, f_mhz in CHANNELS.items():
        off = f_mhz * 1e6 - fc
        if abs(off) + span > usable:
            continue
        y, rate = sl(off, span + 0.15e6, out_fs)
        tsf = TS * rate / FS                       # symbol period in slice samples
        lag = nw                                   # useful part in slice samples = the per-symbol FFT length
        rho = nw / (N * rate / FS)                 # slice bins per HDZero subcarrier (1 to within NW_TOL / 1024)
        cols = lambda k, s=0: (np.round(np.asarray(k) * rho).astype(int) + s) % nw
        n_sym = int((len(y) - lag - 2 * TS) // tsf)
        if n_sym < MIN_SYMBOLS:
            continue
        ph, acc, contrast = _fold_phase(y, lag, tsf, n_sym)
        cfo_frac = float(np.angle(acc)) / (2 * np.pi * lag / rate)
        y = y * np.exp(-2j * np.pi * cfo_frac * np.arange(len(y)) / rate)
        k = np.arange(int((len(y) - ph - CP - nw) // tsf))
        exact = ph + CP // 2 + k * tsf                               # window centred in the CP
        starts = np.round(exact).astype(int)
        S = np.fft.fft(y[starts[:, None] + np.arange(nw)[None, :]], axis=1)
        # the rounded start is up to half a sample off the fractional grid: 1.3 rad on the edge carriers
        S *= np.exp(2j * np.pi * np.outer(exact - starts, np.fft.fftfreq(nw) * nw) / nw)
        smax = int(np.ceil(MAX_PPM * 1e-6 * f_mhz * 1e6 / SCS)) + 1
        best = None
        for s in range(-smax, smax + 1):
            cp_ = _coherence(S, cols(CONTINUAL, s)).mean()
            ref = _coherence(S, cols(_REF, s)).mean()
            if best is None or cp_ > best[1]:
                best = (s, float(cp_), float(ref))
        shift, coh, ref = best
        ok = coh >= PILOT_MIN and coh >= PILOT_CONTRAST * ref
        if log is not None:
            log(dict(channel=ch, fold_contrast=round(contrast, 2), pilot_coherence=round(coh, 3),
                     data_coherence=round(ref, 3), shift=shift, accepted=bool(ok)))
        if not ok:
            continue
        cfo = shift * rate / nw + cfo_frac
        p_band = float(np.mean(np.abs(y) ** 2))
        blocks = _preambles(S, lambda k: cols(k, shift))
        t_sym = tsf / rate
        t0 = (ph + CP) / rate                          # first CP-aligned symbol
        pre_ms = [round((t0 + a * t_sym) * 1e3, 3) for a, _ in blocks]
        frame_syms = None
        if len(blocks) >= 2:
            d = np.diff([a for a, _ in blocks])
            d = d[(d >= FRAME_SYMS[0]) & (d <= FRAME_SYMS[1])]
            if len(d):
                frame_syms = int(np.median(d))
        frame_ms = None if frame_syms is None else frame_syms * t_sym * 1e3
        out.append(dict(channel=ch, channel_mhz=float(f_mhz), off=off, cfo_hz=cfo, ppm=cfo / (f_mhz * 1e6) * 1e6,
                        lo=off + cfo - span, hi=off + cfo + span, pilot_coherence=coh, data_coherence=ref,
                        fold_contrast=contrast, p_band=p_band, half_bw=span, n_symbols=len(S), preamble_ms=pre_ms,
                        preamble_symbols=[b - a + 1 for a, b in blocks], frame_symbols=frame_syms,
                        frame_ms=frame_ms, t0=0.0, t1=n_samples / fs))
    # one link per band: an adjacent channel's slice overlaps a strong link (not with the 37 MHz raster, but
    # E1 / F1 / F2 / F4 sit inside the Raceband channels)
    out.sort(key=lambda d: -d["pilot_coherence"])
    keep = []
    for d in out:
        if all(abs(d["off"] + d["cfo_hz"] - k["off"] - k["cfo_hz"]) > 2 * span * 0.5 for k in keep):
            keep.append(d)
    return keep


def synth(n_sym, rng, pilot_gain=4 / 3, preamble=(30, 8), scattered_base=(12, 4, 28, 20)):
    """Tests / vectors: an HDZero-like link at 27 MSps (the measured layout: continual pilots, scattered pilots every
    32 carriers shifting 8 per symbol, a block of identical preamble symbols; QPSK stands in for the JSCC data).
    Pilot values are a fixed +-pilot_gain pattern (the real values are not known)."""
    used = np.arange(-HALF, HALF + 1)
    pval = pilot_gain * np.where((used * 7) % 3 == 0, -1.0, 1.0)
    qpsk = lambda n: (rng.choice([-1, 1], n) + 1j * rng.choice([-1, 1], n)) / np.sqrt(2)
    pre = qpsk(used.size)
    out = []
    for n in range(n_sym):
        v = pre.copy() if preamble[0] <= n < preamble[0] + preamble[1] else qpsk(used.size)
        sp = (used - scattered_base[n % 4]) % 32 == 0
        v[sp] = pval[sp]
        cp = np.isin(used, CONTINUAL)
        v[cp] = pval[cp]
        F = np.zeros(N, complex)
        F[used % N] = v
        s = np.fft.ifft(F) * np.sqrt(N)
        out.append(np.concatenate((s[-CP:], s)))
    return np.concatenate(out)
