"""802.11 OFDM at a non-standard clock, at any frequency, in either spectrum orientation.

Many links are a stock 802.11 chipset run on a different sample clock and frequency-converted out of
the Wi-Fi bands: the 5 and 10 MHz channels of 802.11 itself (quarter and half clocking), vendor data
links (Microhard at 11 MHz in 2025-2110 MHz, DJI Wi-Fi links at 5 MHz), and up/down-converters
that mirror the spectrum. The waveform is still 802.11: 64-point OFDM, a cyclic prefix of a quarter
of the symbol, the same preamble and codes. Only the clock (so every duration and bandwidth), the
centre frequency and possibly the sign of frequency differ. So instead of a new receiver this
module measures the numerology and then runs the standard receiver on a time-scaled capture.

    find_bands(x, fs, fc)              -> (bands, unions of adjacent fragments): bursty wideband
                                          bands with sharp, repeatable edges
    numerology(x, fs, fc, band)        -> OFDM useful-symbol time, guard fraction, clock (or None)
    scan(x, fs, fc)                    -> [Link]: detector state + PPDUs in real units

Steps
  1. Bands: a ~60 kHz spectrogram; each frame's runs of bins at least 10 dB over a per-bin noise floor
     (20th percentile over time, morphologically opened so a busy bin cannot raise its own floor);
     runs of 1.2-45 MHz whose edges recur in many frames are candidate bursts of one emitter.
  2. Numerology: the band is cut out and gated to its active frames. A CP-OFDM signal correlates with
     itself one useful-symbol time T_u later (the cyclic prefix); the lag of that peak gives T_u, and the
     periodicity of the CP correlation gives the full symbol T_s. 802.11 has T_s/T_u = 1.25 (long
     guard) and 64 subcarriers, so the clock is 64/T_u = 80/T_s. The occupied width at -10 dB from the
     in-band level must fit the clock (52-56 of 64 subcarriers).
  3. Decode: the capture is treated as if sampled at fs * 20 MHz / clock, which makes the link a
     standard 20 MHz Wi-Fi channel, and receivers.wifi.ofdm_rx runs on it as usual (L-STF, L-LTF,
     L-SIG parity, HT-SIG CRC). The clock is tried at the round values (0.5 MHz grid) within 1.5 % of
     the measurement, then the measurement itself; both the capture and its complex conjugate
     (mirrored spectrum) are tried; the clock and orientation kept are the ones that yield the most
     CRC/parity-checked PPDUs. Payloads then go through
     payload.decode_nonht_payload / decode_ht_payload unchanged, because their state is in the
     time-scaled domain.

A standard 20 MHz channel is left to the standard scan (WifiDetector's channel plan); this module
reports 20 MHz only when it sits off the channel grid.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.ndimage import maximum_filter1d, minimum_filter1d
from scipy.signal import firwin, oaconvolve

from . import ofdm_rx as W

TU_RANGE_S = (2.4e-6, 26e-6)          # clocks of 2.5-26.7 MHz
MIN_WIDTH_HZ, MAX_WIDTH_HZ = 1.2e6, 45e6


@dataclass
class Link:
    center_hz: float
    clock_hz: float
    inverted: bool
    numerology: dict
    packets: list = field(default_factory=list)      # ofdm_rx.Packet, converted to real units
    n_rejects: int = 0

    @property
    def scale(self) -> float:
        return 20e6 / self.clock_hz


# ----------------------------------------------------------------------------------------- bands
def _spectrogram(x, fs, res_hz=60e3):
    n = 1 << int(round(np.log2(fs / res_hz)))
    n = max(64, min(n, 4096))
    m = len(x) // n
    X = np.fft.fftshift(np.fft.fft(x[:m * n].reshape(m, n) * np.hanning(n), axis=1), axes=1)
    return np.abs(X) ** 2, n


def _runs(mask):
    d = np.diff(np.concatenate([[0], mask.astype(np.int8), [0]]))
    return np.flatnonzero(d == 1), np.flatnonzero(d == -1)


def find_bands(x, fs, fc, thresh_db=10.0, min_frames=None):
    """Bursty emitters with sharp, repeatable edges -> ([dict(f_lo, f_hi, n_frames, frames)], unions) (Hz, RF)."""
    P, n = _spectrogram(x, fs)
    bin_hz = fs / n
    floor = np.percentile(P, 20, axis=0)
    k = max(3, int(3e6 / bin_hz))
    floor = maximum_filter1d(minimum_filter1d(floor, k, mode="nearest"), k, mode="nearest")
    occ = P > floor * 10 ** (thresh_db / 10)
    occ[:, :max(1, n // 50)] = False                     # analog band edges of the capture
    occ[:, -max(1, n // 50):] = False
    # bridge 1-2 bin holes (the DC subcarrier, a deep fade)
    occ = occ | (np.roll(occ, 1, 1) & np.roll(occ, -1, 1)) | (np.roll(occ, 2, 1) & np.roll(occ, -2, 1))
    lo_w, hi_w = MIN_WIDTH_HZ / bin_hz, MAX_WIDTH_HZ / bin_hz
    edges = {}
    for t in range(occ.shape[0]):
        a, b = _runs(occ[t])
        for s, e in zip(a, b):
            if lo_w <= e - s <= hi_w:
                edges.setdefault((s // 2, e // 2), []).append((t, s, e))
    if min_frames is None:
        min_frames = max(4, int(150e-6 * fs / n))         # >= 150 us of activity
    # merge neighbouring edge cells (within 3 bins) into clusters
    cells = sorted(edges.items(), key=lambda kv: -len(kv[1]))
    used, bands = set(), []
    for key, hits in cells:
        if key in used:
            continue
        group = []
        for k2, h2 in cells:
            if k2 not in used and abs(k2[0] - key[0]) <= 2 and abs(k2[1] - key[1]) <= 2:
                group += h2
                used.add(k2)
        if len(group) < min_frames:
            continue
        s = int(np.median([g[1] for g in group]))
        e = int(np.median([g[2] for g in group]))
        frames = np.zeros(occ.shape[0], bool)
        frames[[g[0] for g in group]] = True
        bands.append(dict(f_lo=fc - fs / 2 + s * bin_hz, f_hi=fc - fs / 2 + e * bin_hz, n_frames=len(group),
                          frames=frames, frame_len=n))
    bands.sort(key=lambda b: -b["n_frames"])
    # A burst whose spectrum dips (a fade, a narrow interferer) splits into adjacent runs that are active
    # in the same frames: offer their union too, after the bands themselves.
    unions = []
    for i, a in enumerate(bands):
        for b in bands[i + 1:]:
            lo, hi = sorted((a, b), key=lambda z: z["f_lo"])
            both = np.sum(a["frames"] & b["frames"])
            if (-0.3e6 < hi["f_lo"] - lo["f_hi"] < 0.6e6 and hi["f_hi"] - lo["f_lo"] <= MAX_WIDTH_HZ
                    and both >= 0.5 * min(a["n_frames"], b["n_frames"])):
                fr = a["frames"] | b["frames"]
                unions.append(dict(f_lo=lo["f_lo"], f_hi=hi["f_hi"], n_frames=int(fr.sum()), frames=fr,
                                   frame_len=n, merged=True))
    unions.sort(key=lambda b: -b["n_frames"])
    return bands, unions


# ----------------------------------------------------------------------------------------- numerology
def _snap(clock):
    g = round(clock / 0.5e6) * 0.5e6
    return g if g and abs(g - clock) / clock < 5e-3 else round(clock, -3)


def numerology(x, fs, fc, band):
    """CP-OFDM useful-symbol time, guard fraction and 802.11 clock of one band, or None."""
    cen = (band["f_lo"] + band["f_hi"]) / 2
    bw = band["f_hi"] - band["f_lo"]
    dec = max(1, int(fs // (1.6 * bw)))
    fsd = fs / dec
    h = firwin(129, 0.55 * bw, fs=fs)
    nn = np.arange(len(x))
    y = oaconvolve(x * np.exp(-2j * np.pi * (cen - fc) / fs * nn), h, "same")[::dec]
    # Gate to the band's active frames, closed over gaps of up to 3 frames: a gate that flickers frame
    # by frame would imprint the spectrogram frame rate on the CP correlation below.
    fr_on = band["frames"].astype(float)
    fr_on = np.convolve(fr_on, np.ones(7), "same") > 0
    gate = np.repeat(fr_on, band["frame_len"])[:len(x)][::dec]
    y = y[:len(gate)] * gate[:len(y)]
    e = np.sum(np.abs(y) ** 2)
    if e <= 0:
        return None
    # Occupied band at -10 dB from the band's own in-band level (a strong burst's skirts stand far
    # above the noise floor, so the floor-relative edges of find_bands overstate the width).
    nseg = 256
    on = y[np.abs(y) > 0]
    if len(on) >= 4 * nseg:
        segs = on[:len(on) // nseg * nseg].reshape(-1, nseg) * np.hanning(nseg)
        psd = np.fft.fftshift(np.mean(np.abs(np.fft.fft(segs, axis=1)) ** 2, 0))
        f = np.fft.fftshift(np.fft.fftfreq(nseg, 1 / fsd))
        inner = np.abs(f) < 0.3 * bw
        lvl = np.median(psd[inner])
        hot = psd > lvl / 10
        c0 = int(np.argmin(np.abs(f)))
        a, b = c0, c0
        while a > 0 and (hot[a - 1] or (a > 1 and hot[a - 2])):
            a -= 1
        while b < nseg - 1 and (hot[b + 1] or (b < nseg - 2 and hot[b + 2])):
            b += 1
        lo, hi = f[a] - 0.5 * fsd / nseg, f[b] + 0.5 * fsd / nseg
        if hi - lo > 0.3 * bw:
            cen, bw = cen + (lo + hi) / 2, hi - lo
    lmin, lmax = int(TU_RANGE_S[0] * fsd), int(np.ceil(TU_RANGE_S[1] * fsd))
    lmax = min(lmax, len(y) // 4)
    if lmax <= lmin + 4:
        return None
    nf = 1 << int(np.ceil(np.log2(2 * len(y))))
    F = np.fft.fft(y, nf)
    r = np.abs(np.fft.ifft(np.abs(F) ** 2)[:lmax + 2]) / e
    seg = r[lmin:lmax + 1]
    i = int(np.argmax(seg)) + lmin
    peak, base = r[i], float(np.median(seg))
    if peak < 0.06 or peak < 5 * base:
        return None
    a, b, c = r[i - 1], r[i], r[i + 1]
    lag = i + 0.5 * (a - c) / (a - 2 * b + c) if (a - 2 * b + c) else i
    tu = lag / fsd
    # symbol period from the periodicity of the CP correlation (window = a quarter symbol)
    L = int(round(lag))
    cc = y[:-L] * np.conj(y[L:])
    w = max(2, L // 4)
    m = np.abs(np.convolve(cc, np.ones(w), "valid"))
    m = m - m.mean()
    nf = 1 << int(np.ceil(np.log2(len(m) * 4)))
    S = np.abs(np.fft.rfft(m, nf))
    fr = np.fft.rfftfreq(nf, 1 / fsd)
    sel = (fr > 0.6 / tu) & (fr < 0.97 / tu)
    frame_rate = fs / band["frame_len"]
    harm = np.round(fr / frame_rate)
    sel &= ~((harm >= 1) & (np.abs(fr - harm * frame_rate) < 0.01 * frame_rate))
    if not sel.any():
        return None
    j = np.flatnonzero(sel)[int(np.argmax(S[sel]))]
    line = S[j] / (np.median(S[sel]) + 1e-20)
    a, b, c = S[j - 1], S[j], S[j + 1]
    fj = fr[j] + (0.5 * (a - c) / (a - 2 * b + c) if (a - 2 * b + c) else 0.0) * (fr[1] - fr[0])
    ts = 1 / fj
    gi = ts / tu - 1
    out = dict(center_hz=float(cen), occupied_hz=float(bw), tu_s=float(tu), ts_s=float(ts), guard=round(float(gi), 3),
               cp_corr=round(float(peak), 3), cp_corr_over_median=round(float(peak / max(base, 1e-9)), 1),
               symbol_line_over_median=round(float(line), 1))
    if not (0.2 <= gi <= 0.3 and line > 3):
        out["reject"] = f"guard {gi:.3f} of the symbol (802.11 long guard is 0.25)"
        return out
    clock = 80 / ts if abs(gi - 0.25) < 0.02 else 64 / tu
    out["clock_hz"] = float(_snap(clock))
    out["clock_measured_hz"] = float(clock)
    ratio = bw / out["clock_hz"]
    out["occupied_over_clock"] = round(float(ratio), 3)
    if not 0.7 <= ratio <= 1.0:
        out["reject"] = f"occupied width {bw / 1e6:.2f} MHz does not fit 52-56 of 64 subcarriers"
    return out


# ----------------------------------------------------------------------------------------- decode
def _decode(x, fs, fc, nm, inverted):
    k = 20e6 / nm["clock_hz"]
    fsv = fs * k
    off = (nm["center_hz"] - fc) * k * (-1 if inverted else 1)
    chans = [dict(name="scaled", ch=None, freq_hz=off + d, offset_hz=off + d, usable_half_hz=fsv / 2, partial=False)
             for d in (-1.0e6, 0.0, 1.0e6)]                    # the STF CFO range is +-625 kHz (time-scaled)
    # A little more sensitive than the standard scan (the band is already known to hold 802.11-numerology
    # OFDM); every PPDU still needs a parity-checked, bit-agreeing L-SIG or a CRC-valid HT-SIG.
    det = W.WifiDetector(fsv, 0.0, channels=chans, usable_bw=fsv, stf_thresh=0.6, power_margin_db=2.5,
                         min_plateau=40)
    return det.run(np.conj(x) if inverted else x), det


def _evidence(pk):
    return sum(2 if p.sig_crc_ok else 1 for p in pk)


def _to_real(p, fs, fc, k, inverted, clock):
    """A Packet decoded in the time-scaled domain -> times, frequencies and rates of the real capture."""
    sgn = -1 if inverted else 1
    cen = fc + sgn * p.center_mhz * 1e6 / k
    half = (8.75e6 if p.fmt == "HT" else 8.3e6) / k
    p.extra.update(scaled_center_mhz=p.center_mhz, scaled_cfo_khz=p.cfo_khz, clock_mhz=clock / 1e6,
                   time_scale=round(k, 6), spectrum_inverted=bool(inverted),
                   lsig_rate_code_mbps=p.lsig_rate_mbps)
    p.center_mhz = round(cen / 1e6, 4)
    p.cfo_khz = round(sgn * p.cfo_khz / k, 3)
    p.start_us = round(p.start_sample / fs * 1e6, 3)
    p.duration_us = round(p.duration_us * k, 2)
    p.end_us = round(p.start_us + p.duration_us, 3)
    if p.energy_end_us is not None:
        p.energy_end_us = round(p.energy_end_us * k, 2)
    p.f_lo_mhz = round(max(cen - half, fc - fs / 2) / 1e6, 4)
    p.f_hi_mhz = round(min(cen + half, fc + fs / 2) / 1e6, 4)
    p.bandwidth_mhz = clock / 1e6
    p.channel = f"{cen / 1e6:.2f}MHz@{clock / 1e6:g}MHz"
    return p


def _standard(nm, fc, fs):
    """A 20 MHz channel on the standard grid: the regular scan owns it."""
    if abs(nm["clock_hz"] - 20e6) > 0.06e6:
        return False
    return any(abs(f - nm["center_hz"]) < 1e6 for _, _, f in W._all_channels())


def scan(x, fs, fc, max_bands=12, max_unions=4):
    """Every 802.11-numerology link in the capture that the standard channel scan cannot see."""
    if fs < 4e6:
        return []
    links, seen = [], []
    bands, unions = find_bands(x, fs, fc)
    for band in bands[:max_bands] + unions[:max_unions]:
        nm = numerology(x, fs, fc, band)
        if not nm or nm.get("reject") or "clock_hz" not in nm:
            continue
        if _standard(nm, fc, fs) or (abs(nm["clock_measured_hz"] - 20e6) < 0.3e6 and _standard(
                dict(nm, clock_hz=20e6), fc, fs)):
            continue
        if any(abs(nm["center_hz"] - c) < 0.25 * nm["clock_hz"] and abs(nm["clock_hz"] - k) < 0.01 * k
               for c, k in seen):
            continue
        seen.append((nm["center_hz"], nm["clock_hz"]))
        # The measured clock is good to ~1 % on a few bursts; the decode decides between it and the
        # round values (0.5 MHz grid) within 1.5 %, nearest first.
        meas = nm["clock_measured_hz"]
        grid = [g * 0.5e6 for g in range(int(meas * 0.985 / 0.5e6), int(meas * 1.015 / 0.5e6) + 2)
                if abs(g * 0.5e6 - meas) <= 0.015 * meas]
        clocks = sorted(grid, key=lambda g: abs(g - meas)) + [round(meas, -3)]
        best = None
        for clock in dict.fromkeys(clocks):
            nm["clock_hz"] = clock
            for inv in (False, True):
                pk, det = _decode(x, fs, fc, nm, inv)
                if pk and (best is None or _evidence(pk) > _evidence(best[0])):
                    best = (pk, det, inv, clock)
        if best is None:
            continue
        pk, det, inv, nm["clock_hz"] = best
        k = 20e6 / nm["clock_hz"]
        pk = [_to_real(p, fs, fc, k, inv, nm["clock_hz"]) for p in pk]
        cen = float(np.median([p.center_mhz for p in pk])) * 1e6
        ln = Link(center_hz=cen, clock_hz=nm["clock_hz"], inverted=inv, numerology=nm, packets=pk,
                  n_rejects=len(det.rejects))
        dup = next((o for o in links if abs(o.center_hz - cen) < 0.25 * ln.clock_hz and o.clock_hz == ln.clock_hz), None)
        if dup is None:
            links.append(ln)
        elif _evidence(pk) > _evidence(dup.packets):          # two band candidates of one link: keep the better
            links[links.index(dup)] = ln
    return links
