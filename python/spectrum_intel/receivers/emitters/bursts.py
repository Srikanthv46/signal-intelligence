#!/usr/bin/env python3
"""
assoc_bursts.py -- self-contained burst source for the emitter-association
experiments (drone_assoc.py). Two things:

  extract_bursts(x, fs, fc)   15 kHz-SCS OFDM bursts (DJI OcuSync-class RC
                              uplink hops, DroneID, video) from a capture, as
                              drone_common.Emission records. Uses the
                              spectrum_survey scan/refine/track functions
                              read-only (no detector modules involved).

  synth_scene(radios, ...)    ground-truth scenes: N independent OcuSync-like
                              uplink radios, each with its own crystal (one
                              ppm error applied to carrier AND frame clock),
                              its own frame phase, hopping on a shared channel
                              raster with the 4 + 6 ms uplink sub-pattern.
                              Returns IQ plus the true per-burst (radio, t, f).

  inject_copy(x, fs, fc, ...) real-signal injection: copy real bursts out of a
                              capture, re-time them on an independent frame
                              phase and re-tune them with a different ppm,
                              add back -> a realistic "second drone" made of
                              real waveforms (same model, different radio).

Usage
  python assoc_bursts.py capture.sigmf-meta      # list the bursts found
"""
from __future__ import annotations
import sys
import numpy as np

from spectrum_intel.core import spectral as SP, ofdm15k as O15
from spectrum_intel.core.emission import Emission, load_capture

# ------------------------------------------------------------------ extraction
_LABELS = {"ocusync_uplink_like": ("dji_ocusync", "ocusync_uplink"),
           "droneid_candidate": ("dji_droneid", "droneid"),
           "ocusync_video_like": ("dji_ocusync", "ocusync_video"),
           "ofdm15k_other": ("unknown_ofdm", "ofdm15k_other")}


def extract_bursts(x, fs, fc, usable_half=None, ch=None, return_ch=False):
    """15 kHz OFDM bursts -> list[Emission] (sorted by start). Features carry what the
    scan measures (coarse); drone_assoc.measure_fingerprints refines them."""
    usable_half = usable_half or min(fs, 56e6) / 2
    spurs = SP.find_spurs(x, fs, fc)
    if ch is None:
        ch = SP.Channelizer(x.astype(np.complex128), fs, [(s["offset_hz"], 30e3) for s in spurs])
    hits, stats = O15.ofdm15k_scan(ch, fc, usable_half=usable_half)
    ob = O15.merge_ofdm15k([b for b in (O15.refine_ofdm15k(ch, h) for h in hits) if b is not None])
    ob = O15.merge_ofdm15k(ob + O15.track_hops(ch, fc, ob, stats, len(x) / fs, usable_half=usable_half))
    out = []
    for b in sorted(ob, key=lambda b: b["start_s"]):
        lab, why = O15.label_ofdm15k(b)
        fam, label = _LABELS[lab]
        cen = (fc + b["f_off"]) / 1e6
        out.append(Emission(fam, label, b["start_s"] * 1e6, b["end_s"] * 1e6,
                            cen - b["bw_hz"] / 2e6, cen + b["bw_hz"] / 2e6, cen, None,
                            confidence=min(0.95, 0.5 + 0.5 * b["score"]),
                            features=dict(scs_khz=15.0, bw_mhz=b["bw_hz"] / 1e6, symbol_us=b["period_us"],
                                          burst_us=(b["end_s"] - b["start_s"]) * 1e6, nsym=b["nsym"],
                                          profile=b["profile"], cp_score=b["score"],
                                          tracked=bool(b.get("tracked", False)), source="assoc_bursts")))
    return (out, ch, spurs) if return_ch else out


# ------------------------------------------------------------------ synthesis
FS_SYN = 61.44e6
N_SYN = 4096                              # 15 kHz at 61.44 MSps
CP_SYN = (320, 288, 288, 288, 288, 288, 288)   # LTE normal-CP slot: 7 symbols = 500 us


def zc(u, n):
    k = np.arange(n)
    return np.exp(-1j * np.pi * u * k * (k + 1) / n)


def uplink_burst(rng, n_sc=72, ref_syms=(0, 6), zc_root=25, data=None):
    """One OcuSync-uplink-like 0.5 ms slot at 61.44 MSps, unit in-band power.
    Symbols in ref_syms carry the same Zadoff-Chu sequence (as observed on real
    DJI uplink hops, where symbols 0 and 6 repeat burst to burst); the rest QPSK."""
    k = np.r_[-n_sc // 2:0, 1:n_sc // 2 + 1]                         # DC unused
    ref = zc(zc_root, n_sc + 1)[:n_sc]
    out = []
    for i, cp in enumerate(CP_SYN):
        X = np.zeros(N_SYN, complex)
        X[k % N_SYN] = ref if i in ref_syms else (rng.choice([-1, 1], n_sc) + 1j * rng.choice([-1, 1], n_sc)) / np.sqrt(2)
        s = np.fft.ifft(X) * N_SYN / np.sqrt(n_sc)
        out.append(np.r_[s[-cp:], s])
    return np.concatenate(out)


def place(x, burst, t0_s, f_hz, fc, fs, gain, phase=0.0):
    """Add burst into x starting at fractional time t0_s, carrier f_hz (absolute)."""
    a = int(np.floor(t0_s * fs)); frac = t0_s * fs - a
    pad = 64
    b = np.r_[np.zeros(pad), burst, np.zeros(pad)]
    kk = np.fft.fftfreq(len(b))
    b = np.fft.ifft(np.fft.fft(b) * np.exp(-2j * np.pi * kk * frac))   # delay by frac samples
    s = a - pad
    lo, hi = max(s, 0), min(s + len(b), len(x))
    if hi <= lo:
        return
    n = np.arange(lo, hi)
    x[lo:hi] += gain * b[lo - s:hi - s] * np.exp(1j * (2 * np.pi * (f_hz - fc) / fs * n + phase))


def synth_scene(radios, dur_s=0.0325, fs=FS_SYN, fc=5770e6, rng=None, raster=(5725.5e6, 1e6),
                usable_half=26e6, noise=True):
    """radios: list of dict(ppm, phase_ms, snr_db, [T_ms=10], [offsets_ms=(0,4)], [channels], [zc_root]).
    Each radio's crystal error e (ppm) scales its carrier by (1+e) and its frame clock
    period by 1/(1+e). Channels hop pseudo-randomly over the in-band raster.
    Returns x, truth (list of dict(radio, t_us, f_hz, f_nom_hz))."""
    rng = rng or np.random.default_rng(1)
    n = int(dur_s * fs)
    x = ((rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)).astype(complex) if noise \
        else np.zeros(n, complex)
    f0, step = raster
    chans = f0 + step * np.arange(200)
    chans = chans[np.abs(chans - fc) < usable_half]
    truth = []
    for r_i, r in enumerate(radios):
        e = r["ppm"] * 1e-6
        T = r.get("T_ms", 10.0) * 1e-3 / (1 + e)
        offs = [o * 1e-3 / (1 + e) for o in r.get("offsets_ms", (0.0, 4.0))]
        my = np.array(r.get("channels", chans))
        g = np.sqrt(10 ** (r["snr_db"] / 10) * 1.08e6 / fs)            # in-band SNR over 1.08 MHz
        t_frame = (r["phase_ms"] * 1e-3) % T - T
        while t_frame < dur_s:
            for o in offs:
                t0 = t_frame + o
                if 0 <= t0 and t0 + 0.6e-3 < dur_s:
                    fn = float(rng.choice(my))
                    f = fn * (1 + e)
                    place(x, uplink_burst(rng, zc_root=r.get("zc_root", 25)), t0, f, fc, fs, g,
                          rng.uniform(0, 2 * np.pi))
                    truth.append(dict(radio=r_i, t_us=t0 * 1e6, f_hz=f, f_nom_hz=fn))
            t_frame += T
    truth.sort(key=lambda d: d["t_us"])
    return x.astype(np.complex64), truth


def inject_copy(x, fs, fc, emissions, dt_us, dppm, gain_db=0.0, df_hz=None, pad_us=30.0):
    """Copy every emission's waveform (time box +- pad, its band) and add it back
    dt_us later, with its carrier moved as if from a crystal dppm different:
    f -> f + f*dppm*1e-6 (+ df_hz, e.g. a raster jump to another channel). Time
    positions scale with the same crystal: copies of bursts that were Dt apart
    are Dt*(1-dppm) apart. Returns new x and the list of (orig_idx, new_t_us, new_f_mhz)."""
    y = x.astype(complex).copy()
    X = np.fft.fft(x.astype(complex))
    f = np.fft.fftfreq(len(x), 1 / fs)
    out = []
    tref = emissions[0].start_us if emissions else 0.0
    for i, e in enumerate(emissions):
        sel = np.abs(f - (e.center_mhz * 1e6 - fc)) <= (e.f_hi_mhz - e.f_lo_mhz) * 1e6 * 0.6
        Y = np.zeros_like(X); Y[sel] = X[sel]
        band = np.fft.ifft(Y)
        a = max(0, int((e.start_us - pad_us) * 1e-6 * fs)); b = min(len(x), int((e.end_us + pad_us) * 1e-6 * fs))
        seg = band[a:b] * np.hanning(b - a) ** 0.25                      # soft edges
        shift = (e.center_mhz * 1e6) * dppm * 1e-6 + (df_hz or 0.0)
        t_new = tref + (e.start_us - tref) * (1 - dppm * 1e-6) + dt_us
        t0 = a / fs + (t_new - e.start_us) * 1e-6
        tmp = np.zeros(len(x), complex)
        place(tmp, seg, t0, fc + shift, fc, fs, 10 ** (gain_db / 20))
        y += tmp
        out.append((i, t_new, e.center_mhz + shift / 1e6))
    return y.astype(np.complex64), out


def copy_emissions(emissions, injected):
    """Emission records for the copies made by inject_copy (known time/frequency),
    so association can be tested independently of the burst detector."""
    from dataclasses import replace
    out = []
    for i, t_new, f_new in injected:
        e = emissions[i]
        d = f_new - e.center_mhz
        out.append(replace(e, start_us=t_new, end_us=t_new + (e.end_us - e.start_us), f_lo_mhz=e.f_lo_mhz + d,
                           f_hi_mhz=e.f_hi_mhz + d, center_mhz=f_new, snr_db=None,
                           features={k: v for k, v in e.features.items()
                                     if k in ("scs_khz", "bw_mhz", "symbol_us", "burst_us", "nsym", "profile")}))
    return out


def fresh(emissions):
    """Copies of emissions with only detector-level features (drop refined ones)."""
    from dataclasses import replace
    keep = ("scs_khz", "bw_mhz", "symbol_us", "burst_us", "nsym", "profile")
    return [replace(e, snr_db=None, features={k: v for k, v in e.features.items() if k in keep}) for e in emissions]


if __name__ == "__main__":
    x, fs, fc, meta = load_capture(sys.argv[1])
    for i, e in enumerate(extract_bursts(x, fs, fc)):
        print(f"{i:3d} {e.family:12s} {e.label:15s} {e.start_us:10.1f}-{e.end_us:10.1f} us "
              f"{e.center_mhz:10.4f} MHz bw {e.features['bw_mhz']:.2f}")
