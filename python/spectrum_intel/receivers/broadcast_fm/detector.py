"""FM broadcast (88-108 MHz) carriers, stereo pilot and RDS (research, signal gaps S20;
research/signal_gaps/06_vhf_below_300.md section 8).

    detect(ch) -> list[Emission]      family "broadcast_fm", label "fm_broadcast"

1. Candidates: the 100 kHz raster 87.6-107.9 MHz inside the capture; a channel's power over +-100 kHz >= 10 dB over
   the median channel power density of the band.
2. Per candidate: +-120 kHz at 240 kSps, FM discriminator -> MPX (Hz). FM broadcast: peak deviation (99.5th
   percentile) 15-110 kHz. Stereo: the 19 kHz line of the MPX spectrum >= 10 dB over its +-0.5-2 kHz surroundings.
3. RDS (when the block holds >= 2 blocks, 45 ms): MPX shifted down by 3 x the measured pilot (57 kHz without one),
   +-2.4 kHz, BPSK phase from the squared signal, biphase matched filter at 1187.5 bit/s (bit timing searched over
   one bit), differential decoding, block sync on the offset words, error-free blocks only (rds.py).
Depends on numpy / scipy and spectrum_intel.core only.
"""
from __future__ import annotations

import numpy as np
from scipy import signal as ss

from spectrum_intel.core.dsp import Channelizer
from spectrum_intel.core.emission import Emission

from . import rds as R

BAND = (87.5e6, 108.0e6)
RASTER = [round(87.6 + 0.1 * k, 1) for k in range(204)]
MPX_FS = 240e3
CHAN_DB = 10.0
PILOT_DB = 10.0
RDS_RATE = 1187.5
MIN_RDS_S = 0.045


def _pow(ch: Channelizer, f0: float, half: float) -> float:
    k = ch._bins(f0 - ch.fc, half)
    return float(np.mean(np.abs(ch.X[k]).astype(np.float64) ** 2)) if len(k) else 0.0


def candidates(ch: Channelizer) -> list[tuple[float, float]]:
    """(channel MHz, dB over the band floor) of occupied raster channels."""
    lo, hi = ch.fc - 0.45 * ch.fs, ch.fc + 0.45 * ch.fs
    chans = [f for f in RASTER if lo + 120e3 < f * 1e6 < hi - 120e3]
    if not chans:
        return []
    sel = (np.abs(ch.f) < 0.45 * ch.fs) & (ch.f + ch.fc >= BAND[0]) & (ch.f + ch.fc <= BAND[1])
    if sel.sum() < 64:
        return []
    floor = float(np.median(np.abs(ch.X[sel]).astype(np.float64) ** 2)) / np.log(2)
    out = []
    for f in chans:
        p = _pow(ch, f * 1e6, 100e3)
        db = 10 * np.log10(p / floor + 1e-30)
        if db >= CHAN_DB:
            out.append((f, float(db)))
    # keep local maxima on the raster (a strong station spills into its +-100 kHz neighbours)
    keep = []
    for i, (f, db) in enumerate(out):
        nb = [d for g, d in out if abs(g - f) < 0.15 and g != f]
        if all(db >= d for d in nb):
            keep.append((f, db))
    return keep


def _line_db(spec: np.ndarray, fr: np.ndarray, f0: float) -> float:
    on = np.abs(fr - f0) <= 15
    bg = (np.abs(fr - f0) > 500) & (np.abs(fr - f0) < 2000)
    if not on.any() or not bg.any():
        return 0.0
    return float(10 * np.log10(spec[on].max() / (np.median(spec[bg]) + 1e-30) + 1e-30))


def rds_bits(mpx: np.ndarray, fs: float, f57: float) -> list[int]:
    """Differentially decoded RDS bits from the MPX (Hz) at fs."""
    t = np.arange(len(mpx)) / fs
    z = mpx * np.exp(-2j * np.pi * f57 * t)
    sos = ss.butter(4, 2400, fs=fs, output="sos")
    z = ss.sosfiltfilt(sos, z.real) + 1j * ss.sosfiltfilt(sos, z.imag)
    th = 0.5 * np.angle(np.sum(z ** 2))
    r = (z * np.exp(-1j * th)).real
    spb = fs / RDS_RATE
    nbits = int(len(r) / spb) - 1
    if nbits < 52:
        return []
    half = spb / 2
    best, best_e = 0.0, -1.0
    for off in np.arange(0, spb, spb / 32):              # bit timing: the biphase correlation peak
        e = 0.0
        for k in range(nbits):
            a = int(off + k * spb)
            m = int(off + k * spb + half)
            b = int(off + (k + 1) * spb)
            e += abs(r[a:m].sum() - r[m:b].sum())
        if e > best_e:
            best, best_e = float(off), e
    enc = []
    for k in range(nbits):
        a = int(best + k * spb)
        m = int(best + k * spb + half)
        b = int(best + (k + 1) * spb)
        enc.append(1 if r[a:m].sum() - r[m:b].sum() > 0 else 0)
    return [enc[k] ^ enc[k - 1] for k in range(1, len(enc))]


def measure(ch: Channelizer, f_mhz: float) -> dict | None:
    y, fs = ch.slice(f_mhz * 1e6, 120e3, MPX_FS)
    if len(y) < 1024:
        return None
    d = np.angle(y[1:] * np.conj(y[:-1])) * (fs / (2 * np.pi))
    mpx = d - float(np.mean(d))                      # the carrier offset removed
    dev = float(np.percentile(np.abs(mpx).astype(np.float64), 99.5))
    spec = np.abs(np.fft.rfft(mpx * np.hanning(len(mpx)))) ** 2
    fr = np.fft.rfftfreq(len(mpx), 1 / fs)
    pilot_db = _line_db(spec, fr, 19e3)
    stereo = pilot_db >= PILOT_DB
    f_pilot = 19e3
    if stereo:
        on = np.abs(fr - 19e3) <= 15
        f_pilot = float(fr[on][np.argmax(spec[on])])
    off = float(np.mean(d))
    out = {"carrier_mhz": round((f_mhz * 1e6 + off) / 1e6, 6), "deviation_khz": round(dev / 1e3, 1),
           "stereo": bool(stereo), "pilot_db": round(pilot_db, 1), "carrier_offset_hz": round(off)}
    if len(mpx) / fs >= MIN_RDS_S:
        bits = rds_bits(mpx, fs, 3 * f_pilot)
        blks = R.blocks(bits) if bits else []
        if blks:
            info = R.decode(R.groups(blks))
            info["blocks_ok"] = len(blks)
            info["bits"] = len(bits)
            out["rds"] = info
    return out


def detect(ch: Channelizer) -> list[Emission]:
    if not (ch.fc - ch.fs / 2 < BAND[1] and ch.fc + ch.fs / 2 > BAND[0]):
        return []
    dur_us = ch.L / ch.fs * 1e6
    out = []
    for f, db in candidates(ch):
        m = measure(ch, f)
        if m is None or not (15.0 <= m["deviation_khz"] <= 110.0):
            continue
        rds_ok = bool(m.get("rds", {}).get("pi"))
        conf = 0.7 + (0.1 if m["stereo"] else 0.0) + (0.15 if rds_ok else 0.0)
        feats = {"channel_mhz": f, "level_db": round(db, 1), **m}
        out.append(Emission(family="broadcast_fm", label="fm_broadcast", start_us=0.0, end_us=dur_us,
                            f_lo_mhz=f - 0.1, f_hi_mhz=f + 0.1, center_mhz=m["carrier_mhz"], snr_db=round(db, 1),
                            confidence=round(min(conf, 0.97), 2), decoded=m.get("rds"), features=feats))
    return out
