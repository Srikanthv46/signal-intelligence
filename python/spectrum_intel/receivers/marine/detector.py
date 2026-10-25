"""AIS receiver (research, signal gaps S15; research/signal_gaps/06_vhf_below_300.md section 5).

    detect(ch, ppm=0.0) -> list[Emission]     family "marine", label "ais"

Per AIS channel inside the capture (161.975 / 162.025 MHz, receiver ppm applied): +-10 kHz at 48 kSps (5 samples
per bit), FM discriminator, the slowly varying offset removed (a 32-bit moving average), integrated over one bit
(5-sample boxcar); for each of the 5 sample
phases the sign of the discriminator gives NRZ levels -> NRZI decoding -> HDLC frames (training, flag, de-stuffing,
closing flag) whose CRC-16/X.25 residue is 0xF0B8 (ais.frames); the frames of all phases are merged by start time.
Depends on numpy and spectrum_intel.core only.
"""
from __future__ import annotations

import numpy as np

from spectrum_intel.core.dsp import Channelizer
from spectrum_intel.core.emission import Emission

from . import ais as A

CHANNELS = {"AIS 1 (87B)": 161.975e6, "AIS 2 (88B)": 162.025e6}
RATE = 9600.0
FS = 48e3
HALF_BW = 10e3


def channel_frames(ch: Channelizer, f_meas: float) -> list[dict]:
    y, fs = ch.slice(f_meas, HALF_BW, FS)
    if len(y) < 512:
        return []
    d = np.angle(y[1:] * np.conj(y[:-1])).astype(np.float64)
    sps = fs / RATE
    w = int(round(32 * sps))
    if len(d) > w:                                   # carrier offset: the 32-bit moving average
        d = d - np.convolve(d, np.ones(w) / w, "same")
    d = np.convolve(d, np.ones(int(round(sps))) / round(sps), "same")    # integrate over one bit
    nb = int((len(d) - 1) / sps)
    seen, out = {}, []
    for ph in range(int(round(sps))):
        idx = (ph + np.arange(nb) * sps).astype(np.int64)
        idx = idx[idx < len(d)]
        lv = (d[idx] > 0).astype(np.int8).tolist()
        bits = A.unnrzi(lv)
        for pos, octs in A.frames(bits):
            t = float(idx[pos]) / fs
            key = octs.hex()
            if key in seen and abs(seen[key] - t) < 0.01:
                continue
            seen[key] = t
            nbits = 8 * len(octs)
            out.append({"t_s": round(t, 5), "octets": octs.hex(), "bits": nbits,
                        "end_s": round(t + (nbits + 40) / RATE, 5)})
    out.sort(key=lambda r: r["t_s"])
    return out


def detect(ch: Channelizer, ppm: float = 0.0) -> list[Emission]:
    lo, hi = ch.fc - 0.45 * ch.fs, ch.fc + 0.45 * ch.fs
    out = []
    for name, f in CHANNELS.items():
        f_meas = f / (1 + ppm * 1e-6)
        if not (lo + 15e3 < f_meas < hi - 15e3):
            continue
        for r in channel_frames(ch, f_meas):
            octs = bytes.fromhex(r["octets"])
            msg = A.fields(octs)
            out.append(Emission(family="marine", label="ais", start_us=r["t_s"] * 1e6, end_us=r["end_s"] * 1e6,
                                f_lo_mhz=(f - 12.5e3) / 1e6, f_hi_mhz=(f + 12.5e3) / 1e6, center_mhz=f / 1e6,
                                confidence=0.97, decoded={"crc_ok": True, **msg, "raw_hex": r["octets"]},
                                features={"channel": name, "modulation": "GMSK 9600 bit/s"}))
    return out
