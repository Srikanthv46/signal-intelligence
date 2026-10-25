"""Bell 202 AFSK packet radio (1200 bit/s, mark 1200 Hz / space 2200 Hz) on an FM voice channel: APRS and other
AX.25 UI frames (research, signal gaps S19; research/signal_gaps/06_vhf_below_300.md section 7).

    decode(d, fs) -> list[dict]     d: FM discriminator output (Hz) at fs

Tone energies over one bit (sliding sums of d x e^-j w n), decision mark vs space at 8 bit phases, NRZI decoding,
HDLC frames with a valid CRC-16/X.25 (core.hdlc), AX.25 UI address / control / PID checks, APRS fields. Depends on
numpy and spectrum_intel.core only.
"""
from __future__ import annotations

import numpy as np

from spectrum_intel.core import hdlc as H

BAUD = 1200.0
MARK, SPACE = 1200.0, 2200.0
PHASES = 8
MIN_S = 0.12                       # shortest frame worth searching: ~18 octets + flags at 1200 bit/s


def _tone_energy(d: np.ndarray, fs: float, f: float, w: int) -> np.ndarray:
    """|sum over the last w samples of d e^-j2pi f n / fs|^2 at every sample."""
    t = np.arange(len(d)) / fs
    z = d * np.exp(-2j * np.pi * f * t)
    c = np.concatenate([[0], np.cumsum(z)])
    s = c[w:] - c[:-w]
    return np.abs(s) ** 2


def decode(d: np.ndarray, fs: float) -> list[dict]:
    if len(d) / fs < MIN_S:
        return []
    d = np.asarray(d, np.float64)
    spb = fs / BAUD
    w = int(round(spb))
    soft = _tone_energy(d, fs, MARK, w) - _tone_energy(d, fs, SPACE, w)      # index i = window ending at i + w
    nb = int(len(soft) / spb) - 1
    out, seen = [], set()
    for ph in range(PHASES):
        idx = (ph * spb / PHASES + np.arange(nb) * spb).astype(np.int64)
        idx = idx[idx < len(soft)]
        lv = (soft[idx] > 0).astype(np.int8).tolist()
        for pos, octs in H.deframe(H.unnrzi(lv)):
            if octs in seen:
                continue
            f = H.ax25(octs)
            if f is None:
                continue
            seen.add(octs)
            rec = {"t_s": round(float(idx[pos] + w) / fs, 4), "raw_hex": octs.hex(), "src": f["src"],
                   "dest": f["dest"], "path": f["path"], "info": f["info"]}
            a = H.aprs(f)
            if a:
                rec["aprs"] = a
            out.append(rec)
    out.sort(key=lambda r: r["t_s"])
    return out
