"""Analog FM video (PAL / NTSC composite) -> a picture, one field per block.

A *streaming* demodulator for a VTX that ``detect()`` already found: build it
from the detection, then feed it consecutive (or non-contiguous) IQ blocks;
every block that holds a vertical sync gives one field: luma, and colour when the
VTX sends a colour burst.

    demod = FieldDemodulator.from_fields(signal_fields)    # report decoded.fields
    demod = FieldDemodulator.from_emission(emission)       # detect() output
    frame = demod.process(iq, fs, fc)   # -> VideoFrame(image uint8 [lines, 720], meta, rgb) or None

Per block:
  1. slice the VTX band and FM-discriminate it (Hz) from the low-passed lag-1 product (each
     phase step weighted by its power, so the FM clicks below threshold count little; colour uses
     the plain discriminator, whose subcarrier this filter would attenuate); flip so
     the sync tip is the minimum; low-pass to the luma bandwidth;
  2. line clock, as an H-PLL over the block (blocks need not be contiguous, so each is timed on
     its own): the sync phase from the line-folded waveform at the detector's rate, then every
     line's sync leading edge (placed by the 4.7 us sync box, timed by a short step filter) and a
     robust straight-line fit t_n = t0 + n T; a weak block's fit is weighted against the
     detector's rate (+-PRIOR_PPM);
  3. resample every line's active part onto a fixed grid -> raster;
  4. levels: black = back-porch (blanking) level, white = 99.5th percentile (of the 3x3-median
     field when it is noisy, step 6: FM clicks do not set it);
  5. vertical sync: broad-pulse lines sit at the sync level for most of the
     line; the field is the active lines after a VBI;
  6. noise = robust spread of the flat sync tip / (white - black) (meta "noise"); above NOISY
     (about 5 dB in-band SNR, where FM clicks dominate) the field is 3x3-median filtered.

Colour (``rgb``, NTSC and PAL), on the field's lines only:
  7. the discriminator output low-passed to 6 MHz (keeps the subcarrier) is sampled on a
     line-locked grid of 1024 per line and mixed down with the line-locked subcarrier
     (NTSC 227.5 x fl, PAL 283.7516 x fl), then low-passed to 1 MHz -> complex chroma;
  8. per line, the burst phasor (mean over the burst gate), averaged over the +-3 nearest lines of
     the same parity (as a TV's burst PLL: a weak VTX gives a noisy burst on every line), is the
     phase reference (the gate tried at BURST_SHIFT_US after the standard's, strongest kept: cameras
     differ); the
     phase reference; PAL: the reference is the mean of the line's and the previous line's
     (the +-45 deg swing cancels), the swing's sign gives the V switch, and U, V are averaged
     over two lines (delay line);
  9. levels from the sync depth, as a TV: white 100 IRE (NTSC, black at 7.5 IRE setup) /
     700 mV (PAL) above blanking; saturation from the burst (ACC): the burst stands for
     20 IRE / 150 mV. (The luma image's white, a percentile, is lifted by the chroma.) A picture
     whiter than SYNC_COMPRESSED x that white has a compressed sync (a weak VTX near the FM
     threshold): the depth is then taken from the picture's white;
 10. Y = the line minus the re-modulated chroma; YUV -> RGB; above NOISY, Y is 3x3- and U, V are
     3x5-median filtered;
 11. colour killer: rgb is None unless the field's burst, coherently averaged over the lines
     (lag-2 products, so slow phase drift and the PAL swing do not cancel it), is at least
     BURST_MIN times the same measure 1 MHz below the subcarrier, and the noise is at most
     NOISY_COLOR.

2 M samples at 61.44 MSps (32.5 ms) hold about two fields.
Ported from signal-vault's sigint demodulator; depends on spectrum_intel.core only.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.ndimage import median_filter
from scipy.signal import firwin, oaconvolve

from spectrum_intel.core.dsp import Channelizer

# per standard: active-video start after the sync leading edge and its length, back-porch
# window (black level), lines skipped after the VBI, active lines per field, nominal line rate
STANDARDS = {
    "ntsc": dict(active_us=(10.9, 52.6), porch_us=(8.0, 9.2), skip=12, lines=240, fl=15734.264),
    "pal": dict(active_us=(12.0, 52.0), porch_us=(8.5, 10.0), skip=17, lines=288, fl=15625.0),
}
WIDTH = 720
# colour: subcarrier / line rate, burst gate (us after the sync edge), burst amplitude as a fraction
# of the sync depth (NTSC 20 IRE / 40 IRE; PAL 150 mV / 300 mV), white and black above blanking as
# multiples of the sync depth (NTSC 100 IRE and 7.5 IRE setup / 40 IRE; PAL 700 mV / 300 mV)
CHROMA = {
    "ntsc": dict(k=227.5, burst_us=(5.6, 7.6), burst_per_sync=0.5, white_per_sync=2.5, black_per_sync=0.1875),
    "pal": dict(k=283.7516, burst_us=(5.9, 8.0), burst_per_sync=0.5, white_per_sync=7 / 3, black_per_sync=0.0),
}
LINE_GRID = 1024        # samples per line for the chroma
BURST_MIN = 4.0         # field burst / off-subcarrier amplitude (coherent over the field) for colour
BURST_SMOOTH = 3        # phase reference: mean burst of +-3 same-parity lines (a TV's burst PLL)
BURST_SHIFT_US = (0.0, 0.2, 0.4)   # burst gate offsets tried (us) after the standard gate
PRIOR_PPM = 2.0         # line rate: the detector's measurement, as a prior of this width (ppm)
SYNC_COMPRESSED = 1.3   # picture white / the sync depth's white above which the sync is taken as compressed
NOISY = 0.09            # sync-tip noise / (white - black) above which the field is 3x3-median filtered
NOISY_COLOR = 0.25      # ... above which colour is killed (synthetic scenes: grey is then closer to the truth)


@dataclass
class VideoFrame:
    image: np.ndarray                    # uint8 [lines, WIDTH], 0 = black (luma, chroma included)
    meta: dict = field(default_factory=dict)
    rgb: np.ndarray | None = None        # uint8 [lines, WIDTH, 3], None: no colour burst (or color=False)


def fm_discriminate(y: np.ndarray, fs: float) -> np.ndarray:
    """Instantaneous frequency (Hz) of a complex baseband signal."""
    return np.angle(y[1:] * np.conj(y[:-1])) * (fs / (2 * np.pi))


def fm_discriminate_weighted(y: np.ndarray, fs: float, cutoff_hz: float = 6e6, ntaps: int = 7) -> np.ndarray:
    """Instantaneous frequency (Hz) from the low-passed lag-1 product y[n] y*[n-1]: each sample's
    phase step is weighted by its power, so the low-amplitude samples where noise spins the phase
    (FM clicks below threshold) count little. Same output as fm_discriminate on a clean signal."""
    y = np.asarray(y, np.complex64)
    z = y[1:] * np.conj(y[:-1])
    h = firwin(ntaps, min(cutoff_hz / (fs / 2), 0.99)).astype(np.float32)
    c = ntaps // 2                                    # symmetric FIR, direct in complex64 (fastest)
    zp = np.pad(z, c)
    zf = h[c] * z
    for i in range(c):
        zf += h[i] * (zp[i:i + len(z)] + zp[2 * c - i:2 * c - i + len(z)])
    return np.angle(zf) * (fs / (2 * np.pi))


def lowpass(x: np.ndarray, fs: float, cutoff_hz: float, ntaps: int = 31) -> np.ndarray:
    h = firwin(ntaps, min(cutoff_hz / (fs / 2), 0.99))
    return oaconvolve(x, h, mode="same")


class FieldDemodulator:
    def __init__(self, center_hz: float, line_rate_hz: float, standard: str = "ntsc",
                 sync_hz: float | None = None, bw_hz: float | None = None, polarity: str = "sync_low",
                 color: bool = True):
        self.std = "pal" if str(standard).lower() == "pal" else "ntsc"
        self.p = STANDARDS[self.std]
        self.center_hz = float(center_hz)
        self.sync_hz = float(sync_hz if sync_hz is not None else center_hz)
        self.half_bw = max((bw_hz or 12e6) / 2, abs(self.sync_hz - self.center_hz)) + 1.5e6
        self.rate0 = self.rate = float(line_rate_hz)   # rate0: the detector's measurement (sub-ppm)
        self.polarity = -1.0 if polarity == "sync_high" else 1.0
        self.color = bool(color)
        self._timing = None                            # set by raster(), used by chroma()

    @classmethod
    def from_fields(cls, f: dict) -> "FieldDemodulator":
        """From a report signal's decoded.fields (analog_vtx plugin)."""
        bw = f.get("bw_mhz")
        return cls(center_hz=float(f["center_mhz"]) * 1e6, line_rate_hz=float(f["line_rate_hz"]),
                   standard=f.get("standard", "ntsc"),
                   sync_hz=float(f["sync_tip_mhz"]) * 1e6 if f.get("sync_tip_mhz") else None,
                   bw_hz=bw * 1e6 if bw else None, polarity=f.get("polarity", "sync_low"))

    @classmethod
    def from_emission(cls, e) -> "FieldDemodulator":
        """From a detect() Emission."""
        return cls.from_fields({**e.features, "center_mhz": e.features.get("center_mhz", e.center_mhz)})

    def in_band(self, fs: float, fc: float) -> bool:
        return abs(self.center_hz - fc) + 2e6 < fs / 2

    # ------------------------------------------------------------------
    def process(self, iq: np.ndarray, fs: float, fc: float, ch: Channelizer | None = None) -> VideoFrame | None:
        if not self.in_band(fs, fc):
            return None
        hb = min(self.half_bw, fs / 2 - abs(self.center_hz - fc))
        ch = ch or Channelizer(np.asarray(iq, np.complex64), fs, fc)
        y, ofs = ch.slice(self.center_hz, hb, min(fs, max(2.4 * hb, 8e6)))
        d = lowpass(fm_discriminate_weighted(y, ofs) * self.polarity, ofs, 4.8e6, 31).astype(np.float32)
        img, meta = self.raster(d, ofs)
        if img is None:
            return None
        meta.update(standard=self.std, center_mhz=round(self.center_hz / 1e6, 3))
        rgb = None
        if self.color:
            # chroma from the plain discriminator: the weighted one's short filter attenuates the
            # subcarrier against the noise (and chroma is killed where clicks matter, NOISY_COLOR)
            d0 = fm_discriminate(y, ofs) * self.polarity
            rgb, cmeta = self.chroma(lowpass(d0, ofs, 6.0e6, 63).astype(np.float32), ofs)
            meta.update(cmeta)
        return VideoFrame(img, meta, rgb)

    def raster(self, d: np.ndarray, ofs: float):
        T_nom = 1 / self.rate
        if len(d) / ofs < 4 * T_nom:
            return None, {}
        # coarse copy (~4 MSps) for the sync search
        k = max(int(ofs // 4e6), 1)
        dc = d[: len(d) // k * k].reshape(-1, k).mean(1)
        cfs = ofs / k
        nbf = 512
        # sync phase from the fold at the detector's rate (sync = lowest frequency after the flip),
        # then the line clock from every line's own sync edge (an H-PLL over the block)
        T = 1 / self.rate0
        prof = self._fold(dc, cfs, self.rate0, nbf)
        w = max(int(round(4.7e-6 / T * nbf)), 1)
        box = np.convolve(np.r_[prof, prof[:w]], np.ones(w) / w, "valid")[:nbf]
        t0 = int(np.argmin(box)) / nbf * T
        tm = None
        for _ in range(2):
            tm = self._edge_fit(d, ofs, t0, T)
            if tm is None:
                break
            t0, T = tm[0], tm[1]
        if tm is not None:
            t0, T = self._prior_blend(t0, T, tm[2], len(d) / ofs)
        self.rate = 1 / T
        prof = self._fold(dc, cfs, self.rate, nbf)
        ph_us = ((np.arange(nbf) + 0.5) / nbf * T - (t0 - (k - 1) / 2 / ofs)) % T * 1e6
        a, b = self.p["porch_us"]
        sync = float(np.median(prof[(ph_us > 0.8) & (ph_us < 4.0)]))
        black = float(np.median(prof[(ph_us > a) & (ph_us < b)]))

        n_lines = int((len(d) / ofs - t0) / T) - 1
        if n_lines < 20:
            return None, {}
        a_us, len_us = self.p["active_us"]
        u = (a_us + np.arange(WIDTH) * (len_us / WIDTH)) * 1e-6
        t = t0 + np.arange(n_lines)[:, None] * T + u[None, :]
        raster = np.interp(t * ofs, np.arange(len(d)), d)
        # broad-pulse (VBI) lines: most of the line at the sync level
        mid = 0.5 * (sync + black)
        vbi = (raster < mid).mean(1) > 0.3
        rows = self._field_rows(vbi, n_lines)
        # noise: robust spread of the (flat) sync tip; below about 5 dB in-band SNR (noise > NOISY of
        # the black-white span) FM clicks dominate and a 3x3 median gives the better picture
        noise = self._sync_noise(d, ofs, t0, T, n_lines)
        field = raster[rows]
        act = ~vbi[rows] if (~vbi[rows]).any() else np.ones(len(rows), bool)
        white = float(np.percentile(field[act], 99.5))
        if white - black < 1e3:
            return None, {}
        rel = noise / (white - black)
        if rel > NOISY:
            field = median_filter(field, 3, mode="nearest")
            white = float(np.percentile(field[act], 99.5))        # FM clicks do not set white
            if white - black < 1e3:
                return None, {}
            rel = noise / (white - black)
        self._timing = dict(t0=t0, T=T, sync=sync, black=black, white=white, rows=rows, noise=rel)
        img = np.clip((field - black) / (white - black), 0, 1)
        img = (img * 255).astype(np.uint8)
        ppm = (self.rate / self.p["fl"] - 1) * 1e6
        meta = dict(line_rate_hz=round(self.rate, 3), line_rate_ppm=round(ppm, 1),
                    vsync=bool(vbi.any()), lines=int(len(rows)), noise=round(float(rel), 3),
                    denoise=bool(rel > NOISY))
        if tm is not None:
            meta["sync_jitter_us"] = round(tm[3], 3)
        return img, meta

    def _edge_fit(self, d, ofs, t0, T, win_us=3.0):
        """Every line's sync leading edge -> robust straight-line fit t_n = t0 + n T.
        Each edge is first placed by the 4.7 us box minimum (the sync is the lowest level, so a
        bright picture cannot pull it), then timed by a short step filter within +-0.6 us of it
        (its 0.5 us front window stays inside the 1.5 us front porch). Equalising and broad pulses
        start at the line start too. -> (t0, T, slope std error s, edge jitter us) or None."""
        c = np.r_[0.0, np.cumsum(d, dtype=np.float64)]
        sw, a, b = (int(round(v * 1e-6 * ofs)) for v in (4.7, 0.5, 2.0))
        w, r = int(round(win_us * 1e-6 * ofs)), int(round(0.6e-6 * ofs))
        n_lines = int((len(d) / ofs - t0) / T) - 1
        k = np.arange(max(n_lines, 0))
        cen = np.round((t0 + k * T) * ofs).astype(np.int64)
        ok = (cen - w - r - a >= 0) & (cen + w + r + sw < len(d))
        k, cen = k[ok], cen[ok]
        if len(k) < 20:
            return None
        i = cen[:, None] + np.arange(-w, w + 1)[None, :]
        j1 = i[np.arange(len(k)), np.argmin(c[i + sw] - c[i], 1)]
        i2 = j1[:, None] + np.arange(-r, r + 1)[None, :]
        m = (c[i2] - c[i2 - a]) / a - (c[i2 + b] - c[i2]) / b
        j = np.clip(np.argmax(m, 1), 1, m.shape[1] - 2)
        n = np.arange(len(j))
        y0, y1, y2 = m[n, j - 1], m[n, j], m[n, j + 1]
        den = y0 - 2 * y1 + y2
        frac = np.where(den < 0, 0.5 * (y0 - y2) / np.where(den < 0, den, -1.0), 0.0)
        te = (i2[n, j] + frac) / ofs
        keep = np.ones(len(te), bool)
        for _ in range(4):                            # VBI half-line pulses, click-hit lines: outliers
            p = np.polyfit(k[keep], te[keep], 1)
            res = te - np.polyval(p, k)
            s = 1.4826 * float(np.median(np.abs(res[keep])))
            keep = np.abs(res) < max(3 * s, 0.05e-6)
            if keep.sum() < 20:
                return None
        kk = k[keep]
        se = float(np.std(res[keep]) / np.sqrt(max(np.sum((kk - kk.mean()) ** 2), 1.0)))
        return float(p[1]), float(p[0]), se, s * 1e6

    def _prior_blend(self, t0, T, se, dur):
        """Weigh the block's line period against the detector's (rate0, +-PRIOR_PPM): a weak block
        times its edges loosely. A fit that disagrees by more than 4 combined sigmas is kept as is
        (the VTX's clock moved since detection). The sync phase is held at mid-block."""
        T0 = 1 / self.rate0
        sp = PRIOR_PPM * 1e-6 * T0
        if abs(T - T0) > 4 * np.hypot(se, sp):
            return t0, T
        Tn = (T / se ** 2 + T0 / sp ** 2) / (1 / se ** 2 + 1 / sp ** 2) if se > 0 else T
        km = (dur / 2 - t0) / T
        return t0 + km * (T - Tn), Tn

    def _sync_noise(self, d, ofs, t0, T, n_lines):
        """Robust std (Hz) of the sync tip (1.0-3.8 us after the edge), each line's own level removed."""
        u = np.linspace(1.0e-6, 3.8e-6, 24)
        s = np.interp((t0 + np.arange(n_lines)[:, None] * T + u[None, :]) * ofs, np.arange(len(d)), d)
        return 1.4826 * float(np.median(np.abs(s - np.median(s, 1, keepdims=True))))

    def chroma(self, d: np.ndarray, ofs: float):
        """Colour picture of the field raster() just selected, from the discriminator output
        low-passed to 6 MHz -> (uint8 [lines, WIDTH, 3] or None, meta)."""
        tm = self._timing
        if tm is None:
            return None, {}
        C = CHROMA[self.std]
        T, rows = tm["T"], tm["rows"]
        depth = tm["black"] - tm["sync"]                                  # luma "black" = blanking (porch)
        # a picture brighter than SYNC_COMPRESSED x the white the sync depth allows (chroma peaks reach
        # ~1.3x) means a compressed sync (weak VTX near the FM threshold, or the VTX itself):
        # take the depth from the picture's white instead
        wspan = C["white_per_sync"] - C["black_per_sync"]
        compressed = depth > 0 and tm["white"] - tm["black"] > SYNC_COMPRESSED * wspan * depth
        if compressed:
            depth = (tm["white"] - tm["black"]) / wspan
        black = tm["black"] + C["black_per_sync"] * depth
        span = (C["white_per_sync"] - C["black_per_sync"]) * depth
        lines = np.r_[rows[0] - 1, rows] if rows[0] > 0 else np.r_[rows[0], rows]   # + the line before (PAL)
        tau = np.arange(LINE_GRID) * (T / LINE_GRID)
        tl = tm["t0"] + lines[:, None] * T + tau[None, :]
        ok = tl[-1, -1] * ofs < len(d) - 1
        if not ok or depth <= 0:
            return None, dict(color=False)
        if tm["noise"] > NOISY_COLOR:    # chroma is noise there: grey beats colour (colour killer)
            return None, dict(color=False, color_killer="noise")
        R = np.interp(tl * ofs, np.arange(len(d)), d)                     # Hz
        lo = np.exp(-2j * np.pi * (C["k"] / T) * tl)
        h = firwin(95, 1.0e6, fs=LINE_GRID / T)
        Cx = oaconvolve(2 * R * lo, h[None, :], mode="same", axes=1)        # complex chroma, Hz
        tau_us = tau * 1e6
        # burst phasor from the gate alone (the 1 MHz chroma filter is ~6 us long: it would pull in
        # the sync edge and the picture); reference: the same gate 1 MHz below (audio is above).
        # Cameras place the burst up to ~0.4 us later than the standard (relative to the 50 % sync
        # edge): the gate is tried at BURST_SHIFT_US and the strongest field burst kept, as a TV's
        # generous burst gate
        best = None
        for sh in BURST_SHIFT_US:
            bg = (tau_us >= C["burst_us"][0] + sh) & (tau_us < C["burst_us"][1] + sh)
            Rg = R[:, bg] - R[:, bg].mean(1, keepdims=True)
            zb = (2 * Rg * lo[:, bg]).mean(1)
            lo_off = np.exp(-2j * np.pi * (C["k"] / T - 1.0e6) * tl[:, bg])
            zo = (2 * Rg * lo_off).mean(1)
            # field burst: |mean zb[l] conj zb[l-2]| is the burst power with the noise averaged out
            pb = np.abs(np.mean(zb[2:] * np.conj(zb[:-2])))
            po = np.abs(np.mean(zo[2:] * np.conj(zo[:-2])))
            # floored at its noise expectation (|mean of n noise products| ~ power / sqrt(n)): a
            # reference small by chance would make any gate look like a burst
            po = max(po, float(np.mean(np.abs(zo) ** 2)) / np.sqrt(max(len(zo) - 2, 1)))
            field = float(np.sqrt(pb / max(po, 1e-18)))
            if best is None or field > best[0]:
                best = (field, sh, zb, zo)
        field, sh, zb, zo = best
        ratio = np.abs(zb) / np.maximum(np.abs(zo), 1e-9)
        meta = dict(color=False, burst_ratio=round(float(np.median(ratio[1:])), 2), burst_field=round(field, 2),
                    burst_gate_us=sh, sync_compressed=bool(compressed))
        if field < BURST_MIN:
            return None, meta
        # phase reference: the burst averaged over same-parity neighbours
        zs = np.empty_like(zb)
        for par in (0, 1):
            z = zb[par::2]
            c = np.r_[0, np.cumsum(z)]
            i = np.arange(len(z))
            a, b = np.maximum(i - BURST_SMOOTH, 0), np.minimum(i + BURST_SMOOTH + 1, len(z))
            zs[par::2] = (c[b] - c[a]) / (b - a)
        # burst sin(ph + theta) mixes to e^{j(theta - pi/2)}: NTSC theta = pi, PAL pi +- pi/4
        zn = zs / np.maximum(np.abs(zs), 1e-9)
        if self.std == "pal":
            ref = zn[1:] + zn[:-1]
            sw = np.where(np.angle(zn[1:] * np.conj(ref)) >= 0, 1.0, -1.0)
        else:
            ref = zn[1:]
            sw = np.ones(len(rows))
        rot = np.exp(-1j * (np.angle(ref) - np.pi / 2))
        # ACC: the burst is burst_per_sync x the sync depth; chroma to luma units (white - black = 1)
        acc = C["burst_per_sync"] * depth / float(np.median(np.abs(zs[1:]))) / span
        Cn = Cx[1:] * (rot * acc)[:, None]
        U = -np.imag(Cn)
        V = np.real(Cn) * sw[:, None]
        if self.std == "pal":                                             # delay line
            Cp = Cx[:-1] * (np.r_[rot[:1], rot[:-1]] * acc)[:, None]
            swp = np.r_[-sw[:1], sw[:-1]]
            U = 0.5 * (U - np.imag(Cp))
            V = 0.5 * (V + np.real(Cp) * swp[:, None])
        Y = (R[1:] - np.real(Cx[1:] * np.conj(lo[1:])) - black) / span
        a_us, len_us = self.p["active_us"]
        u = a_us + np.arange(WIDTH) * (len_us / WIDTH)
        j = np.interp(u, tau_us, np.arange(LINE_GRID))                    # active part onto WIDTH columns
        j0 = np.minimum(j.astype(np.int64), LINE_GRID - 2); fr = j - j0
        def grid(A):
            return A[:, j0] * (1 - fr) + A[:, j0 + 1] * fr
        Y, U, V = grid(Y), grid(U), grid(V)
        if tm["noise"] > NOISY:          # below the FM threshold: clicks in Y, noise in the chroma band
            Y = median_filter(Y, 3, mode="nearest")
            U = median_filter(U, (3, 5), mode="nearest")
            V = median_filter(V, (3, 5), mode="nearest")
        rgb = np.stack([Y + 1.140 * V, Y - 0.395 * U - 0.581 * V, Y + 2.032 * U], -1)
        meta.update(color=True)
        return (np.clip(rgb, 0, 1) * 255).astype(np.uint8), meta

    def _field_rows(self, vbi: np.ndarray, n_lines: int) -> np.ndarray:
        """Rows of one field: the active lines after a vertical-sync run (or,
        if the block ends too soon after it, the lines before the next one)."""
        need, skip = self.p["lines"], self.p["skip"]
        runs = []
        i = 0
        while i < n_lines:
            if vbi[i]:
                j = i
                while j + 1 < n_lines and vbi[j + 1]:
                    j += 1
                if j - i >= 1:                         # two or more broad-pulse lines
                    runs.append((i, j))
                i = j + 1
            else:
                i += 1
        for s, e in runs:
            if e + skip + need <= n_lines:
                return np.arange(e + skip, e + skip + need)
        for s, e in runs:
            if s - need - 2 >= 0:
                return np.arange(s - need - 2, s - 2)
        return np.arange(min(need, n_lines))         # no VBI: rolling picture

    @staticmethod
    def _fold(d, fs, rate, nbf):
        ph = (np.arange(len(d)) / fs * rate) % 1.0
        b = (ph * nbf).astype(np.int64)
        return np.bincount(b, d, nbf) / np.maximum(np.bincount(b, minlength=nbf), 1)
