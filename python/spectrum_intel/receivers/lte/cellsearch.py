"""LTE downlink cell search: PSS fold over 5 ms, ZC ambiguity NMS, SSS (FDD/TDD, CP) -> PCI,
receiver ppm search, occupied bandwidth from energy."""
from __future__ import annotations

import numpy as np
import scipy.fft as sfft


LTE_BANDS = {"b28_dl": "B28", "b5_dl": "B5", "b8_dl": "B8", "b3_dl": "B3", "b1_dl": "B1", "b40": "B40",
             "b41": "B41"}
OFS = 3.84e6
N = 256                                  # 15 kHz useful symbol at 3.84 MSps
CP_N, CP0_N, CP_E = 18, 20, 64           # normal CP (other / first symbol of slot), extended
ROOTS = (25, 29, 34)
HALF = int(5e-3 * OFS)
Z_MIN = 15.0
BWS = (1.4e6, 3e6, 5e6, 10e6, 15e6, 20e6)
SSS_SMOOTH = 15                          # channel-estimate smoothing, subcarriers (225 kHz) per side of DC
OCC = {1.4e6: 1.08e6, 3e6: 2.7e6, 5e6: 4.5e6, 10e6: 9e6, 15e6: 13.5e6, 20e6: 18e6}


def _phy():
    from spectrum_intel.core import seq3gpp as V
    return V


_TPL = None
_TPL_F: dict[int, np.ndarray] = {}


def _templates_f(nf):
    """conj(FFT_nf(templates)), cached per transform length."""
    T = _TPL_F.get(nf)
    if T is None:
        T = _TPL_F[nf] = np.conj(sfft.fft(_templates(), nf, axis=1))
    return T


def _templates():
    global _TPL
    if _TPL is None:
        V = _phy()
        sc = V._SC62
        t = np.arange(N) / OFS
        out = []
        for u in ROOTS:
            x = (V.lte_pss_freq(u)[None, :] * np.exp(2j * np.pi * sc[None, :] * 15e3 * t[:, None])).sum(1)
            out.append(x / np.linalg.norm(x))
        _TPL = np.array(out)
    return _TPL


class _Slice:
    """One 3.84 MSps slice; frequency offsets are integer-bin rotations of its FFT."""

    def __init__(self, ch, f_center):
        self.f_center = f_center
        self.y, self.ofs = ch.slice(f_center, 1.0e6 + 60e3, OFS)
        L = len(self.y)
        self.L = L
        self.nf = sfft.next_fast_len(L + N)          # >= L + N - 1: circular = linear correlation
        self.Y = sfft.fft(self.y.astype(np.complex128), self.nf)
        self.T = _templates_f(self.nf)
        c = np.concatenate(([0.0], np.cumsum(np.abs(self.y.astype(np.complex128)) ** 2)))
        self.e = c[N:L] - c[:L - N] + 1e-30          # energy of y[n:n+N], n < L - N
        self.k = (L - N) // HALF

    def score(self, df, r):
        sh = int(round(df / self.ofs * self.nf))
        c = sfft.ifft(np.roll(self.Y, -sh) * self.T[r])[:self.L - N]
        m = np.abs(c) ** 2 / self.e
        fold = m[:self.k * HALF].reshape(self.k, HALF).sum(0)
        k = int(np.argmax(fold))
        med = np.median(fold)
        mad = 1.4826 * np.median(np.abs(fold - med)) + 1e-12
        far = np.abs((np.arange(len(fold)) - k + HALF // 2) % HALF - HALF // 2) > N
        self.iso = float((fold[k] - med) / (fold[far].max() - med + 1e-12))
        return float((fold[k] - med) / mad), k, m, sh * self.ofs / self.nf

    def best(self, offsets):
        b = None
        for df in offsets:
            for r in range(3):
                z, k, m, dfq = self.score(df, r)
                if b is None or z > b["z"]:
                    b = dict(z=z, df=dfq, nid2=r, phase=k, m=m, iso=self.iso)
        return b

    def shifted(self, df):
        return self.y * np.exp(-2j * np.pi * df * np.arange(self.L) / self.ofs)


def _instances(m, phase):
    """(half-frame index, sample) of every PSS actually present."""
    pk = np.max(m[phase::HALF])
    return [(j, p) for j, p in enumerate(range(phase, len(m), HALF)) if m[max(0, p - 2):p + 3].max() > 0.3 * pk]


def _smooth_half(h, n=SSS_SMOOTH):
    """Moving average over n subcarriers, normalised at the edges (no zero-padding bias)."""
    k = np.ones(n)
    return np.convolve(h, k, "same") / np.convolve(np.ones(len(h)), k, "same")


def sss_decode(y, nid2, inst):
    """Try FDD/TDD x normal/extended CP placements; return the most decisive.

    Per half-frame the SSS is combined with the PSS channel estimate by maximum-ratio combining, the estimate
    smoothed over SSS_SMOOTH subcarriers on each side of DC: a per-subcarrier zero-forcing equaliser amplifies
    noise wherever one PSS subcarrier is weak, and the channel is nearly flat over 225 kHz."""
    V = _phy()
    sc = V._SC62 % N
    d = V.lte_pss_freq(ROOTS[nid2])
    bank = V._SSS_BANK[nid2]                                    # [sub, nid1, 62]
    gaps = {("FDD", "normal"): N + CP_N, ("TDD", "normal"): 3 * (N + CP_N) + (CP0_N - CP_N),
            ("FDD", "extended"): N + CP_E, ("TDD", "extended"): 3 * (N + CP_E)}
    res = []
    for (dup, cp), gap in gaps.items():
        acc = np.zeros((2, 168))
        used = 0
        for k, p in inst:                 # k = half-frame index: SF0 / SF5 alternate with it
            s = p - gap
            if s < 0 or p + N > len(y):
                continue
            Yp = np.fft.fft(y[p:p + N])[sc]
            Ys = np.fft.fft(y[s:s + N])[sc]
            H = Yp / d
            H = np.concatenate([_smooth_half(H[:31]), _smooth_half(H[31:])])
            E = Ys * np.conj(H)
            E = E / (np.sqrt(np.mean(np.abs(E) ** 2)) + 1e-30)
            v = np.real(np.einsum("snk,k->sn", bank, E))
            for par in (0, 1):
                acc[par] += v[0 if (k + par) % 2 == 0 else 1]
            used += 1
        if not used:
            continue
        flat = np.sort(acc.ravel())
        par, nid1 = np.unravel_index(np.argmax(acc), acc.shape)
        med = np.median(flat)
        noise = 1.4826 * np.median(np.abs(flat - med)) + 1e-12
        res.append(dict(duplex=dup, cp=cp, nid1=int(nid1), parity=int(par), n=used,
                        z=float((flat[-1] - med) / noise), z2=float((flat[-2] - med) / noise)))
    return max(res, key=lambda r: r["z"] / max(r["z2"], 1e-3)) if res else None


def sss_decisive(s) -> bool:
    """Many of the 336 hypotheses share one of the two interleaved m-sequences
    with the true one, so the runner-up sits near half the winner: decide on
    the ratio, not on a z gap."""
    return s is not None and s["z"] >= 4.5 and s["z"] / max(s["z2"], 1e-3) >= 1.35


def wide_ppm_search(ch, centres, span=500e3, step=5e3):
    """Strongest PSS over the given centres with a wide offset search -> ppm."""
    best = None
    for c in centres:
        sl = _Slice(ch, c)
        b = sl.best(np.arange(-span, span + step / 2, step))
        if best is None or b["z"] > best["z"]:
            best = dict(b, centre=c)
    if best is None or best["z"] < 30:
        return None
    # refine: the 5 kHz grid alone leaves +-2.5 kHz (+-2.7 ppm at 900 MHz)
    sl = _Slice(ch, best["centre"])
    fine = sl.best(np.arange(best["df"] - 3e3, best["df"] + 3.01e3, 250.0))
    if fine["z"] >= best["z"]:
        best.update(z=fine["z"], df=fine["df"])
    f = best["centre"] + best["df"]
    raster = round(f / 100e3) * 100e3
    return {"ppm": (f - raster) / raster * 1e6, "z": best["z"], "carrier_hz": f, "raster_hz": raster}


FLAT_DB = 8.0          # an LTE carrier's long-term (reference-signal) profile is flat; adjacent carriers and gaps are not


def bandwidth_from_energy(det, fcar):
    """Largest 3GPP BW whose occupied band around fcar is always-on energy with a flat profile (10-90th percentile
    spread <= FLAT_DB): adjacent always-on carriers (GSM beside a 5 MHz LTE) do not widen the claim."""
    C, fl = det.coarse, det.floor_coarse_db
    prof = 10 * np.log10(np.percentile(C.P, 20, axis=0) + 1e-30) + 6.5 - fl
    best = None
    for bw in BWS:
        h = OCC[bw] / 2
        sel = (C.f > fcar - h) & (C.f < fcar + h)
        if sel.sum() < 3 or C.f[sel][0] - C.f[0] < 0 or C.f[-1] - C.f[sel][-1] < 0:
            break
        if np.mean(prof[sel] > 3) >= 0.85 and np.percentile(prof[sel], 90) - np.percentile(prof[sel], 10) <= FLAT_DB:
            best = bw
        else:
            break
    return best


def raster_offsets(f_hz, ppm_known):
    """PSS frequency hypotheses around a ppm-corrected raster carrier. With the receiver error known, the
    residual is within ~2.5 ppm (drift); below 3 kHz (f < 1.2 GHz) the 0 Hz hypothesis alone loses < 0.3 dB."""
    return (0.0,) if ppm_known and 2.5e-6 * f_hz < 3e3 else (-5e3, 0.0, 5e3)


def raster_grid(sig):
    """EARFCN raster carriers whose PSS (+-0.55 MHz) lies inside the energy signal."""
    lo, hi = sig.f0 + 0.55e6, sig.f1 - 0.55e6
    if hi < lo:
        lo = hi = round(sig.fc / 100e3) * 100e3
    return np.arange(np.ceil(lo / 100e3) * 100e3, hi + 1, 100e3)


def search_cells(ctx, sig, ppm, done=None, ppm_known=False, require_sss=False):
    """All PSS/SSS cells on the raster inside an energy signal. done: rasters already searched for another
    (overlapping) energy signal of this capture; skipped here and updated. require_sss: report only cells
    whose SSS is decisive (no PSS-only cells)."""
    rasters = raster_grid(sig)
    hits = []
    for r in rasters:
        if done is not None:
            key = int(round(r / 100e3))
            if key in done:
                continue
            done.add(key)
        guess = r * (1 + ppm * 1e-6)
        sl = _Slice(ctx.ch, guess)
        b = sl.best(raster_offsets(r, ppm_known))
        if b["z"] >= Z_MIN:
            hits.append((r, sl, b))
    # non-max suppression: one cell per (carrier, N_ID2) within +-300 kHz
    hits.sort(key=lambda h: -h[2]["z"])
    keep = []
    for h in hits:
        if any(abs(h[0] - k[0]) < 300e3 and h[2]["nid2"] == k[2]["nid2"] for k in keep):
            continue
        keep.append(h)
    cells = []
    for r, sl, b in keep:
        fine = sl.best(np.arange(b["df"] - 3e3, b["df"] + 3.01e3, 500.0))
        if fine["nid2"] != b["nid2"]:
            fine = b
        y = sl.shifted(fine["df"])
        inst = _instances(fine["m"], fine["phase"])
        s = sss_decode(y, fine["nid2"], inst) if len(inst) >= 2 else None
        fcar = sl.f_center + fine["df"]
        sss_ok = sss_decisive(s)
        # a loaded OFDM band leaves heavy-tailed fold peaks: without a decisive SSS
        # only a very strong, isolated PSS counts
        if not sss_ok and (require_sss or not (fine["z"] >= 40 and fine.get("iso", 0) >= 3)):
            continue
        cells.append(dict(raster_hz=r, carrier_hz=fcar, pss_z=fine["z"], pss_isolation=fine.get("iso"),
                          nid2=fine["nid2"], inst=inst, sss=s, sss_ok=sss_ok, err_hz=fcar - r,
                          phase=fine["phase"]))
    # Zadoff-Chu time-frequency ambiguity: a PSS offset by an integer number of
    # subcarriers (e.g. 300 kHz = 20) still correlates at ~the same timing. Drop a
    # PSS-only hit when a stronger / SSS-resolved cell with the same N_ID2 and
    # timing (+-1 symbol) lies within 20 MHz.
    cells.sort(key=lambda c: (not c["sss_ok"], -c["pss_z"]))
    out = []
    for c in cells:
        if not c["sss_ok"] and any(o["nid2"] == c["nid2"] and abs(o["carrier_hz"] - c["carrier_hz"]) < 20e6 and
                                   abs((o["phase"] - c["phase"] + HALF // 2) % HALF - HALF // 2) <= N
                                   for o in out):
            continue
        out.append(c)
    return out
