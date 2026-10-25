#!/usr/bin/env python3
"""
det_digital_vtx.py -- non-DJI digital video / telemetry links.

    detect(x, fs, fc, ctx) -> list[Emission]

Two detectors share one FFT of the capture:

1. family "herelink"  (CubePilot Herelink, Pinecone-S1 SoC, "LTE over Wi-Fi")
   Measured on testdata/Herelink_20251125-110453 (see HERELINK below): the
   air unit transmits an LTE-Release-8 *downlink* PHY almost verbatim --
   15 kHz SCS, 2048-pt FFT @30.72 MSps, normal CP (160/144), 100 RB = 18 MHz,
   LTE PSS (Zadoff-Chu root 25/29/34) + SSS in the FDD positions (symbols 5/6
   of slot 0 of subframes 0 and 5), cell-specific reference signals -- used in
   a 5 ms TDD rhythm: subframes 0,3,4,5,8,9 downlink (TDD UL/DL config 2 shape),
   1,2,6,7 silent at the air unit.
   Detection = PSS matched filter (66.7 us coherent, the full symbol energy of a
   power-boosted symbol) folded non-coherently over the 5 ms half-frame, searched
   blind over every 10 kHz carrier hypothesis in the band. Per half-frame that is
   ~31 dB of processing gain over the 18 MHz in-band SNR, so the link is found
   ~10 dB below where the CP / energy detectors can still see individual bursts.
   SSS then gives the physical cell id (a stable per-pairing identifier), the
   10 ms frame phase, and the slot grid that per-subframe classification uses.

2. family "digital_vtx" (any CP-OFDM video / data link, blind.py): a coherent
   lag-Tu autocorrelation scan of the whole band (time-AGC'd, so Wi-Fi bursts do
   not dominate) -> symbol-rate line / per-window fold confirmation -> Tu, CP,
   symbol period -> CP cross-spectrum band -> bursts on local symbol grids (a
   link that re-times mid-capture keeps its bursts) -> TDD frame completion.
   No vendor knowledge: a DJI O4 link 5 dB below the noise, Herelink and an
   unknown drone are found the same way. Label "ofdm_link_<scs>kHz"; the
   numerology table (KNOWN_NUMEROLOGIES) says what it resembles and whether that
   is a drone family, infrastructure or unknown. Tu < 10 us (802.11 family) is
   left to the Wi-Fi receivers.

ctx keys used (all optional): "channelizer" (Channelizer of this x -- its FFT is
reused), "spurs" (spectrum_survey.find_spurs output), "generic" (bool, default
True), "report_wifi" (default True), "herelink" (default True).
"""
from __future__ import annotations
import numpy as np
from scipy import ndimage, fft as sfft

from spectrum_intel.core import common as C
from spectrum_intel.core.emission import Emission
from spectrum_intel.core.parallel import pmap
from spectrum_intel.core.spectral import find_spurs
from spectrum_intel.core import bursts as BU
from spectrum_intel.core import bands_india as BI
from spectrum_intel.core.seq3gpp import _SC62, _SSS_BANK, lte_pss_freq
from . import blind as BL
from . import streams as ST
from . import catalog as CAT
from . import hdzero as HZ
# ============================================================== numerology table
# tu_us = useful symbol, cp_us = (normal) cyclic prefix. "verified" means measured on a capture
# in testdata/ (or tmpdata/) by this library; everything else is public info only. role: what a
# blind link of this numerology in an ISM band most likely is (drone | infrastructure | wifi).
KNOWN_NUMEROLOGIES = [
    dict(name="wifi_he", tu_us=12.8, cp_us=3.2, scs_khz=78.125, bw_mhz=(4.4, 19.1, 38.3, 78.1), role="wifi",
         status="verified as 802.11a at a 5 MHz clock (D_20251119-183758, DJI ad-hoc link: wifi_ofdm_scaled)",
         source="IEEE 802.11ax data symbols (GI 0.8/1.6/3.2) / 802.11a quarter clock (GI 3.2)"),
    dict(name="herelink", tu_us=66.67, cp_us=4.69, scs_khz=15.0, bw_mhz=(18.0, 9.0), role="drone",
         status="verified 20 MHz mode (Herelink_20251125, blind: 18.0 MHz, 5 ms TDD); 10 MHz mode unverified",
         source="this module's measurement; CubePilot spec sheet lists 20/10 MHz BW"),
    dict(name="lte_like", tu_us=66.67, cp_us=4.69, scs_khz=15.0, bw_mhz=(1.08, 4.5, 9.0, 18.0), role="drone",
         status="verified: DJI OcuSync 2 video 10/20 MHz (Mini 2, test_signal 2484 MHz), Herelink",
         source="3GPP 36.211 normal CP; in an ISM band: DJI OcuSync / O3, Herelink, Autel SkyLink"),
    dict(name="dji_like", tu_us=66.67, cp_us=4.69, scs_khz=15.0, bw_mhz=(1.4, 10.0, 20.0, 40.0), role="drone",
         status="15 kHz numerology seen on DJI captures; the dji receiver owns classification",
         source="OcuSync / O3 public SDR analyses (LTE-like)"),
    dict(name="autel_skylink_like", tu_us=66.67, cp_us=4.69, scs_khz=15.0, bw_mhz=(9.0, 1.1), role="drone",
         status="unverified (no sample)",
         source="vendor C-UAS white paper: 600 sc / ~9 MHz downlink, 72 sc / 1.1 MHz uplink, 15 kHz SCS"),
    dict(name="dji_o4_like", tu_us=33.33, cp_us=2.34, scs_khz=30.0, bw_mhz=(9.0, 18.0, 36.0), role="drone",
         status="verified: DJI O4 video 10 MHz (Mini 4 Pro) and 20 MHz (D_20260609-151042), 2.5 ms TDD",
         source="this library's measurement (NR-like 30 kHz, CP 72/1024, comb-4 repeated reference)"),
    dict(name="lightbridge_like", tu_us=91.43, cp_us=11.43, scs_khz=10.94, bw_mhz=(9.2,), role="drone",
         status="seen once: DJI_Phantom_4 testdata, 9.6 MHz at 2426.6 MHz (true), CP 1/8, bursty; identity "
                "not verified (Lightbridge 2 PHY is not public)",
         source="802.16e (WiMAX) OFDMA 10 MHz profile: 1024 FFT at 11.2 MSps; Phantom 4 uses Lightbridge 2"),
    dict(name="walksnail_like", tu_us=91.43, cp_us=11.43, scs_khz=10.94, bw_mhz=(18.9,), role="drone",
         status="unverified (no sample)",
         source="waveshark issue #226: 802.16e OFDMA 20 MHz profile, 22.4 MSps, 2048 FFT, CP 1/8, 1729 used sc"),
    dict(name="pmp450_like", tu_us=22.86, cp_us=1.43, scs_khz=43.75, bw_mhz=(18.5,), role="infrastructure",
         status="seen: tmpdata D_20260515-192148 / -192835, 5765.0 MHz, 2.5 ms TDD, Tu/3 preamble",
         source="Cambium PMP/PTP 450 (512 FFT at 22.4 MSps, CP 1/16): fixed-wireless ISP link, not a drone"),
    dict(name="hdzero_like", tu_us=37.926, cp_us=1.185, scs_khz=26.367, bw_mhz=(23.44,), role="drone",
         status="verified 27 MHz mode: HDZero Freestyle V2 (LIVE_HDZERO_5658MHz; the hdzero path names it by its "
                "pilots); 17 MHz narrow mode not measured",
         source="this library's measurement: 1024 FFT at 27 MSps, CP 32, 889 used subcarriers, 26 continual pilots"),
]


def classify_numerology(tu_us, cp_us=None, bw_mhz=None):
    """Table entries within 1.5 % (+0.1 us) of Tu. Several share 15 kHz (all returned so other
    plugins can decide; Herelink needs its PSS); entries sharing a Tu are told apart by the
    nearest listed bandwidth when bw_mhz is given. -> (name, hits)."""
    hits = [n["name"] for n in KNOWN_NUMEROLOGIES
            if n["tu_us"] and abs(tu_us - n["tu_us"]) < 0.015 * n["tu_us"] + 0.1]
    if not hits:
        return "unknown", []
    if "lte_like" in hits:
        return "lte_like", [h for h in hits if h != "herelink"]   # herelink needs its PSS
    if bw_mhz is not None and len(hits) > 1:
        bws = {n["name"]: n["bw_mhz"] for n in KNOWN_NUMEROLOGIES}
        return min(hits, key=lambda h: (min(abs(bw_mhz - b) for b in bws[h]), hits.index(h))), hits
    return hits[0], hits


def _pss_time(u, ofs):
    """Useful part (66.67 us, no CP) of the PSS symbol sampled at ofs."""
    n = int(round(ofs / 15e3))
    t = np.arange(n) / ofs
    return (lte_pss_freq(u)[None, :] * np.exp(2j * np.pi * _SC62[None, :] * 15e3 * t[:, None])).sum(1)


# Herelink, as measured (all times relative to subframe start, 30.72 MSps units)
HERELINK = dict(
    scs_hz=15e3, tu_us=2048 / 30.72, cp_us=144 / 30.72, cp0_us=160 / 30.72, symbol_us=2192 / 30.72,
    half_frame_ms=5.0, frame_ms=10.0, subframe_ms=1.0,
    pss_symbol=6, sss_symbol=5,                    # slot 0 of subframes 0 and 5 (FDD placement)
    pss_offset_us=(160 + 2048 + 5 * 2192 + 144) / 30.72,       # subframe start -> PSS useful part
    dl_subframes=(0, 3, 4, 5, 8, 9), idle_subframes=(1, 2, 6, 7),
    bw_mhz=18.0, n_rb=100,
)
_ROOTS = {25: 0, 29: 1, 34: 2}


# ============================================================== helpers
def _fft_of(x, ctx):
    ch = ctx.get("channelizer")
    if ch is not None and getattr(ch, "L", None) == len(x):
        return ch.X
    return np.fft.fft(x.astype(np.complex128))


def _slice(X, fs, f_off, half_bw, out_fs):
    """Channelizer.slice on a bare FFT (brick-wall), complex128."""
    L = len(X)
    M = int(round(L * out_fs / fs))
    k0 = int(round(f_off / fs * L))
    kk = int(half_bw / fs * L)
    m = np.arange(-kk, kk + 1)
    Y = np.zeros(M, complex)
    Y[m % M] = X[(k0 + m) % L]
    return np.fft.ifft(Y) * (M / L), M * fs / L


class _SliceCache:
    """The blind search asks for the same slice (centre, half band, rate) several times (scan -> confirm -> segment,
    links -> hops): keep the last few, read-only (PLAN P1; identical results)."""

    def __init__(self, X, fs, size=6):
        self.X, self.fs, self.size, self.d = X, fs, size, {}

    def __call__(self, f, h, r):
        k = (float(f), float(h), float(r))
        v = self.d.pop(k, None)
        if v is None:
            y, o = _slice(self.X, self.fs, f, h, r)
            y.flags.writeable = False
            v = (y, o)
            if len(self.d) >= self.size:
                self.d.pop(next(iter(self.d)))
        self.d[k] = v
        return v


def _notch(X, fs, spurs, half=30e3):
    if not spurs:
        return X
    X = X.copy()
    L = len(X)
    f = np.fft.fftfreq(L, 1 / fs)
    for s in spurs:
        X[np.abs(f - s["offset_hz"]) <= half] = 0
    return X


# ============================================================== Herelink: PSS fold search
def pss_fold_search(X, fs, usable_half, step_hz=10e3, period_s=5e-3, batch=48, roots=(25, 29, 34), clip=24.0,
                    f_range=None):
    """Blind PSS search over carrier hypotheses f in [-usable_half, usable_half].
    For every hypothesis: energy-normalised matched-filter output q(n) (~Exp(1)
    in noise, also inside strong non-PSS bursts), clipped, folded modulo the 5 ms
    half-frame (mean over frames), +-1 sample dilated (clock drift). The fold
    peak is normalised by max(row median, global median): a row whose median
    collapses (a CW or empty notch inside the 1.9 MHz slice correlates with
    nothing) must not inflate its own peak. Returns Z[root][f], phase T[root][f]."""
    L = len(X)
    M = int(round(L * 1.92e6 / fs)); ofs = M * fs / L
    P = period_s * ofs
    Pint = int(np.ceil(P))
    integral = abs(P - round(P)) < 1e-9
    ph = np.floor(np.mod(np.arange(M), P)).astype(np.int64)
    cnt = np.bincount(ph, minlength=Pint).astype(float)
    cnt[cnt == 0] = 1
    tpl = {u: _pss_time(u, ofs) for u in roots}
    nt = len(tpl[roots[0]])
    # spectra kept in fftshift order: the IFFT of a circularly shifted spectrum is
    # the true signal times a unit-modulus ramp, and only |.|^2 is used below
    Pc = {u: np.fft.fftshift(np.conj(np.fft.fft(tpl[u], M))).astype(np.complex64) for u in roots}
    Ep = float(np.sum(np.abs(tpl[roots[0]]) ** 2))
    sb = max(1, int(round(step_hz / fs * L)))                     # integer-bin step -> strided gather
    kb = np.arange(int(-usable_half / fs * L), int(usable_half / fs * L) + 1, sb)
    if f_range is not None:                                       # only the hypotheses inside (lo, hi) Hz offsets
        kb = kb[(kb * fs / L >= f_range[0]) & (kb * fs / L <= f_range[1])]
    offs = kb * fs / L
    if not len(offs):
        return dict(offs=offs, Z={u: np.zeros(0) for u in roots}, T={u: np.zeros(0, int) for u in roots}, ofs=ofs,
                    M=M, P=P, nt=nt, global_median=0.0)
    Xs = np.fft.fftshift(X).astype(np.complex64)
    Xs = np.concatenate([Xs, Xs[:M]])
    win = np.lib.stride_tricks.sliding_window_view(Xs, M)
    lo = L // 2 - M // 2                                          # window start for offset 0
    pk = {u: np.zeros(len(offs)) for u in roots}
    md = {u: np.zeros(len(offs)) for u in roots}
    T = {u: np.zeros(len(offs), int) for u in roots}
    K = int(M // P)

    def one_batch(b0):                                            # batches are independent: threads
        kk = kb[b0:b0 + batch]
        B = len(kk)
        Yb = np.ascontiguousarray(win[(lo + kk) % L])
        e = np.abs(sfft.ifft(Yb, axis=1, workers=-1)) ** 2
        # sliding energy over the template length (circular), floored so empty
        # (notched / out-of-band) stretches cannot blow the ratio up
        cs = np.cumsum(np.concatenate([e, e[:, :nt]], 1), 1, dtype=np.float64)
        Ey = (cs[:, nt - 1:nt - 1 + M] - np.concatenate([np.zeros((B, 1)), cs[:, :M - 1]], 1)).astype(np.float32)
        Ey += 1e-3 * Ey.mean(1, keepdims=True)
        inv = (nt / Ep) / Ey
        for u in roots:
            r = sfft.ifft(Yb * Pc[u][None, :], axis=1, workers=-1)
            # clip: one strong chirp-like interferer instance (GFSK looks like a ZC)
            # must not carry the fold; a real PSS train contributes in every frame
            q = np.minimum((r.real ** 2 + r.imag ** 2) * inv, clip)
            if integral:                                          # fold by reshape (fast path)
                Pi = int(round(P))
                Q = q[:, :K * Pi].reshape(B, K, Pi).sum(1)
                Q[:, :M - K * Pi] += q[:, K * Pi:]
            else:
                idx = (np.arange(B)[:, None] * Pint + ph[None, :]).ravel()
                Q = np.bincount(idx, weights=q.ravel(), minlength=B * Pint).reshape(B, Pint)
            Q = Q / cnt[None, :Q.shape[1]]
            Q = np.maximum(Q, np.maximum(np.roll(Q, 1, 1), np.roll(Q, -1, 1)))
            i = np.argmax(Q, 1)
            pk[u][b0:b0 + B] = Q[np.arange(B), i]
            md[u][b0:b0 + B] = np.median(Q, 1)
            T[u][b0:b0 + B] = i

    pmap(one_batch, range(0, len(offs), batch))
    G = float(np.median(np.concatenate(list(md.values()))))
    Z = {u: pk[u] / np.maximum(md[u], G) for u in roots}
    return dict(offs=offs, Z=Z, T=T, ofs=ofs, M=M, P=P, nt=nt, global_median=G)


def _pss_candidates(sr, thr, max_cands=4, ridge=400e3):
    """Hypotheses above thr, strongest first, one per +-ridge (any root: at high
    SNR the ZC delay-Doppler ambiguity lights up neighbouring carriers and the
    other roots too; the refinement resolves the exact carrier)."""
    rows = []
    for u, z in sr["Z"].items():
        for i in np.where(z > thr)[0]:
            rows.append((float(z[i]), u, i))
    rows.sort(reverse=True)
    out = []
    for z, u, i in rows:
        f = float(sr["offs"][i])
        if any(abs(f - o["f_off"]) <= ridge for o in out):
            continue
        out.append(dict(z=z, root=u, f_off=f, tau=int(sr["T"][u][i])))
        if len(out) >= max_cands:
            break
    return out


# ============================================================== Herelink: refinement
def _corr_q(seg, tpl):
    """Normalised correlation q (Exp(1) in noise) of tpl against every lag of seg."""
    nt = len(tpl)
    r = np.correlate(seg, tpl, "valid")
    e = np.convolve(np.abs(seg) ** 2, np.ones(nt), "valid")
    return nt * np.abs(r) ** 2 / (np.sum(np.abs(tpl) ** 2) * e + 1e-30), r


def _carrier_search(X, fs, f0, t_inst_s, roots=(25, 29, 34), n_sc=28, win_s=45e-6):
    """Resolve the exact carrier and root around a fold hit. The fold hit can sit
    on an integer-subcarrier ambiguity of the ZC (strong response at a shifted
    lag) -- so test f0 + k*15 kHz + d for k in +-n_sc, d in a 2.5 kHz grid, and
    keep the hypothesis with the largest summed per-instance peak q."""
    y, ofs = _slice(X, fs, f0, 0.96e6 + n_sc * 15e3, 3.84e6)
    n = np.arange(len(y))
    best = (-1.0, None, None)
    tpls = {u: _pss_time(u, ofs) for u in roots}
    w = int(win_s * ofs)
    wins = [(int(round(t * ofs)) - w, int(round(t * ofs)) + w + len(tpls[roots[0]])) for t in t_inst_s]
    wins = [(a, b) for a, b in wins if a >= 0 and b <= len(y)]
    if not wins:
        return f0, roots[0], 0.0
    for d in np.arange(-7.5e3, 7.5e3, 2.5e3):
        for k in range(-n_sc, n_sc + 1):
            df = k * 15e3 + d
            for u in roots:
                tp = tpls[u] * np.exp(2j * np.pi * df * np.arange(len(tpls[u])) / ofs)
                s = 0.0
                for a, b in wins:
                    q, _ = _corr_q(y[a:b], tp)
                    s += q.max()
                if s > best[0]:
                    best = (s, f0 + df, u)
    # fine: +-1.25 kHz in 250 Hz steps at the winner
    s0, f1, u = best
    for d in np.arange(-1.25e3, 1.26e3, 250.0):
        tp = tpls[u] * np.exp(2j * np.pi * (f1 + d - f0) * np.arange(len(tpls[u])) / ofs)
        s = sum(_corr_q(y[a:b], tp)[0].max() for a, b in wins)
        if s > best[0]:
            best = (s, f1 + d, u)
    return best[1], best[2], best[0]


def _pss_purity(X, fs, f_c, u, t_s):
    """median q / expected q for a pure PSS at the measured window SNR (1.92 MSps,
    centre 1.92 MHz). Window SNR s = window energy over the slice's noise floor
    (20th percentile of the sliding energy over the whole capture) minus 1."""
    y, ofs = _slice(X, fs, f_c, 0.96e6, 1.92e6)
    tpl = _pss_time(u, ofs); nt = len(tpl)
    e = np.convolve(np.abs(y) ** 2, np.ones(nt), "valid")
    n0 = np.percentile(e, 20)
    qs, ss = [], []
    for t in t_s:
        a = int(round(t * ofs)) - 3
        if a < 0 or a + nt + 6 > len(y):
            continue
        q, _ = _corr_q(y[a:a + nt + 6], tpl)
        i = int(np.argmax(q))
        qs.append(q[i]); ss.append(max(e[a + i] / n0 - 1, 1e-3))
    if not qs:
        return 1.0, 0.0
    s_med = float(np.median(ss))
    return float(np.median(qs) / (nt * s_med / (1 + s_med))), s_med


def _instances(y, ofs, tpl, t_guess, win):
    """Peak of the normalised correlation of tpl around each guessed start
    (sample units of y), with parabolic sub-sample timing."""
    nt = len(tpl)
    out = []
    for g in t_guess:
        a = int(round(g)) - win
        if a < 0 or a + 2 * win + nt > len(y):
            out.append(None); continue
        q, r = _corr_q(y[a:a + 2 * win + nt], tpl)
        i = int(np.argmax(q))
        di = 0.0
        if 0 < i < len(q) - 1:
            l, c, rr = np.abs(r[i - 1:i + 2])
            den = l - 2 * c + rr
            di = 0.5 * (l - rr) / den if den != 0 else 0.0
        out.append(dict(t=a + i + di, q=float(q[i]), r=r[i]))
    return out


def _sss_decode(y, ofs, nid2, pss_useful_idx):
    """Equalise SSS (the symbol just before PSS: FDD placement) with the PSS
    channel and test all 168 x 2 hypotheses, accumulated over every half-frame
    (SSS alternates between the subframe-0 and subframe-5 sequences)."""
    N = int(round(ofs / 15e3)); sym = int(round(2192 / 2048 * N))
    d = lte_pss_freq(list(_ROOTS)[nid2])
    bank = _SSS_BANK[nid2]                                        # [sub, nid1, 62]
    acc = np.zeros((2, 168))                                      # [parity of instance 0, nid1]
    used = 0
    for k, p in enumerate(pss_useful_idx):
        p = int(round(p)); s = p - sym
        if s < 0 or p + N > len(y):
            continue
        Yp = np.fft.fft(y[p:p + N])[_SC62 % N]
        Ys = np.fft.fft(y[s:s + N])[_SC62 % N]
        H = Yp / d
        E = Ys * np.conj(H) / (np.abs(H) ** 2 + np.mean(np.abs(H) ** 2) * 0.05)
        sc = np.real(np.einsum("snk,k->sn", bank, E))             # [sub, nid1]
        for par in (0, 1):                                        # instance k is sf0 if (k+par) even
            acc[par] += sc[0 if (k + par) % 2 == 0 else 1]
        used += 1
    if used == 0:
        return None
    flat = np.sort(acc.ravel())
    par, nid1 = np.unravel_index(np.argmax(acc), acc.shape)
    med = np.median(flat)
    noise = 1.4826 * np.median(np.abs(flat - med)) + 1e-12
    return dict(nid1=int(nid1), pci=int(3 * nid1 + nid2), parity=int(par), z=float((flat[-1] - med) / noise),
                z_second=float((flat[-2] - med) / noise), n=used)


def _cfo_cp(y, ofs, sym_starts):
    """Fractional CFO (+-7.5 kHz) from the CP correlation phase of the given
    symbols (sample index of each CP start, normal-CP length used)."""
    N = int(round(ofs / 15e3)); cp = int(round(144 / 2048 * N))
    acc = 0j; per = []
    for s in sym_starts:
        s = int(round(s))
        if s < 0 or s + N + cp > len(y):
            continue
        v = np.vdot(y[s + N:s + N + cp], y[s:s + cp])             # sum y[n] y*[n+N]
        acc += v; per.append(v)
    if not per:
        return None, []
    return float(-np.angle(acc) / (2 * np.pi * N / ofs)), per


def _sym_starts(t_sf, ofs, n_slots=2):
    """CP start (samples) of the 14 symbols of the subframe starting at t_sf (s)."""
    u = ofs / 30.72e6
    a = t_sf * ofs
    return [a + u * (sl * 15360 + (0 if s == 0 else 2208 + (s - 1) * 2192)) for sl in range(n_slots) for s in range(7)]


HERELINK_RF = (2370e6, 2520e6)    # Herelink (Pinecone S1) is 2.4 GHz only; +-20 MHz for frequency-label errors


def herelink_detect(x, fs, fc, X, ctx, info):
    usable = min(fs, 56e6) / 2 - 0.6e6
    # latency (PLAN P1): the PSS fold search costs ~5000 hypotheses x 4 IFFTs per block; Herelink exists only at
    # 2.4 GHz, so only carrier hypotheses inside HERELINK_RF are searched (none at all elsewhere)
    lo, hi = max(-usable, HERELINK_RF[0] - fc), min(usable, HERELINK_RF[1] - fc)
    if lo > hi:
        info["pss_search"] = dict(skipped="capture outside the Herelink band (2.4 GHz)")
        return []
    sr = pss_fold_search(X, fs, usable, f_range=None if (lo <= -usable and hi >= usable) else (lo, hi))
    if not len(sr["offs"]):
        return []
    thr = ctx.get("pss_threshold", PSS_FOLD_THRESHOLD)
    info["pss_search"] = dict(best_z={u: round(float(np.max(z)), 2) for u, z in sr["Z"].items()},
                              threshold=thr, hypotheses=int(len(sr["offs"]) * len(sr["Z"])),
                              global_median=round(sr["global_median"], 3))
    out, accepted = [], []
    for c in _pss_candidates(sr, thr):
        if any(abs(c["f_off"] - a) < 9.5e6 for a, _ in accepted):
            info.setdefault("rejected", []).append(dict(c, why="inside an accepted Herelink band"))
            continue
        # receiver IQ imbalance mirrors a strong carrier to -f (conjugated: root 25 <-> 34)
        if any(abs(c["f_off"] + a) < 0.3e6 and c["z"] < 0.8 * z for a, z in accepted):
            info.setdefault("rejected", []).append(dict(c, why="IQ image of an accepted Herelink carrier"))
            continue
        em, sub = _herelink_refine(x, fs, fc, X, c, sr)
        info.setdefault("links", []).append(sub)
        if sub["accepted"]:
            accepted.append((sub["f_off"], c["z"]))
            out += em
    return out


def hdzero_detect(x, fs, fc, X, ctx, info):
    """HDZero 27 MHz-mode video links on the HDZero channels inside the capture (hdzero.py): continual-pilot
    coherence after the CP-phase CFO; one Emission per link over the whole capture (the link is continuous)."""
    log = info.setdefault("hdzero_log", []).append
    links = HZ.hdzero_links(lambda f, h, r: _slice(X, fs, f, h, r), fs, fc, min(fs, 56e6) / 2, len(x), log=log)
    if not links:
        return []
    noise = C.noise_floor_psd(x)["floor_power_per_sample"]
    out = []
    for lk in links:
        p_noise = noise * 2 * lk["half_bw"] / fs
        snr = float(10 * np.log10(max(lk["p_band"] / max(p_noise, 1e-30) - 1, 1e-3)))
        lo, hi = (fc + lk["lo"]) / 1e6, (fc + lk["hi"]) / 1e6
        cm = lk["channel_mhz"] + lk["cfo_hz"] / 1e6
        fr = lk["frame_ms"]
        feats = dict(channel=lk["channel"], channel_mhz=lk["channel_mhz"], carrier_mhz=round(cm, 4),
                     cfo_hz=round(lk["cfo_hz"], 1), cfo_ppm=round(lk["ppm"], 2), bw_mhz=round(hi - lo, 2),
                     bandwidth_mode="27 MHz (normal)", scs_khz=round(HZ.SCS / 1e3, 3), tu_us=round(HZ.N / HZ.FS * 1e6, 3),
                     cp_us=round(HZ.CP / HZ.FS * 1e6, 3), symbol_us=round(HZ.TS / HZ.FS * 1e6, 3),
                     fft_size=HZ.N, sample_rate_msps=HZ.FS / 1e6, used_subcarriers=2 * HZ.HALF + 1,
                     continual_pilots=len(HZ.CONTINUAL), pilot_coherence=round(lk["pilot_coherence"], 3),
                     data_coherence=round(lk["data_coherence"], 3), cp_fold_contrast=round(lk["fold_contrast"], 1),
                     n_symbols=lk["n_symbols"], preamble_ms=lk["preamble_ms"],
                     frame_symbols=lk["frame_symbols"], frame_period_ms=None if fr is None else round(fr, 3),
                     frame_rate_fps=None if fr is None else round(1e3 / fr, 2), link_snr_db=round(snr, 1))
        out.append(Emission("hdzero", "hdzero_video", round(lk["t0"] * 1e6, 1), round(lk["t1"] * 1e6, 1),
                            round(lo, 3), round(hi, 3), round((lo + hi) / 2, 3), round(snr, 1), 0.9, features=feats))
    return out


PSS_FOLD_THRESHOLD = 7.0         # set from the null captures (see test_det_digital_vtx.py)
PSS_PURITY_MIN = 0.5              # median q / pure-PSS q (see _herelink_refine)
PSS_INSTANCE_Q = 6.0             # per-instance q at the fitted 5 ms grid (Exp(1) in noise)


def _herelink_refine(x, fs, fc, X, cand, sr):
    ofs1, P = sr["ofs"], sr["P"]
    dur_s = len(x) / fs
    nhf = int(dur_s / 5e-3) + 2
    t_inst = [(cand["tau"] + k * P) / ofs1 for k in range(nhf)]
    t_inst = [t for t in t_inst if t + 70e-6 < dur_s]
    sub = dict(fold_z=round(cand["z"], 2), fold_root=cand["root"], accepted=False)
    # ---- exact carrier and root, then 30.72 MSps working slice
    f_c, u, s_sum = _carrier_search(X, fs, cand["f_off"], t_inst)
    y, ofs = _slice(X, fs, f_c, 9.3e6, 30.72e6)
    tpl = _pss_time(u, ofs)
    inst = _instances(y, ofs, tpl, [t * ofs for t in t_inst], win=int(45e-6 * ofs))
    ok = [(k, d) for k, d in enumerate(inst) if d is not None]
    if len(ok) < 2:
        sub["reject_reason"] = "fewer than 2 PSS instances in the capture"
        return [], sub
    # robust 5 ms grid: phase = median of (t_k - k*P); instances on the grid count
    Pn = 5e-3 * ofs
    phase = float(np.median([d["t"] - k * Pn for k, d in ok if d["q"] > PSS_INSTANCE_Q] or
                            [d["t"] - k * Pn for k, d in ok]))
    good = [(k, d) for k, d in ok if abs(d["t"] - (phase + k * Pn)) < 0.6e-6 * ofs and d["q"] > PSS_INSTANCE_Q]
    sub.update(root=u, nid2=_ROOTS[u], carrier_mhz=round((fc + f_c) / 1e6, 4), pss_instances=len(good),
               pss_total=len(ok), pss_q=[round(d["q"], 1) for _, d in ok])
    feats = {}
    if len(good) >= 3:
        k = np.array([k for k, _ in good], float); t = np.array([d["t"] for _, d in good]) / ofs
        Tp, t0 = np.polyfit(k, t, 1)
        feats["timing_rms_us"] = round(float(np.sqrt(np.mean((t - (t0 + Tp * k)) ** 2))) * 1e6, 3)
        feats["frame_period_ms"] = round(Tp * 1e3, 6)
        feats["clock_ppm"] = round((HERELINK["half_frame_ms"] * 1e-3 / Tp - 1) * 1e6, 2)
    else:
        Tp, t0 = 5e-3, phase / ofs
    pss_idx = [(t0 + kk * Tp) * ofs for kk in range(nhf + 1)]
    pss_idx = [p for p in pss_idx if 2192 < p and p + 2048 < len(y)]
    # ---- SSS -> cell id and 10 ms frame phase
    sss = _sss_decode(y, ofs, _ROOTS[u], pss_idx)
    sss_ok = sss is not None and sss["z"] > 8 and sss["z"] - sss["z_second"] > 3
    sub["sss"] = sss
    # ---- purity: a real PSS symbol owns its 1.08 MHz (power-boosted, nothing else
    # in the centre 6 RB of that symbol), so its normalised correlation tracks the
    # window SNR: q ~ nt*s/(1+s). Another 15 kHz system's sync (DJI OcuSync has a
    # 5 ms-periodic ZC symbol) correlates only partially: q << nt*s/(1+s) at high s.
    pur, s_med = _pss_purity(X, fs, f_c, u, [d["t"] / ofs for _, d in good] or [d["t"] / ofs for _, d in ok])
    sub.update(purity=round(pur, 3), pss_window_snr=round(s_med, 2))
    # Real Herelink keeps purity >= 0.62 down to its sensitivity limit (AWGN sweep, d = 69 dB,
    # s = 0.19); DJI OcuSync's slot ZC reaches 0.10-0.35 (0.351 on the 10 MHz Mini 2 video in
    # D_20260916-171218 at s = 0.4) -- 0.5 splits them with margin on both sides
    if s_med > 0.15 and pur < PSS_PURITY_MIN:
        sub["reject_reason"] = f"PSS purity {pur:.3f} at window SNR {10 * np.log10(s_med):.1f} dB: partial ZC match " \
                               "(another 15 kHz system's sync), not an LTE PSS"
        return [], sub
    # ---- acceptance: SSS decoded in the FDD position, or >= 3 on-grid PSS instances
    if not (sss_ok or len(good) >= 3):
        sub["reject_reason"] = f"no SSS and only {len(good)} on-grid PSS instances"
        return [], sub
    sub["accepted"] = True
    sub["f_off"] = f_c
    t_sfp = t0 - HERELINK["pss_offset_us"] * 1e-6                # start of a PSS-bearing subframe
    if sss_ok:
        t_sf0 = t_sfp if sss["parity"] == 0 else t_sfp - 5e-3
        feats["frame_phase_us"] = round(float(np.mod(t_sf0, 10e-3)) * 1e6, 2)
        feats.update(pci=sss["pci"], nid1=sss["nid1"])
    else:
        t_sf0 = t_sfp
        feats["frame_phase_us"] = round(float(np.mod(t_sf0, 5e-3)) * 1e6, 2)
        feats["frame_phase_ref"] = "5 ms (SSS not decoded: subframe 0/5 ambiguous)"
    feats.update(nid2=_ROOTS[u], pss_root=u)
    # ---- subframe grid (index j mod 5 relative to a PSS-bearing subframe)
    sf_len = 1e-3 * Tp / 5e-3                                    # follow the measured clock
    j0 = int(np.floor(-t_sf0 / sf_len)) - 1
    sfs = [(j % 5, t_sf0 + j * sf_len) for j in range(j0, j0 + int(dur_s / sf_len) + 4)]
    sfs = [(j, t) for j, t in sfs if t + sf_len > 0 and t < dur_s]
    bw = _herelink_bw(y, ofs, [t for j, t in sfs if j in (3, 4) and 0 <= t < dur_s - 1e-3])
    yb, _ = _slice(X, fs, f_c, bw / 2 * 1e6, 30.72e6)
    # per-subframe in-band power; noise reference = the quietest subframe positions
    # (trimmed mean over 10 us blocks: a few-us switching transient at a TDD edge
    # is 40 dB above the floor and would otherwise own a silent subframe)
    pw = []
    kb10 = int(10e-6 * ofs)
    for j, t in sfs:
        a, b = int(max(t, 0) * ofs), int(min(t + sf_len, dur_s) * ofs)
        if b - a < 0.2e-3 * ofs:
            pw.append(np.nan); continue
        e = np.sort((np.abs(yb[a:a + (b - a) // kb10 * kb10]) ** 2).reshape(-1, kb10).mean(1))
        pw.append(float(np.mean(e[:max(1, int(len(e) * 0.95))])))
    pw = np.array(pw); jj = np.array([j for j, _ in sfs])
    full = np.array([(t >= 0 and t + sf_len <= dur_s) for _, t in sfs])
    prof = np.array([np.nanmedian(pw[(jj == j) & full]) if np.any((jj == j) & full) else np.nan for j in range(5)])
    pmin = np.nanmin(prof)
    dof = sf_len * bw * 1e6
    # silent positions = the measured Herelink frame (relative to a PSS subframe:
    # +1, +2 = sf 1,2 / 6,7). The measured profile is reported, not trusted: at low
    # SNR every position looks idle, and a TDD-edge transient can make one look busy
    idle_pos = [1, 2]
    ref = pw[np.isin(jj, idle_pos) & full]
    n0 = float(np.median(ref)) if ref.size else float(pmin)
    sig = n0 / np.sqrt(dof)
    if ref.size >= 4:                                            # interference in the reference widens it
        sig = max(sig, 1.4826 * float(np.median(np.abs(ref - np.median(ref)))))
    zsf = (pw - n0) / sig
    dl_pos = [0, 3, 4]
    sub.update(subframe_profile_db=[round(10 * np.log10(max(p / n0, 1e-3)), 2) for p in prof],
               dl_positions=dl_pos, idle_positions=idle_pos, bw_mhz=bw)
    # ---- fine CFO: CP phase over DL symbols (fractional part around f_c)
    starts = sum((_sym_starts(t, ofs) for j, t in sfs if j in dl_pos and 0 <= t and t + sf_len < dur_s), [])
    cfo_frac, _ = _cfo_cp(y, ofs, starts)
    if cfo_frac is not None and abs(cfo_frac) < 3e3:             # PSS search already within +-1.25 kHz
        f_c = f_c + cfo_frac
    f_abs = fc + f_c
    nominal = round(f_abs / 100e3) * 100e3                       # LTE 100 kHz raster (assumed)
    feats.update(cfo_hz=round(f_abs - nominal, 1), cfo_ppm=round((f_abs - nominal) / nominal * 1e6, 2),
                 nominal_carrier_mhz=nominal / 1e6, carrier_mhz=round(f_abs / 1e6, 5),
                 cfo_ref="nearest 100 kHz raster (Herelink channel plan not public)",
                 symbol_us=round(HERELINK["symbol_us"], 3), tu_us=round(HERELINK["tu_us"], 3),
                 cp_us=round(HERELINK["cp_us"], 3), scs_khz=15.0, bw_mhz=bw, numerology="lte_like",
                 modulation="OFDM (LTE Rel-8 downlink PHY)", pss_instances=len(good))
    feats.setdefault("frame_period_ms", HERELINK["half_frame_ms"])
    # ---- emissions
    out = []
    cen = f_abs / 1e6
    lo, hi = cen - bw / 2, cen + bw / 2
    snr_of = lambda p: float(10 * np.log10(max(p - n0, n0 * 1e-3) / n0))
    conf_link = float(np.clip(0.55 + 0.03 * cand["z"] + (0.2 if sss_ok else 0), 0, 0.99))
    grid_q = {int(round((d["t"] / ofs - t0) / Tp)): d for _, d in good}
    for i, (j, t) in enumerate(sfs):                             # control: PSS/SSS-bearing subframes
        if j != 0:
            continue
        k = int(round((t + HERELINK["pss_offset_us"] * 1e-6 - t0) / Tp))
        d = grid_q.get(k)
        if d is None and not (np.isfinite(zsf[i]) and zsf[i] > 6):
            continue
        f = dict(feats)
        f.update(burst_us=round(sf_len * 1e6, 1), pss_q=round(d["q"], 1) if d else None,
                 subframe=(0 if (k + sss["parity"]) % 2 == 0 else 5) if sss_ok else "0|5",
                 content="PSS+SSS+CRS (+PBCH in subframe 0), no PDSCH observed")
        if d is not None:                                        # per-burst CFO from this PSS (vs nominal)
            f["pss_t_us"] = round(d["t"] / ofs * 1e6, 3)
        out.append(Emission("herelink", "control", round(max(t, 0) * 1e6, 1), round(min(t + sf_len, dur_s) * 1e6, 1),
                            round(lo, 3), round(hi, 3), round(cen, 4),
                            round(snr_of(pw[i]), 1) if np.isfinite(pw[i]) else None,
                            conf_link if d is not None else 0.5, features=f))
    # video = runs of active non-sync DL subframes; uplink = activity in idle positions
    act = [(i, j, t) for i, (j, t) in enumerate(sfs) if j != 0 and np.isfinite(zsf[i]) and zsf[i] > 6]
    runs, cur = [], []
    for a in act:
        if cur and (a[0] != cur[-1][0] + 1 or (a[1] in idle_pos) != (cur[-1][1] in idle_pos)):
            runs.append(cur); cur = []
        cur.append(a)
    if cur:
        runs.append(cur)
    for r in runs:
        is_ul = r[0][1] in idle_pos
        t_a, t_b = _edges(yb, ofs, r[0][2], r[-1][2] + sf_len, n0, dur_s, pad=0.15e-3 if is_ul else 0.03e-3)
        # never run into a neighbouring subframe of the other kind
        t_a = max(t_a, r[0][2] - (150e-6 if is_ul else 12e-6))   # UpPTS sits at the end of the S subframe
        t_b = min(t_b, r[-1][2] + sf_len + 3e-6)
        p = float(np.mean([pw[a[0]] for a in r]))
        if is_ul:
            flo, fhi = _occupied(yb, ofs, t_a, t_b, f_c, fc, n0, lo, hi)
            a_, b_ = int(t_a * ofs), int(t_b * ofs)
            blk = (np.abs(yb[a_:a_ + (b_ - a_) // kb10 * kb10]) ** 2).reshape(-1, kb10).mean(1)
            duty = float(np.mean(blk > 1.5 * n0)) if blk.size else 0.0
            if t_b - t_a < 300e-6 or fhi - flo < 3.0 or duty < 0.6:
                continue                                          # short / narrow: not this link's uplink
            f = dict(frame_period_ms=feats.get("frame_period_ms"), frame_phase_us=feats.get("frame_phase_us"),
                     pci_of_downlink=feats.get("pci"), uplink_unverified=True, duty=round(duty, 2),
                     note="TDD-aligned emission in the air unit's silent subframes (UpPTS+U of an LTE TDD "
                          "config-2 frame); presumed ground-unit uplink, PHY not characterised")
        else:
            f = dict(feats)
            cb, _ = _cfo_cp(y, ofs, sum((_sym_starts(a[2], ofs) for a in r), []))
            if cb is not None and abs(cb) < 3e3:                 # per-burst carrier (fingerprint stability)
                f["cfo_burst_hz"] = round(f_abs - (cfo_frac or 0.0) + cb - nominal, 1)
            flo, fhi = lo, hi
        f.update(burst_us=round((t_b - t_a) * 1e6, 1), subframes=[int(a[1]) for a in r],
                 truncated=bool(t_a <= 1e-6 or t_b >= dur_s - 1e-6))
        out.append(Emission("herelink", "uplink" if is_ul else "video", round(t_a * 1e6, 1), round(t_b * 1e6, 1),
                            round(flo, 3), round(fhi, 3), round((flo + fhi) / 2, 4), round(snr_of(p), 1),
                            0.4 if is_ul else conf_link, features=f))
    # link solid but bursts below per-burst sensitivity: report the expected video
    # subframes (measured Herelink frame: sf 3-4 / 8-9), flagged as predicted
    if not any(e.label == "video" for e in out):
        for j, t in sfs:
            if j == 3 and t >= 0 and t + 2 * sf_len <= dur_s:
                f = dict(feats)
                f.update(burst_us=round(2 * sf_len * 1e6, 1), predicted=True,
                         note="below per-burst sensitivity; position from PSS timing + measured frame")
                out.append(Emission("herelink", "video", round(t * 1e6, 1), round((t + 2 * sf_len) * 1e6, 1),
                                    round(lo, 3), round(hi, 3), round(cen, 4), None, 0.3, features=f))
    sub["n0"] = n0
    return out, sub


def _herelink_bw(y, ofs, t_list):
    """20 vs 10 MHz mode from the PSD shape of the data subframes: power in the
    5-8.8 MHz ring vs the inner 4 MHz (both above the 9.6-12 MHz out-of-band level)."""
    segs = [y[int(t * ofs):int(t * ofs) + int(0.9e-3 * ofs)] for t in t_list[:6]]
    segs = [s for s in segs if len(s) >= 4096]
    if not segs:
        return HERELINK["bw_mhz"]
    nf = 2048
    S = np.mean([np.mean(np.abs(np.fft.fft(s[:len(s) // nf * nf].reshape(-1, nf), axis=1)) ** 2, 0) for s in segs], 0)
    f = np.fft.fftfreq(nf, 1 / ofs)
    inner = np.mean(S[np.abs(f) < 4e6]); ring = np.mean(S[(np.abs(f) > 5e6) & (np.abs(f) < 8.8e6)])
    edge = np.mean(S[(np.abs(f) > 9.6e6) & (np.abs(f) < 12e6)])
    return 18.0 if (ring - edge) > 0.3 * (inner - edge) else 9.0


def _edges(yb, ofs, t_a, t_b, n0, dur_s, pad):
    """Snap a subframe run to the actual burst edges at 10 us resolution when the
    burst is strong enough to see them; otherwise keep subframe boundaries."""
    a = int(max(t_a - pad, 0) * ofs); b = int(min(t_b + pad, dur_s) * ofs)
    k = int(10e-6 * ofs)
    if b - a < 4 * k:
        return max(t_a, 0), min(t_b, dur_s)
    e = np.abs(yb[a:b - (b - a) % k]) ** 2
    e = e.reshape(-1, k).mean(1)
    on = e > 4 * n0                                              # ~5 dB above floor
    if on.sum() < 3:
        return max(t_a, 0), min(t_b, dur_s)
    w = np.where(on)[0]
    return (a + w[0] * k) / ofs, min((a + (w[-1] + 1) * k) / ofs, dur_s)


def _occupied(yb, ofs, t_a, t_b, f_c, fc, n0, lo, hi):
    """Occupied band of a (weak) burst: PSD bins > 3 dB over the floor."""
    a, b = int(t_a * ofs), int(t_b * ofs)
    nf = 1024
    s = yb[a:b]
    if len(s) < 4 * nf:
        return lo, hi
    S = np.mean(np.abs(np.fft.fftshift(np.fft.fft(s[:len(s) // nf * nf].reshape(-1, nf) * np.hanning(nf), axis=1),
                                       axes=1)) ** 2, 0) / np.sum(np.hanning(nf) ** 2)
    S = ndimage.uniform_filter1d(S, 9)
    f = (np.arange(nf) - nf // 2) / nf * ofs
    m = S > 2 * n0
    if m.sum() < 5:
        return lo, hi
    w = np.where(m)[0]
    return (fc + f_c + f[w[0]]) / 1e6, (fc + f_c + f[w[-1]]) / 1e6


# ============================================================== generic: blind CP-OFDM links
def generic_detect(x, fs, fc, X, ctx, info, exclude=()):
    """Every CP-OFDM link with Tu >= 10 us, whatever its make (blind.py): one Emission per burst,
    family digital_vtx, label ofdm_link_<SCS>kHz; the numerology table names what it resembles. Then
    the hopping narrowband CP-OFDM links (RC uplinks, blind_hops): label ofdm_hop_<SCS>kHz.
    A link over more than half of whose band an `exclude` emission of the same Tu lies (Herelink, found by its
    PSS: 66.67 us; HDZero, found by its pilots: 37.93 us) is left to it."""
    usable = min(fs, 56e6) / 2
    noise = C.noise_floor_psd(x)["floor_power_per_sample"]
    log = info.setdefault("generic_log", []).append
    tf = ctx.get("tfpower")
    if tf is None:
        tf = ctx["tfpower"] = BU.TFPower(x, fs)
    # a candidate inside a licensed cellular band, and whose IQ mirror is too, is the operator's carrier: the plugin
    # drops it, so it is not confirmed here (PLAN G8: the 1.6-2.33 GHz search covers the B3 / B1 bands)
    def cellular(off):
        return BI.cellular_at(fc + off) and BI.cellular_at(fc - off)

    sl = _SliceCache(X, fs)
    links = BL.blind_links(sl, fs, fc, usable, noise, log=log, tf=tf,
                           keep_slice_snr=ST.SNR_MIN, skip=cellular)
    out = []
    info["generic_links"] = []
    for lk in links:
        lo, hi = (fc + lk["lo"]) / 1e6, (fc + lk["hi"]) / 1e6
        tu_us, ts_us = lk["tu"] * 1e6, lk["ts"] * 1e6
        if any(abs(tu_us / e.features.get("tu_us", 66.667) - 1) < 0.03
               and min(hi, e.f_hi_mhz) - max(lo, e.f_lo_mhz) > 0.5 * (hi - lo) for e in exclude):
            continue
        cp_us = ts_us - tu_us
        scs = 1e3 / tu_us
        name, matches = classify_numerology(tu_us, cp_us, lk["bw"] / 1e6)
        role = next((n["role"] for n in KNOWN_NUMEROLOGIES if n["name"] == name), "unknown")
        feats = dict(numerology=name, matches=matches, role=role, tu_us=round(tu_us, 3), symbol_us=round(ts_us, 3),
                     cp_us=round(cp_us, 3), cp_fraction=round(cp_us / tu_us, 4), scs_khz=round(scs, 3),
                     bw_mhz=round(lk["bw"] / 1e6, 2), carrier_mhz=round((fc + lk["cen"]) / 1e6, 3),
                     n_bursts=len(lk["runs"]), duty_cycle=round(lk["duty"], 3),
                     frame_period_ms=None if lk["frame_s"] is None else round(lk["frame_s"] * 1e3, 3),
                     frame_contrast=round(lk["frame_contrast"], 1), symbol_grids=lk["n_grids"],
                     line_z=round(lk["line_z"], 1), strong_windows=lk["win_z"], acf_z=round(lk["acf_z"], 1),
                     link_snr_db=round(lk["snr_db"], 1), cfo_frac_hz=round(lk["cfo_frac_hz"], 1),
                     cfo_note="CP phase: carrier offset modulo the subcarrier spacing, vs band centre",
                     band_context=BI.uav_link_band(fc + lk["cen"]))
        # transmit streams per subcarrier (1 TX vs 2-TX STBC/SFBC), from the link's slice at N / Tu samples/s
        yb, ofs2 = lk.pop("slice", (None, None))
        feats["tx_antennas"] = ST.tx_streams(yb, ofs2, lk["tu"], lk["ts"], lk["bw"], lk["runs"], lk["snr_db"],
                                             len(x) / fs)
        # research catalog: which system / model the measured PHY most resembles
        feats["system_candidates"] = CAT.candidates(tu_us, cp_us / tu_us, lk["bw"] / 1e6, (fc + lk["cen"]) / 1e6,
                                                    feats["tx_antennas"].get("tx_antennas"),
                                                    None if lk["frame_s"] is None else lk["frame_s"] * 1e3)
        info["generic_links"].append({k: v for k, v in feats.items() if k != "cfo_note"})
        conf = 0.45 if role == "unknown" else 0.6           # inferred from the PHY alone: "-like" / unknown
        lab = f"ofdm_link_{scs:.0f}kHz" if abs(scs - round(scs)) < 0.05 else f"ofdm_link_{scs:.2f}kHz"
        for t0, t1, nsym in lk["runs"]:
            f = dict(feats, burst_us=round((t1 - t0) * 1e6, 1), n_symbols=int(nsym))
            out.append(Emission("digital_vtx", lab, round(t0 * 1e6, 1), round(t1 * 1e6, 1), round(lo, 3), round(hi, 3),
                                round((lo + hi) / 2, 3), round(lk["snr_db"], 1), conf, features=f))
    # hopping narrowband CP-OFDM links (RC uplinks), searched at the numerologies above and the drone
    # families': one Emission per hop, family digital_vtx, label ofdm_hop_<SCS>kHz
    info["generic_hops"] = []
    for hk in BL.blind_hops(sl, fs, fc, usable, links, log=log, skip=cellular):
        tu_us, ts_us = hk["tu"] * 1e6, hk["ts"] * 1e6
        cp_us = ts_us - tu_us
        scs = 1e3 / tu_us
        name, matches = classify_numerology(tu_us, cp_us, hk["bw"] / 1e6)
        role = next((n["role"] for n in KNOWN_NUMEROLOGIES if n["name"] == name), "unknown")
        feats = dict(numerology=name, matches=matches, role=role, tu_us=round(tu_us, 3), symbol_us=round(ts_us, 3),
                     cp_us=round(cp_us, 3), cp_fraction=round(cp_us / tu_us, 4), scs_khz=round(scs, 3),
                     hop_bw_mhz=round(hk["bw"] / 1e6, 2), hop_us=round(hk["dur_s"] * 1e6, 1), n_hops=len(hk["hops"]),
                     n_channels=hk["n_channels"],
                     hop_period_ms=None if hk["frame_s"] is None else round(hk["frame_s"] * 1e3, 3),
                     paired_link_mhz=[round((fc + c) / 1e6, 3) for c in hk["paired"]],
                     band_context=BI.uav_link_band(fc + hk["hops"][0]["channel"]) if hk["hops"] else None)
        info["generic_hops"].append(dict(feats))
        conf = 0.6 if role in ("drone", "infrastructure") else 0.45
        lab = f"ofdm_hop_{scs:.0f}kHz" if abs(scs - round(scs)) < 0.05 else f"ofdm_hop_{scs:.2f}kHz"
        for h in hk["hops"]:
            lo, hi = (fc + h["lo"]) / 1e6, (fc + h["hi"]) / 1e6
            snr = 10 * np.log10(max(10 ** (h["contrast"] / 10) - 1, 1e-3))
            f = dict(feats, channel_mhz=round((fc + h["channel"]) / 1e6, 2), bw_mhz=round(hi - lo, 2),
                     burst_us=round((h["t1"] - h["t0"]) * 1e6, 1), n_symbols=int(round((h["t1"] - h["t0"]) / hk["ts"])),
                     fold_z=round(h["z"], 1), contrast_db=round(h["contrast"], 1), tracked=bool(h["tracked"]))
            out.append(Emission("digital_vtx", lab, round(h["t0"] * 1e6, 1), round(h["t1"] * 1e6, 1), round(lo, 3),
                                round(hi, 3), round((lo + hi) / 2, 3), round(snr, 1), conf, features=f))
    return out


# ============================================================== entry point
def detect(x, fs, fc, ctx=None):
    ctx = ctx if ctx is not None else {}
    info = ctx.setdefault("digital_vtx_info", {})
    x = np.asarray(x)
    spurs = ctx.get("spurs")
    if spurs is None:
        spurs = find_spurs(x, fs, fc)
    X = _notch(_fft_of(x, ctx), fs, spurs)
    out = []
    if ctx.get("herelink", True):
        out += herelink_detect(x, fs, fc, X, ctx, info)
    if ctx.get("hdzero", True):
        out += hdzero_detect(x, fs, fc, X, ctx, info)
    if ctx.get("generic", True):
        excl = [e for e in out if e.family in ("herelink", "hdzero")]
        out += generic_detect(x, fs, fc, X, ctx, info, exclude=excl)
    out.sort(key=lambda e: (e.start_us, e.center_mhz))
    return out


if __name__ == "__main__":
    import sys, json, time
    from collections import Counter
    t = time.time()
    x, fs, fc, meta = C.load_capture(sys.argv[1])
    ctx = {}
    em = detect(x, fs, fc, ctx)
    print(f"{len(em)} emissions in {time.time() - t:.1f}s: {dict(Counter((e.family, e.label) for e in em))}")
    for e in em:
        print(f"  {e.family:11s} {e.label:24s} {e.start_us:9.1f}-{e.end_us:9.1f} us  {e.f_lo_mhz:.2f}-{e.f_hi_mhz:.2f} MHz"
              f"  snr {e.snr_db}  conf {e.confidence:.2f}")
    print(json.dumps({k: v for k, v in ctx["digital_vtx_info"].items()}, default=str, indent=1)[:4000])
    if em:
        print(json.dumps(em[0].features, default=str, indent=1))
