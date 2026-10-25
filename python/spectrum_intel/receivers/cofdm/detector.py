"""
detector.py -- DVB-T-framed COFDM links (drone / UAV video downlinks, DVB-T broadcast) from IQ.

Research (research/drone_transmitters/06_dvbt_cofdm.md; ETSI EN 300 744), synthetic validation only.

1. candidates: the block spectrum (Welch) above floor + 2 dB in contiguous runs whose width is 0.951 B (+-10 %) for a
   channel bandwidth B in BANDWIDTHS (8 / 7 / 6 / 5 MHz standard; 4 / 2.5 / 2 MHz scaled COFDM modules), and present
   in >= 80 % of 1 ms windows (DVB-T is continuous, unlike TDD video links);
2. numerology: the band sliced at the elementary rate 8 B / 7 MHz; the CP metric (lag N) folded on the symbol period
   N (1 + g) for 2K / 8K x GI 1/4 .. 1/32; the best fold peak (>= PEAK_MIN, >= CONTRAST x its median) wins;
3. carrier: fractional offset from the CP phase; integer offset from the 4/3-boosted continual pilots (power summed
   over the symbols at the 45 / 177 pilot carriers, every shift within +-N/8);
4. TPS: DBPSK across consecutive symbols on the 17 / 68 TPS carriers (sum of Re{z_l z*_(l-1)}: the channel cancels),
   sync word (or its inverse) -> frame start, BCH(67,53) check -> mode, guard interval, constellation, hierarchy,
   HP / LP code rate, frame number, cell-ID byte. One TPS frame needs 68 symbols (2K at 8 MHz: 17 ms).
"""
from __future__ import annotations

import numpy as np
from scipy import signal

from spectrum_intel.core.emission import Emission

from . import dvbt2_p1 as P1
from . import tables as T

BANDWIDTHS = (8.0, 7.0, 6.0, 5.0, 4.0, 2.5, 2.0)
ON_RATIO = 1.6                    # +2 dB over the 20th-percentile floor (Welch, ~256 segments: ~6 % spread)
OCC = 0.951                       # occupied / channel bandwidth ((K-1) / Tu)
PEAK_MIN, CONTRAST = 0.2, 3.0
MIN_SYMS = 8
PILOT_DB, PILOT_MARGIN = 1.5, 1.2   # continual-pilot shift: >= 1.5 dB over the median shift and >= 1.2 x the
                                     # runner-up (shifts by 3k partly align the scattered-pilot raster)


def _slice(X, fs, f_off, half_bw, out_fs):
    """Brick-wall slice of the capture FFT X around f_off, at ~out_fs; returns (y, exact rate)."""
    L = len(X)
    M = int(round(L * out_fs / fs))
    k0 = int(round(f_off / fs * L))
    kk = int(half_bw / fs * L)
    m = np.arange(-kk, kk + 1)
    Y = np.zeros(M, complex)
    Y[m % M] = X[(k0 + m) % L]
    return np.fft.ifft(Y) * (M / L), M * fs / L


def candidates(x, fs, fc, nper=1 << 14):
    """[(centre offset Hz, B MHz, width Hz, level dB over floor)] of flat, continuous blocks."""
    f, P = signal.welch(x, fs, nperseg=min(nper, len(x)), return_onesided=False)          # ratios only: scaling-free
    f, P = np.fft.fftshift(f), np.fft.fftshift(P)
    fl = np.percentile(P, 20)
    on = P > ON_RATIO * fl
    on = np.convolve(on.astype(float), np.ones(5), "same") > 0                  # bridge 1-2 bin dips
    df = f[1] - f[0]
    out = []
    edges = np.flatnonzero(np.diff(np.r_[0, on.astype(int), 0]))
    for a, b in zip(edges[0::2], edges[1::2]):
        w = (b - a) * df
        for B in BANDWIDTHS:
            if abs(w - OCC * B * 1e6) <= 0.1 * OCC * B * 1e6:
                lev = 10 * np.log10(np.median(P[a:b]) / fl)
                out.append(((f[a] + f[b - 1]) / 2, B, w, float(lev)))
                break
    return out


def _cp_fold(y, N, Ng):
    """Normalised CP metric folded on the symbol period: (peak, contrast, phase, complex sum at the peak, nsym)."""
    L = N + Ng
    nsym = (len(y) - N) // L
    if nsym < MIN_SYMS:
        return None
    c = y[:-N] * np.conj(y[N:])
    e = (np.abs(y[:-N]) ** 2 + np.abs(y[N:]) ** 2) / 2
    cs = np.concatenate(([0], np.cumsum(c)))
    es = np.concatenate(([0], np.cumsum(e)))
    mc = cs[Ng:] - cs[:-Ng]
    me = es[Ng:] - es[:-Ng]
    n = nsym * L
    if len(mc) < n:
        nsym -= 1
        n = nsym * L
    F = mc[:n].reshape(nsym, L).sum(0)
    E = me[:n].reshape(nsym, L).sum(0) + 1e-30
    prof = np.abs(F) / E
    p = int(np.argmax(prof))
    return float(prof[p]), float(prof[p] / (np.median(prof) + 1e-12)), p, F[p], nsym


def numerology(y, ofs, B):
    """Best (mode, gi, peak, contrast, phase, csum, nsym, N, Ng) of the 2K / 8K x GI hypotheses, or None."""
    best = None
    for mode, M in T.MODES.items():
        N = M["N"]
        for gi, d in T.GUARD.items():
            r = _cp_fold(y, N, N // d)
            if r is None:
                continue
            peak, con, ph, csum, nsym = r
            if peak >= PEAK_MIN and con >= CONTRAST and (best is None or peak > best[2]):
                best = (mode, gi, peak, con, ph, csum, nsym, N, N // d)
    return best


def tps_decode(bits):
    """DBPSK bits b[l] = s_l of consecutive symbols -> decoded TPS dict (first frame whose sync word matches and
    whose BCH(67,53) passes, else the first sync match with bch_ok False), or None."""
    sync = np.array(T.SYNC, np.uint8)
    res = None
    for p in range(1, len(bits) - 66):
        w = bits[p:p + 16]
        if not (np.array_equal(w, sync) or np.array_equal(w, 1 - sync)):
            continue
        s = [0] + list(int(b) for b in bits[p:p + 67])             # s0..s67 (s0 unknown, reference only)
        info, par = s[1:54], s[54:68]
        ok = T.bch_parity(info) == par

        def v(a, n):
            return int("".join(map(str, s[a:a + n])), 2)
        r = dict(bch_ok=ok, frame_number=v(23, 2) + 1, length_indicator=T.LENGTH.get(v(17, 6), f"{v(17, 6):06b}"),
                 constellation=T.CONSTELLATION[v(25, 2)], hierarchy=T.HIERARCHY.get(v(27, 3), f"{v(27, 3):03b}"),
                 code_rate_hp=T.CODE_RATE.get(v(30, 3), "reserved"), code_rate_lp=T.CODE_RATE.get(v(33, 3), "reserved"),
                 guard_interval=T.GI_TPS[v(36, 2)], mode=T.MODE_TPS[v(38, 2)], frame_start_symbol=p - 1)
        if v(17, 6) == 0b011111:
            r["cell_id_byte"] = f"0x{v(40, 8):02X}" + (" (high)" if r["frame_number"] in (1, 3) else " (low)")
        if ok:
            return r
        res = res or r
    return res


def analyse(x, fs, fc, X, cand):
    """One candidate -> features dict, or None."""
    off, B, width, lev = cand
    clk = 8e6 * B / 7
    y, ofs = _slice(X, fs, off, 0.5 * B * 1e6, clk)
    nm = numerology(y, ofs, B)
    if nm is None:
        return None
    mode, gi, peak, con, ph, csum, nsym, N, Ng = nm
    M = T.MODES[mode]
    K = M["K"]
    L = N + Ng
    eps = -np.angle(csum) / (2 * np.pi * N)                       # cycles / sample: y ~ exp(j 2 pi eps n)
    yc = y * np.exp(-2j * np.pi * eps * np.arange(len(y)))
    starts = ph + np.arange(nsym) * L
    starts = starts[starts + Ng + N <= len(yc)]
    off_w = Ng // 2
    W = starts[:, None] + off_w + np.arange(N)[None, :]
    Z = np.fft.fft(yc[W], axis=1)
    kc = (K - 1) // 2
    cp = np.array(M["cp"]) - kc
    Pm = np.mean(np.abs(Z) ** 2, 0)
    S = N // 8
    shifts = np.arange(-S, S + 1)
    sc = np.array([Pm[(cp + s) % N].mean() for s in shifts]) / (Pm.mean() + 1e-30)
    j = int(np.argmax(sc))
    si = int(shifts[j])
    pilot_db = float(10 * np.log10(sc[j] / (np.median(sc) + 1e-30)))   # 4/3 amplitude: +2.5 dB at the right shift
    mad = 1.4826 * np.median(np.abs(sc - np.median(sc))) + 1e-12
    pilot_z = float((sc[j] - np.median(sc)) / mad)
    pilot_margin = float(sc[j] / (np.max(np.delete(sc, j)) + 1e-30))
    scs = ofs / N
    feats = dict(standard="DVB-T (EN 300 744) framing", channel_bw_mhz=B, occupied_mhz=round(width / 1e6, 3),
                 mode=mode, guard_interval=gi, fft_size=N, elementary_rate_mhz=round(ofs / 1e6, 6),
                 scs_hz=round(scs, 2), tu_us=round(N / ofs * 1e6, 3), symbol_us=round(L / ofs * 1e6, 3),
                 cp_peak=round(peak, 3), cp_contrast=round(con, 1), symbols=int(len(starts)),
                 carrier_mhz=round((fc + off + si * scs + eps * ofs) / 1e6, 6), integer_offset_carriers=si,
                 continual_pilot_db=round(pilot_db, 1), continual_pilot_z=round(pilot_z, 1),
                 continual_pilot_margin=round(pilot_margin, 2), level_db=round(lev, 1))
    if pilot_db < PILOT_DB or pilot_margin < PILOT_MARGIN:
        # DVB-T framing always carries the continual pilots: without them a CP-fold peak at a DVB-T period is
        # some other OFDM / periodic signal (real captures: DJI OcuSync, 2.4 GHz ISM traffic), not reported
        return None
    tps = np.array(M["tps"]) - kc
    zt = Z[:, (tps + si) % N]
    d = np.sum(np.real(zt[1:] * np.conj(zt[:-1])), 1)
    bits = np.r_[0, (d < 0).astype(np.uint8)]                       # bits[l] = s_l relative to symbol l-1
    r = tps_decode(bits)
    feats["tps"] = r
    if r is None:
        feats["tps_note"] = (f"no TPS sync word in {len(bits)} symbols (a TPS frame is 68 symbols, "
                             f"{68 * L / ofs * 1e3:.1f} ms)")
    return feats


P1_THR = 0.35                     # P1 correlation: noise / DVB-T / 1K CP-OFDM <= 0.11, a P1 at -3 dB ~0.55
# a real P1 decodes its 127-bit S1 and 256-bit S2 codewords with margins >= 0.74 even at -6 dB; correlation peaks from
# ADS-B pulses, analog video and clipped captures passed P1_THR but decoded with a margin <= 0.13 on S1 or S2
# (regression 2026-10-22): both margins >= P1_MARGIN and a defined S1 are required
P1_MARGIN = 0.4


def analyse_t2(fs, fc, X, cand):
    """DVB-T2 P1 preamble symbols in one candidate band (research, PLAN D5; dvbt2_p1) -> features dict or None.
    A T2 frame lasts up to 250 ms, so a 34 ms block may hold no P1: then nothing is reported for T2."""
    off, B, width, lev = cand
    y, ofs = _slice(X, fs, off, 0.5 * B * 1e6, 8e6 * B / 7)
    found = [f for f in P1.find(y, P1_THR) if f[2]["s1_margin"] >= P1_MARGIN and f[2]["s2_margin"] >= P1_MARGIN
             and f[2]["preamble"] != "reserved"]
    if not found:
        return None
    s0, m, d = found[0]
    scs = ofs / P1.N1
    return dict(standard="DVB-T2 (EN 302 755) P1 preamble", channel_bw_mhz=B, occupied_mhz=round(width / 1e6, 3),
                elementary_rate_mhz=round(ofs / 1e6, 6), p1_count=len(found), p1_metric=round(m, 3),
                p1_start_us=round(s0 / ofs * 1e6, 1), p1=d, carrier_mhz=round((fc + off + d["carrier_offset"] * scs) / 1e6, 6),
                level_db=round(lev, 1))


def detect(x, fs, fc, ctx=None):
    x = np.asarray(x)
    X = np.fft.fft(x.astype(np.complex128))
    out = []
    for cand in candidates(x, fs, fc):
        t2 = analyse_t2(fs, fc, X, cand)
        if t2:
            lo = (fc + cand[0] - cand[2] / 2) / 1e6
            hi = (fc + cand[0] + cand[2] / 2) / 1e6
            ok = t2["p1"]["s1_margin"] > 0.2 and t2["p1"]["s2_margin"] > 0.2
            out.append(Emission("cofdm", "dvbt2" if ok else "dvbt2_like", 0.0, round(len(x) / fs * 1e6, 1),
                                round(lo, 3), round(hi, 3), round((lo + hi) / 2, 3), t2["level_db"],
                                0.95 if ok else 0.7, None, t2))
            continue
        f = analyse(x, fs, fc, X, cand)
        if f is None:
            continue
        lo = (fc + cand[0] - cand[2] / 2) / 1e6
        hi = (fc + cand[0] + cand[2] / 2) / 1e6
        tps = f.get("tps")
        conf = 0.95 if (tps or {}).get("bch_ok") else 0.75 if tps else 0.6
        lab = "dvbt" if (tps or {}).get("bch_ok") else "dvbt_like"
        out.append(Emission("cofdm", lab, 0.0, round(len(x) / fs * 1e6, 1), round(lo, 3), round(hi, 3),
                            round((lo + hi) / 2, 3), f["level_db"], conf, None, f))
    return out
