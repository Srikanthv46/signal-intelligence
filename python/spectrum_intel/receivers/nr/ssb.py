"""5G NR SSB receiver: GSCN search, PSS/SSS (differential), DMRS i_SSB, PBCH LLRs, polar decode
(py3gpp), CRC-24C, MIB fields."""
from __future__ import annotations

from functools import lru_cache

import numpy as np
import scipy.fft as sfft

from spectrum_intel.core.seq3gpp import gold, qpsk


# band key -> (name, SSB SCS)
NR_BANDS = {"n78": ("n78", 30e3), "b28_dl": ("n28", 15e3), "b40": ("n40", 30e3), "b41": ("n41", 30e3),
            "b3_dl": ("n3", 15e3), "b1_dl": ("n1", 15e3), "b5_dl": ("n5", 15e3), "b8_dl": ("n8", 15e3)}
NFFT = 256


def _bands_for(key):
    from spectrum_intel.core.bands_india import BANDS
    return [(b[0], b[1]) for b in BANDS if b[2] == key]


def _gscn(lo, hi):
    out = []
    if hi > 3000e6:
        n0, n1 = int(np.ceil((max(lo, 3000e6) - 3000e6) / 1.44e6)), int((hi - 3000e6) / 1.44e6)
        out += [3000e6 + n * 1.44e6 for n in range(n0, n1 + 1)]
    if lo < 3000e6:
        for n in range(int(lo / 1.2e6), int(min(hi, 3000e6) / 1.2e6) + 1):
            for m in (50e3, 150e3, 250e3):
                f = n * 1.2e6 + m
                if lo <= f <= hi:
                    out.append(f)
    return out


@lru_cache(maxsize=None)
def _pss_seq(nid2):
    from py3gpp.nrPSS import nrPSS
    return _frozen(np.asarray(nrPSS(nid2), float))


@lru_cache(maxsize=1)
def _sss_mseqs():
    from py3gpp.helper import _calc_m_seq
    c = [1, 0, 0, 0, 0, 0, 0]
    return _calc_m_seq(7, c, [0, 4]), _calc_m_seq(7, c, [0, 1])


def _sss_seq(pci):
    """py3gpp.nrSSS(pci): the two m-sequences are the same for every PCI, only rotated."""
    from py3gpp.nrSSS import nrSSSm0, nrSSSm1
    m0, m1 = _sss_mseqs()
    return np.asarray((1 - 2 * np.roll(m0, -nrSSSm0(pci))) * (1 - 2 * np.roll(m1, -nrSSSm1(pci))), float)


def _frozen(a):
    a.flags.writeable = False
    return a


# The SSB search asks for the same Gold sequences (per PCI and i_SSB) thousands of times
# per capture, and py3gpp builds each bit by bit in Python: use core.seq3gpp.gold (bit-exact,
# vectorised) and cache. A Gold sequence of a given c_init is prefix-stable, so one long
# run per c_init serves every length.
_GOLD: dict = {}


def _gold(c_init, n):
    g = _GOLD.get(c_init)
    if g is None or len(g) < n:
        g = _GOLD[c_init] = _frozen(gold(c_init, max(n, 8 * 864)).astype(np.int32))
    return g[:n]


def _pbch_prbs(pci, v, n):
    """py3gpp.nrPBCHPRBS(pci, v, n), cached."""
    return _gold(pci, n + v * n)[v * n:][:n]


@lru_cache(maxsize=None)
def _pbch_dmrs(pci, ib):
    """py3gpp.nrPBCHDMRS(pci, ib): QPSK of the Gold sequence for c_init (38.211 7.4.1.4.1), cached."""
    c_init = 2 ** 11 * (ib + 1) * (pci // 4 + 1) + 2 ** 6 * (ib + 1) + (pci % 4)
    return _frozen(qpsk(gold(c_init, 288)))


@lru_cache(maxsize=None)
def _pss_templates(scs, nf):
    """Conjugate spectra of the three time-domain PSS templates (energy-normalised)."""
    k = np.arange(127) - 64                         # PSS on SSB subcarriers 56..182 -> -64..62 around k=120
    t = np.arange(NFFT) / (NFFT * scs)
    T = []
    for nid2 in range(3):
        x = (_pss_seq(nid2)[None, :] * np.exp(2j * np.pi * k[None, :] * scs * t[:, None])).sum(1)
        T.append(np.conj(np.fft.fft(x / np.linalg.norm(x), nf)))
    return _frozen(np.array(T))


class _SSBSlice:
    def __init__(self, ch, f_ss, scs):
        self.scs, self.ofs = scs, NFFT * scs
        self.cp = int(round(144 / 2048 * NFFT))
        self.y, self.ofs = ch.slice(f_ss, 0.55 * 240 * scs, self.ofs)
        self.L = len(self.y)
        # circular correlation over nf >= L: the lags used (< L - NFFT) never wrap, so no zero padding is needed
        # (half the transform length of the old next power of two above L + NFFT; PLAN P1)
        self.nf = sfft.next_fast_len(self.L, False)
        self.Y = np.fft.fft(self.y, self.nf)
        self.e = np.convolve(np.abs(self.y) ** 2, np.ones(NFFT), "valid")[:self.L - NFFT] + 1e-30
        self.T = _pss_templates(scs, self.nf) if self.ofs == NFFT * scs else self._templates()

    def _templates(self):
        k = np.arange(127) - 64
        t = np.arange(NFFT) / self.ofs
        T = []
        for nid2 in range(3):
            x = (_pss_seq(nid2)[None, :] * np.exp(2j * np.pi * k[None, :] * self.scs * t[:, None])).sum(1)
            T.append(np.conj(np.fft.fft(x / np.linalg.norm(x), self.nf)))
        return np.array(T)

    def best(self, offsets):
        b = None
        n = self.L - NFFT
        for df in offsets:
            sh = int(round(df / self.ofs * self.nf))
            Ys = np.roll(self.Y, -sh)
            M = np.abs(sfft.ifft(Ys[None, :] * self.T, axis=1, workers=-1)[:, :n]) ** 2 / self.e
            for nid2 in range(3):
                m = M[nid2]
                k = int(np.argmax(m))
                med = np.median(m)
                mad = 1.4826 * np.median(np.abs(m - med)) + 1e-12
                z = (m[k] - med) / mad
                if b is None or z > b["z"]:
                    b = dict(z=float(z), nid2=nid2, pos=k, df=sh * self.ofs / self.nf, metric=float(m[k]), m=m)
        return b

    def grid(self, df, pos):
        """4 SSB symbols x 240 subcarriers (frequency-domain, SSB k=0..239)."""
        y = self.y * np.exp(-2j * np.pi * df * np.arange(self.L) / self.ofs)
        G = np.zeros((4, 240), complex)
        for l in range(4):
            s = pos + l * (NFFT + self.cp)
            if s + NFFT > self.L:
                return None
            X = np.fft.fft(y[s:s + NFFT])
            k = np.arange(240) - 120
            G[l] = X[k % NFFT]
        return G


def _dmrs_k(pci):
    v = pci % 4
    k13 = np.arange(v, 240, 4)
    k2 = np.r_[np.arange(v, 48, 4), np.arange(192 + v, 240, 4)]
    return k13, k2


def _pbch_k(pci):
    k13, k2 = _dmrs_k(pci)
    d13 = np.setdiff1d(np.arange(240), k13)
    d2 = np.setdiff1d(np.r_[np.arange(48), np.arange(192, 240)], k2)
    return d13, d2


def _smooth(h, k=7):
    w = np.ones(k) / k
    return np.convolve(h.real, w, "same") + 1j * np.convolve(h.imag, w, "same")


def _peaks(m, frac=0.2, sep=NFFT, max_n=24):
    """PSS instances: local maxima above frac*max, >= one symbol apart."""
    order = np.argsort(m)[::-1]
    top = m[order[0]]
    out = []
    for i in order[:20000]:
        if m[i] < frac * top or len(out) >= max_n:
            break
        if all(abs(i - j) > sep for j in out):
            out.append(int(i))
    return sorted(out)


_SSS_BANK = {}


def _sss_bank(nid2):
    if nid2 not in _SSS_BANK:
        S = np.array([_sss_seq(3 * n1 + nid2) for n1 in range(336)])
        _SSS_BANK[nid2] = S[:, :-1] * S[:, 1:]
    return _SSS_BANK[nid2]


def _sss_scores(G, nid2):
    """Differential SSS metric Re sum_k Y_k Y*_{k+1} s_k s_{k+1}: needs no channel
    estimate. Synchronised TDD cells transmit SSBs simultaneously on the same
    GSCN, so a PSS-based channel estimate is a sum of several cells' channels
    and coherent SSS detection fails; the differential metric does not."""
    Y = G[2, 56:183]
    return np.real(_sss_bank(nid2) @ (Y[:-1] * np.conj(Y[1:])))


def _pbch_llr(G, pci):
    """i_SSB (DMRS) and descrambled-ready PBCH LLRs (864) for one SSB."""
    k13, k2 = _dmrs_k(pci)
    rx = np.r_[G[1, k13], G[2, k2], G[3, k13]]
    best = None
    for ib in range(8):
        d = _pbch_dmrs(pci, ib)
        h = rx * np.conj(d)
        met = np.abs(np.sum(h[1:] * np.conj(h[:-1])))
        if best is None or met > best[0]:
            best = (met, ib, h)
    h = best[2]
    kk = np.arange(240)
    h1, h3 = _smooth(h[:60], 5), _smooth(h[84:], 5)
    H1 = np.interp(kk, k13, h1.real) + 1j * np.interp(kk, k13, h1.imag)
    H3 = np.interp(kk, k13, h3.real) + 1j * np.interp(kk, k13, h3.imag)
    H2 = 0.5 * (H1 + H3)
    d13, d2 = _pbch_k(pci)
    Hd = np.r_[H1[d13], H2[d2], H3[d13]]
    rxd = np.r_[G[1, d13], G[2, d2], G[3, d13]]
    z = rxd * np.conj(Hd)                       # MRC-weighted (|H|^2 scaling)
    llr = np.empty(864)
    llr[0::2], llr[1::2] = z.real, z.imag
    return best[1], llr / (np.mean(np.abs(Hd) ** 2) + 1e-12)


def decode_ssb_multi(sl, b):
    """Combine every SSB of the detected cell in the capture: SSS scores summed
    over all PSS instances; PBCH LLRs (each descrambled with its own i_SSB)
    summed per half-frame, then polar decode + CRC."""
    shifted = sl.y * np.exp(-2j * np.pi * b["df"] * np.arange(sl.L) / sl.ofs)
    sh = int(round(b["df"] / sl.ofs * sl.nf))
    m = np.abs(np.fft.ifft(np.roll(sl.Y, -sh) * sl.T[b["nid2"]])[:sl.L - NFFT]) ** 2 / sl.e
    inst = _peaks(m)
    grids = []
    for p in inst:
        G = sl.grid(b["df"], p)
        if G is not None:
            grids.append((p, G))
    info = {"n_ssb_instances": len(grids)}
    if not grids:
        return None, dict(info, crc_ok=False, sss_z=0, sss_z2=0, nid1=None, ibar_ssb=None)
    # (N_ID2, N_ID1) chosen jointly: the strongest PSS may belong to a different
    # co-channel cell than the strongest SSS
    cells = []
    for n2 in range(3):
        acc = sum(_sss_scores(G, n2) for _, G in grids)
        o = np.argsort(acc)[::-1]
        med = np.median(acc)
        mad = 1.4826 * np.median(np.abs(acc - med)) + 1e-12
        cells.append(dict(nid2=n2, nid1=int(o[0]), z=float((acc[o[0]] - med) / mad), z2=float((acc[o[1]] - med) / mad)))
    cells.sort(key=lambda c: -c["z"])
    best_c = cells[0]
    nid1 = best_c["nid1"]
    info.update(nid1=nid1, nid2=best_c["nid2"], sss_z=best_c["z"], sss_z2=best_c["z2"],
                other_cells=[{"pci": 3 * c["nid1"] + c["nid2"], "sss_z": round(c["z"], 1)} for c in cells[1:]
                             if c["z"] >= 4.5 and c["z"] / max(c["z2"], 1e-3) >= 1.3])
    pci = 3 * nid1 + best_c["nid2"]
    per = []
    for p, G in grids:
        ib, llr = _pbch_llr(G, pci)
        per.append((p, ib, llr))
    info["ssb_indices"] = [ib for _, ib, _ in per]
    info["ibar_ssb"] = per[int(np.argmax([m[p] for p, _, _ in per]))][1]
    half = int(5e-3 * sl.ofs)
    groups = {}
    for p, ib, llr in per:
        groups.setdefault(p // half, []).append((ib, llr))
    tried = 0
    for lmax in (8, 4):
        for key in sorted(groups, key=lambda k: -len(groups[k])):
            soft = np.zeros(864)
            for ib, llr in groups[key]:
                c = _pbch_prbs(pci, ib % lmax, 864).astype(float)
                soft += llr * (1 - 2 * c)
            soft = soft / (np.std(soft) + 1e-12) * 3
            for sign in (1.0, -1.0):
                tried += 1
                try:
                    out = _bch_decode(sign * soft, groups[key][0][0] % 8, pci)
                except Exception:                                  # noqa: BLE001
                    continue
                if len(out) < 6:
                    continue
                _, crc, payload, sfn_lsb, hrf, msb = out
                if int(np.asarray(crc).ravel()[0]) == 0:
                    info.update(crc_ok=True, lmax=lmax, n_combined=len(groups[key]),
                                mib=_parse_mib(payload, sfn_lsb, hrf, msb))
                    return pci, info
    info.update(crc_ok=False, decode_attempts=tried)
    return pci, info


_G = [16, 23, 18, 17, 8, 30, 10, 6, 24, 7, 0, 5, 3, 2, 1, 4, 9, 11, 12, 13, 14, 15, 19, 20, 21, 22, 25, 26, 27, 28,
      29, 31]


def _bch_decode(softbits, lssb, ncellid):
    """BCH decoding for L_max = 8 (TS 38.212 7.1): rate recovery, polar decoding and CRC-24C
    (py3gpp), then the payload descrambling (7.1.2) and de-interleaving (7.1.1) with the
    cached PBCH scrambling sequence. Returns what py3gpp.nrBCHDecode(softbits, 8, lssb,
    ncellid) returns (checked equal in tests); py3gpp would regenerate a 3200-bit Gold
    sequence on every call."""
    from py3gpp.nrCRCDecode import nrCRCDecode
    from py3gpp.nrPolarDecode import nrPolarDecode
    from py3gpp.nrRateRecoverPolar import nrRateRecoverPolar
    K, E, L, N = 32 + 24, 864, 8, 512
    matched = nrRateRecoverPolar(softbits, K, N, False, False)
    decoded = nrPolarDecode(matched, K, E, L)
    scrblk, crc_result = nrCRCDecode(decoded, "24C")
    scrblk = scrblk[:, 0]
    v = 2 * scrblk[_G[7]] + scrblk[_G[8]]           # 3rd and 2nd LSB of the SFN (TS 38.212 7.1.2)
    A, M = 32, 32 - 3                               # L_max 8
    scrambling_seq = _pbch_prbs(ncellid, 0, len(scrblk) * 100)[v * M:][:A]
    final = np.zeros(A, "int")
    j = 0
    for i in range(A):
        if i in (_G[10], _G[8], _G[7]):
            final[i] = 0
        else:
            final[i] = scrambling_seq[j]
            j += 1
    a = np.bitwise_xor(final, scrblk)
    j_sfn, j_other = 0, 14                          # de-interleaving (TS 38.212 7.1.1)
    payload = np.empty(24, "int")
    for i in range(24):
        if 1 <= i < 7:
            payload[i] = a[_G[j_sfn]]
            j_sfn += 1
        else:
            payload[i] = a[_G[j_other]]
            j_other += 1
    lsbotfsfn = np.array([a[_G[j_sfn]], a[_G[j_sfn + 1]], a[_G[j_sfn + 2]], a[_G[j_sfn + 3]]])
    return scrblk, crc_result, payload, lsbotfsfn, a[_G[10]], 0


def _bits(b, i, n):
    return int("".join(str(int(x)) for x in b[i:i + n]), 2)


def _parse_mib(payload, sfn_lsb, hrf, msb):
    b = [int(x) for x in np.asarray(payload).ravel()]
    sfn_msb = _bits(b, 1, 6)
    lsb = _bits([int(x) for x in np.asarray(sfn_lsb).ravel()], 0, 4) if np.size(sfn_lsb) else 0
    k_ssb_lsb = _bits(b, 8, 4)
    k_msb = int(np.asarray(msb).ravel()[0]) if np.size(msb) else 0
    pdcch = _bits(b, 13, 8)
    return {"message": "BCCH-BCH mib" if b[0] == 0 else "messageClassExtension",
            "systemFrameNumber": (sfn_msb << 4) | lsb,
            "subCarrierSpacingCommon_khz": 15 if b[7] == 0 else 30,
            "ssb_SubcarrierOffset": k_ssb_lsb + 16 * k_msb,
            "dmrs_TypeA_Position": 2 if b[12] == 0 else 3,
            "pdcch_ConfigSIB1": {"controlResourceSetZero": pdcch >> 4, "searchSpaceZero": pdcch & 15},
            "cellBarred": "barred" if b[21] == 0 else "notBarred",
            "intraFreqReselection": "allowed" if b[22] == 0 else "notAllowed",
            "half_frame": int(np.asarray(hrf).ravel()[0]) if np.size(hrf) else None,
            "sib1_scheduled": k_ssb_lsb + 16 * k_msb <= 23}
