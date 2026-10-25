"""802.11b DSSS/CCK receiver (IEEE 802.11-2020 Clauses 15/16).

Per 22 MHz channel:
  detect   |Barker matched filter|^2 summed over 16 symbols at each of the 22 sample phases (2 samples/chip);
           DSSS is where the best phase stands out (peak / mean > 2.8, noise < 2.6), at any level; a run is
           split where its best phase moves (a second PPDU). No fixed windows: a PPDU is never cut.
  sync     4 samples/chip: early/late timing loop on the matched filter (fractional, interpolated);
           DBPSK with a sliding CFO estimate; descramble x^7+x^4+1; SFD (long 0xF3A0 / short 0x05CF,
           <= 1 bit error) after a SYNC run.
  train    the scrambled preamble is rebuilt exactly (best of the 128 scrambler states), which gives a
           data-aided CFO, phase and SNR; the PLCP header is demodulated coherently (phase-tracked), with
           a differential fallback, and checked with its CRC-16.
  PSDU     1/2 Mb/s: coherent Barker symbols on the tracked grid. 5.5/11 Mb/s CCK: chips on a grid whose
           rate follows the CFO (locked clocks), a 9-tap LS equaliser trained on the preamble + header
           chips, then ML codeword + phi1 decisions with phase tracking. Descramble, CRC-32 FCS, MAC parse.
  dedup    a PPDU also correlates on neighbouring channels; the strongest decode within +-3 channels is
           kept, and a beacon's DS parameter set names the channel.
"""
from __future__ import annotations

import zlib
from itertools import combinations

import numpy as np
from scipy.ndimage import uniform_filter
from scipy.signal import istft, oaconvolve, stft

from spectrum_intel.core.bands_india import WIFI_24

FS = 22e6                                   # detection rate: 2 samples per chip
FS_DEC = 44e6                               # decoding rate: 4 samples per chip
SPS = 2
SYM = 11 * SPS                              # detection samples per 1 us symbol
BARKER = np.array([1, -1, 1, 1, -1, 1, 1, 1, -1, -1, -1], float)
SFD_LONG = np.array([(0xF3A0 >> i) & 1 for i in range(16)], np.uint8)
SFD_SHORT = np.array([(0x05CF >> i) & 1 for i in range(16)], np.uint8)
RATES = {0x0A: 1.0, 0x14: 2.0, 0x37: 5.5, 0x6E: 11.0}
CCK_QPSK = {(0, 0): 0.0, (0, 1): np.pi / 2, (1, 0): np.pi, (1, 1): 3 * np.pi / 2}
DQPSK_BITS = [(0, 0), (0, 1), (1, 1), (1, 0)]          # phase change k*pi/2 -> dibit

DET_K = 16                                  # symbols averaged by the detector
DET_RATIO = 2.8                             # Barker energy peak / mean over sample phases (noise max 2.6)
DET_MIN_SYM = 32                            # flagged symbols for a candidate (a PPDU has >= 96 Barker symbols)
BARKER_REPORT_RHO = 0.2                     # coherence needed to report undecoded energy as Barker DSSS
EQ_TAPS = 9
FIX_CANDIDATES = 20                         # least reliable decisions offered to the FCS corrector
FIX_MAX_WEIGHT = 4                          # symbol corrections at most (~2e4 patterns: P(false FCS) ~ 5e-6)


def _descramble(bits):
    bits = np.asarray(bits, np.uint8)
    out = bits.copy()
    out[4:] ^= bits[:-4]
    out[7:] ^= bits[:-7]
    out[:7] = 0
    for i in range(min(7, len(bits))):                  # unknown history: the first 7 are meaningless
        out[i] = bits[i]
    return out


def _crc16_ok(bits32, crc_bits):
    crc = 0xFFFF
    for b in bits32:
        fb = ((crc >> 15) & 1) ^ int(b)
        crc = (crc << 1) & 0xFFFF
        if fb:
            crc ^= 0x1021
    crc ^= 0xFFFF
    return [(crc >> (15 - i)) & 1 for i in range(16)] == [int(b) for b in crc_bits]


def _bits_to_int(bits):
    return int(sum(int(b) << i for i, b in enumerate(bits)))


def _bytes(bits):
    bits = np.asarray(bits[:len(bits) // 8 * 8], np.uint8).reshape(-1, 8)
    return bytes((bits << np.arange(8, dtype=np.uint8)).sum(1).astype(np.uint8).tolist())


def _cck_codebook(n_bits):
    """Codewords with phi1 = 0, keyed by data bits d2.. (5.5: 2 bits, 11: 6 bits)."""
    book = []
    if n_bits == 6:
        for d in range(64):
            b = [(d >> i) & 1 for i in range(6)]
            p2, p3, p4 = CCK_QPSK[(b[0], b[1])], CCK_QPSK[(b[2], b[3])], CCK_QPSK[(b[4], b[5])]
            book.append((b, p2, p3, p4))
    else:
        for d in range(4):
            b = [(d >> i) & 1 for i in range(2)]
            p2, p3, p4 = b[0] * np.pi + np.pi / 2, 0.0, b[1] * np.pi
            book.append((b, p2, p3, p4))
    C = np.array([[np.exp(1j * (p2 + p3 + p4)), np.exp(1j * (p3 + p4)), np.exp(1j * (p2 + p4)),
                   -np.exp(1j * p4), np.exp(1j * (p2 + p3)), np.exp(1j * p3), -np.exp(1j * p2), 1.0]
                  for _, p2, p3, p4 in book])
    return C, np.array([b for b, *_ in book], np.uint8)


CCK11, CCK11_BITS = _cck_codebook(6)
CCK55, CCK55_BITS = _cck_codebook(2)


# ------------------------------------------------------------------------------------------- detection
def _candidates(y2):
    """[(sample0, sample1, peak coherence)] at the 2-sps detection rate.

    Per symbol m, P[m, phase] = |Barker MF|^2 summed over the last DET_K symbols. DSSS concentrates it at
    one sample phase: peak / mean over the 22 phases > DET_RATIO (noise: < 2.6 in 180 ms), whatever the
    level or the interference that adds to every phase. A run is split where its best phase moves and
    stays moved (a second PPDU on another timing), so each PPDU gets its own acquisition."""
    mf = np.repeat(BARKER[::-1], SPS)
    c = oaconvolve(y2, mf, mode="full")[len(mf) - 1:]
    nsym = len(c) // SYM
    if nsym < DET_MIN_SYM:
        return []
    P = (np.abs(c[:nsym * SYM]).astype(np.float64) ** 2).reshape(nsym, SYM)
    E = (np.abs(y2[:nsym * SYM]).astype(np.float64) ** 2).reshape(nsym, SYM).sum(1)
    Pk = np.cumsum(P, 0)
    Pk[DET_K:] -= Pk[:-DET_K].copy()
    Ek = np.cumsum(E)
    Ek[DET_K:] -= Ek[:-DET_K].copy()
    Ek[:-1] += E[1:]                                  # the late phases of symbol m reach into m + 1
    peak = Pk.max(1)
    ratio = peak / (Pk.mean(1) + 1e-30)
    rho = peak / (SYM * Ek + 1e-30)                   # Barker coherence: ~0.05 noise, ~0.6-0.8 clean PPDU
    best = Pk.argmax(1)
    flag = ratio > DET_RATIO
    idx = np.where(flag)[0]
    if not len(idx):
        return []
    runs, start, prev = [], idx[0], idx[0]
    for m in idx[1:]:
        if m - prev > 3:
            runs.append((start, prev))
            start = m
        prev = m
    runs.append((start, prev))
    out = []
    for a, b in runs:
        cuts, ph, moved = [a], best[a], 0
        for m in range(a + 1, b + 1):
            d = abs(int(best[m]) - int(ph)) % SYM
            if min(d, SYM - d) > 3:
                moved += 1
                if moved >= 6:                        # a new timing has taken over: new PPDU
                    cuts.append(m - 5)
                    ph, moved = best[m], 0
            else:
                moved = 0
        cuts.append(b + 1)
        for u, v in zip(cuts[:-1], cuts[1:]):
            if v - u >= DET_MIN_SYM // 2 or (u == a and v - u >= 8):
                out.append((max(0, (u - DET_K) * SYM), v * SYM, float(rho[u:v].max())))
    return out


# ------------------------------------------------------------------------------------------- tracking
def _interp(x, t):
    """Linear interpolation of a complex sequence at fractional indices t."""
    t = np.clip(t, 0, len(x) - 1.001)
    i = t.astype(int)
    f = t - i
    return x[i] * (1 - f) + x[i + 1] * f


class _Chan:
    """One 22 MHz channel at 4 samples/chip, with its Barker matched filter."""

    def __init__(self, y4, ofs4, fc):
        self.y = y4
        self.fc = fc
        self.fs = ofs4
        self.T = ofs4 / 1e6                                     # samples per 1 us symbol
        mf = np.repeat(BARKER[::-1], 4)
        self.c = oaconvolve(y4, mf, mode="full")[len(mf) - 1:]  # c[n]: Barker starting at n
        self.ac = np.abs(self.c)
        self.ym = oaconvolve(y4, np.ones(4) / 4, mode="full")[3:]  # chip MF, centred at n + 1.5

    def initial_timing(self, n0):
        """(tau, amplitude): best of the 44 sample phases over 16 symbols, starting at the first strong symbol."""
        L = int(round(self.T))
        seg = self.ac[n0:n0 + 16 * L]
        if len(seg) < 16 * L:
            return None, 0.0
        F = seg.reshape(16, L)
        ph = int(np.argmax((F ** 2).sum(0)))
        A = float(np.sort(F[:, ph])[-4:].mean())                    # peak amplitude of the strong symbols
        k = int(np.argmax(F[:, ph] > 0.4 * A))                       # first symbol of the PPDU
        tau = n0 + k * self.T + ph
        i = int(tau)
        if 1 <= i < len(self.ac) - 1:                                # parabolic refinement
            l, o, r = self.ac[i - 1], self.ac[i], self.ac[i + 1]
            den = l - 2 * o + r
            if den < 0:
                tau += float(np.clip(0.5 * (l - r) / den, -0.5, 0.5))
        return tau, A

    def track(self, tau, nsym, T=None, A=None):
        """Barker symbols from fractional sample tau: early/late loop (2nd order), updated only on symbols
        with at least 40% of the PPDU amplitude A, so noise and CCK chips cannot walk it. Returns (s, taus, T)."""
        T = self.T if T is None else T
        a = self.ac
        lim = len(self.c) - 6
        if A is None:
            i = int(tau)
            A = float(a[max(0, i - 2):i + 3].max()) if 0 <= i < len(a) else 0.0
        taus, s = [], []
        for _ in range(nsym):
            if tau < 4 or tau > lim:
                break
            i = int(tau)
            f = tau - i
            o = a[i] * (1 - f) + a[i + 1] * f
            if o < 0.6 * A:
                j = i - 3 + int(np.argmax(a[i - 3:i + 4]))
                if a[j] > 1.5 * o and a[j] > 0.5 * A:              # off the main lobe: re-acquire the peak
                    tau, i, f = float(j), j, 0.0
                    o = a[j]
            taus.append(tau)
            s.append(self.c[i] * (1 - f) + self.c[i + 1] * f)
            err = 0.0
            if o > 0.4 * A:
                e = a[i - 1] * (1 - f) + a[i] * f
                lt = a[i + 1] * (1 - f) + a[i + 2] * f
                err = max(-1.0, min(1.0, 3.0 * (lt - e) / (lt + e + 1e-12)))
                A = 0.98 * A + 0.02 * o
            T += 0.002 * err
            tau += T + 0.3 * err
        return np.array(s), np.array(taus), T


# ------------------------------------------------------------------------------------------- demod helpers
def _dbpsk_bits(s):
    """Differential DBPSK with a sliding (31-symbol) CFO estimate."""
    d = s[1:] * np.conj(s[:-1])
    d2 = d ** 2 / (np.abs(d) + 1e-30)
    w = np.ones(31)
    rot = np.angle(np.convolve(d2, w, "same")) / 2
    return (np.abs(np.angle(d * np.exp(-1j * rot))) > np.pi / 2).astype(np.uint8)


def _find_sfds(bits):
    """(bit index just after the SFD, preamble) with <= 1 SFD bit error after a SYNC run."""
    out = []
    if len(bits) < 40:
        return out
    W = np.lib.stride_tricks.sliding_window_view(bits, 16)
    ones = np.convolve(bits, np.ones(16, int), "valid")
    for sfd, pre in ((SFD_LONG, "long"), (SFD_SHORT, "short")):
        err = (W != sfd).sum(1)
        for k in np.where(err <= 1)[0]:
            if k < 16:
                continue
            run = ones[k - 16]
            if (pre == "long" and run >= 13) or (pre == "short" and run <= 3):
                if err[k] == 0 or (k + 1 < len(err) and err[k + 1] > 0 and err[k - 1] > 0):
                    out.append((int(k + 16), pre))
    return sorted(set(out))


def _rebuild_scrambled(raw, desc):
    """Scrambled bits s with s_n = desc_n ^ s_{n-4} ^ s_{n-7} that best match raw (all 128 start states)."""
    states = ((np.arange(128)[:, None] >> np.arange(7)) & 1).astype(np.uint8)   # [:, k] = s_{-1-k}
    hist = [states[:, k] for k in range(7)][::-1]                               # s_{-7} .. s_{-1}
    out = np.zeros((128, len(raw)), np.uint8)
    for n in range(len(raw)):
        s = desc[n] ^ hist[-4] ^ hist[-7]
        out[:, n] = s
        hist = hist[1:] + [s]
    agree = (out == raw[None, :]).sum(1)
    b = int(np.argmax(agree))
    return out[b], int(agree[b])


def _coherent(z, grid, ph0, q_prev, gain=0.08, resid=False):
    """Phase-tracked decisions on a k*grid constellation. z is already derotated to the reference.
    Returns (quadrant indices in units of pi/2, final phase[, residual v * exp(-j(phase + decision)))]."""
    step = int(round(grid / (np.pi / 2)))
    q = np.zeros(len(z), int)
    r = np.zeros(len(z), complex)
    ph = ph0
    for k, v in enumerate(z):
        a = np.angle(v * np.exp(-1j * ph))
        qq = int(np.round(a / grid)) * step
        q[k] = qq % 4
        r[k] = v * np.exp(-1j * (ph + qq * np.pi / 2))
        ph += gain * np.angle(r[k])
    return (q, ph, r) if resid else (q, ph)


def _dibit_flips(dq_old, dq_new, at):
    """Raw-bit positions that change when a DQPSK phase step changes from dq_old to dq_new."""
    a, b = DQPSK_BITS[dq_old % 4], DQPSK_BITS[dq_new % 4]
    return [at + i for i in range(2) if a[i] != b[i]]


def _qpsk_alternatives(q, q_prev, r, rate_bits):
    """(margin, raw-bit flips) for the nearest wrong decision of every coherent Barker symbol."""
    alts = []
    n = len(q)
    qq = np.concatenate([[q_prev], q])
    for k in range(n):
        v = r[k]
        if rate_bits == 1:                                              # DBPSK: q -> q + 2 flips bits k, k+1
            alts.append((abs(v.real), [k] + ([k + 1] if k + 1 < n else [])))
            continue
        for d, m in ((1, (v.real - v.imag) / np.sqrt(2)), (-1, (v.real + v.imag) / np.sqrt(2))):
            dq0, dq1 = qq[k + 1] - qq[k], (qq[k + 2] - qq[k + 1]) if k + 1 < n else None
            fl = _dibit_flips(dq0, dq0 + d, 2 * k)
            if dq1 is not None:
                fl += _dibit_flips(dq1, dq1 - d, 2 * k + 2)
            alts.append((abs(m), fl))
    return alts


def _fcs_syndrome(d):
    """CRC-32 of the MPDU body XOR its FCS field; 0 when the FCS passes. d: data bits (LSB first)."""
    b = _bytes(d)
    return zlib.crc32(b[:-4]) ^ int.from_bytes(b[-4:], "little"), b


def fcs_ok(mpdu: bytes) -> bool:
    return len(mpdu) > 4 and zlib.crc32(mpdu[:-4]) == int.from_bytes(mpdu[-4:], "little")


def _fcs_correct(d, alts):
    """Chase-style repair: find <= FIX_MAX_WEIGHT of the least reliable symbol decisions whose flips make
    the FCS pass. CRC-32 and the descrambler are linear, so each flip has a fixed syndrome; singles and
    pairs are matched against each other (meet in the middle). Returns (bits, n_flipped) or None."""
    nbits = len(d)
    S0, _ = _fcs_syndrome(d)
    if S0 == 0:
        return d, 0
    zero = zlib.crc32(bytes(nbits // 8 - 4))
    cand = sorted(alts, key=lambda a: a[0])[:FIX_CANDIDATES]
    errs, syns = [], []
    for _, flips in cand:
        e = np.zeros(nbits + 8, np.uint8)
        for p in flips:                                                 # raw flip -> data bits p, p+4, p+7
            for q in (p, p + 4, p + 7):
                e[q] ^= 1
        e = e[:nbits]
        b = _bytes(e)
        errs.append(e)
        syns.append(zlib.crc32(b[:-4]) ^ zero ^ int.from_bytes(b[-4:], "little"))
    groups = [(i,) for i in range(len(cand))] + list(combinations(range(len(cand)), 2))
    table = {}
    for g in groups:
        sg = 0
        for i in g:
            sg ^= syns[i]
        table.setdefault(sg, g)
    best = None
    for g in groups:
        sg = 0
        for i in g:
            sg ^= syns[i]
        for h in ([()] if sg == S0 else []) + ([table[S0 ^ sg]] if (S0 ^ sg) in table else []):
            combo = set(g) | set(h)
            if len(combo) == len(g) + len(h) and len(combo) <= FIX_MAX_WEIGHT:
                if best is None or len(combo) < len(best):
                    best = combo
    if best is None:
        return None
    out = d.copy()
    for i in best:
        out ^= errs[i]
    return (out, len(best)) if _fcs_syndrome(out)[0] == 0 else None


def _q_to_bits(q, q_prev, rate_bits):
    """Differential bits from absolute quadrant decisions (rate_bits 1: DBPSK, 2: DQPSK)."""
    dq = np.diff(np.concatenate([[q_prev], q])) % 4
    if rate_bits == 1:
        return (dq == 2).astype(np.uint8)
    return np.array([b for v in dq for b in DQPSK_BITS[v]], np.uint8)


def _diff_bits(s, rot, rate_bits):
    """Non-coherent differential detection (fallback)."""
    d = np.angle(s[1:] * np.conj(s[:-1]) * np.exp(-1j * rot))
    if rate_bits == 1:
        return (np.abs(d) > np.pi / 2).astype(np.uint8)
    q = np.round(d / (np.pi / 2)).astype(int) % 4
    return np.array([b for v in q for b in DQPSK_BITS[v]], np.uint8)


def psdu_bytes(rate, length_us, service):
    n = int(np.floor(length_us * rate / 8 + 1e-9))
    if rate == 11.0 and service & 0x80:
        n -= 1
    return n


# ------------------------------------------------------------------------------------------- PPDU decode
class _Ppdu:
    """Decode one PPDU whose SFD ends at bit index `e` of the tracked symbol run (s, taus)."""

    def __init__(self, ch: _Chan, s, taus, T, e, pre, amp):
        self.ch, self.s, self.taus, self.T, self.e, self.pre, self.amp = ch, s, taus, T, e, pre, amp

    def run(self):
        ch, s, e, pre = self.ch, self.s, self.e, self.pre
        npre = 144 if pre == "long" else 72
        p0 = max(0, e - npre + 1)                                        # first SYNC symbol
        desc = np.concatenate([np.full(npre - 16, 1 if pre == "long" else 0, np.uint8),
                               SFD_LONG if pre == "long" else SFD_SHORT])[npre - (e - p0):]
        raw = _dbpsk_bits(s[p0:e + 1])
        if len(raw) < e - p0 or e - p0 < 24:
            return None
        sc, agree = _rebuild_scrambled(raw[:e - p0], desc)
        # data-aided reference: symbol p0 has phase 0; symbol n+1 = symbol n + pi * scrambled bit n
        a = np.concatenate([[1.0], np.cumprod(1.0 - 2.0 * sc)])
        z = s[p0:e + 1] * a
        rot = float(np.angle(np.sum(z[1:] * np.conj(z[:-1]))))
        zr = z * np.exp(-1j * rot * np.arange(len(z)))
        mu = np.mean(zr)
        th0 = float(np.angle(mu))
        snr_sym = float(np.abs(mu) ** 2 / (np.var(zr) + 1e-30))
        self.rot, self.th0, self.p0 = rot, th0, p0
        self.q_pre = ((1 - a) * 1).astype(int) % 4                       # +1 -> 0, -1 -> 2
        rec = {"preamble": pre, "sync_bits_agree": round(agree / (e - p0), 3),
               "cfo_hz": round(rot / (2 * np.pi) * 1e6), "snr_db": round(10 * np.log10(snr_sym + 1e-12) - 10.4, 1),
               "start_s": (self.taus[e] - (npre - 1) * self.T) / ch.fs, "snr_db_symbol": round(10 * np.log10(snr_sym + 1e-12), 1)}
        # ---- PLCP header
        nh = 48 if pre == "long" else 24
        hs = s[e + 1:e + 1 + nh]
        if len(hs) < nh:
            return None
        k = np.arange(e + 1 - p0, e + 1 - p0 + nh)
        zh = hs * np.exp(-1j * (th0 + rot * k))
        hb_bits = 1 if pre == "long" else 2
        q_last = int(self.q_pre[-1])
        tries = []
        qh, ph = _coherent(zh, np.pi if hb_bits == 1 else np.pi / 2, 0.0, q_last)
        tries.append(("coherent", _q_to_bits(qh, q_last, hb_bits), qh, ph))
        tries.append(("differential", _diff_bits(np.concatenate([s[e:e + 1], hs]), rot, hb_bits), None, 0.0))
        for how, hraw, qh, ph in tries:
            hb = _descramble(np.concatenate([sc, hraw]))[len(sc):len(sc) + 48]
            if _crc16_ok(hb[:32], hb[32:48]):
                break
        else:
            hb = _descramble(np.concatenate([sc, tries[0][1]]))[len(sc):len(sc) + 48]
            rec.update(plcp_crc_ok=False, signal=_bits_to_int(hb[0:8]), rate_mbps=RATES.get(_bits_to_int(hb[0:8])),
                       service=_bits_to_int(hb[8:16]), length_us=_bits_to_int(hb[16:32]))
            return rec
        signal, service, length = _bits_to_int(hb[0:8]), _bits_to_int(hb[8:16]), _bits_to_int(hb[16:32])
        rate = RATES.get(signal)
        rec.update(plcp_crc_ok=True, signal=signal, rate_mbps=rate, service=service, length_us=length,
                   locked_clocks=bool(service & 0x04), header_demod=how)
        if rate is None:
            return rec
        if qh is None:                                                    # differential header: rebuild quadrants
            dq = (np.array([DQPSK_BITS.index((hraw[2 * i], hraw[2 * i + 1])) for i in range(24)])
                  if hb_bits == 2 else 2 * hraw.astype(int))
            qh = (q_last + np.cumsum(dq)) % 4
        self.q_hdr, self.ph = qh, ph
        self.prefix = np.concatenate([sc, hraw])
        nbytes = psdu_bytes(rate, length, service)
        rec["psdu_bytes"] = nbytes
        if nbytes < 10 or nbytes > 2400:
            return rec
        mpdu, how_p, nfix = self.psdu(rate, nbytes, e + 1 + nh, service)
        if mpdu is None:
            rec["truncated"] = True
            return rec
        from . import payload as WP
        ok = fcs_ok(mpdu)
        rec["fcs_ok"] = bool(ok)
        rec["psdu_demod"] = how_p
        rec["symbols_corrected"] = nfix
        if ok:
            rec["mac"] = WP.parse_mac(mpdu[:-4])
            rec["mpdu_hex"] = mpdu[:-4][:512].hex()
        return rec

    # ---- PSDU
    def psdu(self, rate, nbytes, k_psdu, service):
        """(mpdu with FCS, method, symbols corrected) -- the first method whose FCS passes, else a corrected
        decode of the methods that report decision reliabilities, else the first method's bytes."""
        nbits = nbytes * 8
        tried = []
        for how, fn in self._psdu_methods(rate):
            got = fn(rate, nbits, k_psdu)
            if got is None:
                continue
            bits, alts = got
            d = _descramble(np.concatenate([self.prefix, bits]))[len(self.prefix):len(self.prefix) + nbits]
            S, mpdu = _fcs_syndrome(d)
            if S == 0:
                return mpdu, how, 0
            tried.append((how, d, alts, mpdu))
        for how, d, alts, _ in tried:
            if alts:
                fixed = _fcs_correct(d, alts)
                if fixed is not None:
                    return _bytes(fixed[0]), how + f", FCS-guided correction of {fixed[1]} symbol(s)", fixed[1]
        return (tried[0][3], tried[0][0], 0) if tried else (None, None, 0)

    def _psdu_methods(self, rate):
        if rate in (1.0, 2.0):
            return [("coherent Barker", self._barker_coherent), ("differential Barker", self._barker_diff)]
        return [("CCK, LS equaliser", lambda r, n, k: self._cck(r, n, k, EQ_TAPS)),
                ("CCK, 1-tap", lambda r, n, k: self._cck(r, n, k, 1)),
                ("CCK, LS equaliser, fitted clock", lambda r, n, k: self._cck(r, n, k, EQ_TAPS, fit_clock=True))]

    def _barker_syms(self, rate, nbits, k_psdu):
        nsym = nbits // (1 if rate == 1.0 else 2)
        s, taus = self.s, self.taus
        have = len(s) - k_psdu
        if have < nsym:                                                   # continue tracking past the run
            more, _, _ = self.ch.track(taus[-1] + self.T, nsym - have + 1, self.T, self.amp)
            s = np.concatenate([s, more])
        seg = s[k_psdu:k_psdu + nsym]
        return seg if len(seg) == nsym else None

    def _barker_coherent(self, rate, nbits, k_psdu):
        seg = self._barker_syms(rate, nbits, k_psdu)
        if seg is None:
            return None
        k = np.arange(k_psdu - self.p0, k_psdu - self.p0 + len(seg))
        z = seg * np.exp(-1j * (self.th0 + self.rot * k))
        nb = 1 if rate == 1.0 else 2
        q, _, r = _coherent(z, np.pi if rate == 1.0 else np.pi / 2, self.ph, int(self.q_hdr[-1]), resid=True)
        return _q_to_bits(q, int(self.q_hdr[-1]), nb), _qpsk_alternatives(q, int(self.q_hdr[-1]), r, nb)

    def _barker_diff(self, rate, nbits, k_psdu):
        seg = self._barker_syms(rate, nbits, k_psdu)
        if seg is None:
            return None
        return _diff_bits(np.concatenate([self.s[k_psdu - 1:k_psdu], seg]), self.rot, 1 if rate == 1.0 else 2), []

    def _cck(self, rate, nbits, k_psdu, ntaps, fit_clock=False):
        ch = self.ch
        bps = 8 if rate == 11.0 else 4
        nsym = nbits // bps
        p0 = self.p0
        ks = np.arange(p0, k_psdu)
        taus = self.taus[p0:k_psdu]
        if len(taus) < 24:
            return None
        # symbol clock: locked clocks -> the chip rate follows the carrier offset
        T = ch.T / (1 + self.rot / (2 * np.pi) * 1e6 / ch.fc) if not fit_clock else np.polyfit(ks, taus, 1)[0]
        tau0 = float(np.mean(taus - T * (ks - p0)))                       # fitted time of symbol p0
        Tc = T / 11
        w_rad = self.rot / T                                              # carrier rotation per sample
        # training chips: preamble (+-1) and header (quadrants) symbols, known after the CRC
        qs = np.concatenate([self.q_pre[:len(ks) - len(self.q_hdr)], self.q_hdr])
        chips_known = (np.exp(1j * np.pi / 2 * qs)[:, None] * BARKER[None, :]).ravel()
        t_tr = tau0 + (np.arange(len(chips_known)) + 0.375) * Tc
        t0_psdu = tau0 + (k_psdu - p0) * T
        t_ps = t0_psdu + (np.arange(nsym * 8) + 0.375) * Tc
        if t_ps[-1] + 8 * Tc >= len(ch.ym):
            return None
        half = ntaps // 2

        def chips_at(t):
            tt = t[:, None] + Tc * np.arange(-half, half + 1)[None, :]
            return _interp(ch.ym, tt - 1.5) * np.exp(-1j * (w_rad * tt + self.th0 - w_rad * tau0))

        A = chips_at(t_tr)
        lam = 1e-3 * np.real(np.vdot(A[:, half], A[:, half]))
        w = np.linalg.solve(A.conj().T @ A + lam * np.eye(ntaps), A.conj().T @ chips_known)
        E = (chips_at(t_ps) @ w).reshape(nsym, 8)
        book, bb = (CCK11, CCK11_BITS) if rate == 11.0 else (CCK55, CCK55_BITS)
        corr = E @ np.conj(book).T / 8                                    # (nsym, ncode)
        out = np.zeros((nsym, bps), np.uint8)
        q0 = int(self.q_hdr[-1])
        q_prev, ph = q0, 0.0
        rot4 = np.exp(-1j * np.pi / 2 * np.arange(4))
        qs = np.zeros(nsym, int)
        second = []                                                       # (margin, codeword, quadrant)
        for i in range(nsym):
            v = corr[i] * np.exp(-1j * ph)
            M = np.real(v[:, None] * rot4[None, :])                       # metric of (codeword, phi1 quadrant)
            flat = np.argsort(M, axis=None)[-2:]
            m, q = np.unravel_index(flat[-1], M.shape)
            m2, q2 = np.unravel_index(flat[-2], M.shape)
            second.append((float(M[m, q] - M[m2, q2]), int(m2), int(q2)))
            a = np.angle(v[m] * np.exp(-1j * q * np.pi / 2))
            ph += 0.05 * a
            dq = (q - q_prev - (2 if i % 2 else 0)) % 4
            out[i, :2] = DQPSK_BITS[dq]
            out[i, 2:] = bb[m]
            qs[i] = q
            q_prev = q
        alts = []
        for i, (margin, m2, q2) in enumerate(second):                    # flips for the runner-up decision
            qp = qs[i - 1] if i else q0
            new = np.zeros(bps, np.uint8)
            new[:2] = DQPSK_BITS[(q2 - qp - (2 if i % 2 else 0)) % 4]
            new[2:] = bb[m2]
            fl = [i * bps + j for j in np.where(new != out[i])[0]]
            if i + 1 < nsym:
                nd = DQPSK_BITS[(qs[i + 1] - q2 - (2 if (i + 1) % 2 else 0)) % 4]
                fl += [(i + 1) * bps + j for j in range(2) if nd[j] != out[i + 1, j]]
            alts.append((margin, fl))
        return out.ravel(), alts


# ------------------------------------------------------------------------------------------- interference
def _excise(y, fs, ref, nper=128):
    """Time-frequency Wiener excision: cells of the STFT (2.9 us x 344 kHz) above the PPDU's own spectrum,
    measured on its clean preamble `ref` (sample range), are scaled down to it. DSSS survives losing part
    of its band; a narrowband or partial-band interferer (BLE, an adjacent-channel OFDM edge) does not."""
    nov = nper * 3 // 4
    _, _, X = stft(y, fs, nperseg=nper, noverlap=nov, return_onesided=False)
    P = np.abs(X) ** 2
    hop = nper - nov
    k0, k1 = ref[0] // hop + 1, max(ref[0] // hop + 5, ref[1] // hop)
    S = uniform_filter(np.median(P[:, k0:k1], axis=1), 5)
    W = np.minimum(1.0, S[:, None] / (uniform_filter(P, size=(3, 5)) + 1e-30))
    _, yc = istft(X * W, fs, nperseg=nper, noverlap=nov, input_onesided=False)
    out = y.copy()
    n = min(len(y), len(yc))
    out[:n] = yc[:n]
    return out


def _retry_excised(ch: _Chan, r):
    """Re-decode a PPDU whose header passed but FCS failed, after interference excision."""
    fs = ch.fs
    a = max(0, int((r["start_s"] - 60e-6) * fs))
    b = min(len(ch.y), int((r["end_s"] + 60e-6) * fs))
    pre_us = 140e-6 if r["preamble"] == "long" else 70e-6
    ref = (int((r["start_s"] - a / fs) * fs), int((r["start_s"] - a / fs + pre_us) * fs))
    seg = _Chan(_excise(ch.y[a:b], fs, ref), fs, ch.fc)
    for g in _decode_region(seg, 0, b - a, r["channel"], ch.fc):
        if g.get("fcs_ok") and abs(g["start_s"] + a / fs - r["start_s"]) < 5e-6:
            g["start_s"] += a / fs
            g["end_s"] += a / fs
            g["psdu_demod"] += ", after STFT interference excision"
            return g
    return None


# ------------------------------------------------------------------------------------------- scan
def _sfd_hits(ch: _Chan, n0, n1):
    """SFD-first acquisition: Barker symbols on a fixed 1 us grid at each of the ~44 sample phases (clock drift
    over a preamble is < 0.3 sample), DBPSK, descramble, SFD search. Hits of one PPDU at neighbouring phases
    are merged. Returns [(sample of the last SFD symbol, preamble, amplitude)], one per PPDU."""
    L = int(round(ch.T))
    nsym = int((n1 - n0) / ch.T) + 8
    k = np.arange(nsym)
    hits = []
    for ph in range(L):
        idx = (n0 + ph + k * ch.T).astype(int)
        idx = idx[idx < len(ch.c) - 1]
        if len(idx) < 60:
            break
        for e, pre in _find_sfds(_descramble(_dbpsk_bits(ch.c[idx]))):
            amp = float(np.median(ch.ac[idx[max(0, e - 48):e + 1]]))
            hits.append((int(idx[e]), pre, amp))
    hits.sort()
    out = []
    for h in hits:                                  # main lobe and sidelobes (+-0.7 us) of one SFD: keep the peak
        if out and h[0] - out[-1][0] <= 30 and h[1] == out[-1][1]:
            if h[2] > out[-1][2]:
                out[-1] = h
        else:
            out.append(h)
    return out


def _decode_region(ch: _Chan, n0, n1, chan, fc_ch, skip=lambda t0, amp: False):
    """Every PPDU whose SFD lies in [n0, n1): tracked from its own SYNC start, phase and amplitude."""
    out = []
    for n_e, pre, amp in _sfd_hits(ch, n0, n1):
        npre = 144 if pre == "long" else 72
        t_start = (n_e - (npre - 1) * ch.T) / ch.fs
        if amp <= 0 or skip(t_start, amp):
            continue
        tau = max(4.0, n_e - (npre - 1) * ch.T)
        i = int(tau)
        j = i - 3 + int(np.argmax(ch.ac[i - 3:i + 4])) if i >= 3 else i     # on the main lobe
        tau = float(j)
        e_guess = int(round((n_e - tau) / ch.T))
        s, taus, T = ch.track(tau, int((max(n1, n_e) - tau) / ch.T) + 80, A=amp)
        if len(s) < e_guess + 30:
            continue
        sf = [e for e, p in _find_sfds(_descramble(_dbpsk_bits(s))) if p == pre and abs(e - e_guess) <= 3]
        if not sf:
            continue
        r = _Ppdu(ch, s, taus, T, sf[0], pre, amp).run()
        if r is None:
            continue
        r.update(channel=chan, ch_freq_mhz=round(fc_ch / 1e6, 3), barker_amplitude=round(amp, 1))
        hdr_us = 192 if pre == "long" else 96
        r["end_s"] = r["start_s"] + hdr_us * 1e-6 + (r["length_us"] * 1e-6 if r["plcp_crc_ok"] else 0.0)
        out.append(r)
    return out


def scan(ctx):
    if "wifi_dsss" in ctx.cache:
        return ctx.cache["wifi_dsss"]
    cap = ctx.cap
    lo, hi = cap.fc - 0.47 * cap.fs, cap.fc + 0.47 * cap.fs
    cands = []
    for chan, f in WIFI_24.items():
        if f - 9e6 < lo or f + 9e6 > hi:
            continue
        y2, ofs2 = ctx.ch.slice(f, 11e6, FS)
        for a, b, rho in _candidates(y2):
            cands.append((rho, chan, f, a / ofs2, b / ofs2))
    cands.sort(key=lambda c: -c[0])                     # strongest first: the true channel wins dedup
    chans, pk = {}, []

    def dup(chan, t0):
        return [i for i, r in enumerate(pk) if r is not None and abs(r["channel"] - chan) <= 3
                and abs(r["start_s"] - t0) < 4e-6]

    def rank(r):
        return bool(r.get("fcs_ok")), bool(r["plcp_crc_ok"]), r["barker_amplitude"]

    for rho, chan, f, t0, t1 in cands:
        # a region inside the payload of a PPDU decoded (FCS-valid) on this channel is that PPDU: CSMA means
        # nothing else starts there that could be decoded. Neighbouring channels are always searched.
        if any(r is not None and r["plcp_crc_ok"] and r["channel"] == chan and r.get("fcs_ok")
               and t0 >= r["start_s"] + 200e-6 and t1 <= r["end_s"] + 10e-6 for r in pk):
            continue
        if chan not in chans:
            y4, ofs4 = ctx.ch.slice(f, 11e6, FS_DEC)
            chans[chan] = _Chan(y4, ofs4, f)
        ch = chans[chan]

        def skip(ts, amp, chan=chan):                   # a PPDU also correlates on neighbouring channels:
            return any(pk[i]["barker_amplitude"] * 1.2 >= amp for i in dup(chan, ts))   # the strongest wins

        for r in _decode_region(ch, int(t0 * ch.fs), int(t1 * ch.fs), chan, f, skip):
            r["coherence"] = round(rho, 3)
            r["start_s"] = float(r["start_s"])
            r["end_s"] = float(r["end_s"])
            same = dup(chan, r["start_s"])
            if same and max(rank(pk[i]) for i in same) >= rank(r):
                continue
            for i in same:
                pk[i] = None
            pk.append(r)
    pk = [r for r in pk if r is not None]
    for i, r in enumerate(pk):                         # last resort for header-ok / FCS-bad PPDUs
        if r["plcp_crc_ok"] and r.get("fcs_ok") is False and not r.get("truncated") and r["snr_db"] > -6:
            g = _retry_excised(chans[r.get("decoded_on_channel", r["channel"])], r)
            if g is not None:
                g["coherence"] = r["coherence"]
                pk[i] = g
    pk.sort(key=lambda r: r["start_s"])
    out = []
    for r in pk:
        ds = (r.get("mac") or {}).get("ds_channel")
        r["seen_on_channels"] = [r["channel"]]
        if ds and ds in WIFI_24:
            r["channel_source"] = "DS parameter set (beacon)"
            r["decoded_on_channel"] = r["channel"]
            r["channel"], r["ch_freq_mhz"] = ds, WIFI_24[ds] / 1e6
        else:
            r["channel_source"] = "strongest Barker coherence among channels within +-15 MHz"
        out.append(r)
    # Barker-coherent regions that no decoded PPDU explains (preamble lost, before the capture, or not 802.11)
    spans = [(r["decoded_on_channel"] if "decoded_on_channel" in r else r["channel"], r["start_s"], r["end_s"])
             for r in out if r["plcp_crc_ok"]]
    regions = []
    for rho, chan, f, t0, t1 in sorted(cands, key=lambda c: c[3]):
        if rho < BARKER_REPORT_RHO or any(abs(c - chan) <= 3 and t0 < e + 20e-6 and t1 > s0 - 20e-6
                                          for c, s0, e in spans):
            continue
        if any(abs(r["channel"] - chan) <= 3 and t0 < r["t1"] and t1 > r["t0"] for r in regions):
            r = next(r for r in regions if abs(r["channel"] - chan) <= 3 and t0 < r["t1"] and t1 > r["t0"])
            if rho > r["coherence"]:
                r.update(channel=chan, coherence=round(rho, 3), ch_freq_mhz=round(f / 1e6, 3))
            r["t0"], r["t1"] = min(r["t0"], t0), max(r["t1"], t1)
            continue
        regions.append({"t0": t0, "t1": t1, "channel": chan, "ch_freq_mhz": round(f / 1e6, 3),
                        "coherence": round(rho, 3)})
    ctx.cache["wifi_dsss"] = out
    ctx.cache["wifi_dsss_candidates"] = len(cands)
    ctx.cache["wifi_dsss_unexplained"] = regions
    return out
