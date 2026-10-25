#!/usr/bin/env python3
"""
wifi_detect.py -- detect and localize 802.11 OFDM PPDUs (a/g/n/ac/ax legacy
preamble) in a wideband complex IQ capture (e.g. SigMF, 61.44 MSps / 56 MHz).

Why a receiver, not an energy detector
--------------------------------------
Energy cannot tell Wi-Fi from other OFDM (e.g. drone video links) and cannot
split back-to-back frames (SIFS = 16 us). Every 802.11 OFDM PPDU starts with
the same legacy preamble, so we detect and DECODE that:

  L-STF 8 us : 0.8 us periodic        -> Schmidl-Cox detector, coarse CFO
  L-LTF 8 us : 2 known 3.2 us symbols -> fine timing, fine CFO, channel AND
                                         per-subcarrier noise (interference-aware)
  L-SIG 4 us : RATE + LENGTH           -> exact PPDU end; parity/reserved/tail/
                                         rate checks + re-encode bit agreement
  next syms  : QBPSK -> HT, BPSK+QBPSK -> VHT, L-SIG repeat -> HE;
               HT-SIG / VHT-SIG-A decoded and CRC-8 checked (MCS, bandwidth)

Robustness details that mattered on real data
  * soft bits weighted by per-subcarrier SNR measured from the L-LTF pair, and
    a reliability-weighted, outlier-rejecting pilot phase estimate -- a single
    narrowband interferer on one pilot otherwise rotates the whole L-SIG
  * free-ending Viterbi so the 6 tail bits are a real check
  * acceptance needs independent evidence: strong L-SIG bit agreement, or a
    passing SIG CRC, or moderate agreement + energy envelope ending where
    L-SIG says
  * partially covered edge channels decoded with out-of-band subcarriers erased
  * optional blind scan (1 MHz grid) finds Wi-Fi on NON-standard centres

Outputs per PPDU: start/end (capture samples and us), true centre frequency,
occupied band, format, L-SIG rate/length, MCS & bandwidth (when SIG CRC ok),
SNR, CFO, evidence. Separately: Wi-Fi activity that could not be decoded
(preamble found but L-SIG destroyed, e.g. collisions; or PPDU already in
progress at capture start), confirmed by preamble/OFDM-numerology evidence.

Usage
  python wifi_detect.py capture.sigmf-meta [--scan] [--out prefix]
  -> prefix.csv / .json (PPDUs), prefix_undecoded.csv, prefix.sigmf-meta
     (annotations), prefix.png (spectrogram overlay)
"""
from __future__ import annotations
import json, math, os
from dataclasses import dataclass, asdict, field
from fractions import Fraction
import numpy as np
from scipy.signal import resample_poly

FS_WIFI = 20e6

from spectrum_intel.core.emission import load_sigmf  # noqa: F401  (re-exported)

# ============================================================== constants
# L-LTF frequency-domain sequence, subcarriers -26..26
_LTF = np.array([1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 0,
                 1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1],
                float)
_STF = np.zeros(53, complex)
for k, v in zip([-24, -20, -16, -12, -8, -4, 4, 8, 12, 16, 20, 24],
                [1, -1, 1, -1, -1, 1, -1, -1, 1, 1, 1, 1]):
    _STF[k + 26] = v * (1 + 1j)
_STF *= np.sqrt(13 / 6)

SC = np.arange(-26, 27)
PILOT_SC = np.array([-21, -7, 7, 21])
PILOT_BASE = np.array([1, 1, 1, -1], float)
DATA_SC = np.array([k for k in SC if k != 0 and k not in PILOT_SC])      # 48 data subcarriers
USED_SC = np.array([k for k in SC if k != 0])                             # 52 used subcarriers

# HT (802.11n) 20 MHz uses 4 more subcarriers than legacy: -28..28, pilots in the same places.
HT_SC = np.arange(-28, 29)
HT_PILOT_SC = PILOT_SC
HT_DATA_SC = np.array([k for k in HT_SC if k != 0 and k not in HT_PILOT_SC])   # 52 data subcarriers
HT_USED_SC = np.array([k for k in HT_SC if k != 0])                            # 56 used subcarriers


def _bins(sc):
    return np.mod(sc, 64)


def _freq_to_time(sym_by_sc):
    X = np.zeros(64, complex)
    X[_bins(SC)] = sym_by_sc
    return np.fft.ifft(X) * 64 / np.sqrt(52)


LTF_TIME = _freq_to_time(_LTF)                  # one 3.2 us L-LTF symbol (64 samples)
STF_TIME = _freq_to_time(_STF)[:16]             # one 0.8 us L-STF period

# HT-LTF, 20 MHz (802.11-2020 Eq. 19-24): the legacy L-LTF sequence with {1,1} and {-1,-1} added on
# the outside, so it also trains the four subcarriers (+-27, +-28) that the legacy preamble never touches.
_HTLTF = np.concatenate([[1.0, 1.0], _LTF, [-1.0, -1.0]])        # subcarriers -28..28


def _pilot_polarity():
    """127-long pilot polarity sequence p_n (scrambler x^7+x^4+1, all-ones seed)."""
    s, out = [1] * 7, []
    for _ in range(127):
        b = s[6] ^ s[3]
        out.append(1 - 2 * b)
        s = [b] + s[:6]
    return np.array(out, float)


PILOT_POL = _pilot_polarity()

# RATE field (R1..R4) -> (Mbps, data bits per OFDM symbol)
RATES = {(1, 1, 0, 1): (6, 24), (1, 1, 1, 1): (9, 36), (0, 1, 0, 1): (12, 48), (0, 1, 1, 1): (18, 72),
         (1, 0, 0, 1): (24, 96), (1, 0, 1, 1): (36, 144), (0, 0, 0, 1): (48, 192), (0, 0, 1, 1): (54, 216)}


# ============================================================== coding
def _interleave_idx(ncbps=48):
    k = np.arange(ncbps)
    return (ncbps // 16) * (k % 16) + k // 16        # BPSK: second permutation is identity


_INTLV = _interleave_idx()

_G0, _G1 = 0o133, 0o171


def _parity(v):
    return bin(v).count("1") & 1


# Viterbi trellis: state = previous 6 input bits (newest at bit 5)
_OUT = np.zeros((64, 2, 2))
_NEXT = np.zeros((64, 2), int)
for s in range(64):
    for b in range(2):
        reg = (b << 6) | s
        _OUT[s, b] = [2 * _parity(reg & _G0) - 1, 2 * _parity(reg & _G1) - 1]
        _NEXT[s, b] = reg >> 1


def conv_encode(bits):
    s, out = 0, []
    for b in bits:
        reg = (b << 6) | s
        out += [_parity(reg & _G0), _parity(reg & _G1)]
        s = reg >> 1
    return np.array(out, int)


# State = the previous 6 input bits, so the trellis has exactly two predecessors per state and both
# reach it on the SAME input bit: dst = (b << 5) | (s >> 1), hence s in {2*(dst & 31), that + 1} and
# b = dst >> 5. Precomputing that turns each step into two 64-wide adds and a compare -- no sort, and
# the decoded bit needs no storage at all. Identical output to the general formulation.
_P0 = 2 * (np.arange(64) & 31)
_P1 = _P0 + 1
_BIT = np.arange(64) >> 5
_OUT0 = _OUT[_P0, _BIT]                            # (64, 2) branch symbols from predecessor 0
_OUT1 = _OUT[_P1, _BIT]


def viterbi(soft, terminated=False):
    """Soft Viterbi for the 802.11 K=7 r1/2 code. soft>0 means bit 1.
    Starts in state 0. terminated=False traces back from the best end state,
    so the 6 tail bits are genuinely decoded and "tail == 0" is a real check
    (tracing back from state 0 would force them to zero)."""
    soft = np.asarray(soft, float).reshape(-1, 2)
    n = len(soft)
    metric = np.full(64, -1e18)
    metric[0] = 0
    took = np.zeros((n, 64), np.uint8)              # which predecessor survived, per state per step
    for t in range(n):
        y = soft[t]
        m0 = metric[_P0] + _OUT0 @ y
        m1 = metric[_P1] + _OUT1 @ y
        win = m1 > m0
        metric = np.where(win, m1, m0)
        took[t] = win
    s = 0 if terminated else int(np.argmax(metric))
    bits = np.zeros(n, int)
    for t in range(n - 1, -1, -1):
        bits[t] = s >> 5                            # the input bit that led into state s
        s = 2 * (s & 31) + int(took[t, s])
    return bits


def crc8_sig(bits):
    """CRC-8 (x^8+x^2+x+1, init 0xFF, complemented) used by HT-SIG and VHT-SIG-A."""
    crc = 0xFF
    for b in bits:
        fb = ((crc >> 7) & 1) ^ int(b)
        crc = (crc << 1) & 0xFF
        if fb:
            crc ^= 0x07
    crc = ~crc & 0xFF
    return [(crc >> (7 - i)) & 1 for i in range(8)]


def _bits_to_int(bits):
    return int(sum(int(b) << i for i, b in enumerate(bits)))


# ============================================================== HE-SIG-A (802.11ax-2021 27.3.11.7)
# Two BPSK r1/2 symbols on 52 data subcarriers (-28..28, the HT set), each interleaved like 20 MHz HT
# (N_COL 13). The outer tones +-27 / +-28 are trained by L-SIG / RL-SIG, which carry -1 -1 -1 +1 there.
HE_EXT_SC = np.array([-28, -27, 27, 28])
HE_EXT_VAL = np.array([-1.0, -1.0, -1.0, 1.0])
_INTLV52 = 4 * (np.arange(52) % 13) + np.arange(52) // 13
HE_BW = {0: 20, 1: 40, 2: 80, 3: 160}
HE_GI_LTF = {0: "1x LTF + 0.8 us GI", 1: "2x LTF + 0.8 us GI", 2: "2x LTF + 1.6 us GI", 3: "4x LTF + 3.2 us GI"}


def he_txop_us(v):
    """TXOP field (7 bits): 127 = no duration; B0 = granularity (8 / 128 us), B1-B6 = the scaled value."""
    if v == 127:
        return None
    return 8 * (v >> 1) if not v & 1 else 512 + 128 * (v >> 1)


def parse_he_sig_a(bits, lsig_length):
    """52 decoded HE-SIG-A bits (A1 B0-B25, A2 B0-B25) -> fields, or None when CRC-4 / tail fail.
    L-SIG LENGTH mod 3 = 2 marks an HE MU PPDU; 1 an HE SU / ER SU / TB PPDU (B0 = 1 SU, 0 TB)."""
    b = [int(v) for v in bits]
    if crc8_sig(b[:42])[:4] != b[42:46] or any(b[46:52]):
        return None
    a1, a2 = b[:26], b[26:]
    g = lambda v, i, j: _bits_to_int(v[i:j])                    # noqa: E731
    out = {"txop_us": he_txop_us(g(a2, 0, 7))}
    if lsig_length % 3 == 2:
        bwi = g(a1, 15, 18)
        out.update(he_format="MU", uplink=bool(a1[0]), sigb_mcs=g(a1, 1, 4), sigb_dcm=a1[4],
                   bss_color=g(a1, 5, 11), spatial_reuse=g(a1, 11, 15),
                   bw=HE_BW.get(bwi, 80 if bwi in (4, 5) else 160), punctured=bwi >= 4,
                   n_sigb_symbols_or_users=g(a1, 18, 22), sigb_compression=a1[22], gi_ltf=HE_GI_LTF[g(a1, 23, 25)],
                   doppler=a1[25], n_he_ltf=g(a2, 8, 11), ldpc_extra=a2[11], stbc=a2[12])
    elif a1[0]:
        out.update(he_format="SU", beam_change=a1[1], uplink=bool(a1[2]), mcs=g(a1, 3, 7), dcm=a1[7],
                   bss_color=g(a1, 8, 14), spatial_reuse=g(a1, 15, 19), bw=HE_BW[g(a1, 19, 21)],
                   gi_ltf=HE_GI_LTF[g(a1, 21, 23)], nsts=g(a1, 23, 26) + 1, ldpc=a2[7], ldpc_extra=a2[8],
                   stbc=a2[9], beamformed=a2[10], doppler=a2[15])
    else:
        out.update(he_format="TB", uplink=True, bss_color=g(a1, 1, 7), bw=HE_BW[g(a1, 24, 26)])
    return out


# ============================================================== channels
def _all_channels():
    c = [(f"2g{ch}", ch, 2484e6 if ch == 14 else 2407e6 + 5e6 * ch) for ch in range(1, 15)]
    c += [(f"5g{ch}", ch, 5000e6 + 5e6 * ch) for ch in [*range(36, 65, 4), *range(100, 145, 4), *range(149, 178, 4)]]
    c += [(f"6g{ch}", ch, 5950e6 + 5e6 * ch) for ch in range(1, 234, 4)]
    return c


def channel_plan(fc_hz, fs_hz, usable_bw_hz, allow_partial=True, max_missing_hz=6e6):
    """20 MHz Wi-Fi channels inside the capture. With allow_partial, channels
    whose edge falls outside the usable band by up to max_missing_hz are kept;
    the missing subcarriers are erased and the K=7 code corrects around them."""
    cands = _all_channels()
    half = min(usable_bw_hz, fs_hz) / 2
    spill = max_missing_hz if allow_partial else 0.0
    return [dict(name=n, ch=c, freq_hz=f, offset_hz=f - fc_hz, usable_half_hz=half,
                 partial=abs(f - fc_hz) + 8.4e6 > half)
            for n, c, f in cands if abs(f - fc_hz) + 8.4e6 <= half + spill]


def scan_plan(fc_hz, fs_hz, usable_bw_hz, step_hz=1e6, max_missing_hz=6e6):
    """Blind scan: candidate 20 MHz centres every step_hz across the capture.
    The STF CFO estimator is unambiguous to +/-625 kHz, so a 1 MHz grid
    reaches a Wi-Fi signal at ANY centre frequency, standard channel or not."""
    std = channel_plan(fc_hz, fs_hz, usable_bw_hz, True, max_missing_hz)
    half = min(usable_bw_hz, fs_hz) / 2
    lim = half + max_missing_hz - 8.4e6
    out = []
    for off in np.arange(-np.floor(lim / step_hz), np.floor(lim / step_hz) + 1) * step_hz:
        f = fc_hz + off
        near = next((c for c in std if abs(c["freq_hz"] - f) < step_hz / 2), None)
        out.append(dict(name=near["name"] if near else f"{f/1e6:.1f}MHz", ch=near["ch"] if near else None,
                        freq_hz=f, offset_hz=off, usable_half_hz=half,
                        partial=abs(off) + 8.4e6 > half, standard=near is not None))
    return out


# Off-grid trigger: 802.11 links away from the standard channels (fixed-wireless backhaul at 4880 MHz, telemetry
# radios) are found by the 3.2 us CP correlation of 20 MHz slices every 5 MHz, with the PPDUs the standard scan
# decoded blanked; the blind 1 MHz scan then runs only around them.
OFFGRID_STEP_HZ = 5e6
OFFGRID_MAX_SLICE_HZ = 20e6                    # slice width: the largest power-of-two bin count within 20 MHz
OFFGRID_REF_LAGS_S = (2.1e-6, 2.3e-6, 2.5e-6, 2.7e-6, 3.7e-6, 3.9e-6, 4.1e-6, 4.3e-6)   # off the 0.8 us STF period
OFFGRID_RATIO = 4.0                            # CP lag correlation over the median of the reference lags
OFFGRID_MIN_CORR = 0.02
OFFGRID_MIN_RUN = 3                            # a 20 MHz signal lights >= 3 adjacent 5 MHz slices
OFFGRID_SPAN_HZ = 5e6                          # blind centres: run midpoint +- 5 MHz
# windowed fallback (PLAN G13): a weak, low-duty link (Mimosa 4880 MHz at ~6.6 dB, 20 % duty) is diluted in the
# whole-capture correlation; the same metric per 250 us window finds its PPDUs (real: ratio 6.5-12, noise << 6)
OFFGRID_WIN_S = 250e-6
OFFGRID_WIN_MIN_CORR = 0.06                    # per-window CP correlation (noise ~0.016 at 3840 samples)
OFFGRID_WIN_RATIO = 6.0
OFFGRID_WIN_MIN = 2                            # windows that pass, per slice
OFFGRID_WIN_GATE = 2.0                         # cheap trigger: whole-capture ratio (Mimosa slices 2.5-4.8, noise ~1)


def offgrid_plan(x, fs, fc, known, usable_bw=60e6, X=None):
    """Blind-scan channels (scan_plan entries, off the standard grid) around slices (every 5 MHz) whose 3.2 us
    cyclic-prefix correlation stands out once the known PPDUs' time spans are blanked. X: the block's FFT if the
    caller has it (core.dsp.Channelizer.X). The slice is a power-of-two number of bins (<= 20 MHz, 15.36 MHz for
    2^21 samples at 61.44 MSps): a fast inverse FFT, and wide enough for the CP test on a 20 MHz signal."""
    plan = [c for c in scan_plan(fc, fs, usable_bw) if not c["standard"]]
    n = len(x)
    df = fs / n
    M = 1 << int(np.floor(np.log2(OFFGRID_MAX_SLICE_HZ / df)))
    if not plan or M * df < 10e6:
        return []
    if X is None:
        X = np.fft.fft(np.asarray(x, np.complex64))
    half = M // 2
    ofs = M * df
    L = int(round(3.2e-6 * ofs))
    refs = [int(round(t * ofs)) for t in OFFGRID_REF_LAGS_S]
    lim = max(abs(c["offset_hz"]) for c in plan)
    offs = np.arange(-np.floor(lim / OFFGRID_STEP_HZ), np.floor(lim / OFFGRID_STEP_HZ) + 1) * OFFGRID_STEP_HZ
    hit = []
    for off in offs:
        idx = (int(round(off / df)) + np.arange(-half, half)) % n
        y = np.fft.ifft(np.fft.ifftshift(X[idx])).astype(np.complex128)    # correlations in double
        for p in known:
            if p.f_hi_mhz * 1e6 > fc + off - ofs / 2 and p.f_lo_mhz * 1e6 < fc + off + ofs / 2:
                y[int(p.start_sample / fs * ofs):int(p.end_sample / fs * ofs) + 1] = 0
        e = float(np.sum(np.abs(y) ** 2)) + 1e-30
        r = [abs(np.vdot(y[:-lag], y[lag:])) / e for lag in (L - 1, L, L + 1, *refs)]
        pk, ref = max(r[:3]), float(np.median(r[3:]))
        hit.append(pk >= OFFGRID_MIN_CORR and pk >= OFFGRID_RATIO * ref or
                   pk >= OFFGRID_WIN_GATE * ref and _offgrid_windows(y, ofs, L, refs))
    centres, k = [], 0
    while k < len(offs):
        if not hit[k]:
            k += 1
            continue
        j = k
        while j + 1 < len(offs) and hit[j + 1]:
            j += 1
        if j - k + 1 >= OFFGRID_MIN_RUN:
            centres.append(fc + 0.5 * (offs[k] + offs[j]))
        k = j + 1
    return [c for c in plan if any(abs(c["freq_hz"] - m) <= OFFGRID_SPAN_HZ + 1 for m in centres)]


def _offgrid_windows(y, ofs, L, refs):
    """Windowed CP test of one off-grid slice: >= OFFGRID_WIN_MIN windows of OFFGRID_WIN_S whose lag-3.2 us
    correlation passes OFFGRID_WIN_MIN_CORR and OFFGRID_WIN_RATIO x the median of the reference lags."""
    wn = int(OFFGRID_WIN_S * ofs)
    nw = (len(y) - max(L + 1, *refs)) // wn
    if nw < OFFGRID_WIN_MIN:
        return False
    m = nw * wn
    a = np.conj(y[:m])
    ew = (np.abs(y[:m]) ** 2).reshape(nw, wn).sum(1) + 1e-30
    rs = np.array([np.abs((a * y[lag:lag + m]).reshape(nw, wn).sum(1)) / ew for lag in (L - 1, L, L + 1, *refs)])
    pk, ref = rs[:3].max(0), np.median(rs[3:], 0)
    return int(np.count_nonzero((pk >= OFFGRID_WIN_MIN_CORR) & (pk >= OFFGRID_WIN_RATIO * ref))) >= OFFGRID_WIN_MIN


def merge_packets(known, extra, tol_us=2.0):
    """known + the extra (off-grid) PPDUs that are not the same PPDU again (start within tol_us, overlapping band);
    among duplicates in extra the highest STF metric is kept. The known PPDUs keep their order."""
    out, added = list(known), []
    for p in sorted(extra, key=lambda p: -p.stf_metric):
        if not any(abs(p.start_us - q.start_us) <= tol_us and p.f_lo_mhz < q.f_hi_mhz and q.f_lo_mhz < p.f_hi_mhz
                   for q in out + added):
            added.append(p)
    return out + sorted(added, key=lambda p: (p.start_sample, p.ch_freq_mhz))


def wide_center(ch_name, bw_mhz):
    """Center frequency (Hz) of the 40/80/160 MHz channel containing a 5 GHz primary."""
    if not ch_name.startswith("5g") or bw_mhz <= 20:
        return None
    ch = int(ch_name[2:])
    groups = {40: [(36, 40), (44, 48), (52, 56), (60, 64), (100, 104), (108, 112), (116, 120), (124, 128),
                   (132, 136), (140, 144), (149, 153), (157, 161), (165, 169), (173, 177)],
              80: [(36, 48), (52, 64), (100, 112), (116, 128), (132, 144), (149, 161), (165, 177)],
              160: [(36, 64), (100, 128), (149, 177)]}
    for lo, hi in groups.get(bw_mhz, []):
        if lo <= ch <= hi:
            return 5000e6 + 5e6 * (lo + hi) / 2
    return None


# ============================================================== records
@dataclass
class Packet:
    channel: str                   # standard channel name, or "<MHz>MHz" if off-grid
    ch_freq_mhz: float             # decode centre used
    center_mhz: float              # estimated true centre = decode centre + CFO
    start_sample: int              # at the capture sample rate
    end_sample: int
    start_us: float
    duration_us: float
    end_us: float
    fmt: str                       # NON_HT | HT | VHT | HE  ("HT?"/"VHT?": SIG CRC failed)
    lsig_rate_mbps: int
    lsig_length: int
    bandwidth_mhz: int
    f_lo_mhz: float                # occupied span (absolute RF)
    f_hi_mhz: float
    snr_db: float
    cfo_khz: float
    stf_metric: float
    mcs: int | None = None
    sig_crc_ok: bool | None = None
    guess: str = ""                # coarse frame-type hint (non-HT only)
    truncated: bool = False        # PPDU runs past end of capture
    energy_end_us: float | None = None
    extra: dict = field(default_factory=dict)


# ============================================================== receiver
class WifiDetector:
    def __init__(self, fs, fc, usable_bw=60e6, stf_thresh=0.75, min_plateau=50,
                 power_margin_db=4.0, channels=None, min_bit_agreement=0.85, strong_agreement=0.90,
                 end_tol_us=5.0, scan_step_hz=None):
        self.fs, self.fc = fs, fc
        self.stf_thresh, self.min_plateau = stf_thresh, min_plateau
        # fraction of received L-SIG coded bits consistent with the decoded word
        self.min_bit_agreement = min_bit_agreement
        self.rejects = []                           # preambles found but not decodable
        self.strong_agreement, self.end_tol_us = strong_agreement, end_tol_us
        self.power_margin = 10 ** (power_margin_db / 10)
        # usable_bw: subcarriers outside it are erased. 60 MHz (0.98 fs on the Pluto) keeps the weak edge
        # subcarriers, which the per-subcarrier SNR from the two L-LTF symbols down-weights instead; with
        # 56 MHz a channel at the capture edge lost ~16 of 52 subcarriers and its long PPDUs failed FCS.
        if channels is None:
            channels = (scan_plan(fc, fs, usable_bw, scan_step_hz) if scan_step_hz
                        else channel_plan(fc, fs, usable_bw))
        self.channels = channels
        r = Fraction(FS_WIFI / fs).limit_denominator(4096)
        self.up, self.down = r.numerator, r.denominator
        self.ratio = fs / FS_WIFI                      # capture samples per 20 MSps sample

    # ---------------------------------------------------------- per channel
    def baseband(self, x, offset_hz):
        n = np.arange(len(x))
        y = x * np.exp(-2j * np.pi * offset_hz / self.fs * n).astype(np.complex64)
        return resample_poly(y, self.up, self.down).astype(np.complex128)

    @staticmethod
    def stf_metric(r, L=48, D=16):
        c = r[:-D] * np.conj(r[D:])
        e = np.abs(r[D:]) ** 2
        P = np.convolve(c, np.ones(L), "valid")
        R = np.convolve(e, np.ones(L), "valid")
        return P, R

    def detect_channel(self, x20, chan):
        P, R = self.stf_metric(x20)
        M = np.abs(P) / (R + 1e-20)
        pwr = R / 48
        floor = np.percentile(pwr, 10)
        ok = (M > self.stf_thresh) & (pwr > floor * self.power_margin)
        d = np.diff(np.concatenate([[0], ok.astype(np.int8), [0]]))
        starts, ends = np.where(d == 1)[0], np.where(d == -1)[0]

        pkts, busy_until = [], -1
        for s0, e0 in zip(starts, ends):
            if e0 - s0 < self.min_plateau or s0 < busy_until:
                continue
            pk = self._decode_at(x20, s0, e0, P, M, floor, chan)
            if pk is not None:
                pkts.append(pk)
                busy_until = pk["_end20"] - 8            # allow a new PPDU right after this one
        return pkts

    # ---------------------------------------------------------- preamble decode
    def _decode_at(self, r, s0, e0, P, M, noise_pwr, chan):
        n = len(r)
        cfo_c = -np.angle(np.mean(P[s0:e0])) / (2 * np.pi * 16) * FS_WIFI    # coarse CFO
        self._cur = dict(freq_hz=chan["freq_hz"], cfo=cfo_c, stf=float(M[s0:e0].mean()),
                         plateau=int(e0 - s0), s0=int(s0))
        # L-LTF fine timing: search where T1 should be (~ STF start + 192)
        lo, hi = s0 + 100, min(s0 + 330, n - 700)
        if hi <= lo:
            return None
        seg = r[lo:hi + 128] * np.exp(-2j * np.pi * cfo_c / FS_WIFI * np.arange(lo, hi + 128))
        xc = np.abs(np.correlate(seg, LTF_TIME, "valid"))
        if len(xc) < 65:
            return None
        pair = xc[:-64] + xc[64:]
        t1 = lo + int(np.argmax(pair))
        start = t1 - 192
        if start < 0:
            return None
        # fine CFO from the two LTF symbols (64 apart)
        dph = cfo_c
        a = r[t1:t1 + 64] * np.exp(-2j * np.pi * cfo_c / FS_WIFI * np.arange(t1, t1 + 64))
        b = r[t1 + 64:t1 + 128] * np.exp(-2j * np.pi * cfo_c / FS_WIFI * np.arange(t1 + 64, t1 + 128))
        cfo = cfo_c - np.angle(np.sum(b * np.conj(a))) / (2 * np.pi * 64) * FS_WIFI

        if start + 560 + 16 > n:
            return None
        rr = r[start:start + 800] * np.exp(-2j * np.pi * cfo / FS_WIFI * np.arange(start, start + min(800, n - start)))
        back = 3                                             # FFT window a few samples into the CP

        def fft_at(off):
            w = rr[off - back:off - back + 64]
            return np.fft.fft(w) / 64 * np.sqrt(52)

        Y1, Y2 = fft_at(192), fft_at(256)
        L = np.zeros(64); L[_bins(SC)] = _LTF
        sc_f = chan["offset_hz"] + USED_SC * 312.5e3
        good_sc = USED_SC[np.abs(sc_f) <= chan.get("usable_half_hz", 1e12) - 0.4e6]
        used = _bins(good_sc)
        H = np.full(64, np.inf, complex); H[used] = (Y1[used] + Y2[used]) / 2 / L[used]
        # A real L-LTF spreads energy over all 52 subcarriers; a CW tone (which also
        # fools the STF detector -- it is periodic at every lag) puts it in 1-2 bins.
        pw = np.abs((Y1[used] + Y2[used]) / 2) ** 2 + 1e-20
        ltf_flat = float(np.exp(np.mean(np.log(pw))) / np.mean(pw))
        self._cur["ltf_flatness"] = ltf_flat
        # Per-subcarrier noise from the two identical L-LTF symbols. Narrowband
        # interference (a CW spur, a neighbour's skirt) shows up here, so it is
        # measured -- not assumed away. Light smoothing, floored for stability.
        nk = np.zeros(64); nk[used] = np.abs(Y1[used] - Y2[used]) ** 2 / 2
        order = np.argsort(good_sc)
        nsm = np.convolve(nk[used][order], np.ones(3) / 3, "same")
        nk_s = np.zeros(64); nk_s[used[order]] = np.maximum(nsm, np.median(nsm) * 0.5)
        noise = float(np.median(nk[used]))
        sig = np.mean(np.abs(H[used]) ** 2)
        snr_db = 10 * np.log10(max(sig / max(noise, 1e-20), 1e-3))
        W_rel = np.zeros(64); W_rel[used] = np.abs(H[used]) ** 2 / nk_s[used]     # per-subcarrier SNR
        W_rel[used] = np.minimum(W_rel[used], np.median(W_rel[used]) * 4)          # cap over-confident bins

        def equalize(off, n_sym):
            Y = fft_at(off)
            Z = np.zeros(64, complex); Z[used] = Y[used] / H[used]   # erased bins stay 0
            pb = _bins(PILOT_SC)
            pil = Z[pb] / (PILOT_BASE * PILOT_POL[n_sym % 127])
            w = W_rel[pb].copy()
            live = w > 0
            if live.sum() >= 3:                                      # drop the worst-fitting pilot
                ph0 = np.angle(np.sum(w * pil))
                err = np.abs(pil * np.exp(-1j * ph0) - 1)
                err[~live] = -1
                w[np.argmax(err)] = 0 if err.max() > 0.8 else w[np.argmax(err)]
            ph = np.angle(np.sum(w * pil))                           # reliability-weighted CPE
            return Z * np.exp(-1j * ph)

        def data_soft(Z, qbpsk=False):
            d = Z[_bins(DATA_SC)]
            v = (d.imag if qbpsk else d.real) * W_rel[_bins(DATA_SC)]   # LLR ~ SNR_k * y
            return v[_INTLV]                                 # deinterleave: coded bit k was on subcarrier i(k)

        def he_sig_soft():
            """104 soft bits of HE-SIG-A1 / A2 (FFT at 496 / 576, pilot polarity p_2 / p_3)."""
            pb, ext = _bins(PILOT_SC), _bins(HE_EXT_SC)

            def cpe(Y, n_sym):
                pil = Y[pb] / H[pb] / (PILOT_BASE * PILOT_POL[n_sym])
                return np.angle(np.sum(W_rel[pb] * np.nan_to_num(pil)))

            Hf, Wf = H.copy(), W_rel.copy()
            Hf[ext] = sum(fft_at(o)[ext] * np.exp(-1j * cpe(fft_at(o), k)) / HE_EXT_VAL
                          for o, k in ((336, 0), (416, 1))) / 2
            edge_ok = np.isin(_bins(np.array([-26, 26])), used).all()
            Wf[ext] = np.median(W_rel[used]) if edge_ok else 0.0
            db = _bins(HT_DATA_SC)
            soft = []
            for o, k in ((496, 2), (576, 3)):
                Y = fft_at(o)
                Z = Y[db] / Hf[db] * np.exp(-1j * cpe(Y, k))
                soft.append((np.nan_to_num(Z.real) * Wf[db])[_INTLV52])
            return np.concatenate(soft)

        # ---- L-SIG
        Zs = equalize(336, 0)
        soft0 = data_soft(Zs)
        lbits = viterbi(soft0)
        rate_bits = tuple(lbits[0:4])
        reserved, length, parity, tail = lbits[4], _bits_to_int(lbits[5:17]), lbits[17], lbits[18:24]
        why = [n for n, bad in (("rate", rate_bits not in RATES), ("reserved", reserved != 0),
                                ("parity", (lbits[:17].sum() & 1) != parity), ("tail", tail.any()),
                                ("zero_length", length == 0)) if bad]
        if why:
            self.rejects.append(dict(self._cur, start20=start, reason="lsig:" + "+".join(why)))
            return None
        dl = Zs[_bins(DATA_SC)][np.isin(_bins(DATA_SC), used)]
        lsig_evm_db = 10 * np.log10(np.mean((dl.real - np.sign(dl.real)) ** 2 + dl.imag ** 2) + 1e-12)
        live = soft0 != 0                                        # erased subcarriers carry 0
        agree = float(np.mean((2 * conv_encode(lbits) - 1)[live] * np.sign(soft0[live]) > 0))
        if agree < self.min_bit_agreement:
            self.rejects.append(dict(self._cur, start20=start, reason=f"lsig_bit_agreement={agree:.2f}"))
            return None
        mbps, ndbps = RATES[rate_bits]
        dur_us = 20 + 4 * math.ceil((16 + 8 * length + 6) / ndbps)

        # ---- format classification (only 6 Mbps L-SIG can front HT/VHT/HE)
        fmt, mcs, bw, crc_ok, extra = "NON_HT", None, 20, None, {}
        if mbps == 6 and start + 560 + 16 <= n:
            Z1, Z2 = equalize(416, 1), equalize(496, 2)
            live = np.isin(_bins(DATA_SC), used)
            d0, d1, d2 = (Z[_bins(DATA_SC)][live] for Z in (Zs, Z1, Z2))
            q1 = np.mean(d1.imag ** 2) / (np.mean(d1.real ** 2) + 1e-12)
            q2 = np.mean(d2.imag ** 2) / (np.mean(d2.real ** 2) + 1e-12)
            rl = np.abs(np.vdot(d0, d1)) / (np.linalg.norm(d0) * np.linalg.norm(d1) + 1e-12)
            if rl > 0.85 and (length % 3) != 0:
                fmt = "HE"
                extra["rl_sig_corr"] = round(float(rl), 3)
                if start + 640 + 16 <= n:
                    he = parse_he_sig_a(viterbi(he_sig_soft()), length)
                    crc_ok = he is not None
                    if he:
                        mcs, bw = he.get("mcs"), he["bw"]
                        extra.update(he)
            elif q1 > 2.0 and length % 3 == 0:
                bits = viterbi(np.concatenate([data_soft(Z1, True), data_soft(Z2, True)]))
                crc_ok = crc8_sig(bits[:34]) == list(bits[34:42]) and not bits[42:48].any()
                fmt = "HT" if crc_ok else "HT?"
                if crc_ok:
                    mcs = _bits_to_int(bits[0:7]); bw = 40 if bits[7] else 20
                    extra.update(ht_length=_bits_to_int(bits[8:24]), aggregation=int(bits[27]),
                                 smoothing=int(bits[24]),
                                 short_gi=int(bits[31]), ldpc=int(bits[30]),
                                 stbc=_bits_to_int(bits[28:30]), n_ess=_bits_to_int(bits[32:34]),
                                 nsts_su=mcs // 8 + 1)
            elif q2 > 2.0 and q1 < 0.5 and length % 3 == 0:
                bits = viterbi(np.concatenate([data_soft(Z1, False), data_soft(Z2, True)]))
                crc_ok = crc8_sig(bits[:34]) == list(bits[34:42]) and not bits[42:48].any()
                fmt = "VHT" if crc_ok else "VHT?"
                if crc_ok:
                    bw = [20, 40, 80, 160][_bits_to_int(bits[0:2])]
                    mcs = _bits_to_int(bits[28:32])
                    gid = _bits_to_int(bits[4:10])
                    extra.update(nsts_su=_bits_to_int(bits[10:13]) + 1, short_gi=int(bits[24]),
                                 ldpc=int(bits[26]), stbc=int(bits[3]), group_id=gid,
                                 beamformed=int(bits[32]), txop_ps_not_allowed=int(bits[22]))
                    # Group ID 0 / 63 is single-user; then B13-B21 are the Partial AID, which names
                    # the receiver even when the payload (80 MHz, LDPC, MIMO) cannot be decoded.
                    if gid in (0, 63):
                        extra["partial_aid"] = _bits_to_int(bits[13:22])
                    else:
                        extra["mu_mimo"] = True                     # 1-62: a multi-user group

        guess = ""
        if fmt == "NON_HT":
            guess = {14: "ACK/CTS", 20: "RTS", 32: "BlockAck"}.get(length, "")
            if not guess:
                guess = "control?" if length < 40 else "data/mgmt"

        end20 = start + int(round(dur_us * FS_WIFI / 1e6))
        ltf_q = float(np.abs(np.vdot(Y1[used], Y2[used])) / (np.linalg.norm(Y1[used]) * np.linalg.norm(Y2[used]) + 1e-12))
        extra["ltf_coherence"] = round(ltf_q, 3)
        extra["ltf_flatness"] = round(ltf_flat, 3)
        extra["lsig_evm_db"] = round(float(lsig_evm_db), 1)
        extra["lsig_bit_agreement"] = round(agree, 3)
        extra["erased_subcarriers"] = int(52 - len(used))
        rx = dict(start=start, cfo=cfo, back=back, used=used, H=H, W_rel=W_rel,
                  good_sc=good_sc, r=r, nk=nk_s, noise=noise)  # state a payload decoder needs
        return dict(_start20=start, _end20=end20, dur_us=dur_us, fmt=fmt, mbps=mbps, length=length,
                    bw=bw, mcs=mcs, crc_ok=crc_ok, snr=snr_db, cfo=cfo, stf=float(M[s0:e0].mean()),
                    guess=guess, extra=extra, noise_pwr=noise_pwr, _rx=rx)

    # ---------------------------------------------------------- energy check
    def _energy_end(self, x20, start20, end20, noise_pwr, gap_us=4.0):
        """Independent check of where the burst stops: the first >= gap_us stretch
        after the preamble where smoothed power stays near the noise floor.
        (A SIFS-separated ACK/BlockAck 16 us later is then NOT counted as part of it.)"""
        k = 20                                              # 1 us smoothing
        stop = min(len(x20), end20 + 800)
        a = start20 + 400                                   # past the 20 us legacy preamble
        if stop - a < 2 * k:
            return None
        p = np.convolve(np.abs(x20[a:stop]) ** 2, np.ones(k) / k, "same")
        low = (p < noise_pwr * 4).astype(np.int8)
        g = int(gap_us * 20)
        run = np.convolve(low, np.ones(g, int), "valid")
        idx = np.where(run == g)[0]
        return None if idx.size == 0 else (a + idx[0]) / FS_WIFI * 1e6

    # ---------------------------------------------------------- driver
    def run(self, x):
        out = []
        n_cap = len(x)
        for chan in self.channels:
            x20 = self.baseband(x, chan["offset_hz"])
            for d in self.detect_channel(x20, chan):
                s = int(round(d["_start20"] * self.ratio))
                e = int(round(d["_end20"] * self.ratio))
                trunc = e > n_cap
                bw = d["bw"]
                cen = wide_center(chan["name"], bw) or (chan["freq_hz"] + d["cfo"])
                half = {20: 8.3e6, 40: 18.1e6, 80: 38.1e6, 160: 78.1e6}[bw]
                ch_lo, ch_hi = self.fc - self.fs / 2, self.fc + self.fs / 2
                out.append(Packet(
                    channel=chan["name"], ch_freq_mhz=chan["freq_hz"] / 1e6,
                    center_mhz=round((chan["freq_hz"] + d["cfo"]) / 1e6, 3),
                    start_sample=s, end_sample=min(e, n_cap),
                    start_us=round(s / self.fs * 1e6, 3), duration_us=float(d["dur_us"]),
                    end_us=round(e / self.fs * 1e6, 3), fmt=d["fmt"],
                    lsig_rate_mbps=d["mbps"], lsig_length=d["length"], bandwidth_mhz=bw,
                    f_lo_mhz=round(max(cen - half, ch_lo) / 1e6, 3),
                    f_hi_mhz=round(min(cen + half, ch_hi) / 1e6, 3),
                    snr_db=round(float(d["snr"]), 1), cfo_khz=round(float(d["cfo"]) / 1e3, 2),
                    stf_metric=round(d["stf"], 3), mcs=d["mcs"], sig_crc_ok=d["crc_ok"],
                    guess=d["guess"], truncated=trunc,
                    energy_end_us=(lambda v: None if v is None else round(v, 2))(
                        self._energy_end(x20, d["_start20"], d["_end20"], d["noise_pwr"])),
                    extra=d["extra"]))
                # spectrum_intel: keep the receiver state so wifi_payload can decode the MPDU
                # from this same pass (receivers.wifi.plugin_ofdm pops it before anything is serialised)
                out[-1].extra["_rx"] = (d.get("_rx"), chan["offset_hz"], chan.get("usable_half_hz", 1e12))
        acc = []
        for p in out:
            if self._accept(p):
                acc.append(p)
            else:
                self.rejects.append(dict(freq_hz=p.ch_freq_mhz * 1e6, cfo=p.cfo_khz * 1e3, stf=p.stf_metric,
                                         plateau=0, s0=0, start20=int(round(p.start_us * 20)),
                                         reason="weak_lsig_evidence"))
        out = acc
        out = self._merge_grid(out)
        out = self._dedupe(out)
        out.sort(key=lambda p: (p.start_sample, p.ch_freq_mhz))
        return out

    def _accept(self, p):
        """Independent lines of evidence that an L-SIG decode is real:
          * the received coded bits agree strongly with the decoded word, or
          * the HT-SIG / VHT-SIG-A CRC passed, or
          * moderate agreement AND the burst's measured energy ends where
            L-SIG says it should (within end_tol_us)."""
        ag = p.extra["lsig_bit_agreement"]
        end_ok = (p.truncated or (p.energy_end_us is not None
                                  and abs(p.energy_end_us - p.end_us) <= self.end_tol_us))
        p.extra["energy_end_consistent"] = bool(end_ok)
        if ag >= self.strong_agreement or p.sig_crc_ok is True:
            p.extra["evidence"] = "sig_crc" if p.sig_crc_ok else "lsig_strong"
            return True
        if ag >= self.min_bit_agreement and end_ok:
            p.extra["evidence"] = "lsig+energy"
            return True
        return False

    @staticmethod
    def _merge_grid(pkts, tol_us=1.0, tol_mhz=2.5):
        """The same PPDU decoded from neighbouring scan centres: keep the decode
        whose centre is closest to the emitter (smallest |CFO|); re-label it with
        the standard channel name when the true centre lands on one."""
        keep = []
        for p in sorted(pkts, key=lambda p: abs(p.cfo_khz)):
            if any(abs(k.start_us - p.start_us) < tol_us and abs(k.center_mhz - p.center_mhz) < tol_mhz
                   for k in keep):
                continue
            keep.append(p)
        std = [(n, f) for n, _, f in _all_channels()]
        for p in keep:
            name, f = min(std, key=lambda nf: abs(nf[1] / 1e6 - p.center_mhz))
            p.extra["off_grid"] = bool(abs(f / 1e6 - p.center_mhz) > 1.0)
            if not p.extra["off_grid"]:
                p.channel = name
        return keep

    @staticmethod
    def _dedupe(pkts, tol_us=1.0):
        """A 40/80/160 MHz PPDU repeats its legacy preamble in every 20 MHz
        sub-channel, so it is found once per sub-channel. Detections that start
        at the same instant on channels >= 20 MHz apart are one PPDU: keep the
        cleanest decode (SIG CRC pass, then lowest L-SIG EVM). If they
        disagree on L-SIG content the weaker one was a mis-decode."""
        def score(p):
            return (p.sig_crc_ok is True, -p.extra.get("lsig_evm_db", 0), p.snr_db)
        keep = []
        for p in sorted(pkts, key=score, reverse=True):
            dup = next((k for k in keep if abs(k.start_us - p.start_us) < tol_us
                        and abs(k.ch_freq_mhz - p.ch_freq_mhz) >= 20), None)
            if dup is None:
                keep.append(p); continue
            agree = k_agree = (dup.lsig_length == p.lsig_length and dup.lsig_rate_mbps == p.lsig_rate_mbps)
            dup.extra.setdefault("also_seen_on", []).append(p.channel + ("" if agree else "(lsig mismatch)"))
            if agree and dup.bandwidth_mhz == 20 and dup.sig_crc_ok is None:
                dup.extra["bw_note"] = "legacy preamble duplicated on >1 sub-channel (>=40 MHz or non-HT dup)"
        return keep


# ============================================================== undecoded Wi-Fi
@dataclass
class Burst:
    """Wi-Fi activity that could not be decoded into a PPDU record."""
    label: str                 # wifi_preamble_undecodable | wifi_truncated_at_capture_start
    center_mhz: float
    start_sample: int
    end_sample: int
    start_us: float
    duration_us: float
    f_lo_mhz: float
    f_hi_mhz: float
    stf_metric: float | None
    cp64_score: float          # normalised CP correlation at the 3.2 us Wi-Fi symbol lag
    cp_other_score: float      # same at non-Wi-Fi lags (baseline)
    reason: str = ""
    seen_on: int = 1           # number of scan centres that saw this preamble


def cp_score(seg, lag, L=8):
    """95th percentile of |sum x[n] x*[n+lag]| / sqrt(E1 E2) over sliding windows
    of L samples. Bounded to [0,1]. Wi-Fi OFDM (64-sample symbol + 16/8 GI at
    20 MSps) peaks at lag 64 when the window sits on a cyclic prefix."""
    if len(seg) < lag + 4 * L:
        return 0.0
    c = seg[:-lag] * np.conj(seg[lag:])
    e1 = np.abs(seg[:-lag]) ** 2; e2 = np.abs(seg[lag:]) ** 2
    w = np.ones(L)
    P = np.abs(np.convolve(c, w, "valid"))
    R = np.sqrt(np.convolve(e1, w, "valid") * np.convolve(e2, w, "valid"))
    return float(np.percentile(P / (R + 1e-20), 95))


def numerology(seg):
    s64 = cp_score(seg, 64)
    other = float(np.median([cp_score(seg, l) for l in (37, 51, 77, 93, 111)]))
    return s64, other


def _extent(x20, start20, noise_pwr, gap_us=4.0, max_us=5500.0):
    """From start20, where does the burst stop (first >= gap_us near the floor)?"""
    k, g = 20, int(gap_us * 20)
    stop = min(len(x20), start20 + int(max_us * 20))
    p = np.convolve(np.abs(x20[start20:stop]) ** 2, np.ones(k) / k, "same")
    low = (p < noise_pwr * 4).astype(np.int8)
    if len(low) < g:
        return stop
    run = np.convolve(low, np.ones(g, int), "valid")
    idx = np.where(run[40:] == g)[0]                    # ignore the first 2 us
    return start20 + 40 + int(idx[0]) if idx.size else stop


def undecoded_wifi(det, x, pkts, same_tol_us=2.0, group_tol_us=5.0, min_ltf_flatness=0.35):
    """(1) Preambles the STF detector found but that could not be decoded
    (collisions, fades, heavy interference) -- clustered across scan centres.
    (2) Wi-Fi OFDM already in progress at capture start (no preamble in file),
    confirmed by numerology on channels where Wi-Fi was decoded."""
    out = []
    # ---- (1) undecodable preambles
    # Only at centres where a Wi-Fi signal can be grid-aligned: standard channels
    # and emitter centres discovered by decoding. Misaligned scan centres trigger
    # on every nearby preamble but can never decode it, so they add only copies.
    anchors = [c["freq_hz"] for c in channel_plan(det.fc, det.fs, 60e6)]
    anchors += [p.center_mhz * 1e6 for p in pkts if p.extra.get("off_grid")]
    cands = []
    # No standard channel in the capture and nothing decoded: there is no
    # centre where a Wi-Fi preamble could be grid-aligned, so skip this pass.
    for r in (det.rejects if anchors else []):
        if min(abs(r["freq_hz"] - a) for a in anchors) > 0.6e6:
            continue
        t = r["start20"] / 20; f = r["freq_hz"] + r["cfo"]
        # The STF stays periodic under ANY frequency shift, so one PPDU triggers every
        # scan centre whose 20 MHz window overlaps it (with a meaningless CFO there).
        # A trigger is explained if a decoded PPDU started at the same instant and
        # overlaps this decode window.
        g = r["freq_hz"]
        overl = [p for p in pkts if min(p.f_hi_mhz * 1e6, g + 10e6) - max(p.f_lo_mhz * 1e6, g - 10e6) > 0]
        if any(abs(p.start_us - t) < same_tol_us for p in overl):
            continue
        # Inside a decoded PPDU's data on an overlapping band: a narrowband interferer
        # (periodic at every lag) can mimic the STF there. Not reported. (A genuine
        # mid-PPDU collision preamble is lost from this list; it was undecodable anyway.)
        if any(p.start_us - same_tol_us < t < p.end_us for p in overl):
            continue
        if r.get("ltf_flatness", 1.0) < min_ltf_flatness:          # energy in a few bins: tone, not L-LTF
            continue
        cands.append(dict(r, t=t, f=f))
    cands.sort(key=lambda c: c["t"])
    # single-linkage: a wideband event triggers many scan centres; chain them
    groups = []
    for c in cands:
        hits = [g for g in groups if any(abs(m["t"] - c["t"]) < group_tol_us and abs(m["f"] - c["f"]) < 12e6
                                         for m in g)]
        if not hits:
            groups.append([c]); continue
        base = hits[0]; base.append(c)
        for g in hits[1:]:
            base.extend(g); groups.remove(g)
    cache = {}
    for g in groups:
        best = max(g, key=lambda c: c["stf"])
        off = best["freq_hz"] - det.fc
        if off not in cache:
            x20 = det.baseband(x, off)
            cache[off] = (x20, np.percentile(np.convolve(np.abs(x20) ** 2, np.ones(48) / 48, "valid"), 10))
        x20, floor = cache[off]
        s20 = max(0, best["start20"])
        e20 = _extent(x20, s20, floor)
        s64, oth = numerology(x20[s20 + 400:e20] if e20 - s20 > 600 else x20[s20:e20])
        f = best["f"]
        flo = min(min(c["f"] for c in g) - 8.3e6, f - 8.3e6); fhi = max(max(c["f"] for c in g) + 8.3e6, f + 8.3e6)
        flo = max(flo, det.fc - det.fs / 2); fhi = min(fhi, det.fc + det.fs / 2)
        out.append(Burst(label="wifi_preamble_undecodable", center_mhz=round(f / 1e6, 3),
                         start_sample=int(round(s20 * det.ratio)), end_sample=int(round(e20 * det.ratio)),
                         start_us=round(s20 / 20, 2), duration_us=round((e20 - s20) / 20, 2),
                         f_lo_mhz=round(flo / 1e6, 3), f_hi_mhz=round(fhi / 1e6, 3),
                         stf_metric=round(best["stf"], 3), cp64_score=round(s64, 3), cp_other_score=round(oth, 3),
                         reason=sorted({c["reason"] for c in g})[0], seen_on=len(g)))
    # ---- (2) truncated at capture start
    centres = sorted({round(p.center_mhz) for p in pkts})
    for fm in centres:
        off = fm * 1e6 - det.fc
        if abs(off) + 8.4e6 > det.fs / 2:
            continue
        if any(p.start_us < 25 and abs(p.center_mhz - fm) < 3 for p in pkts):
            continue
        x20 = det.baseband(x, off)
        floor = np.percentile(np.convolve(np.abs(x20) ** 2, np.ones(48) / 48, "valid"), 10)
        if np.mean(np.abs(x20[40:200]) ** 2) < floor * 4:
            continue
        e20 = _extent(x20, 0, floor)
        s64, oth = numerology(x20[:e20])
        if e20 < 320 or s64 < 0.75 or s64 - oth < 0.15:
            continue
        out.append(Burst(label="wifi_truncated_at_capture_start", center_mhz=float(fm), start_sample=0,
                         end_sample=int(round(e20 * det.ratio)), start_us=0.0, duration_us=round(e20 / 20, 2),
                         f_lo_mhz=round(fm - 8.3, 3), f_hi_mhz=round(fm + 8.3, 3), stf_metric=None,
                         cp64_score=round(s64, 3), cp_other_score=round(oth, 3),
                         reason="PPDU began before the capture; no preamble in file"))
    # overlapping truncated records on neighbouring centres -> keep strongest numerology
    out.sort(key=lambda b: (b.start_us, -b.cp64_score))
    merged = []
    for b in out:
        if any(m.label == b.label == "wifi_truncated_at_capture_start" and abs(m.center_mhz - b.center_mhz) < 12
               for m in merged):
            continue
        merged.append(b)
    return merged


def write_sigmf_annotations(meta, pkts, out_path, bursts=()):
    meta = json.loads(json.dumps(meta))
    meta["annotations"] = [{
        "core:sample_start": p.start_sample,
        "core:sample_count": p.end_sample - p.start_sample,
        "core:freq_lower_edge": p.f_lo_mhz * 1e6,
        "core:freq_upper_edge": p.f_hi_mhz * 1e6,
        "core:label": f"wifi_{p.fmt}",
        "core:comment": f"{p.channel} {p.bandwidth_mhz}MHz L-SIG {p.lsig_rate_mbps}Mbps len={p.lsig_length}"
                        + (f" MCS{p.mcs}" if p.mcs is not None else "") + f" SNR {p.snr_db}dB",
        "custom:wifi_format": p.fmt, "custom:snr_db": p.snr_db,
    } for p in pkts] + [{
        "core:sample_start": b.start_sample, "core:sample_count": b.end_sample - b.start_sample,
        "core:freq_lower_edge": b.f_lo_mhz * 1e6, "core:freq_upper_edge": b.f_hi_mhz * 1e6,
        "core:label": b.label, "core:comment": b.reason,
    } for b in bursts]
    meta["annotations"].sort(key=lambda a: a["core:sample_start"])
    json.dump(meta, open(out_path, "w"), indent=2)


def plot(x, fs, fc, pkts, path, nfft=1024, bursts=()):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    n = len(x) // nfft * nfft
    S = np.fft.fftshift(np.fft.fft(x[:n].reshape(-1, nfft) * np.hanning(nfft), axis=1), axes=1)
    S = 10 * np.log10(np.abs(S) ** 2 + 1e-12).T
    t_ms = n / fs * 1e3
    fig, ax = plt.subplots(figsize=(16, 7), dpi=110)
    ax.imshow(S, aspect="auto", origin="lower", cmap="viridis",
              extent=[0, t_ms, (fc - fs / 2) / 1e6, (fc + fs / 2) / 1e6],
              vmin=np.percentile(S, 5), vmax=np.percentile(S, 99.7))
    col = {"NON_HT": "#ff3b3b", "HT": "#ffd400", "VHT": "#ff00ff", "HE": "#00e5ff",
           "HT?": "#ffffff", "VHT?": "#ffffff"}
    for p in pkts:
        ax.add_patch(Rectangle((p.start_us / 1e3, p.f_lo_mhz), p.duration_us / 1e3,
                               p.f_hi_mhz - p.f_lo_mhz, fill=False, lw=1.1, ec=col[p.fmt]))
    for b in bursts:
        if True:
            ax.add_patch(Rectangle((b.start_us / 1e3, b.f_lo_mhz), b.duration_us / 1e3,
                                   b.f_hi_mhz - b.f_lo_mhz, fill=False, lw=1.1, ec="#ffa500", ls="--"))
    for f, c in col.items():
        ax.plot([], [], color=c, label=f)
    ax.plot([], [], color="#ffa500", ls="--", label="Wi-Fi, not decodable")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.6)
    ax.set_xlabel("time (ms)"); ax.set_ylabel("frequency (MHz)")
    ax.set_title(f"Wi-Fi PPDUs decoded: {len(pkts)}   |   Wi-Fi activity not decodable: {len(bursts)}")
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


# ============================================================== CLI
if __name__ == "__main__":
    import argparse, csv
    ap = argparse.ArgumentParser(description="Detect and localize Wi-Fi PPDUs in a SigMF capture")
    ap.add_argument("meta", help=".sigmf-meta path")
    ap.add_argument("--data", help=".sigmf-data path (default: alongside meta)")
    ap.add_argument("--out", default="wifi_packets")
    ap.add_argument("--usable_bw_mhz", type=float, default=56)
    ap.add_argument("--stf_thresh", type=float, default=0.75)
    ap.add_argument("--scan", action="store_true",
                    help="blind scan every 1 MHz (finds Wi-Fi on non-standard centre frequencies)")
    ap.add_argument("--no_plot", action="store_true")
    a = ap.parse_args()

    x, fs, fc, meta = load_sigmf(a.meta, a.data)
    det = WifiDetector(fs, fc, usable_bw=a.usable_bw_mhz * 1e6, stf_thresh=a.stf_thresh,
                       scan_step_hz=1e6 if a.scan else None)
    print(f"{len(x)} samples @ {fs/1e6:g} MSps, fc {fc/1e6:g} MHz, "
          f"{len(det.channels)} decode centres" + ("" if a.scan else f": {[c['name'] for c in det.channels]}"))
    pkts = det.run(x)

    rows = [asdict(p) for p in pkts]
    for r in rows:
        r["extra"] = json.dumps(r["extra"])
    with open(a.out + ".csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ["channel"])
        w.writeheader(); w.writerows(rows)
    json.dump([asdict(p) for p in pkts], open(a.out + ".json", "w"), indent=1)
    bursts = undecoded_wifi(det, x, pkts)
    with open(a.out + "_undecoded.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(asdict(bursts[0]).keys()) if bursts else ["label"])
        w.writeheader(); w.writerows(asdict(b) for b in bursts)
    write_sigmf_annotations(meta, pkts, a.out + ".sigmf-meta", bursts)
    if not a.no_plot:
        plot(x, fs, fc, pkts, a.out + ".png", bursts=bursts)
    from collections import Counter
    print(f"{len(pkts)} PPDUs:", dict(Counter(p.fmt for p in pkts)))
    print("undecoded Wi-Fi:", dict(Counter(b.label for b in bursts)))
