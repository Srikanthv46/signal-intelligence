"""Reference 802.11 OFDM transmitter (IEEE 802.11-2020 Clauses 17 and 19), for validating the receiver.

    ppdu(mpdu, mcs=0)                 -> complex baseband at 20 MSps, one HT (802.11n) PPDU
    ppdu(mpdu, mcs=7, ldpc=True, stbc=True) -> the same, LDPC-coded, Alamouti STBC from two antennas
    ppdu(mpdu, fmt="NON_HT", rate=6)  -> one non-HT (802.11a/g) PPDU
    ampdu([mpdu, ...])                -> an A-MPDU payload, delimiters and FCS included
    scene(ppdus, n, fs, fc)           -> a capture holding PPDUs at chosen channels, SNRs and offsets
    timescale(x, clock_hz, fs)        -> a 20 MSps PPDU played at another OFDM clock (5/10 MHz channels,
                                         or a vendor clock such as 11 MHz), resampled to fs

Only what one antenna can receive is generated: 20 MHz, one spatial stream (optionally STBC from two
transmit antennas, with the standard's cyclic shifts and a flat channel per antenna), BCC or LDPC. The
constellation and the codes are taken from the receiver's own tables (payload.pam_llr,
ofdm_rx.conv_encode, ldpc), so transmitter and receiver cannot disagree about a convention. The
conventions themselves are pinned separately: the QAM mapping by a hand-written table from the
standard (tests), and LDPC, STBC and the rotating HT pilots by CRC-valid decodes of real 802.11n
frames (testdata/D_20260211-161341).
"""
from __future__ import annotations

import numpy as np

from . import ofdm_rx as W
from . import payload as P

FS = W.FS_WIFI                      # 20 MSps
SERVICE_BITS = 16
TAIL_BITS = 6


# ------------------------------------------------------------------------------------ bits and bytes
def _bits(value, n):
    return [(value >> i) & 1 for i in range(n)]


def _byte_bits(data):
    return [b for byte in bytes(data) for b in _bits(byte, 8)]


def scramble(bits, seed=0x5D):
    reg = [(seed >> i) & 1 for i in range(7)]
    out = []
    for b in bits:
        s = reg[6] ^ reg[3]
        out.append(b ^ s)
        reg = [s] + reg[:6]
    return out


def puncture(coded, rate):
    """Inverse of payload.depuncture: drop the bits the standard's pattern steals."""
    p, q = rate
    pat = {(1, 2): [1, 1], (3, 4): [1, 1, 1, 0, 0, 1], (2, 3): [1, 1, 1, 0],
           (5, 6): [1, 1, 1, 0, 0, 1, 1, 0, 0, 1]}[(p, q)]
    mask = np.tile(np.array(pat, bool), int(np.ceil(len(coded) / len(pat))))[:len(coded)]
    return np.asarray(coded)[mask]


# ------------------------------------------------------------------------------------ constellation
def _axis_levels(per):
    """bit tuple -> PAM level, obtained by inverting the receiver's own demapper."""
    lv = np.arange(-(2 ** per - 1), 2 ** per, 2).astype(float)
    out = {}
    for x in lv:
        out[tuple(1 if v > 0 else 0 for v in P.pam_llr(x, per))] = x
    assert len(out) == len(lv), "demapper is not a bijection on its own constellation"
    return out


_LEVELS = {2: _axis_levels(2), 3: _axis_levels(3)}


def map_bits(bits, nbpsc):
    """Coded bits of one OFDM symbol -> complex constellation points (Gray, unit average power)."""
    b = np.asarray(bits, int).reshape(-1, nbpsc)
    k = P.KMOD[nbpsc]
    if nbpsc == 1:
        return (2 * b[:, 0] - 1) * k + 0j
    if nbpsc == 2:
        return ((2 * b[:, 0] - 1) + 1j * (2 * b[:, 1] - 1)) * k
    per = nbpsc // 2
    tab = _LEVELS[per]
    i = np.array([tab[tuple(r[:per])] for r in b])          # first half of the bits -> I
    q = np.array([tab[tuple(r[per:])] for r in b])          # second half -> Q
    return (i + 1j * q) * k


# ------------------------------------------------------------------------------------ OFDM symbols
def _symbol(data_sc, data_vals, pilot_sc, pilot_vals, cp, csd=0):
    """One OFDM symbol. csd: cyclic shift in 20 MSps samples (802.11 T_CS, e.g. -4 for -200 ns),
    applied to the 64-sample body before the cyclic prefix is taken from it."""
    X = np.zeros(64, complex)
    X[W._bins(data_sc)] = data_vals
    X[W._bins(pilot_sc)] = pilot_vals
    t = np.roll(np.fft.ifft(X) * 64 / np.sqrt(52), csd)
    return np.concatenate([t[-cp:], t])


def _interleave(bits, ncbps, nbpsc, ncol):
    j, _ = P.deinterleave_indices(ncbps, nbpsc, ncol)
    out = np.zeros(ncbps, int)
    out[j] = np.asarray(bits, int)                 # j: coded index -> subcarrier-bit index
    return out


def _legacy_preamble(csd=0):
    stf = np.tile(np.roll(W.STF_TIME, csd), 10)                       # 10 x 0.8 us
    lt = np.roll(W.LTF_TIME, csd)
    ltf = np.concatenate([lt[-32:], lt, lt])                          # GI2 + T1 + T2
    return np.concatenate([stf, ltf])


def _lsig(rate_mbps, length, csd=0):
    rb = next(k for k, v in W.RATES.items() if v[0] == rate_mbps)
    bits = list(rb) + [0] + _bits(length, 12)
    bits.append(sum(bits) & 1)                                        # even parity over the first 17
    bits += [0] * 6
    coded = W.conv_encode(np.array(bits))
    return _symbol(W.DATA_SC, map_bits(_interleave(coded, 48, 1, 16), 1),
                   W.PILOT_SC, W.PILOT_BASE * W.PILOT_POL[0], 16, csd)


def _htsig(mcs, length, short_gi, aggregation, n_ess=0, stbc=0, ldpc=0, csd=0):
    b = _bits(mcs, 7) + [0] + _bits(length, 16)                       # MCS, 20 MHz, HT length
    # smoothing, not sounding, reserved (1), aggregation, STBC (2 bits), FEC coding, short GI, N_ESS
    b += [1, 1, 1, aggregation] + _bits(stbc, 2) + [ldpc, short_gi] + _bits(n_ess, 2)
    crc = W.crc8_sig(b[:34])
    b = b[:34] + [int(x) for x in crc] + [0] * 6
    coded = W.conv_encode(np.array(b))
    out = []
    for k in range(2):
        # Only the data subcarriers rotate to QBPSK; the pilots stay BPSK, which is what lets a
        # receiver measure the common phase error and still see the 90 degrees that mark HT.
        sym = map_bits(_interleave(coded[48 * k:48 * (k + 1)], 48, 1, 16), 1) * 1j
        out.append(_symbol(W.DATA_SC, sym, W.PILOT_SC, W.PILOT_BASE * W.PILOT_POL[k + 1], 16, csd))
    return np.concatenate(out)


def _htstf(csd=0):
    s = np.zeros(57, complex)
    for k, v in zip([-24, -20, -16, -12, -8, -4, 4, 8, 12, 16, 20, 24],
                    [1, -1, 1, -1, -1, 1, -1, -1, 1, 1, 1, 1]):
        s[k + 28] = v * (1 + 1j)
    s *= np.sqrt(13 / 6)
    return _symbol(W.HT_SC, s, np.array([], int), np.array([]), 16, csd)


def _htltf(sign=1, csd=0):
    return _symbol(W.HT_SC, sign * W._HTLTF, np.array([], int), np.array([]), 16, csd)


# ------------------------------------------------------------------------------------ A-MPDU
def ampdu(mpdus, pad_last=True):
    """MPDUs (without FCS) -> an A-MPDU payload: delimiter + MPDU + FCS + 4-octet padding each."""
    out = bytearray()
    for m in mpdus:
        full = bytes(m) + P.crc32_802(bytes(m)).to_bytes(4, "little")
        o = (len(full) << 4) & 0xFFFF                                  # reserved(4) | length(12)
        d = bytes([o & 0xFF, o >> 8])
        out += d + bytes([P._crc8_delim(d), 0x4E]) + full
        if pad_last or m is not mpdus[-1]:
            out += b"\x00" * ((-len(out)) % 4)
    return bytes(out)


# ------------------------------------------------------------------------------------ PPDU
LEGACY_CSD = (0, -4)          # 802.11-2020 Table 19-9: 2 transmit chains, 0 and -200 ns
HT_CSD = (0, -8)              # Table 19-10: space-time streams 1 and 2, 0 and -400 ns


def ppdu(payload, mcs=0, fmt="HT", rate_mbps=6, short_gi=False, aggregation=False, scramble_seed=0x5D,
         ldpc=False, stbc=False, tx_gains=(1.0, 0.6 * np.exp(1j * 2.0))):
    """One PPDU as complex baseband at 20 MSps.

    HT: `payload` is an A-MPDU (aggregation=True) or one MPDU without FCS (the FCS is appended);
    `ldpc` selects LDPC coding, `stbc` Alamouti STBC from two antennas whose flat channels to the
    receiver are `tx_gains` (the returned signal is what one receive antenna sees).
    NON_HT: `payload` is one MPDU without FCS and `rate_mbps` selects the 802.11a/g rate.
    """
    from . import ldpc as LD
    if fmt == "HT":
        nbpsc, coderate = P.HT_MCS[mcs]
        ncbps, ncol, dsc = 52 * nbpsc, 13, W.HT_DATA_SC
    else:
        nbpsc, coderate = P.NONHT_MCS[rate_mbps]
        ncbps, ncol, dsc = 48 * nbpsc, 16, W.DATA_SC
        ldpc = stbc = False
    ndbps = ncbps * coderate[0] // coderate[1]
    nsts = 2 if stbc else 1

    psdu = bytes(payload) if (fmt == "HT" and aggregation) else \
        bytes(payload) + P.crc32_802(bytes(payload)).to_bytes(4, "little")
    length = len(psdu)

    if ldpc:
        LP = LD.params(length, ncbps, coderate, nsts)
        n_sym = LP["nsym"]
        bits = scramble([0] * SERVICE_BITS + _byte_bits(psdu), scramble_seed)
        coded = LD.encode_stream(bits, LP, coderate)
    else:
        n_sym = nsts * int(np.ceil((SERVICE_BITS + 8 * length + TAIL_BITS) / (nsts * ndbps)))
        bits = [0] * SERVICE_BITS + _byte_bits(psdu)
        bits += [0] * (n_sym * ndbps - len(bits))
        bits = scramble(bits, scramble_seed)
        tail0 = SERVICE_BITS + 8 * length
        bits[tail0:tail0 + TAIL_BITS] = [0] * TAIL_BITS                # tail is zeroed after scrambling
        coded = puncture(W.conv_encode(np.array(bits, int)), coderate)
    cp = 8 if (fmt == "HT" and short_gi) else 16

    # constellation points of every data symbol (one space-time stream)
    pts = []
    for k in range(n_sym):
        cb = coded[ncbps * k:ncbps * (k + 1)]
        pts.append(map_bits(cb if ldpc else _interleave(cb, ncbps, nbpsc, ncol), nbpsc))

    if fmt != "HT":
        data = [_symbol(dsc, pts[k], W.PILOT_SC, W.PILOT_BASE * W.PILOT_POL[(k + 1) % 127], cp)
                for k in range(n_sym)]
        return np.concatenate([_legacy_preamble(), _lsig(rate_mbps, length)] + data)

    txtime = 20 + 8 + 4 + 4 * nsts + (3.6 if short_gi else 4.0) * n_sym
    l_len = 3 * int(np.ceil((txtime - 20) / 4)) - 3
    streams = []
    for i in range(nsts):
        if stbc:                                                       # Alamouti over symbol pairs
            st = []
            for m in range(n_sym // 2):
                d0, d1 = pts[2 * m], pts[2 * m + 1]
                st += [d0, d1] if i == 0 else [-np.conj(d1), np.conj(d0)]
        else:
            st = pts
        lc, hc = LEGACY_CSD[i] if nsts > 1 else 0, HT_CSD[i] if nsts > 1 else 0
        head = [_legacy_preamble(lc), _lsig(6, l_len, lc),
                _htsig(mcs, length, int(short_gi), int(aggregation), stbc=int(stbc), ldpc=int(ldpc), csd=lc),
                _htstf(hc)] + [_htltf(P.HT_LTF_P[i, n] if nsts > 1 else 1, hc) for n in range(nsts)]
        data = [_symbol(dsc, st[k], W.PILOT_SC, P.ht_pilots(k, nsts)[i], cp, hc) for k in range(n_sym)]
        streams.append(np.concatenate(head + data) / np.sqrt(nsts))
    if nsts == 1:
        return streams[0]
    return tx_gains[0] * streams[0] + tx_gains[1] * streams[1]


def he_sig_a_bits(fmt="SU", uplink=False, mcs=5, bss_color=14, bw=0, gi_ltf=1, nsts=1, ldpc=1, txop=127,
                  sigb_mcs=0):
    """52 HE-SIG-A bits (802.11ax-2021 Tables 27-18 / 27-20 / 27-21) with CRC-4 and tail. bw: the 2-bit
    (SU / TB) or 3-bit (MU) code; txop: the raw 7-bit field."""
    if fmt == "SU":
        a1 = [1, 0, int(uplink)] + _bits(mcs, 4) + [0] + _bits(bss_color, 6) + [1] + _bits(0, 4) + _bits(bw, 2) \
            + _bits(gi_ltf, 2) + _bits(nsts - 1, 3)
        a2 = _bits(txop, 7) + [ldpc, 0, 0, 0] + _bits(0, 2) + [0, 1, 0]
    elif fmt == "MU":
        a1 = [int(uplink)] + _bits(sigb_mcs, 3) + [0] + _bits(bss_color, 6) + _bits(0, 4) + _bits(bw, 3) \
            + _bits(1, 4) + [0] + _bits(gi_ltf, 2) + [0]
        a2 = _bits(txop, 7) + [1] + _bits(1, 3) + [0, 0] + _bits(0, 2) + [0]
    else:                                                     # TB
        a1 = [0] + _bits(bss_color, 6) + _bits(0, 16) + [1] + _bits(bw, 2)
        a2 = _bits(txop, 7) + [1] * 9
    b = a1 + a2
    assert len(b) == 42
    return b + [int(x) for x in W.crc8_sig(b)[:4]] + [0] * 6


def he_preamble(sig_bits, mu=False, n_filler=12, seed=0):
    """An HE PPDU's legacy part, RL-SIG and HE-SIG-A, then `n_filler` 4 us filler symbols (random QPSK on 56
    tones) standing in for HE-STF / HE-LTF / data, at 20 MSps. L-SIG / RL-SIG carry -1 -1 -1 +1 on +-27 / +-28;
    LENGTH = 3 N - 3 - m (m = 1 MU, 2 otherwise), N = the 4 us symbols after L-SIG."""
    n = 3 + n_filler
    length = 3 * n - 3 - (1 if mu else 2)
    rb = next(k for k, v in W.RATES.items() if v[0] == 6)
    lb = list(rb) + [0] + _bits(length, 12)
    lb.append(sum(lb) & 1)
    lcoded = map_bits(_interleave(W.conv_encode(np.array(lb + [0] * 6)), 48, 1, 16), 1)
    dsc = np.r_[W.DATA_SC, W.HE_EXT_SC]
    sig = [_symbol(dsc, np.r_[lcoded, W.HE_EXT_VAL], W.PILOT_SC, W.PILOT_BASE * W.PILOT_POL[k], 16)
           for k in range(2)]                                # L-SIG, RL-SIG
    coded = W.conv_encode(np.array(sig_bits))
    for k in range(2):
        sig.append(_symbol(W.HT_DATA_SC, map_bits(_interleave(coded[52 * k:52 * (k + 1)], 52, 1, 13), 1),
                           W.PILOT_SC, W.PILOT_BASE * W.PILOT_POL[k + 2], 16))
    rng = np.random.default_rng(seed)
    fill = [_symbol(W.HT_USED_SC, (rng.choice([-1, 1], 56) + 1j * rng.choice([-1, 1], 56)) / np.sqrt(2),
                    np.array([], int), [], 16) for _ in range(n_filler)]
    return np.concatenate([_legacy_preamble()] + sig + fill)


def timescale(x, clock_hz, fs=FS):
    """A PPDU generated at 20 MSps, transmitted with an OFDM clock of clock_hz (every duration
    stretched by 20 MHz / clock_hz, every bandwidth shrunk by the same factor) and resampled to fs.
    This is how 802.11 5/10 MHz channels (half/quarter clocking) and vendor links such as an 11 MHz
    clock are built: the same baseband, a different sample clock."""
    from fractions import Fraction
    from scipy.signal import resample_poly
    r = Fraction(fs / clock_hz).limit_denominator(4096)
    return resample_poly(np.asarray(x), r.numerator, r.denominator)


# ------------------------------------------------------------------------------------ channel
def scene(ppdus, n, fs=61.44e6, fc=2437e6, noise=1.0, seed=0):
    """A capture holding OFDM PPDUs. Each is a dict: samples (from ppdu()), f_hz, t_s, snr_db, and
    optional cfo_hz. `snr_db` is measured in the 20 MHz channel."""
    from scipy.signal import resample_poly
    from fractions import Fraction

    rng = np.random.default_rng(seed)
    x = (rng.normal(size=n) + 1j * rng.normal(size=n)) * np.sqrt(noise / 2)
    r = Fraction(fs / FS).limit_denominator(4096)
    for p in ppdus:
        y = resample_poly(p["samples"], r.numerator, r.denominator)
        i0 = int(p["t_s"] * fs)
        m = min(len(y), n - i0)
        if m <= 0:
            continue
        y = y[:m]
        a = np.sqrt(10 ** (p["snr_db"] / 10) * noise * 20e6 / fs / (np.mean(np.abs(y) ** 2) + 1e-30))
        f = p["f_hz"] - fc + p.get("cfo_hz", 0.0)
        x[i0:i0 + m] += a * y * np.exp(2j * np.pi * f * (np.arange(m) + i0) / fs
                                       + 1j * rng.uniform(0, 2 * np.pi))
    return x.astype(np.complex64)
