#!/usr/bin/env python3
"""
wifi_payload.py -- decode the PAYLOAD and MAC header of 802.11 packets found by
wifi_detect.py.

Scope and why
-------------
A single-antenna capture can fully decode a packet only when it is:
  * 20 MHz, single spatial stream, BCC (convolutional) coded.
That is EVERY non-HT (802.11a/g) packet, which in this capture is the bulk of
traffic -- and, importantly, it is where the MAC-layer control and management
frames live (ACK, RTS/CTS, BlockAck, and any legacy data/mgmt frames). HT/VHT
data frames here are 80 MHz and/or 2-3 spatial streams and/or LDPC, so their
payload is not recoverable from one RX chain; their MAC identity, when present,
tends to appear in the associated non-HT control frames anyway.

The MAC header is NOT encrypted on WPA2/WPA3 -- only the frame body is. So from
a decoded frame we recover: frame type/subtype, the MAC addresses (who talks to
whom, on which BSS), sequence/retry/power bits, the NAV duration, and -- for
beacons/probes, which are entirely unencrypted -- the SSID and capabilities.

The non-HT data-plane chain implemented here (IEEE 802.11-2020 Clause 17):
  equalize each data OFDM symbol (reuse detector's H, per-SC weights, pilots)
  -> soft demap (BPSK/QPSK/16-QAM/64-QAM, Gray)
  -> per-symbol deinterleave (depends on N_BPSC/N_CBPS, NOT the L-SIG interleaver)
  -> depuncture to rate 1/2 (1/2, 3/4, 2/3 puncturing patterns)
  -> soft Viterbi (K=7)
  -> descramble (7-bit LFSR x^7+x^4+1; seed recovered from first 7 SERVICE bits)
  -> strip 16 SERVICE bits + 6 tail + pad -> MPDU bytes
  -> check FCS (CRC-32) -> parse MAC header
"""
from __future__ import annotations
import json
import zlib

import numpy as np

from . import mac80211 as M
from . import ofdm_rx as W
FS = W.FS_WIFI

# ---- MCS table for non-HT: rate(Mbps) -> (N_BPSC, coding rate as (p,q) for p/q)
NONHT_MCS = {6: (1, (1, 2)), 9: (1, (3, 4)), 12: (2, (1, 2)), 18: (2, (3, 4)),
             24: (4, (1, 2)), 36: (4, (3, 4)), 48: (6, (2, 3)), 54: (6, (3, 4))}


# ---------------------------------------------------------------- interleaver
def deinterleave_indices(ncbps, nbpsc, ncol=16):
    """Return perm such that coded_bit[k] = received_bit[perm[k]] for one symbol.
    802.11 has two permutations; we invert their composition. ncol is N_COL: 16 for a non-HT
    20 MHz symbol (48 data subcarriers), 13 for an HT 20 MHz symbol (52). The third HT
    permutation (frequency rotation) applies only to a second spatial stream, so at one
    stream it is the identity and is not needed here."""
    s = max(nbpsc // 2, 1)
    k = np.arange(ncbps)
    i = (ncbps // ncol) * (k % ncol) + (k // ncol)                   # first permutation
    j = s * (i // s) + (i + ncbps - (ncol * i // ncbps)) % s         # second permutation
    # j maps coded index k -> subcarrier-bit index. We need the inverse for gather.
    inv = np.empty(ncbps, int)
    inv[j] = k
    return j, inv


# ---------------------------------------------------------------- demapper (soft)
KMOD = {1: 1, 2: 1 / np.sqrt(2), 4: 1 / np.sqrt(10), 6: 1 / np.sqrt(42)}


def pam_llr(x, bits):
    """Gray-coded PAM LLRs per axis (min-distance approximation). Also defines the constellation:
    ofdm_tx inverts this, so transmitter and receiver cannot disagree about the bit mapping."""
    if bits == 2:                                                   # 16-QAM axis: b0,b1
        return [x, 2 - np.abs(x)]
    return [x, 4 - np.abs(x), 2 - np.abs(4 - np.abs(x))]            # 64-QAM axis: b0,b1,b2


def demap_soft(sym, nbpsc, noise=1.0):
    """Approximate per-bit LLRs for Gray-coded BPSK/QPSK/16/64-QAM.
    Returns N_BPSC soft values per symbol (sign = bit estimate, |.| ~ confidence)."""
    Kmod = KMOD[nbpsc]
    I = sym.real / Kmod; Q = sym.imag / Kmod
    if nbpsc == 1:
        return np.stack([I], -1)
    if nbpsc == 2:
        return np.stack([I, Q], -1)
    per = nbpsc // 2
    # 802.11-2020 Figure 17-10/17-11: a subcarrier's N_BPSC bits split in half, the FIRST half
    # choosing the I level and the SECOND half the Q level -- (b0,b1 | b2,b3) for 16-QAM, not
    # I and Q interleaved. Interleaving them decodes BPSK and QPSK correctly (one bit per axis)
    # but silently destroys every 16-QAM and 64-QAM frame: 24/36/48/54 Mb/s and HT MCS 3-7.
    return np.stack(pam_llr(I, per) + pam_llr(Q, per), -1)


# ---------------------------------------------------------------- depuncture
def depuncture(soft, rate):
    """Insert zero-LLR erasures where the transmitter punctured, giving a rate-1/2
    stream the Viterbi decoder expects. Patterns from 802.11 Table 17-9..11."""
    p, q = rate
    if (p, q) == (1, 2):
        return soft
    if (p, q) == (3, 4):
        pat = [1, 1, 1, 0, 0, 1]                                    # keep A0 B0 A1 -- -- B2
    elif (p, q) == (2, 3):
        pat = [1, 1, 1, 0]                                          # keep A0 B0 A1 --
    elif (p, q) == (5, 6):                                          # 802.11n MCS7: 10 coded bits -> 6
        pat = [1, 1, 1, 0, 0, 1, 1, 0, 0, 1]                        # keep A0 B0 A1 B2 A3 B4
    else:
        return soft
    pat = np.array(pat, bool)
    nin = len(soft)
    reps = int(np.ceil(nin / pat.sum()))
    mask = np.tile(pat, reps)
    keep = np.where(mask)[0]
    out = np.zeros(len(mask), soft.dtype)
    out[keep[:nin]] = soft[:nin]
    # trim to a whole number of output pairs
    return out[:(len(out) // 2) * 2]


# ---------------------------------------------------------------- scrambler
def scrambler_seq(seed7, n):
    s = list(seed7)
    out = []
    for _ in range(n):
        fb = s[6] ^ s[3]
        out.append(fb)
        s = [fb] + s[:6]
    return np.array(out, np.uint8)


def descramble(bits):
    """Recover the 7-bit seed from the first 7 scrambled SERVICE bits (which are
    known zeros pre-scrambling), then descramble the whole stream."""
    b = bits.astype(np.uint8)
    # SERVICE[0:7] are 0 before scrambling, so the scrambler output there == the
    # scrambling sequence. Invert the LFSR to get the seed.
    # scrambler out x_n = seed-driven; first 7 outputs equal s reversed via recurrence.
    # Standard trick: seed bits = descrambling of first 7 known-zero bits.
    s = [0] * 7
    seq = []
    # Solve: for known input 0, scrambled = seq. seq[n] = s[6]^s[3]; then shift in seq[n].
    # Given scrambled bits c[0..6] with input 0: c[n] = state feedback; state updates with c[n].
    st = [0] * 7
    # initialize state from c: state after producing c[n] shifts c[n] in.
    c = b[:7]
    # Recover initial state by running the register backwards is fiddly; instead
    # brute-force the 127 non-zero seeds (cheap, unambiguous).
    best = None
    for seed in range(1, 128):
        s7 = [(seed >> i) & 1 for i in range(7)]
        seq = scrambler_seq(s7, 7)
        if np.array_equal(seq, c):
            best = s7; break
    if best is None:
        best = [1, 0, 1, 1, 1, 0, 1]                                # fallback (spec default)
    full = scrambler_seq(best, len(b))
    return b ^ full, best


# ---------------------------------------------------------------- FCS
def crc32_802(data):
    """IEEE 802.11 FCS (CRC-32, reflected, as zlib)."""
    return zlib.crc32(bytes(bytearray(int(b) & 0xFF for b in data)))


# ---------------------------------------------------------------- MAC parse (mac80211.py)
FTYPE, SUBTYPE = M.FTYPE, M.SUBTYPE
vendor = M.vendor


def parse_mac(mpdu):
    """Dissect an 802.11 MAC frame (mpdu excludes the FCS). See mac80211.parse."""
    return M.parse(mpdu)


# ---------------------------------------------------------------- data-plane decode
def decode_nonht_payload(rx, rate_mbps, length, chan_offset_hz, usable_half_hz):
    """Reuse the detector's per-packet state to demodulate and decode all data
    OFDM symbols of a non-HT PPDU. Returns (mpdu_bytes, fcs_ok, seed, info)."""
    r, start, cfo, back = rx["r"], rx["start"], rx["cfo"], rx["back"]
    H, W_rel, good_sc = rx["H"], rx["W_rel"], rx["good_sc"]
    used = W._bins(good_sc)
    nbpsc, coderate = NONHT_MCS[rate_mbps]
    ncbps = 48 * nbpsc
    ndbps = ncbps * coderate[0] // coderate[1]
    n_sym = int(np.ceil((16 + 8 * length + 6) / ndbps))            # DATA OFDM symbols

    # data subcarriers actually usable (erasure-aware)
    data_bins = W._bins(W.DATA_SC)
    live_data = np.isin(data_bins, used)
    # interleaver perm on the FULL 48-subcarrier grid; erased carriers -> 0 LLR later
    j_perm, inv = deinterleave_indices(ncbps, nbpsc)

    # DATA symbols start 400 samples after `start`:
    #   L-STF(160) + L-LTF(160) + L-SIG(80) = 400. Symbol k occupies
    #   [400+80k, 400+80k+80); its 16-sample cyclic prefix is first, so the
    #   64-sample body is [400+80k+16, 400+80k+80). fft_at(o) reads
    #   seg[o-back : o-back+64], matching the detector's `back`-into-CP window,
    #   so the body-aligned call is o = 400+80k+16+back.
    seg = r[start:start + 400 + 80 * n_sym + 80]
    nn = np.arange(len(seg))
    seg = seg * np.exp(-2j * np.pi * cfo / FS * (start + nn))

    def fft_at(off):
        w = seg[off - back:off - back + 64]
        if len(w) < 64:
            w = np.pad(w, (0, 64 - len(w)))
        return np.fft.fft(w) / 64 * np.sqrt(52)

    pb = W._bins(W.PILOT_SC)
    pbf = W.PILOT_SC                                                # pilot subcarrier indices (-21..21)

    # ---- equalise each DATA symbol, then remove its common phase error (CPE) from the pilots.
    # A residual carrier offset rotates every subcarrier by an angle that grows with the symbol
    # index; the per-symbol CPE removes it. A sampling-clock offset would add a phase slope across
    # subcarriers, but at Wi-Fi tolerances it stays below ~0.15 rad over a 2 ms PPDU (measured on
    # the 2026-10-04 sweep), and fitting it from 4 pilots per symbol cost more than it gained.
    Zsym = []
    for k in range(n_sym):
        Y = fft_at(400 + 80 * k + 16)                              # fft_at applies `back` itself
        Z = np.zeros(64, complex); Z[used] = Y[used] / H[used]
        Zsym.append(Z)
    all_soft = []
    for k in range(n_sym):
        Z = Zsym[k]
        pil = Z[pb] / (W.PILOT_BASE * W.PILOT_POL[(k + 1) % 127])
        w = W_rel[pb].copy(); live = w > 0
        if live.sum() >= 3:                                        # drop one outlier pilot
            ph0 = np.angle(np.sum(w * pil))
            err = np.abs(pil * np.exp(-1j * ph0) - 1); err[~live] = -1
            if err.max() > 0.8:
                w[np.argmax(err)] = 0
        ph = np.angle(np.sum(w * pil))                             # reliability-weighted CPE
        Z = Z * np.exp(-1j * ph)
        d = Z[data_bins]
        conf = W_rel[data_bins]
        # soft demap per data subcarrier -> (48, nbpsc) LLRs, scaled by SC confidence
        llr = demap_soft(d, nbpsc) * conf[:, None]
        llr[~live_data] = 0.0                                       # erased subcarriers = no info
        # arrange into coded-bit order for the symbol: subcarrier-major, bit-minor
        sym_bits = llr.reshape(-1)                                  # length ncbps, order (sc0b0,sc0b1,...)
        # deinterleave: received subcarrier-bit index -> coded index
        coded = np.zeros(ncbps, llr.dtype)
        coded[inv] = sym_bits
        all_soft.append(coded)
    soft = np.concatenate(all_soft) if all_soft else np.zeros(0)

    # depuncture to rate 1/2, then Viterbi
    soft12 = depuncture(soft, coderate)
    databits = W.viterbi(soft12)                                   # includes SERVICE+MPDU+tail+pad
    # descramble
    descr, seed = descramble(databits)
    payload = descr[16:16 + 8 * length]                            # drop 16 SERVICE bits, take LENGTH bytes
    if len(payload) < 8 * length:
        payload = np.pad(payload, (0, 8 * length - len(payload)))
    mpdu = np.packbits(payload.reshape(-1, 8)[:, ::-1])            # LSB-first per octet
    fcs_ok = False
    info = {}
    if length >= 4:
        body, fcs = mpdu[:-4], mpdu[-4:]
        calc = crc32_802(body)
        rx_fcs = int(fcs[0]) | (int(fcs[1]) << 8) | (int(fcs[2]) << 16) | (int(fcs[3]) << 24)
        fcs_ok = (calc == rx_fcs)
        info = parse_mac(body)
    # SERVICE bits 7..15 are zero before scrambling: descrambled to zero they vouch for the seed even
    # when the payload fails its FCS (links.attribute uses seeds as a transmitter fingerprint)
    info["_service_ok"] = bool(len(descr) >= 16 and not descr[7:16].any())
    return mpdu, fcs_ok, seed, info, n_sym


# ---------------------------------------------------------------- VHT (802.11ac) Partial AID
# An 80 MHz VHT payload is wider than a Pluto can capture, but VHT-SIG-A sits on the primary 20 MHz
# and names the receiver (802.11-2020 26.17.x). That attributes otherwise opaque data traffic.
def _octets(mac):
    return [int(x, 16) for x in mac.split(":")]


def vht_uplink_paid(bssid):
    """Partial AID a station puts in an uplink VHT PPDU (Group ID 0): dec(BSSID[39:47])."""
    o = _octets(bssid)
    return ((o[4] >> 7) | (o[5] << 1)) & 0x1FF


def vht_downlink_aid(paid, bssid):
    """Association ID (low 9 bits) of the station an AP addresses with Group ID 63:
    PARTIAL_AID = (AID[0:8] + (BSSID[44:47] xor BSSID[40:43]) * 2^5) mod 2^9, solved for AID."""
    o5 = _octets(bssid)[5]
    return (paid - (((o5 >> 4) ^ (o5 & 0xF)) & 0xF) * 32) % 512


# ---------------------------------------------------------------- 802.11n (HT), one spatial stream
# MCS 0-7 is the single-stream set: one antenna is enough, and it carries most 2.4 GHz phone traffic.
# MCS >= 8 is 2+ spatial streams, which one antenna cannot separate: declined, not guessed.
HT_MCS = {0: (1, (1, 2)), 1: (2, (1, 2)), 2: (2, (3, 4)), 3: (4, (1, 2)), 4: (4, (3, 4)),
          5: (6, (2, 3)), 6: (6, (3, 4)), 7: (6, (5, 6))}
HT_MOD = {1: "BPSK", 2: "QPSK", 4: "16-QAM", 6: "64-QAM"}
HT_MAX_LENGTH = 8192          # bytes; bounds decode work per PPDU (a 4 ms A-MPDU at MCS7 is ~6.5 kB)


def ht_rate_mbps(mcs, short_gi=False):
    nbpsc, (p, q) = HT_MCS[mcs]
    return round(52 * nbpsc * p / q / (3.6 if short_gi else 4.0), 1)


def ht_decodable(mcs, ldpc=0, stbc=0, n_ess=0, bw_mhz=20):
    """Can one antenna decode this HT PPDU? Returns (bool, reason).

    One spatial stream (MCS 0-7) is the limit of a single receive antenna. Within it, BCC and LDPC
    coding are both decoded, and so is Alamouti STBC (STBC=1: one stream sent as two space-time
    streams from two antennas), which a single antenna combines with full diversity."""
    if mcs is None or mcs not in HT_MCS:
        return False, f"MCS {mcs}: 2+ spatial streams (needs as many antennas)"
    if bw_mhz != 20:
        return False, f"{bw_mhz} MHz HT payload (only 20 MHz implemented)"
    if stbc > 1:
        return False, f"space-time block coding with STBC={stbc} (N_STS > 2 for one stream)"
    if n_ess:
        return False, "extension spatial streams"
    return True, ""


def _crc8_delim(two):
    """CRC-8 over an A-MPDU delimiter's first 2 octets (same polynomial as HT-SIG)."""
    bits = [(two[i] >> b) & 1 for i in range(2) for b in range(8)]      # transmission order, LSB first
    c = W.crc8_sig(bits)
    return sum(int(c[i]) << i for i in range(8))


def deaggregate(data):
    """A-MPDU -> [(mpdu_bytes, delimiter_crc_ok)].

    Each subframe is a 4-octet delimiter -- reserved(4) | MPDU length(12) | CRC-8 | signature 0x4E --
    then the MPDU, then padding to a 4-octet boundary; a zero length is pure padding. On a corrupted
    delimiter we resynchronise on the next signature byte rather than abandoning the rest of the A-MPDU.
    """
    data = bytes(data)
    out, i = [], 0
    while i + 4 <= len(data):
        d = data[i:i + 4]
        if d[3] != 0x4E:
            nxt = data.find(b"\x4e", i + 4)
            if nxt < 0:
                break
            i = nxt - 3                                   # >= i+1, so this always makes progress
            continue
        ln = (((d[1] << 8) | d[0]) >> 4) & 0xFFF
        ok = _crc8_delim(d[:2]) == d[2]
        i += 4
        if ln:
            out.append((data[i:i + ln], ok))
            i += ln
            i += (-i) % 4
    return out


# HT pilots (802.11-2020 Eq. 19-58 and Table 19-19): pilot j of data symbol n on space-time stream i is
# PSI[i][(n + j) mod 4] * p_{n+3}. The pattern ROTATES by one position every symbol; a fixed
# [1, 1, 1, -1] is right only on every fourth symbol.
HT_PILOT_PSI = {1: np.array([[1, 1, 1, -1]]), 2: np.array([[1, 1, -1, -1], [1, -1, -1, 1]])}
HT_LTF_P = np.array([[1, -1], [1, 1]])        # HT-LTF n of space-time stream i is multiplied by P[i, n]


def ht_pilots(n, nsts):
    """(nsts, 4) pilot values of HT data symbol n (subcarriers -21, -7, 7, 21)."""
    psi = HT_PILOT_PSI[nsts]
    return psi[:, (n + np.arange(4)) % 4] * W.PILOT_POL[(n + 3) % 127]


def _levels(per):
    """PAM levels and their bit labels, from pam_llr (the single definition of the mapping)."""
    lv = np.arange(-(2 ** per - 1), 2 ** per, 2).astype(float)
    lab = np.array([[1 if v > 0 else 0 for v in pam_llr(x, per)] for x in lv])
    return lv, lab


_PAM = {2: _levels(2), 3: _levels(3)}


def llr_exact(sym, nbpsc, var):
    """Log-MAP LLRs, log P(0)/P(1), of unit-power Gray QAM points. sym (n,), var (n,): complex noise
    variance of each point. Returns (n, nbpsc) in transmit bit order (first half I, second half Q)."""
    k = KMOD[nbpsc]
    v = np.maximum(np.asarray(var, float), 1e-9)
    if nbpsc == 1:
        return (-4 * sym.real * k / v)[:, None]
    if nbpsc == 2:
        return np.stack([-4 * sym.real * k / v, -4 * sym.imag * k / v], -1)
    per = nbpsc // 2
    lv, lab = _PAM[per]
    s2 = v / 2 / k ** 2                                         # per-axis variance in PAM units
    out = []
    for axis in (sym.real, sym.imag):
        d = -(axis[:, None] / k - lv[None, :]) ** 2 / (2 * s2[:, None])
        for b in range(per):
            l0 = np.logaddexp.reduce(np.where(lab[None, :, b] == 0, d, -np.inf), axis=1)
            l1 = np.logaddexp.reduce(np.where(lab[None, :, b] == 1, d, -np.inf), axis=1)
            out.append(l0 - l1)
    return np.stack(out, -1)


def smooth_channel(h, sc, w):
    """Average a channel estimate over 2w+1 neighbouring subcarriers after removing its linear phase
    (the delay: timing offset plus any cyclic shift), then put the phase back. h (64,), sc the used
    subcarrier indices (-28..28). A single HT-LTF gives a noisy estimate; on a channel that is flat
    over w*312.5 kHz this buys up to 10*log10(2w+1) dB of estimation SNR."""
    if w <= 0:
        return h
    sc = np.sort(np.asarray(sc))
    b = W._bins(sc)
    v = h[b]
    step = np.diff(sc)
    tau = np.angle(np.sum(v[1:] * np.conj(v[:-1]) * (step == 1)))       # phase per subcarrier
    g = v * np.exp(-1j * tau * sc)
    c = np.convolve(np.concatenate([np.full(w, g[0]), g, np.full(w, g[-1])]), np.ones(2 * w + 1) / (2 * w + 1),
                    "valid")
    out = h.copy()
    out[b] = c * np.exp(1j * tau * sc)
    return out


def _track(Y, pred, psc):
    """Common phase and phase slope (timing drift) of one symbol from its pilots: Y, pred (4,)."""
    z = Y * np.conj(pred)
    a0 = np.angle(z.sum())
    d = np.angle(np.exp(1j * (np.angle(z) - a0)))
    w = np.abs(pred) ** 2
    b = float(np.sum(w * d * psc) / (np.sum(w * psc ** 2) + 1e-12))
    return a0, b


def decode_ht_payload(rx, mcs, length, short_gi=False, aggregation=False, n_ltf=1, n_ess=0,
                      chan_offset_hz=0.0, usable_half_hz=1e12, stbc=0, ldpc=0, smoothing=1):
    """Demodulate and decode the HT-DATA of a 20 MHz, single-stream 802.11n PPDU: BCC or LDPC, with or
    without Alamouti STBC.

    Differences from the non-HT chain, all of them structural (802.11-2020 Clause 19):
      * 52 data subcarriers over -28..28, not 48 over -26..26, so N_CBPS = 52*N_BPSC; BCC uses the
        interleaver with N_COL = 13, LDPC maps coded bits to subcarriers in order (no interleaver);
      * the channel is re-estimated from the HT-LTFs, the only fields that train the four extra
        subcarriers. With STBC there are two space-time streams (two transmit antennas, the second
        with a -400 ns cyclic shift) and two HT-LTFs; P_HTLTF = [[1, -1], [1, 1]] separates them;
      * STBC: each pair of data symbols carries (d0, d1) on stream 1 and (-d1*, d0*) on stream 2, so
        one antenna recovers both with Alamouti combining;
      * pilots start at p_3 and their pattern rotates by one position per symbol (HT_PILOT_PSI);
        each symbol's common phase and phase slope (sampling-clock drift, ~10 ppm on real links) are
        measured on the pilots and smoothed over neighbouring symbols;
      * a short guard interval makes the symbol 72 samples instead of 80; MCS7 BCC is punctured 5/6;
      * LDPC: codeword count, shortening, puncturing and repetition from the PSDU length
        (ldpc.params), soft log-MAP LLRs, normalised min-sum decoding.
      * the payload is usually an A-MPDU, so it is de-aggregated and each MPDU checked on its own FCS.

    Returns a dict: n_sym, aggregated, subframes [{mpdu_hex, fcs_ok, mac}], n_fcs_ok, coding, stbc,
    scrambler_seed (int 1..127) and service_ok (the 9 reserved SERVICE bits descramble to zero, which
    vouches for the seed even when the payload fails), and for LDPC n_codewords / n_codewords_ok.
    """
    from . import ldpc as LD
    ok, why = ht_decodable(mcs, ldpc, stbc, n_ess)
    if not ok:
        raise ValueError(why)
    if length > HT_MAX_LENGTH:
        raise ValueError(f"HT length {length} B over the {HT_MAX_LENGTH} B decode budget")
    r, start, cfo, back = rx["r"], rx["start"], rx["cfo"], rx["back"]
    nbpsc, coderate = HT_MCS[mcs]
    ncbps = 52 * nbpsc
    ndbps = ncbps * coderate[0] // coderate[1]
    nsts = 2 if stbc else 1
    m_stbc = 2 if stbc else 1
    if ldpc:
        LP = LD.params(length, ncbps, coderate, m_stbc)
        n_sym = LP["nsym"]
    else:
        n_sym = m_stbc * int(np.ceil((16 + 8 * length + 6) / (m_stbc * ndbps)))
    n_ltf = max(n_ltf, nsts)

    # L-STF 160 | L-LTF 160 | L-SIG 80 | HT-SIG 160 | HT-STF 80 | HT-LTF 80 each | data
    ltf1 = 640
    data0 = ltf1 + 80 * (n_ltf + n_ess)
    sym, gi = (72, 8) if short_gi else (80, 16)
    need = data0 + sym * n_sym + 96
    seg = r[start:start + need]
    if len(seg) < data0 + 64:
        raise ValueError("PPDU truncated before HT-LTF")
    seg = seg * np.exp(-2j * np.pi * cfo / W.FS_WIFI * (start + np.arange(len(seg))))

    def fft_at(off):
        w = seg[off - back:off - back + 64]
        if len(w) < 64:
            w = np.pad(w, (0, 64 - len(w)))
        return np.fft.fft(w) / 64 * np.sqrt(52)

    # ---- channel per space-time stream from the HT-LTFs, on the 56 HT subcarriers inside the usable band
    sc_f = chan_offset_hz + W.HT_USED_SC * 312.5e3
    good = W.HT_USED_SC[np.abs(sc_f) <= usable_half_hz - 0.4e6]
    used = W._bins(good)
    Lh = np.zeros(64)
    Lh[W._bins(W.HT_SC)] = W._HTLTF
    Yl = [fft_at(ltf1 + 80 * i + 16) for i in range(nsts)]
    H = np.zeros((nsts, 64), complex)
    for i in range(nsts):
        acc = sum(HT_LTF_P[i, n] * Yl[n] for n in range(nsts)) if nsts > 1 else Yl[0]
        H[i, used] = acc[used] / (nsts * Lh[used])
    # Per-subcarrier noise shape from the L-LTF pair; the four HT-only carriers take the median.
    nk_leg = rx.get("nk")
    base = float(np.median(nk_leg[rx["used"]])) if nk_leg is not None else 1.0
    nk = np.full(64, base)
    if nk_leg is not None:
        have = nk_leg[used] > 0
        nk[used] = np.where(have, nk_leg[used], base)

    data_bins = W._bins(W.HT_DATA_SC)
    live_data = np.isin(data_bins, used)
    pb = W._bins(W.HT_PILOT_SC)
    psc = W.HT_PILOT_SC.astype(float)
    kk = np.where(np.arange(64) < 32, np.arange(64), np.arange(64) - 64)
    Ys = [fft_at(data0 + sym * k + gi) for k in range(n_sym)]

    def tracked(Hs):
        """Pilot tracking: common phase + slope per symbol, the slope smoothed (it drifts slowly).
        Returns the corrected symbols, their pilot predictions and the mean pilot residual."""
        pred = [(ht_pilots(k, nsts) * Hs[:, pb]).sum(0) for k in range(n_sym)]
        slope = np.array([_track(Ys[k][pb], pred[k], psc)[1] for k in range(n_sym)])
        if n_sym >= 5:
            w = min(9, n_sym | 1)
            slope = np.convolve(np.pad(slope, w // 2, mode="edge"), np.ones(w) / w, "valid")
        Zc, resid = [], []
        for k in range(n_sym):
            Yk = Ys[k] * np.exp(-1j * slope[k] * kk)
            Yk = Yk * np.exp(-1j * np.angle(np.sum(Yk[pb] * np.conj(pred[k]))))
            resid.append(np.abs(Yk[pb] - pred[k]) ** 2)
            Zc.append(Yk)
        return Zc, float(np.mean(resid)) if resid else base

    # Noise shape across subcarriers: the L-LTF pair gives one noisy sample per subcarrier, so it is
    # smoothed over 9 and clipped to [0.5, 4] -- enough to de-weight a narrowband interferer without
    # letting estimation noise scatter the LLR weights (the raw shape lost CRC-valid real MCS7 frames).
    nshape = np.convolve(np.pad(nk, 4, mode="wrap"), np.ones(9) / 9, "valid")
    nshape = np.clip(nshape / max(np.median(nshape[used]), 1e-20), 0.5, 4.0)

    def combine(Hs, Zc, lvl):
        """Equalise (one stream) or Alamouti-combine (STBC) -> per data symbol: points, noise variance."""
        noise = lvl * nshape
        pts, var = [], []
        if nsts == 1:
            h = np.where(live_data, Hs[0, data_bins], 1)
            g = np.where(live_data, np.abs(h) ** 2, np.inf)
            for k in range(n_sym):
                pts.append(np.where(live_data, Zc[k][data_bins] / h, 0))
                var.append(noise[data_bins] / g)
        else:
            h1, h2 = Hs[0, data_bins], Hs[1, data_bins]
            g = np.where(live_data, np.abs(h1) ** 2 + np.abs(h2) ** 2, np.inf)
            gs = np.where(live_data, g, 1)
            for m in range(n_sym // 2):
                ya, yb = Zc[2 * m][data_bins], Zc[2 * m + 1][data_bins]
                for s_ in ((np.conj(h1) * ya + h2 * np.conj(yb)) / gs, (np.conj(h1) * yb - h2 * np.conj(ya)) / gs):
                    pts.append(np.where(live_data, s_, 0))
                    var.append(noise[data_bins] / g)
        return pts, var

    # Channel smoothing, unless the transmitter says not to (HT-SIG smoothing bit 0: beamformed).
    # Each width is a candidate receiver; for BCC the one whose estimate best predicts the pilots wins.
    cands = []
    for w in ((0, 1, 2, 3) if smoothing else (0,)):
        Hs = np.array([smooth_channel(H[i], good, w) for i in range(nsts)])
        Zc, lvl = tracked(Hs)
        cands.append((lvl, w, Hs, Zc))
    cands.sort(key=lambda c: c[0])

    out = {"n_sym": n_sym, "aggregated": bool(aggregation), "subframes": [], "coding": "LDPC" if ldpc else "BCC",
           "stbc": bool(stbc)}
    if ldpc:
        # Every candidate receiver decodes; a codeword is taken from the first that satisfies its parity
        # checks (a false pass of 162..972 checks is negligible). Codewords that pass are then known
        # exactly, so their symbols re-train the channel over the whole PPDU (data-aided estimation),
        # and the rest are decoded again -- at the edge of MCS7 this recovers frames one HT-LTF cannot.
        offs = LD.offsets(LP)
        known, firsts = {}, []
        trials = [(lvl, Hs, Zc) for lvl, _, Hs, Zc in cands]
        for rnd in range(8):
            new = 0
            for lvl, Hs, Zc in trials:
                pts, var = combine(Hs, Zc, lvl)
                llr = np.concatenate([np.where(live_data[:, None], llr_exact(p_, nbpsc, v_), 0.0).reshape(-1)
                                      for p_, v_ in zip(pts, var)])
                _, cw_ok, _, hard = LD.decode_stream(llr, LP, coderate, full=True)
                firsts.append((llr[:LP["K"] - LD._split(LP["nshrt"], LP["ncw"])[0]] < 0).astype(np.uint8))
                for i in np.flatnonzero(cw_ok):
                    if i not in known:
                        known[i], new = hard[i], new + 1
                if len(known) == LP["ncw"]:
                    break
            if len(known) == LP["ncw"] or not new or rnd == 7:
                break
            Hda = _data_aided_channel(known, LP, offs, ncbps, nbpsc, n_sym, nsts, cands[0][3],
                                      cands[0][2], data_bins, pb, Lh, Yl, used)
            Zc2, lvl2 = tracked(Hda)
            trials = [(lvl2, Hda, Zc2)]
        lay = LD._layout(LP)
        info = []
        for i, (k_i, _, _) in enumerate(lay):
            if i in known:
                info.append(known[i][:k_i])
            elif i == 0:
                # undecoded first codeword: its systematic bits as received, from the receiver variant whose
                # SERVICE field descrambles cleanly, so the scrambler seed stays usable for attribution
                info.append(next((f for f in firsts if not descramble(f[:16])[0][7:16].any()), firsts[0]))
            else:                                   # undecoded codeword: its systematic bits, as received
                lo = offs[i]
                info.append((llr[lo:lo + k_i] < 0).astype(np.uint8))
        out["n_codewords"], out["n_codewords_ok"] = int(LP["ncw"]), len(known)
        databits = np.concatenate(info)
    else:
        lvl, _, Hs, Zc = cands[0]
        pts, var = combine(Hs, Zc, lvl)
        _, inv = deinterleave_indices(ncbps, nbpsc, ncol=13)
        all_soft = []
        for p_, v_ in zip(pts, var):
            llr = -llr_exact(p_, nbpsc, v_)                          # Viterbi wants +ve = bit 1
            llr[~live_data] = 0.0
            coded = np.zeros(ncbps, llr.dtype)
            coded[inv] = llr.reshape(-1)
            all_soft.append(coded)
        soft = np.concatenate(all_soft) if all_soft else np.zeros(0)
        databits = W.viterbi(depuncture(soft, coderate))
    descr, seed = descramble(databits)
    out["scrambler_seed"] = int(sum(int(b) << i for i, b in enumerate(seed)))
    out["service_ok"] = bool(len(descr) >= 16 and not descr[7:16].any())
    payload = descr[16:16 + 8 * length]                                 # drop the 16 SERVICE bits
    if len(payload) < 8 * length:
        payload = np.pad(payload, (0, 8 * length - len(payload)))
    data = bytes(np.packbits(payload.reshape(-1, 8)[:, ::-1]))

    parts = deaggregate(data) if aggregation else [(data, True)]
    for mpdu, delim_ok in parts:
        if len(mpdu) < 4:
            continue
        body, fcs = mpdu[:-4], mpdu[-4:]
        ok = crc32_802(body) == int.from_bytes(fcs, "little")
        sub = {"fcs_ok": ok, "length_bytes": len(mpdu), "mpdu_hex": body.hex()}
        if aggregation:
            sub["delimiter_crc_ok"] = delim_ok
        if ok:
            sub["mac"] = parse_mac(body)
        out["subframes"].append(sub)
    out["n_subframes"] = len(out["subframes"])
    out["n_fcs_ok"] = sum(1 for s in out["subframes"] if s["fcs_ok"])
    return out


def _qam_points(bits, nbpsc):
    """Coded bits (n, nbpsc) -> unit-power constellation points (inverse of llr_exact's labelling)."""
    k = KMOD[nbpsc]
    if nbpsc == 1:
        return (2 * bits[:, 0] - 1) * k + 0j
    if nbpsc == 2:
        return ((2 * bits[:, 0] - 1) + 1j * (2 * bits[:, 1] - 1)) * k
    per = nbpsc // 2
    lv, lab = _PAM[per]
    code = {tuple(r): v for r, v in zip(lab, lv)}
    i = np.array([code[tuple(r)] for r in bits[:, :per]])
    q = np.array([code[tuple(r)] for r in bits[:, per:]])
    return (i + 1j * q) * k


def _data_aided_channel(known, LP, offs, ncbps, nbpsc, n_sym, nsts, Zc, Hs, data_bins, pb, Lh, Yl, used):
    """Least-squares channel per subcarrier and space-time stream from every data symbol whose coded
    bits are known (codewords that passed their parity checks), every pilot, and the HT-LTFs.
    STBC: a known pair (d0, d1) gives ya = h1 d0 - h2 d1*, yb = h1 d1 + h2 d0*."""
    from . import ldpc as LD
    nb = n_sym * ncbps
    bits = np.zeros(nb, np.uint8)
    have = np.zeros(nb, bool)
    for i, cw in known.items():
        tx = LD.transmitted(cw, LP, i)
        lo = offs[i]
        hi = min(lo + len(tx), nb)
        bits[lo:hi], have[lo:hi] = tx[:hi - lo], True
    nsc = ncbps // nbpsc
    bits = bits.reshape(n_sym, nsc, nbpsc)
    have = have.reshape(n_sym, nsc, nbpsc).all(-1)
    d = np.zeros((n_sym, nsc), complex)
    d[have] = _qam_points(bits[have], nbpsc)
    Y = np.array([z[data_bins] for z in Zc])                           # (n_sym, 52) tracked symbols
    out = Hs.copy()
    # normal equations G h = c per subcarrier; the HT-LTFs enter as rows too (a prior worth one symbol)
    G = np.zeros((nsc, nsts, nsts), complex)
    c = np.zeros((nsc, nsts), complex)
    Lk = Lh[data_bins]
    for n in range(nsts):
        a = np.array([HT_LTF_P[i, n] if nsts > 1 else 1 for i in range(nsts)], float)
        G += np.outer(a, a)[None]
        c += (Yl[n][data_bins] / np.where(Lk != 0, Lk, 1))[:, None] * a[None, :]
    if nsts == 1:
        G[:, 0, 0] += np.sum(np.abs(d) ** 2 * have, 0)
        c[:, 0] += np.sum(Y * np.conj(d) * have, 0)
    else:
        for m in range(n_sym // 2):
            ok = have[2 * m] & have[2 * m + 1]
            d0, d1 = d[2 * m], d[2 * m + 1]
            for y, a1, a2 in ((Y[2 * m], d0, -np.conj(d1)), (Y[2 * m + 1], d1, np.conj(d0))):
                A = np.stack([a1, a2], -1) * ok[:, None]
                G += np.conj(A)[:, :, None] * A[:, None, :]
                c += np.conj(A) * (y * ok)[:, None]
    sol = np.linalg.solve(G + 1e-9 * np.eye(nsts)[None], c[..., None])[..., 0]   # (nsc, nsts)
    live = np.isin(data_bins, used)
    for i in range(nsts):
        out[i, data_bins[live]] = sol[live, i]
    # pilots: known on every symbol
    Yp = np.array([z[pb] for z in Zc])
    for j, b in enumerate(pb):
        rows = [ht_pilots(k, nsts)[:, j] for k in range(n_sym)]
        vals = list(Yp[:, j])
        for n in range(nsts):
            rows.append(np.array([HT_LTF_P[i, n] if nsts > 1 else 1 for i in range(nsts)], float))
            vals.append(Yl[n][b] / Lh[b])
        out[:, b] = np.linalg.lstsq(np.array(rows, float), np.array(vals), rcond=None)[0]
    return out


# ---------------------------------------------------------------- driver
def run(meta_path, data_path=None, out="wifi_decoded", scan=False):
    x, fs, fc, meta = W.load_sigmf(meta_path, data_path)
    det = W.WifiDetector(fs, fc, scan_step_hz=1e6 if scan else None)
    # We need per-packet rx state, so re-run detection but keep _rx.
    results = []
    n_cap = len(x)
    for chan in det.channels:
        x20 = det.baseband(x, chan["offset_hz"])
        for d in det.detect_channel(x20, chan):
            rec = {"channel": chan["name"], "center_mhz": round((chan["freq_hz"] + d["cfo"]) / 1e6, 3),
                   "start_us": round(d["_rx"]["start"] / FS * 1e6, 3), "fmt": d["fmt"],
                   "rate_mbps": d["mbps"], "length": d["length"], "snr_db": round(float(d["snr"]), 1)}
            if d["fmt"] == "NON_HT":
                try:
                    mpdu, ok, seed, info, nsym = decode_nonht_payload(
                        d["_rx"], d["mbps"], d["length"], chan["offset_hz"],
                        chan.get("usable_half_hz", 1e12))
                    rec["fcs_ok"] = ok
                    rec["mpdu_hex"] = bytes(bytearray(int(v) for v in mpdu)).hex()
                    rec.update({k: v for k, v in info.items()})
                except Exception as e:
                    rec["decode_error"] = str(e)
            else:
                rec["note"] = f"{d['fmt']} {d['bw']}MHz " + \
                    f"{d['extra'].get('nsts_su','?')}ss -- payload needs multi-antenna/LDPC RX"
            results.append(rec)
    # de-duplicate wide PPDUs seen on multiple sub-channels (same start, close centre)
    results.sort(key=lambda r: (r["start_us"], r["center_mhz"]))
    dedup = []
    for r in results:
        if any(abs(r["start_us"] - k["start_us"]) < 1.0 and abs(r["center_mhz"] - k["center_mhz"]) >= 20
               for k in dedup):
            continue
        dedup.append(r)
    def clean(o):
        import numpy as _np
        if isinstance(o, dict): return {k: clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)): return [clean(v) for v in o]
        if isinstance(o, (_np.integer,)): return int(o)
        if isinstance(o, (_np.floating,)): return float(o)
        if isinstance(o, (_np.bool_,)): return bool(o)
        return o
    dedup = clean(dedup)
    json.dump(dedup, open(out + ".json", "w"), indent=1)
    return dedup


def write_csv(res, path):
    import csv
    cols = ["start_us", "channel", "center_mhz", "fmt", "rate_mbps", "length", "snr_db",
            "fcs_ok", "type", "subtype", "addr1_ra", "addr2_ta", "addr1", "sa", "da",
            "bssid", "ssid", "retry", "duration_us", "mpdu_hex", "note"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in res:
            w.writerow(r)


def summary(res):
    from collections import Counter
    nh = [r for r in res if r["fmt"] == "NON_HT"]
    ok = [r for r in nh if r.get("fcs_ok")]
    print(f"\n{len(res)} PPDUs | {len(nh)} non-HT decodable | {len(ok)} FCS-verified\n")
    if not ok:
        return
    print("Fully decoded, FCS-verified frames:")
    print("-" * 78)
    for r in sorted(ok, key=lambda r: r["start_us"]):
        ra = r.get("addr1_ra") or r.get("addr1") or "-"
        rav = r.get("addr1_ra_vendor") or r.get("addr1_vendor") or ""
        line = f"  {r['start_us']:8.1f}us  {r.get('subtype',''):11s} {r['center_mhz']:7.1f}MHz"
        line += f"  RA {ra}"
        if rav and "unknown" not in rav:
            line += f" [{rav}]"
        if r.get("addr2_ta"):
            tv = r.get("addr2_ta_vendor", "")
            line += f"  TA {r['addr2_ta']}" + (f" [{tv}]" if tv and 'unknown' not in tv else "")
        if r.get("subtype") == "beacon":
            line += f"\n              SSID={r.get('ssid')!r}  BSSID={r.get('bssid')}"
            bv = r.get("bssid_vendor", "")
            if bv:
                line += f" [{bv}]"
        print(line)
    print("-" * 78)
    # device inventory
    devs = {}
    for r in ok:
        for k in ("addr1_ra", "addr1", "addr2_ta", "sa", "bssid"):
            m = r.get(k)
            if m and m != "ff:ff:ff:ff:ff:ff":
                v = r.get(k + "_vendor") or vendor(m)
                if devs.get(m) in (None, "", "unknown OUI"):
                    devs[m] = v
    print(f"\n{len(devs)} distinct devices:")
    for m, v in sorted(devs.items(), key=lambda x: (x[1] or "zz", x[0])):
        print(f"  {m}  {v}")
    ssids = [r["ssid"] for r in ok if r.get("subtype") == "beacon" and r.get("ssid") and r["ssid"] != "<broadcast>"]
    if ssids:
        print(f"\nNetwork names (SSIDs) seen: {', '.join(repr(s) for s in set(ssids))}")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Decode 802.11 non-HT payloads and MAC frames")
    ap.add_argument("meta"); ap.add_argument("--data"); ap.add_argument("--out", default="wifi_decoded")
    ap.add_argument("--scan", action="store_true")
    a = ap.parse_args()
    res = run(a.meta, a.data, a.out, a.scan)
    write_csv(res, a.out + ".csv")
    summary(res)
