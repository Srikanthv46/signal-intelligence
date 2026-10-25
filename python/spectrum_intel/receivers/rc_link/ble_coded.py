"""
ble_coded.py -- BLE LE Coded PHY (Long Range, S=8 / S=2) for Remote ID over Bluetooth 5 (research, PLAN T2b;
research/drone_transmitters/09_ble_coded_droneid_legacy.md; bench bench/ble_coded). Synthetic validation only.

Packet (Core 5.x Vol 6 Part B 2.2 / 3.3), 1 Msym/s GFSK h 0.5:
  preamble 80 symbols (00111100 x 10, uncoded) | FEC block 1 (S=8): AA 32 b + CI 2 b + TERM1 3 b |
  FEC block 2 (S=8 if CI = 00, S=2 if CI = 01): PDU + CRC-24 + TERM2 3 b.
FEC: rate 1/2, K = 4, G0 = 1+x+x^2+x^3 (octal 17), G1 = 1+x^2+x^3 (octal 13), a0 first, state 0, 3-bit tail.
S=8 pattern mapper: 0 -> 0011, 1 -> 1100 (S=2: identity). Whitening and CRC as on the 1M PHY, before the FEC.
Receiver: the coded preamble + coded access address (376 known symbols) by soft differential correlation, the
coherent GFSK detector (gfsk_coherent) for soft symbols, S=8 / S=2 de-mapping and a soft K=4 Viterbi.
Bench: a Remote ID advert decodes 10/10 at 0 dB in-band SNR (1 MHz), ~7 dB below the 1M PHY.
"""
from __future__ import annotations

import numpy as np

from . import gfsk_coherent as GC

AA_ADV = 0x8E89BED6
PREAMBLE = [0, 0, 1, 1, 1, 1, 0, 0] * 10
PAT = {0: [0, 0, 1, 1], 1: [1, 1, 0, 0]}


def conv_encode(bits):
    s0 = s1 = s2 = 0
    out = []
    for b in bits:
        out += [b ^ s0 ^ s1 ^ s2, b ^ s1 ^ s2]
        s0, s1, s2 = b, s0, s1
    return out


def map_s8(coded):
    return [x for c in coded for x in PAT[c]]


def coded_symbols(aa_bits, ci, pdu_crc_whitened):
    """Air symbols of one Coded PHY packet (preamble, block 1 at S=8, block 2 at S=8 (ci 0) or S=2 (ci 1))."""
    b1 = map_s8(conv_encode(list(aa_bits) + [ci & 1, (ci >> 1) & 1] + [0, 0, 0]))
    c2 = conv_encode(list(pdu_crc_whitened) + [0, 0, 0])
    return PREAMBLE + b1 + (map_s8(c2) if ci == 0 else c2)


def _transitions():
    """For each next state ns and the dropped bit s2: (predecessor state, coded bit a0, a1, input bit)."""
    pre, a0, a1 = np.zeros((8, 2), int), np.zeros((8, 2), int), np.zeros((8, 2), int)
    for ns in range(8):
        b, s0, s1 = ns >> 2, (ns >> 1) & 1, ns & 1
        for s2 in (0, 1):
            pre[ns, s2] = (s0 << 2) | (s1 << 1) | s2
            a0[ns, s2] = b ^ s0 ^ s1 ^ s2
            a1[ns, s2] = b ^ s1 ^ s2
    return pre, 2 * a0 - 1, 2 * a1 - 1


_PRE, _A0, _A1 = _transitions()


def viterbi(llr):
    """Soft Viterbi (K=4, 8 states, start in 0, terminated): llr (n, 2), positive -> 1. Returns n decoded bits."""
    n = len(llr)
    pm = np.full(8, 1e18); pm[0] = 0.0
    surv = np.zeros((n, 8), np.int8)
    for t in range(n):
        m = pm[_PRE] - _A0 * llr[t, 0] - _A1 * llr[t, 1]       # (8 next states, 2 predecessors)
        k = np.argmin(m, axis=1)
        surv[t] = k
        pm = m[np.arange(8), k]
    st = 0
    out = np.zeros(n, np.uint8)
    for t in range(n - 1, -1, -1):
        out[t] = st >> 2
        st = _PRE[st, surv[t, st]]
    return out


def known_symbols(aa_bits, ci=0):
    return np.array(PREAMBLE + map_s8(conv_encode(list(aa_bits) + [ci & 1, (ci >> 1) & 1] + [0, 0, 0])), np.uint8)


def find_coded(p, aa_bits, top=2, min_ncorr=0.4):
    """Start samples of coded packets: soft correlation of the one-bit differential phase with the 376 known symbols
    (CI = 00 and 01 differ only in block 1's tail; the S=8 template is used for both)."""
    return GC.find_access_address(p, known_symbols(aa_bits, 0), top=top, min_ncorr=min_ncorr)


def _soft(bits, rel):
    return (2.0 * bits - 1) * np.minimum(rel, np.r_[rel[0], rel[:-1]])      # alpha_k from b_k and b_(k-1)


def decode(y, s0, aa_bits, max_pdu_octets=257):
    """Coded packet whose preamble starts at sample s0 of y (8 MSps) -> dict(aa_ok, ci, block-2 bits (PDU + CRC,
    still whitened), h) or None."""
    n1 = 80 + 4 * 2 * 37
    nsym = n1 + 4 * 2 * (8 * max_pdu_octets + 24 + 3)
    dm = GC.demodulate(y, s0, nsym)
    if dm is None or len(dm[0]) < n1 + 2 * 43:
        return None
    soft = _soft(dm[0], dm[1])
    p8 = np.array([1.0, 1.0, -1.0, -1.0])                                    # +: pattern of a 1
    llr1 = (soft[80:n1].reshape(74, 4) @ p8).reshape(37, 2)
    blk1 = viterbi(llr1)
    aa = blk1[:32]
    ci = int(blk1[32]) | (int(blk1[33]) << 1)
    aa_err = int(np.count_nonzero(aa != np.asarray(aa_bits)))
    rest = soft[n1:]
    if ci == 0:
        nb = len(rest) // 8
        llr2 = (rest[:8 * nb].reshape(2 * nb, 4) @ p8).reshape(nb, 2)
    elif ci == 1:
        nb = len(rest) // 2
        llr2 = rest[:2 * nb].reshape(nb, 2)
    else:
        return dict(aa_bit_errors=aa_err, ci=ci, h=dm[2])
    return dict(aa_bit_errors=aa_err, ci=ci, bits=viterbi(llr2), h=dm[2])


def parse_extended(by):
    """Common extended advertising payload (ADV_EXT_IND / AUX_ADV_IND, Core Vol 6 Part B 2.3.4): ext header length +
    AdvMode, flags, AdvA, TargetA, CTEInfo, ADI, AuxPtr, SyncInfo, TxPower, ACAD; the rest is AdvData."""
    if not by:
        return {}
    hl, mode = by[0] & 0x3F, by[0] >> 6
    out = dict(adv_mode=mode)
    if hl == 0:
        out["adv_data"] = bytes(by[1:])
        return out
    flags = by[1] if len(by) > 1 else 0
    i, end = 2, 1 + hl
    if flags & 0x01 and i + 6 <= end:
        out["adv_addr"] = ":".join(f"{v:02X}" for v in by[i:i + 6][::-1]); i += 6
    if flags & 0x02 and i + 6 <= end:
        out["target_addr"] = ":".join(f"{v:02X}" for v in by[i:i + 6][::-1]); i += 6
    if flags & 0x04 and i + 1 <= end:
        i += 1
    if flags & 0x08 and i + 2 <= end:
        v = by[i] | (by[i + 1] << 8)
        out["adi"] = dict(did=v & 0xFFF, sid=v >> 12); i += 2
    if flags & 0x10 and i + 3 <= end:
        v = by[i] | (by[i + 1] << 8) | (by[i + 2] << 16)
        units = (v >> 7) & 1
        out["aux_ptr"] = dict(channel=v & 0x3F, ca=(v >> 6) & 1, offset_us=((v >> 8) & 0x1FFF) * (300 if units else 30),
                              phy={0: "1M", 1: "2M", 2: "Coded"}.get((v >> 21) & 7, "?"))
        i += 3
    if flags & 0x20 and i + 18 <= end:
        i += 18
    if flags & 0x40 and i + 1 <= end:
        t = by[i]
        out["tx_power_dbm"] = t - 256 if t > 127 else t; i += 1
    out["adv_data"] = bytes(by[end:])
    return out
