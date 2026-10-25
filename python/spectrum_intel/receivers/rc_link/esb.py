"""2 Mb/s GFSK packet links of the Enhanced ShockBurst family (Nordic nRF24 / nRF52 ESB and its clones: 2.4 GHz
wireless mice, keyboards, game controllers, toy RC links).

On air: an alternating 0101 preamble, an address / sync word, header + payload, CRC-16/CCITT (poly 0x1021,
MSB first). The CRC init, the CRC span and whatever the vendor puts in front of the CRC (address bytes, scrambling
with a fixed pattern) differ per product, so a packet is not checked alone. The packets of one link (same sync
word) are checked together instead: CRC(m, I) = CRC0(m) ^ T_E(I) for a message m of E bits (CRC0 = init 0, T_E =
E zero-input register steps), so K = CRC0(m[:E]) ^ m[E:E+16] equals T_E(I) for every valid packet. The CRC end E
is the one where the most packets of distinct content share one K; I = T_E^-1(K). A wrong E gives every packet
an independent random K (2^-16 each), except next to the true span: at E_true - 1, K takes two values that differ by
0x8000 | poly >> 1 = 0x8810 (the feedback bit of the missing step), so about half the packets agree there and more
when the CRC's last bit is unreliable; such an E is skipped. A span needs at least half its group's packets; ties go
to the shortest (a fixed bit after the CRC makes E_true + 1 constant too: the bench link).

Bench 2026-10-17 (testdata/live/LIVE_ESB2M_2440MHz): a 1 kHz link on 2426 / 2450 MHz, two nodes 145 us apart,
94 us packets, sync 0x31A9D9.. / 0x31A9D8.., CRC end at bit 124, init 0xDD19 at the sync word: 86 of 88 packets
pass (the two failures are bursts cut by Wi-Fi). Depends on numpy / scipy only.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

RATE = 2e6                 # bit/s
HALF_BW = 1.25e6           # slice half-width: GFSK 2 Mb/s, h ~0.5 (99 % bandwidth ~2 MHz) + carrier offset
OUT_FS = 16e6              # 8 samples per bit
MIN_PREAMBLE = 14          # alternating bits: a 16-bit preamble (nRF52 2 Mbit) measures 14-15 (ramp-up, GFSK ISI)
PREAMBLE_START = 32        # ... starting within the first 32 bits (the bench link sends ~20 lead-in bits):
                           # noise passes with p ~32 x 2^-13
SYNC_BITS = 16             # packets of one link share their first 16 bits after the preamble
CRC_POLY = 0x1021
MIN_SPAN = 24              # CRC span >= address
MIN_VERIFIED = 3           # distinct packets on one (E, K): chance ~ n_E * C(n, 3) / 2^32
SHORT_BY_ONE = 0x8000 | (CRC_POLY >> 1)   # K(E_true - 1) classes differ by this


def _two_means(v, it=20):
    a, b = np.percentile(v, 25), np.percentile(v, 75)
    for _ in range(it):
        m = np.abs(v - a) < np.abs(v - b)
        if m.all() or (~m).all():
            break
        a, b = v[m].mean(), v[~m].mean()
    return min(a, b), max(a, b)


def slice_packet(y, ofs):
    """Bits of one burst (y: complex baseband at ofs, centred within +-0.3 MHz of the carrier) from the sync word
    (the alternating preamble's last bit) to the end of the burst. Returns dict(bits, preamble_bits,
    dev_hz, f_off_hz: the carrier relative to y's centre, from the discriminator's two levels) or None when there
    is no >= MIN_PREAMBLE-bit preamble starting in the first PREAMBLE_START bits."""
    sps = ofs / RATE
    if len(y) < 48 * sps:
        return None
    ifr = np.angle(y[1:] * np.conj(y[:-1]))
    pw = np.abs(y[1:]) ** 2
    on = np.flatnonzero(pw > 0.25 * np.percentile(pw, 90))
    if len(on) < 48 * sps:
        return None
    a, b = int(on[0]), min(len(ifr) - 1, int(on[-1] + 4 * sps))   # the CRC's last bits may ride the ramp-down
    ifs = ndimage.uniform_filter1d(ifr, max(1, int(sps / 2)))
    lo, hi = _two_means(ifs[a:b + 1])
    thr = (lo + hi) / 2
    best = None
    for ph in range(int(np.ceil(sps))):                    # bit timing: the phase with the widest eye
        idx = (a + ph + np.arange(int((b - a - ph) / sps)) * sps).astype(int)
        v = ifs[idx] - thr
        sc = float(np.mean(np.abs(v)))
        if best is None or sc > best[0]:
            best = (sc, v)
    v = best[1]
    bits = (v > 0).astype(np.uint8)
    alt = bits[1:] != bits[:-1]                             # alt[j]: bits j and j + 1 differ
    start = None
    for s0 in range(min(PREAMBLE_START, len(alt))):
        if s0 > 0 and alt[s0 - 1]:
            continue                                        # inside a run, not its first bit
        k = s0
        while k < len(alt) and alt[k]:
            k += 1
        if k - s0 + 1 >= MIN_PREAMBLE:                      # bits s0..k alternate
            start = k                                       # the run's last bit begins the sync word
            break
    if start is None or len(bits) - start < MIN_SPAN + 16:
        return None
    return dict(bits=bits[start:], preamble_bits=start, dev_hz=float((hi - lo) / 2 * ofs / (2 * np.pi)),
                f_off_hz=float(thr * ofs / (2 * np.pi)))


def _crc0_prefix(bits):
    """out[k] = CRC-16/CCITT register (init 0) after bits[:k]."""
    out = np.zeros(len(bits) + 1, np.int64)
    c = 0
    for k, bt in enumerate(bits):
        fb = ((c >> 15) & 1) ^ int(bt)
        c = ((c << 1) & 0xFFFF) ^ (CRC_POLY if fb else 0)
        out[k + 1] = c
    return out


def _field16(bits):
    """out[k] = the 16 bits starting at k, MSB first (k <= len - 16)."""
    n = len(bits) - 15
    out = np.zeros(max(n, 0), np.int64)
    for j in range(16):
        out += bits[j:j + n].astype(np.int64) << (15 - j)
    return out


def _unstep(c, n):
    """n inverse zero-input register steps (the poly's x^0 term makes a step invertible)."""
    for _ in range(n):
        c = (((c ^ CRC_POLY) >> 1) | 0x8000) if c & 1 else c >> 1
    return c


def _hex(bits):
    n = (len(bits) + 7) // 8 * 8
    p = np.zeros(n, np.uint8)
    p[:len(bits)] = bits
    return bytes(np.packbits(p)).hex().upper()


def verify_links(frames, starts_us):
    """frames: list of bit arrays (slice_packet()["bits"]), starts_us their burst starts. Returns one dict per
    frame: crc_ok, sync_hex (first 32 bits: address-like, one per node), the node's frame period and, for a link
    whose CRC was found, the CRC span / init and the frame / CRC fields."""
    out = [dict(crc_ok=False, sync_hex=_hex(f[:32])) for f in frames]
    nodes = {}
    for i, o in enumerate(out):
        nodes.setdefault(o["sync_hex"], []).append(starts_us[i])
    for o in out:
        T = frame_period_us(nodes[o["sync_hex"]])
        if T is not None:
            o["node_frame_period_us"] = T
    groups = {}
    for i, f in enumerate(frames):
        groups.setdefault(bytes(f[:SYNC_BITS]), []).append(i)
    for idx in groups.values():
        if len(idx) < MIN_VERIFIED:
            continue
        pre = {i: _crc0_prefix(frames[i]) for i in idx}
        fld = {i: _field16(frames[i]) for i in idx}
        n_max = max(len(frames[i]) for i in idx)
        best = None                                         # (n_distinct, E, K)
        for E in range(MIN_SPAN, n_max - 15):
            seen = {}
            for i in idx:
                f = frames[i]
                if len(f) < E + 16:
                    continue
                K = int(pre[i][E] ^ fld[i][E])
                seen.setdefault(K, set()).add(bytes(f[:E + 16]))
            if not seen:
                continue
            K, s = max(seen.items(), key=lambda kv: (len(kv[1]), -kv[0]))
            if 4 * len(seen.get(K ^ SHORT_BY_ONE, ())) >= len(s):
                continue                                    # one bit short of the true span
            if len(s) >= MIN_VERIFIED and 2 * len(s) >= len(idx) and (best is None or len(s) > best[0]):
                best = (len(s), E, K)
        if best is None:
            continue
        _, E, K = best
        init = _unstep(K, E)
        n_ok = sum(1 for i in idx if len(frames[i]) >= E + 16 and int(pre[i][E] ^ fld[i][E]) == K)
        for i in idx:
            f = frames[i]
            ok = len(f) >= E + 16 and int(pre[i][E] ^ fld[i][E]) == K
            out[i].update(crc_ok=ok, crc="CRC-16/CCITT (0x1021)", crc_span_bits=E, crc_init_hex=f"{init:04X}",
                          link_packets_crc_ok=n_ok)
            if ok:
                out[i].update(frame_hex=_hex(f[:E]), crc_hex=f"{int(fld[i][E]):04X}")
    return out


def frame_period_us(starts_us):
    """Median spacing of a link's packets, folded onto the shortest repeat: gaps of k x T (missed packets) count
    as T. None for < 3 packets."""
    t = np.sort(np.asarray(starts_us, float))
    if len(t) < 3:
        return None
    d = np.diff(t)
    d = d[d > 20.0]
    if not len(d):
        return None
    T = float(np.percentile(d, 10))
    k = np.maximum(np.round(d / T), 1)
    return float(np.median(d / k))
