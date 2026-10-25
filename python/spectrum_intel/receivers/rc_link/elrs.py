"""
elrs.py -- ExpressLRS over-the-air packets (research, PLAN R2c; research/drone_transmitters/07_elrs_details.md, from
ExpressLRS OTA.h / OTA.cpp / crc.cpp, GPL-3, read as reference only). Synthetic validation only.

OTA4 (8 bytes, sub-GHz LoRa / 2.4 GHz modes without "8ch"): byte 0 = type (bits 1:0) + crcHigh (bits 7:2), bytes 1-6
payload, byte 7 = crcLow; CRC-14 (0x2E57, MSB-first, non-reflected, no final XOR) over bytes 0-6 with crcHigh zeroed.
OTA8 (13 bytes): type in byte 0 bits 1:0, CRC-16 (0x3D65) over bytes 0-10, little-endian in bytes 11-12.
CRC init = ((UID4 << 8) | UID5) ^ (OTA_VERSION << 8), XOR the nonce except for SYNC packets; binding: OTA_VERSION
(4.x) / 0 (3.x). A SYNC packet carries UID4 and UID5 itself (3.x also UID3), so it verifies itself: the CRC computed
with its own UID bytes. Model match (UID5 ^= ~modelId & 0x3F) is not searched (64 x the false-pass rate of a 14-bit
CRC); such SYNC packets stay unverified.
SYNC body: 4.x: fhssIndex, nonce, rfRateEnum, {switchEncMode:1, newTlmRatio:3, geminiMode:1, otaProtocol:2}, UID4,
UID5; 3.x: fhssIndex, nonce, {switchEncMode:1, newTlmRatio:3, rateIndex:4}, UID3, UID4, UID5.
"""
from __future__ import annotations

import numpy as np

PACKET_TYPES = {0: "RCDATA / LINKSTATS", 1: "DATA", 2: "SYNC", 3: "?"}
CRC14_POLY, CRC16_POLY = 0x2E57, 0x3D65


def _table(bits, poly):
    t = []
    for i in range(256):
        c = i << (bits - 8)
        for _ in range(8):
            c = ((c << 1) ^ poly) if c & (1 << (bits - 1)) else (c << 1)
        t.append(c & ((1 << bits) - 1))
    return t


_T14, _T16 = _table(14, CRC14_POLY), _table(16, CRC16_POLY)


def crc(data, init, bits):
    t = _T14 if bits == 14 else _T16
    c = init
    for b in data:
        c = (c << 8) ^ t[((c >> (bits - 8)) ^ b) & 0xFF]
    return c & ((1 << bits) - 1)


def _sync_fields(body, ver):
    """OTA_Sync_s (6 bytes) for OTA version 4 (4.x) or 3 (3.x)."""
    f = dict(fhss_index=body[0], nonce=body[1])
    if ver >= 4:
        b = body[3]
        f.update(rate_enum=body[2], switch_enc_mode=b & 1, tlm_ratio_code=(b >> 1) & 7, gemini_mode=(b >> 4) & 1,
                 ota_protocol=(b >> 5) & 3, uid4=body[4], uid5=body[5])
    else:
        b = body[2]
        f.update(switch_enc_mode=b & 1, tlm_ratio_code=(b >> 1) & 7, rate_index=b >> 4, uid3=body[3], uid4=body[4],
                 uid5=body[5])
    return f


def parse(payload: bytes):
    """OTA4 / OTA8 packet bytes -> dict (packet type; for SYNC the fields and whether its own UID bytes verify the
    CRC, which OTA version and whether in binding)."""
    p = bytes(payload)
    if len(p) not in (8, 13):
        return None
    typ = p[0] & 3
    out = dict(protocol="ExpressLRS", ota=("OTA4" if len(p) == 8 else "OTA8"), packet_type=PACKET_TYPES[typ],
               payload_hex=p.hex(), crc_ok=None)
    if typ != 2:
        return out
    if len(p) == 8:
        rx = ((p[0] >> 2) << 8) | p[7]
        data, bits = bytes([p[0] & 3]) + p[1:7], 14
    else:
        rx = p[11] | (p[12] << 8)
        data, bits = p[:11], 16
    body = p[1:7]
    for ver in (4, 3):
        f = _sync_fields(body, ver)
        inits = [(((f["uid4"] << 8) | f["uid5"]) ^ (ver << 8), False), (ver if ver >= 4 else 0, True)]
        for init, binding in inits:
            if crc(data, init, bits) == rx:
                out.update(crc_ok=True, ota_version=ver, binding=binding, sync=f)
                return out
    out["crc_ok"] = False
    out["note"] = "SYNC CRC not verified with its own UID bytes (model match on, other OTA version, or bit errors)"
    return out


def build_sync(fhss_index, nonce, rate, tlm_code, uid, ver=4, binding=False, ota8=False):
    """Encoder for tests: a SYNC packet as the TX builds it (uid = 6 UID bytes)."""
    if ver >= 4:
        body = bytes([fhss_index, nonce, rate, (tlm_code & 7) << 1, uid[4], uid[5]])
    else:
        body = bytes([fhss_index, nonce, ((tlm_code & 7) << 1) | ((rate & 15) << 4), uid[3], uid[4], uid[5]])
    init = (ver if ver >= 4 else 0) if binding else ((uid[4] << 8) | uid[5]) ^ (ver << 8)
    if not ota8:
        c = crc(bytes([2]) + body, init, 14)
        return bytes([2 | ((c >> 8) << 2)]) + body + bytes([c & 0xFF])
    pkt = bytes([2]) + body + bytes(4)
    c = crc(pkt, init, 16)
    return pkt + bytes([c & 0xFF, c >> 8])


# ---------------------------------------------------------------- FHSS hop sequence (research, PLAN R2d; file 07 §4)
# seed = UID[2..5] ^ OTA_VERSION in the low byte; LCG seed = (214013 seed + 2531011) mod 2^31, output seed >> 16;
# 256-entry sequence in blocks of freq_count: each block starts on the sync channel (freq_count // 2), the other
# entries are swapped with a random entry of the same block (rng % (count - 1) + 1). Facts from the published
# description, implemented independently.
FHSS_SEQUENCE_LEN = 256


def fhss_seed(uid2, uid3, uid4, uid5, ver=4):
    return (uid2 << 24) | (uid3 << 16) | (uid4 << 8) | (uid5 ^ ver)


def fhss_sequence(seed, count):
    """The 256-entry (or shorter, a multiple of count) channel-index sequence."""
    sync = count // 2
    n = (FHSS_SEQUENCE_LEN // count) * count
    seq = [sync if i % count == 0 else 0 if i % count == sync else i % count for i in range(n)]
    s = seed & 0xFFFFFFFF
    for i in range(n):
        if i % count:
            s = (214013 * s + 2531011) % (1 << 31)
            r = (s >> 16) % (count - 1) + 1
            j = (i // count) * count + r
            seq[i], seq[j] = seq[j], seq[i]
    return seq


def _lcg_jump(k):
    """(a, c) with LCG^k(s) = (a s + c) mod 2^31."""
    a, c, A, C = 1, 0, 214013, 2531011
    while k:
        if k & 1:
            a, c = (A * a) % (1 << 31), (A * c + C) % (1 << 31)
        A, C = (A * A) % (1 << 31), (A * C + C) % (1 << 31)
        k >>= 1
    return a, c


def uid23_candidates(uid4, uid5, count, observations, ver=4, uid3=None):
    """UID2 / UID3 pairs whose hop sequence puts every observed packet on its channel. observations: [(sequence
    index, channel index)], e.g. from a SYNC packet's fhssIndex plus the hops counted from its nonce. uid3 known
    (3.x SYNC) narrows the search to 256 seeds. Vectorised over the seeds; only the blocks holding an observation are
    simulated (the LCG jumps to each block's start)."""
    n_seq = (FHSS_SEQUENCE_LEN // count) * count
    u3 = np.arange(256) if uid3 is None else np.array([uid3])
    u2 = np.arange(128)                  # the LCG drops seed bit 31: UID2 bit 7 never changes the sequence
    U2, U3 = np.meshgrid(u2, u3, indexing="ij")
    seeds = ((U2.ravel().astype(np.int64) << 24) | (U3.ravel().astype(np.int64) << 16) | (uid4 << 8) | (uid5 ^ ver))
    ok = np.ones(len(seeds), bool)
    sync = count // 2
    blocks = {}
    for i, ch in observations:
        i %= n_seq
        blocks.setdefault(i // count, []).append((i % count, ch))
    base = np.array([sync if j == 0 else 0 if j == sync else j for j in range(count)], np.int16)
    for b, obs in blocks.items():
        a, c = _lcg_jump(b * (count - 1))
        st = (a * (seeds % (1 << 31)) + c) % (1 << 31)
        seq = np.tile(base, (len(seeds), 1))
        rows = np.arange(len(seeds))
        for j in range(1, count):
            st = (214013 * st + 2531011) % (1 << 31)
            r = (st >> 16) % (count - 1) + 1
            tmp = seq[rows, j].copy()
            seq[rows, j] = seq[rows, r]
            seq[rows, r] = tmp
        for j, ch in obs:
            ok &= seq[:, j] == ch
    return sorted((int(x >> 24) & 0xFF | hi, int(x >> 16) & 0xFF) for x in seeds[ok] for hi in (0, 0x80))


def recover_hops(sync, t_sync_us, sync_ch, interval_us, hop_every, count, packets, ver=None):
    """Hop-sequence recovery from one decoded SYNC packet and other packets of the link (research, PLAN R2d).
    sync: the parse() dict of a CRC-verified SYNC (fhss_index, nonce, uid4, uid5, uid3 on 3.x, ota_version), seen on
    channel index sync_ch; packets: [(t_us, channel index)] of the link's other packets. Each packet's nonce is counted from the SYNC by the
    packet interval; the sequence index advances when (nonce + 1) % hop_every == 0. Returns the UID2 / UID3
    candidates and, when they agree on the sequence, the next channels; None without an observation."""
    f = sync.get("sync") or sync
    ver = ver if ver is not None else sync.get("ota_version", 4)
    n0, i0 = f["nonce"], f["fhss_index"]
    seq_len = (FHSS_SEQUENCE_LEN // count) * count

    def index_at(dn):
        hops = sum(1 for n in range(n0, n0 + dn) if (n + 1) % hop_every == 0) if dn >= 0 else \
            -sum(1 for n in range(n0 + dn, n0) if (n + 1) % hop_every == 0)
        return (i0 + hops) % seq_len

    obs = [(i0, sync_ch)]
    for t, ch in packets:
        dn = round((t - t_sync_us) / interval_us)
        if dn != 0:
            obs.append((index_at(dn), ch))
    if len(obs) < 2:
        return None
    c = uid23_candidates(f["uid4"], f["uid5"], count, obs, ver, f.get("uid3"))
    out = dict(observations=len(obs), uid23_candidates=len(c), uid2_uid3=[f"{a:02X}{b:02X}" for a, b in c[:8]])
    seqs = {tuple(fhss_sequence(fhss_seed(a, b, f["uid4"], f["uid5"], ver), count)) for a, b in c}
    if len(seqs) == 1:
        s = seqs.pop()
        out["next_channels"] = [s[(i0 + k) % seq_len] for k in range(1, 9)]
    return out


# ---------------------------------------------------------------- LR1121 GFSK modes (research, PLAN G3; file 14 §3)
# 2.4 GHz FSK 1000 Hz / 500 Hz DVDA / 250 Hz DVDA (OTA4) and 900 MHz FSK 1000 Hz 8ch (OTA8): 300 kb/s, +-100 kHz,
# no pulse shaping, 16-bit 0x55 preamble, 16-bit sync word = UID5, UID4, fixed length, radio CRC off, whitening on
# (LR11xx default seed 0x1FF: the bit order is not published, so the CC2500-style PN9 bytes and their bit-reversed
# form are both tried). A SYNC packet verifies itself (its CRC with its own UID bytes) and its UID bytes must equal
# the sync word: the decode is reported only then.
GFSK_RATE, GFSK_DEV = 300e3, 100e3
GFSK_WHITENING = ("pn9_msb", "pn9_lsb", "none")


def _pn9_bytes(n):
    st, out = 0x1FF, []
    for _ in range(n):
        out.append(st & 0xFF)
        for _ in range(8):
            st = (st >> 1) | (((st ^ (st >> 5)) & 1) << 8)
    return out


_PN9 = _pn9_bytes(16)
_REV8 = [int(f"{v:08b}"[::-1], 2) for v in range(256)]


def gfsk_whiten(b: bytes, variant: str) -> bytes:
    if variant == "none":
        return bytes(b)
    w = _PN9 if variant == "pn9_msb" else [_REV8[v] for v in _PN9]
    return bytes(v ^ w[i] for i, v in enumerate(b))


def gfsk_frame_bits(uid, packet: bytes, variant="pn9_msb") -> np.ndarray:
    """Synthesis / tests: preamble 0x55 0x55, sync UID5 UID4, whitened packet (MSB first)."""
    by = bytes([0x55, 0x55, uid[5], uid[4]]) + gfsk_whiten(packet, variant)
    return np.array([(v >> (7 - i)) & 1 for v in by for i in range(8)], np.uint8)


def _bits_to_bytes(bits):
    n = len(bits) // 8
    return bytes(int("".join(map(str, bits[8 * i:8 * i + 8])), 2) for i in range(n))


def gfsk_try_decode(y, ofs, rate=GFSK_RATE):
    """One 300 kb/s FSK burst at ~8 samples/bit -> dict (verified SYNC packet: UID4 / UID5, OTA version, whitening)
    or None."""
    from . import sik as SIK
    for ph, bits in SIK.demod_bits(y, ofs, rate):
        for inv in (False, True):
            b = (1 - bits) if inv else bits
            d = np.flatnonzero(np.diff(b.astype(np.int8)) == 0)
            edges = np.r_[-1, d, len(b) - 1]
            for a0, a1 in zip(edges[:-1], edges[1:]):
                if a1 - a0 < 12:                                    # >= 12 bits of the 16-bit preamble
                    continue
                for s0 in range(int(a1) - 1, int(a1) + 3):          # the sync starts where the alternation ends
                    for L in (8, 13):
                        seg = b[s0:s0 + 16 + 8 * L]
                        if len(seg) < 16 + 8 * L:
                            continue
                        raw = _bits_to_bytes(seg)
                        for wv in GFSK_WHITENING:
                            p = parse(gfsk_whiten(raw[2:], wv))
                            if not p or p["packet_type"] != "SYNC" or not p["crc_ok"]:
                                continue
                            f = p["sync"]
                            if (raw[0], raw[1]) != (f["uid5"], f["uid4"]):
                                continue
                            return dict(protocol="ExpressLRS GFSK (LR1121)", sync_word=f"0x{raw[0]:02X}{raw[1]:02X}",
                                        uid4=f["uid4"], uid5=f["uid5"], whitening=wv, ota=p["ota"],
                                        ota_version=p["ota_version"], binding=p["binding"], sync=f,
                                        polarity_inverted=bool(inv), crc_ok=True)
    return None
