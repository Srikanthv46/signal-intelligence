#!/usr/bin/env python3
"""
rc_link_synth.py -- synthetic RC / telemetry control-link signals for testing
det_rc_link.py (no real ELRS / Crossfire capture exists yet).

What it makes (complex baseband at any fs, default 61.44 MSps):
  * LoRa packets (Semtech CSS): N preamble upchirps, 2 sync-word symbols,
    2.25 downchirps (SFD), payload symbols; optional IQ inversion (ELRS sets it
    from UID[5] bit 0), carrier offset and chip-clock ppm error (one crystal ->
    the same ppm on carrier and chip clock, as on SX127x/SX128x).
  * GFSK bursts (BT-shaped 2-FSK) with preamble 0101.., for FHSS links
    (FrSky-D16-like, ELRS FLRC / LR1121-GFSK-like, BLE-like confusers).
  * Whole links: ExpressLRS packet timing (interval per air rate, hop every
    FHSShopInterval packets, the real ELRS hop-sequence builder + LCG rng,
    MCU-timer ppm) and generic GFSK FHSS links.
  * embed(): add the bursts to a real capture at a chosen in-band SNR,
    measured against that capture's own per-bin noise floor.

Sources: ExpressLRS master src/src/common.cpp (air-rate table),
src/lib/FHSS/FHSS.cpp + random.cpp (domains, sequence builder, rng),
src/lib/SX1280Driver/SX1280.cpp (implicit header, CRC off, IQ inversion).
Exact over-the-air preamble symbol counts and SF5/SF6 payload symbol counts
are NOT verified against hardware (see PAYLOAD_SYMS note).
"""
from __future__ import annotations
import numpy as np

FS = 61.44e6

# ---------------------------------------------------------------- ExpressLRS tables
# (name, band, sf, bw_hz, cr_den(4/x), preamble, payload_bytes, hop_interval_pkts, interval_us)
# from ExpressLRS_AirRateConfig (SX128x / SX127x builds), master branch Oct 2026.
ELRS_LORA_MODES = {
    "2g4_500hz":     ("2G4", 5, 812.5e3, 6, 12, 8, 4, 2000),
    "2g4_333hz_8ch": ("2G4", 5, 812.5e3, 8, 12, 13, 4, 3003),
    "2g4_250hz":     ("2G4", 6, 812.5e3, 8, 14, 8, 4, 4000),
    "2g4_150hz":     ("2G4", 7, 812.5e3, 8, 12, 8, 4, 6666),
    "2g4_100hz_8ch": ("2G4", 7, 812.5e3, 8, 12, 13, 4, 10000),
    "2g4_50hz":      ("2G4", 8, 812.5e3, 8, 12, 8, 2, 20000),
    "900_200hz":     ("900", 6, 500e3, 7, 8, 8, 4, 5000),
    "900_100hz_8ch": ("900", 6, 500e3, 8, 8, 13, 4, 10000),
    "900_100hz":     ("900", 7, 500e3, 7, 8, 8, 4, 10000),
    "900_50hz":      ("900", 8, 500e3, 7, 10, 8, 4, 20000),
    "900_25hz":      ("900", 9, 500e3, 7, 10, 8, 2, 40000),
}
# FLRC (SX1280): 650 kb/s, BW 0.6 MHz, BT 1.0, CR 1/2, 32-bit preamble, 32-bit sync, 3-byte CRC
ELRS_FLRC = dict(bitrate=650e3, bt=1.0, dev=162.5e3, pre_bits=32, sync_bits=32, payload_bytes=8,
                 crc_bytes=3, cr=2, hop=2, interval_us=1000)   # dev: h~0.5 assumed (unverified)

ELRS_DOMAINS = {   # name: (f_start, f_stop, n_channels)   FHSS.cpp
    "ISM2G4": (2400.4e6, 2479.4e6, 80), "FCC915": (903.5e6, 926.9e6, 40),
    "AU915": (915.5e6, 926.9e6, 20), "EU868": (863.275e6, 869.575e6, 13),
    "IN866": (865.375e6, 866.95e6, 4), "AU433": (433.42e6, 434.42e6, 3),
    "EU433": (433.1e6, 434.45e6, 3), "US433": (433.25e6, 438.0e6, 8),
    "US433W": (423.5e6, 438.0e6, 20), "TH920": (920.5e6, 924.7e6, 8),
}


def elrs_channels(domain="ISM2G4"):
    f0, f1, n = ELRS_DOMAINS[domain]
    return f0 + np.arange(n) * (f1 - f0) / (n - 1)


class _ElrsRng:
    """random.cpp: 32-bit LCG (a=214013, c=2531011) mod 2^31, output seed>>16."""
    def __init__(self, seed): self.s = seed & 0xFFFFFFFF
    def rng(self):
        self.s = ((214013 * self.s + 2531011) & 0xFFFFFFFF) % 2147483648
        return self.s >> 16
    def rngN(self, mx): return self.rng() % mx


def elrs_hop_sequence(seed, n_ch, seq_len=256):
    """FHSSrandomiseFHSSsequenceBuild(): sync channel (n/2) at every block start,
    other entries shuffled within their block of n_ch."""
    count = (seq_len // n_ch) * n_ch
    sync = n_ch // 2
    seq = np.empty(count, int)
    for i in range(count):
        seq[i] = sync if i % n_ch == 0 else (0 if i % n_ch == sync else i % n_ch)
    r = _ElrsRng(seed)
    for i in range(count):
        if i % n_ch:
            off = (i // n_ch) * n_ch
            j = off + r.rngN(n_ch - 1) + 1
            seq[i], seq[j] = seq[j], seq[i]
    return seq


def lora_sf5_implicit_symbols(payload, cr):
    """Raw payload symbols (relative to the preamble, from the SFD on) of an implicit-header
    SF5 SX126x packet with CRC: the inverse of lora_phy.decode_sf5_implicit (2 extra symbols,
    first block CR 4/8 at full SF, then blocks of 4+cr symbols). Feed to lora_packet(5, ...)."""
    from .lora_phy import whitening_sequence, crc16, SF5_EXTRA_SYMS
    sf, N = 5, 32
    pay = np.frombuffer(bytes(payload), np.uint8)
    L = len(pay)
    c = crc16(pay[:-2]) ^ int(pay[-1]) ^ (int(pay[-2]) << 8) if L >= 2 else crc16(b"") ^ int(pay[-1])
    data = list(pay ^ whitening_sequence(L)) + [c & 0xFF, c >> 8]
    nib = [v for b in data for v in (b & 0xF, b >> 4)]
    out = [1] * SF5_EXTRA_SYMS
    for k, cw_len in enumerate([8] + [4 + cr] * ((max(len(nib) - sf, 0) + sf - 1) // sf)):
        blk = (nib[k * sf:(k + 1) * sf] + [0] * sf)[:sf]
        cws = []
        for n in blk:
            d = [(n >> i) & 1 for i in range(4)]                       # cw[0..3] = nibble bits 0..3
            p = [d[0] ^ d[1] ^ d[2], d[1] ^ d[2] ^ d[3], d[0] ^ d[1] ^ d[3]]
            cw = d + (p if cw_len >= 7 else p[:2] if cw_len == 6 else [d[0] ^ d[1] ^ d[2] ^ d[3]])
            if cw_len == 8:
                cw.append(sum(cw) % 2)
            cws.append(cw[:cw_len])
        for i in range(cw_len):                                        # diagonal interleave (inverse)
            g = 0
            for j in range(sf):
                g = (g << 1) | cws[(i - j - 1) % sf][i]
            v, s = g, g >> 1
            while s:
                v ^= s; s >>= 1                                       # inverse Gray
            out.append((v + 1) % N)
    return out


def lora_implicit_symbols(payload, sf, cr):
    """Raw payload symbols of an implicit-header SX127x packet without LoRa CRC (SF >= 7, no LDRO): first block CR
    4/8 at reduced rate (SF-2 nibbles, 2 LSBs 0), then blocks of 4+cr symbols at full rate (inverse of
    lora_phy.decode_implicit)."""
    from .lora_phy import whitening_sequence
    N = 2 ** sf
    pay = np.frombuffer(bytes(payload), np.uint8)
    data = [int(b) for b in (pay ^ whitening_sequence(len(pay)))]
    nib = [v for b in data for v in (b & 0xF, b >> 4)]
    out, k = [], 0
    blocks = [(8, sf - 2)]
    rest = max(len(nib) - (sf - 2), 0)
    blocks += [(4 + cr, sf)] * ((rest + sf - 1) // sf)
    for cw_len, sf_app in blocks:
        blk = (nib[k:k + sf_app] + [0] * sf_app)[:sf_app]
        k += sf_app
        cws = []
        for n in blk:
            d = [(n >> i) & 1 for i in range(4)]
            p = [d[0] ^ d[1] ^ d[2], d[1] ^ d[2] ^ d[3], d[0] ^ d[1] ^ d[3]]
            cw = d + (p if cw_len >= 7 else p[:2] if cw_len == 6 else [d[0] ^ d[1] ^ d[2] ^ d[3]])
            if cw_len == 8:
                cw.append(sum(cw) % 2)
            cws.append(cw[:cw_len])
        for i in range(cw_len):
            g = 0
            for j in range(sf_app):
                g = (g << 1) | cws[(i - j - 1) % sf_app][i]
            v, s = g, g >> 1
            while s:
                v ^= s; s >>= 1
            if sf_app < sf:
                v <<= 2
            out.append((v + 1) % N)
    return out


def lora_payload_symbols(sf, pl_bytes, cr_den, crc=False, implicit=True, ldro=False):
    """Semtech SX127x time-on-air formula (n_payload symbols). PAYLOAD_SYMS note:
    SX128x SF5/SF6 use a slightly different formula -- unverified; good to +-a few symbols."""
    num = 8 * pl_bytes - 4 * sf + 28 + 16 * crc - 20 * implicit
    den = 4 * (sf - 2 * ldro)
    return 8 + max(int(np.ceil(num / den)) * cr_den, 0)


# ---------------------------------------------------------------- waveforms
def _ramp(n_total, fs, ramp_us=2.0):
    k = max(1, int(ramp_us * 1e-6 * fs))
    w = np.ones(n_total)
    r = 0.5 - 0.5 * np.cos(np.pi * np.arange(k) / k)
    w[:k] = r; w[-k:] = r[::-1]
    return w


def lora_packet(sf, bw, n_pre=12, payload=None, n_payload=20, sync_word=0x12, iq_invert=False,
                cfo_hz=0.0, clock_ppm=0.0, fs=FS, rng=None, phase=None):
    """One LoRa packet as complex baseband at fs (unit amplitude). Built by
    integrating instantaneous frequency in continuous time so chip-clock ppm,
    CFO and fractional timing are exact. Returns (iq, info)."""
    rng = rng or np.random.default_rng()
    N = 2 ** sf
    bw_tx = bw * (1 + clock_ppm * 1e-6)                  # same crystal -> chip clock scales
    T = N / bw_tx
    if payload is None:
        payload = rng.integers(0, N, n_payload)
    sw = [((sync_word >> 4) & 0xF) * 8, (sync_word & 0xF) * 8]     # nibble x 8 at every SF (live SF5 / SF9)
    # symbol list: (kind, value, duration_fraction)
    syms = [("up", 0, 1.0)] * n_pre + [("up", s, 1.0) for s in sw] + [("down", 0, 1.0)] * 2 + \
           [("down", 0, 0.25)] + [("up", int(s), 1.0) for s in payload]
    pieces = []
    t_acc = 0.0
    for kind, v, frac in syms:
        ts = t_acc
        n0 = int(np.ceil(ts * fs - 1e-9)); t_acc += frac * T; n1 = int(np.ceil(t_acc * fs - 1e-9))
        t = np.arange(n0, n1) / fs - ts                   # exact time since symbol start
        if kind == "up":
            f = -bw_tx / 2 + bw_tx * np.mod(v / N + t / T, 1.0)
        else:
            f = bw_tx / 2 - bw_tx * (t / T)
        pieces.append(f)
    f = np.concatenate(pieces) + 0.0
    if iq_invert:
        f = -f
    f = f + cfo_hz
    ph = 2 * np.pi * np.cumsum(f) / fs
    ph0 = rng.uniform(0, 2 * np.pi) if phase is None else phase
    iq = np.exp(1j * (ph + ph0)) * _ramp(len(f), fs)
    info = dict(modulation="lora", sf=sf, bw_hz=bw, n_pre=n_pre, n_payload=len(payload),
                iq_invert=iq_invert, cfo_hz=cfo_hz, clock_ppm=clock_ppm,
                dur_us=len(iq) / fs * 1e6, sym_us=T * 1e6)
    return iq.astype(np.complex64), info


def gfsk_burst(bitrate, dev_hz, n_bits=None, bits=None, bt=0.5, pre_bits=32, cfo_hz=0.0,
               clock_ppm=0.0, fs=FS, rng=None):
    """BT-shaped binary FSK with alternating preamble. clock_ppm scales the bit
    clock (same crystal also moves the carrier: give cfo_hz = ppm * fc)."""
    rng = rng or np.random.default_rng()
    if bits is None:
        bits = np.r_[np.arange(pre_bits) % 2, rng.integers(0, 2, n_bits or 200)]
    br = bitrate * (1 + clock_ppm * 1e-6)
    sps = fs / br
    n = int(np.ceil(len(bits) * sps))
    idx = np.minimum((np.arange(n) / sps).astype(int), len(bits) - 1)
    nrz = 2.0 * np.asarray(bits)[idx] - 1
    # Gaussian pulse: sigma (in samples) = sqrt(ln2)/(2 pi BT) * sps
    sig = np.sqrt(np.log(2)) / (2 * np.pi * bt) * sps
    k = np.arange(-int(4 * sig) - 1, int(4 * sig) + 2)
    g = np.exp(-k ** 2 / (2 * sig ** 2)); g /= g.sum()
    fdev = np.convolve(nrz, g, "same") * dev_hz + cfo_hz
    ph = 2 * np.pi * np.cumsum(fdev) / fs + rng.uniform(0, 2 * np.pi)
    iq = np.exp(1j * ph) * _ramp(n, fs)
    return iq.astype(np.complex64), dict(modulation="gfsk", bitrate=bitrate, dev_hz=dev_hz, bt=bt,
                                         n_bits=len(bits), cfo_hz=cfo_hz, clock_ppm=clock_ppm,
                                         dur_us=n / fs * 1e6)


# ---------------------------------------------------------------- Bluetooth BR/EDR (research, signal gaps S7)
def bredr_bits(lap, n_payload=240, rng=None):
    """Access code (preamble + sync word, bit 0 first + trailer), an 18 x 3 header and a payload (random bits)."""
    from .bredr import syncword
    rng = rng or np.random.default_rng()
    sw = syncword(lap)
    s = [(sw >> k) & 1 for k in range(64)]
    pre = [s[0], 1 - s[0], s[0], 1 - s[0]][::-1]          # alternating into the first sync bit
    tr = [1 - s[63], s[63], 1 - s[63], s[63]]
    hdr = [b for b in rng.integers(0, 2, 18).tolist() for _ in range(3)]
    return np.array(pre + s + tr + hdr + rng.integers(0, 2, n_payload).tolist(), np.uint8)


def bredr_packet(lap, h=0.32, n_payload=240, fs=FS, rng=None, cfo_hz=0.0):
    """One BR packet: 1 Mb/s GFSK, BT 0.5, modulation index h."""
    return gfsk_burst(1e6, h * 1e6 / 2, bits=bredr_bits(lap, n_payload, rng), bt=0.5, pre_bits=0, cfo_hz=cfo_hz,
                      fs=fs, rng=rng)


# ---------------------------------------------------------------- Bluetooth LE advertising (Core spec Vol 6 Part B)
BLE_AA = 0x8E89BED6


def _ble_crc_bits(pdu_bits, init=0x555555):
    reg = [(init >> i) & 1 for i in range(24)]
    for b in pdu_bits:
        fb = reg[23] ^ b
        reg = [fb] + reg[:23]
        for p in (1, 3, 4, 6, 9, 10):
            reg[p] ^= fb
    return reg[::-1]


def _ble_whiten(bits, ch):
    reg = [1] + [(ch >> (5 - i)) & 1 for i in range(6)]
    out = []
    for b in bits:
        o = reg[6]
        out.append(b ^ o)
        reg = [o] + reg[:6]
        reg[4] ^= o
    return out


def ble_adv_bits(adv_addr: bytes, adv_data: bytes, pdu_type=2, ch=37, tx_random=True):
    """Air bits (LSB first) of a legacy advertising PDU: preamble, access address, whitened PDU + CRC-24."""
    def lsb(by):
        return [(v >> i) & 1 for v in by for i in range(8)]
    payload = bytes(adv_addr) + bytes(adv_data)
    hdr = bytes([pdu_type | (0x40 if tx_random else 0), len(payload)])
    pdu = lsb(hdr) + lsb(payload)
    aa = [(BLE_AA >> i) & 1 for i in range(32)]
    pre = [aa[0] ^ (i % 2 == 0) for i in range(8)]          # alternating, last bit differs from AA bit 0
    return pre + aa + _ble_whiten(pdu + _ble_crc_bits(pdu), ch)


def ble_adv_burst(adv_addr: bytes, adv_data: bytes, ch=37, cfo_hz=0.0, fs=FS, rng=None):
    """GFSK 1 Mb/s, h 0.5 (deviation 250 kHz), BT 0.5 advertising packet."""
    return gfsk_burst(1e6, 250e3, bits=np.array(ble_adv_bits(adv_addr, adv_data, ch=ch), np.uint8), bt=0.5,
                      cfo_hz=cfo_hz, fs=fs, rng=rng)


def ext_adv_payload(adv_addr: bytes, adv_data: bytes = b"", did=0x123, sid=1, aux=None, mode=0) -> bytes:
    """Common extended advertising payload: ext header (AdvA, ADI, optional AuxPtr = (channel, offset in 30 us units,
    phy 2 = Coded)) + AdvData (Core Vol 6 Part B 2.3.4)."""
    flags = 0x01 | 0x08 | (0x10 if aux else 0)
    eh = bytes([flags]) + bytes(adv_addr) + bytes([(did & 0xFF), ((did >> 8) & 0x0F) | (sid << 4)])
    if aux:
        ch, off, phy = aux
        v = (ch & 0x3F) | ((off & 0x1FFF) << 8) | ((phy & 7) << 21)
        eh += bytes([v & 0xFF, (v >> 8) & 0xFF, (v >> 16) & 0xFF])
    return bytes([len(eh) | (mode << 6)]) + eh + bytes(adv_data)


def ble_coded_burst(payload: bytes, pdu_type=7, ch=37, ci=0, cfo_hz=0.0, fs=FS, rng=None, tx_random=True):
    """LE Coded PHY packet (S=8 if ci 0, S=2 if ci 1) of a PDU with header (pdu_type, len(payload))."""
    from .ble_coded import coded_symbols
    def lsb(by):
        return [(v >> i) & 1 for v in by for i in range(8)]
    hdr = bytes([pdu_type | (0x40 if tx_random else 0), len(payload)])
    pdu = lsb(hdr) + lsb(payload)
    aa = [(BLE_AA >> i) & 1 for i in range(32)]
    sym = np.array(coded_symbols(aa, ci, _ble_whiten(pdu + _ble_crc_bits(pdu), ch)), np.uint8)
    return gfsk_burst(1e6, 250e3, bits=sym, bt=0.5, cfo_hz=cfo_hz, fs=fs, rng=rng)


def crossfire_burst(data: bytes, cfo_hz=0.0, fs=FS, rng=None, preamble_bytes=5):
    """TBS Crossfire 150 Hz FSK packet (research): 0xAA preamble, sync 0x3E6B3E, PN9-whitened data (MSB first),
    2-FSK 85.1 kb/s +-42.48 kHz (plain FSK: BT high)."""
    from .sik import WHITEN
    wd = bytes(v ^ w for v, w in zip(data, WHITEN["pn9_msb_first"][:len(data)]))
    by = bytes([0xAA] * preamble_bytes) + bytes([0x3E, 0x6B, 0x3E]) + wd
    bits = np.array([(v >> (7 - i)) & 1 for v in by for i in range(8)], np.uint8)
    return gfsk_burst(85.1e3, 42.48e3, bits=bits, bt=2.0, cfo_hz=cfo_hz, fs=fs, rng=rng)


FRSKY_PHY = {"D8": (31.04e3, 31.7e3, [0xD3, 0x91]), "D16 v1": (70e3, 57.1e3, [0xD3, 0x91, 0xD3, 0x91]),
             "D16 v2": (76.97e3, 57.1e3, [0xD3, 0x91, 0xD3, 0x91]), "D16 LBT": (99.98e3, 69.8e3, [0xD3, 0x91, 0xD3, 0x91])}


def frsky_burst(variant: str, payload: bytes, cfo_hz=0.0, fs=FS, rng=None, preamble_bytes=4, tail_bits=8):
    """FrSky D8 / D16 CC2500 packet (research, file 10): 0xAA preamble, sync, then cc2500.d8_packet / d16_packet
    (length, payload, CRCs; D16 v2 adds the hardware CRC); GFSK BT 0.5 (LBT: 2-FSK, BT high). payload excludes the
    length byte and the CRCs; tail_bits model the PA ramp-down after the last byte."""
    from .cc2500 import d8_packet, d16_packet
    br, dev, sync = FRSKY_PHY[variant]
    pkt = d8_packet(payload) if variant == "D8" else d16_packet(payload, v2=variant == "D16 v2")
    by = bytes([0xAA] * preamble_bytes + sync) + pkt
    bits = np.r_[np.array([(v >> (7 - i)) & 1 for v in by for i in range(8)], np.uint8), np.zeros(tail_bits, np.uint8)]
    return gfsk_burst(br, dev, bits=bits, bt=2.0 if variant == "D16 LBT" else 0.5, cfo_hz=cfo_hz, fs=fs, rng=rng)


def afhds2a_burst(payload: bytes, cfo_hz=0.0, fs=FS, rng=None, tail_bits=8, **conv):
    """FlySky AFHDS 2A A7105 packet (research, file 10): 4-byte preamble (0101... before the ID's leading 0), ID code
    0x5475C52A, 38-byte payload + CRC-16 (7,4)-Hamming coded per nibble (afhds2a.encode_packet conventions); plain
    2-FSK 500 kb/s +-186 kHz."""
    from .afhds2a import ID_CODE, encode_packet
    pre = np.arange(32) % 2                                              # 0101...: the ID starts with a 0
    idb = [(ID_CODE >> (31 - i)) & 1 for i in range(32)]
    bits = np.r_[pre, idb, encode_packet(payload, **conv), np.zeros(tail_bits)].astype(np.uint8)
    return gfsk_burst(500e3, 186e3, bits=bits, bt=2.0, cfo_hz=cfo_hz, fs=fs, rng=rng)


def bayang_payload(kind="data", tx=b"\x11\x22\x33\x44\x55", sticks=(512, 512, 0, 512)) -> bytes:
    """A 15-byte Bayang payload (research, file 10): bind (0xA4: TX address, hop channels) or data (0xA5: sticks,
    flags); byte 14 = 8-bit sum of bytes 0-13."""
    if kind == "bind":
        h1 = (tx[3] & 0x1F) + 0x10
        p = [0xA4, *tx, 0, h1, h1 + 0x20, h1 + 0x40, tx[0], tx[1], tx[2], 0x0A]
    else:
        st = []
        for v in sticks:
            st += [((v >> 8) & 0x03) | 0x7C, v & 0xFF]
        p = [0xA5, 0xFA, 0x00, 0x00, *st, tx[2], 0x0A]
    return bytes(p + [sum(p) & 0xFF])


def xn297_burst(addr: bytes, payload: bytes, cfo_hz=0.0, fs=FS, rng=None, tail_bits=8):
    """XN297 packet (research, file 10): nRF-style 0x55 preamble, prefix 0x710F55, scrambled address / payload, CRC
    (length constant 0: unknown); GFSK 1 Mb/s, +-160 kHz [U], BT 0.5."""
    from .xn297 import PREFIX, encode
    by = bytes([0x55, PREFIX >> 16, (PREFIX >> 8) & 0xFF, PREFIX & 0xFF]) + encode(addr, payload)
    bits = np.r_[np.array([(v >> (7 - i)) & 1 for v in by for i in range(8)], np.uint8), np.zeros(tail_bits, np.uint8)]
    return gfsk_burst(1e6, 160e3, bits=bits, bt=0.5, cfo_hz=cfo_hz, fs=fs, rng=rng)


def oqpsk_burst(psdu: bytes, cfo_hz=0.0, fs=FS, rng=None):
    """IEEE 802.15.4 2.4 GHz frame (research, file 10): SHR + PHR + psdu (MPDU incl. FCS) as 2 Mchip/s O-QPSK with
    half-sine chips, evaluated at any sample rate fs (even chips on I, odd on Q one chip later)."""
    from .ieee802154 import chips_for
    rng = rng or np.random.default_rng()
    c = 2.0 * chips_for(psdu) - 1
    tc = 1 / 2e6
    t = np.arange(int((len(c) + 2) * tc * fs)) / fs
    k = (t // tc).astype(int)                          # chip interval index
    i = np.zeros(len(t))
    q = np.zeros(len(t))
    for parity, out in ((0, i), (1, q)):
        start = k - ((k - parity) % 2)                 # the chip of this rail whose 2-chip pulse covers t
        ok = (start >= 0) & (start < len(c))
        ph = (t - start * tc) / (2 * tc)
        out[ok] = c[start[ok]] * np.sin(np.pi * ph[ok])
    y = (i + 1j * q) * np.exp(1j * (2 * np.pi * cfo_hz * t + rng.uniform(0, 2 * np.pi)))
    return y.astype(np.complex64), dict(modulation="oqpsk", chip_rate=2e6, n_chips=len(c), dur_us=t[-1] * 1e6)


def flarm_burst(radio_id: int, addr_type=0x20, cfo_hz=0.0, fs=FS, rng=None, n_payload=22):
    """FLARM-like beacon (research, file 10): Manchester preamble, sync F5 31 FA B6, 24-bit ID (LE), address type,
    then random (obfuscated-looking) payload bytes; 100 kchip/s 2-FSK +-50 kHz."""
    from .flarm import SYNC, manchester
    rng = rng or np.random.default_rng()
    by = bytes([SYNC >> 24, (SYNC >> 16) & 255, (SYNC >> 8) & 255, SYNC & 255, radio_id & 255, (radio_id >> 8) & 255,
                radio_id >> 16, addr_type]) + bytes(rng.integers(0, 256, n_payload).tolist())
    bits = [(v >> (7 - i)) & 1 for v in by for i in range(8)]
    chips = np.r_[np.arange(16) % 2, manchester(bits), np.zeros(8)].astype(np.uint8)
    return gfsk_burst(100e3, 50e3, bits=chips, bt=2.0, cfo_hz=cfo_hz, fs=fs, rng=rng)


def remote_id_adv_data(odid_body: bytes) -> bytes:
    """AD structure carrying ASTM F3411 over BLE: Service Data (0x16), UUID 0xFFFA, app code 0x0D, body."""
    v = bytes([0xFA, 0xFF, 0x0D]) + bytes(odid_body)
    return bytes([len(v) + 1, 0x16]) + v


# ---------------------------------------------------------------- MAVLink over SiK (research: 05_sik_mavlink_details.md)
def mavlink_frame(msgid: int, payload: bytes, seq=0, sysid=1, compid=1, v2=False) -> bytes:
    """A MAVLink v1 / v2 frame with its X.25 CRC (CRC_EXTRA from core.mavlink.MSGS); v2 payloads zero-trimmed."""
    from spectrum_intel.core import mavlink as MAV
    payload = bytes(payload)
    if v2:
        while len(payload) > 1 and payload[-1] == 0:
            payload = payload[:-1]
        body = bytes([len(payload), 0, 0, seq, sysid, compid, msgid & 0xFF, (msgid >> 8) & 0xFF, msgid >> 16]) + payload
        stx = 0xFD
    else:
        body = bytes([len(payload), seq, sysid, compid, msgid]) + payload
        stx = 0xFE
    crc = MAV.x25(bytes([MAV.MSGS[msgid][1]]), MAV.x25(body))
    return bytes([stx]) + body + bytes([crc & 0xFF, crc >> 8])


def sik_packet_bits(netid: int, data: bytes, whitening="pn9_lsb_first", crc="ibm_msb_init0000", crc_hi_first=True,
                    preamble_nibbles=16):
    """Air bits (MSB first) of a SiK packet: 1010 preamble, sync 0x2DD4, whitened [netid hi, lo, len, data, CRC]."""
    from . import sik as SIK
    body = bytes([netid >> 8, netid & 0xFF, len(data)]) + bytes(data)
    c = SIK.crc16_variants(body)[crc]
    body += bytes([c >> 8, c & 0xFF]) if crc_hi_first else bytes([c & 0xFF, c >> 8])
    w = np.frombuffer(body, np.uint8) ^ SIK.WHITEN[whitening][:len(body)]
    pre = [1, 0] * (2 * preamble_nibbles)
    sync = [(SIK.SYNC >> (15 - i)) & 1 for i in range(16)]
    return np.array(pre + sync + list(np.unpackbits(w)), np.uint8)


def sik_burst(netid: int, data: bytes, rate=64e3, dev_hz=63.75e3, cfo_hz=0.0, fs=FS, rng=None, **kw):
    """GFSK SiK packet (Si443x GFSK, BT 0.5, deviation = 625 Hz x reg 0x72)."""
    return gfsk_burst(rate, dev_hz, bits=sik_packet_bits(netid, data, **kw), bt=0.5, cfo_hz=cfo_hz, fs=fs, rng=rng)


# ---------------------------------------------------------------- links
def elrs_link(mode="2g4_500hz", dur_s=0.0325, uid_seed=None, crystal_ppm=None, timer_ppm=None,
              iq_invert=None, domain="ISM2G4", t0_us=None, rng=None, fs=FS, f_window=None):
    """ExpressLRS LoRa uplink: one packet every interval_us (MCU timer with its own
    ppm), channel from the ELRS hop sequence, advancing every FHSShopInterval
    packets. Returns list of bursts: dict(t0_us, f_hz, iq, info, channel, nonce).
    f_window=(f_lo, f_hi): drop packets whose channel is outside (capture band)."""
    rng = rng or np.random.default_rng()
    band, sf, bw, cr, pre, pl, hop, ival = ELRS_LORA_MODES[mode]
    uid_seed = int(rng.integers(0, 2 ** 31)) if uid_seed is None else uid_seed
    crystal_ppm = rng.uniform(-20, 20) if crystal_ppm is None else crystal_ppm
    timer_ppm = rng.uniform(-30, 30) if timer_ppm is None else timer_ppm
    iq_invert = bool(rng.integers(0, 2)) if iq_invert is None else iq_invert
    chans = elrs_channels(domain)
    seq = elrs_hop_sequence(uid_seed, len(chans))
    n_pay = lora_payload_symbols(sf, pl, cr)
    ptr = int(rng.integers(0, len(seq)))
    nonce = int(rng.integers(0, 256))
    t = rng.uniform(0, ival) if t0_us is None else t0_us
    T = ival * (1 + timer_ppm * 1e-6)
    out = []
    while t < dur_s * 1e6:
        f = chans[seq[ptr % len(seq)]]
        cfo = crystal_ppm * 1e-6 * f
        iq, info = lora_packet(sf, bw, pre, n_payload=n_pay, iq_invert=iq_invert, cfo_hz=cfo,
                               clock_ppm=crystal_ppm, fs=fs, rng=rng)
        info.update(mode=mode, crystal_ppm=crystal_ppm, timer_ppm=timer_ppm, uid_seed=uid_seed)
        if (f_window is None or (f_window[0] <= f - bw / 2 and f + bw / 2 <= f_window[1])) and \
                t + info["dur_us"] < dur_s * 1e6:
            out.append(dict(t0_us=t, f_hz=f, iq=iq, info=info, channel=int(seq[ptr % len(seq)]),
                            nonce=nonce, bw_hz=bw))
        nonce += 1
        if nonce % hop == 0:
            ptr += 1
        t += T
    return out


def gfsk_fhss_link(bitrate=70.0e3, dev_hz=57e3, n_bytes=30, interval_us=9000.0, channels=None,
                   bt=0.5, pre_bits=32, dur_s=0.0325, crystal_ppm=None, timer_ppm=None,
                   t0_us=None, rng=None, fs=FS, f_window=None, name="frsky_d16_like"):
    """Generic GFSK FHSS link: one burst per interval on a pseudo-random channel
    from `channels` (Hz). Default = FrSky ACCST D16-like (CC2500 70 kBd, +-57 kHz,
    30-byte packets, 9 ms, 47 hops on a ~300 kHz raster)."""
    rng = rng or np.random.default_rng()
    if channels is None:
        channels = 2400e6 + 0.2996e6 * np.sort(rng.choice(np.arange(1, 235), 47, replace=False))
    crystal_ppm = rng.uniform(-20, 20) if crystal_ppm is None else crystal_ppm
    timer_ppm = rng.uniform(-30, 30) if timer_ppm is None else timer_ppm
    order = rng.permutation(len(channels))
    t = rng.uniform(0, interval_us) if t0_us is None else t0_us
    out, k = [], 0
    while t < dur_s * 1e6:
        f = channels[order[k % len(order)]]
        iq, info = gfsk_burst(bitrate, dev_hz, n_bits=8 * n_bytes, bt=bt, pre_bits=pre_bits,
                              cfo_hz=crystal_ppm * 1e-6 * f, clock_ppm=crystal_ppm, fs=fs, rng=rng)
        info.update(link=name, crystal_ppm=crystal_ppm, timer_ppm=timer_ppm)
        bw = 2 * (dev_hz + bitrate / 2)
        if (f_window is None or (f_window[0] <= f - bw and f + bw <= f_window[1])) and \
                t + info["dur_us"] < dur_s * 1e6:
            out.append(dict(t0_us=t, f_hz=f, iq=iq, info=info, channel=int(order[k % len(order)]),
                            bw_hz=bw))
        k += 1
        t += interval_us * (1 + timer_ppm * 1e-6)
    return out


# ---------------------------------------------------------------- embedding
def floor_density(x, nfft=1024):
    """Per-bin noise floor (same estimator as drone_common.noise_floor_psd)."""
    from spectrum_intel.core import common as C
    return C.noise_floor_psd(x, nfft)["floor_bin"]


def embed(x, fs, fc, bursts, snr_db, floor_bin=None, nfft=1024):
    """Add bursts to capture x so each has in-band SNR snr_db (signal power over
    the capture's per-bin noise floor integrated over the signal bandwidth --
    the definition drone_common.measure_snr uses). snr_db may be a scalar or a
    per-burst list. Returns (y, bursts) with t1_us / snr_target filled in."""
    fl = floor_density(x, nfft) if floor_bin is None else floor_bin
    fbin = (np.arange(nfft) - nfft // 2) / nfft * fs + fc
    y = x.astype(np.complex64).copy()
    n = len(x)
    snrs = np.broadcast_to(np.asarray(snr_db, float), (len(bursts),))
    for b, s in zip(bursts, snrs):
        bw = b.get("bw_hz", 1e6)
        sel = np.abs(fbin - b["f_hz"]) <= bw / 2
        n0 = fl[sel].mean() if sel.any() else np.median(fl)
        # floor_bin = per-bin power of a 1-sample-power white noise -> noise power in bw = n0 * bw/fs * nfft/nfft
        p_noise = n0 * bw / fs
        a = int(round(b["t0_us"] * 1e-6 * fs))
        iq = b["iq"]
        m = min(len(iq), n - a)
        if m <= 0:
            continue
        g = np.sqrt(10 ** (s / 10) * p_noise / np.mean(np.abs(iq) ** 2))
        tt = np.arange(a, a + m)
        y[a:a + m] += (g * iq[:m] * np.exp(2j * np.pi * (b["f_hz"] - fc) / fs * tt)).astype(np.complex64)
        b["t1_us"] = b["t0_us"] + m / fs * 1e6
        b["snr_target"] = float(s)
    return y, bursts


def noise_capture(n=2_000_000, rng=None):
    rng = rng or np.random.default_rng()
    return ((rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)).astype(np.complex64)


def radiolink_burst(fifo: bytes, sync=0xABCD, cfo_hz=0.0, fs=FS, rng=None, tail_bits=8):
    """RadioLink CC2500 packet (research, PLAN G4): 4 x 0xAA, the sync word twice, PN9-whitened length + payload +
    CRC (radiolink.air_frame); GFSK 37.99 kBd +-38.09 kHz, BT 0.5."""
    from .radiolink import RATE, DEV, air_frame
    by = bytes([0xAA] * 4) + air_frame(sync, fifo)
    bits = np.r_[np.array([(v >> (7 - i)) & 1 for v in by for i in range(8)], np.uint8), np.zeros(tail_bits, np.uint8)]
    return gfsk_burst(RATE, DEV, bits=bits, bt=0.5, cfo_hz=cfo_hz, fs=fs, rng=rng)


def a7105_burst(id_code: int, payload: bytes, rate=500e3, crc=None, cfo_hz=0.0, fs=FS, rng=None, tail_bits=8):
    """A7105 packet without FEC (AFHDS 1, Hubsan; research, PLAN G5): 32-bit alternating preamble flowing into the ID
    code, the ID code, the payload, an optional CRC-16 (poly, init) big-endian; 2-FSK +-186 kHz."""
    pre = (np.arange(32) + ((id_code >> 31) & 1) + 1) % 2              # 0101... ending opposite the ID's first bit
    idb = [(id_code >> (31 - i)) & 1 for i in range(32)]
    body = bytes(payload)
    if crc:
        from .afhds2a import crc16
        c = crc16(body, *crc)
        body += bytes([c >> 8, c & 0xFF])
    pb = [(v >> (7 - i)) & 1 for v in body for i in range(8)]
    bits = np.r_[pre, idb, pb, np.zeros(tail_bits)].astype(np.uint8)
    return gfsk_burst(rate, 186e3, bits=bits, bt=2.0, cfo_hz=cfo_hz, fs=fs, rng=rng)


def _bits_of(by: bytes):
    return np.array([(v >> (7 - i)) & 1 for v in by for i in range(8)], np.uint8)


def hott_burst(data: bytes, sync=0xD391, cfo_hz=0.0, fs=FS, rng=None, tail_bits=8):
    """Graupner HoTT CC2500 packet (research, PLAN G12): 4 x 0xAA, the sync twice, the FEC stream of the whitened
    data + CRC (hott.air_bits); MSK 249.939 kBd (CPFSK h = 0.5, no Gaussian filter)."""
    from .hott import RATE, DEV, air_bits
    bits = np.r_[_bits_of(bytes([0xAA] * 4)), air_bits(sync, data), np.zeros(tail_bits, np.uint8)]
    return gfsk_burst(RATE, DEV, bits=bits, bt=5.0, cfo_hz=cfo_hz, fs=fs, rng=rng)


def hitec_burst(fifo: bytes, sync: int, cfo_hz=0.0, fs=FS, rng=None, tail_bits=8):
    """Hitec Optima CC2500 packet (research, PLAN G12): 24 x 0xAA, the sync (= TX ID bytes 3, 2) twice, whitened
    FIFO + CRC (hitec.hitec_air); MSK 249.939 kBd."""
    from .hitec import RATE, DEV, hitec_air
    bits = np.r_[_bits_of(bytes([0xAA] * 24) + hitec_air(sync, fifo)), np.zeros(tail_bits, np.uint8)]
    return gfsk_burst(RATE, DEV, bits=bits, bt=5.0, cfo_hz=cfo_hz, fs=fs, rng=rng)


def sfhss_burst(packet: bytes, crc_init=0xFFFF, cfo_hz=0.0, fs=FS, rng=None, tail_bits=8):
    """Futaba S-FHSS CC2500 packet (research, PLAN G12): 4 x 0xAA, D391 twice, 13 bytes + CRC (CC2400 mode);
    2-FSK 128.143 kBd +-38.086 kHz."""
    from .hitec import SF_RATE, SF_DEV, sfhss_air
    bits = np.r_[_bits_of(bytes([0xAA] * 4) + sfhss_air(packet, crc_init)), np.zeros(tail_bits, np.uint8)]
    return gfsk_burst(SF_RATE, SF_DEV, bits=bits, bt=5.0, cfo_hz=cfo_hz, fs=fs, rng=rng)


def adsl_burst(pkt: bytes, cfo_hz=0.0, fs=FS, rng=None):
    """ADS-L SRD-860 M-band packet (research, PLAN G10): chip preamble, Manchester(F5 72 4B + length + header + data +
    CRC-24); 100 kchip/s 2-GFSK +-50 kHz."""
    from .adsl import air_bits
    chips = np.r_[np.arange(16) % 2, air_bits(pkt), np.zeros(8)].astype(np.uint8)
    return gfsk_burst(100e3, 50e3, bits=chips, bt=0.5, cfo_hz=cfo_hz, fs=fs, rng=rng)


def ogn_burst(pkt26: bytes, cfo_hz=0.0, fs=FS, rng=None):
    """OGN tracker packet (research, PLAN G10): chip preamble, Manchester(0A F3 65 6C + 26 bytes); FLARM PHY."""
    from .ogn import air_chips
    chips = np.r_[np.arange(16) % 2, air_chips(pkt26), np.zeros(8)].astype(np.uint8)
    return gfsk_burst(100e3, 50e3, bits=chips, bt=2.0, cfo_hz=cfo_hz, fs=fs, rng=rng)
