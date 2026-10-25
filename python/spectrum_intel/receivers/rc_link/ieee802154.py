"""
ieee802154.py -- IEEE 802.15.4 2.4 GHz O-QPSK frames (XBee 802.15.4 / Zigbee telemetry; research, PLAN T4;
research/drone_transmitters/10_cc2500_a7105_xn297_flarm_xbee.md section 6). Synthetic validation only.

PHY: 2 Mchip/s O-QPSK with half-sine chips (MSK-equivalent), 16 x 32-chip symbols (4 bits each, low nibble first),
channels 11-26 at 2405 + 5 (k - 11) MHz. Frame: 8 zero symbols, SFD 0xA7, PHR (length 1-127), MPDU, FCS CRC-16 KERMIT
(reflected 0x8408, init 0) over the MPDU, sent little-endian.
Receiver: non-coherent. The FM discriminator at the chip rate gives one frequency sign per chip pair; each symbol's
expected 31 inner signs (its own chips only) are templates made by modulating the symbol through the same shaping and
discriminator (no analytic sign rule assumed). Symbols are decided by template correlation; the preamble + SFD fix the
timing; a CRC pass is the label gate. MAC header: frame type, PAN IDs, addresses; a Digi OUI (00 13 A2) in an extended
address names the vendor.
Nothing is transmitted or controlled.
"""
from __future__ import annotations

import numpy as np

CHIP_RATE = 2e6
CHIPS = ["11011001110000110101001000101110", "11101101100111000011010100100010", "00101110110110011100001101010010",
         "00100010111011011001110000110101", "01010010001011101101100111000011", "00110101001000101110110110011100",
         "11000011010100100010111011011001", "10011100001101010010001011101101", "10001100100101100000011101111011",
         "10111000110010010110000001110111", "01111011100011001001011000000111", "01110111101110001100100101100000",
         "00000111011110111000110010010110", "01100000011101111011100011001001", "10010110000001110111101110001100",
         "11001001011000000111011110111000"]
CHIP_ARR = np.array([[int(c) for c in s] for s in CHIPS], np.int8)
FRAME_TYPES = {0: "beacon", 1: "data", 2: "ack", 3: "mac command"}
DIGI_OUI = "0013a2"


def channel(f_hz, tol_hz=150e3):
    k = int(round((f_hz - 2405e6) / 5e6)) + 11
    return k if 11 <= k <= 26 and abs(f_hz - (2405e6 + 5e6 * (k - 11))) <= tol_hz else None


def crc16_kermit(data: bytes) -> int:
    c = 0
    for b in data:
        c ^= b
        for _ in range(8):
            c = (c >> 1) ^ 0x8408 if c & 1 else c >> 1
    return c


def chips_for(psdu: bytes) -> np.ndarray:
    """Chips of SHR + PHR + PSDU (psdu = MPDU incl. FCS)."""
    by = bytes(4) + bytes([0xA7, len(psdu)]) + bytes(psdu)
    syms = [s for b in by for s in (b & 15, b >> 4)]
    return CHIP_ARR[syms].reshape(-1)


def modulate(chips, sps=4):
    """O-QPSK, half-sine pulses of 2 chips, even chips on I, odd on Q delayed by one chip; sps samples per chip."""
    c = 2.0 * np.asarray(chips, float) - 1
    n = len(c) * sps + 2 * sps
    t = np.arange(2 * sps) / (2 * sps)
    p = np.sin(np.pi * t)
    i = np.zeros(n)
    q = np.zeros(n)
    for k in range(0, len(c) - 1, 2):
        i[k * sps:(k + 2) * sps] += c[k] * p
        q[(k + 1) * sps:(k + 3) * sps] += c[k + 1] * p
    return (i + 1j * q).astype(np.complex64)


def _disc(y, sps):
    """Discriminator averaged over each chip: one value per chip boundary-to-boundary interval."""
    d = np.angle(y[1:] * np.conj(y[:-1]))
    return np.convolve(d, np.ones(sps), "valid")


_TPL = None


def _templates(sps=4):
    """[16, 31] expected discriminator signs inside each symbol (chip intervals 1..31), by modulation."""
    global _TPL
    if _TPL is None or _TPL[0] != sps:
        rows = []
        for s in range(16):
            ch = np.r_[CHIP_ARR[0], CHIP_ARR[s], CHIP_ARR[0]]               # neighbours: any symbol works inside
            d = _disc(modulate(ch, sps), sps)
            base = 32 * sps + sps                                        # start of chip interval 1 of the symbol
            v = d[base + np.arange(31) * sps - 1]
            rows.append(np.sign(v))
        _TPL = (sps, np.array(rows))
    return _TPL[1]


def demod(y, ofs, max_bytes=133):
    """One burst at 4 samples / chip (ofs = 8 MHz) centred on the carrier -> dict or None."""
    sps = int(round(ofs / CHIP_RATE))
    if sps != 4:
        return None
    T = _templates(sps)
    d = _disc(y, sps)
    n_sym = (len(d) - sps) // (32 * sps) - 1
    if n_sym < 12:
        return None
    # timing: the symbol-0 template against every offset in the first ~4 symbols (the preamble); keep the best
    span = max(1, min(len(d) - 32 * sps * 10, 32 * sps * 4))
    idx = np.arange(span)[:, None] + np.arange(31)[None, :] * sps
    v = np.sign(d[idx]) @ T[0]
    off0 = int(np.argmax(v))
    if v[off0] < 20:
        return None

    def sym_at(k):
        idx = off0 + k * 32 * sps + np.arange(31) * sps
        if idx[-1] >= len(d):
            return None, 0.0
        sc = T @ d[idx]
        j = int(np.argmax(sc))
        return j, float(sc[j])

    # find SFD (symbols 7, A) after a run of zero symbols, within the first 12 symbols
    k = 0
    while k < 12:
        a, _ = sym_at(k)
        b, _ = sym_at(k + 1)
        if a == 7 and b == 0xA:
            break
        k += 1
    else:
        return None
    k += 2
    syms = []
    plen = None
    while True:
        s, _ = sym_at(k)
        if s is None:
            break
        syms.append(s)
        k += 1
        if len(syms) == 2:
            plen = (syms[0] | syms[1] << 4) & 0x7F
            if plen < 5:
                return None
        if plen is not None and len(syms) >= 2 + 2 * plen:
            break
    if plen is None or len(syms) < 2 + 2 * plen:
        return None
    by = bytes((syms[2 + 2 * i] | syms[3 + 2 * i] << 4) for i in range(plen))
    fcs = by[-2] | by[-1] << 8
    ok = crc16_kermit(by[:-2]) == fcs
    out = dict(protocol="IEEE 802.15.4 (O-QPSK 2.4 GHz)", length=plen, crc_ok=bool(ok), mpdu_hex=by.hex())
    if ok:
        out.update(parse_mac(by[:-2]))
    return out


def parse_mac(m: bytes) -> dict:
    """Frame control, sequence number, PAN IDs and addresses (little-endian on air)."""
    if len(m) < 3:
        return {}
    fcf = m[0] | m[1] << 8
    ftype, pan_comp = fcf & 7, (fcf >> 6) & 1
    dmode, smode = (fcf >> 10) & 3, (fcf >> 14) & 3
    out = dict(frame_type=FRAME_TYPES.get(ftype, ftype), seq=m[2], security=bool(fcf >> 3 & 1),
               ack_request=bool(fcf >> 5 & 1))
    i = 3

    def take(n):
        nonlocal i
        v = m[i:i + n][::-1].hex() if i + n <= len(m) else None
        i += n
        return v

    if dmode:
        out["dst_pan"] = take(2)
        out["dst_addr"] = take(2 if dmode == 2 else 8)
    if smode:
        if not (pan_comp and dmode):
            out["src_pan"] = take(2)
        out["src_addr"] = take(2 if smode == 2 else 8)
    out["payload_bytes"] = max(0, len(m) - i)
    ext = [a for a in (out.get("src_addr"), out.get("dst_addr")) if a and len(a) == 16]
    if any(a.startswith(DIGI_OUI) for a in ext):
        out["vendor"] = "Digi International (XBee)"
    return out
