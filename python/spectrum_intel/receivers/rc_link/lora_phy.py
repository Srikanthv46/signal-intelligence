"""LoRa payload decoding from hard symbol values (explicit-header packets).

    decode(raw_syms, sf, ldro=None) -> dict | None

raw_syms are the de-chirped peak positions of the payload symbols, relative to
the preamble tone (what detector._walk_packet measures), in transmit order and
starting with the first symbol after the 2.25-symbol SFD.

The chain is the Semtech SX127x/SX126x one as implemented by gr-lora_sdr
(Tapparel et al., EPFL), which was used as the reference:

    symbol - 1 (mod 2^SF)                       fft_demod
    header block / LDRO: keep SF-2 bits (/4)    reduced rate: the 2 LSBs are a parity bit and 0
    Gray: g = v ^ (v >> 1)                      gray_mapping
    diagonal deinterleave, block of 4+CR syms   deinterleaver
    Hamming (8,4)/(7,4) syndrome, (6,4)/(5,4)   hamming_dec
    header: length, CR, CRC flag, 5-bit chksum  header_decoder
    bytes = low nibble first, XOR whitening     dewhitening (LFSR x^8+x^6+x^5+x^4+1, seed 0xFF)
    CRC-16/CCITT(0) over len-2 bytes XOR the    crc_verif
    last two bytes == the two trailing bytes

The first block (8 symbols) is always coded at CR 4/8 and reduced rate; it holds
the 5 header nibbles and the first SF-7 payload nibbles. Implicit-header packets
(ExpressLRS, most RC links) carry no length/CR and are not decoded here.

SF5 (SX126x), decode_sf5_implicit: measured on the bench (testdata/live/LIVE_LORA_SF5_890MHz,
2026-10-19). The SFD is followed by 2 extra symbols, and the first block is CR 4/8 at the full
SF (no reduced rate: 4*SF bits, 5 nibbles), as the SX126x time-on-air term
ceil((8PL + 16CRC - 4SF + 20H) / 4SF) and its 6.25 (not 4.25) symbols of overhead imply.
With an implicit header these 5 nibbles are payload; length and CR are found by the CRC.
Explicit-header SF5 and SF6 are not verified and not decoded.
"""
from __future__ import annotations

import numpy as np

from spectrum_intel.core import mavlink


def whitening_sequence(n):
    """LoRa whitening bytes: 0xFF, 0xFE, 0xFC, 0xF8, 0xF0, 0xE1, 0xC2, 0x85, ..."""
    out, r = [], 0xFF
    for _ in range(n):
        out.append(r)
        fb = ((r >> 7) ^ (r >> 5) ^ (r >> 4) ^ (r >> 3)) & 1
        r = ((r << 1) | fb) & 0xFF
    return np.array(out, np.uint8)


def crc16(data):
    """CRC-16/CCITT, polynomial 0x1021, initial value 0, MSB first, no final XOR."""
    crc = 0
    for b in data:
        for j in range(7, -1, -1):
            fb = ((crc >> 15) & 1) ^ ((int(b) >> j) & 1)
            crc = (crc << 1) & 0xFFFF
            if fb:
                crc ^= 0x1021
    return crc


def _bits(v, n):
    return [(int(v) >> (n - 1 - i)) & 1 for i in range(n)]          # MSB first


def _deinterleave(block, sf_app, cw_len):
    """cw_len symbols of sf_app bits -> sf_app codewords of cw_len bits (MSB first)."""
    inter = [_bits(s, sf_app) for s in block]
    de = [[0] * cw_len for _ in range(sf_app)]
    for i in range(cw_len):
        for j in range(sf_app):
            de[(i - j - 1) % sf_app][i] = inter[i][j]
    return de


def _hamming(cw, cw_len):
    """One codeword (list of cw_len bits, MSB first) -> (nibble, corrected, parity_error)."""
    d = [cw[3], cw[2], cw[1], cw[0]]                                 # data bits, MSB first
    corrected = bad = False
    if cw_len >= 7:
        s0 = cw[0] ^ cw[1] ^ cw[2] ^ cw[4]
        s1 = cw[1] ^ cw[2] ^ cw[3] ^ cw[5]
        s2 = cw[0] ^ cw[1] ^ cw[3] ^ cw[6]
        syn = s0 | (s1 << 1) | (s2 << 2)
        even = sum(cw) % 2 == 0
        if syn and (cw_len == 7 or not even):                          # (8,4): a single error flips overall parity
            flip = {5: 3, 7: 2, 3: 1, 6: 0}.get(syn)
            if flip is not None:
                d[flip] ^= 1
                corrected = True
            else:
                bad = True                                           # syndrome hits a parity bit
        elif syn:
            bad = True                                               # (8,4): double error, detected only
    elif cw_len == 6:
        bad = bool(cw[4] ^ cw[0] ^ cw[1] ^ cw[2]) or bool(cw[5] ^ cw[1] ^ cw[2] ^ cw[3])
    elif cw_len == 5:
        bad = bool(sum(cw) % 2)
    return (d[0] << 3) | (d[1] << 2) | (d[2] << 1) | d[3], corrected, bad


def header_checksum(n0, n1, n2):
    b = lambda v, k: (v >> k) & 1
    c4 = b(n0, 3) ^ b(n0, 2) ^ b(n0, 1) ^ b(n0, 0)
    c3 = b(n0, 3) ^ b(n1, 3) ^ b(n1, 2) ^ b(n1, 1) ^ b(n2, 0)
    c2 = b(n0, 2) ^ b(n1, 3) ^ b(n1, 0) ^ b(n2, 3) ^ b(n2, 1)
    c1 = b(n0, 1) ^ b(n1, 2) ^ b(n1, 0) ^ b(n2, 2) ^ b(n2, 1) ^ b(n2, 0)
    c0 = b(n0, 0) ^ b(n1, 1) ^ b(n2, 3) ^ b(n2, 2) ^ b(n2, 1) ^ b(n2, 0)
    return (c4 << 4) | (c3 << 3) | (c2 << 2) | (c1 << 1) | c0


def _block_nibbles(syms, sf, sf_app, cw_len):
    N = 2 ** sf
    v = (np.asarray(syms, np.int64) - 1) % N
    if sf_app < sf:
        v = v >> 2                                                   # drop parity bit and the 0 bit
    g = v ^ (v >> 1)
    out, n_corr, n_bad = [], 0, 0
    for cw in _deinterleave(g, sf_app, cw_len):
        nib, c, b = _hamming(cw, cw_len)
        out.append(nib); n_corr += c; n_bad += b
    return out, n_corr, n_bad


def n_symbols(sf, length, cr, crc, ldro):
    """Payload symbols after the SFD for an explicit-header packet (SX127x datasheet formula)."""
    num = 8 * length - 4 * sf + 28 + 16 * crc
    return 8 + max(int(np.ceil(num / (4 * (sf - 2 * ldro)))) * (cr + 4), 0)


def default_ldro(sf, bw):
    """Semtech rule: low-data-rate optimisation on when a symbol lasts more than 16 ms."""
    return 2 ** sf / bw > 16e-3


def decode(raw_syms, sf, bw=125e3, ldro=None):
    """Explicit-header decode. Returns None when the header block is missing or its
    checksum fails (then the packet is implicit-header or too noisy to say)."""
    if ldro is None:
        ldro = default_ldro(sf, bw)
    if sf < 7 or len(raw_syms) < 8:
        return None                                                  # SF5/6 (SX126x/SX128x) differ; not handled
    hdr, c_h, b_h = _block_nibbles(raw_syms[:8], sf, sf - 2, 8)
    n0, n1, n2, n3, n4 = hdr[:5]
    if ((n3 & 1) << 4 | n4) != header_checksum(n0, n1, n2):
        return None
    length, crc_on, cr = (n0 << 4) | n1, n2 & 1, n2 >> 1
    if not 1 <= cr <= 4 or length == 0:
        return None
    nib = list(hdr[5:])
    corr, bad = c_h, b_h
    need = 2 * length + 4 * crc_on
    sf_app = sf - 2 if ldro else sf
    cw_len = 4 + cr
    j = 8
    while len(nib) < need and j + cw_len <= len(raw_syms):
        n, c, b = _block_nibbles(raw_syms[j:j + cw_len], sf, sf_app, cw_len)
        nib += n; corr += c; bad += b
        j += cw_len
    out = dict(header_ok=True, length=length, cr=f"4/{4 + cr}", crc_present=bool(crc_on), ldro=bool(ldro),
               n_symbols=n_symbols(sf, length, cr, crc_on, ldro), symbols_used=j)
    if len(nib) < 2 * length:
        out.update(complete=False, reason=f"packet ends after {len(nib) // 2} of {length} bytes")
        return out
    raw = np.array([(nib[2 * i + 1] << 4) | nib[2 * i] for i in range(length)], np.uint8)
    pay = raw ^ whitening_sequence(length)
    out.update(complete=True, payload_hex=bytes(pay).hex(), fec_corrected=int(corr), fec_detected_errors=int(bad))
    if crc_on:
        if len(nib) < need:
            out.update(crc_ok=None, reason="CRC bytes missing")
            return out
        rx = (nib[2 * length + 1] << 4 | nib[2 * length]) | ((nib[2 * length + 3] << 4 | nib[2 * length + 2]) << 8)
        if length >= 2:
            c = crc16(pay[:-2]) ^ int(pay[-1]) ^ (int(pay[-2]) << 8)
        else:
            c = crc16(b"") ^ int(pay[-1])
        out["crc_ok"] = bool(c == rx)
    else:
        out["crc_ok"] = None
    return out


def n_symbols_implicit(sf, length, cr, crc, ldro):
    """Payload symbols after the SFD of an implicit-header SX127x packet (datasheet formula with IH = 1)."""
    num = 8 * length - 4 * sf + 8 + 16 * crc
    return 8 + max(int(np.ceil(num / (4 * (sf - 2 * ldro)))) * (cr + 4), 0)


def decode_implicit(raw_syms, sf, cr, length, ldro=False):
    """Implicit-header payload of known length and coding rate (SX127x, SF >= 6; ExpressLRS sub-GHz uses this with its
    own CRC and the LoRa CRC off): the first block is CR 4/8 at reduced rate (SF-2 bits / symbol) and holds the first
    SF-2 payload nibbles (no header), then blocks of 4+CR symbols. Returns dict(payload bytes, FEC counts) or None.
    Validated against gr-lora_sdr implicit-header packets (tests/receivers/gen_lora_vectors.py, PLAN R2c; SF6 PLAN R2e:
    the SX127x SF6 first block is also reduced-rate, research file 12 §6). SX126x SF5 / 6 framing differs: not here."""
    if sf < 6 or len(raw_syms) < 8 or not 1 <= cr <= 4:
        return None
    nib, corr, bad = _block_nibbles(raw_syms[:8], sf, sf - 2, 8)
    sf_app = sf - 2 if ldro else sf
    cw_len = 4 + cr
    j = 8
    while len(nib) < 2 * length and j + cw_len <= len(raw_syms):
        n, c, b = _block_nibbles(raw_syms[j:j + cw_len], sf, sf_app, cw_len)
        nib += n; corr += c; bad += b
        j += cw_len
    if len(nib) < 2 * length:
        return None
    raw = np.array([(nib[2 * i + 1] << 4) | nib[2 * i] for i in range(length)], np.uint8)
    pay = raw ^ whitening_sequence(length)
    return dict(implicit_header=True, length=length, cr=f"4/{cw_len}", payload=bytes(pay), fec_corrected=int(corr),
                fec_detected_errors=int(bad), symbols_used=j)


SF5_EXTRA_SYMS = 2          # symbols between the SFD and the first block (SX126x SF5/6)


def n_symbols_sf5(length, cr, crc=True, sf=5):
    """Symbols after the SFD of an implicit-header SF5 packet (SX126x formula, measured)."""
    return SF5_EXTRA_SYMS + 8 + max(int(np.ceil((8 * length + 16 * crc - 4 * sf) / (4 * sf))), 0) * (cr + 4)


def decode_sf5_implicit(raw_syms, n_pay, sf=5):
    """Implicit-header SF5 packet with a payload CRC-16. Every (CR, length) whose symbol count
    fits the measured burst (n_pay - 8 .. n_pay + 8: the walk ends early when the last block is
    all-zero nibbles, i.e. identical symbols) is tried; only a CRC pass without an uncorrectable
    FEC error is accepted (~36 hypotheses: a false pass ~5e-4). Passes with different payloads:
    None (ambiguous). Same payload, different CR (a last block of zero nibbles is a codeword at
    every CR): the symbol count nearest n_pay, with cr_ambiguous. None otherwise: no CRC (nothing
    verifies the length) or not decodable."""
    if sf != 5 or len(raw_syms) < SF5_EXTRA_SYMS + 8:
        return None
    first, c1, b1 = _block_nibbles(raw_syms[SF5_EXTRA_SYMS:SF5_EXTRA_SYMS + 8], sf, sf, 8)
    if b1:
        return None
    passes = []
    for cr in (1, 2, 3, 4):
        cw_len = 4 + cr
        nib, corr, bad, j = list(first), c1, 0, SF5_EXTRA_SYMS + 8
        for length in range(1, 256):
            ns = n_symbols_sf5(length, cr, True, sf)
            if ns > min(n_pay + 8, len(raw_syms)):
                break
            if ns < n_pay - 8:
                continue
            while j < ns:                                            # extend the nibble stream to ns symbols
                n, c, b = _block_nibbles(raw_syms[j:j + cw_len], sf, sf, cw_len)
                nib += n; corr += c; bad += b; j += cw_len
            if bad or len(nib) < 2 * length + 4:
                continue
            raw = np.array([(nib[2 * i + 1] << 4) | nib[2 * i] for i in range(length + 2)], np.uint8)
            pay = raw[:length] ^ whitening_sequence(length)
            rx = int(raw[length]) | int(raw[length + 1]) << 8
            c = crc16(pay[:-2]) ^ int(pay[-1]) ^ (int(pay[-2]) << 8) if length >= 2 else crc16(b"") ^ int(pay[-1])
            if c == rx:
                passes.append(dict(header_ok=False, implicit_header=True, length=length, cr=f"4/{cw_len}",
                                   crc_present=True, ldro=False, n_symbols=ns, symbols_used=ns, complete=True,
                                   payload_hex=bytes(pay).hex(), fec_corrected=int(corr), fec_detected_errors=0,
                                   crc_ok=True))
    if not passes or len({p["payload_hex"] for p in passes}) > 1:
        return None
    best = min(passes, key=lambda p: (abs(p["n_symbols"] - n_pay), p["fec_corrected"]))
    best["cr_ambiguous"] = len({p["cr"] for p in passes}) > 1
    return best


def sync_word(sync_vals):
    """Sync-word byte from the two sync symbols (each nibble is sent as 8 x nibble), or None."""
    if len(sync_vals) != 2:
        return None
    nib = []
    for v in sync_vals:
        q = int(round(v / 8))
        if abs(v - 8 * q) > 1 or not 0 <= q < 16:
            return None
        nib.append(q)
    return (nib[0] << 4) | nib[1]


# well-known sync words; the byte is a network separator, not a guarantee of the protocol
SYNC_WORDS = {0x12: "private (Semtech default, e.g. RadioLib / RadioHead / point-to-point)",
              0x34: "LoRaWAN public network", 0x2B: "Meshtastic", 0xF1: "FANET"}

# FANET (research, PLAN G10; research/drone_transmitters/16_uat_adsl_fanet_ogn.md section C): manufacturer IDs
FANET_MFR = {0x01: "Skytraxx", 0x03: "BitBroker", 0x04: "AirWhere", 0x05: "Windline", 0x06: "Burnair", 0x07: "SoftRF",
             0x08: "GXAirCom", 0x09: "Airtribune", 0x0A: "FLARM", 0x0B: "FlyBeeper", 0x0C: "Leaf Vario", 0x0D: "LYNKX",
             0x0E: "Syride", 0x0F: "FlyBrief", 0x10: "alfapilot", 0x11: "FANET+ (FLARM + FANET device)",
             0x20: "XC Tracer", 0xBA: "Burnair", 0xCB: "Cloudbuddy", 0xE0: "OGN Tracker", 0xE4: "4aviation",
             0xFA: "various", 0xFB: "Espressif base station", 0xFC: "unregistered", 0xFD: "unregistered"}
FANET_TYPES = {0: "ACK", 1: "tracking", 2: "name", 3: "message", 4: "service", 5: "landmarks", 6: "remote config",
               7: "ground tracking", 8: "HW info", 9: "thermal", 10: "HW info v2"}
FANET_AIRCRAFT = ("other", "paraglider", "hang glider", "balloon", "glider", "powered aircraft", "helicopter", "UAV")


def _s24(b):
    v = int.from_bytes(bytes(b), "little")
    return v - (1 << 24) if v & 0x800000 else v


def _fanet_scaled(v, mult):
    """7-bit value with bit 7 = x mult scale."""
    return (v & 0x7F) * (mult if v & 0x80 else 1)


def fanet_parse(payload, strict=False):
    """FANET MAC header + payload -> dict, or None. strict (sync word 0x12, shared with other private links): only a
    tracking packet of the exact length with a known manufacturer and a plausible position is accepted."""
    p = bytes(payload)
    if len(p) < 4:
        return None
    ext, fwd, typ, mfr = bool(p[0] & 0x80), bool(p[0] & 0x40), p[0] & 0x3F, p[1]
    dev = int.from_bytes(p[2:4], "little")
    if typ not in FANET_TYPES or mfr in (0x00, 0xFF) or dev in (0x0000, 0xFFFF):
        return None
    k = 4
    out = dict(protocol="FANET", message_type=FANET_TYPES[typ], forward=fwd, manufacturer_id=f"0x{mfr:02X}",
               manufacturer=FANET_MFR.get(mfr, "unknown"), device_id=f"{mfr:02X}:{dev:04X}")
    if ext:
        if len(p) < 5:
            return None
        e = p[4]
        k = 5
        out.update(ack=(e >> 6) & 3, unicast=bool(e & 0x20), signed=bool(e & 0x10))
        if e & 0x20:
            if len(p) < k + 3:
                return None
            out["destination"] = f"{p[k]:02X}:{int.from_bytes(p[k + 1:k + 3], 'little'):04X}"
            k += 3
        if e & 0x10:
            k += 4
    body = p[k:]
    if strict and (typ != 1 or mfr not in FANET_MFR or ext or len(body) not in (11, 12, 13)):
        return None
    if typ in (1, 7, 4):
        off = 1 if typ == 4 else 0
        if len(body) < off + 6:
            return None
        lat, lon = _s24(body[off:off + 3]) / 93206, _s24(body[off + 3:off + 6]) / 46603
        if abs(lat) > 90 or abs(lon) > 180:
            return None
        out.update(lat=round(lat, 5), lon=round(lon, 5))
    if typ == 1:
        if len(body) < 11:
            return None
        w = body[6] | body[7] << 8
        alt = (w & 0x7FF) * (4 if w & 0x800 else 1)
        clb = body[9] & 0x7F
        clb = clb - 128 if clb & 0x40 else clb
        out.update(aircraft_type=FANET_AIRCRAFT[(w >> 12) & 7], uav=((w >> 12) & 7) == 7,
                   online_tracking=bool(w & 0x8000), alt_m=alt, speed_kmh=_fanet_scaled(body[8], 5) / 2,
                   climb_ms=round(clb * (5 if body[9] & 0x80 else 1) / 10, 1), heading_deg=round(body[10] * 360 / 256, 1))
    elif typ == 7 and len(body) >= 7:
        out.update(ground_type=body[6] >> 4, online_tracking=bool(body[6] & 1))
    elif typ == 2:
        out["name"] = body.decode("utf-8", "replace")
    return out


def _le(b):
    return int.from_bytes(bytes(b), "little")


def identify(sw, payload):
    """Plain-text header fields of protocols known by their sync word. Only framing is
    read; LoRaWAN FRMPayload and Meshtastic data are encrypted and left as they are."""
    p = bytes(payload)
    if sw == 0x34 and len(p) >= 12:
        mtype, major = p[0] >> 5, p[0] & 3
        names = ["join_request", "join_accept", "unconfirmed_data_up", "unconfirmed_data_down",
                 "confirmed_data_up", "confirmed_data_down", "rejoin_request", "proprietary"]
        if major != 0:
            return None
        out = dict(protocol="LoRaWAN", mtype=names[mtype])
        if mtype == 0 and len(p) == 23:
            out.update(join_eui=p[1:9][::-1].hex(), dev_eui=p[9:17][::-1].hex(), dev_nonce=_le(p[17:19]))
        elif mtype in (2, 3, 4, 5):
            fopts = p[5] & 0x0F
            out.update(dev_addr=f"{_le(p[1:5]):08x}", adr=bool(p[5] & 0x80), ack=bool(p[5] & 0x20),
                       fcnt=_le(p[6:8]), fopts_len=fopts)
            if len(p) > 8 + fopts + 4:
                out["fport"] = p[8 + fopts]
        return out
    if sw in (0xF1, 0x12):
        f = fanet_parse(p, strict=sw == 0x12)
        if f:
            return f
    if sw == 0x2B and len(p) >= 16:
        to, frm, pid, fl = _le(p[0:4]), _le(p[4:8]), _le(p[8:12]), p[12]
        return dict(protocol="Meshtastic", to=("broadcast" if to == 0xFFFFFFFF else f"!{to:08x}"),
                    sender=f"!{frm:08x}", packet_id=pid, hop_limit=fl & 7, hop_start=fl >> 5,
                    want_ack=bool(fl & 8), via_mqtt=bool(fl & 16), channel_hash=p[13])
    return None


def mavlink_in(payload):
    """CRC-valid MAVLink frames in a LoRa payload -> {mavlink_frames, sysid, messages, ...} or None. A frame split
    over two LoRa packets (small module packet sizes) is not reassembled."""
    fr = [f for f in mavlink.scan(bytes(payload)) if f.get("crc_ok")]
    return {"mavlink_frames": len(fr), **mavlink.summary(fr)} if fr else None


def ascii_view(payload):
    return "".join(chr(b) if 32 <= b < 127 else "." for b in bytes(payload))
