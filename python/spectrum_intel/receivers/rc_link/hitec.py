"""
hitec.py -- Hitec Optima / Minima 2.4 GHz RC link (CC2500 MSK) and Futaba S-FHSS (CC2500 2-FSK): passive
identification and decode (research, PLAN G12; research/drone_transmitters/17_hott_hitec_r9.md sections 0, 2, 3).
Synthetic validation only.

Hitec: MSK 249.939 kBd, f = 2405.59985 + 0.299927 CHANNR MHz (CHANNR even, 0..140; 21 per TX, one per 22.5 ms frame,
4 identical packets per frame), 24-byte preamble, sync = TX ID bytes id3 id2 sent twice (bind: 0xD391 twice),
variable length (byte 0), PN9 whitening, CRC-16 0x8005 / 0xFFFF over length + data.
- data (length 0x1A): TX ID (3 bytes), 9 channels u16 BE, 0xEE, repeat flag 0x80 / 0x40 / 0x20 / 0x10;
- bind (0x16): TX ID, bind phase 0x72.. (Minima 0x82 / 0x83), hop table / 2;
- receiver telemetry (0x0C): TX ID, frame number, data (frame 0x00 / 0x11: RX battery = word / 28 V).
The sync word must equal the TX ID's bytes 3, 2 (data and telemetry): a second, independent check besides the CRC.

S-FHSS: 2-FSK 128.143 kBd +-38.086 kHz, f = 2403.9987 + 1.49963 k MHz (k = 0..29), 4-byte preamble, sync 0xD391
(16 or 32 bits), 13 bytes, no whitening, CRC-16 0x8005 over the 13 bytes in CC2400 mode (init not documented: 0xFFFF
and 0x0000 are tried; the one that passed is reported). Fields: header 0x80..0x82, 16-bit TX ID, model-hint bytes 3-4,
hop index k (byte 5, checked against the measured channel), 4 x 12-bit channels, hop code, command bits.
A packet is reported only when the CRC passes (and, for Hitec, the sync / TX ID check). Nothing is transmitted.
"""
from __future__ import annotations

import numpy as np

from . import cc2500 as CC
from . import radiolink as RL
from . import sik as SIK

RATE = 249.939e3
DEV = RATE / 4
F0, STEP = 2405.59985e6, 0.299927e6
BIND_SYNC = 0xD391
LENGTHS = {0x1A: "data", 0x16: "bind", 0x0C: "telemetry", 0x0A: "bind reply"}

SF_RATE, SF_DEV = 128.143e3, 38.086e3
SF_F0, SF_STEP = 2403.9987e6, 1.49963e6
SF_LEN = 13
SF_MODEL_HINTS = {bytes.fromhex("4207"): "T8J", bytes.fromhex("0f09"): "T10J", bytes.fromhex("9a06"): "TM-FH"}


# ---------------------------------------------------------------- Hitec
def hitec_frame(tx_id: bytes, channels, repeat=0) -> bytes:
    """Synthesis / tests: the FIFO content (length + 26 bytes) of a data packet."""
    p = bytearray(27)
    p[0] = 0x1A
    p[1:4] = bytes(tx_id)
    for i, v in enumerate(list(channels)[:9]):
        p[4 + 2 * i:6 + 2 * i] = int(v).to_bytes(2, "big")
    p[22], p[23] = 0xEE, (0x80, 0x40, 0x20, 0x10)[repeat % 4]
    return bytes(p)


def hitec_air(sync, fifo: bytes) -> bytes:
    """Synthesis / tests: bytes after the preamble: sync twice, whitened FIFO + CRC."""
    c = CC.crc16_cc2500(fifo)
    return bytes([sync >> 8, sync & 0xFF] * 2) + RL.whiten(fifo + bytes([c >> 8, c & 0xFF]))


def _hitec_parse(sync, raw):
    L = raw[0]
    if L not in LENGTHS or len(raw) < L + 3:
        return None
    fifo = raw[:L + 1]
    if CC.crc16_cc2500(fifo) != (raw[L + 1] << 8 | raw[L + 2]):
        return None
    kind = LENGTHS[L]
    tid = fifo[1:4]
    if sync != BIND_SYNC and sync != (tid[2] << 8 | tid[1]):
        return None                                           # the sync word is the TX ID's bytes 3, 2
    out = dict(protocol="Hitec Optima / Minima", packet_name=kind, tx_id=tid.hex(), length=L, crc_ok=True,
               sync_matches_tx_id=sync != BIND_SYNC, packet_hex=fifo.hex())
    if kind == "data":
        out.update(channels=[fifo[4 + 2 * i] << 8 | fifo[5 + 2 * i] for i in range(9)],
                   repeat_flag=f"0x{fifo[23]:02X}", marker_ok=fifo[22] == 0xEE)
    elif kind == "bind":
        ph = fifo[4]
        out.update(bind_phase=f"0x{ph:02X}", variant="Minima" if ph & 0xF0 == 0x80 else "Optima")
        if ph in (0x72, 0x82):
            out["hop_channels_part"] = [2 * v for v in fifo[5:19]]
        elif ph in (0x73, 0x83):
            out["hop_channels_part"] = [2 * v for v in fifo[5:12]]
    elif kind == "telemetry":
        out.update(protocol="Hitec receiver telemetry", frame=f"0x{fifo[5]:02X}")
        if fifo[5] in (0x00, 0x11):
            out["rx_battery_v"] = round((fifo[6] << 8 | fifo[7]) / 28, 2)
    return out


def hitec_try_decode(y, ofs, f_hz=None, rate=RATE):
    """One burst at ~8 samples/symbol centred on the carrier -> dict or None. The sync word (= TX ID) is not known
    in advance: any 16-bit word sent twice after the preamble is a candidate; the CRC and the ID check decide."""
    for ph, bits in SIK.demod_bits_mid(y, ofs, rate):
        for inv in (False, True):
            b = 1 - bits if inv else bits
            for pos in _doubled_words(b):
                w = int("".join(map(str, b[pos:pos + 16])), 2)
                raw = RL.whiten(CC._bytes(b[pos + 32:pos + 32 + 8 * (0x1A + 3)]))
                if not raw:
                    continue
                d = _hitec_parse(w, raw)
                if d is None:
                    continue
                d.update(sync=f"0x{w:04X}", polarity_inverted=bool(inv))
                if f_hz is not None:
                    n = (f_hz - F0) / STEP
                    d["channr"] = int(round(n))
                return d
    return None


def _doubled_words(b, min_pre=16):
    """Candidate sync positions: a 16-bit word sent twice, starting near the end of an alternating preamble run of
    >= min_pre bits (the sync's first bits may continue the alternation). The preamble itself repeats every 16 bits,
    so a fully alternating window is not a candidate."""
    if len(b) < min_pre + 32:
        return []
    alt = np.flatnonzero(np.diff(b.astype(np.int8)) == 0)             # breaks of the alternation
    edges = np.r_[-1, alt, len(b) - 1]
    out = []
    for a0, a1 in zip(edges[:-1], edges[1:]):
        if a1 - a0 < min_pre:
            continue
        for p in range(max(int(a0) + 1, int(a1) - 14), min(int(a1) + 3, len(b) - 32)):
            w = b[p:p + 32]
            if np.array_equal(w[:16], w[16:]) and np.any(np.diff(w.astype(np.int8)) == 0):
                out.append(p)
    return out


# ---------------------------------------------------------------- Futaba S-FHSS
def crc16_cc2400(data: bytes, init=0xFFFF) -> int:
    """CRC-16 x^16 + x^15 + x^2 + 1, MSB first, over the data (CC2400-compatible mode; init per caller)."""
    c = init
    for v in data:
        for i in range(8):
            bit = (v >> (7 - i)) & 1
            c = ((c << 1) ^ 0x8005) & 0xFFFF if ((c >> 15) & 1) ^ bit else (c << 1) & 0xFFFF
    return c


def sfhss_packet(tx_id: int, k: int, channels, code: int, command=0x02, header=0x81, model=b"\0\0") -> bytes:
    """Synthesis / tests: the 13-byte packet."""
    c = [int(v) & 0xFFF for v in list(channels)[:4]]
    p = bytearray(13)
    p[0], p[1], p[2] = header, tx_id >> 8, tx_id & 0xFF
    p[3:5] = bytes(model)
    p[5] = (k << 3) & 0xF8 | (c[0] >> 9) & 7
    p[6] = (c[0] >> 1) & 0xFF
    p[7] = (c[0] & 1) << 7 | (c[1] >> 5) & 0x7F
    p[8] = (c[1] & 0x1F) << 3 | (c[2] >> 9) & 7
    p[9] = (c[2] >> 1) & 0xFF
    p[10] = (c[2] & 1) << 7 | (c[3] >> 5) & 0x7F
    p[11] = (c[3] & 0x1F) << 3 | (code >> 2) & 7
    p[12] = (code & 3) << 6 | command & 0x3F
    return bytes(p)


def sfhss_air(p: bytes, init=0xFFFF) -> bytes:
    c = crc16_cc2400(p, init)
    return bytes([0xD3, 0x91, 0xD3, 0x91]) + p + bytes([c >> 8, c & 0xFF])


def _sfhss_parse(p: bytes, crc: int):
    init = next((i for i in (0xFFFF, 0x0000) if crc16_cc2400(p, i) == crc), None)
    if init is None or p[0] not in (0x80, 0x81, 0x82):
        return None
    k = p[5] >> 3
    if k > 29:
        return None
    ch = [(p[5] & 7) << 9 | p[6] << 1 | p[7] >> 7, (p[7] & 0x7F) << 5 | p[8] >> 3,
          (p[8] & 7) << 9 | p[9] << 1 | p[10] >> 7, (p[10] & 0x7F) << 5 | p[11] >> 3]
    code = (p[11] & 7) << 2 | p[12] >> 6
    cmd = p[12] & 0x3F
    hint = SF_MODEL_HINTS.get(p[3:5])
    return dict(protocol="Futaba S-FHSS", packet_name="data", header=f"0x{p[0]:02X}", tx_id=p[1:3].hex(),
                model_bytes=p[3:5].hex(), model_hint=hint, hop_index=k, hop_code=code, channels=ch,
                command=f"0x{cmd:02X}", data_half=2 if cmd & 1 else 1, failsafe=bool(cmd & 4),
                channel_bank="CH5-8" if cmd & 8 else "CH1-4", crc_ok=True, crc_init=f"0x{init:04X}",
                packet_hex=p.hex())


def sfhss_try_decode(y, ofs, f_hz=None, rate=SF_RATE):
    """One burst at ~8 samples/bit centred on the carrier -> dict or None (CRC-verified; hop index checked)."""
    for ph, bits in SIK.demod_bits_mid(y, ofs, rate):
        for sync in ([0xD3, 0x91, 0xD3, 0x91], [0xD3, 0x91]):
            for pos, inv, err in CC._find_sync(bits, sync):
                b = bits[pos:pos + 8 * (SF_LEN + 2)]
                if len(b) < 8 * (SF_LEN + 2):
                    continue
                if inv:
                    b = 1 - b
                raw = CC._bytes(b)
                d = _sfhss_parse(raw[:SF_LEN], raw[SF_LEN] << 8 | raw[SF_LEN + 1])
                if d is None:
                    continue
                d.update(sync_bit_errors=err, polarity_inverted=bool(inv), sync_bits=8 * len(sync))
                if f_hz is not None:
                    kk = (f_hz - SF_F0) / SF_STEP
                    d["channel"] = int(round(kk))
                    d["hop_check"] = "match" if abs(kk - d["hop_index"]) <= 0.3 else "mismatch"
                return d
    return None
