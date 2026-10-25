"""
hubsan.py -- Hubsan toy / hobby drones (H107 / H107 Plus / H301 / H501 family, A7105): passive identification and
decode (research, PLAN G5; research/drone_transmitters/14_market_gap_survey.md section 3). Synthetic validation only.

PHY (A7105): 2-FSK 100 kb/s, +-186 kHz, no Gaussian filter, one fixed channel per transmitter out of 12 (2410-2465 MHz
in 5 MHz steps; the H107 sends every 5th packet 17.5 MHz higher), 10 ms period. On air: 4-byte alternating preamble,
the 32-bit ID code (bind: 0x55201041 / 0xAA201041; data: a per-session code from the bind handshake), 16 payload
bytes, CRC-16 in the data phase (convention unverified, not used). The payload's own checksum makes the 16 bytes sum
to 0 mod 256; the model is told by its constant bytes (H107: byte 10 = 0x64, bytes 11-14 = DB 04 26 79; H301:
0x18 / 0x5C / 0xF6; H107 Plus: 0x19 / 0x5C / 0x49; H501: flight-mode byte 0x41-0x43).
The ID code is read where the preamble ends (it is not known in advance), so a packet is reported only when the sum
and a known packet header verify (and, for data packets, a model's constant bytes). Nothing is transmitted.
"""
from __future__ import annotations

import numpy as np

from . import sik as SIK

RATE, DEV = 100e3, 186e3
BIND_IDS = {0x55201041: "bind (normal)", 0xAA201041: "bind (plus)"}
CHANNELS = tuple(range(0x14, 0x83, 10))             # A7105 channel numbers: f = 2400.001 + 0.5 n MHz
F0, STEP = 2400.000977e6, 0.5e6
LEN = 16
HEADERS = {0x20: "data", 0x40: "vtx", 0xE0: "telemetry", 0xE1: "telemetry", 0xE8: "keep-alive", 0xE9: "keep-alive",
           0x01: "bind 1", 0x03: "bind 3", 0x09: "bind 9"}


def checksum(p: bytes) -> int:
    return (256 - (sum(p[:15]) & 0xFF)) & 0xFF


def packet(header, body15: bytes) -> bytes:
    """Synthesis / tests: header + 14 bytes + checksum."""
    p = bytes([header]) + bytes(body15[:14]).ljust(14, b"\0")
    return p + bytes([checksum(p)])


def model_of(p: bytes):
    """Model candidates from a data packet's constant bytes."""
    if p[10] == 0x64 and p[11:15] == bytes.fromhex("db042679"):
        return ["H107 / H107L / H107C / H107D"]
    if p[10] == 0x18 and p[12] == 0x5C and p[14] == 0xF6:
        return ["H301"]
    if p[10] == 0x19 and p[12] == 0x5C and p[14] == 0x49:
        return ["H107 Plus"]
    if p[11] in (0x41, 0x42, 0x43) and p[3] == p[5] == p[7] == 0:
        return ["H501 / H122D / H502 class"]
    return []


def _alternating_runs(bits, min_len=24):
    d = np.flatnonzero(np.diff(bits.astype(np.int8)) == 0)          # breaks of the alternation
    edges = np.r_[-1, d, len(bits) - 1]
    return [(int(a) + 1, int(b) + 1) for a, b in zip(edges[:-1], edges[1:]) if b - a >= min_len]


def _bytes(bits):
    n = len(bits) // 8
    return bytes(int("".join(map(str, bits[8 * i:8 * i + 8])), 2) for i in range(n))


def try_decode(y, ofs, f_hz=None, rate=RATE):
    """One burst at ~8 samples/bit centred on the carrier -> dict or None."""
    for ph, bits in SIK.demod_bits_mid(y, ofs, rate):
        for inv in (False, True):
            b = 1 - bits if inv else bits
            for r0, r1 in _alternating_runs(b, 16):
                # the 32-bit ID code follows the 32-bit preamble; the burst window can start inside the preamble and an
                # ID whose first byte alternates (0x55 / 0xAA) extends the run, so every start from 8 bits into the
                # run to its end is tried (the checksum, header and model bytes decide)
                for s in range(r0 + 8, r1 + 1):
                    seg = b[s:s + 32 + 8 * LEN]
                    if len(seg) < 32 + 8 * LEN:
                        break
                    idc = int("".join(map(str, seg[:32])), 2)
                    p = _bytes(seg[32:])
                    if p[0] not in HEADERS or sum(p) & 0xFF:
                        continue
                    kind = HEADERS[p[0]]
                    models = model_of(p) if p[0] == 0x20 else []
                    if p[0] == 0x20 and not models:
                        continue                                    # data packets must show a model's constants
                    if kind.startswith("bind") and idc not in BIND_IDS:
                        continue
                    out = dict(protocol="Hubsan (A7105)", id_code=f"0x{idc:08X}",
                               id_kind=BIND_IDS.get(idc, "session"), packet_type=f"0x{p[0]:02X}", packet_name=kind,
                               checksum_ok=True, polarity_inverted=bool(inv), packet_hex=p.hex(), models=models)
                    if p[0] == 0x20:
                        out["sticks"] = dict(throttle=p[2], rudder=255 - p[4], elevator=255 - p[6], aileron=p[8])
                    if kind.startswith("bind"):
                        out["tx_id"] = p[2:6].hex()
                    if kind == "telemetry":
                        out["battery_v"] = round(p[13] / 10, 1)
                    if f_hz is not None:
                        n = int(round((f_hz - F0) / STEP))
                        out["channel"] = n
                        out["channel_check"] = "match" if n in CHANNELS or n - 0x23 in CHANNELS else "mismatch"
                    return out
    return None
