"""
mlrs.py -- mLRS (olliw42) RC / telemetry link frames: passive identification (research, PLAN R3 / T6;
research/drone_transmitters/12_rc_links_more.md section 1). Synthetic validation only.

Air modes (LoRa, preamble 12, IQ never inverted, implicit header, 91-byte frames): 2.4 GHz SF5 / 6 / 7 at 812.5 kHz
(50 / 31 / 19 Hz), sub-GHz SF5 (31 Hz, SX126x), SF6 SX126x framing (19 Hz) and SF6 SX127x framing (19 Hz "7x").
Frame: bytes 0-1 FrameSyncWord (the per-binding 16-bit ID, also the FSK sync word), 2-6 status, 7-88 RC data /
payload, 89-90 CRC-16/MCRF4XX (reflected 0x1021, init 0xFFFF) over bytes 0-88. The byte order of the 16-bit fields and
whether an uplink's running CRC skips its embedded crc1 bytes are not published, so those variants are tried and the
one that passed is named. Only the SX127x-framed SF6 mode can be decoded here (SX126x SF5 / 6 and SX1280 long
interleaving are not); the others are identified by their air mode ("-like").
"""
from __future__ import annotations

FRAME_LEN = 91
# name: (band, sf, bw, cr_den, preamble, frame bytes, period us, decodable here)
MLRS_LORA = {
    "mlrs_50hz_2g4": ("2G4", 5, 812.5e3, 5, 12, FRAME_LEN, 20000, False),
    "mlrs_31hz_2g4": ("2G4", 6, 812.5e3, 5, 12, FRAME_LEN, 32000, False),
    "mlrs_19hz_2g4": ("2G4", 7, 812.5e3, 5, 12, FRAME_LEN, 53000, False),
    "mlrs_31hz_900": ("900", 5, 500e3, 5, 12, FRAME_LEN, 32000, False),
    "mlrs_19hz_900": ("900", 6, 500e3, 5, 12, FRAME_LEN, 53000, False),
    "mlrs_19hz_7x": ("900", 6, 500e3, 5, 12, FRAME_LEN, 53000, True),
}
_CRC1_POS = (13, 14)                   # uplink V1: header 7 B + 6 B of RC data, then crc1 (research file 12 §1.2)


def crc16_mcrf4xx(data: bytes, init=0xFFFF) -> int:
    c = init
    for b in data:
        c ^= b
        for _ in range(8):
            c = (c >> 1) ^ 0x8408 if c & 1 else c >> 1
    return c


def parse(frame: bytes):
    """A 91-byte frame -> dict (FrameSyncWord, CRC convention) when one CRC convention verifies it, else None."""
    if len(frame) != FRAME_LEN:
        return None
    body = frame[:89]
    variants = (("all bytes 0-88", body), ("crc1 bytes skipped", body[:_CRC1_POS[0]] + body[_CRC1_POS[1] + 1:]))
    for cov, data in variants:
        c = crc16_mcrf4xx(data)
        for order in ("little", "big"):
            if c == int.from_bytes(frame[89:91], order):
                return dict(protocol="mLRS", crc_ok=True, frame_sync_word=f"0x{int.from_bytes(frame[0:2], order):04X}",
                            frame_bytes=FRAME_LEN, convention=dict(crc_coverage=cov, byte_order=order),
                            status_hex=frame[2:7].hex())
    return None


def build_frame(sync_word: int, body: bytes) -> bytes:
    """Synthesis / tests: sync word (LE) + 87 body bytes + CRC (LE, over bytes 0-88)."""
    f = sync_word.to_bytes(2, "little") + bytes(body[:87]).ljust(87, b"\0")
    return f + crc16_mcrf4xx(f).to_bytes(2, "little")
