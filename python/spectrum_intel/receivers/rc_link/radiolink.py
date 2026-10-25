"""
radiolink.py -- RadioLink 2.4 GHz RC links (CC2500; AT9S / AT10 / T8FB / T8S air, RC4GS / RC6GS / RC8X surface,
DumboRC, RC4G): passive identification and decode (research, PLAN G4; research/drone_transmitters/14_market_gap_survey.md
section 3). Synthetic validation only.

PHY (CC2500, 26 MHz crystal): GFSK 37.99 kBd, +-38.09 kHz, channels f = 2401.99994 MHz + n x 399.902 kHz; preamble
4 x 0xAA, the 16-bit sync word sent twice (SYNC_MODE 30/32), then PN9-whitened [length, payload, CRC-16]
(0x8005, init 0xFFFF over length + payload, before whitening). Sync words: 0xABCD surface / air, 0xBADC DumboRC,
0xABA5 RC4G.
- Surface / air / DumboRC data (length 33): header (0x01 surface, 0x21 air; | 0x02 telemetry request), TX ID (4 bytes),
  16 channels x 11 bits LSB first, hop counter `pseudo` (LE), 2 zero bytes, replaced-slot index, checksum = sum of
  bytes 1..32. 16 hop channels from the TX ID (two LCG shuffles of 12 i + inc); the packet goes out on
  hop[pseudo & 15] (one slot may be replaced by a random channel 192 + inc, signalled in the packet). 19.9 ms period.
- Receiver telemetry (length 15): type, TX ID, echoed header, RSSI, battery; checksum = sum of bytes 1..14.
- RC4G (sync 0xABA5, length 15): 3-byte ID, 4 channels, next channel; CRC only. 14.46 ms period.
A packet is reported only when the CC2500 CRC passes (and, for the 33-byte data packet, the checksum); the channel the
ID's hop table predicts for the packet is compared with the measured frequency. Nothing is transmitted or controlled.
"""
from __future__ import annotations

from . import cc2500 as CC
from . import sik as SIK

RATE, DEV = 37.987e3, 38.09e3
F0, STEP = 2401.99994e6, 399.902e3
SYNCS = {0xABCD: "surface / air", 0xBADC: "DumboRC", 0xABA5: "RC4G"}
PERIOD_US = {"data": 19900, "RC4G": 14460}
MAX_LEN = 0x21


def pn9(n):
    """CC2500 PN9 data-whitening bytes (x^9 + x^5 + 1, seed 0x1FF): FF E1 1D 9A ED 85 ..."""
    st, out = 0x1FF, []
    for _ in range(n):
        out.append(st & 0xFF)
        for _ in range(8):
            st = (st >> 1) | (((st ^ (st >> 5)) & 1) << 8)
    return out


_PN9 = pn9(64)


def whiten(b: bytes) -> bytes:
    return bytes(v ^ _PN9[i] for i, v in enumerate(b))


def _prng(r):
    return (0xA5E2A705 * r + 0x754DB79B) & 0xFFFFFFFF


def _shuffle(hf, seed):
    r2, r1 = seed, seed
    for _ in range(31):
        r1 = _prng(r1)
    for i in range(16):
        k = (((r2 >> 16) + r2 + (r1 >> 16) + r1) & 0xFF) & 0x0F
        r2, r1 = _prng(r2), _prng(r1)
        hf[k], hf[i] = hf[i], hf[k]


def hop_table(tx_id: bytes):
    """The 16 hop channel numbers of a TX ID (4 bytes, as sent), and the replacement channel."""
    inc = 3 * (tx_id[0] & 3)
    hf = [12 * i + inc for i in range(16)]
    for a, b in ((0, 1), (2, 3)):
        _shuffle(hf, ((tx_id[a] + 256 * tx_id[b]) * 0xF65EF9F9 + 0x2EDDF6CA) & 0xFFFFFFFF)
    return hf, 192 + (6 if inc == 9 else inc)


def channel_mhz(n):
    return (F0 + n * STEP) / 1e6


def data_packet(header, tx_id, channels, pseudo, slot=16):
    """Synthesis / tests: the 34-byte FIFO content (length + 33 bytes) of a surface / air data packet."""
    v = 0
    for i, c in enumerate(channels[:16]):
        v |= (int(c) & 0x7FF) << (11 * i)
    p = bytes([header]) + bytes(tx_id) + v.to_bytes(22, "little") + pseudo.to_bytes(2, "little") + b"\0\0" + \
        bytes([slot])
    return bytes([len(p) + 1]) + p + bytes([sum(p) & 0xFF])


def air_frame(sync, fifo: bytes) -> bytes:
    """Synthesis / tests: bytes after the preamble: sync twice, whitened FIFO + CRC."""
    c = CC.crc16_cc2500(fifo)
    return bytes([sync >> 8, sync & 0xFF] * 2) + whiten(fifo + bytes([c >> 8, c & 0xFF]))


def _parse(sync, raw):
    """Dewhitened bytes after the sync -> dict when the CRC (and checksum) verify, else None."""
    L = raw[0]
    if L not in (0x21, 0x0F) or len(raw) < L + 3:
        return None
    fifo = raw[:L + 1]
    if CC.crc16_cc2500(fifo) != (raw[L + 1] << 8 | raw[L + 2]):
        return None
    sub = SYNCS[sync]
    out = dict(sync=f"0x{sync:04X}", length=L, crc_ok=True, packet_hex=fifo.hex())
    if sync == 0xABA5:
        if L != 0x0F:
            return None
        ch = [fifo[4 + k] | (((fifo[8 + k // 2] >> (4 * (k % 2))) & 0x0F) << 8) for k in range(4)]
        out.update(protocol="RadioLink RC4G", sub_protocol="RC4G", packet_name="data", tx_id=fifo[1:4].hex(),
                   channels=ch, next_channel=fifo[15], checksum_ok=None)
        return out
    if L == 0x21:
        if sum(fifo[1:33]) & 0xFF != fifo[33]:
            return None
        hdr = fifo[1]
        if hdr & ~0x02 not in (0x01, 0x21):
            return None
        v = int.from_bytes(fifo[6:28], "little")
        pseudo = fifo[28] | fifo[29] << 8
        kind = "air" if hdr & 0x20 else ("DumboRC" if sync == 0xBADC else "surface")
        out.update(protocol=f"RadioLink ({kind})", sub_protocol=kind, packet_name="data", header=f"0x{hdr:02X}",
                   telemetry_request=bool(hdr & 0x02), tx_id=fifo[2:6].hex(),
                   channels=[(v >> (11 * i)) & 0x7FF for i in range(16)], pseudo=pseudo, hop_slot=pseudo & 0x0F,
                   replaced_slot=fifo[32], checksum_ok=True)
        hf, rep = hop_table(fifo[2:6])
        if fifo[32] < 16:
            hf[fifo[32]] = rep
        out["hop_channels"] = hf
        out["expected_channel"] = hf[pseudo & 0x0F]
        return out
    if sum(fifo[1:15]) & 0xFF != fifo[15]:
        return None
    out.update(protocol=f"RadioLink receiver telemetry ({sub})", sub_protocol=sub, packet_name="telemetry",
               telemetry_type=fifo[1], tx_id=fifo[2:6].hex(), echoed_header=f"0x{fifo[6]:02X}",
               rssi_raw=fifo[7], checksum_ok=True)
    return out


def _find_sync(bits, sync):
    """[(pos after the doubled sync, inverted, bit errors)] for the 32-bit doubled sync word (<= 2 errors)."""
    return CC._find_sync(bits, [sync >> 8, sync & 0xFF] * 2)


def try_decode(y, ofs, f_hz=None, rate=RATE):
    """One burst at ~8 samples/bit centred on the FSK carrier -> dict or None (CRC-verified only). f_hz: the burst's
    carrier, compared with the hop table's prediction."""
    for ph, bits in SIK.demod_bits(y, ofs, rate):
        for sync in SYNCS:
            for pos, inv, err in _find_sync(bits, sync):
                b = bits[pos:pos + 8 * (MAX_LEN + 3)]
                if inv:
                    b = 1 - b
                raw = whiten(CC._bytes(b))
                d = _parse(sync, raw)
                if d is None:
                    continue
                d.update(sync_bit_errors=err, polarity_inverted=bool(inv))
                if f_hz is not None:
                    n = (f_hz - F0) / STEP
                    d["channel"] = int(round(n))
                    if "expected_channel" in d:
                        d["hop_check"] = "match" if abs(n - d["expected_channel"]) <= 1.5 else "mismatch"
                return d
    return None
