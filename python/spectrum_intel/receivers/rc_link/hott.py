"""
hott.py -- Graupner HoTT 2.4 GHz RC link (CC2500 MSK + FEC): passive identification and decode (research, PLAN G12;
research/drone_transmitters/17_hott_hitec_r9.md sections 0, 1, 5). Synthetic validation only.

PHY: MSK (CPFSK h = 0.5) 249.939 kBd, channels f = 2404.0264 + 1.014038 n MHz (n = 0..74), 4-byte preamble, sync
0xD391 sent twice (receiver telemetry: 0x2C6E twice in "Sync" mode), then the CC2500 FEC stream: [PN9-whitened data +
CRC-16 (0x8005, init 0xFFFF, over the data)] + trellis terminator 0x0B (x2 for an even count), rate-1/2 K = 4
convolutional code (g1 = 1 + D^2 + D^3 sent first, g0 = 1 + D + D^2 + D^3), 4 x 4 symbol interleaver per 32 coded
bits (TI DN504; the clean-room coder reproduces its test vector).
- TX packet (50 bytes, every 10 ms, 3.71 ms on air): hop-table tag (LE, one of 16 tables), hop index, channel-group
  byte (0x00 / 0x01; 0x40 / 0x41 failsafe), 12 channels u16 LE in 0.5 us, sensor byte, mode byte, 9 x 0xFF, 0x07,
  TX ID (5 bytes; bind: FA x 5), RX ID (5 bytes).
- RX telemetry (22 bytes): TX ID, RX ID, warnings, page, page data (page 0: RX voltage, temperature, RSSI, LQI).
A packet is reported only when the CRC passes. MSK polarity (the CC2500 inverts it) is tried both ways. Nothing is
transmitted or controlled.
"""
from __future__ import annotations

import numpy as np

from . import cc2500 as CC
from . import radiolink as RL
from . import sik as SIK

RATE = 249.939e3
DEV = RATE / 4
F0, STEP = 2404.0264e6, 1.014038e6
SYNC_TX, SYNC_TLM = 0xD391, 0x2C6E
LEN_TX, LEN_TLM = 50, 22
TABLE_TAGS = (0xC06B, 0xC34A, 0xDB24, 0x8E09, 0x272E, 0x217F, 0x155B, 0xEDE8, 0x1D31, 0x0986, 0x56F7, 0x6454,
              0xC42D, 0x01D2, 0xC253, 0x1180)
# trellis: symbol of (state = x[n-3] x[n-2] x[n-1], input x[n]) = T[8 x[n-3] + 4 x[n-2] + 2 x[n-1] + x[n]] (DN504)
_T = (0, 3, 1, 2, 3, 0, 2, 1, 3, 0, 2, 1, 0, 3, 1, 2)


def coded_len(n_data: int) -> int:
    """Coded bits on air for n_data bytes (data + 2 CRC + terminator to an even byte count)."""
    n = n_data + 2
    return 2 * 8 * (n + (1 if n % 2 else 2))


def fec_encode(data: bytes) -> np.ndarray:
    """Bytes (whitened data + CRC) -> interleaved coded bits (terminator added)."""
    b = bytes(data) + (b"\x0b" if len(data) % 2 else b"\x0b\x0b")
    st, sym = 0, []
    for v in b:
        for i in range(8):
            x = (v >> (7 - i)) & 1
            sym.append(_T[(st << 1) | x])
            st = ((st << 1) | x) & 7
    out = []
    for k in range(0, len(sym), 16):
        s = sym[k:k + 16]
        for j in range(16):
            v = s[4 * (3 - (j % 4)) + (3 - j // 4)]
            out += [(v >> 1) & 1, v & 1]
    return np.array(out, np.uint8)


def fec_decode(bits: np.ndarray, n_bytes: int):
    """Interleaved coded bits -> (n_bytes decoded bytes, path metric). Hard-decision Viterbi from state 0."""
    nb = len(bits) // 32
    sym = np.empty(16 * nb, np.int64)
    for k in range(nb):
        blk = bits[32 * k:32 * k + 32]
        rx = blk[0::2] * 2 + blk[1::2]
        for j in range(16):
            sym[16 * k + 4 * (3 - (j % 4)) + (3 - j // 4)] = rx[j]
    n = min(len(sym), 8 * n_bytes + 3)
    INF = 1 << 30
    pm = np.full(8, INF, np.int64)
    pm[0] = 0
    prev = np.zeros((n, 8), np.int8)
    # branch: from state s with input x -> next (s << 1 | x) & 7, symbol _T[(s << 1) | x]
    T = np.array(_T, np.int64)
    pc = np.array([0, 1, 1, 2], np.int64)                        # popcount of a 2-bit symbol difference
    s_from = np.arange(16) >> 1
    nxt = np.arange(16) & 7
    for t in range(n):
        cand = pm[s_from] + pc[T ^ sym[t]]
        new = np.full(8, INF, np.int64)
        arg = np.zeros(8, np.int8)
        for k in range(16):                                      # 2 candidates per next state
            if cand[k] < new[nxt[k]]:
                new[nxt[k]], arg[nxt[k]] = cand[k], k
        pm = new
        prev[t] = arg
    st = int(np.argmin(pm))
    metric = int(pm[st])
    out = np.zeros(n, np.uint8)
    for t in range(n - 1, -1, -1):
        k = int(prev[t, st])
        out[t] = k & 1
        st = k >> 1
    return CC._bytes(out[:8 * n_bytes]), metric


def tx_packet(table, index, channels_halfus, tx_id: bytes, rx_id: bytes, group=0x00, sensor=0x80, mode=0x02) -> bytes:
    """Synthesis / tests: the 50-byte TX packet."""
    p = bytearray(LEN_TX)
    p[0:2] = TABLE_TAGS[table].to_bytes(2, "little")
    p[2], p[3] = index, group
    for i, v in enumerate(list(channels_halfus)[:12]):
        p[4 + 2 * i:6 + 2 * i] = int(v).to_bytes(2, "little")
    p[28], p[29] = sensor, mode
    p[30:39] = b"\xff" * 9
    p[39] = 0x07
    p[40:45], p[45:50] = bytes(tx_id), bytes(rx_id)
    return bytes(p)


def air_bits(sync, data: bytes) -> np.ndarray:
    """Synthesis / tests: bits after the preamble: sync twice, then the FEC stream of whitened data + CRC."""
    c = CC.crc16_cc2500(data)
    s = [(sync >> (15 - i)) & 1 for i in range(16)] * 2
    return np.concatenate([np.array(s, np.uint8), fec_encode(RL.whiten(bytes(data) + bytes([c >> 8, c & 0xFF])))])


def _parse_tx(p: bytes):
    tag = p[0] | p[1] << 8
    tid, rid = p[40:45], p[45:50]
    bind = tid == b"\xfa" * 5
    grp = p[3]
    ch = [p[4 + 2 * i] | p[5 + 2 * i] << 8 for i in range(12)]
    out = dict(protocol="Graupner HoTT", packet_name="bind" if bind else ("failsafe" if grp & 0x40 else "data"),
               hop_table_tag=f"0x{tag:04X}", hop_table=TABLE_TAGS.index(tag) if tag in TABLE_TAGS else None,
               hop_index=p[2], channel_group="CH1-8 + CH13-16" if grp & 1 else "CH1-8 + CH9-12",
               channels_us=[round(v / 2, 1) for v in ch], sensor_byte=f"0x{p[28]:02X}", mode_byte=f"0x{p[29]:02X}",
               tx_id=(rid if bind else tid).hex(), rx_id=None if bind else rid.hex(), crc_ok=True)
    if p[30:39] != b"\xff" * 9 or p[39] != 0x07:
        out["note"] = "constant bytes 30-39 differ from the documented FF x 9, 07"
    return out


def _parse_tlm(p: bytes):
    out = dict(protocol="Graupner HoTT receiver telemetry", packet_name="telemetry", tx_id=p[0:5].hex(),
               rx_id=p[5:10].hex(), warnings=p[10], page=p[11], crc_ok=True)
    if p[11] == 0:
        out.update(rx_voltage_v=round(p[13] / 10, 1), temperature_c=p[14] - 20, rx_rssi_raw=p[15], lqi_pct=p[16])
    return out


def try_decode(y, ofs, f_hz=None, rate=RATE):
    """One burst at ~8 samples/symbol centred on the carrier -> dict or None (CRC-verified only)."""
    for ph, bits in SIK.demod_bits_mid(y, ofs, rate):
        for sync, n, parse in ((SYNC_TX, LEN_TX, _parse_tx), (SYNC_TLM, LEN_TLM, _parse_tlm)):
            nc = coded_len(n)
            for pos, inv, err in CC._find_sync(bits, [sync >> 8, sync & 0xFF] * 2):
                b = bits[pos:pos + nc]
                if len(b) < nc:
                    continue
                if inv:
                    b = 1 - b
                dec, metric = fec_decode(b, n + 2)
                raw = RL.whiten(dec)
                if CC.crc16_cc2500(raw[:n]) != (raw[n] << 8 | raw[n + 1]):
                    continue
                d = parse(raw[:n])
                d.update(sync=f"0x{sync:04X}", sync_bit_errors=err, polarity_inverted=bool(inv),
                         fec_corrected_bits=metric)
                if f_hz is not None:
                    k = (f_hz - F0) / STEP
                    d["channel"] = int(round(k))
                return d
    return None
