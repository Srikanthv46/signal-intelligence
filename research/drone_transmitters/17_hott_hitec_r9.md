# 17. Graupner HoTT, Hitec Optima / Minima, Futaba S-FHSS, FrSky R9 (ACCST): decoder-level facts, plus CC2500 FEC

Research only (2026-10-22). This file adds detail to `12_rc_links_more.md` §2–3: register values, bit order, FEC, hop
tables and the steps of a passive decoder. Where the two files disagree, this file is the correction (see "Corrections"
at the end).

**Licence and sources.**
- The protocol facts come from DIY-Multiprotocol-TX-Module (MPM, pascallanger, **GPL-3**). It was read only as a
  specification / oracle: no code was copied. Re-implement from the facts below.
- Local copies: `src_gpl_reference/mpm/` (`HOTT_cc2500.ino`, `Hitec_cc2500.ino`, `Futaba_cc2500.ino` (= S-FHSS),
  `FrSkyR9_sx1276.ino`, `FrSkyDVX_common.ino`, `Protocols_Details.md`).
- Also read online (not copied into the repo): `SX1276_SPI.ino`, `iface_sx1276.h`, `Convert.ino`.
- Chip behaviour comes from TI and Semtech documents (CC2500 SWRS040C, DN504 SWRA113A, CC2400 SWRS042A, SX1276).
- Tags:
  - **[code]** = read in MPM;
  - **[spec]** = chip datasheet;
  - **[calc]** = computed here from register values with the datasheet formulas (26 MHz crystal for CC2500; 32 MHz /
    2^19 = 61.03515625 Hz FRF step for SX1276);
  - **[ver]** = checked here with a clean-room script against a published test vector;
  - **[U]** = uncertain or unverified.

**What MPM is.** MPM is a reverse-engineered transmitter that interoperates with the genuine receivers. A decoder built
from it therefore decodes what genuine receivers accept. Where the genuine transmitter may behave differently (ID
generation, hop-table choice), the item is marked [U].

---

## 0. CC2500 behaviour shared by HoTT, Hitec and S-FHSS [spec]

| Item | Fact | Source |
|---|---|---|
| Frequency | f = FREQ[23:0] · 26 MHz / 2^16 + CHANNR · Δch, with Δch = 26 MHz / 2^18 · (256 + CHANSPC_M) · 2^CHANSPC_E | CC2500 SWRS040C |
| Data rate | R = (256 + DRATE_M) · 2^DRATE_E · 26 MHz / 2^28 | idem |
| 2-FSK / GFSK deviation | f_dev = 26 MHz / 2^17 · (8 + DEVIATION_M) · 2^DEVIATION_E; '0' = −dev, '1' = +dev (Table 23) | idem p.33 |
| MSK | the whole burst (preamble, sync, payload) is MSK. "Identical to offset QPSK with half-sine shaping (data coding may differ)". "The MSK modulation format implemented in CC2500 **inverts the sync word and data** compared to e.g. signal generators". DEVIATN.DEVIATION_M sets the fraction of the symbol used for the phase change (7 = whole symbol); the mapping is only in SmartRF Studio | idem p.33–34, reg 0x15 |
| SYNC_MODE | 001 = 15/16, 010 = 16/16, **011 = 30/32: the 16-bit SYNC1:SYNC0 is sent twice** | idem Table 24 |
| Preamble | NUM_PREAMBLE (MDMCFG1[6:4]) → {2, 3, 4, 6, 8, 12, 16, 24} bytes of 1010… | idem |
| Whitening | PN9 (x^9 + x^5 + 1), register reset to all ones at the start of every packet. It covers everything after the sync word (length byte, data and CRC) and is applied **before** the FEC / interleaver stage | idem §15.1, §15.5 |
| CRC (CC2400_EN = 0) | CRC-16 poly 0x8005, init 0xFFFF, MSB first, over the length byte (if any) and the data, computed on **un-whitened** data. Sent high byte first | DN504 `culCalcCRC`; our `rc_link/cc2500.py` |
| CRC (CC2400_EN = 1) | polynomial x^16 + x^15 + x^2 + 1 (= 0x8005), over the data field only. Whitening is impossible in this mode. **Init value not stated [U]**: try 0xFFFF, then 0x0000 | CC2400 SWRS042A §23.2; CC2500 §15.4.2 |
| FEC | fixed-length mode only. Applies to data + CRC. Details in §5 | CC2500 §18; DN504 |

**MSK demodulation note.** Demodulate MSK as CPFSK with h = 0.5 (deviation = R/4): an FM discriminator plus symbol
timing.
- Because of the datasheet's inversion note, **try both bit polarities**.
- If neither polarity shows the sync word, try differential decoding (b_k ⊕ b_k−1), because of "data coding may
  differ" [U]. Settle it with the first real capture.

---

## 1. Graupner HoTT (`HOTT_cc2500.ino`, MPM protocol 57)

### 1.1 PHY
| Parameter | Register (value) | Result | Tag |
|---|---|---|---|
| Base frequency | FREQ2/1/0 = 0x5C 0x6C **0x6E** (`HOTT_FREQ0_VAL`, +`HOTT_COARSE` = 0; the comment "B9" in the file is stale) | **2403.01239 MHz** | [code][calc] |
| Channel spacing | MDMCFG1[1:0] = 3, MDMCFG0 = 0xAA | 338.0127 kHz | [calc] |
| Channel use | CHANNR = 3 · (n + 1), n = 0…74 | **f_n = 2404.0264 + 1.014038 · n MHz** (2404.026 to 2479.065 MHz, 75 channels) | [code][calc] |
| Modulation | MDMCFG2 = 0x73 → MOD_FORMAT 7 = **MSK**, DC filter off, no Manchester | MSK, h = 0.5 (±62.5 kHz) | [code][spec] |
| Data rate | MDMCFG4 = 0x2D, MDMCFG3 = 0x3B (E = 13, M = 59) | **249.939 kbaud** (coded); 124.97 kbit/s net | [calc] |
| MSK phase fraction | DEVIATN = 0x47 → M = 7 = whole symbol (true MSK) | | [spec] |
| RX filter (MPM) | MDMCFG4[7:4] = 2 | 541.7 kHz | [calc] |
| Sync mode | MDMCFG2[2:0] = 3 → 30/32 | on air: **D3 91 D3 91** | [code][spec] |
| Sync word, TX→RX | SYNC1 = 0xD3, SYNC0 = 0x91 (rewritten before every TX) | D391 | [code] |
| Sync word, RX→TX telemetry | **0x2C6E** when the TX's byte 29 has any of bits 7..3 set ("Sync" sub-protocol with telemetry); otherwise D391. Note 0x2C6E = bitwise NOT of 0xD391 | | [code] |
| Preamble | MDMCFG1 = 0xA3 → NUM_PREAMBLE = 2 | 4 bytes | [calc] |
| FEC | MDMCFG1[7] = 1 | **on**: K = 4, rate 1/2, 4×4 interleave (§5) | [code][spec] |
| Whitening | PKTCTRL0 = 0x44 → WHITE_DATA = 1 | on | [code] |
| CRC | PKTCTRL0[2] = 1, CC2400_EN = 0 | HW CRC-16 0x8005 / 0xFFFF | [code][spec] |
| Length mode | PKTCTRL0[1:0] = 0 (fixed); PKTLEN = 0x32 | **TX 50 B**; RX telemetry PKTLEN = **22 B** | [code] |
| Address check | PKTCTRL1 = 0x04 (APPEND_STATUS, no address check) | none | [code] |
| Air time | TX: 4 + 4 + 2·(50 + 2 + 2) = 116 B = 928 bit | **3.713 ms** | [calc] |
| | RX: 4 + 4 + 2·(22 + 2 + 2) = 60 B = 480 bit | 1.920 ms | [calc] |

### 1.2 Timing (10 ms frame, one hop per frame) [code]
| t (µs, from frame start) | Event |
|---|---|
| 0 | tune to the frame's channel, LBT receive (CCA), build the packet |
| 1095 | STX: the TX packet (3.71 ms on air). If CCA was busy, transmit at reduced power (LBT) |
| 4975 | switch to RX, same channel. Sync = 2C6E or D391 (see above), PKTLEN 22 |
| 9000 | read telemetry; next frame at 10000 |

So the TX burst repeats every **10.0 ms**, with the receiver's telemetry burst 1.9 ms long in the 4.0 ms window after
it, on the **same** channel.

### 1.3 Hopping [code]
- There are 16 fixed hop tables of 75 entries. Each is a permutation of 0…74 (checked here).
- MPM picks the table at random (`num_ch = random % 16`). It is **not** derived from the TX ID.
- Bytes 0–1 of every TX packet carry the table's 16-bit tag (`HOTT_hop_val`, little-endian), so the table is known
  from any decoded packet.
- Sequence: the frame k channel is `table[k mod 75]`, i.e. f = 2404.0264 + 1.014038 · table[i] MHz.
- **Byte 2 = the index i of the channel this packet is sent on.** MPM tunes to `table[i]`, then writes i into byte 2,
  then increments. The source comment says "next frequency", but the code order gives the current one; verify on a
  capture [U].
- Whether genuine Graupner TXs use only these 16 tables is unknown [U]. A table tag outside the list would mean a new
  table.

| Table | Tag (u16) | Bytes 0,1 | Channel sequence n (i = 0…74) |
|---|---|---|---|
| 0 | 0xC06B | 6B C0 | 48, 37, 16, 62, 9, 50, 42, 22, 68, 0, 55, 35, 21, 74, 1, 56, 31, 20, 70, 11, 45, 32, 24, 71, 8, 54, 38, 26, 61, 13, 53, 30, 15, 65, 7, 52, 34, 28, 60, 3, 47, 39, 18, 69, 2, 49, 44, 23, 72, 5, 51, 43, 19, 64, 12, 46, 33, 17, 67, 6, 58, 36, 29, 73, 14, 57, 41, 25, 63, 4, 59, 40, 27, 66, 10 |
| 1 | 0xC34A | 4A C3 | 50, 23, 5, 34, 67, 53, 22, 12, 39, 62, 51, 21, 10, 33, 63, 59, 16, 1, 43, 66, 49, 19, 8, 30, 71, 47, 24, 2, 35, 68, 45, 25, 14, 41, 74, 55, 18, 4, 32, 61, 54, 17, 11, 31, 72, 52, 28, 6, 38, 65, 46, 15, 9, 40, 60, 48, 26, 3, 37, 70, 58, 29, 0, 36, 64, 56, 20, 7, 42, 69, 57, 27, 13, 44, 73 |
| 2 | 0xDB24 | 24 DB | 73, 51, 39, 18, 9, 64, 56, 34, 16, 12, 66, 58, 36, 25, 11, 61, 47, 40, 15, 8, 71, 50, 43, 20, 6, 62, 54, 42, 19, 3, 63, 46, 44, 29, 14, 72, 49, 33, 22, 5, 69, 57, 30, 21, 10, 70, 45, 35, 26, 7, 65, 59, 31, 28, 1, 67, 48, 32, 24, 0, 60, 55, 41, 17, 2, 74, 52, 38, 27, 4, 68, 53, 37, 23, 13 |
| 3 | 0x8E09 | 09 8E | 52, 60, 40, 21, 14, 50, 72, 41, 23, 13, 59, 61, 39, 16, 6, 58, 66, 33, 17, 5, 55, 64, 43, 20, 12, 54, 74, 35, 29, 3, 46, 63, 37, 22, 10, 48, 65, 31, 27, 9, 49, 73, 38, 24, 11, 56, 70, 32, 15, 1, 51, 71, 44, 18, 8, 45, 67, 36, 25, 7, 57, 62, 34, 28, 2, 53, 69, 42, 19, 4, 47, 68, 30, 26, 0 |
| 4 | 0x272E | 2E 27 | 50, 16, 34, 6, 71, 51, 24, 40, 7, 68, 57, 27, 33, 14, 70, 55, 26, 30, 5, 74, 47, 28, 44, 11, 67, 49, 15, 32, 9, 61, 52, 22, 37, 13, 66, 59, 18, 42, 3, 62, 46, 29, 31, 12, 60, 48, 19, 38, 1, 72, 58, 17, 36, 4, 64, 53, 21, 39, 0, 63, 56, 20, 41, 2, 65, 45, 25, 35, 10, 69, 54, 23, 43, 8, 73 |
| 5 | 0x217F | 7F 21 | 55, 38, 12, 62, 23, 52, 44, 3, 66, 18, 54, 36, 10, 74, 16, 56, 42, 9, 70, 17, 58, 33, 5, 69, 20, 50, 40, 1, 63, 24, 53, 37, 13, 65, 15, 48, 34, 4, 61, 22, 57, 31, 6, 64, 26, 46, 35, 11, 72, 21, 47, 30, 7, 68, 29, 45, 32, 8, 60, 19, 49, 43, 2, 67, 27, 51, 39, 0, 71, 28, 59, 41, 14, 73, 25 |
| 6 | 0x155B | 5B 15 | 70, 32, 18, 10, 58, 69, 38, 22, 2, 54, 67, 36, 19, 12, 57, 62, 34, 20, 14, 52, 63, 41, 15, 3, 51, 73, 42, 28, 6, 48, 60, 43, 29, 5, 45, 64, 31, 17, 4, 56, 65, 35, 26, 13, 53, 61, 37, 23, 1, 49, 68, 40, 16, 9, 47, 71, 39, 25, 7, 50, 66, 33, 24, 8, 59, 72, 44, 27, 11, 46, 74, 30, 21, 0, 55 |
| 7 | 0xEDE8 | E8 ED | 6, 45, 71, 27, 44, 10, 46, 74, 22, 32, 0, 55, 69, 21, 33, 4, 50, 66, 18, 38, 7, 57, 62, 19, 36, 1, 48, 70, 20, 40, 8, 47, 68, 15, 43, 2, 58, 61, 26, 42, 3, 56, 72, 23, 34, 14, 54, 67, 16, 37, 5, 59, 64, 24, 30, 12, 52, 65, 25, 39, 13, 49, 73, 17, 31, 9, 53, 60, 28, 35, 11, 51, 63, 29, 41 |
| 8 | 0x1D31 | 31 1D | 31, 65, 50, 20, 13, 37, 66, 45, 23, 5, 32, 69, 54, 19, 7, 39, 74, 52, 27, 1, 42, 64, 53, 22, 4, 43, 70, 58, 16, 3, 40, 71, 57, 17, 0, 35, 63, 56, 18, 9, 44, 72, 51, 21, 6, 33, 67, 46, 25, 11, 30, 73, 55, 15, 8, 36, 62, 48, 24, 10, 41, 60, 49, 29, 14, 34, 61, 47, 26, 2, 38, 68, 59, 28, 12 |
| 9 | 0x0986 | 86 09 | 67, 22, 49, 36, 13, 64, 28, 57, 37, 6, 65, 29, 46, 39, 3, 70, 26, 45, 35, 1, 62, 24, 58, 34, 10, 68, 19, 53, 33, 4, 66, 21, 52, 31, 7, 74, 18, 47, 32, 5, 61, 16, 51, 38, 8, 72, 23, 55, 30, 12, 73, 17, 59, 44, 0, 60, 15, 50, 43, 14, 63, 27, 48, 42, 11, 71, 20, 54, 41, 9, 69, 25, 56, 40, 2 |
| 10 | 0x56F7 | F7 56 | 19, 38, 14, 66, 57, 18, 44, 7, 74, 48, 23, 30, 6, 71, 58, 26, 32, 5, 61, 46, 20, 34, 0, 68, 45, 24, 36, 1, 70, 50, 27, 33, 10, 63, 52, 16, 42, 9, 65, 51, 15, 41, 11, 64, 53, 22, 37, 3, 60, 56, 28, 35, 4, 67, 49, 17, 39, 13, 69, 54, 25, 43, 2, 73, 55, 21, 31, 8, 62, 47, 29, 40, 12, 72, 59 |
| 11 | 0x6454 | 54 64 | 4, 52, 64, 28, 44, 14, 46, 74, 16, 32, 11, 50, 68, 27, 36, 0, 47, 70, 26, 34, 13, 57, 61, 18, 38, 6, 56, 62, 19, 40, 5, 58, 67, 17, 31, 12, 54, 63, 22, 33, 3, 53, 72, 21, 41, 10, 48, 66, 15, 35, 7, 45, 60, 20, 37, 9, 51, 69, 25, 42, 2, 59, 71, 24, 39, 1, 55, 65, 23, 30, 8, 49, 73, 29, 43 |
| 12 | 0xC42D | 2D C4 | 44, 66, 19, 1, 56, 35, 62, 20, 4, 54, 39, 70, 24, 5, 55, 31, 74, 26, 12, 58, 32, 60, 17, 10, 45, 37, 63, 22, 3, 50, 33, 64, 16, 7, 51, 34, 61, 21, 8, 48, 38, 68, 29, 0, 46, 36, 72, 28, 14, 49, 42, 69, 25, 6, 57, 43, 65, 18, 2, 52, 30, 71, 23, 13, 47, 41, 67, 15, 9, 53, 40, 73, 27, 11, 59 |
| 13 | 0x01D2 | D2 01 | 12, 16, 36, 46, 69, 6, 20, 44, 58, 62, 11, 19, 34, 48, 71, 1, 18, 42, 50, 74, 3, 25, 31, 47, 65, 0, 24, 33, 45, 72, 2, 23, 35, 56, 64, 10, 22, 38, 49, 63, 7, 26, 37, 51, 70, 14, 21, 30, 53, 67, 5, 15, 40, 52, 66, 9, 17, 39, 55, 60, 13, 27, 41, 54, 73, 4, 28, 32, 57, 61, 8, 29, 43, 59, 68 |
| 14 | 0xC253 | 53 C2 | 63, 42, 18, 2, 57, 71, 34, 22, 10, 48, 67, 36, 25, 4, 46, 60, 31, 28, 6, 47, 74, 37, 15, 0, 55, 65, 32, 24, 12, 56, 66, 40, 27, 14, 52, 62, 38, 19, 3, 50, 73, 33, 29, 11, 53, 61, 35, 16, 7, 58, 72, 41, 26, 5, 59, 69, 30, 20, 9, 51, 68, 44, 23, 1, 49, 70, 39, 17, 8, 54, 64, 43, 21, 13, 45 |
| 15 | 0x1180 | 80 11 | 52, 1, 71, 17, 36, 47, 7, 64, 26, 32, 53, 5, 60, 20, 42, 57, 2, 66, 18, 34, 56, 4, 63, 24, 35, 46, 13, 72, 22, 30, 48, 0, 67, 21, 39, 50, 3, 74, 16, 31, 59, 14, 61, 23, 37, 45, 6, 65, 19, 44, 51, 11, 62, 27, 41, 55, 9, 68, 15, 38, 58, 8, 70, 29, 40, 54, 10, 69, 28, 33, 49, 12, 73, 25, 43 |

(The tables were extracted by script from `src_gpl_reference/mpm/HOTT_cc2500.ino` and checked to be permutations.)

### 1.4 TX packet: 50 bytes after FEC decoding and de-whitening, then CRC [code]
| Byte | Content |
|---|---|
| 0–1 | hop-table tag, little-endian (table above) |
| 2 | hop index i (0…74) of this packet's channel (see 1.3 [U]) |
| 3 | 0x00 = channels 1–8 + 9–12; 0x01 = channels 1–8 + 13–16 (alternates every packet). Failsafe packets: 0x40 = failsafe for CH1–12, 0x41 = CH13–24 |
| 4–27 | 12 channels × u16 **little-endian**. Value = pulse width in **0.5 µs**: v = ch·5/4 + 1720 for ch 0…2047, i.e. 1720…4278 = 860…2139 µs (−125…+125 %); 3000 = 1500 µs. Slots 0–7 = CH1–8; slots 8–11 = CH9–12 or CH13–16 according to byte 3. In failsafe packets: bit 15 = hold, bit 14 = position flag with the value in the low bits |
| 28 | 0x80 = no sensor. With telemetry: 0x89…0x8E = sensor being polled (sensor id 9…14). In text-menu mode: key codes 0xXF none, 0xXD down, 0xXB up, 0xX9 enter, 0xXE right, 0xX7 left |
| 29 | 0x02 normal; 0x01 text-config menu; "Sync" mode: ((seq + 1) << 3) ∣ 2 (binary sensor) or ∣ 1 (text). Bind: 0x02, then 0x1A / 22 / 2A / 0A / 12 as the RX replies (source comment) |
| 30–38 | 0xFF × 9 |
| 39 | 0x07 (constant, meaning unknown) |
| 40–44 | **TX ID (5 B)**. During bind: 0xFA × 5 |
| 45–49 | **RX ID (5 B)**, learned at bind. During bind: the TX ID |
| (50–51) | HW CRC-16 over bytes 0–49, big-endian, then whitened, then FEC with the trellis terminator |

Example TX IDs in the source: 7C 94 00 0D 50 and EA 4D 00 01 50. Example telemetry TX IDs: D4 20 F2 E6 F6.

### 1.5 RX telemetry packet: 22 bytes (same coding: FEC + whitening + CRC) [code]
| Byte | Content |
|---|---|
| 0–4 | TX ID (must match the TX packet's bytes 40–44) |
| 5–9 | RX ID |
| 10 | warnings |
| 11 | page number (sensors 0–4; config menu 0–0x12) |
| 12–21 | page data. Page 0 (RX): 13 = RX voltage ×10, 14 = temperature + 20 °C, 15 = RX RSSI (CC2500 format), 16 = LQI %, 17 = minimum voltage ×10, 18–19 = longest lost-packet time in ms (LE), 20 = events. Sensor page 1: byte 12 = sensor id << 4 (0x90…0xE0). Text pages: 13–21 = 9 ASCII characters (bit 7 = highlight) |

Bind reply examples (source comments): GR-16 `D4 20 F2 E6 F6 31 BD 01 00 90 00 FF 03 00 9E 1B 00 00 00 00 00 00`;
GR-12L `D4 20 F2 E6 F6 6E EE 01 00 B1 …`; Vector `D4 20 F2 E6 F6 00 00 3A 01 A1 …`. Bytes 0–4 are the TX ID and
5–9 the RX ID.

### 1.6 Decoder recipe (HoTT)
1. **Find the burst.** A 3.7 ms MSK burst about 0.5 MHz wide every 10 ms, on the 1.014 MHz grid from 2404.026 MHz.
   It is followed on the same channel by a 1.9 ms burst (telemetry). Demodulate as CPFSK with h = 0.5, 249.939 kBd.
2. **Sync.** Find 32 bits of 1010 preamble, then D391D391 (30/32 tolerance) in either polarity. For a telemetry burst,
   also try 2C6E2C6E.
3. **Collect the coded bits.** Take 2 · (50 + 2 + 2) · 8 = 864 coded bits (telemetry: 2 · 26 · 8 = 416).
4. **De-interleave and decode.** De-interleave per 32-bit block, Viterbi-decode (§5) → 54 bytes (telemetry 26).
5. **De-whiten and check.** De-whiten the first 52 bytes (PN9, seed all ones) → 50 data + 2 CRC. Check CRC-16 0x8005 /
   0xFFFF.
6. **Fix the polarity.** The wrong polarity makes the Viterbi path metric large: the all-inverted coded stream is not a
   codeword, because g0 has an even number of taps. Use the polarity that gives the small metric and a passing CRC.
7. **Report:**
   - TX ID (bytes 40–44) and RX ID (45–49);
   - hop table (bytes 0–1) and index (byte 2): the whole future hop sequence follows from them;
   - channel values, the bind flag (bytes 40–44 = FA × 5) and the failsafe flag (byte 3 = 0x4x);
   - for a telemetry packet: the TX / RX IDs, the page, RX voltage and RSSI.

---

## 2. Hitec Optima / Minima (`Hitec_cc2500.ino`, MPM protocol 39; sub-protocols OPT_FW, OPT_HUB, MINIMA)

### 2.1 PHY
| Parameter | Register (value) | Result | Tag |
|---|---|---|---|
| Base frequency | FREQ = 0x5C 0x85 0xE8 | **2405.59985 MHz** | [calc] |
| Channel spacing | CHANSPC_E = 3 (MDMCFG1 = 0x73), CHANSPC_M = 0x7A | 299.927 kHz | [calc] |
| Channel use | CHANNR **even**, 0…140 (data); bind CHANNR = 0, 10, …, 130 (14 channels) | f = 2405.59985 + 0.299927 · CHANNR MHz: data channels 2405.600…2447.590 MHz (71 possible, 0.6 MHz apart); bind 2405.600…2444.590 MHz | [code][calc] |
| Modulation | MDMCFG2 = 0x73 → **MSK**, sync mode 3 (30/32) | h = 0.5 | [code][spec] |
| Data rate | MDMCFG4 = 0x3D, MDMCFG3 = 0x3B | **249.939 kbaud** | [calc] |
| MSK phase fraction | DEVIATN = 0x01 → M = 1 (short phase transition; exact SmartRF meaning [U]) | still ±90° per symbol | [spec] |
| RX filter (MPM) | 0x3D → 464.3 kHz | | [calc] |
| Preamble | MDMCFG1 = 0x73 → NUM_PREAMBLE = 7 | **24 bytes** (a strong detection trigger) | [calc] |
| FEC | MDMCFG1[7] = 0 | off | [code] |
| Whitening | PKTCTRL0 = 0x45 → WHITE_DATA = 1 | on | [code] |
| CRC | CRC_EN = 1, CC2400_EN = 0 | HW CRC-16 0x8005 / 0xFFFF over the length byte + data | [code][spec] |
| Length mode | LENGTH_CONFIG = 1 (**variable**, byte 0 = length); PKTLEN = 0xFF (maximum) | | [code] |
| Sync word, bind | SYNC1:SYNC0 = D3:91 | on air D3 91 D3 91 | [code] |
| **Sync word, normal** | SYNC1 = TXID[3], SYNC0 = TXID[2] (the same bytes as packet bytes 3 and 2) | on air **id3 id2 id3 id2**, e.g. 6A 03 6A 03 for TX ID 00 03 6A | [code] |
| Air time, data packet | 24 + 4 + 27 + 2 = 57 B = 456 bit | **1.824 ms** | [calc] |
| Air time, bind packet | 24 + 4 + 23 + 2 = 53 B | 1.696 ms | [calc] |

### 2.2 Timing (22.5 ms frame; one channel per frame, 4 identical packets) [code]
| t (µs) | Event |
|---|---|
| 0 | fast-tune to channel `hop[k]`, build the packet |
| 462, 3198, 5934, 8670 | 4 TX packets (2.736 ms apart). They are identical except byte 23 (0x80, 0x40, 0x20, 0x10) |
| 11406 | RX on the same channel (window 4.636 ms; the receiver answers in it) |
| 16042 | read telemetry; next frame at 22500 |

### 2.3 Hopping
- **21 channels per TX**, used in table order, one per 22.5 ms frame → a period of **472.5 ms** [code].
- MPM generates the table from its 32-bit protocol ID, not from what a genuine Hitec TX does:
  - LCG r ← r · 0x0019660D + 0x3C6EF35F (32-bit), seeded with the ID;
  - candidate = ((r >> 8) mod 141) & 0xFE;
  - reject duplicates and keep at most 8 channels in each third (≤ 47, 48–93, ≥ 94) [code].
- **A genuine TX's algorithm is unknown [U]**; MPM's own dump table (TX ID 00 03 6A) is
  `00 3A 4A 32 0C 58 2A 10 26 20 08 60 68 70 78 80 88 56 5E 66 6E` (CHANNR values).
- The table is **sent in clear at bind** as CHANNR/2 values:
  - bind phase 0x72: entries 0–13 in bytes 5–18;
  - 0x73: entries 14–20 in bytes 5–11;
  - Minima: alternating 0x82 / 0x83 with the same content.
- **Passive approach:** the sync word already gives the ID. Learn the 21 channels by observing (a burst of 4 packets
  per channel every 22.5 ms), or read them from bind packets.

### 2.4 TX packet (after de-whitening; byte 0 = length, CRC follows the last byte) [code]
| Byte | Normal (length 0x1A → 27 B) | Bind (length 0x16 → 23 B) |
|---|---|---|
| 0 | 0x1A | 0x16 |
| 1–3 | **TX ID** id1 id2 id3 (sync = id3 id2) | TX ID |
| 4–21 | 9 channels × u16 **big-endian**. Linear: −100 % → 0x1B87 (7047), +100 % → 0x3905 (14597), centre 0x2A46 (10822); extended limits beyond. Pulse-width scale unknown [U] | 4: bind phase (Optima 0x72, 0x73, 0x74, …, 0x7B; Minima phase + 0x10, i.e. 0x82 / 0x83); 5–18: hop table / 2, or 0x55 at 7–8 (phase 0x74), or 0x2A 0x46 0x2A 0x47 0x2A at 14–18 (phase 0x7B) |
| 19 | (channel data) | sequence 0x08, 0x04, 0x02, 0x01 across the 4 repeats |
| 20–21 | (channel data) | F5 marker (as bytes 24–25 below, +1: 0xE0 / 0xDC) or 00 00 |
| 22 | **0xEE** (constant) | 0xEE |
| 23 | sequence 0x80, 0x40, 0x20, 0x10 (repeat 1…4 in the frame) | – |
| 24–25 | alternate frames: 00 00, or **F5 DF** (F5 DB on every 9th such frame; MPM's counter wraps at 59) | – |
| 26 | 0x00 | – |
| + 2 | HW CRC-16 (big-endian) over bytes 0…L | idem |

### 2.5 RX → TX telemetry (same sync = TX ID) [code]
- Data reply:
  - length byte 0x0C (13 B + CRC);
  - bytes 1–3 = TX ID; byte 4 = 0; byte 5 = frame number (0x00, 0x11…0x1C); byte 12 = 0;
  - CRC passes.
- Frame contents (from the MPM notes):
  - 0x00 / 0x11: RX battery = (H·256 + L) / 28 V;
  - 0x12 / 0x13: GPS latitude / longitude;
  - 0x14: speed, altitude, temperature 1;
  - 0x15: fuel, RPM;
  - 0x16: date and time;
  - 0x17: heading, satellites;
  - 0x18: voltage / current;
  - 0x19: servo currents;
  - 0x1A: airspeed;
  - 0x1B: variometer.
- Example: `0C 1C A1 2B 00 00 00 00 00 00 00 8D 00` → TX ID 1C A1 2B, RX battery 0x008D / 28 = 5.03 V.
- Bind reply: length 0x0A; byte 4 = 0x7X; bytes 5–9 = 05 06 07 08 09.

### 2.6 Decoder recipe (Hitec)
1. Trigger on ≥ 16 bytes of 1010 MSK preamble on the 0.6 MHz grid from 2405.6 MHz.
2. Read the next 32 bits; they must be the pattern A B A B (a 16-bit word repeated). Bind uses D391D391.
3. De-whiten from the length byte. Byte 0 ∈ {0x1A, 0x16, 0x0C, 0x0A}. Check the CRC over 1 + L bytes.
4. Consistency check: bytes 3, 2 equal the sync word.
5. Report the TX ID (3 B), the 9 channels, the repeat index (byte 23), bind (and the hop table from bind packets), and
   telemetry frames.
6. The 4-repeat burst group every 22.5 ms is an additional signature.

---

## 3. Futaba S-FHSS (`Futaba_cc2500.ino`, MPM protocol 21; mirrors deviation `sfhss_cc2500.c`)

| Parameter | Value | Tag |
|---|---|---|
| Base / grid | FREQ 0x5C 0x4E 0xC4 = 2399.99963 MHz; spacing (E 3, M 0x3B) 249.939 kHz; CHANNR = 16 + 6k, k = 0…29 → **f_k = 2403.9987 + 1.49963 · k MHz** (to 2447.488 MHz) | [code][calc] |
| Modulation | MDMCFG2 = 0x83 → **2-FSK** (no Gaussian), sync mode 3 | [code] |
| Rate / deviation | 0x7C / 0x43 → **128.143 kbaud**; DEVIATN 0x44 → **±38.086 kHz**; RX filter 232 kHz | [calc] |
| Preamble / sync | 4 B; **D391 sent twice** (30/32 per the register). The file comment says "15/16", and RFENG counts a 16-bit sync, so accept 16 or 32 bits [U] | [code][U] |
| FEC / whitening | off / **off** | [code] |
| CRC | PKTCTRL0 = 0x0C → CRC on with **CC2400_EN = 1** → poly 0x8005, over the 13 data bytes; init [U] (try 0xFFFF, then 0x0000) | [code][spec][U] |
| Length | fixed **13 B** | [code] |
| Air time | 4 + 4 + 13 + 2 = 23 B = 184 bit → 1.436 ms (1.311 ms with a 16-bit sync) | [calc] |
| Timing | 2 packets per channel, **1.625 ms** apart (DATA1, DATA2), then hop; hop period **6.8 ms** | [code] |
| Hop | k ← k + code + 2; if k > 29: { if k < 31: k ← k + code + 2 }; k ← k − 31. code ∈ 0…27, fixed per TX (MPM: random) | [code] |

**Packet (13 B)** [code]:
| Byte | Bits |
|---|---|
| 0 | 0x81 (0x80 seen on Orange; TM-FH 0x82) |
| 1–2 | **16-bit TX ID** (MPM limits it to runs ≤ 6 identical bits) |
| 3–4 | 0x00 0x00 in MPM. Seen on genuine TXs: T8J 42 07, T10J 0F 09, TM-FH 9A 06, so these bytes **hint at the TX model** |
| 5 | k << 3 ∣ ch0[11:9] (k = current channel, 0…29) |
| 6 | ch0[8:1] |
| 7 | ch0[0] << 7 ∣ ch1[11:5] |
| 8 | ch1[4:0] << 3 ∣ ch2[11:9] |
| 9 | ch2[8:1] |
| 10 | ch2[0] << 7 ∣ ch3[11:5] |
| 11 | ch3[4:0] << 3 ∣ code[4:2] |
| 12 | code[1:0] << 6 ∣ command |

- **Command** bits:
  - bit 0 = DATA1 / DATA2;
  - bit 1 = 1 normally (meaning unknown);
  - bit 2 = failsafe data;
  - bit 3 = CH5–8 (else CH1–4).
- **Channel values:** 12-bit, 1020…2020, centre 1520, "values grow down" (−100 % = 2020).
- **Channel order** per model: T8J sends CH1–4 / CH5–8 on DATA1 / DATA2 with the order alternating between frames; XK,
  T10J and TM-FH send CH1–4 on DATA1 and CH5–8 on DATA2.
- **Verification:** the CRC (init to be settled), plus the hop index in byte 5 matching the measured channel, plus the
  hop code from bytes 11–12, which predicts the next channel.
- **No bind:** the RX latches onto the first TX it hears.

---

## 4. FrSky R9 ACCST (`FrSkyR9_sx1276.ino`, MPM protocol 65)

### 4.1 PHY [code][spec]
| Parameter | Value |
|---|---|
| Modulation | LoRa **SF6**, **BW 500 kHz**, **CR 4/5**, LDRO off (ModemConfig3 bit 3 = 0), IQ normal |
| Header | **implicit** (ModemConfig1 bit 0 = 1; SF6 requires it). The decoder must know the length: 26 (FCC / flex) or 14 (EU) for the uplink, 13 for telemetry |
| Payload CRC (chip) | **off** (ModemConfig2 bit 2 = 0). The CRC is in software, inside the payload |
| Preamble | **9 symbols** (+ 2 sync + 2.25 SFD) |
| Sync word | not written → SX1276 reset value **0x12** (RegSyncWord 0x39) |
| Symbol time | 2^6 / 500 kHz = **128 µs** |
| Air time [calc] | 26 B: 8 + ⌈(8·26 − 24 + 28 − 20)/24⌉ · 5 = 48 payload symbols → (13.25 + 48) · 128 µs = **7.84 ms**. 14 B / 13 B: 28 symbols → **5.28 ms** |
| Period | **20 ms**: 0.46 ms settle, TX, (telemetry builds: TX slot 11.14 ms, RX 7.4 ms, 1 ms); without telemetry 20 ms |
| Power | PA_BOOST, 17 − (15 − option) dBm (≤ 17 dBm) |

The existing `rc_link/lora_phy.decode_implicit(sf=6, cr=1, length=26)` path applies directly.

### 4.2 Frequency grids (FRF = num · 256; f = num · 15.625 kHz exactly) [code][calc]
| Sub-protocol | num | Grid | Channels | Hop |
|---|---|---|---|---|
| R9_FCC (4; 8CH = 6) | 0xE200 + 32·n | **904.0 + 0.5·n MHz** | n = 0…42 (43; 904.0–925.0) | n ← (n + chanskip) mod 43 |
| R9_915 flex (0; 8CH = 2) | 0xE4C0 + 32·m | **915.0 + 0.5·m MHz** | n = 0…28 (29 slots); m = n for n < 27, m = chanskip (n = 27), chanskip + 1 (n = 28); 915.0–928.0 | n ← (n + chanskip) mod 29 |
| R9_868 flex (1; 8CH = 3) | 0xD700 + 32·m | **860.0 + 0.5·m MHz** (860.0–873.0); **binding uses the 915 flex grid** | as 915 flex | idem |
| R9_EU (5; 8CH = 7) | 0xD7D0 + 16·n | **863.25 + 0.25·n MHz** | n = 0…18 (19; 863.25–867.75) | n ← (n + chanskip) mod 19 |

- **chanskip** = 1…24, chosen at random per session by MPM. It is sent in every packet.
- With EU, chanskip = 19 would never hop (MPM quirk [U]).
- The index is updated **before** each transmission, and byte 3 (FCC / flex) holds the index of the channel being used.

### 4.3 Uplink packet, FCC / flex: 26 B [code]
| Byte | Content |
|---|---|
| 0–2 | **TX ID** (rx_tx_addr[1..3]) |
| 3 | hop index n of this packet (check: the frequency must equal the grid value for n) |
| 4 | chanskip (1…24) |
| 5 | RX number (model match, 0…63) |
| 6 | flags: 0x00 normal; 0x10, 0x12, …, 0x1E failsafe (low nibble = channel / 2); 0x20 range check. **Bind**: 0x01 ∣ 0x20 (868) ∣ 0x40 (telemetry off) ∣ 0x80 (CH9–16) |
| 7 | 0x00 |
| 8–19 | 8 channels × 12 bits in 4 triplets (a, b): byte0 = a[7:0]; byte1 = a[11:8] ∣ b[3:0] << 4; byte2 = b[11:4]. PXX value 64…1984 (1024 = 1500 µs; µs = 1500 + (v − 1024)·2/3); bit 11 set = channels 9–16. In 16-channel mode, packets alternate CH1–8 / CH9–16 |
| 20 | (RX_seq << 4) ∣ TX_seq (S.Port flow control) |
| 21 | S.Port byte count (0–2) |
| 22–23 | S.Port data |
| 24–25 | **CRC-16, little-endian**, over bytes 0–23 |

**CRC-16 (FrSkyX_crc)** [code]:
- c = 0; per byte b: c = ((c << 8) & 0xFFFF) ^ T[(c >> 8) ^ b].
- T is the 256-entry table of **reflected** CRC-16 (poly 0x8408, the KERMIT table, T[1] = 0x1189) used in this
  non-reflected shift.
- MPM builds T from a 16-entry low-nibble table: T[v] = Tlo[v & 15] ^ 0x1081 · (v >> 4). This integer product equals
  the carry-less one because the bits do not overlap. Tlo = 0000 1189 2312 329B 4624 57AD 6536 74BF 8C48 9DC1 AF5A
  BED3 CA6C DBE5 E97E F8F7.
- It is **the same function** as our FrSky D16 CRC (`rc_link/cc2500.py: crc16_d16`).

### 4.4 Uplink packet, EU: 14 B [code] (MPM's channel packing is marked TODO → layout partly [U])
| Byte | Content |
|---|---|
| 0–2 | TX ID |
| 3 | chanskip |
| 4 | RX number |
| 5 | flags (bind: 0x01 ∣ 0x10 (16 ch) ∣ 0x40 (telemetry off) ∣ 0x80 (CH9–16)) |
| 6–11 | 4 channels (packing not implemented in MPM [U]) |
| 12 | (RX_seq << 4) ∣ 0x08 |
| 13 | **CRC-8**: poly 0x07, init 0xFF, MSB first, no final XOR, over bytes 0–12 |

### 4.5 Telemetry from the receiver: 13 B, same LoRa settings, inside the 7.4 ms window after the uplink [code]
| Byte | Content |
|---|---|
| 0 | bit 7 = 1: RX RSSI (value << 1); otherwise RX battery voltage ((value << 1) + 1) |
| 1–8 | the FrSky X (D16) telemetry body: sequence and S.Port bytes, as D16 telemetry bytes 5–12 [U detail] |
| 9–10 | TX ID bytes 0, 1 (= rx_tx_addr[1], [2]) |
| 11–12 | CRC-16 LE, the same function seeded with **id0 + 256 · id1** (rx_tx_addr[1] + rx_tx_addr[2]·256), over bytes 0–10 |

### 4.6 Decoder recipe (R9)
1. Look for LoRa SF6 / BW500 chirps (128 µs symbols) with a 9-symbol preamble every 20 ms on 904–925 (0.5 MHz grid),
   915–928, 860–873 or 863.25–867.75 (0.25 MHz grid).
2. Check the sync word 0x12. Decode the implicit header with CR 4/5 for lengths 26, 14 and 13.
3. Accept a packet when its CRC passes (CRC-16 LE for 26 B; CRC-8 for 14 B; the seeded CRC-16 for 13 B telemetry).
4. For FCC / flex: check byte 3 against the measured frequency, then predict the next channel from byte 4.
5. Report the TX ID, RX number, channels, bind / failsafe / range-check flags, the plan (FCC / flex915 / flex868 / EU,
   from the frequency) and 8 vs 16 channels (bit 11 alternating).
6. **Stock-firmware caveat:** ACCESS R9 firmware is different and unpublished [U]. ExpressLRS on R9M hardware is a
   separate protocol (file 07).

---

## 5. CC2500 / CC1101 FEC: convolutional code and interleaver (for HoTT) [spec][ver]

Sources:
- DN504 (SWRA113A) §4 (encoder figure, trellis) and §6 (reference C code and test vector);
- CC2500 datasheet §18.

All of the following was re-implemented clean-room in Python and **reproduces the DN504 test vector bit-exactly**
[ver]. The script is in this session's scratchpad (`chk/fec.py`); it is not part of the repo.

**Encoder:** non-recursive, rate 1/2, constraint length L = 4 (memory 3), d_free = 6. Let x[n] be the input bit (data
taken MSB first). Then:
- **first output bit (g1) = x[n] ⊕ x[n−2] ⊕ x[n−3]**, i.e. 1 + D² + D³;
- **second output bit (g0) = x[n] ⊕ x[n−1] ⊕ x[n−2] ⊕ x[n−3]**, i.e. 1 + D + D² + D³.

Equivalently, in DN504's table form: symbol = T[8·x[n−3] + 4·x[n−2] + 2·x[n−1] + x[n]] with
T = {0,3,1,2, 3,0,2,1, 3,0,2,1, 0,3,1,2}, where bit 1 of the symbol (g1) is sent first.

**Encoder state:** starts at 0 at the start of each packet's data field. The preamble and sync are not coded.

**Input to the encoder:**
- the whitened data + CRC (whitening happens before FEC);
- then the **trellis terminator**: one byte **0x0B** if the number of data + CRC bytes is odd, two bytes **0x0B 0x0B**
  if it is even. The total is therefore even (a multiple of 2 bytes = 32 coded bits = one interleaver block).
- The leading 000 of 0x0B terminates the trellis; the rest fills the block and keeps transitions on air.
- Whether the terminator is whitened is not stated. It does not matter: it is discarded after decoding.

**Interleaver:** a 4×4 matrix of 2-bit coded symbols, one block = 16 symbols = 4 coded bytes.
- Write the 16 symbols s[0..15] (in coder output order) row-wise: M[r][c] = s[4r + c].
- Transmitted symbol j (j = 0..15, each sent MSB first): **out[j] = M[3 − (j mod 4)][3 − ⌊j/4⌋]**, i.e. the columns
  are read from right to left, each column bottom to top.
- De-interleave in the receiver with s[4·(3 − (j mod 4)) + (3 − ⌊j/4⌋)] = rx[j].

**Test vector (DN504 §6) [ver]:**
| Stage | Bytes |
|---|---|
| Packet (length + data) | 03 01 02 03 |
| + CRC-16 0x8005 / 0xFFFF | 03 01 02 03 30 3A |
| + terminator | 03 01 02 03 30 3A 0B 0B |
| Encoder output | 00 0E 8C 03 7C 0D F0 0E 82 8C 0E 5E F0 D1 8C D1 |
| Interleaver output (on air) | C8 3C 00 20 84 CF 33 31 A2 FC 40 4A 44 30 47 EF |

**Viterbi:**
- 8 states; start in state 0. Use the Hamming branch metric, or soft metrics from the discriminator; the chip itself
  uses 8-bit soft values.
- Trace back from the best (ideally zero) state after all coded bits, or after 3 terminator bits.
- The clean-room check decoded the vector exactly, and also with 3 coded-bit errors spread over different blocks.
- The coded length on air is 2 · (N_data + N_crc + N_term) bytes. HoTT TX: 2 · (50 + 2 + 2) = 108 B; HoTT
  telemetry: 2 · (22 + 2 + 2) = 52 B.

**Effective rate:** HoTT's 249.94 kBd on air carries 124.97 kbit/s of data.

---

## Summary: implementable now vs uncertain

| Protocol | Implementable now (passive, verifying) | Uncertain [U] |
|---|---|---|
| **HoTT** | Full chain: MSK h = 0.5 at 249.94 kBd on the 75-channel 1.014 MHz grid → sync D391 ×2 → 4×4 de-interleave → K = 4 Viterbi (DN504 vector reproduced) → PN9 de-whitening → CRC-16 0x8005. TX / RX IDs (5 B each), the 16 hop tables (identified by bytes 0–1) and the hop index, 12 channel values, telemetry pages | MSK bit polarity and any differential coding (try both); whether byte 2 is the current or the next index; whether genuine TXs use only these 16 tables; the telemetry sync 2C6E applies only in "Sync" mode |
| **Hitec** | MSK 249.94 kBd, 24 B preamble, sync = TX ID (id3 id2 ×2), variable length, whitening, CRC-16. TX ID, 9 channels, the 4-repeat structure, bind tables, telemetry frames | The genuine TX's hop-table generation (learn it from observation or bind); the µs scale of the channel values; the DEVIATN = 1 MSK transition shape; MSK polarity |
| **S-FHSS** | 2-FSK 128.14 kBd, ±38.1 kHz, sync D391, 13 B fixed, no whitening. 16-bit TX ID, hop index + hop code (next channel predictable), 8 channels; bytes 3–4 hint at the model | The **CRC init** in CC2400 mode (poly 0x8005 known; try 0xFFFF / 0x0000); 16- vs 32-bit sync |
| **R9 ACCST** | LoRa SF6 / BW500 / CR 4/5, implicit header, preamble 9, sync 0x12, 26 B with CRC-16 LE (= our D16 CRC), TX ID, hop index + chanskip (next channel predictable), exact grids, RX telemetry with seeded CRC | The EU 14 B channel packing (TODO in MPM); the telemetry body bytes 1–8 in detail; stock ACCESS R9 (unpublished) |

## Corrections to earlier files
- **R9 frequencies** (`12_rc_links_more.md` §2 used 61 Hz per FRF LSB, MPM's debug approximation). The real step is
  32 MHz / 2^19 = 61.035 Hz, so num·256 LSBs = num · 15.625 kHz exactly:
  - FCC 904.0–925.0 MHz (not 903.479);
  - flex 915.0–928.0 and 860.0–873.0 MHz;
  - EU 863.25–867.75 MHz;
  - spacing exactly 500 / 250 kHz.
- **HoTT FREQ0** is 0x6E (the code value), not 0xB9 (a stale comment, used in `src_gpl_reference/cc2500_register_decoder.py`).
  Base 2403.0124 MHz, which `12_rc_links_more.md` already has right.
- **S-FHSS CRC:** CC2400_EN = 1 (PKTCTRL0 = 0x0C). The poly is still 0x8005, but the init is not established, and
  whitening is impossible in this mode. File 12 lists "HW CRC-16 0x8005 / FFFF"; the init there is [U].

## Sources
- MPM (GPL-3, oracle only), https://github.com/pascallanger/DIY-Multiprotocol-TX-Module/tree/master/Multiprotocol :
  - `HOTT_cc2500.ino`
  - `Hitec_cc2500.ino`
  - `Futaba_cc2500.ino` (S-FHSS)
  - `FrSkyR9_sx1276.ino`
  - `FrSkyDVX_common.ino` (FrSkyX_crc, FrSkyX_channels, S.Port sequence)
  - `Convert.ino` (channel scaling: `convert_channel_16b_nolimit`, `FrSkyX_scaleForPXX`; CHANNEL_MIN_100 = 204,
    CHANNEL_MAX_100 = 1844)
  - `SX1276_SPI.ino` and `iface_sx1276.h` (ModemConfig bit packing, BW / CR enums)
  - `Protocols_Details.md` (sub-protocols, tested receivers)
- TI CC2500 datasheet SWRS040C: https://www.ti.com/lit/ds/symlink/cc2500.pdf (§15 packet handling / whitening / CRC,
  §16 modulation and MSK note, §18 FEC / interleaving, register descriptions)
- TI DN504 "FEC Implementation" SWRA113A: https://www.ti.com/lit/an/swra113a/swra113a.pdf (encoder, trellis, terminator,
  reference code, test vector)
- TI CC2400 datasheet SWRS042A: https://www.ti.com/lit/ds/symlink/cc2400.pdf (§23.2 CRC polynomial x^16 + x^15 + x^2 + 1)
- Semtech SX1276/77/78/79 datasheet: https://www.semtech.com/products/wireless-rf/lora-connect/sx1276 (time-on-air
  formula, RegSyncWord 0x39 reset 0x12, SF6 implicit-header requirement)
- RFENG, "Futaba S-FHSS Protocol Overview": https://rfengfpv.wordpress.com/category/futaba/ (16-bit sync claim)
