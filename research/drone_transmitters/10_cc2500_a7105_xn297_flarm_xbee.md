# 10. CC2500 / A7105 / XN297 / FLARM / 802.15.4 (XBee) links: on-air formats for passive detection

Researched 2026-10-22 from the web only. Every fact has a source. `[U]` marks uncertain items. `[calc]` marks values
computed here from register values using the vendor datasheet formulas. Each constant table carries a licence tag:
- `[DS]`: vendor datasheet or standard; usable in the C++ library.
- `[GPL-only]`: found only in GPL-3 Multiprotocol / Deviation / SoftRF code. The behaviour can be described and
  re-implemented from this description, but the code must not be copied. Treat these values as reverse-engineered
  facts and verify them on real captures before trusting them.

Main sources (cited below by tag):
- **[CC2500]** TI CC2500 datasheet SWRS040C: https://www.ti.com/lit/ds/symlink/cc2500.pdf. Pages: §15.1–15.4
  p.29–32, register details p.62–67.
- **[DN502]** TI DN502, CRC implementation for CC1100/CC2500: https://ti.com/lit/swra111
- **[MPM-D]** Multiprotocol (GPL-3) `FrSkyD_cc2500.ino`:
  https://github.com/pascallanger/DIY-Multiprotocol-TX-Module/blob/master/Multiprotocol/FrSkyD_cc2500.ino
- **[MPM-X]** Multiprotocol `FrSkyX_cc2500.ino`:
  https://github.com/pascallanger/DIY-Multiprotocol-TX-Module/blob/master/Multiprotocol/FrSkyX_cc2500.ino
- **[MPM-C]** Multiprotocol `FrSkyDVX_common.ino` (register tables, CRC, hop tables):
  https://github.com/pascallanger/DIY-Multiprotocol-TX-Module/blob/master/Multiprotocol/FrSkyDVX_common.ino
- **[MPM-RX]** Multiprotocol `FrSky_Rx_cc2500.ino`:
  https://github.com/pascallanger/DIY-Multiprotocol-TX-Module/blob/master/Multiprotocol/FrSky_Rx_cc2500.ino
- **[A7105]** AMICCOM A7105 datasheet v1.1 (Feb 2010): https://www.g7smy.co.uk/img/rf600tx/a7105_datasheet.pdf.
  Pages: §9.2.14–9.2.22 p.19–21, §9.2.32–34 p.24–25, ch.14 p.56–57, ch.16 p.65–67.
- **[MPM-A2]** Multiprotocol `AFHDS2A_a7105.ino`:
  https://github.com/pascallanger/DIY-Multiprotocol-TX-Module/blob/master/Multiprotocol/AFHDS2A_a7105.ino
- **[MPM-A7]** Multiprotocol `A7105_SPI.ino` (A7105 register tables and ID):
  https://github.com/pascallanger/DIY-Multiprotocol-TX-Module/blob/master/Multiprotocol/A7105_SPI.ino
- **[MPM-XN]** Multiprotocol `XN297_EMU.ino`:
  https://github.com/pascallanger/DIY-Multiprotocol-TX-Module/blob/master/Multiprotocol/XN297_EMU.ino
- **[DEV-NRF]** Deviation (GPL-3) `nrf24l01.c`:
  https://github.com/DeviationTX/deviation/blob/master/src/protocol/spi/nrf24l01.c
- **[MPM-BY]** Multiprotocol `Bayang_nrf24l01.ino`:
  https://github.com/pascallanger/DIY-Multiprotocol-TX-Module/blob/master/Multiprotocol/Bayang_nrf24l01.ino
- **[DEVF]** DeviationTX forum threads on XN297 sniffing:
  https://deviationtx.com/forum/protocol-development/8969-capturing-data-from-cheap-drone-with-the-hackrf and
  https://www.deviationtx.com/forum/protocol-development/5655-cx-10-new-red-rf-sniffing
- **[OXF]** B. Wang et al., "On the Security of the FLARM Collision Warning System", ASIACCS 2022:
  https://www.cydcampus.admin.ch/dam/de/sd-web/qkQH--RLJBPN/On%20the%20Security%20of%20the%20FLARM%20Collision%20Warning%20System.pdf
  (p.3 Tables 2–3, §2.2–2.3)
- **[SRF-L]** SoftRF (GPL-3) `Legacy.h` / `Legacy.cpp`:
  https://github.com/lyusupov/SoftRF/blob/master/software/firmware/source/SoftRF/src/protocol/radio/Legacy.h
- **[SRF-FP]** SoftRF `freqplan.h`:
  https://github.com/lyusupov/SoftRF/blob/master/software/firmware/source/libraries/OGN/freqplan.h
- **[OGN-FAQ]** https://wiki.glidernet.org/ogn-receiver-faq
- **[EASA]** EASA ADS-L SRD-860 specification (ED Decision 2022/024/R): https://www.easa.europa.eu/en/downloads/137503/en
- **[CC2420]** TI CC2420 datasheet SWRS041c: https://www.ti.com/lit/ds/symlink/cc2420.pdf. Pages: §10 p.24–25,
  §16 p.36–38.
- **[XB9]** Digi XBee-PRO 900HP datasheet:
  https://lv.mouser.com/pdfdocs/Digi_International_XbeePro_900HP_Datasheet.pdf; user guide 90002173:
  https://digi.com/resources/documentation/digidocs/PDFs/90002173.pdf
- **[XB8]** Digi XBee-PRO 868 manual 90001020: https://www.sos.sk/a_info/resource/c/90001020_B.pdf

---

## 0. CC2500 PHY primitives (shared by FrSky D8/D16) [DS]

Source: [CC2500]. With f_xosc = 26 MHz:
- **Data rate:** R = (256+DRATE_M)·2^DRATE_E / 2^28 · f_xosc. DRATE_E = MDMCFG4[3:0], DRATE_M = MDMCFG3 (p.66).
- **Deviation:** f_dev = f_xosc/2^17 · (8+DEV_M)·2^DEV_E, from DEVIATN[6:4] = E and [2:0] = M (p.33).
- **Carrier frequency:** f_c = f_xosc/2^16 · FREQ[23:0] (p.65). The channel is f_c + CHANNR·Δch, with
  Δch = f_xosc/2^18·(256+CHANSPC_M)·2^CHANSPC_E (CHANSPC_M = MDMCFG0, CHANSPC_E = MDMCFG1[1:0]).
- **Channel filter bandwidth:** BW = f_xosc / (8·(4+CHANBW_M)·2^CHANBW_E) (p.66).
- **MDMCFG2** (p.67): bits [6:4] = MOD_FORMAT, where 000 is 2-FSK, 001 is GFSK and 111 is MSK. Bit 3 =
  MANCHESTER_EN. Bits [2:0] = SYNC_MODE: 001 means 15/16 bits of a 16-bit sync are needed; 011 means 30/32, a
  repeated 32-bit sync, so SYNC1 SYNC0 SYNC1 SYNC0 is sent.
- **GFSK shaping:** Gaussian, BT = 1 (p.33). A '0' is sent at −dev and a '1' at +dev (Table 23).
- **Sync word reset value:** SYNC1 = 0xD3, SYNC0 = 0x91 (p.62). Neither FrSky table writes SYNC1/0 ([MPM-C],
  [MPM-RX]), so FrSky uses **0xD391**, sent MSB first [U: MSB-first is the CC1101-family convention; the bit order
  is not stated on these pages].
- **Preamble:** alternating 1010…; the byte count comes from MDMCFG1.NUM_PREAMBLE (p.30). For MDMCFG1 = 0x23,
  NUM_PREAMBLE = 010, which is 4 bytes [U: the table on datasheet p.68 was not re-read].
- **Packet:** preamble, then sync, then the optional length byte, the optional address byte, the data, and the
  optional CRC-16 (Fig. 11, p.29).
  - Variable-length mode (PKTCTRL0.LENGTH_CONFIG = 01): the first byte after sync is the length, and it excludes
    itself and the CRC (p.30).
  - Address check (PKTCTRL1.ADR_CHK ≠ 0): the byte after the length is compared with ADDR (p.31, p.63).
- **Hardware CRC (PKTCTRL0.CRC_EN = 1, CC2400_EN = 0):** CRC-16 x^16+x^15+x^2+1 (0x8005) with the register
  initialised to 0xFFFF. It covers every byte after the sync word: the length, address and data ([DN502], [CC2500]
  p.32).
- **Whitening (PKTCTRL0.WHITE_DATA, bit 6):** PN9 with the register reset to all 1s. It covers everything after the
  sync (p.29). FrSky does not use it, because its PKTCTRL0 values (0x05 or 0x01) clear bit 6.

## 1. FrSky D8 ("D" / ACCST D8, CC2500)

### PHY
Register values are [GPL-only] (from [MPM-C] and [MPM-RX]); they are decoded with the [DS] formulas above.

| Reg | Value | Meaning [calc] |
|---|---|---|
| FREQ2/1/0 | 5C 76 27 | base 2404.00 MHz |
| MDMCFG4/3 | AA 39 | BW 135 kHz; **R = 31.04 kbit/s** |
| MDMCFG2 | 11 | **GFSK** (BT 1), no Manchester, 15/16 of the 16-bit sync |
| MDMCFG1/0 | 23 7A | 4-byte preamble [U]; Δch = **299.93 kHz** |
| DEVIATN | 42 | **±31.7 kHz** |
| PKTCTRL1 / PKTCTRL0 | 04 (TX) or 05 (RX/addr) / 05 | variable length, **HW CRC on**, no whitening |
| ADDR | rx_tx_addr[3] (0x03 when binding) | address byte = TX-ID low byte |
| PKTLEN | 0x19 | max length |

The task brief guessed 250 kb/s; that is wrong. The FrSky D8 register values give about 31 kbit/s GFSK.

### Frame
On air: `pre(4) | D3 91 | len=0x11 | 17 bytes | CRC16`. Total 18 + 2 bytes of payload and CRC, plus 6 bytes of
preamble and sync, which is about 208 bits, or **about 6.7 ms on air** [calc].

**Data packet** ([MPM-D]) [GPL-only]:

| Byte | Content |
|---|---|
| [0] | 0x11 |
| [1] | ID_lo (also the CC2500 address byte) |
| [2] | ID_hi |
| [3] | frame counter 0–187 |
| [4] | telemetry counter (or 0) |
| [5] | rx_tx_addr[1] (constant; 0x01 for Multi's non-clone ID) |
| [6..9] | low 8 bits of channels 1–4 |
| [10..11] | high nibbles of channels 1–4 (ch1 low nibble, ch2 high nibble, …) |
| [12..15], [16..17] | the same for channels 5–8 (12-bit values) |

**Bind packet** [GPL-only]: `11 03 01 ID_lo ID_hi idx h0..h4 00×6 rx_tx_addr[1]`. Bind uses ADDR 0x03 and channel 0.
Ten bind packets carry the 47-entry hop table, five entries each ([MPM-D]).

### CRC
Only the CC2500 hardware CRC is used: 0x8005, init 0xFFFF, over the length and the payload [DS]. The receiver also
checks CRC_OK and the 3 ID bytes ([MPM-RX]).

### Timing and hopping
- **Timing:** frames go out about every 9 ms. The TX cycle is 36 ms: three TX frames plus a 1.3 ms telemetry listen
  window ([MPM-D]).
- **Hop channel:** `hop[counter % 47]`; counter = (counter+1) % 188 ([MPM-D]) [GPL-only].
- **Hop table** (from the TX ID; [MPM-C] `Frsky_init_hop`) [GPL-only]:
  1. ch = ID0 & 7. spacing = ID1, then adjusted: if < 2, add 2; if > 0xE9, subtract 0xE7; if it is a multiple of
     47, add 1.
  2. ch_i = (ch + i·spacing) mod 235.
  3. The values 0x00, 0x5A and 0xDC are stored +1.
  4. 47 entries.
- **Multi's own D TX file** uses a fixed (i·0x1E) % 0xEB table, with i = 3, 23 and 47 bumped by 1 ([MPM-D]).
- The original TX sends its table in the bind packets, so **a passive receiver cannot rely on deriving the table.**
  Use a CHANNR grid of 0–234 at 299.93 kHz (2404.0–2474.5 MHz) and follow the sequence empirically.

### Passive identity
- **TX ID:** 16 bits in clear (bytes 1–2), plus byte 5.
- **Frame counter:** byte 3.
- **Channel values:** 8 × 12 bit.
- **Hop spacing:** the ID-specific step can be measured from consecutive packets.

## 2. FrSky D16 (ACCST "X", v1 and v2), FCC and EU-LBT

### PHY
Register values are [GPL-only] ([MPM-C], [MPM-RX]); decoded [calc].

| Mode | FREQ | MDMCFG4/3 | Rate | MDMCFG2 | DEVIATN | PKTCTRL0 |
|---|---|---|---|---|---|---|
| FCC v1 | 5C7627 = 2404.0 MHz | 7B/61 | **70.0 kb/s**, BW 232 kHz | 0x13: **GFSK**, **32-bit sync D391D391** (30/32) | 0x51 = ±57.1 kHz | 0x01: **no HW CRC** |
| LBT (EU) v1 | 5C8000 = 2405.0 MHz | 7B/F8 | **99.98 kb/s** | 0x03: **2-FSK**, 32-bit sync | 0x53 = ±69.8 kHz | 0x01 |
| FCC v2 | as FCC | MDMCFG3 = 0x84 | **76.97 kb/s** (Multi's comment: "70K to 77K") | 0x13 | 0x51 | **0x05: HW CRC on** |
| LBT v2 | as LBT | as LBT | 99.98 kb/s | 0x03 | 0x53 | **0x05** |

All modes use MDMCFG1/0 = 23/7A, so Δch = 299.93 kHz.

### Frame
- **Length byte:** 0x1D (30 bytes with the length byte) for FCC and v2; **0x20 (33 bytes) for LBT v1** ([MPM-X],
  [MPM-RX]).
- **Data packet** [GPL-only] ([MPM-X], [MPM-C]):

| Byte | Content |
|---|---|
| [1] | ID_lo |
| [2] | ID_hi |
| [3] | rx_tx_addr[1] (0x02 for Multi's non-clone ID) |
| [4] | (chanskip & 3) << 6 \| hop_index (0–46) |
| [5] | chanskip >> 2 |
| [6] | RX number |
| [7] | failsafe flag |
| [8] | (channel-bank flags) [U] |
| [9..20] | 8 channels × 12 bit, packed in pairs: b0 = c1, b1 = (c1>>8 & 0x0F) \| (c2<<4), b2 = c2>>4 |
| [21] | sequence byte (RX seq in the high nibble; TX seq) |
| [22] | S.Port byte count |
| [23..] | S.Port uplink |
| last 2 | software CRC, big-endian |

- Channels 1–8 and 9–16 alternate between frames (chan_start toggles with 0x08) ([MPM-C]).
- **Software CRC** [GPL-only] ([MPM-C], [MPM-RX]):
  - It is computed over packet[3 … len−?]: the TX computes `crc(&p[3], len−4)` and the RX checks `crc(&p[3], L−5)`.
    It therefore covers bytes 3 up to the last payload byte, and **excludes the length byte and the two ID bytes**.
  - The table starts 0x0000, 0x1189, 0x2312, 0x329B, 0x4624, 0x57AD, 0x6536, 0x74BF. That is the reflected-0x1021
    (KERMIT) table, but it is applied MSB-first: `crc = (crc<<8) ^ T[(crc>>8) ^ byte]`. T is built from a 16-entry
    nibble table: T[v] = Tshort[v & 15] ^ 0x1081·(v>>4). The seed is 0.
  - **The brief's "TX-ID-derived seed in v2" is not supported.** Multi calls the same `FrSkyX_crc(…, init=0)` for v1
    and v2, and no CRC code references the TX ID ([MPM-C], [MPM-X]) [U for real FrSky v2 firmware].
- **v2 adds** ([MPM-RX]): the CC2500 hardware CRC (0x8005/0xFFFF) appended after the software CRC.
- **v2 bind:** bytes 3…len−2 XORed with **0xA7** ([MPM-X], [MPM-RX]) [GPL-only].
- **Bind packet v1:** `1D 03 01 ID_lo ID_hi idx h0..h4 rx_tx_addr[1] RXnum` and the flag pattern 55 AA 5A A5.

### Hopping
- **Rule:** idx = (idx + chanskip) mod 47 after every frame. chanskip is 1–46, random, and is sent in bytes 4–5
  ([MPM-X]).
- **v1:** the TX sends its 47-entry table in bind packets, with the same ID-derived generator as D8 ([MPM-C]).
- **v2 derives the table from the TX ID** ([MPM-C] `FrSkyX2_init_hop`) [GPL-only]:
  1. id = ID_hi·256 + ID_lo. inc = id % 46 + 1, bumped to 13 or 36 if it is 12 or 35. offset = id % 5.
  2. ch_i = 5·((inc·i) % 47) + offset, for i = 0..46.
  3. LBT adds +2 to channels {0, 1, 43, 44, 87, 88, 129, 130, 173, 174} and +3 to {216, 217, 218}.
  4. FCC adds +2 to channels {3, 4, 46, 47, 90, 91, 133, 134, 176, 177, 220, 221}.
  5. Entry 47 = 0 (the bind channel).
- **Consequence:** in v2 the visible channel set is mostly congruent to (id % 5) mod 5, and the spacing between
  hops follows inc·chanskip. **The TX ID can be cross-checked against the observed channel set.**

### Timing
- **Period:** 9 ms. FCC: TX, then 5.2 ms, RX 3.4 ms, 0.4 ms. LBT adds a 0.4 ms RSSI/LBT check ([MPM-X]).
- **Telemetry slot:** the RX answers with a 17-byte frame (length 0x0E? [U]) about 5.2 ms after the TX frame, on the
  same channel ([MPM-X]).
- **Air time:** about 36 bytes, which is about 4.1 ms at 70 kb/s [calc].

### Passive identity
TX ID (16 bit plus byte 3), RX number (byte 6), chanskip (11 bits split over bytes 4–5), hop index, failsafe and
channel values; and v1 against v2 (from the data rate on FCC: 70 vs 77 kb/s; HW CRC present or absent).

## 3. FlySky AFHDS 2A (A7105)

### A7105 primitives [DS] ([A7105])
- **Packet (FIFO mode):** `preamble (1–4 B) | ID code (2 or 4 B) | payload | CRC (2 B)`.
  - Preamble: 0101… if the first ID bit is 0, otherwise 1010…. ID byte 0 is recommended to be 0x5? or 0xA?
    (§16.1, p.65).
  - **CRC:** "CCITT-16" over the **payload only** (preamble and ID excluded). The datasheet writes the polynomial as
    x^16+x^15+x^2+1, which is the 0x8005 form, contradicting the "CCITT" label [U: the polynomial and init need
    capture verification] (§16.2, p.66).
- **FEC:** (7,4) Hamming per nibble. It covers the payload and the CRC.
- **Whitening:** XOR with PN7, seed WS[6:0] in register 0x21.
- **Code Register I (0x1F):** bit 6 MCS, bit 5 WHTS, bit 4 FECS, bit 3 CRCS, bit 2 IDL (1 = 4-byte ID), bits [1:0]
  PML (11 = 4-byte preamble) (p.24–25).
- **Data rate:** = F_SYCK / 32 / (SDR+1), with SDR in register 0x0E. F_SYCK = F_MCLK / {1, 2, 2, 4} by CSC (register
  0x0D) (p.19–20).
- **LO:** f_LO = f_LO_BASE + CHN·F_CHSP, with F_CHSP = F_XTAL·(DBL+1)/4/(CHR+1); 500 kHz is recommended. With a
  16 MHz crystal, DBL = 1, BIP = 0x4B and BFP = 0x0002, the base is about 2400.001 MHz (ch.14, p.56–57).
- **Modulation:** register 0x14 bit 3 FS enables a Gaussian filter (BT 0.7). Deviation comes from register 0x15
  FD[4:0]: FD = 0x0B at F_PFD 32 MHz gives **186 kHz** for rates above 50 kb/s (p.21).

### AFHDS 2A register values and their decoding
Register values are [GPL-only] ([MPM-A7]); decoded [calc] with the formulas above.
- **0x0D = 0x05:** CSC = 01, so F_SYCK = 16 MHz. **0x0E = 0x00 gives 500 kb/s.**
- **0x14 = 0x16:** FS = 0, so **plain 2-FSK, no Gaussian filter** [U: some sources call it GFSK].
- **0x15 = 0x2B:** FD = 0x0B, so **±186 kHz**.
- **0x10 = 0x9E:** DBL = 1, CHR = 0xF, so the channel step is **500 kHz**.
- **0x11/0x12/0x13 = 4B 00 02:** base 2400.001 MHz.
- **0x1F = 0x1F:** FEC **on**, CRC **on**, **whitening off**, 4-byte ID, 4-byte preamble.
- **0x20 = 0x1E:** ID error tolerance ETH = 3 bits; PMD = 8-bit preamble detect.
- **ID code:** `A7105_WriteID(0x5475C52A)`, so the ID bytes 0x54 0x75 0xC5 0x2A are sent in that byte order
  ([MPM-A7]) [U: the bit order within each byte is assumed MSB-first].
- **AFHDS 1 (FlySky):** uses the same ID, with 0x1F = 0x0F (CRC on, no FEC) ([MPM-A7]).

### Frame
TX payload 38 bytes; RX telemetry 37 bytes ([MPM-A2]) [GPL-only].
- **On air** [calc]: 32 bits of preamble + 32 bits of ID + (38+2)·2 nibbles·7 bits = 624 bits, which is
  **1.25 ms at 500 kb/s**.
- **Bytes:** [0] type, [1..4] **TX ID (4 B)**, [5..8] **RX ID (4 B; 0xFF×4 while binding)**, [37] = 0.
- **Packet types:**
  - **0x58:** sticks. 14 channels as 12-bit little-endian at [9+2k] and [10+2k] (low nibble). Channels 15–17 sit in
    the spare nibbles; channel 17 = 2000 − 10·LQI. The next hop index is in byte 34 (high nibble) and byte 36.
  - **0x56:** failsafe.
  - **0xAA:** settings: [9] = 0xFD, [10] = 0xFF, [11–12] frame rate 50–400 Hz, [21] 0xDD (SBUS) or 0xDE (IBUS).
  - **0xBB / 0xBC:** bind. Bytes [11..26] carry the 16 hop channels.
  - **RX telemetry** 0xAA / 0xAC uses 4-byte sensor entries from [9] (0x00 RX V, 0xFC RSSI, 0xFE error rate, 0xFF end).

### Hopping and timing [GPL-only] ([MPM-A2])
- **Hop table:** 16 channels, generated from the TX ID with the LCG rnd = rnd·0x0019660D + 0x3C6EF35F (the public
  Numerical-Recipes LCG).
  - ch = band·41 + 1 + ((rnd >> idx) % 41), with the band chosen from the index and ID byte 3.
  - Each channel is 1–164, which is **2400.5–2482 MHz** [calc].
  - Each new channel must be ≥ 5 channels (2.5 MHz) from every other.
- **Rule:** the index advances +1 mod 16 per packet [U: Multi advances every packet].
- **Period:** **3.85 ms** per packet.
- **Other packets:** a settings packet every 1313 packets (about 5 s); failsafe every 1569 packets.
- **Bind:** on A7105 channels 0x0D (2406.5 MHz) and 0x8C (2470.0 MHz), alternating.

### Passive identity
4-byte TX ID and 4-byte RX ID in clear; packet type; frame-rate and output settings (from 0xAA); 14–17 channels; and
the RX telemetry (RX voltage, RSSI).

## 4. XN297 / XN297L toy drones and the Bayang protocol

### PHY
- **Modulation:** GFSK, 2.400–2.483 GHz; 250 kb/s, 1 Mb/s or 2 Mb/s. The crystal must be ±40 ppm for 1/2 Mb/s and
  ±20 ppm for 250 kb/s (Panchip XN297L datasheet:
  https://wmsc.lcsc.com/wmsc/upload/file/pdf/v2/lcsc/1809291411_Panchip-Microele-XN297LCQH_C88024.pdf) [DS].
- **Deviation and BT:** not found [U]. Emulation on the nRF24L01 works, which suggests nRF-like ±160 kHz (1 Mb/s)
  and BT 0.5 [U].
- **Channel:** 2400 + RF_CH MHz (nRF24 convention, used by the emulation) [U].
- **Data rate in use:** Bayang uses 1 Mb/s ([MPM-BY]).

### On-air format against the nRF24L01 [GPL-only] ([MPM-XN], [DEV-NRF], [DEVF])
**XN297 preamble prefix:**
- The emulation loads the nRF TX_ADDR with `55 0F 71 0C 00` (that is, 0xC710F55, 28 bits). It sets SETUP_AW to
  len−2, so the bytes ...**0x0C 0x71 0x0F 0x55** appear on air before the XN297 address.
- The CC2500 path sends `0C 71 0F 55` explicitly.
- Forum reverse engineering gives the XN297 preamble as **0x710F55** ([DEVF]).

**Address (3–5 B):**
- The bytes are sent in **reverse order** of the protocol's `addr[]` array (addr[len−1] first).
- Each byte is XORed with `scramble[i]`; **there is no bit reversal of the address.**
- For RX: rx_addr[i] = addr[i] ^ scramble[len−1−i].

**Payload:**
- Each byte is **bit-reversed** (the XN297 sends LSB-first, while the nRF sends MSB-first), then XORed with
  `scramble[addr_len + i]`.

**Scramble table:**
- The table is 35–37 bytes, starting `E3 B1 4B EA 85 BC E5 66`.
- A forum post generates it with a 7-bit LFSR seeded with 0x71 ([DEVF]) [U: the polynomial is not given; x^7+x^4+1
  is likely].

**CRC-16:**
- Polynomial 0x1021, MSB-first, **init 0xB5D2** ([DEV-NRF]).
- It covers the scrambled address bytes (in on-air order) and the scrambled payload.
- The result is then XORed with `xorout[addr_len−3+payload_len]`.

**xorout tables:**

| Mode | First entries |
|---|---|
| scrambled | 0x0000, 0x3448, 0x9BA7, 0x8BBB, 0x85E1, 0x3E8C, 0x451E, 0x18E6 … |
| unscrambled | 0x0000, 0x3D5F, 0xA6F1, 0x3A23, 0xAA16, 0x1CAF, 0x62B2, 0xE0EB … |
| enhanced scrambled | 0x0000, 0x7EBF, 0x3ECE, 0x07A4 … |

**Enhanced mode (ESB-like):**
- The 9-bit PCF is length(6), PID(2), NO_ACK(1), placed after the address. The payload is then shifted by 2 bits
  (not byte-aligned).
- The CRC covers the address, PCF, payload and 2 extra bits.
- The xorout comes from the enhanced table.

**Implementation note:** `scramble` and `xorout` are data tables, not code. They are [GPL-only] facts but can be
regenerated independently: the scramble table from the LFSR; xorout[k] is the CRC of k zero bytes under scrambling
against without it [U]. **Validate both on a real capture.**

### Bayang [GPL-only] ([MPM-BY])
- **Settings:** XN297, 1 Mb/s, scrambled, CRC on, **5-byte address**, **15-byte payload**.
- **On air:** 4 + 5 + 15 + 2 = 26 bytes, which is about **208 µs** [calc].
- **Bind:** address 00 00 00 00 00 on RF channel 0 (2400 MHz). Channel 10 for the X16_AH and IRDRONE variants.
- **Bind packet:**
  - [0] = 0xA4 (or A1–A3, 0x53 for QX100).
  - [1..5] = **TX address = rx_tx_addr[0..4]**.
  - [6..9] = 4 hop channels.
  - [10..12] = rx_tx_addr[0], [1], [2] (variant-specific constants instead for some variants).
  - [13] = 0x0A.
  - [14] = checksum.
- **Data packet:**
  - [0] = 0xA5 (0xA6 for some variants).
  - [1] = 0xFA.
  - [2] = flags: flip 0x08, RTH 0x01, picture 0x20, video 0x10, headless 0x02.
  - [3] = flags: inverted 0x80, take-off 0x20, emergency stop 0x0C.
  - [4..11] = A, E, T, R as 10-bit big-endian pairs; the high byte carries a trim offset of 0x7C.
  - [12] = rx_tx_addr[2].
  - [13] = 0x0A.
  - **[14] = 8-bit sum of [0..13].**
- **Hopping:** 4 channels: [0] = 0, [1] = (rx_tx_addr[3] & 0x1F) + 0x10, [2] = [1] + 0x20, [3] = [2] + 0x20. The
  channel advances on every packet.
- **Period:** 2 ms per packet; 5 ms with telemetry.
- **Telemetry:** frames start with 0x85.

### Passive identity
- **XN297 address:** the 5-byte data address equals the TX ID, after descrambling.
- **Bayang:** the TX ID also appears in clear in bind packets [1..5] and in data byte [12].
- **Hop set:** the 4-channel set follows from ID byte 3.
- **Sticks and flags:** readable.
- **Protocol family:** the header byte (A4/A5).

## 5. FLARM (868 MHz SRD; India)

### PHY
- **Frequencies:** Europe/Africa/Asia: **868.2 and 868.4 MHz**, one transmission per second on each. 868.2 MHz is
  used at 0.4–0.8 s in the second and 868.4 MHz at 0.4–1.2 s ([OXF] §2.2.1). The Americas, Oceania and Israel use an
  undisclosed hopping scheme ([OXF]).
- **India:** SoftRF's band plan "IN" uses **866.0 MHz, 1 channel, 200 kHz** for the Legacy/OGNTP/ADS-L protocols
  ([SRF-FP]) [U: FLARM's own India band is unconfirmed].
- **Modulation:** 2-FSK/GFSK, **±50 kHz**. The chip rate is **100 kb/s, Manchester-coded, so 50 kb/s of data**
  ([OXF], [SRF-L], [OGN-FAQ]). OGN lists "50 kbps GFSK, 250 kHz BW, 5 ms packets" ([OGN-FAQ]). FLARM is based on the
  nRF905 ([OXF]).
- **Not 38.4 kb/s:** that rate is PilotAware at 869.525 MHz.

### Frame
**Preamble:** nRF905-style, Manchester (`1010…`). SoftRF sends 1 byte of 0x55 ([SRF-L]).

**Sync word:**
- SoftRF lists the sync as the 8 Manchester bytes `55 99 A5 A9 55 66 65 96` ([SRF-L]) [GPL-only].
- Decoding these with 01 → 1 and 10 → 0 [calc] gives **F5 31 FA B6**.
- This matches [OXF]: "the first three bytes of the header statically contain **0x31FAB6**".
- **SoftRF marks the payload "inverted"** (RF_PAYLOAD_INVERTED).

**Header** (after 31 FA B6; [OXF] Tables 2–3):
- **24-bit radio ID, little-endian.**
- **A type byte:** 0x20 = default (manufacturer) ID, 0x10 = ICAO 24-bit address.
- Example: `31:fa:b6 | 8b:08:de | 20 | 7e:47…`.

**Payload:**
- 24 bytes. The v7 header in clear is addr:24, type:4, addr_type:2, urgency:2, …; the encrypted payload is 16
  bytes or more ([SRF-L]).
- The decrypted v6 fields are vs, gps, aircraft_type, lat:19, alt:13, lon:20, ns[4]/ew[4] speed predictions
  ([SRF-L]; [OXF] Table 3).

**CRC:** CRC-16 CCITT, init 0xFFFF, 2 bytes ([SRF-L] `RF_CHECKSUM_TYPE_CCITT_FFFF`) [GPL-only; U: coverage].

**Air time:** SoftRF uses about 6 ms. The TX interval is randomised over 600–1400 ms ([SRF-L]).

### Obfuscation [GPL-only]
- **Order:** XXTEA encryption, then Manchester coding, then FSK ([OXF] Fig. 2).
- **Algorithm:** XXTEA (Corrected Block TEA, delta 0x9E3779B9, 6 rounds in SoftRF; [OXF] says 8).
  - v6 keys come from a static table, the UTC time (time >> 6, with bit 23 choosing the table half), and the
    address.
  - v7 (2017 onward) uses another key set, mixed from the header words and time >> 4 ([SRF-L], [OXF] §2.3).
- **Key history:** keys were changed in 2008, 2015 and 2017, and were each time reverse-engineered or leaked ([OXF]).
- **Licence and legal note:** the key constants appear only in GPL code. **Do not embed them in the C++ library.**
  Decrypting FLARM may also raise legal and ToS issues [U].

**Readable without decryption:** the **24-bit aircraft or device ID**, the **address type** (ICAO or FLARM), and the
v7 header fields (message type, stealth / no-track bits [U]). Also the transmission timing pattern (1 Hz, two
frequencies) and RSSI.

### Related
ADS-L SRD-860 (EASA): 868.2/868.4 MHz, GFSK, 100 kchip/s Manchester = 50 kb/s, sync 0x72 0x4B (Manchester
encoded), XXTEA 6 rounds ([EASA]). It is open and specified, so it is the cleanest path for an interoperable 868
detector.

### Open-source decoders and licences
- **SoftRF** (Legacy/FLARM, OGNTP, ADS-L, FANET, P3I): **GPL-3** [U: licence per the repo]
  (https://github.com/lyusupov/SoftRF).
- **OGN glidernet:** `python-ogn-client` and `ogn-python` are **AGPL-3**. `ogn-rf` shows no licence on GitHub
  (https://github.com/glidernet).
- **rtlsdr-ogn:** the receiver that decodes FLARM; distributed as binaries [U].
- **GXAirCom:** FLARM/FANET/OGN (https://github.com/ElDevelish/GXAirCom), licence not checked [U].
- **Usage:** all of these are test oracles only.

## 6. IEEE 802.15.4 2.4 GHz O-QPSK (XBee 802.15.4, Zigbee telemetry) [DS]

### PHY ([CC2420] §10, p.24–25)
- **Spreading:** each byte becomes 2 symbols, low nibble first. Each 4-bit symbol maps to one of 16 32-chip
  sequences. **2 Mchip/s**, so **62.5 ksym/s and 250 kb/s**.
- **Chip order:** C0 is sent first.
- **Modulation:** O-QPSK with **half-sine** chip shaping, which is MSK-equivalent. Even chips go on I and odd chips
  on Q, with Q offset by Tc = 0.5 µs.
- **Channels:** 11–26, f = 2405 + 5·(k−11) MHz ([CC2420] frequency programming section; IEEE 802.15.4).

### Chip table (C0…C31) [DS: IEEE 802.15.4, reproduced in [CC2420] Table 3]
```
0  11011001110000110101001000101110
1  11101101100111000011010100100010
2  00101110110110011100001101010010
3  00100010111011011001110000110101
4  01010010001011101101100111000011
5  00110101001000101110110110011100
6  11000011010100100010111011011001
7  10011100001101010010001011101101
8  10001100100101100000011101111011
9  10111000110010010110000001110111
10 01111011100011001001011000000111
11 01110111101110001100100101100000
12 00000111011110111000110010010110
13 01100000011101111011100011001001
14 10010110000001110111101110001100
15 11001001011000000111011110111000
```
Symbols 1–7 are 4-chip cyclic shifts of symbol 0. Symbols 8–15 are symbols 0–7 with the odd chips inverted
(conjugated).

### Frame ([CC2420] §16, p.36–38)
- **SHR:** preamble of **4 bytes of 0x00 (8 zero symbols)**, then **SFD 0xA7** (sent as symbols 7, then A).
- **PHR:** 1 byte. Bits 0–6 are the frame length (MPDU bytes including the FCS, at most 127); bit 7 is reserved.
- **MPDU:**
  - **FCF (2 B, LE):** bits 0–2 frame type (0 beacon, 1 data, 2 ACK, 3 MAC command), 3 security, 4 frame pending,
    5 ACK request, 6 PAN-ID compression ("intra-PAN"), 10–11 destination address mode, 12–13 frame version,
    14–15 source address mode. Address modes: 0 = none, 2 = 16-bit short, 3 = 64-bit extended.
  - **Sequence number:** 1 B.
  - **Addressing:** destination PAN, destination address, source PAN (omitted if intra-PAN), source address. All LE.
  - **Payload.**
  - **FCS.**
- **FCS:** **CRC-16 x^16+x^12+x^5+1**, computed over the MPDU (the length byte excluded), with data fed LSB-first
  ([CC2420] p.38, Fig. 20). The IEEE register init is **0x0000** and there is no xorout, so this is CRC-16/KERMIT:
  reflected, poly 0x8408 [U: the init value is from the IEEE standard, not printed on the pages read].
- **ACK frame:** 5 bytes (FCF 0x0002, seq, FCS).

### XBee specifics
- **XBee 802.15.4 and Zigbee:** use this PHY. The 64-bit extended addresses of Digi modules start with the OUI
  **00 13 A2** (SH = 0013A200…) [U: Digi docs convention, not re-verified here].
- **Zigbee:** a NWK header follows the MAC header, with 16-bit NWK addresses and an 8-byte extended PAN ID in
  beacons.
- **XBee-PRO 900HP (902–928 MHz):**
  - **Rate:** 10 kb/s or 200 kb/s RF; sensitivity −101 dBm at 200k and −110 dBm at 10k; up to 24 dBm ([XB9]).
  - **Hopping:** **FHSS**. **64 channels**; the module picks 25 by noise scan, with 8 hop patterns. The hop sync is
    carried in a packet preamble (retailer and forum sources, via [XB9]) [U].
  - **Not published:** the modulation (2-FSK or GFSK?), the channel spacing (≈400 kHz?) and the framing [U]. The
    radio is a SiLabs-based proprietary DigiMesh PHY.
  - **Not relevant for India.** 902–928 MHz is not an Indian ISM band; India has 865–867 MHz.
- **XBee-PRO 868:**
  - **24 kb/s RF**, 10 % duty cycle, single channel, −112 dBm, up to +25 dBm ([XB8]).
  - **Frequency:** 869.525 MHz (g3 sub-band) per a reseller [U].
  - **Modulation:** not published [U].

---

## Implementation notes for a detector

General approach: use the block FFT / channelizer to find narrowband bursts. Then FM-demodulate (arg of
x[n]·x*[n−1]) for the FSK systems, or MSK/O-QPSK-correlate for 802.15.4. Slice at the expected symbol rate,
correlate the sync, and check the CRC. A CRC pass is the label gate. Without one, report "-like".

**1. FrSky D8 (CC2500)**
- **Search:** 2404.0 MHz + k·299.93 kHz (k = 0–234). Burst about 6.7 ms. GFSK ±31.7 kHz at **31.04 kbaud**.
- **Sync:** 4 bytes of 0xAA preamble, then the 16-bit sync **0xD391**.
- **Fields:** length byte 0x11. Bytes 1–2 = TX ID.
- **CRC:** CRC-16 **0x8005, init 0xFFFF, MSB-first**, over the length and 17 bytes [DS].
- **Timing gate:** about 9 ms repetition, a 36 ms cycle, a 47-channel hop. Telemetry replies sit in the gaps.

**2. FrSky D16**
- **Search:** same grid. FCC: GFSK ±57 kHz at 70.0 kbaud (v1) or 76.97 kbaud (v2). LBT: 2-FSK ±70 kHz at 100 kbaud,
  base 2405.0 MHz.
- **Sync:** **32-bit 0xD391D391**. Length 0x1D (or 0x20 for LBT v1).
- **CRC:** the software CRC described in §2 (KERMIT table applied MSB-first, seed 0) over bytes 3…N−2 [GPL-only:
  verify on a capture]. v2 also has the HW CRC 0x8005/0xFFFF after it.
- **Fields:** TX ID bytes 1–2, chanskip bits from bytes 4–5, RX number byte 6.
- **Timing gate:** 9 ms period; hop index +chanskip mod 47. For v2, check ID ↔ channel set (5·((inc·i)%47) + id%5).

**3. AFHDS 2A (A7105)**
- **Search:** 2400.0 + 0.5·ch MHz (ch 1–164). 2-FSK ±186 kHz at **500 kbaud**. Burst 1.25 ms every **3.85 ms**.
- **Sync:** 32 bits of preamble 0101…, then the **ID code 54 75 C5 2A** (allow up to 3 bit errors, like the RX
  setting).
- **Decode:** (7,4) Hamming on each nibble of 40 bytes, giving 38 B payload + CRC. Check the CRC-16 (poly and init
  [U]; try 0x8005 and 0x1021, with init 0x0000 and 0xFFFF, and keep the one that validates on testdata).
- **Fields:** byte 0 ∈ {0x58, 0x56, 0xAA, 0xBB, 0xBC}. TX ID = bytes 1–4; RX ID = bytes 5–8.
- **Timing gate:** 16-channel hop set with ≥ 2.5 MHz separation.

**4. XN297 / Bayang**
- **Search:** 2400 + n MHz. GFSK at 1 Mbaud (also 250k and 2M). Burst about 0.21 ms every 2 ms, over 4 channels
  spaced 32 MHz.
- **Sync:** correlate **…0C 71 0F 55** (MSB-first). Then XOR the next 5 bytes with E3 B1 4B EA 85 and reverse their
  order to get the address. Descramble the payload with the table continuing from index 5, and bit-reverse each
  byte.
- **CRC:** CRC-16 0x1021, init 0xB5D2, over the on-air (scrambled) address and payload, XORed with xorout[17]
  [GPL-only table: regenerate and verify].
- **Bayang checks:** payload[0] ∈ {A4, A5, …}, and the 8-bit sum. Hop channels = (ID3 & 0x1F) + 0x10 + {0, 32, 64}.

**5. FLARM**
- **Search:** 868.2 and 868.4 MHz (India: 866.0 MHz [U]). 2-FSK ±50 kHz at 100 kchip/s Manchester. Burst about
  5–6 ms, about 1/s per frequency.
- **Sync:** correlate the chip pattern `55 99 A5 A9 55 66 65 96` (this is F5 31 FA B6 after Manchester decoding,
  01→1). Try both polarities.
- **CRC:** CRC-16 CCITT init 0xFFFF [GPL-only; U: coverage].
- **Readable fields:** the 24-bit ID (LE) and the address type from the clear header. **Do not decrypt** (key
  material is GPL-only and legally sensitive).

**6. IEEE 802.15.4**
- **Search:** 2405 + 5(k−11) MHz. 2 Mchip/s O-QPSK half-sine; about 2 MHz occupied.
- **Sync:** correlate the 32-chip symbol 0 sequence (MSK-demodulate, then match against the chip table), then
  SFD = symbols 7, A.
- **Frame:** PHR length ≤ 127.
- **CRC:** CRC-16/KERMIT (poly 0x1021 reflected, init 0, LSB-first) over the MPDU [DS].
- **Fields:** FCF type, PAN ID, short or extended addresses. An extended source address starting 00:13:A2 suggests
  an XBee [U].
- **XBee 900HP / 868:** no published PHY. Report as an FHSS-FSK "-like" detection only (energy plus hop statistics).

## Licence status summary

| Item | Status |
|---|---|
| CC2500 formulas, sync 0xD391, HW CRC 0x8005/0xFFFF, PN9 | [DS] |
| A7105 packet / FEC / CRC / whitening / LO formulas | [DS] |
| 802.15.4 chip table, SFD, FCS | [DS] (IEEE, CC2420) |
| FrSky register values, packet layouts, hop generators, D16 software CRC, 0xA7 bind XOR | [GPL-only] (Multiprotocol GPL-3) |
| AFHDS2A register values, ID 0x5475C52A, layout, hop LCG use | [GPL-only] (the LCG constants themselves are public) |
| XN297 scramble / xorout tables, CRC init 0xB5D2, Bayang layout | [GPL-only] (Multiprotocol / Deviation GPL-3; forum reverse engineering) |
| FLARM sync, layout, CRC, XXTEA keys | [GPL-only] (SoftRF GPL-3) plus an academic paper ([OXF], usable as a description). **The keys must not be used.** |
| ADS-L SRD-860 | [DS] (EASA spec) |
| XBee 900HP / 868 PHY details | not published [U] |
