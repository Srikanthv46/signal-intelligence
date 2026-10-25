# 08 — TBS Crossfire, Spektrum DSM2/DSMX, Walkera DEVO (air links for passive detection)

Research date 2026-10-21. Web research only. `[U]` = uncertain or inferred. `[GPL-only]` = the constant is published only in
GPL code, so for the C++ library it must come from a vendor document or be re-derived (for example from captures).
Values marked "computed" are arithmetic on the cited values, not quoted.

Key sources (abbreviations used below):
- **G3** — g3gg0, "[FPV] Analysis of TBS Crossfire, reverse engineering the air link", 2021-08-28 (updated 2025-12-01): https://www.g3gg0.de/default/fpv-analysis-of-tbs-crossfire/
- **G3-INO** — g3gg0 ESP32 Crossfire sniffer (MIT licence): https://github.com/g3gg0/ESP32_CRSFSniffer (`ESP32_CRSFSniffer.ino`, `SX1276.cpp`)
- **G3-BRUTE** — g3gg0 CRC brute-forcer: https://github.com/g3gg0/CRSF_RevEng/blob/master/shell_crcbrute/brute.c (no licence seen [U])
- **HAD** — Hackaday on G3: https://hackaday.com/2022/02/20/reverse-engineering-a-900-mhz-rc-transmitter-and-receiver/
- **SX1276-DS** — Semtech SX1276/77/78/79 datasheet (register meanings): https://www.semtech.com/products/wireless-rf/lora-connect/sx1276
- **CY-DS** — Cypress CYRF6936 datasheet, doc 38-16015 Rev *J (2011-04-18): https://datasheet.octopart.com/CYRF6936-40LTXC-Cypress-Semiconductor-datasheet-14183427.pdf (also https://www.farnell.com/datasheets/1487446.pdf)
- **DEV-DSM** — Deviation `dsm2_cyrf6936.c` (GPL-3+): https://github.com/DeviationTX/deviation/blob/master/src/protocol/dsm2_cyrf6936.c
- **DEV-DEVO** — Deviation `devo_cyrf6936.c` (GPL-3+): https://github.com/DeviationTX/deviation/blob/master/src/protocol/devo_cyrf6936.c
- **MPM-DSM** — DIY-Multiprotocol `DSM_cyrf6936.ino` (GPL-3): https://github.com/pascallanger/DIY-Multiprotocol-TX-Module/blob/master/Multiprotocol/DSM_cyrf6936.ino
- **MPM-CYRF** — DIY-Multiprotocol `CYRF6936_SPI.ino` (GPL-3): https://github.com/pascallanger/DIY-Multiprotocol-TX-Module/blob/master/Multiprotocol/CYRF6936_SPI.ino
- **GRDSMX** — Morin & Cardoso, gr-dsmx-rc (GPL-3+): https://github.com/lscardoso/gr-dsmx-rc (`lib/Despreader_impl.cc`); paper "Hacking the DSMx Drone RC protocol", Proc. GNU Radio Conf. v2 n1: https://pubs.gnuradio.org/index.php/grcon/article/view/97 (TLS certificate expired when fetched; paper content known only from a search snippet)
- **SPM** — Spektrum DSMX page: https://www.spektrumrc.com/spm-bs-dsmx.html

---

## A. TBS Crossfire (868/915 MHz)

### A1. 150 Hz FSK mode — PHY and packet

| Item | Value | Source |
|---|---|---|
| Radio | Semtech SX1272 (TX/RX), FSK mode of the LoRa modem | G3, HAD |
| Bit rate | 85.1 kbit/s ("85.1 kBaud"); SX127x reg BitRate = 0x0178 = 376 → 32 MHz/376 = 85.106 kbit/s (computed) | G3; G3-INO `SX1276.cpp` initRegs |
| Deviation | 42.48 kHz; reg Fdev = 0x02B8 = 696 × 61.035 Hz = 42.48 kHz (computed). Sniffer sketch also uses `cfgFreqDev 42300`, `cfgBitrate 85000` | G3; G3-INO |
| Modulation index | h = 2·42.48/85.1 ≈ 1.0 (computed): 2-FSK, tones at ±42.5 kHz | computed |
| Gaussian shaping | Not published for the TX. The sniffer writes RegPaRamp 0x0A = 0x49 → bits 6:5 = 10 = Gaussian BT 0.5 on SX1276, but that is the sniffer's own config [U] | G3-INO; SX1276-DS |
| Preamble | Sniffer RX config RegPreamble = 0x0005 (5 bytes); preamble polarity bit = 0 (0xAA…) [U: TX preamble length not measured] | G3-INO `SX1276.cpp` |
| Sync word | RegSyncConfig 0x27 = 0x12 → SyncOn = 1, SyncSize field = 2 → **3 bytes: 0x3E 0x6B 0x3E** (regs 0x28-0x2A; 0x2B-0x2F = 0x55 unused). Taken from g3gg0's working sniffer; whether it is the same for all units/binds is not stated [U] | G3-INO `SX1276.cpp`; SX1276-DS (RegSyncConfig) |
| Packet format | RegPacketConfig1 0x30 = 0x48 → fixed length (bit7 = 0), DcFree = 10 = **whitening on** (SX127x CCITT PN9 whitening), CrcOn = 0 (CRC done in software), CrcAutoClearOff = 1, no address filter. RegPacketConfig2 0x31 = 0x40 → packet mode | G3-INO `SX1276.cpp`; SX1276-DS |
| Payload length | Uplink (TX→RX) **23 bytes**, downlink (RX→TX telemetry) **13 bytes** over the air | G3; G3-INO (`currentLength = downlink ? 0x0D : 0x17`) |
| Uplink CRC | CRC-16, poly 0x1021, refin = true, refout = true, xorout 0x0000 (KERMIT-like) with an unknown init; over bytes 0..20, stored little-endian in bytes 21-22 | G3; G3-BRUTE (`brute_crc16`) |
| Downlink CRC | CRC-8, poly 0x07, refin/refout = false, xorout 0x00, unknown init; over bytes 0..11, CRC in byte 12 | G3; G3-BRUTE (`brute_crc8`) |
| "Seeded per hop" | Identical payloads gave different CRCs depending on the hop number; the init "depends on the hop number" and repeats every hop cycle ("salted CRC"). Derivation not solved. G3-BRUTE tried a model where a 32-bit secret is prepended to the data with the hop number added to its low byte (hop numbers 61-66 in the test set); no result published [U] | G3; G3-BRUTE |
| Bind / firmware update | LoRa mode, CRC-8 with init 0 | G3 |

Uplink packet layout (23 bytes; G3, G3-INO header comment and `parseChannels()`):
- byte 0: low nibble = packet type (3 = normal 150 Hz FSK uplink; 1 = "50 Hz mode uplink data?" [U]); bit 5 (0x20) set = this packet carries sticks 8-11 (12-channel mode alternate packet).
- bytes 1-10: eight 10-bit stick values packed LSB-first (80 bits). 8-ch mode: sticks 0-7 every packet; 12-ch mode: sticks 0-3 every packet, 4-7 and 8-11 alternate at half rate.
- byte 11: high nibble = length of an extra (telemetry/serial) payload, low nibble = flags [U].
- bytes 12-20: extra payload (tunnelled CRSF frames; addresses `EA` = transmitter, `EC` = receiver seen) [U].
- bytes 21-22: CRC-16 (little-endian).
- Stick range seen 0xA1…0x2BB (G3-INO `chanPercent()`). Over-the-air resolution 10 bit, serial CRSF 11 bit (G3).
- Example frame (G3-INO): `03 AA FC F7 1F 80 A1 84 12 4A 28 00 00 00 00 00 00 00 00 00 00 32 E2`.

Downlink packet layout (13 bytes; G3-INO):
- byte 0: low nibble type (3 normal, 0x0B carries telemetry payload); high-nibble flags incl. 0x80 = "metadata cycle start" and a flag for "no uplink received for 8 times".
- bytes 1-2: rotating metadata (RSSI/LQ etc.) [U]. bytes 3-11: telemetry payload. byte 12: CRC-8.
- Example: `03 75 64 00 00 00 00 00 00 00 00 00 0F`.

### A2. Frequency plan and timing

Channel grid (G3 table; f = Fmin + n·spacing is consistent with the table — computed):

| Band / mode | Fmin (MHz) | Fmax (MHz) | Channels | Spacing |
|---|---|---|---|---|
| 868 | 860.165 | 885.905 | 50 uplink + 50 downlink (n = 0..99) | 260 kHz |
| 915 | 902.165 | 927.905 | 50 + 50 | 260 kHz |
| 868 RACE | 860.165 | 885.905 | 100 shared | 260 kHz |
| 915 RACE | 902.165 | 927.905 | 100 shared | 260 kHz |
| 868 CE | 863.093 | 868.581 | 50 | 112 kHz |
| 915 AU | 915.165 | 927.905 | 50 | 260 kHz |

- Two hardware setups differed by ~15 kHz in absolute frequency (G3): per-unit crystal offsets of ~17 ppm occur (computed).
- Normal mode: uplink on channel n ∈ 0..49, downlink on **n + 50** (+13.0 MHz at 260 kHz spacing) (G3; G3-INO `chanAbs = chan + 50`). RACE: uplink and downlink on the same channel, hopping over 0..99 (G3).
- Hopping: after each uplink/downlink pair the TX moves to the next channel; sequence length **150 hops** (`hoppingPos %= 150`), so one cycle = 150 × 6.667 ms = 1.0 s (computed); each of the 50 channels appears ~3 times per cycle (150/50, computed; the sniffer stores up to 4 slots per channel) (G3-INO). G3 says the full sequence is swept 150 times, then repeats [U wording].
- Hop sequence generator: **not published**. g3gg0 read it from SPI logs and says it is "most likely only valid for my configuration and my bind key", i.e. probably bind-dependent [U] (G3, HAD). His sniffer learns it empirically: sit on one channel, measure its slot offsets modulo the 150-slot cycle, repeat for all 50 channels (G3-INO states 0-4).
- Frame timing: 150 Hz → one uplink/downlink exchange every **6.666 ms** (300 bursts/s); the downlink comes at most ~2.6 ms after the uplink ends (G3). Sniffer uses `packet_ms = 6666.67 µs × 1.00101` (≈ 6673 µs, crystal-corrected) and re-arms 2650 µs after an uplink for the downlink (G3-INO). The repeat of channel 0 is found by a 365-374 ms gap (G3-INO).
- Air time (computed, assuming 5-byte preamble + 3-byte sync): uplink (5+3+23)·8/85.1k ≈ 2.9 ms, downlink (5+3+13)·8/85.1k ≈ 2.0 ms.

### A3. 50 Hz LoRa mode
- Uses the SX1272 LoRa modem; G3 mentions the 50 Hz mode and a slower LoRa mode but did not analyse them. **SF, BW, CR, preamble, sync word and payload length were not found in any public source** [U]. A forum answer says LoRa is used at 50 Hz because LoRa cannot carry the packet at 150 Hz (https://intofpv.com/t-crsf-lora). TBS switches to 50 Hz on low LQ (TBS KB: https://www.team-blacksheep.com/new-support/kb/4000112664/). Packet type nibble 1 = "50 Hz mode uplink?" (G3-INO) [U]. Must be measured on the bench (blind LoRa SF/BW estimation on a capture).

### A4. Identity recoverable passively (no bind data)
- No address/ID field was identified in the 23/13-byte payloads (G3 found no authentication; no ID documented) [U].
- Usable fingerprints: (1) the per-hop CRC-init vector (16-bit uplink / 8-bit downlink per hop slot): recoverable from a received packet because a CRC is affine in its init (solve or brute force, as G3-BRUTE did), and bind-specific per G3's observation [U that it differs per bind]; (2) the 150-slot hop sequence (bind-dependent [U]); (3) carrier offset in ppm (units differ by ~15 kHz, G3); (4) band/mode (868/915/CE/AU, RACE vs normal) from the grid; (5) 8 vs 12 channel mode (byte 0 bit 5 alternation); (6) telemetry contents (CRSF frames) after de-whitening.

---

## B. Spektrum DSM2 / DSMX (Cypress CYRF6936)

### B1. CYRF6936 physical layer (vendor datasheet CY-DS unless noted)

| Item | Value | Source |
|---|---|---|
| Band | RF range 2.400-2.497 GHz; "98 discrete 1 MHz channels"; CHANNEL_ADR (0x00) 7-bit, default 0x48 (72). f = 2400 + CHANNEL MHz [U: the formula comes from drivers/TRM, not stated in the datasheet text] | CY-DS p4, p14, p17 |
| Modulation | GFSK; GFSK raw mode 1 Mbit/s; DSSS chips at 1 Mchip/s (implied: 32-chip 8DR = 250 kbit/s; framing symbol "32 or 64 µs") | CY-DS p4-6 |
| Deviation | 270 kHz typ (PN pattern 10101010) to 323 kHz typ (pattern 11110000); occupied BW 500 min, 876 typ kHz (−6 dBc) | CY-DS Table 8 p17 |
| Gaussian BT | Not given in the datasheet [U] (0.5 commonly assumed) | CY-DS |
| Data modes | GFSK 1000 kbps; 8DR: 32-chip 250 kbps, 64-chip 125 kbps ("eight bits are encoded in each derived code symbol"); DDR: 32-chip 62.5 kbps, 64-chip 31.25 kbps (2 bits/symbol); SDR 64-chip 15.625 kbps (1 bit/symbol). The length, payload and CRC are sent in the data mode; the SOP is separate | CY-DS "Data Transmission Modes", "Data Rates" p4, p6 |
| TX_CFG_ADR 0x03 | bit5 DATA CODE LENGTH (1 = 64 chip), bits4:3 DATA MODE, bits2:0 PA. Mode encoding from drivers [U]: 0x28 = 64-chip 8DR, 0x38 = 64-chip SDR, 0x08 = 32-chip 8DR, 0x00 = GFSK | CY-DS p14; DEV-DSM; DEV-DEVO; MPM-CYRF |
| Packet | Preamble (n × 16 µs) – SOP1 – SOP2 – Length (1 byte) – payload – CRC16 (2 bytes). SOP = two-symbol marker, each 32 or 64 chips; required in GFSK and 8DR, optional in DDR, not in SDR | CY-DS Figure 2 p5, "Link Layer Modes" p4 |
| SOP code | SOP_CODE_ADR 0x22, 8-byte register file (64 chips); default **0x17FF9E213690C782** | CY-DS p14 note 5 |
| Data code | DATA_CODE_ADR 0x23, 16-byte register file (two 64-chip codes); default printed as 0x02F9939702FA5CE3012BF1DB0132BE6EF (33 hex digits as printed: a datasheet typo [U]) | CY-DS p14 note 6 |
| Preamble | PREAMBLE_ADR 0x24, 3 bytes: 2-byte pattern + count; default **0x333302**; count must be > 4 for DDR, > 8 for SDR. DSM uses **0x333304** (pattern 0x3333, 4 × 16 µs = 64 µs, computed). The GFSK init in MPM uses 0xAAAA00 \| len | CY-DS p14 note 7; DEV-DSM; MPM-CYRF |
| Framing | FRAMING_CFG 0x10: SOP EN, SOP LEN, LEN EN, SOP threshold (5 bits). DSM data uses 0xEA (SOP on, 64-chip SOP, length on, threshold 0x0A) | CY-DS p14; DEV-DSM |
| Length | First byte after the SOP, sent at the payload data rate; EOP = length bytes + 2 CRC bytes | CY-DS p4 |
| CRC16 | "USB CRC polynomial" (x^16+x^15+x^2+1 = 0x8005) with programmable seed (CRC_SEED_LSB/MSB 0x15/0x16, default 0x0000); the RX checks against both the configured and the zero seed. GRDSMX computes it with poly 0x8005 over the length byte + payload. Reflection/xorout are not stated in the datasheet [U] | CY-DS p5, p14; GRDSMX |
| Buffers | 16-byte FIFOs, packets up to 255 bytes (max 40 bytes at < 60 ppm crystal offset; 16 for 64-DDR/SDR) | CY-DS p5, p18 |
| Bit order | PN code bytes are "in order transmitted (LSB 1st)" (comment in GPL code) [GPL-only] | DEV-DSM; MPM-CYRF |

How bits map to chips in 8DR (64-chip, as used by DSM): the 16-byte data code is two 64-chip codes A and B. Each 64-chip symbol is one of 4 base codes {¬A, A, ¬B, B} (2 bits) at one of 64 cyclic rotations (6 bits): byte = 64·index + rotation. GRDSMX decodes this way, accepting ≥ 54/64 chip matches (≤ 10 chip errors). The order of {¬A, A, ¬B, B} and the rotation direction come from GPL code only [GPL-only][U]. For 32-chip 8DR (DEVO): 8 bits per 32-chip symbol by an analogous derivation, not documented publicly [U]. The datasheet says only "eight bits are encoded in each derived code symbol" (CY-DS p4). An Infineon forum post quotes the TRM: "128 derivative PN codes are generated based upon the organization of the contents of the DATA_CODE_ADR file register" (search snippet of https://community.infineon.com/t5/Other-Technologies-General/8DR-support-in-CYRF-modules/td-p/125239; page returned 403) [U].
- SDR (Cypress patent US7653123B1, https://patents.google.com/patent/US7653123): a 1 is sent as the inverted PN code, a 0 as the code.

PN code table status:
- The datasheet contains **only the default SOP and DATA codes** (above). The **5 × 9 × 8-byte DSM `pncodes` table** appears in Deviation, Multiprotocol and gr-dsmx-rc (all GPL) with no source citation → **[GPL-only]** as far as found.
- The **10-row SOP table used by Walkera DEVO** (`DEVO_j6pro_sopcodes` in MPM-CYRF) is commented as following "the Cypress recommended table" in the **WirelessUSB LP/LPstar Technical Reference Manual** (vendor document; MPM notes one byte changed from the Cypress value). The TRM itself was not found online in this research [U]; if obtained, its table can be used directly, and it may also be the origin of the DSM codes [U].
- GPL-free fallback: learn the codes from captures (the paper states the SOP code is found among ~7 candidates and PN codes can be extracted without the ID; search snippet of GRDSMX paper) [U].

### B2. DSM2 / DSMX protocol (reverse engineered; GPL sources, facts only)

Manufacturer ID: 6 bytes read from CYRF MFG_ID_ADR 0x25 (CY-DS p14); the protocol uses id0..id3 (DEV-DSM, MPM-DSM).

Radio settings: data packets in 64-chip **8DR** (TX_CFG 0x28|PA), 64-chip SOP, preamble 0x333304, framing 0xEA, length 16; bind in 64-chip **SDR** (0x38|PA) with hardware CRC disabled (RX_OVERRIDE 0x14, TX_OVERRIDE 0x04), framing 0x4A (DEV-DSM; MPM-DSM comment "64-SDR").

Air time per data packet (computed): 64 µs preamble + 128 µs SOP + 19 bytes × 64 µs (length + 16 + CRC 2) ≈ 1.41 ms.

Bind packet (RF channel **13**, = 2413 MHz if f = 2400 + ch [U]; 16 bytes; MPM-DSM):
- bytes 0-3 = ~id0..~id3; bytes 4-7 = repeat of 0-3; bytes 8-9 = checksum, big-endian: 384 − 0x10 + Σ bytes 0-7; byte 10 = 0x01; byte 11 = number of channels (≤ 12); byte 12 = protocol/rate code; byte 13 = 0x00; bytes 14-15 = running checksum continued over bytes 8-13.
- Byte 12: **0x01** DSM2 1024/22 ms (< 8 ch), **0x02** DSM2 1024/22 ms (≥ 8 ch), **0x12** DSM2 2048/11 ms, **0xA2** DSMX 22 ms, **0xB2** DSMX 11 ms (MPM-DSM; DEV-DSM). (0x23 = DSM2 "SFC" variant, MPM only.) MPM notes the 11/22 ms choice is not fully encoded in this byte for DSMX 1F/2F [U].
- Bind sent every ~10 ms (MPM) / 8.5 ms (DEV) for ~1.8 s or more; bind data code = pncodes[0][8] for both halves (DEV-DSM) [GPL-only].

Data packet (16 bytes; MPM-DSM, DEV-DSM):
- bytes 0-1: **DSMX: id2, id3; DSM2: ~id2, ~id3**.
- bytes 2-15: 7 channel words, big-endian 16-bit. 11-bit format (2048; DSMX and DSM2 11 ms): bits 10:0 value, bits 14:11 channel number, bit 15 set only on the first word of the second ("upper") frame of a > 7-channel set. 10-bit format (DSM2 1024/22 ms): bits 9:0 value, channel number from bit 10. Unused words 0xFFFF. Default throws 1100-1900 µs.
- GRDSMX parses words 1-7 as bit 15 flag, bits 14:11 channel, bits 10:0 value.

Channel selection:
- **DSM2**: 2 RF channels. MPM picks them by scanning for quiet channels (`CYRF_FindBestChannels`, then 2 random distinct choices); Deviation derives them from the mfg ID. For a genuine Spektrum DSM2 TX they are not assumed derivable from the ID [U].
- **DSMX**: 23 channels from the mfg ID (identical in DEV-DSM and GRDSMX; originally by S. Gimaev / A. Alexandrov, Deviation forum https://www.deviationtx.com/forum/protocol-development/1160-dsmx-channels-list-for-frequency-shift-hopping):
  - seed s0 = ~(id0<<24 | id1<<16 | id2<<8 | id3) (32-bit, bitwise inverted); s = s0;
  - repeat: s = s × 0x0019660D + 0x3C6EF35F (mod 2^32); candidate ch = ((s >> 8) mod 0x49) + 3 → 3..75;
  - reject if (ch & 1) ≠ (s0 & 1) (parity tied to the seed's LSB), if already in the list, or if its band is full: 3-27 max 8, 28-51 max 7, 52-75 max 8; stop at 23 channels.
  - Spektrum confirms 23 channels per TX, a pattern calculated from the TX GUID (SPM).

PN code selection (DEV-DSM, MPM-DSM):
- sop_col = (id0 + id1 + id2 + 2) & 7 (MPM flips an ID bit if the result is 0: its own quirk); data_col = 7 − sop_col.
- row = **channel mod 5 (DSM2)**, **(channel − 2) mod 5 (DSMX)**.
- SOP code = pncodes[row][sop_col]; data code = pncodes[row][data_col] ‖ pncodes[row][data_col+1] [GPL-only table].

CRC seed: base = (id0<<8) | id1; the seed alternates every packet between base and ~base (DEV: seed starts as ~base and toggles each `set_sop_data_crc`; MPM: init base, toggled per packet). Which of the CH1/CH2 packets gets the inverted seed differs between the two descriptions [U]; a receiver should test both.

Timing (DEV-DSM, MPM-DSM):
- Two packets per frame on consecutive hop channels (CH1/CH2), **4010 µs apart**; frame period **11 ms** or **22 ms** (22 ms when < 8 channels / DSM2 1024). DSMX advances the channel index mod 23, DSM2 mod 2. Telemetry RX windows (MPM: 1950 µs / 600 µs delays) fill the rest of the frame.
- Telling DSM2 from DSMX on air: bytes 0-1 inverted (DSM2) or not; row formula; distinct channels (2 vs 23); frame rate from the timing; bind byte 12.

### B3. Walkera DEVO (same CYRF6936; DEV-DEVO, MPM-CYRF)
- 32-chip **8DR** (TX_CFG 0x08|PA → 250 kbps), 32-chip SOP, SOP threshold 0x04, 16-byte packets, one packet every **2.4 ms**. Air time (computed): ≈ 19 bytes × 32 µs + SOP 64 µs + preamble ≈ 0.7 ms.
- 3 RF channels found by scanning (`CYRF_FindBestChannels(…, 3, 4, 4, 80)`), **4 packets per channel** (~9.6 ms dwell).
- Byte 0: high nibble = number of channels, low nibble = type: 0x0B-0x0D data (4 channels each, 2 bytes per channel, signs in byte 9), 0x07/0x08 beacon/failsafe, 0x0A bind; telemetry replies start 0x30-0x36. Bytes 10-15: byte 10 = state (0x80 fixed-ID, 0xC0 binding) + countdown 3..0 within the 4-packet group; bytes 11-12 = next two hop channels; bytes 13-15 = 3-byte "fixed ID" (model ID, ≤ 999 999).
- Scrambling: bytes 1-15 XORed with the mfg ID bytes in a repeating 4-byte pattern (not in bind packets).
- SOP: one of a **10-row table** (Cypress TRM recommended codes per the MPM comment), index = low byte of (id0·4 + id1 + id2) mod 10; CRC seed = id0 + (id1 >> 6) + id2 (0 → 1), the same value in both seed bytes [U: from a prose reading of GPL code]. Before bind: default SOP, seed 0.
- Bind: type 0x0A packets alternating with data; mfg ID bytes 0-3 in clear in bytes 6-9; fixed ID XOR id0..id2 in bytes 13-15.

### B4. What a passive receiver can extract
- **DSMX, from one data packet**: id2, id3 in clear (bytes 0-1); id0, id1 from the CRC seed (the CRC is affine in the seed: solve, or try 2^16 with both polarities); (id0+id1+id2+2) & 7 cross-checks the SOP column. With id0..id3 the full 23-channel hop list follows and can be checked against observed channels → **32-bit unit identity** (the mfg-ID part DSM uses), plus frame rate and channel count. The GRDSMX paper reports identifying "the emitter by its manufacturing ID" (search snippet of https://pubs.gnuradio.org/index.php/grcon/article/view/97).
- **DSM2**: ~id2, ~id3 in bytes 0-1; id0/id1 from the CRC seed; hop channels not ID-derived.
- **Bind packets**: the whole id0..id3 (inverted), plus protocol byte (DSM2/DSMX, rate) and channel count.
- **DEVO**: bind packets give the mfg ID + fixed ID; data packets are XOR-scrambled with the mfg ID, but known-plaintext fields (countdown in byte 10, next-hop channels in bytes 11-12, which the receiver observes) reveal mfg-ID bytes [U: inference].
- Demodulation needed: GFSK discriminator over a 1 MHz channel (deviation ~±270-323 kHz) at ≥ 2-4 MS/s, 1 Mchip/s clock recovery, 64-chip SOP correlation (32-chip for DEVO), then 8DR despreading (4 base codes × 64 rotations per byte; 32-chip analogue for DEVO). Without the PN table, SOP/data codes must be learned from captures [U].
- Prior work: J. Andersson "Icarus" (2016) hijacked DSMx by extracting/brute-forcing the pairing secret (https://www.pcworld.com/article/410844/no-need-to-shoot-down-drones-many-of-them-can-now-be-hijacked.html).

---

## Implementation notes for a detector

**Crossfire 150 Hz (FSK)**
1. Channelise 860-886 / 902-928 MHz (863-869 CE) on the 260 kHz (112 kHz CE) grid: f = Fmin + n·Δ (table A2).
2. Detect ~2.9 ms (uplink) and ~2.0 ms (downlink) 2-FSK bursts, ±42.5 kHz, 85.1 kbit/s (~170 kHz occupied), on a 6.667 ms cadence, a new channel each exchange; in normal mode the downlink is 13.0 MHz above the uplink, ≤ 2.6 ms after it.
3. Demodulate (FM discriminator, h ≈ 1), find the 0xAA/0x55 preamble + sync **0x3E 0x6B 0x3E** [U: unit independence], de-whiten (SX127x PN9 x^9+x^5+1, seed 0x1FF — SX1276-DS), read 23 / 13 bytes.
4. Validate: byte 0 low nibble = 3; CRC-16 (0x1021, reflected, xorout 0) over bytes 0-20 vs bytes 21-22 LE with an unknown init → solve the init (linear) and track it per hop slot; CRC-8 (0x07) over bytes 0-11 for the downlink.
5. Identity: per-slot CRC-init vector, learned 150-slot hop sequence, carrier ppm, band/mode, 8/12-ch alternation. Label "TBS Crossfire-like" until the CRC-init model is verified on the bench.
6. 50 Hz LoRa: parameters unknown — detect LoRa chirps on the same grid and estimate SF/BW blind; a bench capture is needed.

**Spektrum DSMX / DSM2**
1. Channelise 2400-2480 MHz in 1 MHz steps (DSMX ch 3-75, bind ch 13). Bursts ≈ 1.41 ms, pairs 4.01 ms apart, every 11 or 22 ms; a new channel each burst (DSMX) or 2 alternating channels (DSM2).
2. GFSK demod (~±270-323 kHz), slice at 1 Mchip/s.
3. Detect preamble 0x3333 (4 × 16 µs), then the SOP: correlate against the 64-chip SOP candidates (needs the PN table [GPL-only] or codes learned from captures).
4. Despread 64-chip 8DR: each 64-chip block → byte = 64·index{¬A, A, ¬B, B} + rotation (≥ 54/64 matches). First byte = length 0x10, then 16 bytes, then CRC16 (0x8005, seed from the ID, both polarities).
5. Classify: DSMX (bytes 0-1 = id2, id3; row (ch−2) mod 5) vs DSM2 (inverted; row ch mod 5); frame rate from timing; 11-bit vs 10-bit channel words.
6. Identity: id0..id3 (B4) → compute the 23-channel list (LCG ×0x0019660D +0x3C6EF35F, (s>>8) mod 0x49 + 3, parity, band caps 8/7/8) and confirm it on later bursts.
7. Constants: SOP default 0x17FF9E213690C782 (vendor), preamble 0x333304, PN table 5×9×8 bytes [GPL-only], LCG constants, CRC poly 0x8005, bind channel 13, CH1→CH2 4010 µs, frame 11/22 ms, bind byte 12 values.

**Walkera DEVO**: same front end; 32-chip 8DR (250 kbps), 32-chip SOP from the 10-row Cypress TRM table (vendor via the TRM, not found online; GPL copy in MPM-CYRF), 16-byte packets every 2.4 ms, 3 channels × 4 packets; descramble with the mfg-ID XOR.

**Licence status**: Crossfire constants come from g3gg0's MIT sniffer and blog (facts; usable). CYRF6936 PHY constants come from the Cypress datasheet (vendor; usable). DSM/DEVO protocol logic (packet layouts, LCG, row/column rules, CRC seeds) is published in GPL-3 code (Deviation, Multiprotocol, gr-dsmx-rc); the algorithms and numeric constants are facts that can be reimplemented from this description, but the **DSM 5×9×8 PN table is GPL-only** as found, and the DEVO 10-row SOP table is attributed to the Cypress TRM (vendor), which was not located.
