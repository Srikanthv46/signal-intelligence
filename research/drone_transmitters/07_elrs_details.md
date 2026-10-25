# 07 — ExpressLRS on-air packet details (SX1280 FLRC/LoRa, SX127x/LR1121 LoRa)

Research date 2026-10-21. Web sources only. `[U]` = uncertain / not confirmed by a primary source.
Purpose: reference for a passive RF-analysis library that detects and labels drone RC links from IQ.
ExpressLRS (ELRS) firmware is GPL-3; the format below is described in prose with exact constants, no code copied.
The fetch tool strips line numbers, so citations give **file + function**, not line numbers.
"master" = ExpressLRS master as fetched 2026-10-21 (4.x line).

Source key:
- [S1] Semtech SX1280/SX1281 datasheet Rev 3.2 (March 2020, DS.SX1280-1.W.APP). https://cdn.components.ru/media/PDF/1531189.pdf
- [S2] src/lib/SX1280Driver/SX1280.cpp
- [S3] src/lib/SX1280Driver/SX1280_Regs.h
- [S4] src/src/common.cpp
- [S5] src/include/common.h
- [S6] src/lib/OTA/OTA.cpp
- [S7] src/lib/OTA/OTA.h
- [S8] src/src/tx_main.cpp
- [S9] src/src/rx_main.cpp
- [S10] src/lib/FHSS/FHSS.cpp
- [S11] src/lib/FHSS/FHSS.h
- [S12] src/lib/FHSS/random.cpp
- [S13] src/lib/CRC/crc.cpp and crc.h
- [S14] src/lib/SX127xDriver/SX127x.cpp + SX127xRegs.h
- [S15] src/lib/LR1121Driver/LR1121.cpp
- [S16] PR #2933 "Update OTA_VERSION_ID to V4" (merged 2024-09-23): github.com/ExpressLRS/ExpressLRS/pull/2933/files
- [S17] 3.5.3 tag equivalents of S4/S8/S7/S14
- [S18] 4.0 release notes mirrors: forum.wearefpv.fr/topic/25870-expresslrs-400 ; intofpv.com/archive/index.php/thread-29408.html
(All [S2]..[S15],[S17] are files under github.com/ExpressLRS/ExpressLRS/blob/<ref>/ )

---

## 1. SX1280/SX1281 FLRC on air (datasheet [S1])

**Modem.** FLRC = "Fast Long Range Communication": coherent demodulation of **GMSK** plus FEC and
interleaving ([S1] §6.3 p.38, §7.3 p.47). GMSK implies **modulation index h = 0.5**; the datasheet names
GMSK and gives no separate h. One bit = one symbol, so Rb is both the bit rate and the symbol rate.

**Raw bit rates (Table 6-4 p.38 / Table 14-31 p.121).** The column Rb is the **raw on-air symbol/bit rate**;
FEC is applied on top, giving the lower effective rate Rb_eff (Table 6-5 p.39). So e.g. FLRC_BR_0_650 = 0.65 Mb/s raw.

| enum (modParam1 value) | Rb raw Mb/s | BW DSB MHz | Rb_eff @CR 1 / 3-4 / 1-2 (Mb/s) |
|---|---|---|---|
| FLRC_BR_1_300_BW_1_2 (0x45) | 1.3 | 1.2 | 1.3 / 0.975 / 0.65 |
| FLRC_BR_1_000_BW_1_2 (0x69) | 1.04 | 1.2 | 1.04 / 0.78 / 0.52 |
| FLRC_BR_0_650_BW_0_6 (0x86) | 0.65 | 0.6 | 0.65 / 0.488 / 0.325 |
| FLRC_BR_0_520_BW_0_6 (0xAA) | 0.52 | 0.6 | 0.52 / 0.39 / 0.26 |
| FLRC_BR_0_325_BW_0_3 (0xC7) | 0.325 | 0.3 | 0.325 / 0.244 / 0.163 |
| FLRC_BR_0_260_BW_0_3 (0xEB) | 0.26 | 0.3 | 0.26 / 0.195 / 0.130 |

BW is the double-sideband modulation bandwidth. RX frequency tolerance (Table 6-6 p.40): +/-150 kHz for the
0.52/0.65/1.04/1.3 rates, +/-75 kHz for 0.26/0.325. Enum byte values confirmed in both [S1] Table 14-31 and
[S3].

**On-air packet structure (Fig 7-6/7-7, p.47-48).** In transmit order:
1. **AGC preamble** — variable, 1 to 4 bytes. Can be 1 byte at 1.3 Mb/s; >=2 bytes for all other rates ([S1] §7.3.2).
2. **21-bit timing-recovery preamble** — fixed, always present. Related to a 21-bit Barker-like pattern used for
   start-of-packet detection ([S1] §7.3.2, errata 16.4 p.151). This is the "fixed internal preamble" that exists
   before any user sync word. Not independently configurable.
3. **Sync word** — 32 bits, of which 31 bits may be configured (1 bit effectively fixed by the modem) ([S1]
   Fig 7-6; §14.3.1 Table 14-35). Length options: `FLRC_SYNC_NOSYNC` (0x00) = no sync word (21-bit preamble only),
   or `FLRC_SYNC_WORD_LEN_P32S` (0x04) = 21-bit preamble + 32-bit sync word. Up to 3 sync words can be matched
   at once (3 correlators); the RX match combination is set by packetParam3 (RX_MATCH_SYNC_WORD_1..1_2_3).
4. **Header** — only in variable-length format (Fig 7-7): 16 bits = 2-bit type field + 7-bit RFU + 7-bit length.
   Fixed-length format has no header ([S1] §7.3.2/§7.3.3). ELRS uses fixed length, so no header on air.
5. **Payload** — 6 to 127 bytes ([S1] Table 14-38; minimum payload 6 bytes).
6. **CRC** — optional, 2 / 3 / 4 bytes ([S1] Table 14-39).
7. **Tail** — fixed 6-bit sequence of trailing zeros (convolutional-code flush) ([S1] §7.3.2).

**Bit order / coding.** Sync word is written MSB-first into registers (see §2). The sync word, header and CRC/
payload are transmitted after **convolutional encoding**; the CRC checksum covers all preceding bits **except the
AGC preamble** ([S1] §7.3.2, "CRC is performed on all of the preceding bits except for the preamble"). The
convolutional encoding span (Fig 7-6/7-7) covers header+payload+CRC+tail — i.e. the sync word is NOT FEC-coded,
it is sent raw and used for correlation; FEC applies from the header/payload onward.

**FEC (coding rate, Table 14-32 p.122).** Convolutional code. Options: `FLRC_CR_1_2` (0x00, rate 1/2),
`FLRC_CR_3_4` (0x02, rate 3/4), `FLRC_CR_1_1` (0x04, rate 1 = no coding). The tail is 6 bits when CR = 1/2 or 3/4,
and 0 bits when CR = 1 ([S1] §7.3.4, n_tail). **The convolutional generator polynomials and constraint length are
NOT published** in the datasheet. [U] The 6-bit flush tail implies constraint length K = 7 (6 memory elements),
consistent with a common K=7 rate-1/2 industry code, but Semtech does not state the polynomials. [U] Any
interleaver pattern for FLRC is also unpublished.

**BT / Gaussian filter (Table 14-33 p.122).** Set by modParam3 (TX pulse shaping only): `BT_DIS` (0x00, no filter),
`BT_1` (0x10, BT=1.0), `BT_0_5` (0x20, BT=0.5).

**Preamble length (AGC) (Table 14-34 p.122).** packetParam1, value = PREAMBLE_LENGTH_n_BITS: 0x00 reserved,
0x10=8, 0x20=12, 0x30=16, 0x40=20, 0x50=24, 0x60=28, 0x70=32 bits. Minimum 8 bits at 1 Mb/s, 16 bits otherwise.

**Whitening (Table 14-41 p.125).** Whitening **cannot be enabled** in FLRC. packetParam7 must always be
`WHITENING = 0x08` (= "whitening disabled"). So FLRC payloads are **not whitened**.

**CRC (Tables 14-39/14-40 p.124-125).** Length 0/2/3/4 bytes. The CRC polynomial is **fixed per CRC length**:
- CRC-16: 0x755B, P16(x)=x16+x14+x13+x12+x10+x8+x6+x4+x3+x+1
- CRC-24: 0x5D6DCB (datasheet), P24 as listed
- CRC-32: 0xF4ACFB13, P32 as listed

The **CRC init/seed** is loaded via register writes ("CrcInit"), MSB at 0x9C8, LSB at 0x9C9 ([S1] Table 14-40).
These are the same `0x9C8/0x9C9` registers the datasheet also labels "CRC MSB/LSB Initial Value — CRC Seed used for
GFSK and FLRC" ([S1] Register Map p.102). [U] The datasheet does not state whether FLRC CRC is reflected or has a
final XOR; the ELRS code never reads back FLRC CRC in software (hardware checks it), so the reflection/xor-out
convention is not determinable from ELRS.

**Sync word registers (Table 14-42 p.125-126).** SyncWord1 = 0x09CF..0x09D2 (MSB..LSB), SyncWord2 = 0x09D4..0x09D7,
SyncWord3 = 0x09D9..0x09DC. Sync-word bit-error tolerance is set at register 0x9CD ("Synch Address Control",
0..15 bit errors; reset 0x80) ([S1] p.123, p.103).

**Errata 16.4 (p.151) — FLRC elevated PER with sync word.** If the 32-bit sync word too closely resembles the
21-bit Barker preamble pattern, PER rises. For CR 1/2 and CR 3/4 the first 16 bits of the sync word must avoid
`0x8C38xxxx` and `0x630Exxxx`. For CR 3/4 the two LSBytes must additionally not be in range 0x0000..0x3EFF.
(This is why ELRS applies a byte-swap workaround, see §2.)

---

## 2. ExpressLRS FLRC mode (2.4 GHz SX1280)

**Air-rate table (common.cpp `ExpressLRS_AirRateConfig` SX1280 block [S4]; struct fields from [S5]).** ELRS has
four FLRC entries, all sharing: bw = `FLRC_BR_0_650_BW_0_6`, sf-field = `FLRC_BT_1`, cr = `FLRC_CR_1_2`,
PreambleLen = 32, TLMinterval = TLM_RATIO_1_128, FHSShopInterval = 2, PayloadLength = `OTA4_PACKET_SIZE` (8).

| index | enum_rate | interval (us) | numOfSends |
|---|---|---|---|
| 0 | RATE_FLRC_2G4_1000HZ | 1000 | 1 |
| 1 | RATE_FLRC_2G4_500HZ | 2000 | 1 |
| 2 | RATE_FLRC_2G4_500HZ_DVDA | 1000 | 2 |
| 3 | RATE_FLRC_2G4_250HZ_DVDA | 1000 | 4 |

So FLRC always runs: **raw 0.65 Mb/s GMSK, BW 0.6 MHz DSB, CR 1/2 (effective 0.325 Mb/s), BT=1, 32-bit AGC
preamble, fixed-length (no header), 8-byte payload, CRC on.** DVDA = "Deja Vu Diversity Aggregation": same packet
sent 2 or 4 times to make a lower effective loop rate; on-air radio config is identical to 1000 Hz.

**Radio config applied (SX1280.cpp `Config` -> `ConfigModParamsFLRC` + `SetPacketParamsFLRC` [S2]):**
- ModParams = {bw, cr, bt} = {FLRC_BR_0_650_BW_0_6, FLRC_CR_1_2, FLRC_BT_1}.
- Packet params: HeaderType = `FLRC_PACKET_FIXED_LENGTH` (0x00); preamble encoded as `((PreambleLen/4)-1)<<4`
  (<8 raised to 8); SyncWordLen = `FLRC_SYNC_WORD_LEN_P32S` (0x04, 32-bit); match = `RX_MATCH_SYNC_WORD_1`
  (0x10, SyncWord1 only); PayloadLength = 8; CRC = `FLRC_CRC_3_BYTE` (0x20, **3-byte / 24-bit CRC**);
  whitening byte = 0x08 (disabled, forced).
- **CRC seed:** ELRS writes the 16-bit `flrcCrcSeed` (= `OtaCrcInitializer`, see §3) to register
  `REG_FLRC_CRC_SEED` = 0x9C8 (MSB), 0x9C9 implied (LSB), high byte then low byte. So ELRS seeds the FLRC CRC with
  the **same 16-bit UID-derived initializer** it uses for the software LoRa CRC. [S2] `SetPacketParamsFLRC`, [S3].

**32-bit FLRC sync word from UID (SX1280.cpp `SetPacketParamsFLRC`, tx_main.cpp `SetRFLinkRate` [S2][S8]).**
The sync word passed to the driver is `OtaGetUidSeed()`:
```
syncWord32 = (UID[2]<<24) | (UID[3]<<16) | (UID[4]<<8) | (UID[5] ^ OTA_VERSION_ID)
```
([S6] `OtaGetUidSeed`; master). It is written **MSB-first** to 0x09CF..0x09D2:
`buf[0]=sync>>24 (UID2)`, `buf[1]=sync>>16 (UID3)`, `buf[2]=sync>>8 (UID4)`, `buf[3]=sync (UID5 ^ OTA_VERSION_ID)`.
So on air the FLRC SyncWord1 = **UID bytes 2,3,4 then (UID5 XOR OTA_VERSION_ID)**.

**Errata workaround (SX1280.cpp `SetPacketParamsFLRC` [S2]).** Because CR = 1/2: if `buf[0],buf[1]` == `0x8C,0x38`
or `0x63,0x0E`, ELRS swaps buf[0]<->buf[1] before writing (per errata 16.4). For CR 3/4 it would also OR 0x80 into
buf[3] if buf[3] <= 0x3E, but ELRS FLRC always uses CR 1/2 so only the swap applies. A passive receiver must apply
the same conditional swap when reconstructing UID[2],UID[3] from an observed sync word.

**FLRC CRC seed value.** `OtaCrcInitializer` is computed from the UID (see §3) and, in binding, forced to
`OTA_VERSION_ID` (master) / `0` (3.5.3). Same value drives both FLRC-hardware-CRC seed and the OTA packet CRC.

---

## 3. ExpressLRS OTA packets and CRC (LoRa and FLRC share these)

Current **OTA_VERSION_ID = 4** (common.h; raised from 3 in PR #2933, merged 2024-09-23 [S16]). It must be <=15.
3.x used 3. 4.x is not cross-compatible with 1.x-3.x [S18].

**IQ inversion rule (SX1280 LoRa, and whole link) — tx_main/rx_main `SetRFLinkRate` [S8][S9]:**
```
invertIQ = InBindingMode || (UID[5] & 0x01)
```
So: IQ is inverted when the **LSB of UID[5] = 1**, and **always inverted while binding**. This is a single
recoverable identity bit for SX1280 LoRa (FLRC GMSK has no IQ-polarity concept; SX127x sub-GHz does not apply IQ
inversion in ELRS — see §note). Binding on 2.4 GHz LoRa therefore always uses inverted IQ.

**SX1280 LoRa sync word.** ELRS does **not** set a LoRa sync word on SX1280 ([S2] `Config` LoRa branch writes no
sync-word register; the SX1280 modem actually has no LoRa sync word — RX never gets SYNCWORD_VALID on LoRa). The
"sync" on SX1280 LoRa is only the preamble + implicit-header framing + CRC. [U] (SX1280 LoRa has a configurable
1-byte sync via regs 0x944/0x945 per [S1] Table 13-1, but ELRS leaves it default.)

**SX127x (900 MHz) LoRa sync word.** ELRS leaves the default LoRa sync word `SX127X_SYNC_WORD = 0x12` (18)
(SX127x.cpp `ConfigLoraDefaults`->`SetSyncWord`, SX127xRegs.h [S14]). The `Config` `syncWord`/`InvertIQ`
parameters are accepted but **not applied** on SX127x — ELRS does NOT invert IQ on sub-GHz SX127x (verified
master and 3.5.3: no write to `SX127X_REG_INVERT_IQ` = 0x33, mask 0x41). Implicit header, CRC off in software
(computed in payload), preamble LSB default 8.

**OTA4 layout — 8 bytes (`OTA4_PACKET_SIZE`=8, OTA_Packet4_s [S7]):**
- Byte 0: `type` (bits 1:0) + `crcHigh` (bits 7:2, the high 6 bits of a 14-bit CRC).
- Bytes 1..6: union — RC data, or the **SYNC** struct `OTA_Sync_s` (6 bytes), or uplink/downlink data.
- Byte 7: `crcLow` (low 8 bits of the CRC).
- `OTA4_CRC_CALC_LEN` = offsetof(crcLow) = **7** bytes covered by the CRC.

**OTA8 layout — 13 bytes (`OTA8_PACKET_SIZE`=13, OTA_Packet8_s [S7]):**
- Byte 0: `packetType` in low 2 bits (+ rc fields / per-variant).
- SYNC variant: byte0 packetType, bytes 1..6 `OTA_Sync_s`, bytes 7..10 free[4], bytes 11..12 `uint16_t crc` (LE).
- `OTA8_CRC_CALC_LEN` = offsetof(crc) = **11** bytes covered by the CRC.

**Packet type field (OTA.h [S7]).** Low 2 bits of byte 0: `PACKET_TYPE_RCDATA`=0b00 (uplink),
`PACKET_TYPE_DATA`=0b01, `PACKET_TYPE_SYNC`=0b10, `PACKET_TYPE_LINKSTATS`=0b00 (downlink). `MODELMATCH_MASK`=0x3F.

**OTA_Sync_s — SYNC packet fields (6 bytes).** NOTE the layout differs between versions:
- **master / 4.x [S7]:** byte0 `fhssIndex`; byte1 `nonce`; byte2 `rfRateEnum` (full 8-bit rate enum);
  byte3 bitfield `switchEncMode:1, newTlmRatio:3, geminiMode:1, otaProtocol:2, free:1`; byte4 `UID4`; byte5 `UID5`.
- **3.5.3 [S17]:** byte0 `fhssIndex`; byte1 `nonce`; byte2 bitfield `switchEncMode:1, newTlmRatio:3, rateIndex:4`;
  byte3 `UID3`; byte4 `UID4`; byte5 `UID5`.
So 3.5.3 SYNC leaks **UID3, UID4, UID5** (3 bytes); 4.x SYNC leaks only **UID4, UID5** (2 bytes) and moved the
rate to a full byte + added geminiMode/otaProtocol. The field that carries identity in SYNC: UID bytes 4,5
(and 3 on 3.5.3).

**SYNC field semantics (tx_main `GenerateSyncPacketData` [S8]):**
- `fhssIndex` = current FHSS sequence pointer `FHSSgetCurrIndex()`.
- `nonce` = `OtaNonce` (8-bit packet counter, wraps).
- rate = current air-rate index/enum.
- `switchEncMode` = configured switch mode; `newTlmRatio` = effective TLM ratio **minus** `TLM_RATIO_NO_TLM`
  (stored as offset); 4.x adds `geminiMode` and `otaProtocol` (= link mode: CRSF/MAVLink).
- `UID4`/`UID5` = UID[4], UID[5], except when model-match is on and not binding: `UID5 ^= (~modelId) & 0x3F`.
  (3.5.3 also sends `UID3 = UID[3]`.)

**CRC algorithm (CRC/crc.cpp class `Crc2Byte` [S13]).** Generic table-driven CRC, **MSB-first, non-reflected,
no final XOR**. Table seed per byte i: `i << (bits-8)`, 8 iterations of "test highbit (bit bits-1); shift left 1;
if set XOR poly". calc: per data byte `idx = ((crc >> (bits-8)) ^ byte) & 0xFF; crc = (crc<<8) ^ table[idx]`,
final `crc & bitmask` (bitmask = (1<<bits)-1). Init value used as supplied (not masked by generator).

**CRC polynomials (OTA.h [S7]).** OTA4 = **14-bit CRC**, `ELRS_CRC14_POLY = 0x2E57` (Koopman 0x372B).
OTA8 = **16-bit CRC**, `ELRS_CRC16_POLY = 0x3D65` (Koopman 0x9EB2).

**CRC init / seed (OTA.cpp [S6]).**
```
OtaCrcInitializer = ((UID[4] << 8) | UID[5]) ^ ((uint16_t)OTA_VERSION_ID << 8)
```
i.e. high byte = `UID[4] ^ OTA_VERSION_ID`, low byte = `UID[5]`. Per packet the actual init is
`OtaCrcInitializer ^ nonceValidator`, where `nonceValidator = (type==PACKET_TYPE_SYNC) ? 0 : OtaNonce`.
So for SYNC packets the CRC init is just `OtaCrcInitializer`; for all other packets the current `OtaNonce` XORs
into the low byte. During binding, `OtaCrcInitializer = OTA_VERSION_ID` (master) / `0` (3.5.3) [S8][S17][S18].

**Bytes the CRC covers & bit placement.**
- OTA4: CRC computed over bytes [0 .. 6] (7 bytes; byte0 includes the type bits and the crcHigh field which the
  **validator zeroes before recomputing**; generator assumes crcHigh=0 when computing). The 14-bit result is split
  `crcHigh = crc>>8` (into byte0 bits 7:2) and `crcLow = crc & 0xFF` (byte 7).
- OTA8: CRC-16 over bytes [0 .. 10] (11 bytes); result stored as `uint16_t crc` little-endian in bytes 11-12.

---

## 4. ExpressLRS FHSS (FHSS.cpp / FHSS.h [S10][S11])

**Domain tables (`fhss_config_t`: start Hz, stop Hz, channel count; sync_channel = count/2).** Channel i centre:
`freq_start + freq_spread * i / FREQ_SPREAD_SCALE`, where
`freq_spread = (stop-start)*FREQ_SPREAD_SCALE/(count-1)` and `FREQ_SPREAD_SCALE = 256` (SX127x/SX128x; = 1 for
LR1121/LR2021/unit-test) [S11]. So channel spacing = (stop-start)/(count-1).

| Domain | start MHz | stop MHz | count | spacing MHz | sync_channel (count/2) |
|---|---|---|---|---|---|
| ISM2G4 (SX128x; "CE_LBT" if EU_CE_2400) | 2400.4 | 2479.4 | 80 | 1.0 | 40 |
| FCC915 | 903.5 | 926.9 | 40 | 0.6 | 20 |
| EU868 | 863.275 | 869.575 | 13 | 0.525 | 6 |
| IN866 | 865.375 | 866.95 | 4 | 0.525 | 2 |
| AU915 | 915.5 | 926.9 | 20 | 0.6 | 10 |
| AU433 | 433.42 | 434.42 | 3 | 0.5 | 1 |
| EU433 | 433.1 | 434.45 | 3 | 0.675 | 1 |
| US433 | 433.25 | 438.0 | 8 | ~0.679 | 4 |
| US433-wide | 423.5 | 438.0 | 20 | ~0.763 | 10 |
| TH920 | 920.5 | 924.7 | 8 | 0.6 | 4 |

(ISM2G4 2.4 GHz is the one used for FLRC and 2G4 LoRa. India = IN866 sub-GHz, or ISM2G4 at 2.4 GHz.)

**sync_channel rule.** `sync_channel = freq_count / 2` (integer division) [S10]. (Pre-PR#2933 it was
`freq_count/2 + 1`; changed to `freq_count/2` in the same PR that bumped OTA_VERSION_ID to 4 [S16]. So the sync
channel index for a given domain depends on firmware major version.)

**Sequence build (`FHSSrandomiseFHSSsequenceBuild` [S10]).** `FHSS_SEQUENCE_LEN = 256` [S11]. Seed =
`OtaGetUidSeed()` (= the same UID[2..5]^OTA_VERSION_ID value as the FLRC sync word) [S9][S6]. Each block of
`freqCount` entries starts with the sync channel; within a block entries are Fisher-Yates-shuffled using the LCG,
leaving the block-start sync entry fixed. `FHSSonSyncChannel()` = current sequence entry == sync_channel.
`FHSSgetInitialFreq()` returns the sync-channel frequency.

**RNG (random.cpp [S12]).** LCG: `seed = (214013*seed + 2531011) mod 2^31`; output = `seed>>16` (0..0x7FFF).
`rngSeed(x)` sets seed=x directly; `rngN(max) = rng() % max`. Because the high bit of the 32-bit seed is dropped
(mod 2^31) and only bits 16..30 are output, the top bit of the 4-byte seed does not affect the sequence — i.e.
seeds that differ only in bit 31 give identical hop sequences [S18].

**Hop interval per rate = `FHSShopInterval` (channels per hop; hop every that many packets).** From the air-rate
tables ([S4]): FLRC 2G4 (all 4 rates) = **2**. 2G4 LoRa: 500/333/250/150/100 Hz = **4**; 50 Hz = **2**.
900 LoRa: 200/100/100_8CH/50 Hz = **4**; 25 Hz = **2**; 50HZ_DVDA = **2**. A hop occurs when
`(OtaNonce+1) % FHSShopInterval == 0` [S8]. `NonceFHSSresult = OtaNonce % FHSShopInterval`.

---

## 5. ExpressLRS 2.4 GHz LoRa coding (long interleaving)

**Confirmed: 2.4 GHz LoRa uses SX1280 long-interleaving coding rates** (common.cpp SX1280 LoRa rows [S4]):
all 2G4 LoRa entries use `SX1280_LORA_CR_LI_4_6` (0x06) or `SX1280_LORA_CR_LI_4_8` (0x07); none use the legacy
CR_4_x. Per-rate: 500 Hz = BW0800/SF5/CR_LI_4_6; 333_8CH = SF5/CR_LI_4_8; 250 Hz = SF6/CR_LI_4_8; 150 Hz =
SF7/CR_LI_4_8; 100_8CH = SF7/CR_LI_4_8; 50 Hz = SF8/CR_LI_4_8. (BW always `LORA_BW_0800` = 0x18 = 812.5 kHz.)
CR enum values from [S3]: CR_LI_4_5=0x05, CR_LI_4_6=0x06, CR_LI_4_8=0x07.

Datasheet [S1] Table 14-49 (p.131) lists LORA_CR_LI_4_5/4_6/4_8 and notes "a new interleaving scheme ... to
increase robustness to burst interference/strong Doppler; complexity kept the same; coding rate signalled in
header; previously only 0x0-0x4 valid. Compatible with implicit header. Scrambling (whitening) occurs with long
interleaving as with legacy." Time-on-air formulas for LI are in [S1] §7.4.4.2 (p.51-52).

**Public documentation of the LI interleaver pattern: essentially none.** The datasheet gives ToA formulas but
**not the bit-level interleaver map**. gr-lora_sdr (EPFL, GPL-3) implements the **legacy** LoRa interleaver only;
its SX1280 issue #143 is unresolved and does not cover CR_LI. [U] No open-source or published reverse-engineered
description of the SX1280 long-interleaver was found. A passive decoder of 2G4 LoRa payload would have to
reverse-engineer LI from captures (preamble + implicit-header + whitening + LI de-interleave + Hamming decode).

SX1280 LoRa misc constants ([S1]/[S3]): BW_0800=0x18 (812.5 kHz), SF5=0x50..SF12=0xC0; IQ: INVERTED=0x00,
NORMAL=0x40 (note inverted==0 in the SX1280 encoding); header EXPLICIT=0x00, IMPLICIT/FIXED=0x80; CRC
ENABLE=0x20/DISABLE=0x00 (ELRS uses CRC_OFF and carries its own OTA CRC in the payload). SF-specific reg writes
to 0x925 (ELRS: 0x1E for SF5/6, 0x37 for SF7/8, 0x32 for SF9+). LoRa 2G4 preamble lengths (symbols) per rate:
12 (500/333/150/100), 14 (250), 12 (50).

---

## What a passive receiver can extract (per mode)

Goal: recover identity bits and/or label the link **without** the binding phrase. The binding phrase is MD5'd
and its first 6 bytes become UID[0..5]; UID[0],UID[1] are not sent over the air (bind MSP sends only UID[2..5],
and SYNC sends only UID[3..5] or UID[4..5]) [S18]. "Identity" recoverable from IQ = the UID bytes that appear
on air, plus the IQ-polarity bit.

- **2.4 GHz FLRC (0.65 Mb/s GMSK, BW 0.6 MHz, CR 1/2, 32-bit sync word):** demod = GMSK (h=0.5) coherent, find
  21-bit preamble, read the 32-bit sync word by correlation. The sync word directly = **UID[2], UID[3], UID[4],
  (UID[5]^OTA_VERSION_ID)** MSB-first (undo the errata byte-swap of bytes 0/1 if they match 0x8C38/0x630E). That
  gives **4 of the 6 UID bytes immediately, with no FEC/decode needed** (sync word is raw, pre-FEC). This is the
  strongest passive identity cue for ELRS. Hop interval 2, spacing 1 MHz on ISM2G4.

- **2.4 GHz LoRa (SF5-8, BW 812.5 kHz, CR_LI, implicit header):** demod = LoRa CSS (standard chirp detection).
  **IQ polarity = 1 bit of identity** (= UID[5]&0x01; inverted means that bit is 1, or the link is binding).
  Detecting/labelling the link needs only CSS preamble + chirp-slope/polarity — no payload decode. Decoding the
  payload (to read the SYNC packet) additionally needs de-whitening + the unpublished LI de-interleaver + Hamming
  decode + the 14/16-bit CRC with UID-derived init — not feasible without reversing LI [U].

- **900 MHz / sub-GHz LoRa (SX127x, SF6-9, BW 500 kHz, CR 4/7 or 4/8, implicit header, sync word 0x12, non-inverted
  IQ):** demod = standard LoRa CSS; decodable with ordinary LoRa receivers (gr-lora_sdr). The implicit-header
  payload is the raw OTA4 (8 B) or OTA8 (13 B) packet. A **SYNC packet** (type bits = 0b10 in byte0) yields, after
  CRC check: **fhssIndex, nonce, rate, switchEncMode, newTlmRatio, and UID[3],UID[4],UID[5]** (3.5.3) or
  **UID[4],UID[5]** (4.x). The CRC-14 init = (UID[4]<<8|UID[5]) ^ (OTA_VERSION_ID<<8), and for SYNC the nonce does
  not XOR in — so once UID[4],UID[5] are read from the SYNC body the CRC is self-consistent and confirms them.
  fhssIndex + nonce let a passive receiver align to the hop sequence.

- **UID bytes 0,1:** never transmitted; not recoverable passively. The remaining unknown for full UID is UID[2]
  (and UID[3] on 4.x LoRa), which on FLRC is in the sync word but on sub-GHz LoRa is not in SYNC.

- **Reconstructing the hop sequence:** seed = UID[2..5]^OTA_VERSION_ID via the known LCG (`214013/2531011 mod 2^31`,
  out>>16). With UID[2..5] known (e.g. from an FLRC sync word) the full 256-entry hop table and all channel
  frequencies are computable; because mod 2^31 drops seed bit 31, the top bit is irrelevant [S18].

---

## 10-line constants summary

1. OTA_VERSION_ID = 4 (master; was 3 in 3.x) — appears XORed into CRC init high byte and UID seed low byte.
2. FLRC on air: GMSK h=0.5, raw 0.65 Mb/s, BW 0.6 MHz DSB, CR 1/2, BT=1, 32-bit AGC preamble, 21-bit fixed Barker preamble, fixed-len (no header), 8-B payload, 3-byte(24-bit) CRC, whitening always OFF.
3. FLRC 32-bit sync word (regs 0x9CF..0x9D2, MSB-first) = UID[2],UID[3],UID[4],(UID[5]^OTA_VERSION_ID); errata swaps bytes0/1 if 0x8C38/0x630E; sent raw (not FEC-coded).
4. FLRC_BR enum bytes: 1.3=0x45,1.04=0x69,0.65=0x86,0.52=0xAA,0.325=0xC7,0.26=0xEB; CR 1/2=0x00,3/4=0x02,1=0x04; BT 1=0x10,0.5=0x20.
5. FLRC CRC polys fixed per length: 16b=0x755B, 24b=0x5D6DCB, 32b=0xF4ACFB13; CRC seed reg 0x9C8/0x9C9 (=OtaCrcInitializer). Conv-code polynomials NOT published [U].
6. OTA4 = 8 B (14-bit CRC poly 0x2E57, crcHigh in byte0 bits7:2, crcLow byte7, covers 7 B); OTA8 = 13 B (16-bit CRC poly 0x3D65, LE in bytes11-12, covers 11 B).
7. CRC init = (UID[4]<<8|UID[5]) ^ (OTA_VERSION_ID<<8); per-packet XOR nonce except SYNC; binding init = OTA_VERSION_ID (master)/0 (3.5.3). Crc2Byte = MSB-first, non-reflected, no xorout.
8. IQ inversion = InBindingMode || (UID[5]&0x01); SX1280 LoRa IQ_INVERTED=0x00/NORMAL=0x40; SX127x sub-GHz LoRa: sync word 0x12, IQ NOT inverted by ELRS.
9. 2G4 LoRa: BW 812.5 kHz, SF5-8, CR_LI_4_6/4_8 (long interleaving, enum 0x06/0x07) — LI interleaver undocumented publicly [U]. 900 LoRa: BW 500 kHz, SF6-9, CR 4/7-4/8.
10. FHSS: seed = UID[2..5]^OTA_VERSION_ID; LCG 214013/2531011 mod 2^31, out>>16; seq len 256; sync_channel = count/2; ISM2G4 2400.4-2479.4 MHz x80 @1 MHz; IN866 865.375-866.95 x4; FCC915 x40; hop interval 2 (FLRC) / 2-4 (LoRa).
