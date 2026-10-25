# 12. More RC links: mLRS, FrSky R9, CC2500 links (S-FHSS, Hitec, HoTT, RadioLink), OpenLRS-ng, AFHDS 3 / FASST / Tracer / Ghost, SX127x SF6 framing

Researched 2026-10-22 from the web only. It extends `03_rc_links.md`, `07_elrs_details.md` and
`10_cc2500_a7105_xn297_flarm_xbee.md` and does not repeat them. Scope: the on-air facts needed to **passively detect
and identify** a link. It does not cover keys, binding procedures, or transmitting.

**Tags:**
- **[spec]:** datasheet, standard or application note. Usable in the C++ library.
- **[vendor]:** manufacturer or retailer claim.
- **[code]:** open-source code. The licence is given for each source. **GPL code is a source of facts only**: re-implement
  from this description and never copy it.
- **[comm]:** forum post or third-party analysis.
- **[calc]:** computed here from cited register values with the cited datasheet formula.
- **[U]:** uncertain, inferred, or unpublished.

## Sources (cited by key)
| Key | URL | Licence |
|---|---|---|
| ML-128 | https://github.com/olliw42/mLRS/blob/main/mLRS/Common/sx-drivers/sx128x_driver.h | GPL-3 |
| ML-126 | https://github.com/olliw42/mLRS/blob/main/mLRS/Common/sx-drivers/sx126x_driver.h | GPL-3 |
| ML-127 | https://github.com/olliw42/mLRS/blob/main/mLRS/Common/sx-drivers/sx127x_driver.h | GPL-3 |
| ML-LR | https://github.com/olliw42/mLRS/blob/main/mLRS/Common/sx-drivers/lr11xx_driver.h | GPL-3 |
| ML-FT | https://github.com/olliw42/mLRS/blob/main/mLRS/Common/frame_types.h | GPL-3 |
| ML-FR | https://github.com/olliw42/mLRS/blob/main/mLRS/Common/frames.h | GPL-3 |
| ML-SET | https://github.com/olliw42/mLRS/blob/main/mLRS/Common/setup.h | GPL-3 |
| ML-CONF | https://github.com/olliw42/mLRS/blob/main/mLRS/Common/common_conf.h | GPL-3 |
| ML-FH | https://github.com/olliw42/mLRS/blob/main/mLRS/Common/fhss.h | GPL-3 |
| ML-FHC | https://github.com/olliw42/mLRS/blob/main/mLRS/Common/fhss.cpp | GPL-3 |
| ML-TOO | https://github.com/olliw42/mLRS/blob/main/mLRS/Common/theory_of_operation.h | GPL-3 |
| ML-RM | https://github.com/olliw42/mLRS/blob/main/README.md | GPL-3 |
| SXL | https://github.com/olliw42/sx12xx-lib/blob/main/src/sx128x.h, `.../src/lr11xx.cpp` | GPL-3 |
| FMAV | https://github.com/olliw42/fastmavlink/blob/master/c_library/lib/fastmavlink_crc.h | [U] |
| MPM-R9 | https://github.com/pascallanger/DIY-Multiprotocol-TX-Module/blob/master/Multiprotocol/FrSkyR9_sx1276.ino | GPL-3 |
| MPM-SX | https://github.com/pascallanger/DIY-Multiprotocol-TX-Module/blob/master/Multiprotocol/SX1276_SPI.ino | GPL-3 |
| MPM-C | https://github.com/pascallanger/DIY-Multiprotocol-TX-Module/blob/master/Multiprotocol/FrSkyDVX_common.ino | GPL-3 |
| MPM-FU | https://github.com/pascallanger/DIY-Multiprotocol-TX-Module/blob/master/Multiprotocol/Futaba_cc2500.ino | GPL-3 |
| MPM-HI | https://github.com/pascallanger/DIY-Multiprotocol-TX-Module/blob/master/Multiprotocol/Hitec_cc2500.ino | GPL-3 |
| MPM-HO | https://github.com/pascallanger/DIY-Multiprotocol-TX-Module/blob/master/Multiprotocol/HOTT_cc2500.ino | GPL-3 |
| MPM-RL | https://github.com/pascallanger/DIY-Multiprotocol-TX-Module/blob/master/Multiprotocol/RadioLink_cc2500.ino | GPL-3 |
| MPM-PD | https://github.com/pascallanger/DIY-Multiprotocol-TX-Module/blob/master/Protocols_Details.md | GPL-3 |
| OL-B | https://github.com/openLRSng/openLRSng/blob/master/openLRSng/binding.h | GPL-3 |
| OL-C | https://github.com/openLRSng/openLRSng/blob/master/openLRSng/common.h | GPL-3 |
| OL-RFM | https://github.com/openLRSng/openLRSng/blob/master/openLRSng/RFM.h | GPL-3 |
| OL-TX | https://github.com/openLRSng/openLRSng/blob/master/openLRSng/TX.h | GPL-3 |
| RH22 | https://github.com/hallard/RadioHead/blob/master/RH_RF22.h | GPL-2 / commercial |
| SX1276 | Semtech SX1276/77/78/79 datasheet Rev. 4 (March 2015): https://cdn-shop.adafruit.com/product-files/3179/sx1276_77_78_79.pdf | [spec] |
| SX1262 | Semtech SX1261/2 datasheet Rev. 1.2 (June 2019): https://cdn.sparkfun.com/assets/6/b/5/1/4/SX1262_datasheet.pdf | [spec] |
| CC2500 | TI CC2500 datasheet SWRS040C: https://www.ti.com/lit/ds/symlink/cc2500.pdf | [spec] |
| DN504 | TI DN504, FEC implementation (SWRA113A): https://www.ti.com/lit/an/swra113a/swra113a.pdf | [spec] |
| RL-LR | https://github.com/jgromes/RadioLib/blob/master/src/modules/LR11x0/LR11x0.cpp | MIT |
| RL-126 | https://github.com/jgromes/RadioLib/blob/master/src/modules/SX126x/SX126x.cpp | MIT |
| GRL | gr-lora_sdr: https://github.com/tapparelj/gr-lora_sdr (README), `lib/interleaver_impl.cc`, `lib/modulate_impl.cc` | GPL-3 |
| GRL62 | https://github.com/rpp0/gr-lora/issues/62 | GPL-3 project, forum |
| SEMF | Semtech developer forum: https://forum.lora-developers.semtech.com/t/sx1272-and-sx1262-lora-sync-word-compatibility/988 and https://forum.lora-developers.semtech.com/t/sx126x-127x-interoperability-at-500khz/1568/6 (seen in a search snippet only; the host did not resolve) | [comm] |
| RFENG | https://rfengfpv.wordpress.com/category/futaba/ ("Futaba S-FHSS Protocol Overview") | [comm] |
| EDGE | https://github.com/EdgeTX/edgetx/blob/main/radio/src/pulses/afhds3.h | GPL-2 |

---

## 1. mLRS (olliw42)

### 1.1 PHY per mode
| Mode / band | Chip | Modulation | Params | Preamble | Header / CRC / IQ | Frame period | Air time (code) | Ref |
|---|---|---|---|---|---|---|---|---|
| 50 Hz, 2.4 GHz | SX128x / LR1121 | LoRa | SF5, BW reg `LORA_BW_800` = **0x18**, CR **LI 4/5** | 12 sym | implicit, chip CRC off, **IQ normal** | 20 ms | 7.892 ms | [code] ML-128, ML-LR, SXL |
| 31 Hz, 2.4 GHz | same | LoRa | SF6, BW 0x18, LI 4/5 | 12 | same | 32 ms | 13.418 ms | same |
| 19 Hz, 2.4 GHz | same | LoRa | SF7, BW 0x18, LI 4/5 | 12 | same | 53 ms | 23.527 ms | same |
| FLRC 111 Hz, 2.4 GHz | SX128x | FLRC 650 kb/s, BW 0.6 MHz, CR 1/2, BT 1 | – | 32-bit AGC preamble | fixed length, chip CRC off, **32-bit sync word** (see 1.3) | 9 ms | 2.383 ms | [code] ML-128 |
| 31 Hz, 868 / 915 / 433 | SX126x / LLCC68 / LR1121 | LoRa | **SF5**, BW 500, CR 4/5 | 12 | implicit, CRC off, IQ normal | 32 ms | 13.2 ms | [code] ML-126, ML-LR |
| 19 Hz, 868 / 915 / 433 | SX126x / LR1121 | LoRa | **SF6 with SX126x framing**, BW 500, CR 4/5 | 12 | same | 53 ms | 22.56 ms | same |
| **19 Hz 7x**, 868 / 915 | SX1276 / LR1121 (SX127x-compatible mode) | LoRa | **SF6 with SX127x framing**, BW 500, CR 4/5, LDRO off | 12 | same; sync word = chip default 0x12 (its `SetSyncWord(0x12)` is commented out) | 53 ms | 22.3 ms | [code] ML-127, ML-LR |
| FSK 50 Hz, 868 / 915 / 433 | SX126x / LR1121 | GFSK BT 1.0, **100 kb/s, fdev 50 kHz**, RX BW 312 kHz | – | 16 bits (SX126x unit = bits [spec] SX1262 §13 PacketParams), detector 8 bits | **16-bit sync = FrameSyncWord**, fixed length, chip CRC off, **whitening on** (PN9 x⁹+x⁵+1, default seed 0x0100 [spec] SX1262 §6.2.3.4 p.46) | 20 ms | 7.6 ms | [code] ML-126, ML-LR |

- **Mode availability ([code] ML-RM):**
  - 50 Hz: 2.4 GHz only.
  - FLRC: 2.4 GHz only.
  - FSK: sub-GHz only.
  - 19 Hz 7x: 915 / 868 only. It is "the only mode" that ExpressLRS 900 hardware and FrSky R9 hardware support.
  - SX127x hardware cannot talk to SX126x hardware.
- The 7x mode on the LR1121 sets register 0xF20414 (bit 18 = 1, bit 23 = 0). The comment calls it "SF6 compatible with
  SX127x, user manual 8.3.1" ([code] SXL `lr11xx.cpp`; the same register is in [code] RL-LR, MIT).
- **Air times check [calc]** (formulas in §6, BW 500 kHz, T_sym = 128 µs at SF6):
  - 91 B, SF6, CR 4/5, preamble 12, SX127x formula: 22.304 ms, which matches the code's 22.3 ms.
  - The same frame with the SX126x formula: 22.56 ms, which matches the code's 22.56 ms.
  - So the two mLRS SF6 modes differ on air by 2 symbols.

### 1.2 Frame (91 bytes, both directions; `FRAME_TX_RX_LEN 91`) [code] ML-CONF, ML-FT
| Bytes | Uplink (Tx → Rx), `tTxFrame` | Downlink (Rx → Tx), `tRxFrame` |
|---|---|---|
| 0–1 | `sync_word` (uint16) = **FrameSyncWord** | same |
| 2–6 | status, 40 bits: seq_no 3, ack 1, frame_type 3, is_32channels 1, antenna 1, rssi_u7 7, fhss_index_band 1, **fhss_index 6**, LQ_serial 7, transmit_antenna 1, spare 2, payload_len 7 | seq_no 3, ack 1, frame_type 3, spare 1, antenna 1, rssi_u7 7, LQ_rc 7, LQ_serial 7, transmit_antenna 1, spare 2, payload_len 7 |
| 7–24 | RC data (18 B): V1 = 6 B + crc1 (2) + 10 B; V2 (32 ch) = 11 B + crc1 (2) + 5 B | payload[82] (bytes 7–88) |
| 25–88 | payload[64] | |
| 89–90 | `crc` (uint16) | `crc` |

- **The bit order inside the status bytes is not given.** C bitfields are packed by the compiler (GCC on ARM:
  LSB-first, little-endian) [U].
- The uint16 fields are probably little-endian, as on the STM32 / ESP32 targets [U].

### 1.3 CRC and identity
- **CRC:** `fmav_crc` = **CRC-16/MCRF4XX**: reflected polynomial 0x1021 (0x8408), init 0xFFFF, no final XOR
  ([code] FMAV).
  - It covers bytes 0–88 (`FRAME_TX_RX_LEN − 2`).
  - It is **not seeded by the bind**: the sync word is just a covered byte ([code] ML-FR).
  - Unauthenticated uplink: the running CRC after the header + RC1 block is also stored in `crc1`, then the CRC
    continues over the rest. Whether the continuation covers the `crc1` bytes themselves is [U]; try both.
- **FrameSyncWord** = CRC-16/MCRF4XX over the 4 bytes of `bind_dblword = u32_from_bindphrase(phrase)`
  ([code] ML-SET). Then:
  - +0x1111·Ortho;
  - +0x5AA5 or +0xA55A for the two dual-band settings;
  - +0x1212·Privacy.
- The same value is:
  - the FHSS seed;
  - the **GFSK sync word**;
  - the first two bytes of **every** frame.
- The **FLRC sync word** = `bind_dblword` (+0x11111111·Ortho) ([code] ML-SET).
- **Passive identity:** a CRC-valid frame yields the 16-bit FrameSyncWord. This is a stable per-binding ID. It also
  yields fhss_index, seq_no, frame_type and LQ / RSSI fields.
- **Encryption:** in privacy mode only the payload / RC data of some frame types are encrypted. `Encrypt(…, 5)` suggests
  the 5-byte status stays clear [U] ([code] ML-FR). Key handling is out of scope.
- **Whitening:** the frame code applies none ([code] ML-FR). LoRa applies the chip's own whitening as usual. FSK uses
  the chip's PN9.

### 1.4 Channel plans and hopping [code] ML-FH, ML-CONF, ML-FHC
| Band | List | Used per binding | Fixed bind channel |
|---|---|---|---|
| 2.4 GHz | 2401–2480 MHz, 1 MHz, 80 ch | 24 (50 Hz), 18 (31 Hz), 12 (19 Hz); with Ortho ≤ 18 | 2447 MHz |
| 915 FCC | 902.4–927.6 MHz, 0.6 MHz, 43 ch | 25 | 913.8 MHz |
| 868 | 863.275–868.0 MHz, 0.525 MHz, 10 ch | 6 | 863.275 MHz |
| 866 IN | 865.375 / 865.900 / 866.425 / 866.950 MHz (the source says "needs to be adjusted") | 3 | 865.375 MHz |
| 433 | 433.360 / 433.920 / 433.480 MHz | 2 | 433.360 MHz |
| 70 cm | 430.4–449.6 MHz, 0.6 MHz, 33 ch | 18 / 12 | 436.4 MHz |

- **PRNG:** the same LCG as ELRS (214013, 2531011, mod 2³¹, output >> 16), seeded with FrameSyncWord.
- **Channel draw:** take the k-th unused channel with `prng() % (N − k)`. Reject:
  - the bind channel;
  - for all bands except 433 / IN866, a channel within ±1 of the previous one.
- **Ortho:** uses every 3rd channel.
- **2.4 GHz exclusion:** an optional Wi-Fi channel exclusion.
- **Header:** the uplink header carries fhss_index (6 bits). With dual band, it alternates which band's index it carries
  ([code] ML-TOO).
- Hop once per frame period, with the Rx reply on the same channel inside the period [U]: not stated in the files read.

### 1.5 mLRS vs ExpressLRS on air
| Cue | mLRS | ELRS | Discriminates? |
|---|---|---|---|
| 2.4 LoRa BW | SX1280 register **0x18** ([code] SXL) | 0x18 (file 07) | **No.** Both are the SX1280 "BW_0800" setting = **812.5 kHz** (file 07, [spec] SX1280 table). The 800 vs 812.5 premise is wrong. |
| IQ polarity | always normal ([code] ML-128 / ML-126 / ML-LR) | inverted when UID[5] bit 0 = 1, or while binding | inverted ⇒ not mLRS |
| Preamble (2.4 LoRa) | 12 at SF5 / 6 / 7 | 12, but 14 at SF6 (250 Hz) | SF6 with 14 ⇒ ELRS |
| Preamble (sub-GHz) | 12 | 8 (10 at SF8 / 9) | **yes** |
| CR | 4/5 (LI 4/5 at 2.4) | 4/7–4/8 sub-GHz, LI 4/6–4/8 at 2.4 | yes, after decode (implicit header) |
| LoRa sync word | chip default (0x12 / 0x1424); none set ([code] ML-127 comment) | 0x12 | no |
| Payload | 91 B | 8 / 13 B | **yes**: burst length |
| Air time / period | 7.9–23.5 ms bursts every 20 / 32 / 53 ms, sub-GHz SF5 / 6 | 1.5–10.8 ms, 2–20 ms period | **yes** |
| Sub-GHz SF6 framing | SX126x (19 Hz) or SX127x (7x) | SX127x | partly (§6) |
| FSK | 100 kb/s, ±50 kHz, 16-bit sync, 91 B, whitened | LR1121 GFSK 300 kb/s, ±100 kHz | yes |
| FLRC | 91 B, 9 ms, CRC off | 8 B, 1 / 2 ms | yes |

---

## 2. FrSky R9 / R9M (as implemented in Multiprotocol, "FrSkyR9"; stock ACCESS R9 firmware: unpublished [U])
| Item | Value | Ref |
|---|---|---|
| PHY | LoRa SF6, BW 500 kHz, CR 4/5, implicit header (`ConfigModem1(..., true)` sets bit 0 of reg 0x1D), chip payload CRC **off** (bit 2 of reg 0x1E = false), SF6 detect-optimise / threshold | [code] MPM-R9, MPM-SX |
| Preamble | 9 symbols | [code] MPM-R9 |
| Sync word / IQ | not written ⇒ SX1276 reset values: **0x12**, IQ normal (RegInvertIQ bit 6 = 0) | [code] MPM-SX + [spec] SX1276 p.114 (0x33, 0x39) |
| Uplink length | **26 B** (FCC / flex); **14 B** (EU / LBT) | [code] MPM-R9 |
| 26 B layout | 0–2 **TX ID** (3 B); 3 hop index; 4 chanskip; 5 RX number; 6 flags (0x00 normal; failsafe 0x10–0x1E; range check 0x20); 7 ?; 8–19 channels (4 × 3 B = 8 channels × 12 bits, alternating CH1–8 / CH9–16); 20 RX/TX sequence; 21 S.Port byte count; 22–23 S.Port data; **24–25 CRC-16 LE** | [code] MPM-R9, MPM-C |
| 26 B CRC | `FrSkyX_crc`: nibble table of reflected 0x1021 entries (0x1189 …) applied MSB-first, seed 0, over bytes 0–23 (the same software CRC as D16 in file 10) | [code] MPM-C |
| 14 B layout | 0–2 TX ID; 3 chanskip; 4 RX number; 5 flags; 6–11 channels (marked TODO in MPM); 12 = (RX_seq << 4) ∣ TX_seq; **13 = CRC-8 poly 0x07, init 0xFF** over bytes 0–12 | [code] MPM-R9 |
| Telemetry | 13 B from the RX in a 7.4 ms window after the uplink; bytes 9–10 = TX ID[1..2]; CRC-16 (same function, seed = TX ID[1] + TX ID[2]·256) LE in bytes 11–12 | [code] MPM-R9 |
| Period | **20 ms**: 0.46 ms settle + 11.14 ms TX slot + 7.4 ms RX + 1 ms | [code] MPM-R9 |
| Air time [calc] | 26 B, preamble 9, SX127x formula: **7.84 ms** (48 payload symbols). 14 B: **5.28 ms** | [spec] SX1276 §4.1.1.7 |

**Channel plans** (FRF = num·256·61 Hz, i.e. 15.625 kHz units [calc]):

| Variant | Plan | Hop rule |
|---|---|---|
| FCC | base 0xE200 → **903.479 MHz**, 43 ch × 499.7 kHz (to ~924.5) | `n = (n + chanskip) % 43`, chanskip 1–24 |
| 915 flex | 0xE4C0 → 914.473 MHz, 29 ch × 499.7 kHz | % 29; indices 27–28 replaced by chanskip, chanskip+1 |
| 868 flex | 0xD700 → 859.505 MHz, 29 ch × 499.7 kHz | same |
| EU | 0xD7D0 → 862.753 MHz, 19 ch × 249.9 kHz | % 19 |

**Passive identity:** the 3-byte TX ID is in the clear in bytes 0–2 of every uplink. The hop index and chanskip in
bytes 3–4 give the next channels directly.

**Versus others:** SF6, BW 500, preamble 9, ~7.8 ms burst every 20 ms. ELRS 200 Hz: preamble 8, 4.38 ms burst. mLRS 7x:
preamble 12, 22.3 ms burst.

---

## 3. Other CC2500 links
Values come from the register tables [code] (GPL-3), decoded with the CC2500 formulas (file 10 §0) [calc].

**Common chip behaviour [spec] (CC2500):**
- **Sync:** `SYNC_MODE` 011 = **30/32: the 16-bit sync word is sent twice** (p.34 Table 24).
- **Hardware CRC:** CRC-16 0x8005, init 0xFFFF, over the length + data (DN504 `culCalcCRC`).
- **Whitening:** PN9, all-ones seed, applied **before** FEC (CC2500 §15.1, p.29–30).
- **FEC:** only in fixed-length mode. Rate ½, K = 4 convolutional code with encoder table {0,3,1,2, 3,0,2,1, 3,0,2,1,
  0,3,1,2}, 4×4 matrix interleaver, trellis terminator 0x0B (one or two bytes) (CC2500 §18 p.37; DN504 §4, §6 code).
- **MSK:** "inverts the sync word and data compared to e.g. signal generators". `DEVIATN` then sets the phase-change
  fraction of the symbol, not a deviation (p.34).

| | Futaba **S-FHSS** | **Hitec** Optima / Minima | Graupner **HoTT** | **RadioLink** (Surface / Air; DumboRC; RC4G) |
|---|---|---|---|---|
| Mod / rate | **2-FSK 128.143 kBd, ±38.09 kHz** | **MSK 249.94 kBd** | **MSK 249.94 kBd + FEC ½** (125 kb/s net) | **GFSK 37.99 kBd, ±38.09 kHz** |
| RX BW | 232 kHz | 464 kHz | 542 kHz | 232 kHz |
| Preamble | 4 B | **24 B** | 4 B | 4 B |
| Sync (sent twice) | **D391** | **TX ID bytes [3],[2]** (bind: D391) | TX **D391**; RX telemetry **2C6E** | **ABCD**; DumboRC **BADC**; RC4G **ABA5** |
| Length mode | fixed **13 B** | variable (length byte); data **27 B** (0x1A + 26), bind 23 B | fixed **50 B** TX / 22 B RX | variable; **34 B** (0x21 + 33); RC4G 16 B |
| Whitening | no | **yes** | **yes** | **yes** |
| CRC | HW CRC-16 0x8005 / FFFF | HW | HW (before FEC) | HW, plus byte 33 = sum of bytes 1–32 mod 256 |
| Grid | base 2399.9996 MHz, CHANNR = 16 + 6k (k = 0–29) × 249.94 kHz → **2404.00 + 1.49963·k MHz** (to 2447.49) | 2405.600 MHz + CHANNR × 299.927 kHz, CHANNR even 0–140 (to 2447.59) | 2403.012 MHz + CHANNR × 338.013 kHz, CHANNR = 3(n+1) → **2404.026 + 1.01404·n MHz**, n = 0–74 | 2402.000 MHz + CHANNR × 399.90 kHz; CHANNR = 12i + inc (i = 0–15, inc ∈ {0, 3, 6, 9}), one slot replaced by 0xC0 + inc |
| Hops | 30 positions; `n += code + 2` (code 0–27, per TX) with wrap | **21 ch from the 32-bit ID** (LCG ×0x19660D + 0x3C6EF35F, `((r >> 8) % 141) & 0xFE`, ≤ 8 per third); bind uses CHANNR 0, 10, … 130 | 75 ch, one of **16 fixed tables** (chosen at random, not from the ID) | 16 ch: 12-ch-spaced set offset by `inc = 3·(ID[0]&3)`, shuffled by a seeded LCG from ID[0..3] |
| Timing | **2 packets per channel, 1.625 ms apart; hop every 6.8 ms** | **22.5 ms**: 4 TX packets 2.736 ms apart, then an RX window | **10 ms**: TX, then a 4.025 ms RX window | **19.9 ms** (RC4G 14.46 ms); RX window 7.4 ms |
| Air time [calc] | 184 bit → 1.44 ms | 27 B data: 456 bit → 1.82 ms | 54 B coded ×2 + 8 B → 928 bit → 3.71 ms | 36 B + 8 B → 352 bit → 9.27 ms |
| Clear IDs | byte 0 = 0x81; **bytes 1–2 = 16-bit TX ID**; byte 5 carries the hop number; bytes 11–12 the 5-bit hop code | **bytes 1–3 = TX ID** (and the sync word) | bytes 0–1 = hop-table value; byte 2 hop index; **bytes 40–44 TX ID, 45–49 RX ID** (5 B each); RX packet bytes 0–4 TX ID, 5–9 RX ID | byte 1 type (0x01 Surface, 0x21 Air, ∣0x02 telemetry request); **bytes 2–5 = 4-byte TX ID**; 28–29 pseudo counter; 32 rf_ch_num |
| Refs | [code] MPM-FU; [comm] RFENG | [code] MPM-HI | [code] MPM-HO | [code] MPM-RL |

**Notes:**
- **S-FHSS sync length conflict:** [comm] RFENG counts 168 bits per packet with a **16-bit** D391 sync (1.311 ms).
  The MPM register MDMCFG2 = 0x83 means 30/32, i.e. **32 bits** [calc]. Accept either on detection; a capture will
  settle it [U].
- **S-FHSS channel width:** RFENG describes 30 channels of 1.5 MHz, 2-FSK, 128 kb/s, no FEC, 4 channels per packet and
  2 packets per frame [comm], consistent with the table.
- **Hitec:** the 30/32 sync is the TX ID (bytes [3],[2] repeated), so the sync word itself identifies the transmitter.
  Identification needs a blind sync search (preamble 24 B of 1010 is a strong trigger).
  - The 0xF5 0xDF / 0xDB marker alternates in bytes 24–25.
  - Byte 22 = 0xEE. Byte 23 = the halving sequence 0x80, 0x40, …
- **HoTT:** FEC means the bits on air are the interleaved, convolutionally coded bits.
  - Decode: de-interleave 4×4 → Viterbi (DN504 table) → de-whiten → CRC.
  - Bytes 0–1 identify the hop table: `HOTT_hop_val[table]`.
  - **LBT:** MCSM1 = 0x00 means always CCA ([code] MPM-HO).
- **RadioLink:** the hop channels are 4.8 MHz apart (12 × 0.4 MHz). The one extra channel sits at 2478.8–2481.2 MHz.

---

## 4. OpenLRS-ng (RFM22B / Si4432)
**PHY.** Si443x packet handler, GFSK. Modem table rows (`modem_params`) [code] OL-B. Rate and deviation are
[calc] (rate = txdr·1 MHz/2¹⁶, or /2²¹ when reg 0x70 bit 5 is set; deviation = fd·625 Hz, consistent with the bps column
and the 625 Hz unit in [code] OL-RFM):

| idx | bps | r_6e/6f → rate | r_72 → deviation | r_70 |
|---|---|---|---|---|
| 0 | 4800 | 0x2752, /2²¹ → 4.80 k | 0x30 → **30 kHz** | 0x2C |
| 1 | 9600 | 0x4EA5 → 9.60 k | 30 kHz | 0x2C |
| 2 (default) | 19200 | 0x9D49 → 19.20 k | 30 kHz | 0x2C |
| 3 | 57600 | 0x0EBF, /2¹⁶ → 57.6 k | 0x2E → **28.75 kHz** | 0x00 |
| 4 | 125000 | 0x2000 → 125 k | 0xC8 → **125 kHz** | 0x00 |

- **Modulation:** r_71 = 0x23, i.e. GFSK via FIFO.
- **Coding:** r_70 bits 0–1 are clear, so **no whitening and no Manchester** [calc from the Si443x register map, U:
  datasheet not fetched (Silicon Labs 403)].

**Packet** (`rfmInit`, [code] OL-RFM):

| Element | Value |
|---|---|
| Preamble | 0x0A = 10 nibbles = **40 bits** (80 with diversity); RH22 confirms nibble units [code] |
| Sync | 2 bytes **0x2D 0xD4** |
| Header | **4-byte header** = `rf_magic` bytes 3…0, checked on all bits. On each hop, header byte 3 = `(rf_magic & 0xFF) ^ hop_index` ([code] OL-C `setHopChannel`) |
| Length | variable: length byte after the header |
| CRC | DACTL = 0x8C: packet handler on, MSB first, CRC on, crc[1:0] = 00 = **CCITT**, computed over header + length + data (CRCDONLY = 0; RH22 describes the CRC as "computed on HEADER, LENGTH and DATA"). Polynomial x¹⁶+x¹²+x⁵+1; **init not confirmed [U]** (try 0x0000 / 0xFFFF) |

**Payload** ([code] OL-TX, OL-C):

| Item | Value |
|---|---|
| Byte 0 | flags: 0x80 / 0x40 sequence bits; low bits 0x00 servo, 0x01 failsafe, 0x38–0x3F serial (1–8 B) |
| Channels | from byte 1, 4 channels in 5 bytes (10-bit), optionally + 1 byte of 4 × 2-bit channels |
| Packet size | `pktsizes[flags & 7]` = {0, 7, 11, 12, 16, 17, 21} B (default flags 0x0A → 11 B) |
| Bind packet | starts with 'b'; the reply starts with 'B' |

**Hopping and timing:**
- **Carrier:** default base 435 / 868 / 915 MHz (limits 413–463, 848–888, 895–935).
- **Channel:** base + hopchannel[i] × spacing, with spacing 5 = 50 kHz (10 kHz step units) [code] OL-B.
- **Hop list:** up to 24 hops (default 22, 10, 19, 34, 49, 41). Hop after every packet.
- **Interval:** `(bytes + 15) × 8.2e6 / bps + 2000 µs` (+ a telemetry term), rounded up to whole ms. Default 11 B at
  19.2 k → ≈ 13.1 ms + telemetry [calc].

**Passive identity:** the 32-bit `rf_magic` is the Si443x header, sent in the clear in every packet (byte 3 XORed with
the hop index). The default `0xDEADFEED` [code] OL-B means an unconfigured unit.

---

## 5. FlySky AFHDS 3, Futaba FASST / FASSTest / T-FHSS, TBS Tracer, ImmersionRC Ghost
| System | Published facts | Status |
|---|---|---|
| **AFHDS 3** | 2.402–2.480 GHz, 140 channels, 16 hops (32 JP / KR) [vendor] (03 ref). PHY modes in the EdgeTX module driver: `CLASSIC_FLCR1_18CH`, `CLASSIC_FLCR6_10CH`, `ROUTINE_FLCR1_18CH`, `ROUTINE_FLCR6_8CH`, **`ROUTINE_LORA_12CH`** (a LoRa mode exists); EMI standards FREE / CE / FCC [code] EDGE (GPL-2). Latency claims: < 10 ms (FRM302), 5 ms (PL18) [vendor] https://hobbyking.com/jp_jp/flysky-frm302-2-4ghz-afhds-3-long-range-transmitter-module-pl18-nv14-compatible-w-ppm-uart-modes.html. RF chip, modulation, rate, sync and packet: **unpublished** (FCC exhibits for FCC IDs 2A2UNNB4PLUS00 / N4ZFGR4B000 not reachable: 403) | **[U] unpublished** |
| **FASST** | Radio IC **Micro Linear ML2724**, not CC2500; not reverse-engineered [comm] https://rfengfpv.wordpress.com/2017/01/10/futaba-s-fhss-protocol-overview/ and https://www.deviationtx.com/forum/protocol-development/7583-futaba-fasst-protokoll. Patent summary: 36 channels 2.025 MHz apart, 36 hop sequences, ~2 ms burst every 8 ms; the FCC filing reportedly says "DSSS" [comm] https://www.flyinggiants.com/forums/showthread.php?t=56339 (secondhand) | **[U]** |
| **FASSTest / T-FHSS** | FASSTest = DSSS + FHSS hybrid, modes 18CH / 14CH / 12CH (12CH faster), incompatible with FASST; T-FHSS = FHSS with telemetry [vendor] https://futabausa.com/product-support/protocols/ and https://futaba.uk/pages/radio-protocols | **[U]** PHY unpublished |
| **TBS Tracer** | 250 Hz, 3 ms latency, 25 / 100 mW [vendor] https://www.getfpv.com/drone-brands/tbs/tbs-tracer-micro-tx.html; FLRC on the same Semtech 2.4 GHz chip family as Ghost; a 50 Hz LoRa option announced for firmware 6.14 [comm] https://intofpv.com/t-are-tbs-about-to-release-lora-firmware-for-tracer | **[U]**. Sync, packet and grid unpublished; if it is SX1280 FLRC, the FLRC detector of file 07 applies ("FLRC-like, non-ELRS" when the 1 / 2 ms ELRS timing is absent) |
| **ImmersionRC Ghost** | Frame rates 222.22 Hz (Pure Race), 166 Hz (Race), 62 Hz (Normal), 15 Hz (Long Range); "Chirp Spread Spectrum + adaptive FHSS"; Atto RX −117 dBm in LR [vendor] https://www.getfpv.com/discontinued/immersionrc-ghost-2-4ghz-combo-kit-1x-tx-3x-rx-3x-antennas.html; 250 Hz Pure Race / Race250 uses MSK instead of LoRa [comm] https://oscarliang.com/immersionrc-ghost/; "LoRa / FLRC" [vendor] https://www.getfpv.com/radios/radio-transmitters/immersionrc-ghost-uber-lite-2-4ghz-rc-transmitter-module.html | **[U]**. Rates are a timing signature; SF / BW / sync are unpublished |

---

## 6. SX127x SF6 implicit-header framing (closes R2e research)
**Q1. Is the first block at SF6 sent at reduced rate (SF−2 bits/symbol, CR 4/8) as at SF ≥ 7?**
**Yes on SX127x**, from these sources:
- **[spec] SX1276 §4.1.1.6 p.29–30:**
  - "The header is transmitted with maximum error correction code (4/8)".
  - "With SF = 6 selected, implicit header mode is the only mode of operation possible."
- **[spec] SX1276 §4.1.1.7 p.31:** n_payload = 8 + max(ceil((8PL − 4SF + 28 + 16CRC − 20IH) / (4(SF − 2DE)))·(CR+4), 0),
  with "SF is the spreading factor (6 to 12)". It has no SF6 exception.
  - With IH = 1 the numerator is 8PL − 4SF + 8 + 16CRC, so the first 8 symbols carry 4(SF−2) payload bits.
  - That is the reduced-rate block: at SF6, 8 symbols × 4 bits at CR 4/8 = **2 payload bytes**.
- **[code] gr-lora_sdr (GPL-3):**
  - `interleaver_impl.cc` uses `sf_app = sf − 2` and `cw_len = 8` for the first `sf − 2` codewords for every SF, with no
    sf < 7 or implicit-header exception.
  - The README lists SF 5–12 "tested with RFM95, SX1276, SX1262", with the footnote "SF 5 and 6 are not compatible with
    SX126x". So its SF6 is the SX127x variant.
- **Air-time checks [calc]:** the SX127x formula reproduces the firmware air-time constants:
  - ELRS 900 200 Hz (SF6, CR 4/7, 8 B, preamble 8): **4.384 ms**, which matches file 03's 4.38 ms;
  - mLRS 7x (91 B): 22.304 ms vs the code's 22.300.

  (ELRS's number comes from a calculation, not a measurement. The bench should confirm with one ELRS 900 200 Hz
  capture.)

**Q2. Other SF6 differences on SX127x:**
- **Receiver only:** RegDetectOptimize = 0x05 and RegDetectionThreshold = 0x0C ([spec] p.27, p.114). The transmitted
  waveform is unaffected.
- **Sync word:** RegSyncWord 0x39 (default 0x12) has no SF6 note ([spec] p.114).
  - gr-lora_sdr sends the sync word as two symbols, nibble·8, after the up-chirps for every SF, then 2.25 down-chirps
    ([code] GRL `modulate_impl.cc`).
  - At SF6, symbol values 8 and 16 (for 0x12) fit within 64 bins.
- **Preamble:** T_preamble = (n + 4.25)·T_sym for all SF ([spec] p.31). There are no extra symbols.
- **LDRO:** mandated only above 16 ms symbol time. At SF6 / BW 500 the symbol is 128 µs, so DE = 0.

**SX126x SF5 / SF6 is a different waveform** (this matters for mLRS SX126x modes and LLCC68 gear):
- **[spec] SX1262 §6.1.1.1 p.37–38:**
  - SF5 / 6 "have been modified slightly … operate in both implicit and explicit mode … incompatible with previous
    device generations. Especially, SF6 on the SX1261/2 will not be backward compatible with the SF6 used on the
    SX1276".
  - 12 preamble symbols are recommended.
- **ToA difference ([code] RL-126, MIT, which implements the datasheet §6.1.4 formula):** for SF5 / 6 the constant
  after the preamble is **6.25** instead of 4.25 symbols, and the first-block term is −4SF + **0** instead of −4SF + 8.
  So:
  - there are **2 extra symbols** before the payload;
  - the first 8-symbol block carries **4·SF** bits, i.e. full SF bits per symbol, **not reduced rate**, still 8
    symbols (CR 4/8 implied).
  - The mLRS constants confirm it: SF6 / 91 B gives 22.56 ms (SX126x) vs 22.30 ms (SX127x) [calc].
  - What the 2 extra symbols contain is **[U]**.
- **Compatibility:**
  - Semtech forum: SF6 is not compatible between SX127x and SX126x; the sync words are 0x12 ↔ 0x1424 and 0x34 ↔ 0x3444
    ([comm] SEMF; [spec] SX1262 register 0x0740 / 0x0741 defaults 0x14 / 0x24).
  - The LR1121 can switch to the SX127x SF6 framing (register 0xF20414 bit 18) ([code] RL-LR, SXL).
- **Older gr-lora (rpp0) SF6 issue:** inconclusive ([comm] GRL62).

**Decoder consequence:** the SX127x-SF6 path of R2e can use the same reduced-rate first block as SF ≥ 7, which is
gr-lora_sdr's behaviour. An SX126x-SF5/6 path would need the "+2 symbols, full-rate first block" variant, which is
unverified [U].

---

## 7. Implementability summary
| Item | Detect + identify from IQ? | Gaps |
|---|---|---|
| mLRS LoRa / FSK / FLRC | **yes**: PHY table, 91 B layout, CRC-16/MCRF4XX, FrameSyncWord in clear, hop LCG | status bit order and uint16 endianness [U]; crc1 coverage [U]; SX126x SF5 / 6 waveform details [U] |
| FrSky R9 (MPM-compatible) | **yes**: SF6 / 500 / 4/5, preamble 9, 26 / 14 B, CRCs, TX ID, grids | stock ACCESS R9 unpublished |
| S-FHSS, Hitec, HoTT, RadioLink | **yes**: full PHY + CRC + IDs (HoTT needs CC2500 FEC from DN504) | S-FHSS sync 16 vs 32 bits [U] |
| OpenLRS-ng | **yes**: rates, deviation, sync 2DD4, 4 B rf_magic header, grid | Si443x CCITT CRC init [U] |
| SF6 framing | **resolved for SX127x** (reduced-rate first block) | SX126x SF5 / 6 extra symbols [U] |
| AFHDS 3, FASST / FASSTest / T-FHSS, Tracer, Ghost | detect-only by timing / grid ("-like") | PHY unpublished |
