# Drone telemetry links and identification broadcasts (research 2026-10-21)

**Tags:**
- **[S]** source code read;
- **[D]** datasheet or standard;
- **[P]** paper;
- **[M]** measured in this repo (testdata / bench);
- **[U]** uncertain.

`src_sik/` holds the SiK firmware sources (ArduPilot/SiK, **BSD-2**, usable in C++): `radio_443x.c`, `radio_446x.c`,
`tdm.c`, `freq_hopping.c`, `golay.c`, `packet.c`, `radio.h`, `board_hm.h`.

## 1a. Telemetry / data links
| Link | Bands | PHY | Rates | Channels / hopping | Frame / sync | Decodable | Reference (licence) |
|---|---|---|---|---|---|---|---|
| **SiK** (3DR / Holybro / mRo; Si1000 = 8051 + Si443x) | 433.05–434.79 (10 ch), 470–471 (10), 868–870 (10), 915–928 (50 ch max) [S] | 2-GFSK (reg 0x71 = 0x23); PN9 whitening on (reg 0x70 bit 0); Manchester optional ≤ 128 kbps [S] | 2, 4, 8, 16, 19, 24, 32, 48, **64 (default)**, 96, 128, 192, 250 kbps [S] | FHSS: spacing = (max − min)/(n + 2); start = min + spacing/2 + random offset; order = shuffle seeded by NETID (default 25), or by the CRC16 of the AES key when encrypted [S] | preamble 64 bits 1010…; sync 2 bytes, Si443x reset value **0x2D 0xD4** (not written by SiK) [S, derived]; header 2 bytes = NETID; length; hardware CRC-16. Golay mode: software header {netid lo, hi, len} + {crc lo, hi, len}, Golay(24,12), own CRC16 [S] | **Yes:** NETID, MAVLink v1/v2, TDM trailer (13-bit window + command / bonus / resend bits), RSSI / radio-status packets. Not when AES is on [S] | ArduPilot/SiK **BSD-2** |
| RFD900x / ux / 868x | 902–928; 868 | GFSK FHSS [U] | 4, 64, 125 (≈ "128"), 250, 500 (+750) kbps [D/U] | ≤ 50 ch, NETID 0–255, P2P / multipoint / async mesh [D] | SiK-like; over-the-air SiK-compatible only in 128 kbps compatibility mode [D, inferred] | MAVLink if SiK-like [U] | binary |
| **mLRS** | 2.4 (2401–2480, 1 MHz; 24 / 18 / 12 used at 50 / 31 / 19 Hz); 915 FCC (0.6 MHz, 25 used); 868 (863.275–869.575, 6 used); **IN866 (3 used)**; 433 (2); 70 cm [S] | LoRa SX128x SF5/6/7 **BW 800 kHz** CR4/5 LI (2.4: 50 / 31 / 19 Hz); SX126x SF5/6 BW 500 kHz (sub-GHz 31 / 19 Hz); FLRC 650 kbps (111 Hz); GFSK 100 kbps Fdev 50 kHz (50 Hz FSK) [S] | 4100 / 2562 / 1547 B/s | FHSS, one hop per frame [S] | LoRa preamble 12, **normal IQ**, fixed length [S] | PHY yes; payload carries compressed MAVLink ("MavlinkX") [D] | olliw42/mLRS GPL-3 (oracle) |
| ExpressLRS telemetry | FCC915 903.5–926.9 (40 ch); EU868 (13); **IN866 865.375–866.95 (4 ch)**; AU915 / 433; ISM2G4 2400.4–2479.4 (80) [S] | LoRa (sub-GHz SX127x / 126x; SX1280 BW 812.5 kHz), FLRC; **inverted IQ** [S/M] | 2.4 "100 Hz Full" = SF7 / 812.5 kHz, 13-byte OTA8, every 10 ms [M] | FHSS; telemetry 1 slot every N packets (1:2 … 1:128) [U] | | PHY; payload blocked by SX1280 CR_LI [M] | ExpressLRS GPL-3 |
| TBS Crossfire | 868 / 915 | LoRa 50 Hz, FSK 150 Hz [U] | | proprietary FHSS | | PHY | – |
| Dragon Link V3 | 433 "ARS band" FHSS 1 W [D] | ? | | FHSS | | PHY | – |
| **LoRa MAVLink modules** (MicoAir …) | sub-GHz | **LoRa SF5 / 500 kHz, implicit header, sync 0x12, preamble 8, CR4/5, CRC, at 890.03 MHz, every 112 ms** [M, LIVE_LORA_SF5_890MHz] | | fixed or hopping | | yes (CRC) [M] | gr-lora_sdr GPL-3 (oracle) |
| **Microhard pDDL / pMDDL** (900 / 1624 / 1800 / 2280 / 2350 / 2450 / 2550 / 4000 / 5824 variants) | 2.402–2.478 (2450) [D] | 802.11n-derived OFDM 2×2, LDPC, AES [D]. **Measured: 802.11n at an 11 MHz clock, mirrored, Tu 5.82 µs, GI 1/4, 9.8 MHz, OUI 00:0f:92** [M, D_20260211-161341] | ≤ 25–27 Mbps; 8 MHz MIMO / 4 MHz SISO [D] | fixed channel | HT PLCP at the scaled clock [M] | PHY + MAC headers; CCMP payload [M] | repo `wifi` |
| Microhard P900 (Holybro) | 902–928 FHSS [D] | ? | 57.6, 115.2, 172, 230, 276 kbps [D] | FHSS | | PHY | – |
| Doodle Labs Mesh Rider | 600 MHz–6 GHz; 902–928 + 2400–2482; Helix M1–M6 1625–2500 [D] | COFDM BPSK–64QAM, AES-256 [D]; reportedly Atheros 802.11 + translator [U] | 80 / 40 / 20 / 12 Mbps at 20 / 10 / 5 / 3 MHz [D] | software width | likely scaled 802.11 [U] | PHY (+ MAC if 802.11) [U] | – |
| Silvus StreamCaster | 300 MHz–6 GHz | MN-MIMO COFDM | | 20 / 10 / 5 (2.5 / 1.25); Lite 1.25 / 2.5 / 5 [D] | | PHY | – |
| Herelink | 2.4 | **LTE Rel-8 DL, 5 ms TDD, PCI 411, +12.6 ppm, 2410.63 MHz** [M] | | | PSS / SSS / CRS | cell-level | repo |
| SIYI HM30 / MK15 / MK32 | 5.8 (MK15E: 2.4) [D] | ? (LTE-like SoC class?) [U] | | | | PHY | – |
| **XBee / DigiMesh 2.4** | 2405–2480, ch 11–26 = 2405 + 5(k − 11) | **802.15.4 O-QPSK, 2 Mchip/s, half-sine, 32-chip DSSS, 250 kbps** [D] | | fixed | SHR 4 × 0x00 + SFD **0xA7**, PHR length | **yes:** MAC frame, 64-bit addresses, payload (AES-CCM if on) | gr-ieee802-15-4 GPL-3 (oracle) |
| XBee-PRO 900HP / HPX | 902–928 | FSK FHSS [D] | 10 / 200 kbps | 64 ch, 25 hop | proprietary | PHY | – |
| LTE modem telemetry | India B1 / B3 / B5 / B8 / B28 / B40 / B41 uplink | SC-FDMA PUSCH / PUCCH / SRS, PRACH | | eNB-scheduled | | no | – |
| Iridium SBD (9602/9603) | 1616–1626.5 | DE-QPSK 25 ksym/s, 41.667 kHz channels, TDMA 90 ms, uplink slots ~8.28 ms [U] | | | unique word | burst type only | gr-iridium / iridium-toolkit |

## 1b. Identification broadcasts
| System | Band | PHY | Timing | Structure | Fields | Reference |
|---|---|---|---|---|---|---|
| **ASTM F3411-22a / EN 4709-002, BT4 legacy** | ch 37 / 38 / 39 | BLE 1M | Location ≥ 1 Hz; Basic ID / System ≤ 3 s | ADV_NONCONN_IND; AD = [len 0x1E][type **0x16**][UUID **0xFFFA** (FA FF)][app code **0x0D**][counter][25-byte message] [S] | all, in the clear | opendroneid-core-c **Apache-2.0**; receiver-android Apache-2.0 |
| F3411, BT5 Long Range | same | LE Coded PHY (S = 8 / S = 2): ADV_EXT_IND on primary channels → AUX_ADV_IND on a data channel | | same header, usually a **Message Pack** (0xF, ≤ 9 × 25 B) | same | same |
| F3411, Wi-Fi Beacon | 2.4 (ch 6) / 5 (ch 149) | 802.11 beacons | **5 Hz** on another channel | vendor IE OUI **FA-0B-BC**, type **0x0D**, counter + Message Pack [S] | same | `wifi.c` |
| F3411, Wi-Fi NAN | ch 6 (opt. 5 GHz) | NAN Sync Beacon + SDF (action frame, OUI type 0x13) | 1 s | Service ID **88 69 19 9D 92 09**; cluster **50-6F-9A-01-00-FF**; DA 51-6F-9A-01-00-00; service control 0x10 [S] | same | same |
| **DJI DroneID** (OcuSync 2/3/O4) | 2.4 (2399.5 / 2414.5 / 2429.5 / 2444.5 / 2459.5), 5.8 (5756.5 / 5776.5 / 5796.5) | 1024 FFT @ **15.36 MSps**, 15 kHz, **601 sc** (~9 MHz); CP 80 on first and last, 72 on the others; **9 symbols, ZC roots 600 / 147 in symbols 4 / 6** (0-based 3 / 5); legacy Mavic 2: 8 symbols (no first, ~576 µs) [P/S] | ~**640 ms** [P] | QPSK; Gold scrambler seed **0x12345678**; LTE turbo + rate matching; 176-byte block, **CRC-24A**; inner **91-byte** payload, **CRC-16/Kermit**; struct `<BBBHH16siihhhhhhQiiiiBB20sH` | serial (16), drone lat/lon (int / 174533 = deg), altitude, height, vN/E/U, yaw, GPS time (ms), **pilot lat/lon**, **home lat/lon**, product type, UUID, sequence | proto17 **MIT**; RUB-SysSec **AGPL** (oracle) |
| DJI Wi-Fi DroneID (Spark, Mavic Air, Mini SE) | Wi-Fi | beacons | | vendor IE OUI 26:37:12 | same | repo `wifi/droneid.py` |
| ADS-B 1090ES (ping200X / ping20Si; pingRX = receive only) | 1090 | PPM 1 Mbps, DF17 / 18 | ~2 Hz | | ICAO, callsign, position, velocity; **emitter category B6 = UAV** | repo `adsb` |
| UAT (ping2020, US only) | 978 | CPFSK 1.041667 Mbps h ≈ 0.6; 36-bit sync; RS FEC; 144 / 272 bits [U] | 1 Hz | | same | flightaware/dump978 **BSD-2** |
| **FLARM-compatible** | EU 868.2 / 868.4; **India 866.0 (SoftRF plan)**; US 902.2–927.8 (65 × 400 kHz) [S] | **2-FSK 100 kbps ±50 kHz, Manchester**; sync = Manchester(F531FAB6) = 55 99 A5 A9 55 66 65 96; 24-byte payload + CRC-CCITT [S] | ~6 ms every 0.6–1.4 s | payload XXTEA-obfuscated (keys in SoftRF GPL); aircraft type 13 = UAV | ID, position, velocity, type | SoftRF GPL-3 (oracle); FTD-12 public description |

## 2. Identification and decoding

### SiK (highest priority)
**Identification cues:**
- **Deviation ≈ bit rate (h ≈ 2):** Δf = 625 Hz × fd. From `radio_443x.c` tables, verify on air:

  | Rate (kbps) | 2 | 4 | 8 | 16 | 19 | 24 | 32 | 48 | 64 | 96 | 128 | 192 / 250 |
  |---|---|---|---|---|---|---|---|---|---|---|---|---|
  | Deviation (kHz) | 1.875 | 3.75 | 8.125 | 16.25 | 18.75 | 23.75 | 31.875 | 48.125 | 63.75 | 96.25 | 128.1 | 158.75 (cap) |
- **Preamble:** ~64 bits of 1010.
- **TDM timing (64 kbps):** 144 µs per byte; transmit window 7140 × 16 µs ≈ **114 ms** (3 packets ≤ 252 B); silence 2 × packet latency ≈ 4.2 ms; channel change per window.
- **Channel spacing:** 250 kHz on 915 (50 channels), 145 kHz on 433 (10 channels). Window capped at 400 ms (FCC) / MAX_WINDOW 131 ms.
- **India:** 433.05–434.79 is licence-exempt. 915 SiK units sit on GSM/LTE-900, so keep SiK labels away from GSM bursts
  (`LIVE_GSM_958MHz_fsk` once gave false "sik like" labels).

**Decoding:**
1. Find sync 0x2DD4 (assumed).
2. De-whiten with PN9.
3. Read NETID (2 B) and length; check CRC-16. In Golay mode, decode the Golay(24,12) header, length and CRC first.
4. Parse MAVLink: 0xFE (v1) / 0xFD (v2), sysid, compid, msgid, CRC-X.25 with the CRC_EXTRA seed per message.
5. The last 2 bytes are the TDM trailer.

**Outputs:**
- NETID (identifies the radio pair);
- MAVLink sysid;
- HEARTBEAT (autopilot, vehicle type);
- GLOBAL_POSITION_INT / GPS_RAW_INT;
- RADIO_STATUS.

The hop sequence can be regenerated from NETID with the LCG `r = r * 1103515245 + 12345` and the SiK shuffle, which
predicts the next hops.

### Other FHSS
- **RFD900x:** SiK-like; test the h ≈ 2 GFSK hypothesis at 125 / 250 / 500 kbps.
- **mLRS vs ELRS:** **2.4 GHz LoRa BW 800 kHz + normal IQ** vs ELRS **812.5 kHz + inverted IQ**. Sub-GHz mLRS: SF5/6 500 kHz, or GFSK 100 kbps / 50 kHz.
- **ELRS telemetry:** an extra packet in the hop slot pattern.
- **Crossfire and Dragon Link:** PHY only, "-like" from band / BW / hop rate.

### OFDM data links
- **Microhard:** scaled-clock 802.11 (done in `wifi_ofdm_scaled`).
- **Doodle Labs:** try scaled 802.11 at 3 / 5 / 10 / 20 MHz equivalents [U].
- **Silvus:** blind CP-OFDM, numerology label only.
- **Herelink:** an LTE cell inside ISM spectrum is itself a Herelink cue.
- **SIYI:** unknown.

### XBee / 802.15.4
2 MHz O-QPSK on the 2405 + 5k raster. Decodable: SFD 0xA7, PHR, MAC addresses, PAN ID.

### LTE UE and Iridium
- **LTE UE:** PRB-quantised SC-FDMA uplink bursts; only "LTE UE uplink" can be stated.
- **Iridium SBD:** DE-QPSK bursts on the 90 ms frame.

### Remote ID message layout (25 bytes; byte 0 = type << 4 | protocol version 0 / 1 / 2)
| Type | Message | Contents |
|---|---|---|
| 0x0 | Basic ID | ID type (serial CTA-2063-A / CAA / UTM UUID / session ID), UA type, 20-char ID; up to 2 messages |
| 0x1 | Location / Vector | status (ground / airborne / emergency / RID failure); height type; direction; speed (with multiplier); vertical speed; lat/lon int32 × 1e-7; pressure and geodetic altitude ((raw × 0.5) − 1000 m); height; accuracies; timestamp (1/10 s since the hour) |
| 0x2 | Authentication | up to 16 pages |
| 0x3 | Self-ID | description type + 23 characters |
| 0x4 | System | operator location type; **operator lat/lon**; area count / radius / ceiling / floor; EU category / class; operator altitude; timestamp (Unix seconds since 2019-01-01) |
| 0x5 | Operator ID | 20 characters |
| 0xF | Message Pack | up to 9 messages |

Data is self-declared and unauthenticated.

**Regions:**
- **EU:** prEN 4709-002 (protocol version 1).
- **Japan:** both IDs plus an authentication signature.
- **China:** **GB 46750-2025** (from 2026-05-01) uses an incompatible flag-bitmask wire format over BT5 / Wi-Fi
  (`libopendroneidcn/` in opendroneid-core-c).
- **India:** Drone Rules 2021 / 2022 mention NPNT and a "real-time tracking beacon" to be notified later. No broadcast
  RID standard is notified, so expect imported F3411 drones [U].

### DJI DroneID decode
1. Equalise on the ZC symbols; QPSK.
2. Descramble (Gold, seed 0x12345678); turbo-decode; check CRC-24A.
3. Parse the 91-byte struct; check CRC-16/Kermit.
4. Serial, drone, **pilot** and **home** positions come out unencrypted on OcuSync 2 (Mini 2, Air 2, Mavic 2) [P].

**Variants:**
- 8-symbol legacy form.
- The local decoder's "c2" variant is the RC uplink (7 symbols, 73 carriers, ZC in 0-based symbols 0 and 6).
- **O4 / newer:** encryption unresolved. CRC-24 pass + CRC-16 fail means "DroneID present, payload opaque".

### ADS-B / UAT / FLARM
- **ADS-B out on drones is rare.** DF17 with emitter category B6 = UAV (TC 1–4) is the drone flag. India uses 1090ES only.
- **FLARM:** India 866.0 MHz, 2-FSK 100 kbps, Manchester sync; type 13 = UAV. De-obfuscation needs XXTEA keys
  (GPL SoftRF; reimplement from the FTD-12 description).

## 3. Gaps
1. SiK sync / deviation on air; RFD900x rates.
2. SIYI, Doodle Labs and Silvus numerologies.
3. Crossfire and Dragon Link PHY.
4. DroneID encryption on Mini 3/4 / O4.
5. India RID status.

**Licences:**
- Permissive (OK in C++): SiK BSD-2, opendroneid-core-c Apache-2.0, proto17 MIT, flightaware dump978 BSD-2.
- Oracle only: mLRS, ELRS, SoftRF, gr-lora_sdr, gr-ieee802-15-4, RUB-SysSec, pyModeS, mutability dump978.

## 4. Sources
**Source code:**
- ArduPilot/SiK (BSD-2): https://github.com/ArduPilot/SiK/tree/master/Firmware/radio
- opendroneid-core-c (Apache-2.0): https://github.com/opendroneid/opendroneid-core-c
- OpenDroneID receiver-android: https://github.com/opendroneid/receiver-android/blob/master/Android/app/src/main/java/org/opendroneid/android/bluetooth/BluetoothScanner.java
- mLRS: https://github.com/olliw42/mLRS
- ExpressLRS FHSS: https://github.com/ExpressLRS/ExpressLRS/blob/master/src/lib/FHSS/FHSS.cpp
- SoftRF: https://github.com/lyusupov/SoftRF
- proto17/dji_droneid: https://github.com/proto17/dji_droneid
- RUB-SysSec/DroneSecurity: https://github.com/RUB-SysSec/DroneSecurity
- flightaware/dump978: https://github.com/flightaware/dump978
- local `../droneid-decoder-reference/src` (helpers.py, droneid_packet.py, enhanced_decoder.py)

**Paper:**
- NDSS 2023: https://www.ndss-symposium.org/wp-content/uploads/2023/02/ndss2023_f217_paper.pdf

**Datasheets and vendor pages:**
- RFD900x: https://cdn.hackaday.io/files/276031180966752/RFD900x%20DataSheet.pdf
- Microhard: https://microhardcorp.com/pMDDL2450.php and https://www.microhardcorp.com/brochures/pMDDL2450.Brochure.Rev.1.3.0.pdf
- Holybro P900: https://docs.px4.io/v1.14/en/telemetry/holybro_microhard_p900_radio.html
- Doodle Labs: https://mouser.com/datasheet/2/895/miniOEM_2450-3446604.pdf
- Silvus: https://silvustechnologies.com/wp-content/uploads/2021/04/StreamCaster-Lite-4200-SL4200-Datasheet-1.pdf
- XBee-PRO 900HP: https://lv.mouser.com/pdfdocs/Digi_International_XbeePro_900HP_Datasheet.pdf
- Dragon Link: https://www.getfpv.com/dragonlink-v3-advanced-slim-complete-system.html
- SIYI forum thread: https://discuss.ardupilot.org/t/notice-a-2-9-mini-android-controller-from-siyi-10-km-range-1080p-180-ms-8-core-cpu/64810

**Secondary:**
- DroneID on newer firmware: https://www.rtl-sdr.com/tag/droneid/
- India drone rules:
  - https://www.drishtiias.com/daily-updates/daily-news-analysis/draft-drone-rules-2021
  - https://ssrana.in/articles/drone-policy-of-india/
  - https://www.unmannedairspace.info/emerging-regulations/indias-new-drone-laws-mandate-real-time-tracking-but-ban-cargo-deliveries
