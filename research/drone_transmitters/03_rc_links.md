# Drone / RC control links: PHY reference for identification from IQ (research 2026-10-21)

**Tags:**
- **[S]** read in source code, or decoded from datasheet registers;
- **[M]** measured (repo or a cited third party);
- **[U]** uncertain or vendor claim.

The best sources (ExpressLRS, Multiprotocol TX module (MPM), mLRS, OpenLRSng) are **GPL-3**: reference / oracle only.
Copies are in `src_gpl_reference/`.

**Chip facts:**
- CC2500: 26 MHz crystal, default sync 0xD391.
- CYRF6936: channels 2400 + n MHz, 1 Mchip/s; 8DR mode = 250 kb/s with 32-chip codes, 125 kb/s with 64-chip.
- A7105: 2400 + 0.5·n MHz.
- nRF24: 2400 + n MHz.

## 1. Main table
| Protocol | Bands | Modulation | Rate | BW / grid | Channels / hop | Frame | Preamble / sync | Decodable / IDs | Ref |
|---|---|---|---|---|---|---|---|---|---|
| **ExpressLRS 2.4 LoRa** (SX1280 / LR1121) [S] | 2400.4–2479.4 | LoRa CSS, implicit header, LoRa CRC off, CR **LI** 4/6 (500 Hz), LI 4/8 otherwise | 500 Hz SF5; 333 Hz SF5 (8 ch); 250 Hz SF6; 150 Hz SF7; 100 Hz SF7 (8 ch); 50 Hz SF8 | BW **812.5 kHz**, 1 MHz grid | 80 ch; hop every 4 packets (2 at 50 Hz) | 2 / 3.003 / 4 / 6.666 / 10 / 20 ms; air time 1.51 / 2.37 / 3.30 / 5.87 / 7.61 / 10.8 ms | preamble 12 (14 at SF6); sync 0x12 (SX1280 0x1424); **IQ inverted when UID[5] bit 0 = 1, and always while binding** | 8 B (OTA4) / 13 B (OTA8); CRC-14 (poly 0x2E57) or CRC-16 (0x3D65), init = (UID[4]<<8 ∣ UID[5]) ^ (OTA_VERSION<<8) ^ nonce. SYNC packet: fhssIndex, nonce, rate enum, TLM ratio, Gemini, UID4 / UID5 | ELRS `common.cpp`, `OTA.h/.cpp`, `FHSS` |
| **ELRS 2.4 FLRC** [S] | same | FLRC BT 1.0, CR 1/2 | 1000 Hz, 500 Hz, DVDA 500 / 250 (packet ×2 / ×4 at 1 ms) | 650 kb/s, 0.6 MHz | 80 ch, hop every 2 packets | 1 / 2 ms, air time 389 µs | preamble 32 bits; **32-bit sync word = UID[2..5] (UID[5] ^ version)**; HW CRC seed = OtaCrcInitializer | **the sync word identifies the binding** | same |
| **ELRS 900 LoRa** [S] | FCC915 903.5–926.9 (40 ch); EU868 863.275–869.575 (13); **IN866 865.375–866.950 (4 ch, 0.525 MHz)**; AU915 (20); EU/AU433 (3), US433 (8) | LoRa CR 4/7 (4/8 on 8-ch rates) | 200 Hz SF6; 100 Hz SF7; 100 Hz 8 ch SF6; 50 Hz SF8; 25 Hz SF9; 50 Hz DVDA (SF6 ×4); LR1121 + 250 Hz SF5, 200 Hz 8 ch SF5 | BW 500 kHz | hop every 4 (2 at 25 Hz); sync channel = count/2 | 5–40 ms; air time 4.38 / 8.77 / 6.69 / 18.6 / 30.0 ms | preamble 8 (10 at SF8/9); sync 0x12; same IQ rule | same | same |
| ELRS LR1121 GFSK [S] | 900 / 2.4 | GFSK 300 kb/s, fdev 100 kHz, BT 1.0 | 1000 Hz (900 8 ch, 2.4), DVDA 500 / 250 | RX BW 467 kHz | | 1 ms, air time 658–690 µs | preamble 16 | | |
| ELRS Gemini / dual band [S] | Gemini: same packet at the same time at hop index **+count/2**. Dual band: 900 SF6 + 2.4 SF7 together (150 Hz, 100 Hz 8 ch) | | | | | | | **two simultaneous copies = strong signature** | `FHSS.h` |
| **mLRS** [S] | 2.4 (24 / 18 / 12 ch at 50 / 31 / 19 Hz); 868 (6); 915 FCC (25); **IN866 (3, 865.375–866.950)**; 433 (2) | LoRa CR LI 4/5 at **BW 800** (SF5 / 6 / 7 = 50 / 31 / 19 Hz); sub-GHz SF5/6 BW 500; **FLRC 111 Hz**; **GFSK 50 Hz (100 kb/s, fdev 50 kHz, whitening)** | | 0.5 / 0.6 MHz grids | | 20 / 32 / 53 / 9 ms | preamble 12 (FSK 16); implicit, chip CRC off; **normal IQ** | **91-byte fixed frame** both ways; air time 7.9–23.5 ms | mLRS `sx-drivers`, `fhss.h` |
| **TBS Crossfire** [M, 3rd party] | 868: 860.165–885.905; 915: 902.165–927.905; 868 CE: 863.093–868.581 (112 kHz grid); AU | **150 Hz = FSK 85.1 kBd, shift 42.48 kHz** (SX1272); 50 Hz = LoRa [U]; bind / firmware update = LoRa | | 260 kHz grid | uplink ch 0–49, downlink +50 channels; RACE: 100 shared; pseudo-random, repeats after 150 pairs | 6.666 ms; downlink ≤ 2.6 ms after the uplink | ? | uplink 23 B + CRC-16/0x1021 (reflected, seeded per hop); downlink 13 B + CRC-8/0x07; 10-bit sticks | g3gg0.de |
| TBS Tracer [U] | 2.4 | SX1280, FLRC per forum | fixed 250 Hz | ~0.6–1.2 MHz | | 4 ms | | | TBS manual |
| ImmersionRC Ghost [U] | 2.4 | SX1280 LoRa "CSS + adaptive FHSS"; pure race maybe FLRC / MSK | 15 / 55 / 160 / 222 (250 / 500) Hz | ~0.8 MHz | | 4.5–67 ms | | | |
| **FrSky D8** (CC2500) [S] | 2404–2474 | GFSK **31.0 kb/s**, dev 31.7 kHz | | **299.927 kHz** grid (CC2500 ch 0–235) | 47 ch, fixed table indexed by a counter (`(i·0x1E)%0xEB`) | 9 ms | 4 B preamble, 0xD391 | length, TX ID, chip CRC | MPM `FrSkyD_cc2500.ino` |
| **FrSky D16 / ACCST v1** [S] | same | GFSK **70.0 kb/s**, dev 57 kHz; LBT/EU **2-FSK 100 kb/s**, dev 70 kHz, base 2405 | | 300 kHz | 47 ch, `next = (n + chanskip) % 47` | 9 ms (FCC 5.2 + 3.4); LBT 4 + 4.2 | 4 B + 32-bit sync | 30 B (FCC) / 35 B (LBT): TX ID (2–3 B), RX num, chanskip + index in byte 4, **software CRC-16 (table 0x1189)** | MPM `FrSkyX_cc2500.ino` |
| **FrSky ACCST v2** [S] | same | GFSK **77.0 kb/s**, chip CRC | | | hops from the ID: inc = id%46+1, ch = 5·((inc·i)%47) + id%5 | 9 ms | | hop table from the 16-bit TX ID | `FrSkyX2_init_hop` |
| FrSky ACCESS, TD / Tandem [U] | 2.4; 900; dual | proprietary | 24 ch | | | 12–23 ms latency | | | – |
| **FrSky R9 (ACCST 900)** [S] | FCC 904–925 (43 ch, 0.5 MHz); 915 flex 915–929 (29); EU 863.25–867.75 (19, 250 kHz); 868 flex 859–873 | **LoRa SF6 BW 500 CR 4/5, implicit, preamble 9** | | 500 kHz | `(n + chanskip) % N` | 20 ms | 0x12 | 26 B: index, TX ID, CRC-16 + CRC-8 | MPM `FrSkyR9_sx1276.ino` |
| **FlySky AFHDS 2A** (A7105) [S] | 2400.5–2482 | GFSK (500 kb/s [U]) | | 0.5 MHz, ch 1–164 | **16 ch from the TX ID (LCG)**, 4 sub-bands, hop every packet | **3.85 ms** | A7105 ID code 0x5475C52A | 38 B TX / 37 B RX, 4 B TX ID + RX ID; telemetry every 5th | MPM `AFHDS2A_a7105.ino` |
| FlySky AFHDS (1) [S] | same | GFSK | | 0.5 MHz | 16 ch from a table row/offset by ID | ~1.5 ms | | TX ID | `FlySky_a7105.ino` |
| FlySky AFHDS 3 [U] | 2402–2480 (140 ch) | ? | | | 16 hops (32 JP/KR) | | | | vendor |
| **Spektrum DSM2 / DSMX** (CYRF6936) [S] | 2403–2475, 1 MHz | **DSSS 8DR, 64-chip → 125 kb/s at 1 Mcps**; bind SDR | | ~1 MHz | DSM2: 2 fixed ch; **DSMX: 23 ch from the 4-byte mfg ID (LCG, 3–75, 8/7/8 per third, parity = mfgID[3])** | **11 / 22 ms, CH1 + CH2 packets ~4 ms apart**; DSM2 2F 16.5 ms | preamble 0x333304; SOP row = ch%5 (DSMX (ch−2)%5), column = (id0+id1+id2+2)&7 | 16 B; CRC-16 seed ~(id0<<8∣id1) on CH1, (id0<<8∣id1) on CH2; byte 12 = type; **mfg ID → hops, CRC seed, SOP = unit identity** | MPM `DSM*.ino`; GRCon "Hacking DSMx" |
| Walkera DEVO (CYRF) [S] | 2.4 | DSSS 8DR 32-chip 250 kb/s | | 1 MHz | 3 "best" channels at power-up | 1.2 ms steps | SOP from mfg ID | payload XOR mfg ID | `Devo_cyrf6936.ino` |
| **Futaba S-FHSS** (CC2500) [S] | 2400–2443.5 | **2-FSK 128.1 kb/s**, dev 38.1 kHz | | 1.5 MHz hop grid | 30 ch: `ch += code+2` (code 0–27) | **6.8 ms**, each packet twice | 4 B, 0xD391 | 13 B, 2-byte TX ID, chip CRC | `Futaba_cc2500.ino` |
| Futaba FASST / FASSTest / T-FHSS [U] | 2.4 | GFSK CC2500-class | ~2 ms hop | | | | | | – |
| **Hitec AFHSS / Minima** (CC2500) [S] | 2405.6–~2447 | **MSK 250 kb/s**, whitening | | 300 kHz grid, even ch 0–140 | 21 ch from the TX ID | **22.5 ms**: 4 data + RX slot | **24 B preamble; sync = TX ID bytes** | 13 B, chip CRC | `Hitec_cc2500.ino` |
| **Graupner HoTT** (CC2500) [S] | ~2403–2428+ | **MSK 250 kb/s + FEC ½ + whitening** | | 338 kHz | 75 ch, one of 16 fixed tables | **10 ms** | 4 B; 0xD391 (bind) / 0x2C6E | 50 B TX / 22 B RX, chip CRC | `HOTT_cc2500.ino` |
| **RadioLink** (CC2500) [S] | 2402 + 0.4·n | GFSK 38 kb/s, dev 38 kHz, whitening | | 400 kHz | 16 ch from ID | **20 ms** (RC4G 14.46 ms) | 0xABCD / 0xBADC / ..A5 | 11-bit channels, chip CRC | `RadioLink_cc2500.ino` |
| JR DMSS [U]; Multiplex M-LINK (CYRF) [S] | 2.4 | DMSS FHSS + DSSS | | | | | | | `MLINK_cyrf6936.ino` |
| **Hubsan** (A7105) [S] | 2410–2465 | GFSK | | single channel, one of 12 (0x14–0x82 × 0.5 MHz) by ID | **no hopping** | ~10 ms | ID = A7105 ID code | 16 B + checksum, session ID | `Hubsan_a7105.ino` |
| **Bayang / toy XN297** [S] | 2.4 | GFSK **1 Mb/s**, XN297 scrambled | | 1 MHz | 4 ch | **2 ms** (5 with telemetry) | **XN297 preamble 0x71 0x0F 0x55**, 5 B scrambled address, CRC-16 (scrambled xorout) | 15 B with TX address | `Bayang_nrf24l01.ino`, `XN297_EMU.ino` |
| WLtoys V2x2 (nRF24) [S] | 2.4 | GFSK 1 Mb/s (MR101 250 kb/s) | | 1 MHz | 16 ch from 4 tables, 2 packets per channel | **4 ms** | fixed address 66 88 68 68 68 | 16 B | `V2X2_nrf24l01.ino` |
| **OpenLRS-ng** (RFM22B / Si4432) [S] | 433 / 868 / 915 | GFSK 4.8 / 9.6 / 19.2 / 57.6 / 125 kb/s (dev 30 / 30 / 30 / 28.75 / 125 kHz) | | 10 kHz × spacing (default 50 kHz) | default 6-ch hop list (bind-time) | rate-dependent | Si443x header = 32-bit **rf_magic** (bind ID) | bind ID | openLRSng `binding.h` |
| Dragon Link V3 [U] | 433 (868 / 915) | FSK FHSS | ~50 Hz | | | | | | – |
| **DJI RC-N1 / RC2 uplink (OcuSync 2/3/O4)** [M, ours] | 2.4 / 5.8 | **15 kHz OFDM, 72 ZC subcarriers**, LTE-like | | **~1.08–1.14 MHz** | hops ~4–6 ms | burst **≈ 513 µs** | ZC | encrypted | repo `dji` |
| Siyi / Skydroid / Herelink uplink [U] | 2.4 / 5.8 | narrow TDD uplink on an OFDM SoC (AR8020/8030, Pinecone) | | | | | | | – |
| BLE gamepads [S] | 2402–2480 | GFSK 1M/2M, h ≈ 0.5 | | 2 MHz | 37 data ch, interval 7.5–15 ms | | AA 0x8E89BED6 (adv) | adv: MAC, name | BLE Core; `rc_link` |

### ExpressLRS details [S]
- **UID** = first 6 bytes of MD5(`-DMY_BINDING_PHRASE="<phrase>"`) (`build_flags.py`).
- **Hop sequence:**
  - seed = UID[2..5] (UID[5] ^ OTA_VERSION);
  - LCG: seed = (214013·seed + 2531011) mod 2³¹, output = seed>>16;
  - 256 entries in blocks of freq_count; each block leader = sync channel (count/2), the rest shuffled.
  - A known phrase therefore reproduces the sequence. The CRC-14 seed confirms UID[4..5]; the FLRC sync word gives UID[2..5].
- **Telemetry** (ratio 1:N, 1:2 … 1:128): the RX transmits in place of one uplink slot on the same channel, i.e. a
  second emitter (different CFO / power) inside the hop.
- **IN866:** 4 channels × 0.525 MHz from 865.375. The grid separates it from EU868.

## 2. Telling them apart from IQ
1. **Modulation family:**
   - **Chirp (LoRa):**
     - BW 812.5 kHz → SX1280 2.4 (ELRS / mLRS [800!] / Ghost).
     - BW 500 sub-GHz SF5–9 → ELRS 900, mLRS, R9 (SF6 CR4/5 preamble 9), Crossfire 50 Hz.
     - Then separate by SF + preamble + CR + air time:
       - ELRS: CR 4/7–4/8 sub-GHz, LI at 2.4; preamble 8 / 10 / 12 / 14.
       - mLRS: CR 4/5, preamble 12, 91 B frame (long air time).
       - R9: SF6 CR 4/5, 26 B every 20 ms.
   - **DSSS:** a 1 MHz, 1 Mcps burst on a 1 MHz grid = CYRF6936 (DSM2 / DSMX / DEVO / M-LINK). DSMX: 23 ch, burst
     pairs ~4 ms apart every 11 / 22 ms; DSM2: 2 frequencies.
   - **GFSK / FSK / MSK:** bit rate + deviation:

     | Bit rate | Candidates |
     |---|---|
     | ~31 kb/s | D8 |
     | 38 kb/s | RadioLink |
     | 70 / 77 / 100 kb/s | D16 v1 / v2 / LBT |
     | 85.1 kBd, shift 42.5 kHz, sub-GHz | Crossfire 150 Hz |
     | 100 kb/s, 50 kHz dev, sub-GHz | mLRS FSK |
     | 128 kb/s 2-FSK | S-FHSS |
     | 250 kb/s MSK | HoTT / Hitec |
     | 300 kb/s, fdev 100 kHz | ELRS LR1121 GFSK |
     | 500 kb/s | FlySky A7105 [U] |
     | 650 kb/s FLRC | ELRS / Tracer [U] / mLRS FLRC |
     | 1 Mb/s | XN297 / nRF24 toys, BLE (h 0.5, AA 0x8E89BED6) |
     | 2 Mb/s | ESB (`esb.py`) |
   - **OFDM:** 15 kHz, ~72 sc, ~0.5 ms hops = DJI RC uplink-like. Wider = video / combined SoC.
2. **Grid** (hop-centre residue; remove the receiver ppm first, since ~20 kHz is a big fraction of a 300 kHz step):

   | Grid | Protocols |
   |---|---|
   | 299.927 kHz from 2404.0 | FrSky D8/D16 |
   | 299.927 kHz from 2405.6 | Hitec |
   | 299.927 kHz from 2405.0 | D16 LBT |
   | 1.5 MHz from 2400 | S-FHSS |
   | 338 kHz | HoTT |
   | 400 kHz from 2402 | RadioLink |
   | 500 kHz | A7105 |
   | 1 MHz at 2400.4 + n | ELRS 2.4 |
   | 1 MHz at 2400 + n | CYRF / nRF24 |
   | 260 kHz from 860.165 / 902.165 | Crossfire |
   | 0.6 MHz from 903.5 | ELRS FCC915 |
   | 0.525 MHz from 863.275 / 865.375 | ELRS EU868 / IN866 |
   | 0.5 MHz from 904 / 915 | R9 |
3. **Timing** (packet period / dwell):

   | Period | Protocols |
   |---|---|
   | 1 / 2 / 4 / 6.67 / 10 / 20 ms, hop every 2–4 packets | ELRS |
   | 2 ms | Bayang |
   | 3.85 ms, hop every packet | AFHDS 2A |
   | 4 ms | V2x2 / Tracer |
   | 6.67 ms, up/down 50 ch apart | Crossfire 150 Hz |
   | 6.8 ms, pairs | S-FHSS |
   | 9 ms | FrSky D8/D16 |
   | 10 ms | HoTT |
   | 11 / 22 ms, ~4 ms pairs | DSM |
   | 20 ms | RadioLink / R9 |
   | 22.5 ms, 4 + RX slot | Hitec |

   A 34 ms block holds 1.5–34 packets; for the 20 ms-class protocols use `StreamAnalyzer` / several snapshots.
4. **Bidirectional structure:**
   - ELRS: same frequency, 1 slot in N.
   - Crossfire: +50 ch, ≤ 2.6 ms later.
   - FrSky / HoTT / Hitec: RX slot in the frame.
5. **Unit identity** (cheapest first):
   1. ELRS FLRC: sync word = UID[2..5].
   2. ELRS LoRa: IQ inversion = UID[5] bit 0; CRC seed → UID[4..5]; SYNC packet → UID4/5 (XOR model-match 0x3F) + hop index; hop sequence from UID[2..5].
   3. DSMX: channel set / CRC seed / SOP column from the 4-byte mfg ID.
   4. D16 v2: hops from the 16-bit TX ID; packets carry ID bytes.
   5. AFHDS 2A / Hubsan / Bayang: TX ID in payload or address.
   6. Hitec: sync = TX ID.
   7. OpenLRS-ng: rf_magic.
   8. DJI / BLE / OFDM: carrier ppm.
6. **Gap:** the ELRS 2.4 payload uses SX1280 LI coding (undocumented), and the standard chain fails (bench 2026-10-20).
   FLRC and sub-GHz LoRa use standard coding, so they are decodable.

## 3. Open (needs captures)
- Crossfire 50 Hz LoRa parameters / sync / IN plan.
- Tracer, Ghost.
- ACCESS / TD.
- FASST / FASSTest / T-FHSS.
- AFHDS 3.
- DMSS.
- Dragon Link.
- Siyi / Skydroid / Herelink uplinks.
- DJI RC2 / O3 / O4 uplinks (only the Mini 2 is measured).
- A7105 rate / deviation.

## 4. Sources
**ExpressLRS (GPL-3)** — https://github.com/ExpressLRS/ExpressLRS/blob/master/src/
- `src/common.cpp`, `include/common.h`
- `lib/OTA/OTA.h`, `lib/OTA/OTA.cpp`
- `lib/FHSS/FHSS.cpp`, `lib/FHSS/FHSS.h`, `lib/FHSS/random.cpp`
- `src/tx_main.cpp`
- `lib/SX1280Driver/SX1280.cpp`, `lib/SX127xDriver/SX127xRegs.h`
- `python/build_flags.py`

**Multiprotocol TX module (GPL-3)** — https://github.com/pascallanger/DIY-Multiprotocol-TX-Module/tree/master/Multiprotocol
- FrSky: `FrSkyDVX_common.ino`, `FrSkyD_cc2500.ino`, `FrSkyX_cc2500.ino`, `FrSkyR9_sx1276.ino`
- CYRF6936: `DSM.ino`, `DSM_cyrf6936.ino`, `Devo_cyrf6936.ino`
- A7105: `AFHDS2A_a7105.ino`, `FlySky_a7105.ino`, `Hubsan_a7105.ino`, `A7105_SPI.ino`
- CC2500: `Futaba_cc2500.ino`, `Hitec_cc2500.ino`, `HOTT_cc2500.ino`, `RadioLink_cc2500.ino`
- nRF24: `Bayang_nrf24l01.ino`, `V2X2_nrf24l01.ino`, `XN297_EMU.ino`
- `Protocols_Details.md`

**Other source code (GPL-3)**
- mLRS — https://github.com/olliw42/mLRS/blob/main/mLRS/Common/: `common_conf.h`, `fhss.h`, `sx-drivers/*.h`
- OpenLRS-ng — https://github.com/openLRSng/openLRSng/blob/master/openLRSng/binding.h

**Other protocols**
- Crossfire:
  - https://www.g3gg0.de/default/fpv-analysis-of-tbs-crossfire/
  - https://www.team-blacksheep.com/new-support/kb/4000112664/
- Tracer:
  - https://www.team-blacksheep.com/media/files/tbs-tracer-manual.pdf
  - https://oscarliang.com/tbs-tracer/
  - https://intofpv.com/t-are-tbs-about-to-release-lora-firmware-for-tracer
- Ghost:
  - https://immersionrc.com/?p=6107
  - https://oscarliang.com/immersionrc-ghost/
- DSMX:
  - https://pubs.gnuradio.org/index.php/grcon/article/view/97
  - https://wiki.paparazziuav.org/wiki/DSM
- CYRF6936 8DR:
  - https://community.infineon.com/t5/Other-Technologies-General/8DR-support-in-CYRF-modules/td-p/125239
  - datasheet https://www.farnell.com/datasheets/97396.pdf
- AFHDS 3: https://rcvisions.com/products/flysky-fgr4-4-channel-afhds3-2-4ghz-receiver-noble-nb4

**Datasheets and regulation**
- TI CC2500 (SWRS040), Semtech SX1276 / SX1280 datasheets
- FCC 15.247: https://www.ecfr.gov/current/title-47/chapter-I/subchapter-A/part-15/subpart-C/subject-group-ECFR2f2e5828339709e/section-15.247
