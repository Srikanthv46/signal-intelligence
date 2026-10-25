# Digital drone video / data links: PHY reference (research 2026-10-21)

**Evidence tags:**
- **[M-lib]**: measured by this library on testdata (`digital_vtx/detector.py` `KNOWN_NUMEROLOGIES` / `HERELINK`, bench notes).
- **[M-3p]**: measured or reverse-engineered by a third party.
- **[D]**: datasheet or regulatory filing.
- **[U]**: unverified, or the researcher's inference.
- "?": unknown.

Source texts in this folder: `src_ndss.txt` (NDSS 2023 DroneID paper), `src_bender.txt` (Bender DroneID preprint),
`src_lb1.txt` (Lightbridge 1 FCC test report).

## 1. DJI

| System | Bands | BW | SCS / FFT / CP | Frame / TDD | Sync | Antennas | Notes | Src |
|---|---|---|---|---|---|---|---|---|
| **DroneID burst** (OcuSync 2/3; O1 Mavic Pro also seen) | 2399.5–2474.5 MHz in 15 MHz steps; 5741.5–5831.5 MHz in 15 MHz steps [M-3p]. proto17: 2399.5/2414.5/2429.5/2444.5/2459.5, 5756.5/5776.5/5796.5 | 600 used + DC, ~9 MHz | 15 kHz, 1024 FFT @ 15.36 MSps, CP 72 (80 on symbols 1 and 9, NDSS; Bender swaps them) | **9 symbols (~644 µs)**; Mavic 2 and older send **8 (576 µs)**. Every ~640 ms (~600 ms per proto17), hopping centres | ZC roots 600 (symbol 4) and 147 (symbol 6), length 601 | – | QPSK, Gold scrambler, turbo code, CRC. Unencrypted on OS2/3. Newer models "may be encrypted" [U] | [1][2][3][4] |
| Lightbridge 1 (Phantom 3 Pro/Adv, Inspire 1) | 2.4 | OFDM downlink, BW ? | ? | Uplink GFSK FHSS: 36 channels, 2405.376 + 2.048(n−1) MHz, 20 dB BW 2.05 MHz, 1.88 ms pulses, ~83 hops/s [D] | ? | MIMO downlink [D] | | [5][6] |
| Lightbridge 2 (P4/P4P v1, Inspire 2, M600) | 2.4 (5.8 later [U]) | 9.2–9.6 MHz [M-lib] | 10.94 kHz, CP 1/8 (802.16e 1024 @ 11.2 MSps) [M-lib, once] | bursty | ? | ? | identity not confirmed | lib |
| OcuSync 1.0 (Mavic Pro, Goggles) | 2.4 only, 2403.5–2477.5 on a 1 MHz raster | 1.4 / 10 / 20 [D] | 15 kHz [U] | ? | ? | ? | | [4][7][8] |
| OcuSync 2.0 (Mini 2/2SE/4K, Air 2, Mavic 2, P4P v2) | 2.4 / 5.8 | 10 / 20 | 15 kHz LTE-like, CP 4.69 / 5.21 µs [M-lib] | ? | ZC | 2T2R-class SoC [U] | Uplink: 15 kHz OFDM hops ~1.1 MHz [M-lib]. Mini 3: 3.5 MHz hops, 0.56 ms on / 5.96 ms [M-3p] | lib, [1][9] |
| **O3 / O3+** (Air 2S, Mini 3/3P, Mavic 3, Avata) | 2.4 / 5.8 | 10 / 20 / 40 | **15 kHz, 2048 @ 30.72 MSps, ±600 bins at 20 MHz, CP 144/160** [M-3p, Air 2S] | ? | ZC length 601 per branch, "dual-comb" | **SFBC (Alamouti across adjacent subcarriers), 2 TX** [M-3p]; Mavic 3 Pro 2T4R [D] | Mavic 3 Pro uplink: 4.78 MHz hops, 1.7 ms / 6 ms [M-3p]; "Pigeon" P1 SoC | [10][11][9] |
| **O4** (Mini 4 Pro, Air 3/3S, Avata 2, Neo/Flip) | 2.4 / 5.1 / 5.8 | 10 / 20 / 40 / 60 | **30 kHz, 1024 @ 30.72 MSps, ±300, CP 72/88**: 13 short + 1 long-CP symbol per 0.5 ms (2192 / 2224 @ 61.44) [M-lib] | **2.5 ms TDD** [M-lib] | ZC root 330 length 601 [M-3p]; comb-2 reference at odd slot positions, same sequence on every O4 unit [M-lib] | 2T4R [D] | x.5 MHz raster [M-lib]. Uplink: 30 kHz, 4.4 MHz × 0.41 ms every 2.5 ms [M-lib]; Avata 2 0.41 / 2 ms [M-3p]. 10 MHz mode = 300-subcarrier slice offset ~150 bins | lib, [10][9][12] |
| O4+ (Mini 5 Pro, Mavic 4 Pro [U]) | | ? | 15 kHz, 2048, ±600, CP 144/160 [M-3p, one repo] | | ZC root 232 length 1201 | | verify | [10] |
| **DJI FPV O1** (Air Unit, Caddx Vista, Polar, Runcam Link, DJI FPV drone) | 5.725–5.850. 25 Mbps channels: 5660, 5695, 5735, 5770, 5805, 5878, 5914 (FCC) + CH8 5839; CE: 5735/5770/5805/5839 | **20 MHz (25 Mbps), 40 MHz (50 Mbps**, channels 1/2/3/8) | ? | ? | ? | ? | DJI FPV Combo uplink: 5 MHz hops, 0.64 ms / 4 ms; DroneID-like 10 MHz; video 40 MHz [M-3p] | [13][14][9] |
| **O3 Air Unit** | TX 5.725–5.850; **2.4 GHz RX only** | 40 MHz: 5794.5; 20 / 10 MHz: 5768.5 / 5804.5 / 5839.5 [D] | as O3 [U] | | | 2 antennas [D] | | [15] |
| **O4 Air Unit / Pro** | 5.170–5.250, 5.725–5.850 | 60 / 40 MHz: 5794.5; 20 / 10 MHz: 5768.5 / 5789.5 / 5814.5 [D] | as O4 [U] | | | **Lite 1T1R, Pro 2T2R** [D] | 100 / 130 Mbps; EIRP < 30 / < 33 dBm | [16] |
| Enterprise: M300, M30/M350/Dock 2 (O3 Ent.), Matrice 4/Dock 3 (O4 Ent.) | 2.4 / 5.1 / 5.8 | ? | O3 / O4 numerology [U] | | | M350 2T4R, Matrice 4 2T4R [D] | optional LTE dongle | [17][18] |
| Wi-Fi / Enhanced Wi-Fi (Spark, Mavic Air 1, Mini/Mini SE) | 2.4 / 5.8 | 802.11 at a reduced clock | **802.11a at a 5 MHz clock: 78.125 kHz, Tu 12.8 µs, GI 3.2 µs, ~4.4 MHz** [M-lib] | CSMA | STF/LTF | ? | DroneID in the beacon vendor IE (OUI 26:37:12, subcmd 0x10/0x11, ~200 ms) | lib, [4][19] |
| Tello | 2.4 | 802.11n 20 MHz AP [U] | 312.5 kHz | | | | | |

## 2. Walksnail / Artosyn AR8030
| System | Bands | BW | SCS / FFT / CP | Frame | Sync | Antennas | Notes | Src |
|---|---|---|---|---|---|---|---|---|
| **Walksnail Avatar HD** (V1/V2/Pro/Mini 1S/Moonlight) | 5.725–5.850 (some 5.1 / 5.6–6); channels like DJI FPV (5695 …) | 20 MHz (25 Mbps); 40 MHz at 50 Mbps [U] | **802.16e 20 MHz profile: 2048 @ 22.4 MSps, 10.94 kHz, CP 1/8, 1729 used (18.9 MHz)** [M-3p] — **matches our bench measurement exactly** | **frame 7.163 ms** (6-bit counter, wraps at 64 frames = 458 ms) | symbols 0–1: BPSK headers, one per antenna; DVB-T-like PRBS pilots | **2 TX, STBC matrix A over symbol pairs** — **matches our bench STBC finding** | 802.16e LDPC n=2304, r 1/2…5/6, 64QAM. Goggle uplink: 938 µs every frame, ~4.7 MHz (carriers +6…+438), BPSK — **matches the bench (0.94 ms, 4.8 MHz, BPSK)** | [20][21] |
| AR8030 family: Walksnail Ascent, BETAFPV ArtLynk, HGLRC Draco, Skyzone SteadyDigital, StartRC VT5, Potensic Atom 2 [U] | 19 centres: 5658–5917 (37 MHz step), 5758/5788/5828, 5340–5620 (40 MHz step) | **5 / 10 / 20 MHz** video, 2.5 MHz return | 802.16e-style; 10 MHz = 1024 @ 11.2 MSps (10.94 kHz) [U]; CP ? | ? | ? | 1T default; 2T STBC / 2T2R MIMO supported | LDPC 576/1152/2304, 25 MCS (BPSK–256QAM) | [22][23] |

Note: our bench frame period was 6.67 ms vs 7.163 ms here. Ours came from energy gaps (~0.7 ms gap + ~6 ms burst), so
check with the 6-bit frame counter in the header before trusting either.

## 3. HDZero / Shark Byte / Connex
**Measured 2026-10-22** (HDZero Freestyle V2, user-confirmed, R1, `testdata/live/LIVE_HDZERO_5658MHz`; the library's
`digital_vtx/hdzero.py`): 1024 FFT at 27 MSps (SCS 26.367 kHz, Tu 37.926 µs), CP 32 = 1/32, 889 used subcarriers
(23.44 MHz at −6 dB, 28 MHz at −20 dB: the "27 MHz" mode), 26 continual pilots incl. DC, scattered pilots every 32
carriers shifting 8 per symbol, a ~17-symbol repeated preamble once per video frame (16.66 ms, 60 fps), continuous,
1 TX. Data carriers: near-constant amplitude per carrier (kurtosis ~1.2), no clean QAM grid (JSCC).

| System | Bands | BW | PHY | Frame | Antennas | Notes | Src |
|---|---|---|---|---|---|---|---|
| **HDZero** (Race V1/V2/V3, Freestyle V1/V2, Whoop / Lite, Nano) | R1–R8 5658…5917 (37 MHz), F2 5760, F4 5800, E1 5705 [U]; L1–L8 5362–5621; CE 5732–5843 | **~27 MHz**; Narrow ~17 MHz (Nano 90 at 540p60) [U] | OFDM with **JSCC**: video samples mapped straight onto constellation points; numerology not public | continuous, locked to the video frame (60/90 fps) [U] | 1 TX | **No uplink.** 25 mW–1 W | [24][25][13] |
| Fat Shark Shark Byte (TX5S.1, TX5M.1, RX5.1) | DJI 25 Mbps channel set [U] | ? | Divimath OFDM (HDZero's predecessor) | ? | ? | discontinued | [26] |
| Connex / Amimon ProSight | 5.8 | ~40 MHz [U] | WHDI JSCC OFDM [U] | ? | MIMO [U] | | – |

## 4. Autel, Parrot, Skydio, others
| System | Bands | PHY | Uplink | Antennas | Src |
|---|---|---|---|---|---|
| **Autel SkyLink 1/2/3** (EVO II v1–v3, Lite+, Nano, Max) | 2.4 / 5.8; 900 MHz (FCC) from v2; 5.1 / 5.7 per one retailer | **15 kHz LTE-like, 600 sc (~9 MHz)** [M-3p]; Nano video ~20 MHz | 72 sc (1.1 MHz); Nano 4 MHz hops / 4 ms | EVO II v3 1T2R [D] | [27][9][28] |
| Parrot Anafi / USA / Ai | 2.4 / 5.8 Wi-Fi a/b/g/n; Ai + 4G LTE UE | standard Wi-Fi | Wi-Fi | ? | [29][9] |
| Skydio X10 / X10D | Connect SL 2.4 + 5.15–5.85; X10D MH 1.8–2.47; 5G | ? (Skydio 2 = Wi-Fi [U]) | ? | ? | [30] |
| Yuneec Typhoon H / H520 | 5.8 Wi-Fi video, 2.4 control [U] | | | | |
| FIMI X8SE | 2.4 [D], "TDMA" / RokLink 3.0 | ? | | | [31] |
| Hubsan Zino (SyncLeas) | 5.725–5.850 [D] | ? | | | [32] |
| Potensic Atom (PixSync 3/4) | 2.4 [D]; Atom 2 = AR803x [U] | | | | [32][23] |
| Holy Stone / Wi-Fi toys | 2.4 / 5.8 phone AP, standard 802.11 [U] | | separate 2.4 GHz GFSK RC [U] | | |

## 5. Third-party links
| System | Bands | PHY | Timing / uplink | Notes | Src |
|---|---|---|---|---|---|
| **Herelink** | 2.4 | **LTE 15 kHz, 100 RB = 18 MHz** [M-lib]; 10 MHz mode [D] | **10 ms / 5 ms half-frame; PSS/SSS roots 25/29/34 at symbols 6/5 of subframes 0/5; DL subframes 0/3/4/5/8/9** [M-lib]. Uplink hops 2.96 MHz, 0.52 ms / 5.16 ms; video 19.1 MHz [M-3p] | Pinecone S1 | lib, [9][33] |
| Siyi HM30 / MK15 / MK32 | HM30: FCC 5745–5825 (vendor 5.1–5.8); MK15/MK32: 2.4 | ? | hops: MK15 3.32 MHz, 1.55 ms / 4.99 ms; MK32 7.23 MHz, 0.47 / 5.01 ms; FT24 4.59 MHz, 1.61 ms | | [34][9] |
| Skydroid H12 / T10 (H16 / T12 similar [U]) | 2.4 | ? | H12 6.06 MHz, 0.25 / 2.97 ms; T10 5.27 MHz, 0.27 / 2.91 ms; YunZhuo H16 10 MHz hops + 20 MHz video | | [9] |
| **wfb-ng / OpenIPC / OpenHD** | 2.4 / 5.8 any channel | 802.11n/ac 20/40 MHz; RTL8812EU injects at **10 MHz** (156.25 kHz, Tu 6.4 µs); 5 MHz doesn't work | broadcast, no ACK, RS FEC (k-of-n, default 8/12), libsodium encryption | RTL8812AU 2T2R | [35][36] |
| RubyFPV | 433/868/915 SiK + 2.4/5.8 802.11 injection | adaptive MCS | | custom licence | [37] |
| Microhard pDDL / pMDDL2450 | 2.402–2.478 (also 900 / 1.8 GHz) | **OFDM 4 / 8 MHz**; pMDDL 2×2 MIMO, LDPC, BPSK–64QAM | | our Microhard capture = 802.11n at an 11 MHz clock, STBC | [38] |
| Doodle Labs Mesh Rider | 900 MHz–6 GHz | COFDM, software-defined BW (802.11-derived [U]) | | MIMO | [39] |
| Silvus StreamCaster | 300 MHz–6 GHz | MN-MIMO OFDM 20/10/5 (2.5/1.25) MHz; Lite 1.25/2.5/5 | | 2×2 / 4×4 | [40] |
| Domo / DTC (SOL8, MeshUltra) | wide | COFDM 1.25 … 20 MHz, ≤ 64QAM | | MIMO | [41] |
| **Cheap "COFDM" UAV modules; IMT / Vislink DVB-T class** | 300–1400 MHz (also to 2200 / 4400) | **DVB-T 2K style, 2 / 2.5 / 4 / 8 MHz, GI 1/4–1/32, QPSK/16QAM.** 8 MHz: 4.464 kHz, Tu 224 µs, 1705 carriers; scaled BWs scale the SCS [U] | continuous | AES on some | [42][43] |

## PHY features that separate systems
1. **SCS / Tu:**
   - **15 kHz (66.67 µs):** OcuSync 1/2/3/O3+, DroneID, Herelink, Autel, O4+?
   - **30 kHz (33.33 µs):** O4.
   - **10.94 kHz (91.43 µs, CP 1/8):** Walksnail / AR8030 / Lightbridge 2-like.
   - **802.11 family:** 312.5 kHz (wfb-ng / OpenHD / Parrot / toys), 156.25 kHz (10 MHz injection), 78.125 kHz (DJI Wi-Fi at 5 MHz).
   - **1–4.5 kHz (Tu 224–900 µs):** DVB-T COFDM.
   - **43.75 kHz:** PMP450 (infrastructure).
2. **CP:**
   - LTE 144/2048 (160 for the first symbol);
   - O4 72/1024 with an 88-sample symbol once per 0.5 ms;
   - DroneID 72 (80 for the first and last symbols);
   - 802.16e 1/8; Wi-Fi 1/4; DVB-T 1/4–1/32.
3. **TDD:**
   - O4 2.5 ms; Herelink 10 / 5 ms;
   - Walksnail 7.163 ms (+ a 938 µs uplink each frame);
   - HDZero, COFDM and analog continuous; Wi-Fi CSMA.
   - DroneID burst every ~640 ms: **9 symbols (OS2/3) vs 8 (Mavic 2 and older)**.
4. **Uplink hops (width / on / period):**
   - O4 4.4–7 MHz / 0.41 ms / 2–2.5 ms;
   - Mini 3 3.5 / 0.56 / 6;
   - Mavic 3 Pro 4.8 / 1.7 / 6;
   - DJI FPV 5 / 0.64 / 4;
   - OS2 1.1 MHz;
   - Herelink 3 / 0.52 / 5.2;
   - Skydroid 5–6 / 0.25 / 2.9;
   - Siyi MK15 3.3 / 1.55 / 5;
   - LB1 GFSK 2.05 MHz / 1.88 ms / ~12 ms.

   RFUAV values are approximate.
5. **Sync:**
   - DroneID: ZC 600 / 147.
   - Herelink: PSS 25/29/34.
   - O4: ZC 330 + a fixed comb-2 reference (system, not unit; unit = crystal ppm).
   - Walksnail: per-antenna BPSK headers + PRBS pilots.
   - O3: dual-comb ZC.
6. **Antennas:**
   - O3: SFBC pairs on adjacent subcarriers.
   - Walksnail / AR8030 2T: STBC over symbol pairs.
   - **O4 Air Unit (1T1R) vs O4 Air Unit Pro (2T2R).**
   - DJI drones: 2T.
   - HDZero: 1T.
7. **Raster and bandwidth:**
   - **Raster:** DJI uses x.5 MHz centres. Raceband 37 MHz steps mean HDZero or AR8030. DJI FPV and Walksnail share 5660/5695/5735/5770/5805/5839/5878/5914.
   - **BW sets:** O3 AU 10/20/40; O4 + 60; O1 20/40; Walksnail 20; HDZero ~27 / 17; AR8030 5/10/20; Microhard 4/8; Silvus 1.25–20.
   - **O3 Air Unit:** TX at 5.8 only, goggles uplink at 2.4.

## Open-source decoders
| Project | Licence | C++ use |
|---|---|---|
| RUB-SysSec/DroneSecurity (DroneID) | AGPL-3.0 | oracle only |
| proto17/dji_droneid | MIT | yes |
| anarkiwi/samples2djidroneid (OS2 DroneID, 15.36 / 30.72 MSps) | Apache-2.0 | yes |
| v0l/waveshark (Walksnail 20 MHz demod, AR8030 notes) | GPL-3.0+ | oracle only |
| svpcom/wfb-ng | GPL-3.0 | oracle only |
| Kismet DJI IE221 Kaitai | GPL-2 [U] | oracle only |
| huangjiansheng321-cloud/dji-drone-ocusync-research | no licence | facts only |
| RubyFPV | custom | ? |

## Unknown
- Lightbridge 2 and OcuSync 1 numerology / TDD (HDZero 27 MHz mode measured 2026-10-22; 17 MHz mode open).
- O3 TDD period.
- Autel, Siyi and Skydroid video numerology.
- Connex.
- Skydio PHY.
- Whether O4 DroneID is encrypted.

## Sources
1. NDSS 2023, Schiller et al.: https://www.ndss-symposium.org/wp-content/uploads/2023-217-paper.pdf
2. https://github.com/RUB-SysSec/DroneSecurity
3. https://github.com/proto17/dji_droneid
4. Bender: https://arxiv.org/pdf/2207.10795
5. https://fccid.io/SS3-201402241/Test-Report/Test-Report-2219514
6. https://dl.djicdn.com/downloads/lightbridge/en/Lightbridge_release_note_en.pdf
7. https://fccid.io/SS3-M1P1607
8. https://developer.dji.com/api-reference/android-api/Components/OcuSyncLink/DJIOcuSyncLink.html
9. RFUAV: https://arxiv.org/html/2503.09033v1
10. https://github.com/huangjiansheng321-cloud/dji-drone-ocusync-research
11. https://www.dji.com/mavic-3-pro/specs
12. https://www.dji.com/mini-4-pro/specs and https://www.dji.com/avata-2/specs
13. https://oscarliang.com/fpv-channels/
14. https://intofpv.com/archive/index.php/thread-23233.html
15. https://www.dji.com/au/support/product/o3-air-unit and https://www.getfpv.com/dji-o3-air-unit.html
16. https://www.dji.com/o4-air-unit/specs
17. https://enterprise.dji.com/matrice-350-rtk/specs
18. https://enterprise.dji.com/matrice-4-series/specs
19. https://github.com/kismetwireless/kismet/blob/master/kaitai_definitions_disabled/dot11_ie_221_dji_droneid.ksy
20. https://github.com/v0l/waveshark/issues/226
21. https://github.com/v0l/waveshark/issues/227
22. https://github.com/v0l/waveshark/issues/228
23. https://fpvwiki.co.uk/caddx-walksnail-ascent-lite-vtx
24. https://docs.hd-zero.com/goggles-operation
25. https://www.hd-zero.com/better-transmiition
26. https://intofpv.com/t-shark-byte-hdzero-hd-fpv-system?pid=153609
27. https://www.aerobotpro.com/news/white-paper-on-uav-interference-technology-based-on-vco-dds-and-sdr-technology-1/ (offline)
28. https://www.autelrobotics.com/productdetail/evo-max-series-old/
29. https://www.parrot.com/en/drones/anafi/technical-specifications
30. https://skydio.com/x10/technical-specs and https://support.skydio.com/hc/en-us/articles/23748914885659
31. https://www.maisondudrone.com/en/fimi-x8-se-2022-our-complete-test/
32. https://www.unmannedtechshop.co.uk/fr/products/hubsan-zino-mini-pro-4k-drone-sub-250g and https://www.maisondudrone.com/en/potensic-atom-full-review-a-real-competitor-for-dji-mini/
33. https://www.getfpv.com/cubepilot-herelink-hd-air-unit-v1-1.html
34. https://fccid.io/2AZBB-HM30 and https://siyi.biz/en/product/image-digital-link/hm30/faq/
35. https://github.com/svpcom/wfb-ng
36. https://docs.openipc.org/use-cases/fpv/net-cards/rtl8812eu/ and https://github.com/svpcom/rtl8812eu
37. https://github.com/RubyFPV/RubyFPV
38. https://www.microhardcorp.com/brochures/pMDDL2450.Brochure.Rev.1.3.0.pdf
39. https://doodlelabs.com/news/now-available-multi-band-mesh-rider-radios-unlicensed-bands
40. https://silvustechnologies.com/wp-content/uploads/2025/10/StreamCaster-4200-SC4200-Plus-Drop-In-Module-Datasheet.pdf and https://silvustechnologies.com/wp-content/uploads/2021/04/StreamCaster-Lite-4200-SL4200-Datasheet-1.pdf
41. https://www.unmannedsystemstechnology.com/2021/07/dtc-unveils-new-software-defined-radio/
42. https://www.asmag.com/product_detail.aspx?id=35973
43. https://www.cnx-software.com/?p=33985
