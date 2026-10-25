# 06 — DVB-T / COFDM drone video links

Sources:
- [S1] ETSI EN 300 744 V1.6.2 (2015-10), DVB-T: https://www.etsi.org/deliver/etsi_en/300700_300799/300744/01.06.02_60/en_300744v010602p.pdf (downloaded and read: clauses 4.4–4.6, Annexes E, F and G).
- [S2] ETSI EN 302 755 V1.4.1 (2015-07), DVB-T2: https://www.etsi.org/deliver/etsi_en/302700_302799/302755/01.04.01_60/en_302755v010401p.pdf (read: clauses 7.2.1, 9.5 and 9.8).
- Vendor sources are cited inline. Lines marked **[derived]** I computed from the standard. Lines marked **[unverified]** are vendor claims or my inference.

## 1. OFDM parameters (DVB-T) [S1 §4.4, Tables 4–5, Annexes E and G]

| Parameter | 2K | 8K |
|---|---|---|
| FFT size (Tu in units of T) | 2048 | 8192 |
| Carriers K | 1705 | 6817 |
| Kmin / Kmax | 0 / 1704 | 0 / 6816 |
| Useful data carriers per symbol | 1512 | 6048 |
| Continual pilots / TPS carriers | 45 / 17 | 177 / 68 |
| GI Δ (1/4, 1/8, 1/16, 1/32) | 512, 256, 128, 64 T | 2048, 1024, 512, 256 T |
| Ts = Tu + Δ | 2560, 2304, 2176, 2112 T | 10240, 9216, 8704, 8448 T |

The carrier index relative to the centre is k' = k − (Kmax+Kmin)/2. The centre carrier (852 in 2K, 3408 in 8K) sits at fc.

| Channel | Elementary period T | Clock 1/T | Tu 2K / 8K | Δf 2K / 8K | (K−1)/Tu occupied |
|---|---|---|---|---|---|
| 8 MHz | 7/64 µs | 64/7 = 9.142857 MHz | 224 / 896 µs | 4464.29 / 1116.07 Hz | 7.61 MHz (7.611607 in 2K, 7.608259 in 8K) |
| 7 MHz | 1/8 µs | 8 MHz | 256 / 1024 µs | 3906.25 / 976.56 Hz | 6.66 MHz |
| 6 MHz | 7/48 µs | 48/7 = 6.857143 MHz | 298.667 / 1194.667 µs | 3348.21 / 837.05 Hz | 5.71 MHz |
| 5 MHz (Annex G) | 7/40 µs | 40/7 = 5.714286 MHz | 358.4 / 1433.6 µs | 2790.18 / 697.55 Hz | 4.75 MHz |

- **Scaling rule.** Annex G says the system "can be tailored changing only the parameter elementary period T … By adjusting the system clock rate the bandwidth and the transmission capacity are modified accordingly". All four standard values fit T = 7/(8·B) µs, so the clock is 8B/7 MHz for channel bandwidth B in MHz. **[derived]**
- **Extrapolated non-standard bandwidths [derived; unverified that products use this rule]:**

  | B | Clock | 2K Δf | 2K occupied |
  |---|---|---|---|
  | 4 MHz | 4.5714 MHz | 2232.14 Hz | 3.80 MHz |
  | 2.5 MHz | 2.8571 MHz | 1395.09 Hz | 2.38 MHz |
  | 2 MHz | 2.2857 MHz | 1116.07 Hz | 1.90 MHz |
  | 1 MHz | 1.1429 MHz | 558.04 Hz | 0.95 MHz |

  Occupied bandwidth is about 0.951·B in every case.
- **Frame and superframe.** A frame is 68 OFDM symbols, numbered l = 0..67. A superframe is 4 frames [S1 §4.4, §4.6.1]. One superframe always carries an integer number of RS(204) packets [S1 §4.7].
- **Frame durations [derived]:**
  - 2K, 8 MHz, GI 1/8: 68 × 252 µs = 17.14 ms, so one 2^21-sample block at 61.44 MSps (34 ms) holds a full TPS frame.
  - 8K, 8 MHz, GI 1/8: 68.5 ms.
  - 2K scaled to 2 MHz, GI 1/8: 68.5 ms. These need more than one block (StreamAnalyzer).
- **4K mode (DVB-H only, Annex F):** K = 3409, Kmax = 3408, 89 continual pilots, 3024 data carriers.

## 2. Pilots [S1 §4.5]

**Continual pilots.** Table 7: 45 in 2K, 177 in 8K. The 8K list begins with the full 2K list. I checked that the lists are strictly increasing and that, together with the scattered and TPS cells, they leave exactly 1512 / 6048 data cells in all 4 symbol phases.

2K (45):
```
0 48 54 87 141 156 192 201 255 279 282 333 432 450 483 525 531 618 636 714 759 765 780 804 873 888 918 939 942 969 984 1050 1101 1107 1110 1137 1140 1146 1206 1269 1323 1377 1491 1683 1704
```
8K (177) = the 2K list above, then:
```
1752 1758 1791 1845 1860 1896 1905 1959 1983 1986 2037 2136 2154 2187 2229 2235 2322 2340 2418 2463 2469 2484 2508 2577 2592 2622 2643 2646 2673 2688 2754 2805 2811 2814 2841 2844 2850 2910 2973 3027 3081 3195 3387 3408 3456 3462 3495 3549 3564 3600 3609 3663 3687 3690 3741 3840 3858 3891 3933 3939 4026 4044 4122 4167 4173 4188 4212 4281 4296 4326 4347 4350 4377 4392 4458 4509 4515 4518 4545 4548 4554 4614 4677 4731 4785 4899 5091 5112 5160 5166 5199 5253 5268 5304 5313 5367 5391 5394 5445 5544 5562 5595 5637 5643 5730 5748 5826 5871 5877 5892 5916 5985 6000 6030 6051 6054 6081 6096 6162 6213 6219 6222 6249 6252 6258 6318 6381 6435 6489 6603 6795 6816
```

**TPS carriers.** Table 8: 17 in 2K, 68 in 8K. The 8K list begins with the 2K list.

2K (17):
```
34 50 209 346 413 569 595 688 790 901 1073 1219 1262 1286 1469 1594 1687
```
8K (68) = the 2K list above, then:
```
1738 1754 1913 2050 2117 2273 2299 2392 2494 2605 2777 2923 2966 2990 3173 3298 3391 3442 3458 3617 3754 3821 3977 4003 4096 4198 4309 4481 4627 4670 4694 4877 5002 5095 5146 5162 5321 5458 5525 5681 5707 5800 5902 6013 6185 6331 6374 6398 6581 6706 6799
```

**Scattered pilots.** In symbol l, the scattered pilots are k = Kmin + 3·(l mod 4) + 12p (p ≥ 0, k ≤ Kmax). The pattern repeats every 4 symbols, and each continual pilot coincides with a scattered pilot every 4th symbol [§4.5.1, §4.5.3].

**Pilot modulation and boost.**
- Continual and scattered pilots: Re{c} = 4/3 · 2(1/2 − w_k), Im{c} = 0. So w_k = 0 gives +4/3 and w_k = 1 gives −4/3, which is E[|c|²] = 16/9 (+2.5 dB) [§4.5.3–4.5.5].
- Data cells are normalised to E = 1. TPS cells use normal power (amplitude 1) [§4.6.3].
- Pilot values do not depend on l or m: they depend only on k.

**PRBS w_k** [§4.5.2, Fig. 10]:
- Polynomial: X^11 + X^2 + 1.
- Initialisation: all ones.
- One output per carrier k = Kmin..Kmax, on every used carrier whether or not it is a pilot. The first output bit belongs to the first active carrier.
- The standard states the sequence starts "1111111111100…".
- A Fibonacci implementation reproduces it: s[n] = s[n−9] ⊕ s[n−11], output = oldest stage. Its first 48 bits are `111111111110000000001100000001111000001100110001`. **[derived; I computed these 48 bits myself and did not check them against a reference decoder such as gr-dtv]**

## 3. TPS [S1 §4.6]

**Modulation.** Each OFDM symbol carries one TPS bit s_l. All TPS carriers in a symbol carry the same differentially encoded bit (DBPSK, real axis).
- In symbol 0 of each frame: Re{c} = 2(1/2 − w_k), which carries s0.
- For l > 0: s_l = 0 keeps Re{c_l} = Re{c_{l−1}}, and s_l = 1 negates it.
- The DBPSK is re-initialised at every TPS block (frame). Bits are sent left-most first.

**68-bit TPS block** (Table 9). There is 1 initialisation bit, 16 sync bits, 37 information bits (31 used, 6 zero) and 14 BCH parity bits. Bits s25–s39 sent in superframe m′ apply to superframe m′+1; all other bits refer to the current superframe.

| Bits | Content | Values |
|---|---|---|
| s0 | Initialisation | Set by w_k (see above) |
| s1–s16 | Sync word | Frames 1 and 3: `0011010111101110`<br>Frames 2 and 4: `1100101000010001` (the bitwise inverse) |
| s17–s22 | Length indicator | `010111` = 23 bits, no cell ID<br>`011111` = 31 bits, cell ID present<br>`100001` = 33 bits, DVB-H (Annex F) |
| s23–s24 | Frame number (Table 10) | 00 = 1, 01 = 2, 10 = 3, 11 = 4 |
| s25–s26 | Constellation (Table 11) | 00 = QPSK, 01 = 16-QAM, 10 = 64-QAM, 11 = reserved |
| s27–s29 | Hierarchy (Table 12) | 000 = non-hierarchical, 001 = α1, 010 = α2, 011 = α4, 1xx = Annex F<br>Annex F redefines s27 as 0 = native / 1 = in-depth interleaver, and s28–s29 as 00 non-hierarchical, 01 α1, 10 α2, 11 α4 |
| s30–s32 | HP code rate (Table 13) | 000 = 1/2, 001 = 2/3, 010 = 3/4, 011 = 5/6, 100 = 7/8, 101–111 = reserved |
| s33–s35 | LP code rate | Same coding as HP; 000 when non-hierarchical |
| s36–s37 | Guard interval (Table 14) | 00 = 1/32, 01 = 1/16, 10 = 1/8, 11 = 1/4 |
| s38–s39 | Mode (Tables 15 / F.9) | 00 = 2K, 01 = 8K, 10 = 4K (Annex F), 11 = reserved |
| s40–s47 | Cell ID (Table 15a) | cell_id b15–b8 in frames 1 and 3, b7–b0 in frames 2 and 4; zero if not provided |
| s48–s53 | DVB-H / reserved | Standard: all 0<br>Annex F: s48 = 1 means time slicing in use, s49 = 1 means MPE-FEC in use (HP in frames 1/3, LP in frames 2/4 if hierarchical); s50–s53 = 0 |
| s54–s67 | BCH parity | See below |

**BCH error protection** [§4.6.2.11]:
- Code: BCH(67,53, t = 2), shortened from systematic BCH(127,113, t = 2).
- It protects s1–s53.
- Generator polynomial: x^14 + x^9 + x^8 + x^6 + x^5 + x^4 + x^2 + x + 1. The standard labels it h(x), but it is the generator.
- Shortening: prepend 60 zero bits to the encoder input, then discard them after encoding.

**Test vectors [derived; my own encoder, not cross-checked]:**
- Frame 1, length 010111, QPSK, non-hierarchical, HP 1/2, LP 000, GI 1/8, 2K, cell 0 → parity `10101011100111`.
- Same with length 011111, 16-QAM, HP 2/3, GI 1/16 → parity `11001100000100`.

**Detection and decoding notes [derived]:**
- TPS bits are identical across all 17 or 68 carriers, so a receiver can sum the carriers after removing the w_k sign.
- The sync word or its inverse aligns the frame.
- The BCH check validates the block.
- Mode, GI, constellation and code rate decode fully without a DVB-T demodulator chain.

## 4. DVB-T2 P1 (to tell T2 from DVB-T) [S2 §7.2.1, §9.8]

**P1 symbol.** It starts every T2 frame and every FEF part.
- It is a 1K OFDM symbol, Tu = 1024T, with a C-A-B structure:
  - C: the first 542T of A, frequency-shifted by f_SH = 1/(1024T);
  - A: 1024T;
  - B: the last 482T of A, frequency-shifted by f_SH.
- Total length T_P1 = 2048T, which is 224 µs at 8 MHz (C 59 µs, A 112 µs, B 53 µs).
- **Why it matters for detection:** the C/A and A/B correlations at lags 1024T and 542T, with the f_SH frequency shift, are a unique detector. They are not a cyclic prefix.

**Carriers and modulation.**
- 384 of the 853 carriers are active, at k = 44..809 of the 1K grid (Table 68). The pattern is about 6.83 MHz wide at 8 MHz.
- Active carriers are boosted by √(853/384) (3.47 dB).
- Modulation chain:
  1. Sequence = CSS_S1 (64) | CSS_S2 (256) | CSS_S1 (64).
  2. DBPSK, starting from +1.
  3. Scrambling with the 384 first bits of the PRBS from §5.2.4, initial state `100111001000110`.

**Signalling fields.**
- S1 is 3 bits (8 sequences of 64 chips). Table 18:

  | S1 | Meaning |
  |---|---|
  | 000 | T2_SISO |
  | 001 | T2_MISO |
  | 010 | Non-T2 (FEF) |
  | 011 | T2-Lite SISO |
  | 100 | T2-Lite MISO |
  | 101–111 | Reserved |

- S2 is 4 bits (16 sequences of 256 chips).
  - S2 field 1 (3 bits) gives the FFT size and GI class. For T2-base (Table 19a):

    | Field 1 | FFT size | Guard intervals |
    |---|---|---|
    | 000 | 2K | any |
    | 001 | 8K | 1/32–1/4 |
    | 010 | 4K | any |
    | 011 | 1K | any |
    | 100 | 16K | any |
    | 101 | 32K | 1/32, 1/16, 1/8 |
    | 110 | 8K | 1/128, 19/256, 19/128 |
    | 111 | 32K | 1/128, 19/256, 19/128 |

  - S2 field 2 (1 bit) is the mixed flag: 0 means all preambles in the transmission are the same type, 1 means preambles of different types are transmitted (Table 20).

**Other differences from DVB-T.**
- T2 adds 1K, 16K and 32K FFT, GI 1/128, 19/256 and 19/128, PP1–PP8 pilot patterns and LDPC/BCH coding.
- T2 bandwidths: 1.7 MHz (T = 71/131 µs), 5, 6, 7 and 8 MHz, and 10 MHz (T = 7/80 µs, professional only) [S2 Table 65].
- T2 frames are separated by at least 10 000T between P1s [S2, FEF clause].
- **Rule:** a P1 means DVB-T2 (or a FEF). DVB-T has no P1, uses only 2K/8K (4K for DVB-H), and its TPS carries the sync word.

## 5. Drone / UAV COFDM video transmitters

| Product or class | Waveform | Bandwidth | Modes | Encryption | Bands | Source |
|---|---|---|---|---|---|---|
| Domo / DTC (now Codan) SOL7 OBTX, 2018 | **DVB-T (standard)** plus proprietary NB and UMVL | DVB-T 6/7/8 MHz<br>NB 2.5 MHz, 1.25 MHz, 625 kHz<br>UMVL 6/7/8 MHz | DVB-T: QPSK/16/64QAM, FEC 1/2–7/8, GI 1/4–1/32, 3.732–31.668 Mbit/s<br>NB/UMVL: BPSK/QPSK/8PSK/16QAM, FEC 1/3 and 2/3, GI 1/8 and 1/16 | ABS 32-bit standard; AES 128/256 licensed (export-controlled) | Swappable bands 200 MHz–8.9 GHz | https://www.broadcastrental.com/api/download-file/32766 |
| DTC SOLO family incl. SOL8SDR (brochure, Nov 2016) | DVB-T (6/7/8 MHz) with MPEG-2 and ABS in the base tier; NB 2.5/1.25 MHz; 625 kHz, UMVL, ultra-low latency in Platinum | As left | Not stated | AES 128/256 "subject to UK Export Control" | Not stated | https://www.top10vpn.com/assets/2021/01/DTC_Video_Brochure_Nov16.pdf |
| Codan acquisition note | NB is "DTC proprietary narrowband COFDM", as distinct from DVB-T | — | — | — | — | https://armadainternational.com/?p=13121 |
| IMT VSTx (Vislink) | "Standard 2k DVB-T COFDM mode"; also 1.25 / 2.5 MHz NBCOFDM (proprietary) | 2K DVB-T; NB 1.25 / 2.5 MHz | — | — | — | https://www.securityinfowatch.com/video-surveillance/video-transmission-equipment/product/10271618/integrated-microwave-technologies-vstx-cofdm-transmitter , search summary of https://www.officer.com/home/product/10827690/[sectionalias] |
| IMT SkymasterTX (UAV downlink) | COFDM (DVB-T per IMT range) **[unverified mode]** | — | up to 8 W | "BCRYPT AES 128/256" | "Wide variety of bands" | https://www.officer.com/command-hq/technology/computers-software/product/10773713/ |
| Vislink AeroLink, 2025 | COFDM DVB-T | 6/7/8 MHz | up to 31.7 Mbit/s, 8 W | AES, FIPS 140-2 | 1.7–2.2, 2.0–2.4, 4.4–5.0, 6.425–6.540 GHz | https://www.vislink.com/wp-content/uploads/2025/06/AeroLink-Datasheet-06-25c.pdf (search summary; the direct PDF download failed) |
| HawkEye I (consumer drone COFDM, 2016) | DVB-T; unencrypted version "can be received by standard DVB-T tuners" | 6/7/8 MHz | FEC 1/2–7/8, GI 1/32–1/4 | Encrypted variant not STB-compatible | 200–2400 MHz | https://www.cnx-software.com/2016/05/23/hawkeye-i-is-a-compact-lightweight-cofdm-dvb-t-video-transmitter-with-up-to-10-km-range/ |
| Cheap Chinese "COFDM video transmitter" modules (asmag / made-in-china listings) | "COFDM 2K", some 2K/8K | **2 / 2.5 / 4 / 8 MHz adjustable** | QPSK/16QAM (64QAM on higher-power units), FEC 1/2–7/8, GI 1/32–1/4, ~1 W (30 dBm) up to 15–25 W | AES 128/256 | 300–900 MHz, 1.2, 1.4, 2.4, 4.4 GHz options; one 300–4400 MHz | https://www.asmag.com/product_detail.aspx?id=35991 , https://www.asmag.com/suppliers/inquiry.aspx?id=32987 , https://www.asmag.com/product_detail.aspx?id=42002 (fetch gave 403; content from search summaries) **[unverified]** |
| "Ranger mini" body-worn / drone COFDM | COFDM | 1 / 2 / 4 / 8 MHz | — | — | — | https://www.securityworldmarket.com/uk/Newsarchive/ranger-mini-digital-video-transmitter1 (search summary only) **[unverified]** |
| iWave UAV downlink (2.3 / 2.4 / 2.5 GHz) | COFDM with **LDPC**: not DVB-T FEC, so proprietary | 4 / 8 MHz | — | AES-128 | 2304–2390, 2402–2482, 2500–2570 MHz | https://iwavecomms.com/12km-tcpipudp-cofdm-uav-video-and-data-downlink-product |

**Observations:**
- **Standard vs proprietary.** Professional units (Domo/DTC, IMT/Vislink) advertise an interoperable DVB-T mode at 6/7/8 MHz, which is the standard 2K framing. Their narrowband modes (625 kHz–2.5 MHz) and UMVL are proprietary COFDM, with different FEC (1/3, 2/3), BPSK/8PSK and GI 1/8 or 1/16 only.
- **Low-cost modules.** The 2/2.5/4/8 MHz modules call themselves "2K COFDM". **[unverified]** Possibly they reuse DVB-T framing with T scaled as in Annex G; if so, Δf = (8B/7 MHz)/2048 and the continual/TPS pilot sets above apply unchanged. This must be confirmed on a capture: measure Δf and look for TPS carriers and the sync word.
- **Typical modes.** 2K, QPSK or 16QAM, GI 1/8 or 1/16 (1/32 also offered), FEC 1/2–2/3 for range. These are vendor menus; no source states the defaults.
- **Encryption.**
  - Domo ABS (32-bit) and AES 128/256; IMT/Vislink AES; consumer units AES-128/256.
  - **[inference]** In DVB-T-mode links the encryption scrambles the transport stream payload. The PHY (pilots, TPS, FEC) stays standard, so detection and TPS decoding still work, and only the video is unreadable. This fits HawkEye, whose encrypted version is not set-top-box compatible.
- **BISS.** BISS-1 (12-hex-digit session key) and BISS-E exist for DVB broadcast contribution (https://en.wikipedia.org/wiki/Basic_Interoperable_Scrambling_System). I found BISS only on Vislink's satellite encoder (DVE6100), not on any COFDM drone downlink. The Soliton and Ultra Electronics searches returned nothing useful on COFDM or BISS. Soliton Zao is, to my knowledge, a cellular-bonding product rather than COFDM **[unverified]**.
- **Frequency ranges seen.**
  - 170/200–900 MHz (UHF), 1.2–1.4 GHz, 1.7–2.5 GHz, 4.4–5.0 GHz, 6.4–6.5 GHz (public-safety downlinks: twelve channels at 6.4–6.5 GHz, licensed). Source for the 6.4 GHz channels: search summary of https://www.eenewseurope.com/en/expanding-the-use-of-aerial-downlinks/.
  - The 170 MHz lower edge in the task brief was not confirmed by any source; the lowest I saw is 200 MHz.
