# 02 — GNSS: one-block acquisition, C/N0, jamming and spoofing indicators

Research note for batch item **S3** (`research/signal_gaps/PLAN.md`). Written 2026-10-25. It covers code generation
(with self-test values), an acquisition design that fits one IQ block, and jamming/spoofing indicators.
Status tags used below:
- **[ICD]**: copied from the primary ICD and checked against it;
- **[verified]**: also reproduced by an independent generator written for this note (scratch Python, see §0);
- **[lit]**: from published literature;
- **[derived]**: computed here;
- **[unconfirmed]**: not checked against a primary source, or an open question.

No code has been added to `python/` or `c++/`. Everything here still needs the bench step (bench → Python → C++).

---------------------------------------------------------------------------------------------------------------------

## 0. Summary and verification done for this note

| Signal | Centre (MHz) | Chip rate (Mcps) | Code | Period | Data | Min. power (dBW) | Codes self-tested here |
|---|---|---|---|---|---|---|---|
| GPS L1 C/A | 1575.42 | 1.023 | Gold 1023 (G1, G2 10-bit) | 1 ms | 50 bps (20 ms bits) | −158.5 (IIR-M and later) | PRN 1–63 first-10-chip octal reproduced, 63/63 match [verified] |
| GPS L5 I5/Q5 | 1176.45 | 10.23 | XA ⊕ XB (13-bit) 10230 | 1 ms | I5: 100 sps + NH10; Q5: pilot + NH20 | −157.9 (IIF), −157.0 (III) | all 126 XB advance↔initial-state pairs reproduced [verified] |
| NavIC L5 SPS | 1176.45 | 1.023 | Gold 1023 (GPS polynomials) | 1 ms | 50 sps (25 bps, r=1/2 K=7) | −159.0 | PRN 1–14 first-10-chip octal reproduced, 14/14 [verified] |
| NavIC S SPS | 2492.028 | 1.023 | Gold 1023 | 1 ms | 50 sps | −162.3 | PRN 1–14, 14/14 [verified] |
| NavIC L1 SPS (NVS-01+) | 1575.42 | 1.023 | IZ4, 10230 chips, data + pilot (+1800-chip overlay) | 10 ms | 100 sps | −155.8 total | not implemented (nonlinear generator) |
| GLONASS L1OF | 1602 + k·0.5625, k = −7…+6 | 0.511 | m-sequence 511, 1+x⁵+x⁹ | 1 ms | 50 bps + 100 Hz meander | −161 [unconfirmed] | generator matches GNSS-SDR's mapping [verified] |
| BeiDou B1I | 1561.098 | 2.046 | Gold 2047 truncated to 2046 (11-bit) | 1 ms | D1 50 bps + NH20 (MEO/IGSO); D2 500 bps (GEO) | −163 | table from ICD v3.0; generator balanced (1023/1024 ones) [verified] |
| Galileo E1-B/C | 1575.42 | 1.023 | memory codes 4092, CBOC(6,1,1/11) | 4 ms | E1-B 250 sps; E1-C pilot + CS25 | −157.25 (E1 total) | codes not embedded (see §6) |

What was run (scratch, not in the repo):
- the GPS C/A generator in both forms, two-tap selector (Table 3-Ia) and G2 delay; the PRN 38–63 first-10-chip
  octals of IS-GPS-200N Table 3-Ib all match (63/63);
- the L5 XB register advanced from all ones by the Table 3-I advance count, which gives the Table 3-I initial XB state
  for all 63 I5 and 63 Q5 PRNs. The same check passes on the SBAS PRNs 120, 121, 124 and 127 of the 2026-01 GPS.gov
  assignment list;
- NavIC: loading the Table 7 initial state with its **first character in stage 10** (the output stage) reproduces
  every first-10-chip octal (28/28). With the first character in stage 1 they do not match;
- BeiDou B1I, using the gnss-sdr register mapping (read for facts only): the codes come out balanced and PRN 1–63
  match ICD Table 4-1.

A wrong turn this note avoids: GNSS-SDR's `gps_l5_signal_replica.cc` takes the XB chip at index
`(advance + n) mod 10230` of a pre-generated XB stream. IS-GPS-705J §3.3.2.2 says XB runs its natural 8191-chip course
and restarts. Per the ICD the chip is `XA[n] ⊕ m8191[(advance + n) mod 8191]`. The two disagree for
n ≥ 10230 − advance: 1.3 % of chips for PRN 1 and 39 % for PRN 37. This is [unconfirmed] as a GNSS-SDR bug, since we
have not tested it on air. Our generator should use a register clocked from the ICD initial state, never that
indexing.

---------------------------------------------------------------------------------------------------------------------

## 1. GPS L1 C/A (IS-GPS-200N)

**Signal.** 1575.42 MHz (154 × 10.23 MHz), BPSK(1), 1.023 Mcps, 1023-chip Gold code, 1 ms period. LNAV data is
50 bps: 20 code periods per bit, and bit edges fall on code epochs. Minimum received power is −158.5 dBW (−128.5 dBm),
0 dBi RHCP, elevation ≥ 5°, for Block IIR-M, IIF and III [ICD]. The maximum is "not expected to exceed −153.0 dBW",
or −150 dBW for programmable-power IIR-M/IIF [ICD, §3.3.1.6.1]. The task's "−130 dBm" is roughly right; −128.5 dBm is
the spec minimum.

Typical C/N0 [lit]:
- 37–45 dB-Hz with a passive patch antenna at mid elevation;
- 45–50 dB-Hz open sky with an active antenna;
- above 52 dB-Hz is rare.

Thermal floor: C/N0 = P + 204 − NF (dBW → dB-Hz), so −158.5 dBW with a 3 dB NF gives about 42.5 dB-Hz.

**Generator** [ICD + verified]:
- G1 = 1 + x³ + x¹⁰ and G2 = 1 + x² + x³ + x⁶ + x⁸ + x⁹ + x¹⁰, 10-stage shift registers, both initialised to all ones
  at every 1 ms epoch. The output is G1 stage 10 ⊕ G2i.
- G2i is either G2 stage-s₁ ⊕ stage-s₂ (two-tap phase selector, PRN 1–37), or G2 delayed by D chips:
  `ca[n] = g1[n] ⊕ g2[(n − D) mod 1023]` (equivalently `np.roll(g2, D)`).
- PRN 38–63 have no two-tap form. Use the delay D, or load the "Initial G2" octal with its **MSB into stage 10** (the
  output stage) and the LSB into stage 1 [verified].
- Feedback with stages 1…10 and output from stage 10: G1 `new = s3 ⊕ s10`; G2
  `new = s2 ⊕ s3 ⊕ s6 ⊕ s8 ⊕ s9 ⊕ s10`.
- "First 10 chips (octal)" is the first 10 output chips as a binary number, first chip = MSB. PRN 1 = 1440 means
  1100100000.
- PRN 34 and PRN 37 have the same code (both 4⊕10, delay 950) [ICD].
- PRN 33–37 are reserved for other uses (pseudolites, test).
- Cross-correlation peaks of C/A Gold codes are at most 65/1023 (−23.9 dB). Any secondary peak above about −20 dB of
  the main peak needs an explanation: another strong PRN, multipath, or spoofing (§9).

| PRN | G2 taps (Table 3-Ia) | G2 delay (chips) | Initial G2 (octal, 3-Ib) | First 10 chips (octal) |
|---|---|---|---|---|
| 1 | 2⊕6 | 5 | – | 1440 |
| 2 | 3⊕7 | 6 | – | 1620 |
| 3 | 4⊕8 | 7 | – | 1710 |
| 4 | 5⊕9 | 8 | – | 1744 |
| 5 | 1⊕9 | 17 | – | 1133 |
| 6 | 2⊕10 | 18 | – | 1455 |
| 7 | 1⊕8 | 139 | – | 1131 |
| 8 | 2⊕9 | 140 | – | 1454 |
| 9 | 3⊕10 | 141 | – | 1626 |
| 10 | 2⊕3 | 251 | – | 1504 |
| 11 | 3⊕4 | 252 | – | 1642 |
| 12 | 5⊕6 | 254 | – | 1750 |
| 13 | 6⊕7 | 255 | – | 1764 |
| 14 | 7⊕8 | 256 | – | 1772 |
| 15 | 8⊕9 | 257 | – | 1775 |
| 16 | 9⊕10 | 258 | – | 1776 |
| 17 | 1⊕4 | 469 | – | 1156 |
| 18 | 2⊕5 | 470 | – | 1467 |
| 19 | 3⊕6 | 471 | – | 1633 |
| 20 | 4⊕7 | 472 | – | 1715 |
| 21 | 5⊕8 | 473 | – | 1746 |
| 22 | 6⊕9 | 474 | – | 1763 |
| 23 | 1⊕3 | 509 | – | 1063 |
| 24 | 4⊕6 | 512 | – | 1706 |
| 25 | 5⊕7 | 513 | – | 1743 |
| 26 | 6⊕8 | 514 | – | 1761 |
| 27 | 7⊕9 | 515 | – | 1770 |
| 28 | 8⊕10 | 516 | – | 1774 |
| 29 | 1⊕6 | 859 | – | 1127 |
| 30 | 2⊕7 | 860 | – | 1453 |
| 31 | 3⊕8 | 861 | – | 1625 |
| 32 | 4⊕9 | 862 | – | 1712 |
| 33 | 5⊕10 | 863 | – | 1745 |
| 34 | 4⊕10 | 950 | – | 1713 |
| 35 | 1⊕7 | 947 | – | 1134 |
| 36 | 2⊕8 | 948 | – | 1456 |
| 37 | 4⊕10 | 950 | – | 1713 |
| 38 | – | 67 | 0017 | 1760 |
| 39 | – | 103 | 0541 | 1236 |
| 40 | – | 91 | 1714 | 0063 |
| 41 | – | 19 | 1151 | 0626 |
| 42 | – | 679 | 1651 | 0126 |
| 43 | – | 225 | 0103 | 1674 |
| 44 | – | 625 | 0543 | 1234 |
| 45 | – | 946 | 1506 | 0271 |
| 46 | – | 638 | 1065 | 0712 |
| 47 | – | 161 | 1564 | 0213 |
| 48 | – | 1001 | 1365 | 0412 |
| 49 | – | 554 | 1541 | 0236 |
| 50 | – | 280 | 1327 | 0450 |
| 51 | – | 710 | 1716 | 0061 |
| 52 | – | 709 | 1635 | 0142 |
| 53 | – | 775 | 1002 | 0775 |
| 54 | – | 864 | 1015 | 0762 |
| 55 | – | 558 | 1666 | 0111 |
| 56 | – | 220 | 0177 | 1600 |
| 57 | – | 397 | 1353 | 0424 |
| 58 | – | 55 | 0426 | 1351 |
| 59 | – | 898 | 0227 | 1550 |
| 60 | – | 759 | 0506 | 1271 |
| 61 | – | 367 | 0336 | 1441 |
| 62 | – | 299 | 1333 | 0444 |
| 63 | – | 1018 | 1745 | 0032 |

Values for PRN 1–37 come from Table 3-Ia (taps and delays). The first-10-chip column is computed and matches the ICD
for 38–63, the octals that were checked in this session. PRN 1–37 octals are the well-known ICD values (1440, 1620,
1710, 1744, 1133, …) and come out of the same generator. Table 3-Ib for PRN 38–63 is from IS-GPS-200N §3.2.1.3
("Expanded Code Phase Assignments (GPS III, GPS IIIF, and subsequent blocks only)").

---------------------------------------------------------------------------------------------------------------------

## 2. GPS L5 (IS-GPS-705J)

**Signal.** 1176.45 MHz (115 × 10.23), QPSK:
- I5 carries data: CNAV at 50 bps, rate-1/2 K=7 convolutional coding → 100 sps, then a 10-chip Neuman-Hofman code at
  1 kHz;
- Q5 is a dataless pilot with a 20-chip NH code.

Each code is 10230 chips at 10.23 Mcps (1 ms). Minimum power per component is −157.9 dBW (IIF) and −157.0 dBW
(III/IIIF); the maximum is −150 dBW [ICD Table 3-V/§3.3.1.6]. The main lobe is ±10.23 MHz, so acquisition needs
≥ 20.46 MSps complex (decimate 61.44 → 20.48 MSps, ÷3).

**Generator** [ICD §3.3.2.2 + verified]:
- XA = 1 + x⁹ + x¹⁰ + x¹² + x¹³, initialised to all ones. It is short-cycled at 8190 chips: reset to all ones after
  8190 chips and again at the 1 ms epoch.
- XB = 1 + x + x³ + x⁴ + x⁶ + x⁷ + x⁸ + x¹² + x¹³. It is initialised to the PRN's initial XB state at each 1 ms epoch,
  then runs freely: its 8191 period wraps inside the millisecond.
- Stages 1–13, output = stage 13. Feedback into stage 1: XA `s9⊕s10⊕s12⊕s13`; XB
  `s1⊕s3⊕s4⊕s6⊕s7⊕s8⊕s12⊕s13`.
- Code chip n = XA(n) ⊕ XB(n), n = 0…10229.
- The initial XB state in the table is written stage 1 … stage 13, left to right. It equals the state of an all-ones
  XB register clocked `advance` times [verified for all 126 + 8 SBAS entries].

**Neuman-Hofman** [ICD]:
- NH10 = 0000110101 on I5;
- NH20 = 00000100110101001110 on Q5. This is the same NH20 as BeiDou D1.

One NH chip per 1 ms code period, first chip aligned with the data-symbol (I5) or 20 ms (Q5) boundary.

**Acquisition implication.** With NH codes, coherent integration beyond 1 ms needs the NH phase (10 or 20 hypotheses).
Non-coherent 1 ms summation ignores the NH code. Pilot+data combined acquisition (|I5|² + |Q5|²) gains about 3 dB.
L5 shares 1164–1188 MHz with DME/TACAN pulses (962–1213 MHz). The library already detects DME, so pulses must be
blanked (zero samples above k·σ) before correlation. Do not report DME as L5 jamming.

| PRN | XB advance I5 | XB advance Q5 | Initial XB I5 (binary) | Initial XB Q5 (binary) |
|---|---|---|---|---|
| 1 | 266 | 1701 | 0101011100100 | 1001011001100 |
| 2 | 365 | 323 | 1100000110101 | 0100011110110 |
| 3 | 804 | 5292 | 0100000001000 | 1111000100011 |
| 4 | 1138 | 2020 | 1011000100110 | 0011101101010 |
| 5 | 1509 | 5429 | 1110111010111 | 0011110110010 |
| 6 | 1559 | 7136 | 0110011111010 | 0101010101001 |
| 7 | 1756 | 1041 | 1010010011111 | 1111110000001 |
| 8 | 2084 | 5947 | 1011110100100 | 0110101101000 |
| 9 | 2170 | 4315 | 1111100101011 | 1011101000011 |
| 10 | 2303 | 148 | 0111111011110 | 0010010000110 |
| 11 | 2527 | 535 | 0000100111010 | 0001000000101 |
| 12 | 2687 | 1939 | 1110011111001 | 0101011000101 |
| 13 | 2930 | 5206 | 0001110011100 | 0100110100101 |
| 14 | 3471 | 5910 | 0100000100111 | 1010000111111 |
| 15 | 3940 | 3595 | 0110101011010 | 1011110001111 |
| 16 | 4132 | 5135 | 0001111001001 | 1101001011111 |
| 17 | 4332 | 6082 | 0100110001111 | 1110011001000 |
| 18 | 4924 | 6990 | 1111000011110 | 1011011100100 |
| 19 | 5343 | 3546 | 1100100011111 | 0011001011011 |
| 20 | 5443 | 1523 | 0110101101101 | 1100001110001 |
| 21 | 5641 | 4548 | 0010000001000 | 0110110010000 |
| 22 | 5816 | 4484 | 1110111101111 | 0010110001110 |
| 23 | 5898 | 1893 | 1000011111110 | 1000101111101 |
| 24 | 5918 | 3961 | 1100010110100 | 0110111110011 |
| 25 | 5955 | 7106 | 1101001101101 | 0100010011011 |
| 26 | 6243 | 5299 | 1010110010110 | 0101010111100 |
| 27 | 6345 | 4660 | 0101011011110 | 1000011111010 |
| 28 | 6477 | 276 | 0111101010110 | 1111101000010 |
| 29 | 6518 | 4389 | 0101111100001 | 0101000100100 |
| 30 | 6875 | 3783 | 1000010110111 | 1000001111001 |
| 31 | 7168 | 1591 | 0001010011110 | 0101111100101 |
| 32 | 7187 | 1601 | 0000010111001 | 1001000101010 |
| 33 | 7329 | 749 | 1101010000001 | 1011001000100 |
| 34 | 7577 | 1387 | 1101111111001 | 1111001000100 |
| 35 | 7720 | 1661 | 1111011011100 | 0110010110011 |
| 36 | 7777 | 3210 | 1001011001000 | 0011110101111 |
| 37 | 8057 | 708 | 0011010010000 | 0010011010001 |
| 38 | 5358 | 4226 | 0101100000110 | 1111110011101 |
| 39 | 3550 | 5604 | 1001001100101 | 0101010011111 |
| 40 | 3412 | 6375 | 1100111001010 | 1000110101010 |
| 41 | 819 | 3056 | 0111011011001 | 0010111100100 |
| 42 | 4608 | 1772 | 0011101101100 | 1011000100000 |
| 43 | 3698 | 3662 | 0011011111010 | 0011001011001 |
| 44 | 962 | 4401 | 1001011010001 | 1000100101000 |
| 45 | 3001 | 5218 | 1001010111111 | 0000001111110 |
| 46 | 4441 | 2838 | 0111000111101 | 0000000010011 |
| 47 | 4937 | 6913 | 0000001000100 | 0101110011110 |
| 48 | 3717 | 1685 | 1000101010001 | 0001001000111 |
| 49 | 4730 | 1194 | 0011010001001 | 0011110000100 |
| 50 | 7291 | 6963 | 1000111110001 | 0100101011100 |
| 51 | 2279 | 5001 | 1011100101001 | 0010100011111 |
| 52 | 7613 | 6694 | 0100101011010 | 1101110011001 |
| 53 | 5723 | 991 | 0000001000010 | 0011111101111 |
| 54 | 7030 | 7489 | 0110001101110 | 1100100110111 |
| 55 | 1475 | 2441 | 0000011001110 | 1001001100110 |
| 56 | 2593 | 639 | 1110111011110 | 0100010011001 |
| 57 | 2904 | 2097 | 0001000010011 | 0000000001011 |
| 58 | 2056 | 2498 | 0000010100001 | 0000001101111 |
| 59 | 2757 | 6470 | 0100001100001 | 0101101101111 |
| 60 | 3756 | 2399 | 0100101001001 | 0100100001101 |
| 61 | 6205 | 242 | 0011110011110 | 1101100101011 |
| 62 | 5053 | 3768 | 1011000110001 | 1010111000100 |
| 63 | 6437 | 1186 | 0101111001011 | 0010001101001 |

PRN 64–119 are reserved; 120–158 are SBAS (GAGAN = PRN 127, 128, 132, with I5/Q5 states in the GPS.gov assignment
list); 159–210 are other RNSS [ICD/GPS.gov 2026-01].

---------------------------------------------------------------------------------------------------------------------

## 3. NavIC (IRNSS) SPS

### 3.1 L5 and S (ISRO IRNSS SIS ICD for SPS v1.1, Aug 2017)

- **Frequencies:** L5 at 1176.45 MHz (24 MHz band, 1164.45–1188.45) and S at 2492.028 MHz (16.5 MHz band,
  2483.5–2500). The S band sits next to 2.4 GHz Wi-Fi channel 14 / ISM edge: check for Wi-Fi leakage.
- **SPS signal:** BPSK(1), 1.023 Mcps, 1023 chips, 1 ms. The RS component is BOC(5,2) data + pilot plus an interplex
  term; not public.
- **Data:** 25 bps, rate-1/2 K=7 convolutional code (G1 = 171₈, G2 = 133₈) → **50 sps**, so 20 code periods per
  symbol. There is no NH code on SPS.
- **Subframe:** 600 symbols = 12 s, 16-bit sync word 0xEB90 followed by 584 symbols.
- **Minimum received power:** L5 −159.0 dBW, S −162.3 dBW. Maximum: L5 −154.0 dBW, S −157.3 dBW [ICD Tables 4 and 5].
- **Generator** [ICD §4.1.1 + verified]: G1 = X¹⁰ + X³ + 1 and G2 = X¹⁰ + X⁹ + X⁸ + X⁶ + X³ + X² + 1, which are the
  GPS C/A polynomials. G1 starts at all ones. G2 starts at the per-PRN state below, so the initial state is the code
  phase. Output = G1(s10) ⊕ G2(s10). **Load the 10-character string with its first character in stage 10** (the first
  G2 output chip). Self-test: the first-10-chip octal column, 28/28 reproduced.

| PRN | L5 initial G2 | L5 first 10 chips (octal) | S initial G2 | S first 10 chips (octal) |
|---|---|---|---|---|
| 1 | 1110100111 | 0130 | 0011101111 | 1420 |
| 2 | 0000100110 | 1731 | 0101111101 | 1202 |
| 3 | 1000110100 | 0713 | 1000110001 | 0716 |
| 4 | 0101110010 | 1215 | 0010101011 | 1524 |
| 5 | 1110110000 | 0117 | 1010010001 | 0556 |
| 6 | 0001101011 | 1624 | 0100101100 | 1323 |
| 7 | 0000010100 | 1753 | 0010001110 | 1561 |
| 8 | 0100110000 | 1317 | 0100100110 | 1331 |
| 9 | 0010011000 | 1547 | 1100001110 | 0361 |
| 10 | 1101100100 | 0233 | 1010111110 | 0501 |
| 11 | 0001001100 | 1663 | 1110010001 | 0156 |
| 12 | 1101111100 | 0203 | 1101101001 | 0226 |
| 13 | 1011010010 | 0455 | 0101000101 | 1272 |
| 14 | 0111101010 | 1025 | 0100001101 | 1362 |

(ISRO-IRNSS-ICD-SPS-1.1 Table 7, © Indian Space Research Organization.)

**Licence/terms:** the ICD may be reproduced only for non-profit purposes with its terms reproduced and the © notice
kept, and may not be altered. The table above is a short factual extract for internal engineering. Before this note
is distributed outside the team, replace it with a pointer to the ICD. Implementing the code from the published
polynomials and states in our own C++ is not a copyright issue [unconfirmed legal view].

**India relevance:** NavIC satellites are GEO/GSO over India and always visible from India. Their Doppler is small,
|f_D| ≲ a few hundred Hz [derived from GSO range rate, unconfirmed]: the Doppler search can be very narrow. Within one
capture at 1176.45 MHz the same block holds GPS L5, NavIC L5 and Galileo E5a.

### 3.2 NavIC L1 SPS (ISRO-NAVIC-ICD-SPS-L1-1.0, Aug 2023)

- **Signal:** 1575.42 MHz, 24 MHz band, SBOC: BOC(1,1) and BOC(6,1) sub-carriers at 1.023 MHz and 6.138 MHz,
  coefficients α = √(6/11), β = √(4/110), γ = √(4/11), η = √(6/110).
- **Components:** data and pilot quadrature-multiplexed with 41.82 % / 58.18 % power.
- **Codes:** primary codes are **IZ4** (interleaved Z4-linear) of 10230 chips (10 ms) at 1.023 Mcps, from a pair of
  coupled nonlinear feedback shift registers. Data and pilot codes are tabulated for PRN 1–64 (ICD Tables 7 and 8).
- **Pilot overlay:** 1800 chips at 100 cps (18 s).
- **Data:** 100 sps.
- **Minimum received power:** data −159.6, pilot −158.2, total −155.8 dBW. Maximum −152 dBW.
- **Transmitters:** NVS-01 (launched May 2023, transmitting L1 since June 2023). NVS-02 (Jan 2025) did not reach its
  operational orbit [unconfirmed status today].

**Verdict:** implementable from the ICD (section 4, figures 6–15). The generator is a nonlinear register pair, more
work than a Gold code, and there is a single satellite. Low priority: do GPS C/A first. The 10 ms primary needs
coherent ≥ 10 ms correlation per hypothesis (FFT size 10230·k). Acquisition with the pilot primary plus an unknown
overlay chip is fine over 34 ms: 3 non-coherent sums.

---------------------------------------------------------------------------------------------------------------------

## 4. GLONASS L1OF (GLONASS ICD Edition 5.1, 2008)

- **Frequencies (FDMA):** f_k = 1602 MHz + k × 0.5625 MHz, k = −7…+6, so 1598.0625–1605.375 MHz. Antipodal satellites
  share a k. Every satellite uses the same code.
- **Ranging code:** 511-chip m-sequence, 0.511 Mcps, 1 ms. Polynomial G(x) = 1 + x⁵ + x⁹, 9-stage register
  initialised to all ones at each 1 ms epoch; **output from stage 7** [ICD; cross-checked with GNSS-SDR's mapping].
  Feedback into stage 1 = s5 ⊕ s9. First 24 chips from the all-ones start: `111111100000111101111100`, with 256 ones
  in 511 [derived].
- **Data:** 50 bps NRZ ⊕ 100 Hz meander (10 ms), giving 10 ms symbols. Each 2 s string ends in a 0.3 s time mark
  (30-chip shortened sequence at 100 sps). The exact bits are [unconfirmed] here: take them from the ICD.
- **Minimum power:** about −161 dBW for L1 C/A [unconfirmed: the ICD figure was not retrieved in this session].
- **Acquisition:** one code, 14 frequency channels. Search each channel: mix by −f_k, then the same 1 ms parallel
  code search with 511 code phases, two samples per chip, and ±5 kHz + ppm Doppler. The satellite identity is the
  channel number k (the slot needs the almanac).
- **Coverage:** a block centred at 1580 MHz with 61.44 MSps covers −30.72…+30.72 MHz, which is 1549.3–1610.7 MHz.
  GPS L1, Galileo E1, BeiDou B1I and GLONASS k = −7…+6 are all inside (band edges 1561.098 − 2.046 and
  1605.375 + 0.511 are fine).

---------------------------------------------------------------------------------------------------------------------

## 5. BeiDou B1I (BDS-SIS-ICD-B1I-3.0, Feb 2019)

- **Signal:** 1561.098 MHz, BPSK(2), 2.046 Mcps, 2046 chips (a 2047 Gold code truncated by its last chip), 1 ms.
  Minimum received power −163 dBW [ICD 4.2.5].
- **Generator** [ICD]: G1(X) = 1+X+X⁷+X⁸+X⁹+X¹⁰+X¹¹ and G2(X) = 1+X+X²+X³+X⁴+X⁵+X⁸+X⁹+X¹¹. Both initial phases are
  **01010101010** (stages 1…11). Output = G1 stage 11 ⊕ (XOR of the G2 stages in the phase assignment). Feedback into
  stage 1: G1 `s1⊕s7⊕s8⊕s9⊕s10⊕s11`; G2 `s1⊕s2⊕s3⊕s4⊕s5⊕s8⊕s9⊕s11`. Reset both every 2046 chips.
- **Phase assignment of G2, Table 4-1:**

| PRN | taps | PRN | taps | PRN | taps | PRN | taps |
|---|---|---|---|---|---|---|---|
| 1 (GEO) | 1⊕3 | 17 | 4⊕5 | 33 | 8⊕10 | 49 | 1⊕5⊕11 |
| 2 (GEO) | 1⊕4 | 18 | 4⊕6 | 34 | 8⊕11 | 50 | 1⊕6⊕9 |
| 3 (GEO) | 1⊕5 | 19 | 4⊕8 | 35 | 9⊕10 | 51 | 1⊕8⊕9 |
| 4 (GEO) | 1⊕6 | 20 | 4⊕9 | 36 | 9⊕11 | 52 | 1⊕9⊕10 |
| 5 (GEO) | 1⊕8 | 21 | 4⊕10 | 37 | 10⊕11 | 53 | 1⊕9⊕11 |
| 6 | 1⊕9 | 22 | 4⊕11 | 38 | 1⊕2⊕7 | 54 | 2⊕3⊕7 |
| 7 | 1⊕10 | 23 | 5⊕6 | 39 | 1⊕3⊕4 | 55 | 2⊕5⊕7 |
| 8 | 1⊕11 | 24 | 5⊕8 | 40 | 1⊕3⊕6 | 56 | 2⊕7⊕9 |
| 9 | 2⊕7 | 25 | 5⊕9 | 41 | 1⊕3⊕8 | 57 | 3⊕4⊕5 |
| 10 | 3⊕4 | 26 | 5⊕10 | 42 | 1⊕3⊕10 | 58 | 3⊕4⊕9 |
| 11 | 3⊕5 | 27 | 5⊕11 | 43 | 1⊕3⊕11 | 59 (GEO) | 3⊕5⊕6 |
| 12 | 3⊕6 | 28 | 6⊕8 | 44 | 1⊕4⊕5 | 60 (GEO) | 3⊕5⊕8 |
| 13 | 3⊕8 | 29 | 6⊕9 | 45 | 1⊕4⊕9 | 61 (GEO) | 3⊕5⊕10 |
| 14 | 3⊕9 | 30 | 6⊕10 | 46 | 1⊕5⊕6 | 62 (GEO) | 3⊕5⊕11 |
| 15 | 3⊕10 | 31 | 6⊕11 | 47 | 1⊕5⊕8 | 63 (GEO) | 3⊕6⊕9 |
| 16 | 3⊕11 | 32 | 8⊕9 | 48 | 1⊕5⊕10 | | |

PRN 6–58 are MEO/IGSO.

- **Self-test** [derived, from our generator, not printed in the ICD]: first 24 chips
  - PRN 1 `011001011011011011001101`
  - PRN 2 `100100100110001000111000`
  - PRN 6 `011010100100010000001001`
  - PRN 38 `110111111011011010101100`
  - PRN 63 `111100011110110100100001`

  Each code has 1023 or 1024 ones. Cross-check these against an independent oracle (e.g. a GNSS-SDR run, or MATLAB)
  before trusting them.
- **Navigation data:**
  - **D1** (MEO/IGSO): 50 bps with a 20-bit NH code at 1 kbps, NH = `00000100110101001110`. Data bits are 20 ms.
  - **D2** (GEO): 500 bps (2 ms bits), no NH. Coherent integration on GEOs is ≤ 2 ms.
  - **Sync:** 11-bit preamble `11100010010`. This is the standard BDS preamble; the ICD figure was not re-read here
    [unconfirmed].
- **Acquisition:** 1 ms coherent, ≥ 2 samples/chip, so ≥ 4.092 MSps (61.44/15 = 4.096 MSps works). In India the BDS
  GEO/IGSO satellites at 80–160° E are well visible; GEO Doppler is near zero.
- **Terms:** the ICD carries "© China Satellite Navigation Office"; Section 1 claims only "the right for final
  explanation". It has no reproduction clause of the ISRO/EU kind [ICD §1].

---------------------------------------------------------------------------------------------------------------------

## 6. Galileo E1 OS (Galileo OS SIS ICD Issue 2.1, Nov 2023)

- **Signal:** 1575.42 MHz, CBOC(6,1,1/11):
  - E1-B is data: I/NAV at 250 sps, 1 symbol per 4 ms code period;
  - E1-C is a pilot with the 25-chip secondary code CS25₁ = `0011100000001010110110010` (hex `380AD90` with 3 padding
    zeros) [ICD Table 20].
  - Primary codes are 4092 chips at 1.023 Mcps (4 ms), 50 codes each for E1-B and E1-C.
  - Minimum received power for E1 is −157.25 dBW [ICD].
- **Memory codes:** Annex C (Table 106) gives each code as 1023 hex symbols with 0 fill zeros. The Issue 2.1 text says
  "the primary codes can be accessed/saved from the attachments panel of the pdf". `pdfdetach -list` on the copy
  downloaded from gsc-europa.eu on 2026-10-25 showed **0 embedded files** [unconfirmed why: attachments stripped by the
  host?]. Options: the GSC's separate code files, if published; the Issue 1.x/2.0 PDFs, whose Annex C printed the hex
  inline [unconfirmed]; or ask the GSC. **Do not** copy `Galileo_E1.h` from GNSS-SDR (GPL-3).
- **Size:** 2 × 50 × 4092 bits = 409,200 bits = **51.2 kB** packed (100 kB as hex text). Fine to embed, as a generated
  `.inc`/`.npz` with a SHA-256 checked in a test.
- **Terms:** copyright © EU 2023. Reproduction is allowed if the "Terms of Use and Disclaimers" and Annex H go with it
  unmodified. Annex H authorises **royalty-free, non-exclusive** use of the OS SIS ICD IPRs, including the "E1 and E5
  Memory Codes", for R&D, manufacturing, commercialisation, distribution, sale, supply and maintenance of Products
  (software included). It is non-transferable and non-licensable, and each user is authorised directly by the EU.

  Practical reading: embedding the codes in our C++ is permitted. Ship a NOTICE pointing to Annex H, and do not
  relicense the codes under our own licence terms [unconfirmed legal view: get a lawyer's read before external
  distribution].
- **Acquisition:** the 4 ms primary gives FFT size 4 × 1023 × k. For a BOC(1,1)-type replica the correlation main
  peak is ±⅓ chip wide with −6 dB side peaks at ±½ chip: use ≥ 4 samples per chip (≥ 4.092 MSps; 8.192 MSps is
  better) or a sideband ("BPSK-like") method. Coherent integration is ≤ 4 ms on E1-B (symbol flips) and on E1-C
  without the CS25 phase. In 34 ms, 8 code periods non-coherent is fine.
- **Priority:** after GPS C/A and BeiDou. Its main value is the multi-constellation spoofing check in §9.

---------------------------------------------------------------------------------------------------------------------

## 7. Acquisition design for one block (GPS L1 C/A; the others are parameter changes)

### 7.1 Front end and constants

- **Production input:** IQ + fs + fc + `receiver_ppm` (the Pluto TCXO is about +7 to +10 ppm).
- **One TCXO for LO and ADC clock.** Sign convention to confirm on the bench [unconfirmed]: with an LO high by ε, a
  carrier at true f appears at baseband offset (f − fc) − ε·fc. Correct the carrier search centre by −ε·fc, i.e.
  −15.75 kHz at +10 ppm at L1. The code rate measured in samples scales by 1/(1+ε): 10 ppm over 34 ms is 0.35 chip
  of code-phase drift. Below ~0.5 chip of drift it only smears non-coherent sums. Correct anyway, by generating the
  replica at fs_eff = fs·(1+ε) (it costs nothing).
- **Doppler window.** Satellite Doppler is about ±4.2 kHz for a static user (±5 kHz design margin) plus ±1 kHz for
  vehicles. Plus the residual receiver error after the `receiver_ppm` correction, about ±1–2 ppm (drift across the
  sweep): ±1.6–3.2 kHz at L1.
  - With `receiver_ppm` known: search ±(5 + 1.6·Δppm) kHz. Use ±6.6 kHz for Δppm = 1, about 27 bins at 500 Hz.
  - Without it: ±(5 + 15.75) kHz at +10 ppm, or 0…+21 kHz one-sided if the sign is known. About 85 bins at 500 Hz.
    Then use the first acquired PRN to fix the common offset, and search the rest within ±5 kHz of it.
- **Decimation.** The C/A main lobe is ±1.023 MHz.
  - 61.44 MSps ÷ 30 = **2.048 MSps**: 2048 samples per ms (a power of two) and 2.002 samples/chip.
  - ÷ 60 = 1.024 MSps: 1024 per ms, 1 sample/chip. Half the cost, but up to about 3.5 dB worst-case (about 1.5 dB
    average) code-phase straddle loss.

  Use an anti-alias FIR to ±1.0 MHz passband: a polyphase decimator, or keep the block FFT, crop bins to ±1.024 MHz
  and IFFT. That crop reuses the pipeline's block FFT, per the latency rule. With 2^21 samples at 61.44 MSps,
  ÷30 gives 69,905 samples = 34 full ms.
- **Lower capture rates give longer blocks.** At 2.5 MSps a 2^21 block is 0.84 s, which allows much longer
  non-coherent sums (§7.4).

### 7.2 Core: parallel code-phase search with folding (all PRNs share the expensive parts)

For each Doppler bin f_d (spacing Δf), with code period N = 2048 samples at 2.048 MSps:
1. **Wipe-off and fold.**
   - Multiply the decimated block by e^{−j2π f_d t}. Use a recursive phasor in double precision, and use the true
     sample times.
   - Cut it into K-ms segments.
   - Fold each segment: sum its K consecutive 1 ms (N-sample) chunks into one N-sample vector. This is coherent
     integration over K ms, done in the time domain before any FFT ("averaging code periods").
   - Cost: one complex MAC per sample per Doppler bin, independent of the number of PRNs.
2. **Forward FFT** of each folded segment: one N-point FFT per (segment, bin), **shared by all PRNs**.
3. **Per PRN:**
   - multiply by the conjugate FFT of the PRN replica, precomputed once per (PRN, fs_eff);
   - **IFFT**;
   - |·|²;
   - accumulate over segments (non-coherent sum, M = segments).
   - This step dominates: PRNs × bins × M IFFTs.
4. **Alternative Doppler shift.** Since the 1 ms FFT bin is 1 kHz, a Doppler shift of an integer number of kHz is a
   circular shift of FFT bins. Compute the FFT once per segment for each of the 1/Δf_kHz fractional offsets
   (2 for 500 Hz), then roll. The wipe-off in step 1 is then needed only for those fractional offsets. This saves the
   forward FFTs; the IFFT count is unchanged.

**Bin spacing:** Δf = 1/(2T_coh). That is 500 Hz at 1 ms (≤ 0.9 dB worst-case loss) and 125 Hz at 4 ms. Code Doppler
within a segment is negligible: 3.2 chips/s at 5 kHz Doppler.

### 7.3 Data-bit transitions and the coherent length

- **C/A:** bits are 20 ms with edges on code epochs, unknown. In a K-ms segment the probability of a sign change is
  ≤ K/40.
  - K = 1–5 ms is safe on average.
  - K = 10 loses a segment's energy 25 % of the time.
  - K = 20 always straddles one edge unless it is aligned.
- **Robust trick for K = 10–20:** split the segment into halves A and B (2 IFFTs) and keep
  max(|a+b|², |a−b|²). This doubles the IFFT count and tolerates one edge.
- **Other signals:** BeiDou D1 has NH20 (1 kHz flips), so non-coherent 1 ms, or 20 NH-phase hypotheses. BDS GEO D2
  allows ≤ 2 ms. Galileo E1-B allows ≤ 4 ms. GLONASS allows ≤ 10 ms (meander).

### 7.4 Sensitivity versus cost (GPS C/A, ±6 kHz search, 2046 code cells, P_FA = 1e-3 per PRN, P_D = 0.9, 2 dB implementation loss)

[derived: noncentral χ²; noise per cell ~χ²₂M]

| T_coh (ms) | M (non-coh) | data used (ms) | Doppler bins | IFFTs per PRN | required C/N0 (dB-Hz) | threshold (peak / mean noise) | 32 PRNs, numpy f64 (ms) N=1024 / 2048 |
|---|---|---|---|---|---|---|---|
| 1 | 4 | 4 | 25 | 100 | 41.0 | 6.5 | 15 / 30 |
| 1 | 8 | 8 | 25 | 200 | 38.7 | 4.3 | 30 / 61 |
| 1 | 16 | 16 | 25 | 400 | 36.5 | 3.0 | 60 / 122 |
| 1 | 34 | 34 | 25 | 850 | 34.3 | 2.2 | 128 / 258 |
| 2 | 17 | 34 | 49 | 833 | 33.4 | 3.0 | 125 / 253 |
| **4** | **4** | **16** | 97 | 388 | **35.3** | 6.8 | 58 / 118 |
| 4 | 8 | 32 | 97 | 776 | 32.9 | 4.5 | 117 / 236 |
| 5 | 6 | 30 | 121 | 726 | 32.9 | 5.3 | 109 / 221 |
| 10 | 1 | 10 | 241 | 241 | 36.5 | 20.0 | 36 / 73 |
| 10 | 3 | 30 | 241 | 723 | 32.4 | 8.6 | 109 / 220 |

How to read the table:
- **Timings** are measured on this Mac (Apple M1, numpy pocketfft, complex128, one thread): 4.6 µs per 1024-point and
  9.5 µs per 2048-point IFFT. C++ pocketfft in complex64, multithreaded (`SPECTRUM_INTEL_THREADS` = 4), should be
  about 3–6× faster [unconfirmed until measured]. The fold and forward FFTs add about 10–20 %.
- **For a fixed IFFT budget, longer coherent integration wins.** Compare 4 × 4 ms at 388 IFFTs (35.3 dB-Hz) with
  16 × 1 ms at 400 IFFTs (36.5 dB-Hz). Coherent-4 ms is the sweet spot for C/A: low bit-edge risk and 125 Hz bins.
- **What the Pluto can see.** With a passive antenna and no LNA, a Pluto probably sees strong satellites at about
  38–44 dB-Hz [unconfirmed: NF of the AD9363 front end at 1.5 GHz plus cable loss]. The "fast gate" (1 ms × 8 or
  4 ms × 4) then finds the strong satellites. The full 34 ms run finds those down to about 33 dB-Hz.
- **Pluto has no bias-tee.** An active GNSS antenna needs an external bias-tee: a bench note, not a library one.

**Recommended parameters (GPS C/A, "milliseconds" budget):**
1. **Fast gate:** 1.024 MSps (÷60), T_coh = 4 ms, M = 4 (16 ms), 125 Hz bins over ±6.6 kHz around −ε·fc, 32 PRNs.
   About 12.4k IFFT-1024, roughly 15 ms in C++ with 4 threads [estimate]. Detects ≥ 36 dB-Hz.
2. **Confirm/refine** only the hits, and only when the gate finds ≥ 1 PRN or when GNSS jamming/spoofing is being
   assessed:
   - 2.048 MSps, T_coh = 4–5 ms, M = 6–8 over the whole block, Doppler ±250 Hz around the hit (5 bins), 3 code
     phases.
   - Then a "full sky" run for the missed PRNs, using the common clock offset learned from the hits to narrow the
     Doppler to ±5 kHz.
3. **Mandatory gate before any of this:** the block must contain 1575.42 ± 1.1 MHz inside the usable band.
   Otherwise skip GNSS entirely. This is a cheap trigger, per the latency rule.

### 7.5 Detection threshold and false alarms

- **Noise statistic.** For a cell without signal, 2·Z/σ² ~ χ²₂M, where Z = Σ|R|² and σ² is the noise power per
  complex correlator output.
- **Estimate σ²** robustly from the whole search grid of that PRN, excluding ±2 chips × ±1 bin around the peak.
  Use the mean, or the median / ln 2 for M = 1 (exponential).
- **Threshold:** γ = χ²₂M⁻¹(1 − P_fa,cell)/(2M) on Z/σ², with P_fa,cell = P_fa,sys / N_cells. N_cells = code
  phases × bins, about 5·10⁴–2·10⁵. Use P_fa,sys = 10⁻³ per PRN per block, so about 3 % of blocks have a false PRN
  over 32 PRNs. Report a satellite only above γ. For "GNSS present", require **≥ 2 PRNs**, or 1 PRN with C/N0 ≥ 40.
  That makes false labels about 10⁻⁴ per block.
- **Second criterion** (used by SoftGNSS-style receivers): peak / second-highest peak outside ±1 chip, in the same
  Doppler bin, ≥ 2.0–2.5 [lit: Borre et al. 2007 use 2.5]. It is less principled than the P_fa threshold. Keep it as
  a sanity check, and as the multipeak detector for spoofing (§9).
- **Cross-correlation guard.** A strong PRN (≥ 50 dB-Hz, or a spoofer) produces cross-correlation peaks at −24 dB in
  other PRNs. Reject a weak hit whose (code phase, Doppler mod 1 kHz) coincides with a strong PRN's, and whose level
  is ≤ (strong − 21 dB).

### 7.6 C/N0 from the correlation peak

Let R be the complex correlator over T_coh, signal amplitude A, noise E|n|² = σ². Then SNR_post = |A|²/σ² =
(C/N0)·T_coh [lit: Kaplan & Hegarty; Van Diggelen]. With non-coherent sums (M segments), Z̄ = Z/M:

```
SNR_post = (Z̄_peak − σ̂²) / σ̂²            # subtract noise bias of |R|^2
C/N0 [dB-Hz] = 10·log10(SNR_post / T_coh) + L_corr
L_corr = straddle loss (code, from parabolic interpolation of 3 code cells) + Doppler scalloping (sinc² over the bin offset) + front-end filter loss (~0.5–1 dB at ±1 MHz)
```

Notes:
- Interpolate the peak (parabola on |R| at −1/0/+1 code cell and on the Doppler axis) to estimate and remove the
  straddle and scalloping losses. Otherwise C/N0 is biased low by up to 3 dB.
- The estimate's standard deviation is about 1/√M of SNR_post in linear terms at high SNR. With M = 8 expect about
  ±1 dB.
- Below about 33 dB-Hz the bias correction dominates: report "< 35 dB-Hz" rather than a number.

### 7.7 Output per block (proposal, additive to the response schema)

`gnss: {system: "GPS L1 C/A", prns: [{prn, cn0_dbhz, doppler_hz (after −ε·fc), code_phase_chips, peak_ratio}],
common_clock_offset_hz, noise_floor_dbfs_hz, in_band_power_rise_db, indicators: {...}}`

Labels:
- "GNSS: GPS L1 C/A (n satellites)";
- also "GNSS: BeiDou B1I (n)", "GNSS: GLONASS L1 (n channels)", "GNSS: NavIC L5 (n)";
- "GNSS jamming-like";
- "GNSS spoofing-like".

All carry `detection` + `decoded.steps` provenance. The label is research-tagged and "-like" except the satellite list.
A satellite list confirmed by ≥ 2 PRNs above the P_fa threshold is a verified detection, not an inference.

---------------------------------------------------------------------------------------------------------------------

## 8. Jamming indicators

1. **In-band power / floor rise (AGC-equivalent).** The library runs at fixed gain, so it has no AGC. Use:
   - J/N = P(1575.42 ± 1.023 MHz) / P_ref, where P_ref is a calibrated thermal floor in dBFS/Hz for the gain setting,
     or the median floor of quiet parts of the same block outside 1559–1610 MHz.
   - Authentic GNSS adds < 0.1 dB, since the signals are 20–30 dB below the noise in 2 MHz.
   - **J/N ≥ 3 dB** is suspicious and **≥ 6 dB** is interference [heuristic: the C/N0 loss is already about 1–3 dB
     there].
   - Effective C/N0 under jamming [lit: Kaplan & Hegarty]: (C/N0)_eff = [1/(C/N0) + (J/C)/(Q·R_c)]⁻¹, with
     R_c = 1.023e6 and Q ≈ 1 (narrowband at a spectral line) / 1.5 (matched-spectrum) / 2 (wideband).
   - Watch out for the Pluto's own self-noise: the floor rises 4–5 dB about 12 ms into each 2^21 buffer (CLAUDE.md).
     Measure the floor rise per 1 ms segment and compare against the block's own early segments, or the jammer flag
     will trip on our own USB noise.
2. **Chirp / swept jammers ("personal privacy devices").** Mitch et al. 2011 measured 18 commercial jammers:
   - every one is a swept tone (linear FM);
   - L1 sweep period 1–27 µs, most of them **9 µs**: 12 of 18 at 9 µs, 2 at 8 µs, 1 at 12 µs, 2 at 26–27 µs, 1 at 1 µs;
   - sweep extent about 4–31 MHz below/above L1, typically < 20 MHz total;
   - some are off-centre by several MHz; one moved > 2.5 MHz in 41 µs;
   - all jammed L1, 6 of 18 also L2, none L5 [lit, Table 1].

   Detection:
   - The library's `jammer` receiver detects swept carriers with sweep rates between max(300 Hz, 8 sweeps/block) and
     row rate/8 = 120 kHz (STFT fs/64 bins, 1.04 µs rows). A 9 µs PPD (111 kHz) is at that limit, and 1 µs (1 MHz)
     is outside it.
   - **Gap to close:** for GNSS-band PPDs, add a short-STFT variant: fs/16 bins = 3.84 MHz bins, 0.26 µs rows, up to
     about 500 kHz sweep rate. Or use instantaneous frequency (phase difference of the 61.44 MSps samples): a PPD is
     typically 20–60 dB above the floor near the device, so IF is clean.
   - Report "GNSS jamming-like (swept, period P µs, span lo–hi MHz)" when the swept band covers 1575.42 MHz (or
     1561.098 / 1602 / 1176.45).
3. **CW / narrowband.** A spectral line within ±1.023 MHz of L1 with J/N_line ≫ 0 dB. The C/A spectrum has lines every
   1 kHz, so a CW tone close to a line is worst (Q ≈ 1).
   - **Exclude own spurs first.** The Pluto makes integer-MHz reference lines, so **1575.000 MHz sits 420 kHz from L1
     inside the main lobe**. It also makes 33.333 MHz multiples (1566.67 and 1600.00 MHz, the latter inside the
     GLONASS band) and k × 61.44 MHz lines. The library already labels these receiver spurs.
   - A CW is reported as jamming-like only if it is not one of those, or if it is ≥ 20 dB stronger than the
     calibrated spur level.
4. **Loss of all satellites.** Zero PRNs with a raised floor (J/N ≥ 6 dB) or a swept/CW interferer in band gives
   "GNSS jamming-like". Zero PRNs alone (antenna indoors, no sky) gives **no** label: never guess silently. With
   multi-block context (`StreamAnalyzer`, or signal-vault replays), PRNs present then absent with a floor rise is
   strong evidence.
5. **Pulsed interference in L5/NavIC L5** (DME/TACAN, radar). Report it as its own signal, not as GNSS jamming. Blank
   pulses before acquisition.

---------------------------------------------------------------------------------------------------------------------

## 9. Spoofing indicators usable from one snapshot

Psiaki & Humphreys (2016) group defences into signal-processing (power, distortion, multiple peaks), drift/
consistency and cryptographic checks. Only the first group and some consistency checks fit a single block.

| # | Indicator | Computation | Threshold (source) | Notes |
|---|---|---|---|---|
| S1 | **Absolute C/N0 too high** | max C/N0 over PRNs | > 55 dB-Hz is implausible for a 0 dBi antenna: −150 dBW max + 204 = 54 dB-Hz with 0 dB NF [derived from IS-GPS-200 max]. Calibrate per setup: with a passive antenna on a Pluto, > 50 dB-Hz is already suspicious [unconfirmed] | A spoofer must out-power the true signals (typically by ≥ 3 dB, often 10+ dB) to capture tracking loops [lit: Humphreys 2008; Psiaki & Humphreys 2016] |
| S2 | **Uniform C/N0** | standard deviation of C/N0 over ≥ 4 PRNs | σ < 1.5–2 dB with mean > 45 dB-Hz [heuristic, unconfirmed]. Authentic skies span about 10–15 dB (elevation, antenna gain) [lit: Jafarnia-Jahromi 2012 §C/N0 monitoring] | A single-antenna spoofer gives near-equal powers unless it simulates elevation (gps-sdr-sim does vary power with elevation — test against it) |
| S3 | **In-band power rise with good tracking** | the J/N of §8.1 with PRNs at high C/N0 and no visible jammer | J/N > 1–2 dB with ≥ 4 PRNs above 45 dB-Hz [lit: Akos 2012 (AGC); Wesson et al. 2018 (power-distortion detector)] | Thermal noise dominates authentic GNSS, so the total power is unchanged by real satellites. A spoofer adds power. Strongest single-snapshot indicator |
| S4 | **Multiple peaks per PRN** | second peak above threshold separated by > 1.5 chips or > 1 Doppler bin from the first, after the cross-correlation guard (§7.5) | in ≥ 2 PRNs [lit: Psiaki & Humphreys 2016 — signal-processing defences; Jafarnia-Jahromi 2012] | Authentic + spoofed (or meaconed) copies. A meaconer's delay is ≥ µs-scale, so the code offset is several chips or more. Multipath gives < 1 chip offsets and much weaker copies |
| S5 | **Doppler inconsistency** | Dopplers after the −ε·fc correction | any \|f_D\| > 5.5 kHz for a static receiver, or ≥ 4 PRNs within ±100 Hz of each other [heuristic] | Weak without position/time. Optional upgrade: an almanac + approximate time + India location predicts visible PRNs and Dopplers (± a few hundred Hz). Large mismatches flag spoofing. Needs almanac input (out of single-block scope) |
| S6 | **Visible-PRN set inconsistency** | PRN count and identities | > 12–14 GPS PRNs at once is implausible (typical 8–12 above 5°), and so is a set that cannot co-exist; needs an almanac [lit: standard] | Cheap check: count only |
| S7 | **Single-constellation sky** | in one block centred at about 1580 MHz (covering GPS, Galileo E1, BeiDou B1I, GLONASS), GPS strong (≥ 4 PRNs > 45 dB-Hz) and **no** BeiDou B1I / GLONASS | flag as supporting evidence [heuristic] | In India 6–10 BeiDou and 5–8 GLONASS satellites are normally visible. Common low-cost spoofers (HackRF + gps-sdr-sim) transmit GPS L1 C/A only. Needs B1I + GLONASS acquisition (§4–5); GLONASS is FDMA, so cheap |
| S8 | **SQM / correlation-peak distortion** | at ≥ 8 samples/chip (8.192 MSps), the ratios Δ = (E−L)/P and R = (E+L)/(2P) at ±0.5 chip | deviation from a clean-capture baseline by > 3σ [lit: Pini et al.; Manfredini et al. 2014; Wesson et al. 2018]. No universal number: calibrate on our own clean captures | Only at high C/N0 (≥ 45 dB-Hz) with ≥ 20 ms integration; also triggers on multipath |
| S9 | **Code-phase / clock pattern** | all code phases mod 1 ms vs predicted (needs time + ephemeris) | — | Not usable from one snapshot without a time reference. Listed for completeness |

**Label rule (proposal):** "GNSS spoofing-like" when **S3** holds, or when ≥ 2 of {S1, S2, S4, S7} hold, with ≥ 4 PRNs
acquired. Confidence 0.6, raised to 0.75 with three indicators. The evidence lists the indicators and their values.
Never label on S5/S6 alone.

---------------------------------------------------------------------------------------------------------------------

## 10. Synthetic test signals (one 2^21 block)

Per satellite i, with t = n/fs (use float64 and integer sample indices; never accumulate t):

```
f_c_i   = f_IF_i + f_D_i - eps*f_RF              # carrier at baseband: IF from tuning offset, Doppler, receiver LO error
code_phase_chips_i(t) = tau0_i + R_c*(1 + f_D_i/f_RF) * t_true(t),   t_true = t*(1+eps)   # receiver clock fast by eps
s_i(t)  = A_i * d_i(t) * c_i[floor(code_phase) mod 1023] * exp(j*(2*pi*f_c_i*t + phi_i))
d_i(t)  = ±1 nav bits, 20 ms, edges at code epochs; random first-bit offset
A_i     = sqrt(10**(CN0_i/10) * sigma2 / fs)     # complex noise variance sigma2 per sample => N0 = sigma2/fs
x       = sum_i s_i + sqrt(sigma2/2)*(randn + j*randn)
```

- Use fractional chip positions for exactness. Sample the code as rectangular chips, or band-limit it with the
  same FIR as the receiver to test straddle losses. At 61.44 MSps a chip is 60.06 samples.
- **Scale:** at 45 dB-Hz and 61.44 MSps the per-sample SNR is −32.9 dB; at 35 dB-Hz it is −42.9 dB.
- **Store as SigMF ci16** like the testdata, with σ at about 300–1000 LSB so the 12-bit quantisation is not the
  limit.
- **Cases to include:**
  - 8 PRNs with C/N0 from 32 to 48 dB-Hz, Doppler ±4.5 kHz, eps = +9.5 ppm (receiver_ppm given and withheld);
  - a bit edge inside a 4 ms fold;
  - a 34/37 duplicate;
  - a strong 55 dB-Hz PRN that leaves −24 dB cross-correlation;
  - a spoofed sky: all PRNs at 50 ± 0.5 dB-Hz, +2 dB in-band power, the authentic sky still present 1.5–10 chips
    away (multipeak);
  - a PPD chirp: 9 µs sawtooth over L1 ± 10 MHz at J/N = 30 dB, and 1 µs;
  - CW at L1 + 3 kHz at J/N = 20 dB;
  - an own-spur line at 1575.000 MHz (must not label);
  - noise only (must not label).
- **Oracles:**
  - **gps-sdr-sim (MIT)** generates a realistic GPS L1 C/A sky from a RINEX ephemeris at a chosen location/time
    (default 2.6 MSps, 8/16-bit IQ). The Doppler/power/code-phase geometry is realistic, which makes it the best
    oracle for S2/S5.
  - **GNSS-SDR (GPL-3)** is an acquisition/tracking oracle only, run as an external program on our synthetic files.
  - **MATLAB Satellite Communications Toolbox** has NavIC and GPS waveform generators (proprietary; oracle only).
  - **Never transmit:** keep it IQ-file only. Radiating GNSS signals is illegal without authorisation (in India under
    the WPC/DoT).

---------------------------------------------------------------------------------------------------------------------

## 11. Implementation recipe (receiver `gnss/`, Python first)

0. **Trigger** (µs): 1575.42 ± 1.1 MHz (or 1561.098 ± 2.1, 1598–1606, 1176.45 ± 1.1/10.3) is inside the block's usable
   band, fs ≥ 2.1 MSps, and the block is long enough: ≥ 4 ms for the fast gate. Otherwise skip.
1. **Interference pre-check** (reuse the block FFT): in-band J/N per 1 ms segment against the floor; CW lines
   (excluding the known own-spur list); hand off to the `jammer` sweep detector (extended range).
2. **Decimate** by cropping the existing block FFT to ±1.024 MHz (C/A) and IFFT, which gives 2.048 MSps. A further
   ÷2 gives 1.024 MSps for the gate. Apply −ε·fc in the same step (bin shift + residual phasor).
3. **Fast gate:** T_coh = 4 ms, M = 4, ±(5 + 1.6·Δppm) kHz, 125 Hz bins, 32 PRNs, threshold at P_fa = 10⁻³ per PRN.
   Replica FFTs cached per (PRN, fs_eff), like the existing "reuse what the pipeline computed" caches.
4. **Refine** each hit: 2.048 MSps, the whole block, parabolic peak interpolation → code phase, Doppler, C/N0 (§7.6).
   Search the rest of the sky using the common offset.
5. **Optional:** BeiDou B1I (2.046 Mcps code; 4.096 MSps crop at 1561.098) and GLONASS (14 channels, 511 chips), for
   the S7 indicator and the multi-GNSS label.
6. **Indicators** (§8, §9) → labels: "GNSS: GPS L1 C/A (n satellites)", "GNSS jamming-like", "GNSS spoofing-like".
   Below threshold the output is `unknown` with the evidence attached.
7. **Tests:**
   - code self-tests: the C/A octal column for 1–63, L5 initial states, NavIC octals;
   - synthetic scenarios of §10;
   - bench: a live L1 capture through an active antenna and bias-tee (testdata/live/GNSS_L1_*). Use the
     `capture_testdata.py` crop with fc = 1580 MHz for multi-GNSS.
8. **C++:** generators and acquisition are plain LFSR + pocketfft (BSD-3), with no GPL dependency. Galileo memory codes
   only after the Annex H decision.

---------------------------------------------------------------------------------------------------------------------

## 12. Open questions / unconfirmed

- The receiver-ppm sign convention on the Pluto (−ε·fc vs +ε·fc) and whether `receiver_ppm` describes the LO or the
  reference: confirm with the LTE-measured ppm on a bench L1 capture.
- GLONASS minimum power (−161 dBW) and time-mark bits: re-read ICD Ed. 5.1.
- The BDS preamble bits: re-read ICD v3.0 §5.
- The Galileo memory-code source, given that the PDF attachments are missing; the legal reading of Annex H for a
  redistributed library.
- The NavIC GSO Doppler range (expected small) and the NVS-01 L1 status today.
- Whether GNSS-SDR's L5 XB indexing is a real bug (§0). It does not affect us if we follow the ICD.
- Actual achievable C/N0 with the Pluto + antenna on the bench: it decides between the fast-gate parameters and the
  full-block search.
- Spoofing thresholds S2/S5/S7 are heuristics. Calibrate on bench captures and on gps-sdr-sim/TEXBAT-like data.
  TEXBAT recordings are public from UT Austin; their licence is [unconfirmed].

---------------------------------------------------------------------------------------------------------------------

## 13. Sources

| # | Source | URL | Licence / terms |
|---|---|---|---|
| 1 | IS-GPS-200N, Navstar GPS Space Segment/Navigation User Segment Interfaces (Tables 3-Ia, 3-Ib; §3.3.1.6 power) | https://www.gps.gov/technical/icwg/IS-GPS-200N.pdf (retrieved via web.archive.org; the gps.gov path returned 404 on 2026-10-25) | US Government work, public release, distribution unlimited |
| 2 | IS-GPS-705J, Navstar GPS Space Segment/User Segment L5 Interfaces (Table 3-I; §3.3.2.2) | https://www.gps.gov/technical/icwg/IS-GPS-705J.pdf (via web.archive.org) | US Government work, public release |
| 3 | GPS.gov, L5 PRN Code Assignments, Jan 2026 edition (SBAS I5/Q5 states) | https://www.gps.gov/sites/default/files/2026-02/L5-PRN-code-assignments-2026-Jan.pdf | US Government, public |
| 4 | ISRO, IRNSS Signal-in-Space ICD for SPS, v1.1 (Aug 2017) | https://www.isro.gov.in/media_isro/pdf/SateliteNavigation/irnss_sps_icd_version1.1-2017.pdf | © ISRO; reproduction only for non-profit purposes, with the terms attached and unaltered |
| 5 | ISRO, NavIC Signal-in-Space ICD for SPS in L1 frequency, v1.0 (Aug 2023) | https://www.isro.gov.in/media_isro/pdf/SateliteNavigation/NavIC_SPS_ICD_L1_final.pdf | © ISRO; same kind of terms |
| 6 | CSNO, BDS SIS ICD Open Service Signal B1I v3.0 (Feb 2019) | http://en.beidou.gov.cn/SYSTEMS/ICD/201902/P020190227702348791891.pdf | © China Satellite Navigation Office; no explicit reuse clause |
| 7 | EU, Galileo OS SIS ICD Issue 2.1 (Nov 2023), Annex C (codes), Annex H (IPR authorisation) | https://www.gsc-europa.eu/sites/default/files/sites/all/files/Galileo_OS_SIS_ICD_v2.1.pdf | © EU 2023; royalty-free, non-exclusive, non-transferable authorisation (Annex H) |
| 8 | Russian Space Systems, GLONASS ICD, General description of FDMA signals, Edition 5.1 (2008) | https://www.glonass-iac.ru/ (ICD downloads; not retrieved this session) | © RSS; terms not checked [unconfirmed] |
| 9 | GNSS-SDR (`GPS_L5.h`, `gps_l5_signal_replica.cc`, `beidou_b1i_signal_replica.cc`, `glonass_l1_signal_replica.cc`, `Galileo_E1.h`) — read for facts and as an oracle only | https://github.com/gnss-sdr/gnss-sdr | **GPL-3.0-or-later**: do not copy into c++/ |
| 10 | gps-sdr-sim (GPS L1 C/A baseband simulator) — oracle | https://github.com/osqzss/gps-sdr-sim | MIT (GitHub licence API, checked 2026-10-25) |
| 11 | RTKLIB — not needed for acquisition; reference for nav-message/positioning if ever needed | https://github.com/tomojitakasu/RTKLIB | BSD-2-Clause plus two additional clauses (readme). rtklibexplorer fork: GitHub reports NOASSERTION, licence.txt says BSD-2 |
| 12 | K. Borre, D. Akos, N. Bertelsen, P. Rinder, S. H. Jensen, *A Software-Defined GPS and Galileo Receiver*, Birkhäuser 2007 (parallel code-phase search, peak-ratio threshold 2.5; SoftGNSS MATLAB code) | — | Book ©; SoftGNSS code GPL-2 [unconfirmed] — oracle only |
| 13 | E. D. Kaplan, C. J. Hegarty (eds.), *Understanding GPS/GNSS: Principles and Applications*, 3rd ed., Artech House 2017 (C/N0, effective C/N0 under interference, Q factor) | — | Book © |
| 14 | F. van Diggelen, *A-GPS: Assisted GPS, GNSS, and SBAS*, Artech House 2009 (coherent/non-coherent sensitivity, bin spacing) | — | Book © |
| 15 | R. H. Mitch, R. C. Dougherty, M. L. Psiaki, S. P. Powell, B. W. O'Hanlon, J. A. Bhatti, T. E. Humphreys, "Signal Characteristics of Civil GPS Jammers," ION GNSS 2011, pp. 1907–1919 | https://radionavlab.ae.utexas.edu/images/stories/files/papers/Signal%20Characteristics%20Civil%20GPS%20Jammers_Mitch.pdf | Paper © authors/ION; facts cited |
| 16 | M. L. Psiaki, T. E. Humphreys, "GNSS Spoofing and Detection," Proc. IEEE 104(6):1258–1270, 2016, doi:10.1109/JPROC.2016.2526658 | https://ieeexplore.ieee.org/document/7470255 | © IEEE; facts cited |
| 17 | T. E. Humphreys, B. M. Ledvina, M. L. Psiaki, B. W. O'Hanlon, P. M. Kintner, "Assessing the Spoofing Threat: Development of a Portable GPS Civilian Spoofer," ION GNSS 2008 | https://radionavlab.ae.utexas.edu/ | © ION |
| 18 | D. M. Akos, "Who's Afraid of the Spoofer? GPS/GNSS Spoofing Detection via Automatic Gain Control (AGC)," NAVIGATION 59(4):281–290, 2012, doi:10.1002/navi.19 | https://doi.org/10.1002/navi.19 | © ION |
| 19 | K. D. Wesson, J. N. Gross, T. E. Humphreys, B. L. Evans, "GNSS Signal Authentication via Power and Distortion Monitoring," IEEE Trans. AES 54(2):739–754, 2018 | https://doi.org/10.1109/TAES.2017.2765258 | © IEEE |
| 20 | A. Jafarnia-Jahromi, A. Broumandan, J. Nielsen, G. Lachapelle, "GPS Vulnerability to Spoofing Threats and a Review of Antispoofing Techniques," Int. J. Navig. Observ. 2012, 127072 | https://doi.org/10.1155/2012/127072 | Open access, Creative Commons Attribution (version [unconfirmed]) |
| 21 | T. E. Humphreys, J. A. Bhatti, D. Shepard, K. Wesson, "The Texas Spoofing Test Battery (TEXBAT)," ION GNSS 2012 | https://radionavlab.ae.utexas.edu/texbat | data terms [unconfirmed] |
| 22 | E. Manfredini, F. Dovis, B. Motella, "Validation of a signal quality monitoring technique over a set of spoofed scenarios," NAVITEC 2014 (SQM metrics) | https://doi.org/10.1109/NAVITEC.2014.7045150 | © IEEE |
| 23 | IISc (P. V. Kumar group), "The IZ4 spreading code family incorporated into the NavIC L1 SPS signal" | https://ece.iisc.ac.in/~pvkece/NavIC.html | web page, facts |
| 24 | MathWorks, NavIC waveform generation example (oracle) | https://www.mathworks.com/help/satcom/ug/navic_waveform_generation.html | proprietary; oracle only |
| 25 | This repo: `python/spectrum_intel/receivers/jammer/README.md` (swept-carrier detector limits), `CLAUDE.md` (Pluto spurs, self-noise, ppm) | local | project |
