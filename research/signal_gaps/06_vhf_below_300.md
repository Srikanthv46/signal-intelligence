# 06 — VHF below 300 MHz (70–300 MHz): band plan, airband AM, ACARS, VDL Mode 2, AIS, marine VHF + DSC, APRS, FM + RDS, VOR / ILS localizer

Research file for PLAN.md items S13 (VHF part), S14–S21. Written 2026-10-25. Nothing here is implemented yet; all
labels produced from this research start as `research` and "-like" unless a CRC / parity / FEC check verifies them.

Conventions used below:
- **[conf]** = confirmed from a primary standard, or from two independent sources that agree, or from source code
  that is known to interoperate with real traffic (read for facts only; see licences in the Sources list).
- **[unconf]** = not confirmed; secondary source only or from memory. Verify before relying on it.
- **[computed]** = test vector computed for this file by
  `scratchpad/vec/vectors.py` / `dsc.py` (not committed; the code is ~80 lines and is reproduced in prose here,
  re-derive in the receiver's test). Every CRC vector was self-checked by recomputing the residue.
- Bit numbering: "MSB-first" = the bit order of the message as printed in the standard; "on air" = the order of
  transmission.

Block arithmetic used throughout (input is 2^21 complex samples):

| Sample rate | 2^21 samples last | One 30 Hz period (33.3 ms) | One AIS slot (26.67 ms) | One RDS group (87.6 ms) |
|---|---|---|---|---|
| 61.44 MSps | 34.13 ms | 1.02 periods | fits (if aligned) | no (0.39 group) |
| 30.72 MSps | 68.3 ms | 2.05 | yes | no (0.78) |
| 2.4 MSps | 873.8 ms | 26 | many | ~10 groups |
| 1.92 MSps | 1092 ms | 32.8 | many | ~12 groups |

At 61.44 MSps the VHF signals here are 5–200 kHz wide, so every receiver should first **channelise** (the existing
`Channelizer` / block FFT) to a narrow rate (e.g. 48–250 kSps). The narrow channel inherits the block duration
(34 ms), so duration — not bandwidth — is the binding limit; a `StreamAnalyzer` dwell is required for anything
longer than ~30 ms.

---

## 1. India band plan 70–300 MHz

Basis: ITU Radio Regulations Art. 5 (Region 3) as adopted by India's National Frequency Allocation Plan 2022
(NFAP-2022, WPC/DoT). The NFAP-2022 IND footnotes for the government land-mobile bands were **not** retrieved
(the NFAP PDF and the 2024 revision minutes on dot.gov.in were not readable here) → those rows are [unconf]. The
aeronautical, maritime, met-sat and broadcast rows follow ICAO / ITU allocations that India uses and are [conf]
at allocation level.

Proposed rows for `core/bands_india.py` (same tuple shape: lo, hi, key, description, expected technologies):

| Lo (MHz) | Hi (MHz) | key | Use in India | Expected signals | Status |
|---|---|---|---|---|---|
| 70 | 74.8 | govt_mobile_70 | Fixed / mobile (government, defence; some legacy police links) | narrowband_fm, pmr | [unconf] |
| 74.8 | 75.2 | marker_75 | Aeronautical marker beacons 75 MHz (ILS OM/MM/IM; many decommissioned) | marker_beacon | alloc [conf], use [unconf] |
| 75.2 | 87.5 | govt_mobile_75 | Fixed / mobile (government); Band II-low not broadcast in India | narrowband_fm, pmr | [unconf] |
| 87.5 | 108 | fm_bcast | FM broadcasting: AIR (Prasar Bharati) + private FM (Phase I–III); stations seen 88.1–107.9 | fm_broadcast, rds | [conf] |
| 108 | 111.975 | vor_ils_loc | VOR (even tenths +/-.05) and ILS localizer (odd tenths and odd tenths +.05) | vor, ils_loc | [conf] |
| 111.975 | 117.975 | vor | VOR (and GBAS VDB 108–117.975) | vor, gbas_vdb | [conf] |
| 117.975 | 137 | air_vhf | Aeronautical mobile (R): ATC AM voice 118–136.975; 121.5 emergency; 123.1 SAR; ACARS 131.x; VDL2 136.725–136.975 | am_voice, acars, vdl2 | [conf] |
| 137 | 138 | metsat_137 | Met-sat / space ops downlinks (NOAA APT 137.1/137.62/137.9125 — NOAA-15/18/19 retired 2025 [unconf]; Meteor-M LRPT 137.1/137.9) | apt, lrpt | [conf] alloc |
| 138 | 144 | govt_mobile_138 | Fixed / mobile (defence) | narrowband_fm, pmr | [unconf] |
| 144 | 146 | ham_2m | Amateur 2 m (India: 144–146). APRS (see §7), FM repeaters, satellites 145.8–146 | aprs, narrowband_fm | [conf] band; APRS freq [unconf] |
| 146 | 156 | govt_mobile_146 | Land mobile: police / state govt / railways / PSU VHF (12.5/25 kHz FM, some DMR) | narrowband_fm, dmr, pmr | [unconf] |
| 156 | 162.05 | marine_vhf | Maritime mobile (RR Appendix 18): ch 16 156.800, DSC ch 70 156.525, AIS 161.975/162.025 | marine_fm, dsc, ais | [conf] |
| 162.05 | 174 | govt_mobile_162 | Land mobile (police, railways, forest, PSU) | narrowband_fm, dmr, pmr | [unconf] |
| 174 | 230 | vhf_band3 | Broadcasting Band III. Doordarshan analog VHF TV transmitters largely shut 2017–2023 [unconf]; no DAB service in India [unconf]; IAFI proposal (2024) to use 174–230 for audio / wireless mics [conf: WG-1 minutes] | analog_tv, wireless_mic | [unconf] |
| 230 | 235 | govt_230 | Fixed / mobile (defence) | — | [unconf] |
| 225 | 400 | mil_air_uhf | Military aeronautical UHF AM (225–400, 243.0 guard); overlaps the 230–300 rows; satcom 240–270 / 290–320 (UFO/MUOS-class) | am_voice, satcom | ITU [conf]; India use [unconf] |
| 235 | 267 | govt_235 | Fixed / mobile / mobile-satellite (defence) | — | [unconf] |
| 267 | 300 | govt_267 | Fixed / mobile (defence); 272–273 space ops | — | [unconf] |

Notes
- The existing table starts at 300 MHz with `govt_mobile_300`; the new rows join it seamlessly at 300.
- Overlap (225–400 military air vs. 230–300 govt rows) is intentional: the military air row is an expectation for
  AM voice, not an exclusive allocation. The table already contains an overlap (865–867 inside 851–894).
- Pluto with the AD9364 mod reaches 70 MHz but the RF front end (balun / matching) is specified 325 MHz–3.8 GHz;
  expect a few dB to >10 dB loss below 200 MHz [unconf]. Strong FM broadcast (−20 dBm near a tower) will
  intermodulate in a wide-open front end: FM-band images/intermod products in 108–137 MHz must not become
  "airband AM" detections (gate on AM vs. FM modulation type).
- Mixer harmonic images (bench/live_lo/FINDINGS.md) also apply: at LO ≈ 100 MHz the 3rd harmonic (300 MHz) and
  5th (500 MHz) fold UHF signals into the VHF block. Same harmonic map logic applies.

---

## 2. Civil air band AM voice (118–136.975 MHz)

| Item | Value | Status |
|---|---|---|
| Band | 117.975–137.000 MHz (voice 118.000–136.975) | [conf] Annex 10 Vol V |
| Modulation | A3E DSB full-carrier AM, audio 300–2500 Hz (25 kHz channels), 350–2500 Hz (8.33 kHz) | [conf] |
| Channel spacing India | 25 kHz (8.33 kHz is mandatory in the ICAO EUR region above FL195; not adopted in India as of the sources seen) | [unconf] |
| Carrier offset (climax) | multi-carrier offset systems ±2.5 / ±5 / ±7.5 kHz in 25 kHz channels | [conf] Annex 10 Vol III |
| Emergency | 121.500 MHz; SAR on-scene 123.100; air-to-air 123.450 (Annex 10 Vol V) | [conf] |
| ELT 121.5 sweep tone | swept AM tone 1600→300 Hz (≥700 Hz range), 2–4 sweeps/s (legacy ELT homing signal) | [conf] |
| Military UHF AM | 225–400 MHz, 25 kHz, guard 243.000 | [conf] |

**8.33 kHz channel naming** [conf] (EUROCONTROL 8.33 kHz channel naming, ICAO Annex 10 Vol V Table 4-1): each
25 kHz block starting at frequency F (F ends in .x00, .x25, .x50, .x75) carries either one 25 kHz channel named F,
or three 8.33 kHz channels at F, F+8.333 kHz, F+16.667 kHz named F+0.005, F+0.010, F+0.015 MHz. Names therefore
end in 00/25/50/75 for 25 kHz and 05/10/15/30/35/40/55/60/65/80/85/90 for 8.33 kHz. Example: 118.0083 MHz is
named "118.010"; 132.0167 MHz is "132.015".

```
def name_to_freq(name_mhz):            # e.g. 118.010 -> 118.008333
    k = round(name_mhz * 1000)          # kHz
    blk = (k // 25) * 25; r = k - blk   # r in {0,5,10,15}
    return blk if r == 0 else blk + {5: 0, 10: 8.3333, 15: 16.6667}[r]   # kHz
```

**India ATC examples** (only public, secondary sources; check AIP India ENR / AD 2.18 at aim-india.aai.aero):
- Mumbai VABB: Tower 118.100, Ground 121.900, Approach 127.900, ATIS 126.400 [unconf, two secondary sources agree].
- Delhi VIDP: sources disagree (Tower 118.100 / 118.250 / 118.750 per runway vs. 119.300); ATIS 126.400 listed
  [unconf].
- Chennai VOMM: not found [unconf].
- ATIS is a continuous AM voice loop; VOLMET similar. Ground stations transmit continuously → a stable AM carrier
  at a 25 kHz-raster frequency in 118–137 MHz with voice-band envelope is "ATIS/VOLMET-like".

**Recipe (AM voice detector, S13)**: channel at 25 kHz raster → measure carrier (narrow line in the block FFT at a
raster frequency, ±10 kHz offset tolerated for ppm and climax) → envelope `|x|` → modulation depth m = (max−min)/
(max+min) on a 20 ms smoothed envelope; AM if a stable carrier is present and the envelope spectrum has energy in
300–3000 Hz while the instantaneous frequency (FM discriminator) is flat (FM/AM discrimination: ratio of envelope
variance to frequency variance). Push-to-talk aircraft bursts are seconds long: one 34 ms block sees a carrier +
a fraction of a syllable; label "am_voice-like" with the raster channel name. Synthetic: carrier + m·s(t) with s =
band-limited noise or a recorded speech clip, m = 0.3–0.9, add ±5 kHz offset and CFO.

---

## 3. ACARS (Plain Old ACARS, ARINC 618)

| Item | Value | Status |
|---|---|---|
| India frequencies | 131.725 MHz and 136.975 MHz assigned to SITA (28 airports, licence renewed to 2021); 131.825 MHz to BCS (4 airports, not renewed after 2014) — DoT statement in TRAI consultation paper 10-12-2022 | [conf] for assignment; current use [unconf] |
| Other common channels | 131.550 (worldwide primary), 131.525, 131.725 (EU SITA), 131.450 (Japan), 130.025/130.450/131.125/136.900 (US) | [conf] secondary |
| Note | 136.975 is also the VDL2 common signalling channel; in India it is probably VDL2 (D8PSK), not ACARS — detect by modulation | [unconf] |
| Channel | 25 kHz AM (A2D): DSB AM carrier modulated by an audio MSK signal | [conf] |
| Audio modulation | MSK, 2400 bit/s, tones 1200 Hz and 2400 Hz (h = 0.5 at 2400 bps: 1200 Hz = ½ cycle/bit, 2400 Hz = 1 cycle/bit; centre 1800 Hz) | [conf] |
| Bit encoding | Frequency encodes the relation between successive bits (MSK = OQPSK with half-sine; acarsdec demodulates coherently at 1800 Hz, alternates I/Q decisions and flips sign every 2 bits — so the decoded I/Q signs are the data). Mapping "tone ↔ bit change" to be verified by loop-back through acarsdec | [unconf] mapping |
| Character | 7-bit ASCII + **odd** parity in bit 8; transmitted **LSB first** (acarsdec shifts bits in from the top: `outbits >>= 1; if bit: outbits |= 0x80`) | [conf] |

**Frame** (characters, each 8 bits with odd parity) [conf: ARINC 618 via acarsdec/airframes docs]:

| # | Field | Chars | Value (hex incl. parity) | Notes |
|---|---|---|---|---|
| 1 | Pre-key | ≥ 16 bits (often 128 bit times) of '1' / carrier | — | length [unconf] |
| 2 | Bit sync | 2 | `+` `*` = 0xAB 0x2A | The prompt's "++ / **" is wrong: it is "+" then "*" (0x2B, 0x2A; 0x2B gets parity bit → 0xAB) |
| 3 | Char sync | 2 | SYN SYN = 0x16 0x16 | acarsdec syncs on SYN SYN SOH and accepts ~SYN (0xE9) for inverted polarity |
| 4 | SOH | 1 | 0x01 | start of heading |
| 5 | Mode | 1 | '2' (Category A, all GS) or '@'…']' / 'a'…'}' (Category B, specific GS) | CRC starts here |
| 6 | Address | 7 | aircraft registration, right-justified, padded with '.' (e.g. ".VT-ANA"); uplinks may carry flight-ID | India regs: VT-xxx |
| 7 | Technical ack | 1 | NAK 0x15 or the block ID being acknowledged ('0'–'9', 'A'–'Z', 'a'–'z') | |
| 8 | Label | 2 | e.g. "H1", "_d" (0x5F 0x7F general response / squitter), "Q0" (link test), "SQ" (squitter, uplink), "5Z", "80" | label tables: ARINC 620 |
| 9 | Block ID (UBI/DBI) | 1 | downlink '0'–'9'; uplink 'A'–'Z','a'–'z'; NUL for general response | |
| 10 | STX | 1 | 0x02 (absent if no text) | |
| 11 | Text | 0–220 | downlink text starts with MSN (4 chars, e.g. "M01A") + flight ID (6 chars, "AI0101") | |
| 12 | Suffix | 1 | ETX 0x83 (end) or ETB 0x97 (more blocks follow) | CRC ends here (inclusive) |
| 13 | BCS | 2 | CRC-16, low byte first | |
| 14 | BCS suffix | 1 | DEL 0x7F | |

**BCS (block check sequence)** [conf: acarsdec `syndrom.h` table starts 0x0000, 0x1189 → reflected 0x8408]:
CRC-16 with polynomial x^16+x^12+x^5+1, reflected (0x8408), init 0x0000, no final XOR, i.e. **CRC-16/KERMIT**
(check("123456789") = 0x2189). Computed over the bytes **including parity bits**, from the Mode character through
ETX/ETB inclusive. Validity: running the same CRC over [mode…ETX, BCS_lo, BCS_hi] gives 0. Odd parity on every
character is a second check; acarsdec corrects up to a few parity-flagged characters by syndrome search over the
CRC — a pure 1-bit error in character k has a known syndrome, so single-bit errors are correctable.

**Worked example** [computed; synthetic, not a recorded message]: link test from VT-ANA, no text:
characters `2 . V T - A N A NAK Q 0 1 ETX` → with odd parity
`32 AE D6 54 AD C1 CE C1 15 51 B0 31 83`, BCS = 0x9902 → transmitted `02 99`, then `7F`.
Full on-air byte sequence after pre-key: `AB 2A 16 16 01 32 AE D6 54 AD C1 CE C1 15 51 B0 31 83 02 99 7F`, each
byte LSB first. CRC over the 13 bytes + `02 99` = 0x0000.

**Durations**: 2400 bps → 3.33 ms/char. Minimum block (no text) = 2+2+1+13+2+1 = 21 chars + pre-key ≈ 70 ms +
pre-key (≥ 7–53 ms) ≈ 80–125 ms; a full 220-char block ≈ 0.8 s. **Never fits a 34 ms block**; a 0.87 s block
(2.4 MSps) holds most single blocks only if it starts before the pre-key. Needs StreamAnalyzer (≥ 1 s dwell on
131.725 / 131.825 / 131.55).

**Recipe**: channelise 25 kHz (sample 12–48 kSps) → AM envelope (remove DC/carrier) → audio band-pass 600–3000 Hz
→ coherent MSK demod at 1800 Hz with bit clock 2400 Hz (acarsdec-style PLL: VCO at 1800 Hz, matched filter
half-sine, decision every bit alternately on I and Q) or non-coherent: quadrature FM discriminator on the audio
(1200 vs 2400 Hz) + clock recovery + differential decoding → bytes LSB-first → search SYN SYN SOH (and inverted)
→ collect until ETX/ETB (abort at > 240 chars) → parity check → BCS check (KERMIT residue 0) → fields: mode, reg,
label, block ID, MSN, flight → label "acars" (verified) else "acars-like". Fields map to the response schema as
`decoded.fields.{registration, flight, label, msn, text}`.

**Synthetic**: build the byte sequence above (+ text), bits LSB first, MSK at 2400 bps between 1200/2400 Hz with
continuous phase (h = 0.5), AM-modulate a carrier at m ≈ 0.6, add noise. Oracle: acarsdec (`-f` on a WAV; its
test.wav ships in the repo) — LGPL-2.0 per the source headers ("GNU Library General Public License version 2"),
test only.

---

## 4. VDL Mode 2

Sources: ICAO Doc 9776 (Manual on VDL Mode 2), ETSI EN 301 841-1 (physical layer); facts cross-checked in dumpvdl2
(GPL-3, read only). All constants in this section are [conf] from dumpvdl2 unless marked.

| Item | Value |
|---|---|
| Frequencies | 136.975 (CSC, common signalling channel, worldwide); 136.725, 136.775, 136.825, 136.875, 136.925 (25 kHz). India: 136.975 assigned to SITA [conf DoT] |
| Access | CSMA (p-persistent), bursts |
| Modulation | D8PSK, Gray-coded, differential: phase increment Δφ = k·π/4; bits ↔ k via graycode[k] = {0:000, 1:001, 2:011, 3:010, 4:110, 5:111, 6:101, 7:100} (dumpvdl2 table `{0,1,3,2,6,7,5,4}`); raised-cosine α = 0.6 [conf EN 301 841] |
| Symbol rate | 10 500 sym/s → 31.5 kbit/s; occupied ~16 kHz in a 25 kHz channel |
| Burst | ramp-up (~5 symbols, [unconf]) → **training sequence 16 symbols** → header (25 bits) → data (AVLC frames, bit-stuffed, HDLC flags) → RS parity octets → ramp-down |
| Training sequence (16 symbols, 3 bits each) | `000 010 011 110 000 001 101 110 001 100 011 111 101 111 100 010` — equivalently cumulative phases 0, 3, −3, 1, 1, 2, 0, 4, −3, 4, −2, 3, 1, −2, −3, 0 (×π/4) [conf: derived from dumpvdl2 `pr_phase` + Gray table; matches EN 301 841] |
| Scrambler | additive, LFSR x^15 + x + 1, seeded at the first header bit with dumpvdl2's register value **0x6959** (shift: `b = (r ^ (r>>14)) & 1; r = (r>>1) | (b<<14); bit ^= b`). The seed's bit order follows that code; EN 301 841 gives it as a 15-bit vector [conf for the code form] |
| Header | 3 reserved bits (0) + **17-bit transmission length** (number of data bits before FEC, sent **LSB first**) + **5-bit header FEC**: (25,20) linear block code, parity-check rows (25-bit, MSB = first header bit): `0000000011111111111110000`, `0011111100001111111101000`, `1100011100110000111100100`, `1101101101010011001100010`, `0110100111100101010100001`; syndrome → single-error (and some double) correction table. Header bits read MSB-first from the descrambled stream |
| FEC | Reed–Solomon (255,249) over GF(2^8), field generator **0x187** (x^8+x^7+x^2+x+1), first consecutive root **α^120**, primitive element step 1, 6 parity octets (libfec `init_rs_char(8, 0x187, 120, 1, 6, 0)`). Shortened: data placed at the start of each 249-octet block, zero-filled. Last (partial) block carries 0 / 2 / 4 / 6 parity octets for last-block length < 3 / < 31 / < 68 / ≥ 68 octets; dropped parity is decoded as erasures |
| Interleaving | octet interleave across blocks: with B blocks, transmitted data octet i belongs to block (i mod B), position (i div B); parity octets interleaved the same way after all data octets |
| Bit order | octets transmitted LSB first (dumpvdl2 `bitstream_read_lsbfirst`) |
| AVLC frame | HDLC-like: flag 0x7E, bit stuffing (0 after five 1s), [dst addr 4 octets][src addr 4 octets][LCF 1][info][FCS 2], flag |
| FCS | CRC-16/X.25 (reflected 0x8408, init 0xFFFF, final XOR 0xFFFF), sent low byte first; good-frame residue (without final XOR) **0xF0B8** |
| Address (4 octets) | bit 0 of every octet = HDLC extension bit (1 only in the last octet of the address field, i.e. src octet 4). Take bits 1–7 of octets 0..3 into a 28-bit value `v = (o0>>1) | o1<<6 | o2<<13 | (o3&0xFE)<<20`, bit-reverse 28 → bits 0–23 = 24-bit address, bits 24–26 = type, bit 27 = C/R (dst) or A/G (src). Type: 1 = aircraft (ICAO 24-bit, same as Mode S/ADS-B), 4 = ground station ICAO-administered, 5 = ground station delegated, 7 = all stations (0, 2, 3, 6 reserved) |
| LCF | HDLC control: I (bit0 = 0), S (bits 1..0 = 01), U (11); XID (U) carries GS info, frequency lists, position (in XID params) |
| ACARS over AVLC | I-frame whose info starts `FF FF 01` → then an ACARS block (mode, reg, ack, label, block ID, STX, text, ETX/ETB, BCS) without the pre-key/sync |
| ATN | otherwise X.25 → CLNP → COTP/ IDRP / CPDLC (ASN.1 PER) — out of scope |

**Minimum to get AVLC addresses (the SIGINT-relevant fields: aircraft ICAO 24-bit, GS ID, direction)**:
1. preamble correlation (16 differential phases) → symbol timing + CFO (dumpvdl2 fits a line to the 16 phase
   errors; CFO up to ~±1.5 kHz at 136.975 MHz with Pluto ppm → widen search);
2. differential D8PSK demod + Gray → bits;
3. descramble from the header start; decode the 25-bit header (syndrome must be 0 or correctable) → data length;
4. read ceil(len/8) data octets; if len ≤ 249 octets there is one RS block and data octets are in order — the
   first AVLC frame starts at the first 0x7E; else de-interleave first;
5. de-stuff, check FCS (0xF0B8 residue) — this alone verifies the frame; RS is optional when FCS passes;
6. parse dst/src addresses (8 octets) + LCF.
Per-frame decode needs all octets of that frame (FCS); addresses are bytes 0–7, so even an FCS-failed frame
yields "-like" addresses.

**Durations**: overhead (ramp + 16 training + 25-bit header) ≈ 2.4 ms. A 40-octet frame (e.g. ground XID / short
I-frame) ≈ 40·8/31.5k + RS 6 octets ≈ 12 ms total → **fits a 34 ms block**. Long uplinks (ACARS 220 chars + AVLC
≈ 250 octets + 12 RS) ≈ 70 ms → needs ≥ 0.1 s. GS "squitter" XIDs repeat every ~2 min on the CSC [unconf].

**Synthetic**: build AVLC frame (addresses, LCF, info `FF FF 01` + ACARS block), FCS, flags, stuffing; prepend
nothing; data length = bit count; header = 3 zeros + length LSB-first + 5 FEC bits (generator: pick the 5 parity
bits so that H·h = 0 — H above is systematic in the last 5 columns: rows end `10000, 01000, 00100, 00010, 00001`);
scramble from header start with the LFSR; append RS parity (encode with the same field/roots) interleaved;
map 3 bits → Gray k → Δφ = kπ/4 (training first), RRC/raised-cosine α=0.6 at 10.5 ksym/s. Oracle: dumpvdl2 on an
IQ file (`--iq-file`, complex S16 or U8 at 105 kSps = 10 SPS... its sample-rate/oversample options are documented
in its README) — GPL-3, test only. The dumpvdl2 repo ships `test/vdl2_model_16b_1050kHz.wav`.

---

## 5. AIS (ITU-R M.1371-5)

| Item | Value | Status |
|---|---|---|
| Channels | AIS 1 = ch 87B **161.975 MHz**, AIS 2 = ch 88B **162.025 MHz** (25 kHz); long-range (msg 27, satellite) ch 75/76 156.775 / 156.825 | [conf] |
| Modulation | GMSK, BT 0.4 (transmit), modulation index 0.5, 9600 bit/s (±50 ppm) | [conf] |
| Access | SOTDMA / ITDMA / CSTDMA (class B "CS"), 2250 slots per minute per channel, slot = 26.667 ms = 256 bits | [conf] |
| Burst structure (one slot, 256 bits) | ramp-up 8 bits → **training 24 bits** `0101…01` (alternating, starts with 0; after NRZI → alternating tones) → **start flag 0x7E** (`01111110`) → data 168 bits (msg 1/2/3/4/18; up to 1008 bits multi-slot) → **FCS 16 bits** → end flag 0x7E → buffer 24 bits (bit stuffing 4 + distance delay 12 + repeater 2 + sync jitter 6) | [conf] |
| Line coding | **NRZI: a 0 bit = level change, 1 = no change**, applied to everything from training to end flag | [conf] |
| Bit stuffing | after five consecutive 1s in data+FCS a 0 is inserted (not in flags) | [conf] |
| Bit order | message bits are defined MSB-first in the field tables, packed into octets in that order; **each octet is transmitted LSB first** (HDLC). So on air, field bits appear bit-reversed within each octet | [conf] (AIS-catcher / gr-aistx behaviour; residue check below confirms) |
| FCS | CRC-16-CCITT x^16+x^12+x^5+1, HDLC convention = **CRC-16/X.25**: reflected (0x8408) over octets as transmitted, init 0xFFFF, final XOR 0xFFFF (ones' complement); sent low byte first, LSB first; check("123456789") = 0x906E; good-frame residue (CRC register over data+FCS, no final XOR) = **0xF0B8** | [conf] [computed] |
| MMSI | 9 digits. MID 419 = India (ships 419xxxxxx). Base station 00MIDxxxx, AtoN 99MIDxxxx, SAR aircraft 111MIDxxx, AIS-SART 970xxyyyy, MOB 972xxyyyy, EPIRB-AIS 974xxyyyy | [conf] |

**Message 1/2/3 (position report class A, 168 bits)** [conf gpsd AIVDM]:

| Bits | Len | Field | Scaling |
|---|---|---|---|
| 0–5 | 6 | type (1, 2, 3) | |
| 6–7 | 2 | repeat | |
| 8–37 | 30 | MMSI | |
| 38–41 | 4 | nav status | 0 under way engine, 1 at anchor, 2 not under command, 3 restricted manoeuvrability, 4 constrained by draught, 5 moored, 6 aground, 7 fishing, 8 sailing, 9–10 reserved (HSC/WIG), 11 towing astern, 12 pushing/towing alongside, 14 AIS-SART/MOB/EPIRB active, 15 not defined |
| 42–49 | 8 | ROT (signed) | ROT°/min = sign·(raw/4.733)^2; ±127 = turning > 5°/30 s without TI; −128 = n/a |
| 50–59 | 10 | SOG | 0.1 kn; 1023 n/a; 1022 = ≥102.2 kn |
| 60 | 1 | position accuracy | 1 = < 10 m (DGNSS) |
| 61–88 | 28 | longitude (signed) | 1/10000 min → deg = raw/600000; 181° (0x6791AC0) n/a |
| 89–115 | 27 | latitude (signed) | raw/600000; 91° (0x3412140) n/a |
| 116–127 | 12 | COG | 0.1°; 3600 n/a |
| 128–136 | 9 | true heading | 511 n/a |
| 137–142 | 6 | UTC second | 60 n/a, 61 manual, 62 dead reckoning, 63 inoperative |
| 143–144 | 2 | manoeuvre indicator | |
| 145–147 | 3 | spare | |
| 148 | 1 | RAIM | |
| 149–167 | 19 | radio status (SOTDMA/ITDMA comm state) | |

**Message 4 (base station report, 168)**: 0–5 type, 6–7 repeat, 8–37 MMSI, 38–51 year (14), 52–55 month, 56–60 day,
61–65 hour, 66–71 minute, 72–77 second, 78 accuracy, 79–106 lon (28), 107–133 lat (27), 134–137 EPFD type, 138–147
spare, 148 RAIM, 149–167 radio. (Msg 11 = same layout, UTC response.) [conf]

**Message 5 (static and voyage, 424 bits, 2 slots)** [conf gpsd]: 0–5 type, 6–7 repeat, 8–37 MMSI, 38–39 AIS
version, 40–69 IMO, 70–111 call sign (7 six-bit chars), 112–231 name (20 chars), 232–239 ship type, 240–248 to bow,
249–257 to stern, 258–263 to port, 264–269 to starboard (m), 270–273 EPFD, 274–277 ETA month, 278–282 ETA day,
283–287 ETA hour, 288–293 ETA minute, 294–301 draught (0.1 m), 302–421 destination (20 chars), 422 DTE, 423 spare.
Many transmitters send 420/422 bits. On air: 8+24+8+424+16+8 + stuffing ≈ 490 bits = 51 ms (> one 34 ms block).

**Message 18 (class B position, 168)**: 0–5, 6–7, 8–37 MMSI, 38–45 reserved, 46–55 SOG, 56 accuracy, 57–84 lon,
85–111 lat, 112–123 COG, 124–132 heading, 133–138 second, 139–140 regional, 141 CS unit, 142 display, 143 DSC,
144 band, 145 msg22, 146 assigned, 147 RAIM, 148–167 radio (20). [conf from memory of gpsd; re-check, unconf]

**Message 19 (extended class B, 312)**: as 18 through 138, 139–142 regional, 143–262 name, 263–270 ship type,
271–279 bow, 280–288 stern, 289–294 port, 295–300 starboard, 301–304 EPFD, 305 RAIM, 306 DTE, 307 assigned,
308–311 spare. [unconf, re-check]

**Message 21 (aid to navigation, 272–360)**: 8–37 MMSI, 38–42 aid type (5; 1 reference point … 31 light vessel/
LANBY), 43–162 name (20 chars), 163 accuracy, 164–191 lon, 192–218 lat, 219–227 bow, 228–236 stern, 237–242 port,
243–248 starboard, 249–252 EPFD, 253–258 second, 259 off-position, 260–267 regional, 268 RAIM, 269 virtual-aid
flag, 270 assigned, 271 spare, 272– name extension (≤ 88 bits). [unconf, re-check]

**Message 24 (class B static, 160/168)**: 8–37 MMSI, 38–39 part number. Part A (0): 40–159 name. Part B (1): 40–47
ship type, 48–65 vendor ID (3 chars), 66–69 unit model, 70–89 serial, 90–131 call sign, 132–140 bow, 141–149 stern,
150–155 port, 156–161 starboard (or 132–161 mothership MMSI if MMSI is 98MIDxxxx), 162–167 spare. [unconf, re-check]

**6-bit text** [conf]: value v (0–63) → ASCII: v < 32 → chr(v+64) ('@','A'..'Z','[','\\',']','^','_'), v ≥ 32 →
chr(v) (space…'?'). '@' (0) is padding; strip trailing '@' and spaces.

**NMEA armoring** (for test vectors and for emitting AIVDM): payload char c → v = ord(c) − 48; if v > 40: v −= 8;
6 bits MSB first. Sentence `!AIVDM,<count>,<num>,<seqid>,<A|B>,<payload>,<fillbits>*<xor checksum hex>`.

**Worked test vector** [conf sentence from gpsd AIVDM doc; decode + on-air bits computed here]:
`!AIVDM,1,1,,B,177KQJ5000G?tO`K>RA1wUbN0TKH,0*5C`
→ type 1, repeat 0, MMSI **477553000**, nav status 5 (moored), ROT 0, SOG 0.0, accuracy 0,
lon **−122.345833**, lat **47.582833**, COG 51.0, heading 181, second 15, manoeuvre 0, RAIM 0, radio 149208.
- 21 payload octets (MSB-first packing): `04 71 db 85 a1 40 00 05 cf f1 fa 1b 3a 24 41 fe 5a 9e 02 46 d8`
- FCS = CRC-16/X.25 = **0xA844** → transmitted octets `44 A8` (each LSB first). Residue over 23 octets = 0xF0B8.
- Bit stuffing inserts 2 zeros.
- On-air pre-NRZI bits (training 24 + flag 8 + stuffed data+FCS 186 + flag 8 = 226 bits):
```
010101010101010101010101 01111110
001000001000111011011011101000011000010100000010000000001010000011110011100011110101111101101100
001011100001001001000001001111101101011010011110010100000001100010000110110010001000010101
01111110
```
  (spaces/line breaks for readability only; the data line is octet-reversed: 0x04 → `00100000`.)
- After NRZI (0 = toggle, start level 1):
```
0011001100110011001100110000000101101010010111100011100001101011101011001010100101010101100101011111011110100000110000001110001010011110101101101101010010000001110011100100000100110101010001011010111000100101101011001100000001
```
Use this as the golden vector for both Python and C++ (`c++/tests/vectors/ais/`).

**Durations**: msg 1/2/3/4/18 burst = 232 bits = **24.2 ms** (26.67 ms slot) → fits a 34 ms block when it falls
inside (a random 34 ms window fully contains a given burst ~29 % of the time; channel load near ports is high, so
detections per block are frequent). Msg 5/19/21/24B are 2–3 slots (≥ 51 ms) → need ≥ 0.1 s; a 0.87 s block holds
~33 slots per channel.

**Recipe**: channelise both 161.975 and 162.025 (they are 50 kHz apart; one 2^21 block at any rate ≥ 1 MSps
centred 162.0 covers both) → 48 kSps (5 SPS) per channel → FM discriminator (arg(x[n]·x*[n−1])) → matched
Gaussian filter → **correlate with the NRZI-encoded training+flag pattern** (32 bits: the 24 alternating training
bits produce a 4800 Hz square tone in the discriminator output; the flag after NRZI gives a unique edge) →
clock phase from the training → slice bits → NRZI decode (bit = 1 if level unchanged) → find 0x7E, de-stuff until
next 0x7E → length check (168 + 16 for single slot; any multiple of 8) → reverse bits in each octet → CRC-16/X.25
residue 0xF0B8 → parse type → fields → label "ais" (verified) with `decoded.fields.{mmsi, msg_type, lat, lon,
sog, cog, heading, nav_status, name, callsign, ship_type, dest}`; MID → flag state (419 India). With CRC failure:
try 1-bit flips (cheap: 200 candidates, use syndrome) — then "ais-like".
Coherent alternatives (Viterbi on GMSK) gain ~2–3 dB; start with the discriminator.

**Synthetic**: payload bits (from an armored sentence) → octets → LSB-first → FCS → stuffing → training + flags →
NRZI → NRZ ±1 → Gaussian filter BT 0.4 (span 3 bits) → FM with h = 0.5 (frequency deviation ±2400 Hz) at 9600 bps
→ ramp (raised-cosine power ramp 8 bits) → add CFO (±(162 MHz × 10 ppm) = ±1.6 kHz), noise, random slot timing.
Oracles: AIS-catcher (GPL-3) and rtl-ais (GPL-2 [unconf]) — test only; libais (Apache-2.0) and pyais (MIT) can
decode NMEA payloads in tests (permissive; still do not vendor into C++ without need).

---

## 6. Marine VHF voice + DSC

**Channel plan (RR Appendix 18)** [conf]:
- Channels 1–28: ship TX = 156.000 + 0.050·n MHz; channels 60–88: ship TX = 156.025 + 0.050·(n − 60) MHz.
  Duplex channels: coast TX = ship TX + 4.600 MHz.
- Key simplex channels: 06 156.300 (inter-ship safety/SAR), 08 156.400, 09 156.450, 10 156.500, 12 156.600 (port),
  13 156.650 (bridge-to-bridge), 14 156.700, 15 / 17 156.750 / 156.850 (on-board, 1 W), **16 156.800 (distress,
  calling)**, 67 156.375, 68 156.425, 69 156.475, **70 156.525 (DSC only)**, 71 156.575, 72 156.625, 73 156.675,
  74 156.725, 77 156.875.
- AIS 87B 161.975, 88B 162.025 (§5). Since WRC-12 some duplex channels are split into 10xx / 20xx simplex variants
  (e.g. 1019 = 156.950 ship-side simplex, 2019 = 161.550 coast-side simplex).
- Voice: F3E/G3E FM, 25 kHz channel, ±5 kHz peak deviation, 6 dB/octave pre-emphasis, 300–3000 Hz audio. In India
  Coast Guard and ports (VTS) use ch 16, 12, 13, 14, 67–71 [unconf].

**DSC (ITU-R M.493-15, VHF class A/B/D)** [conf unless marked]:

| Item | Value |
|---|---|
| Frequency | ch 70, 156.525 MHz |
| Modulation | FM transmitter (G2B), audio subcarrier FSK **1300 Hz = Y = binary 1**, **2100 Hz = B = binary 0** (centre 1700 Hz, ±400 Hz), **1200 baud** |
| Symbol | 10 bits: 7 information bits, **LSB first**, followed by 3 check bits = binary count of **0 (B) bits** in the 7 info bits, MSB first. Symbol values 0–127: 0–99 = two decimal digits, 100–127 = commands |
| Dot pattern | alternating B/Y, 20 bits at VHF (200 bits at MF/HF for some calls) |
| Phasing | DX/RX alternation: DX positions carry symbol **125**; RX positions carry **111, 110, 109, 108, 107, 106** (then 105, 104 while the first message symbols occupy DX) — 12 symbols minimum for phasing [unconf: exact count] |
| Time diversity | every message symbol is sent twice: DX then RX **after 4 other symbols** (RX copy at +5 symbol positions). Receiver combines DX/RX (use whichever passes the 3-bit check) |
| Format specifier (sent twice in DX) | 102 geographic area, **112 distress**, 114 group, 116 all ships, 120 individual, 123 individual automatic |
| Category | 100 routine, 108 safety, 110 urgency, 112 distress |
| Address / self-ID | 5 symbols = 10 digits: the 9-digit MMSI + trailing 0 |
| Nature of distress | 100 fire/explosion, 101 flooding, 102 collision, 103 grounding, 104 listing/capsizing, 105 sinking, 106 disabled and adrift, 107 undesignated, 108 abandoning ship, 109 piracy/armed attack, 110 man overboard, 112 EPIRB emission |
| Position | 5 symbols = 10 digits: quadrant (0 NE, 1 NW, 2 SE, 3 SW), lat DD MM, lon DDD MM; "9999999999" = unknown |
| Time | 2 symbols HH MM (UTC), 8888 = n/a |
| Subsequent comms | 100 = F3E/G3E simplex telephony, 109 = J3E telephony, … |
| EOS | **117** acknowledgement required, **122** acknowledgement (answer), **127** other calls. EOS transmitted 3 times in DX and once in RX [unconf] |
| ECC | one symbol = bitwise XOR of all information symbols from the (first) format specifier through EOS (the format specifier counted once); sent DX and RX |

**Worked example** [computed; synthetic]: distress, MMSI 419000001, nature 106 (disabled and adrift),
position quadrant 0, 12°59'N 077°35'E, 10:30 UTC, simplex telephony:
symbols `112 112 | 41 90 00 00 10 | 106 | 01 25 90 77 35 | 10 30 | 100 | 127` → **ECC = 64**.
10-bit encodings (info LSB first + zero-count MSB first): 0 → `0000000111`, 125 → `1011111001`,
111 → `1111011001`, 112 → `0000111100`, 127 → `1111111000`, 106 → `0101011011`, 64 → `0000001110`.

**Durations**: 1200 baud → 8.33 ms/symbol. Distress call ≈ 20 dot bits + 12 phasing + ~2×18 message + EOS repeats
≈ 55 symbols ≈ **0.47 s**; individual/all-ships calls 0.35–0.5 s. Not decodable in 34 ms; a 0.87 s block holds one
call if aligned → StreamAnalyzer ≥ 1 s. Distress alerts repeat (automatic re-alert every 3.5–4.5 min).

**Recipe**: channel 156.525 ±12.5 kHz → FM demod → audio LPF 3 kHz → tone discriminator (1300 vs 2100 Hz, Goertzel
pair over 1/1200 s or quadrature FM demod of the audio about 1700 Hz) → clock recovery 1200 Hz → bits → search the
phasing pattern (symbols 125 / 111..104 with the 10-bit code, LSB-first) for symbol alignment and DX/RX parity →
de-diversity → check bits → ECC → fields (format, category, MMSIs (self, addressed), nature, position, time) →
label "dsc" (ECC + all symbol checks pass) / "dsc-like". MMSI MID 419 = India.
**Marine voice**: FM voice on an Appendix 18 channel → "marine_vhf_voice-like" with the channel number.

**Synthetic**: symbols → 10-bit words → DX/RX interleave with phasing → dot pattern → continuous-phase FSK audio
1300/2100 Hz at 1200 baud → pre-emphasis → FM ±3 kHz. Oracle: multimon-ng has no DSC; YaDDNet / "DSC Decoder"
are closed; use the M.493 annex examples (not retrieved) [unconf] — self-consistency only.

---

## 7. APRS (AX.25 v2.2 + APRS 1.0.1)

| Item | Value | Status |
|---|---|---|
| India frequency | Commonly cited **144.390 MHz** for terrestrial APRS by Indian amateurs, but not found in an Indian primary source; IARU R3 recommends 144.390 or 144.640 for some countries; satellite APRS 145.825 (ISS, AMSAT-India AISAT-1 digipeater) | [unconf] — make the detector frequency-agnostic over 144–146 MHz |
| RF | NBFM, ±3 kHz (12.5 kHz) or ±5 kHz (25 kHz) | [conf] |
| Modem | **Bell 202 AFSK**, 1200 baud, **mark 1200 Hz, space 2200 Hz** | [conf] |
| Line coding | **NRZI: 0 = tone change, 1 = no change** | [conf] |
| Framing | HDLC flags 0x7E (TXDELAY = many flags, typ. 100–500 ms), bit stuffing after five 1s, octets LSB first | [conf] |
| FCS | **CRC-16/X.25** (same as AIS/AVLC), low byte first; residue 0xF0B8 | [conf] [computed] |
| Frame | dest (7) + src (7) + 0–8 digipeaters (7 each) + control **0x03** (UI) + PID **0xF0** (no L3) + info (≤ 256) + FCS | [conf] |
| Address octet coding | 6 callsign chars ASCII uppercase, space-padded, **each shifted left 1 bit**; 7th octet = `C R R S S S S E` → 0b C 1 1 SSID(4) E: bits 6–5 reserved = 1 (0x60), bits 4–1 SSID, bit 0 = extension (1 on the last address). C bit: dest C=1/src C=0 for command (v2). Digipeater 7th octet bit 7 = H ("has been repeated") | [conf] |
| Dest field | APRS "tocall" (APxxxx identifies the software, e.g. APDW16 = Dire Wolf), or Mic-E encoded latitude | [conf] |

**Worked example** [computed; synthetic]: `VU2ABC-9>APRS:!1259.00N/07735.00E>Test`
AX.25 octets: `82 a0 a4 a6 40 40 e0 | ac aa 64 82 84 86 73 | 03 f0 | 21 31 32 35 39 2e 30 30 4e 2f 30 37 37 33 35 2e 30 30 45 3e 54 65 73 74`
→ FCS = **0xB7A6** → transmitted `A6 B7`; residue 0xF0B8. ('A' = 0x41 → 0x82; SSID 9 last-address octet = 0x60|9<<1|1 = 0x73; dest SSID octet with C bit = 0xE0.)

**APRS position formats** (APRS 1.0.1 ch. 8–10) [conf]:
- Data type identifier (first info char): `!` / `=` position without timestamp (`=` = messaging capable), `/` / `@`
  with timestamp (7 chars: DDHHMMz, DDHHMM/, HHMMSSh), `` ` `` / `'` Mic-E, `:` message, `>` status, `;` object,
  `)` item, `T` telemetry, `_` weather without position, `$` raw NMEA.
- **Uncompressed**: `DDMM.hhN` + symbol table (`/` primary, `\` alternate, or overlay) + `DDDMM.hhE` + symbol
  code; optional `CSE/SPD` (course deg / speed knots), `PHGphgd`, `/A=aaaaaa` altitude feet in comment. Position
  ambiguity: digits replaced by spaces.
- **Compressed** (13 chars): symbol table, `YYYY` (base-91 lat), `XXXX` (base-91 lon), symbol code, `cs`, `T`.
  Each char value = ord − 33. lat = 90 − (y1·91³+y2·91²+y3·91+y4)/380926; lon = −180 + (x1·91³+…)/190463.
  cs: if c ≠ ' ': course = (c−33)·4°, speed = 1.08^(s−33) − 1 knots (or altitude / range by T bits).
- **Mic-E**: destination callsign (6 chars) encodes 6 latitude digits + 3 message bits + N/S + longitude offset
  (+100°) + W/E; info: `` ` ``/`'`, then lon degrees d+28, minutes m+28, hundredths h+28, speed/course SP+28,
  DC+28, SE+28, symbol code, symbol table; then optional altitude / telemetry / device suffix. (Table 10 of APRS
  1.0.1; implement from the spec.)

**Durations**: 1200 bps; a 40-octet position frame = (40+2+flags)·8·~1.03 ≈ 0.31 s + TXDELAY 0.1–0.5 s. Not in 34
ms; a 0.87 s block catches about one in two packets whole → StreamAnalyzer ≥ 2 s.

**Recipe**: channel 25 kHz → FM demod → audio → correlators at 1200 / 2200 Hz (or delay-and-multiply) → clock
1200 → NRZI decode → flag search, de-stuff → LSB-first octets → FCS residue → address parse → APRS info parse
(position, symbol, comment) → label "aprs" (FCS) with `decoded.fields.{src, dest, path, lat, lon, symbol,
comment}`; any AX.25 UI frame with valid FCS but non-APRS info → "ax25". The same HDLC/AFSK chain also decodes
1200-baud packet in other bands (e.g. 145.825 satellite digipeaters).

**Synthetic**: frame above → stuffing → NRZI → continuous-phase AFSK 1200/2200 at 1200 baud (note: the tones are
not phase-coherent with bit boundaries; CPFSK) with 50 flags of TXDELAY → optional pre-emphasis → FM ±3 kHz.
Oracles: direwolf (GPL-2), multimon-ng (GPL-2) — test only.

---

## 8. FM broadcast + RDS (IEC 62106 / EN 50067)

**India FM** [conf for band; raster unconf]:
- 88–108 MHz; carriers on a **100 kHz raster** (e.g. AIR FM Gold 106.4 / 100.1, FM Rainbow 102.6, private 91.1,
  92.7, 93.5, 94.3, 98.3, 104.0 in Delhi/Mumbai [unconf individual stations]); in a city stations are ≥ 800 kHz
  apart (Phase III policy) [unconf].
- ±75 kHz peak deviation, 50 µs pre-emphasis, ~200 kHz occupied (Carson 2·(75+57) ≈ 264 kHz incl. RDS).

**MPX baseband** [conf]: L+R 30 Hz–15 kHz; **pilot 19 kHz ± 2 Hz** (8–10 % deviation, ~6.75–7.5 kHz); L−R DSB-SC
on 38 kHz (23–53 kHz); **RDS at 57 kHz ± 6 Hz** (= 3 × pilot, phase-locked to the pilot's 3rd harmonic,
in-phase or quadrature), injection ~2–4 kHz deviation (typ. 2 kHz); ±2.4 kHz occupied.

**RDS physical layer** [conf]:
- Bit rate **1187.5 bit/s** = 57 000 / 48.
- **Differential encoding** at the transmitter: e_k = d_k ⊕ e_{k−1}; each encoded bit becomes a **biphase**
  (Manchester) symbol: +pulse then −pulse (or inverse) over one bit period, shaped so the overall spectrum is a
  cosine roll-off (zero at 2.4 kHz from the subcarrier). The carrier is suppressed (DSB-SC BPSK).
- Demod: recover 57 kHz (from 3× pilot, or Costas loop), BPSK matched filter at the 2375 Hz biphase symbol rate,
  bit clock 1187.5 Hz, then differential decode d_k = e_k ⊕ e_{k−1} (removes 180° ambiguity).

**Block / group** [conf; offsets and syndromes cross-checked in redsea `block_sync.cc` (MIT)]:
- Group = 4 blocks × 26 bits = **104 bits = 87.6 ms**; block = 16 information bits + 10-bit checkword, MSB first.
- Generator **g(x) = x^10 + x^8 + x^7 + x^5 + x^4 + x^3 + 1** (0x5B9), shortened cyclic (26,16) code.
- Checkword = (d(x)·x^10 mod g(x)) ⊕ offset word.
- Offset words (10 bits): **A = 0x0FC** (0011111100), **B = 0x198** (0110011000), **C = 0x168** (0101101000),
  **C' = 0x350** (1101010000), **D = 0x1B4** (0110110100); E = 0x000 (RBDS MMBS).
- Syndrome decoding: (received 26-bit word) mod g(x) = offset word for an error-free block (simplest form; this is
  what the [computed] example below shows). EN 50067 instead multiplies by its parity-check matrix H, giving the
  "syndromes" A 1111011000, B 1111010100, C 1001011100, C' 1111001100, D 1001011000 (redsea's table). Both are
  fine; pick one and use it in both Python and C++.
- Sync: find a bit position where consecutive 26-bit windows give offsets in sequence A→B→C/C'→D (two consecutive
  matches with correct spacing ⇒ lock). Burst-error correction up to 5 bits per block (syndrome → error-pattern
  table), applied only once locked (redsea does this).

**Group content** [conf]:
- Block 1 always **PI** (16 bits): country code (4) | area coverage (4) | programme reference (8).
  **India: PI country nibble 0x5 with ECC 0xF2** (redsea `tables.cc` country table, from EN 50067 Annex D /
  RDS Forum R08/008_7: ECC F2 list position 5 = "in"). ECC is carried in group 1A, variant 0.
- Block 2: group type (4) | B0 version (1) | TP (1) | PTY (5) | 5 group-specific bits.
- **0A/0B** (basic tuning, PS): block 2 low bits = TA (1), M/S (1), DI segment bit (1), PS segment address c1c0
  (2). 0A: block 3 = two AF codes (AF code n → 87.5 + 0.1·n MHz, n = 1..204; 224–249 = "N AFs follow");
  0B: block 3 = PI (offset C'). Block 4 = 2 PS characters at positions 2·c..2·c+1. **PS (8 chars) needs 4 groups
  ≥ 350 ms** (and stations interleave other groups: typically ≥ 1 s).
- **2A/2B** (RadioText): block 2 low bits = A/B text flag (1), segment address (4). 2A: blocks 3+4 = 4 chars
  (64 chars total, 16 groups); 2B: block 4 = 2 chars (32 chars). 0x0D terminates.
- **4A** (clock-time): block 2 low 2 bits = MJD bits 16–15; block 3 = MJD bits 14–0 (15) + hour bit 4; block 4 =
  hour bits 3–0, minute (6), offset sign (1), offset in half-hours (5). **India +5:30 → sign 0, offset 11.** Sent
  once per minute at the minute edge. MJD → date (EN 50067 Annex G): Y' = ⌊(MJD − 15078.2)/365.25⌋,
  M' = ⌊(MJD − 14956.1 − ⌊Y'·365.25⌋)/30.6001⌋, D = MJD − 14956 − ⌊Y'·365.25⌋ − ⌊M'·30.6001⌋, K = 1 if M' ∈ {14,
  15} else 0, Y = 1900 + Y' + K, M = M' − 1 − 12K.
- PTY (RDS, Europe/India table): 0 none, 1 News, 2 Current affairs, 3 Information, 4 Sport, 5 Education,
  6 Drama, 7 Culture, 8 Science, 9 Varied, 10 Pop, 11 Rock, 12 Easy listening, 13 Light classical, 14 Serious
  classical, 15 Other music, 16 Weather, 17 Finance, 18 Children's, 19 Social affairs, 20 Religion, 21 Phone-in,
  22 Travel, 23 Leisure, 24 Jazz, 25 Country, 26 National music, 27 Oldies, 28 Folk, 29 Documentary, 30 Alarm test,
  31 Alarm.
- RDS use in India: some AIR and private stations transmit PS/RT, many transmit none [unconf] → "fm_broadcast"
  must not require RDS.

**Worked example** [computed]: group 0A, PI **0x5201** (India-style), TP 0, PTY 10, TA 0, M/S 0, DI 1, segment 0,
AF "1 AF follows" + 87.6 MHz, PS chars "AI":

| Block | Offset | Data | Checkword | 26-bit word |
|---|---|---|---|---|
| 1 | A | 0x5201 | 0100010001 | 01010010000000010100010001 |
| 2 | B | 0x0144 | 1000101110 | 00000001010001001000101110 |
| 3 | C | 0xE001 | 0111101110 | 11100000000000010111101110 |
| 4 | D | 0x4149 | 0100000010 | 01000001010010010100000010 |

Each 26-bit word mod g(x) returns its offset word (0x0FC, 0x198, 0x168, 0x1B4) — the receiver's validity test.

**Durations**: FM carrier + stereo pilot: one 34 ms block is plenty (pilot = 19 kHz tone in the demodulated MPX:
SNR-gain over 34 ms is large). RDS: one 34 ms block holds ~40 bits = one full block plus fragments → at most a PI
candidate from a single A-offset hit (weak; 10-bit check ≈ 1/1024 false rate per position × 40 positions ⇒ ~4 %
false; require the PI to repeat across blocks/dwells). **PI reliably: ≥ 2 blocks of one group (~45 ms) — 2^21 at
≤ 30.72 MSps or a dwell; PS: ≥ 0.35 s (best ≥ 1 s); RT: ≥ 1.4–5 s; CT: up to 60 s.** A 0.87 s block (2.4 MSps)
gives PI, PTY, TP, usually PS.

**Recipe**: detect 200 kHz-wide FM carriers on the 100 kHz raster (block FFT) → per carrier: channelise to
~228 kSps (= 4 × 57 kHz) → FM demod → MPX → pilot detection (19 kHz line, Goertzel; "stereo") → RDS: bandpass
57 kHz ±3 kHz, mix down with 3 × pilot phase (or Costas), LPF, biphase matched filter, clock 1187.5 Hz, BPSK
decisions, differential decode → block sync with offsets → error correction → groups → PI/PTY/PS/RT/CT/ECC →
label "fm_broadcast" (always) + "rds" (when ≥ 2 blocks pass in sequence). Country from PI nibble + ECC.

**Synthetic**: encode groups (above) → differential → biphase symbols at 2375 Hz → shape (raised cosine), × cos(2π·
57k·t) phase-locked to the pilot sin(2π·19k·t) → MPX = 0.9·audio(L+R) + 0.09·pilot + 0.03–0.05·RDS (+ L−R DSB-SC)
→ FM ±75 kHz at e.g. 1.92 MSps. Oracle: redsea (**MIT**, verified LICENSE file — permissive, may even be consulted
for C++ structure), gr-rds (GPL-3, test only).

---

## 9. VOR and ILS localizer (ICAO Annex 10 Vol I §3.3 and §3.1.3)

### VOR (conventional CVOR and Doppler DVOR)

| Item | Value | Status |
|---|---|---|
| Frequencies | 108.000–117.950 MHz, 50 kHz; in 108–111.975 only **even tenths and even tenths + 50 kHz** (108.00, 108.05, 108.20, 108.25 …) — odd tenths are ILS LOC | [conf] |
| Carrier | horizontally polarised, A9W composite AM | [conf] |
| 30 Hz AM | depth 28–32 % (space modulation at the receiver) | [conf] |
| 9960 Hz subcarrier | AM depth 28–32 %, **FM by 30 Hz with deviation ±480 Hz (16 ± 1 deviation ratio)** | [conf] |
| CVOR | reference phase = 30 Hz FM on the 9960 Hz subcarrier; variable phase = 30 Hz AM (rotating cardioid). Bearing (magnetic radial FROM the station) = phase lag of variable relative to reference, 0–360° unambiguous | [conf] |
| DVOR | reference = 30 Hz AM (omni carrier); variable = 30 Hz FM produced by the Doppler shift of a sideband (±9960 Hz) commutated around a ~13.5 m (≈ 5 λ) circle at 30 rps. The rotation direction is reversed relative to CVOR (**sense reversal**) so that an unmodified receiver computes the same bearing | [conf] |
| Ident | Morse, **1020 Hz ± 50 Hz**, ~7 words/min, 2–3 letters, AM depth ≤ 10 %, at least once every 30 s (VOR/DME: 3 VOR idents per 1 DME ident) | [conf] |
| Voice | optional 300–3000 Hz AM ≤ 30 % (ATIS on some VORs) | [conf] |
| Accuracy | ±2° (radial, at the station) | [conf] |

**Bearing algorithm**: channel ±12.5 kHz → envelope a(t) = |x(t)| (normalise by its mean) → (i) variable/reference
30 Hz AM: LPF 0–60 Hz of a(t), complex correlate with e^{−j2π·30t} → phasor P_am; (ii) 9960 Hz subcarrier: BPF
9960 ± 700 Hz of a(t) → analytic signal → instantaneous frequency f(t) − 9960 → 30 Hz correlation → P_fm. Bearing
θ = arg(P_am) − arg(P_fm) (CVOR; for DVOR identical by design) mod 360. Calibrate fixed group-delay differences
of the two filter paths with a synthetic VOR at known radials (0, 90, 180, 270).
Validity: subcarrier present (9960 Hz peak), 30 Hz AM and FM at the expected depths (AM 0.28–0.32 nominal,
deviation 400–560 Hz) → "vor" label; ident decode (Morse at 1020 Hz) needs ~5–10 s of dwell → `fields.ident`.
**How much signal**: one 30 Hz period (33.3 ms) gives one bearing estimate — exactly what a 34 ms 61.44 MSps block
holds, but the window must contain an integer number of 30 Hz periods to avoid leakage between the 30 Hz and DC
terms (use exactly 1/30 s = 2 048 000 samples at 61.44 MSps, or a Hann-windowed correlation over 2+ periods).
At moderate SNR one period gives ±2–5° [unconf: simulate]; ≥ 0.2–0.5 s (6–15 periods) for ≤ 1°.
Also note the receiver's own position cannot be resolved from one VOR; bearing is to the station, which is
identified by frequency (+ ident) from a VOR database (AIP India ENR 4.1).

### ILS localizer

| Item | Value | Status |
|---|---|---|
| Frequencies | 108.10–111.95 MHz, **odd tenths and odd tenths + 50 kHz**: 108.10, 108.15, 108.30, 108.35, … 111.90, 111.95 → 40 channels, each paired with a glide-path frequency in 329.15–335.00 MHz (file 05) and a DME channel | [conf] |
| Modulation | carrier AM by **90 Hz and 150 Hz**, each **20 % depth on the course line** (CSB); SBO (sideband-only) antenna pattern makes one tone dominate off course | [conf] |
| DDM | DDM = (m_90 − m_150) (modulation depths as fractions); **0 on course**; **0.155 at the edge of the course sector** (full-scale deflection, sector half-width typically ±2–3°); 90 Hz predominates to the **left** of course as seen from an approaching aircraft (right side → 150 Hz) | [conf] |
| SDM | m_90 + m_150 = 40 % (36–44 %) | [conf] |
| Two-frequency (capture effect) | course and clearance carriers separated by **5–14 kHz** around the channel frequency | [conf] |
| Tones | 90 and 150 Hz phase-locked (both harmonics of 30 Hz) | [conf] |
| Ident | Morse 1020 Hz ± 50 Hz, depth 5–15 %, ~7 wpm, typically 4 letters starting with "I" (e.g. "IBBB") | [conf] ident "I…" convention [unconf for India] |

**DDM from one block**: the LCM of the 90 Hz and 150 Hz periods is 1/30 s → **one 33.3 ms window contains exactly 3
and 5 cycles**, so a 34 ms block measures m_90, m_150 cleanly (correlate the normalised envelope with e^{−j2π·90t}
and e^{−j2π·150t} over exactly 1/30 s). Validity: both tones present, SDM 0.36–0.44, channel on the LOC raster →
"ils_localizer" with `fields.{ddm, sdm, side}`; the DDM tells only which side of the extended runway centre-line the
receiver is (and roughly how far, in angle), with SDM sanity. Two-frequency systems: process both carriers or the
composite envelope.

**Marker beacons (75 MHz, if any remain)**: A2A AM 95 %, OM 400 Hz dashes 2/s, MM 1300 Hz alternate dot-dash,
IM 3000 Hz dots 6/s [conf]. Detection only.

**Synthetic**: VOR: x(t) = [1 + 0.3 cos(2π30t − θ) + 0.3 cos(2π9960t + (480/30) sin(2π30t)) + 0.1·ident(t)]·e^{jΔωt}
(CVOR; θ = radial). LOC: x(t) = 1 + m_90 sin(2π90t) + m_150 sin(2π150t) with m_90 = 0.2 + DDM/2, m_150 = 0.2 −
DDM/2 (+ ident). Add CFO (±1.2 kHz at 110 MHz for ±10 ppm), noise, Pluto DC spur.

---

## Summary table: what fits where

| Signal | Minimum dwell for detection | Minimum dwell for full decode | 34 ms block? | 0.87 s block? | Verifier |
|---|---|---|---|---|---|
| Airband AM voice | 1 block (carrier + envelope) | n/a (no data) | detect | yes | modulation type only ("-like") |
| ACARS | ~80 ms (one short block) | 0.1–0.8 s | **no** | mostly | BCS CRC-16/KERMIT + odd parity |
| VDL Mode 2 | ~3 ms (training) | 10–70 ms per burst | short bursts yes | yes | header FEC + RS + FCS X.25 |
| AIS | 24.2 ms burst | same (msg 5: 51 ms) | msg 1–4/18 yes | yes | FCS X.25 |
| DSC | ~0.5 s | ~0.5 s | **no** | if aligned | 3-bit symbol checks + ECC + DX/RX |
| APRS | 0.3–0.8 s | same | **no** | ~half the packets | FCS X.25 |
| FM broadcast | 1 block | stereo pilot: 1 block | yes | yes | pilot 19 kHz |
| RDS | 45 ms (PI, 2 blocks) | PS 0.35–1 s; RT 1.4–5 s; CT ≤ 60 s | PI candidate only | PI/PTY/PS | (26,16) offset checkwords |
| VOR | 1 block | bearing 33 ms (coarse) / 0.2–0.5 s; ident 5–10 s | coarse bearing | yes | 9960 Hz subcarrier + depths |
| ILS LOC | 1 block | DDM 33.3 ms; ident 5–10 s | **yes (DDM)** | yes | 90/150 Hz + SDM |

---

## Sources (with licences)

Standards / primary documents (no code; implement from these):
1. ITU-R M.1371-5, "Technical characteristics for an automatic identification system using TDMA in the VHF maritime
   mobile band" — https://www.itu.int/rec/R-REC-M.1371 (ITU publication; free download).
2. ITU-R M.493-15, "Digital selective-calling system for use in the maritime mobile service" —
   https://www.itu.int/rec/R-REC-M.493 (ITU; free download).
3. ITU Radio Regulations, Article 5 (Region 3) and Appendix 18 (maritime VHF channels) —
   https://www.itu.int/pub/R-REG-RR (ITU).
4. ICAO Annex 10 Vol I (VOR §3.3, ILS §3.1, markers §3.1.7), Vol III (VDL, ACARS-related), Vol V (frequency
   utilisation, 8.33 kHz naming) — ICAO store (copyright ICAO, purchase).
5. ICAO Doc 9776, Manual on VHF Digital Link (VDL) Mode 2 — ICAO (purchase).
6. ETSI EN 301 841-1, VDL Mode 2 ground equipment, physical layer — https://www.etsi.org (free download).
7. ARINC 618 (air/ground character-oriented protocol) and ARINC 620 (labels / ground-ground) — SAE-ITC/ARINC
   (purchase). Not read here; facts below came via sources 13 and 18.
8. AX.25 Link Access Protocol v2.2 (TAPR, 1998) — https://www.ax25.net/AX25.2.2-Jul%2098-2.pdf (TAPR, free).
9. APRS Protocol Reference 1.0.1 (APRS Working Group, 2000) — http://www.aprs.org/doc/APRS101.PDF (free).
10. IEC 62106 (RDS; current edition 2015+ parts) / CENELEC EN 50067:1998 — IEC webstore (purchase);
    https://webstore.iec.ch/publication/6467 (2009 ed., withdrawn).
11. NFAP-2022, National Frequency Allocation Plan, WPC / DoT India — https://dot.gov.in (government publication);
    WG-1 minutes (174–230 MHz audio proposal) — https://dot.gov.in/sites/default/files/WG1_Minutes%20of%20the%204th%20meeting%20dated%2012th%20Dec%2024.pdf
12. TRAI consultation paper 10-12-2022 quoting DoT on ACARS assignments (SITA 131.725 / 136.975, BCS 131.825) —
    https://trai.gov.in/sites/default/files/2024-09/Consultation_10122022.pdf (government publication).

Code read for protocol facts only (test oracles; never copied into python/ or c++/):
13. acarsdec, Thierry Leconte — https://github.com/TLeconte/acarsdec — source headers: **GNU Library GPL v2
    (LGPL-2.0)**; README says GPL. Treat as copyleft, oracle only. Facts used: SYN/SOH/ETX(0x83)/ETB(0x97)/DEL
    constants, CRC table (reflected 0x8408, init 0), CRC range mode..ETX, LSB-first bits, MSK demod at 1800 Hz.
14. dumpvdl2, Tomasz Lemiech — https://github.com/szpajder/dumpvdl2 — **GPL-3.0** (LICENSE). Facts used: training
    phases, Gray table, LFSR x^15+x+1 seed 0x6959, header H matrix, length field LSB-first, RS (255,249) 0x187 /
    fcr 120, partial-block parity 0/2/4/6, interleave, AVLC addresses and 0xF0B8 residue, ACARS `FF FF 01` prefix.
    Ships test WAV `test/vdl2_model_16b_1050kHz.wav`.
15. redsea, Oona Räisänen — https://github.com/windytan/redsea — **MIT** (LICENSE verified). Facts used: syndromes,
    offset sequence, country table (ECC 0xF2, cc 5 = India).
16. gr-rds, Bastian Bloessl — https://github.com/bastibl/gr-rds — **GPL-3.0** (COPYING). Oracle only.
17. AIS-catcher, jvde-github — https://github.com/jvde-github/AIS-catcher — **GPL-3.0** (README). Oracle only.
18. airframes.io ACARS docs — https://docs.airframes.io/docs/technology/acars/intro (docs; licence not stated).
19. gpsd "AIVDM/AIVDO protocol decoding", Eric S. Raymond — https://gpsd.gitlab.io/gpsd/AIVDM.html (document
    licence not stated on the page; gpsd code is BSD-2). Test sentence and field tables.
20. rtl-ais — https://github.com/dgiardini/rtl-ais — GPL-2 [unconf]. Oracle only.
21. direwolf — https://github.com/wb2osz/direwolf — GPL-2. Oracle for AX.25/APRS.
22. multimon-ng — https://github.com/EliasOenal/multimon-ng — GPL-2. Oracle for AFSK1200.
23. libais (Apache-2.0) https://github.com/schwehr/libais and pyais (MIT) https://github.com/M0r13n/pyais —
    permissive AIS payload decoders for Python tests [licences from memory, unconf].
24. Signal Identification Wiki, APRS / ACARS pages — https://www.sigidwiki.com/wiki/APRS ,
    https://www.sigidwiki.com/wiki/ACARS (CC BY-SA; secondary).
25. Wikipedia, "List of amateur radio frequency bands in India" —
    https://en.wikipedia.org/wiki/List_of_amateur_radio_frequency_bands_in_India (CC BY-SA; secondary).
26. ATC frequency aggregators for VABB / VIDP (learn-atc.com, acukwik.com, metarcentral.com) — secondary,
    conflicting for VIDP; authoritative source is AIP India (https://aim-india.aai.aero).

## Open items / uncertainties
- India APRS frequency (144.390 vs other); India airband 8.33 kHz status; current ACARS use of 131.725 / 131.825.
- NFAP-2022 IND footnotes for 70–88, 138–144, 146–174, 230–300 MHz (government mobile) — rows are placeholders.
- Band III (174–230): residual analog TV / any digital audio in India.
- ACARS: pre-key length and the exact tone↔bit (differential) mapping — settle by loop-back against acarsdec.
- VDL2: ramp-up symbol count; scrambler seed in EN 301 841 notation (dumpvdl2 register form is given).
- AIS msg 18/19/21/24 bit tables were written from memory of gpsd — re-check against M.1371-5 Annex 8 before coding.
- DSC: phasing length (12 vs. more symbols), EOS repetition pattern — check M.493-15 Annex 1 §3–4.
- RDS: which Indian stations carry RDS; PI area/programme allocation in India.
- VOR bearing error vs. SNR for a single 33.3 ms window — simulate before promising per-block bearings.
