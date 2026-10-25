# 01 — Satellite phones and terminals in L-band (SIGINT gap research)

**Batch S, items S1 / S2.** Goal: detect a satellite phone *in use* over a Pluto block (300 MHz–6 GHz, India).
Sat phones need a DoT licence in India, so **handset uplinks are the high-value evidence** (see §6, legal). Downlinks
(sat → handset) confirm that a terminal is being paged / served nearby. Typical input: 2^21 cx samples, 34 ms at
61.44 MSps; min 65 536 samples. All L-band systems here fit comfortably in one block's *bandwidth* (each is ≤ 1.23 MHz),
but **34 ms is only ~0.4 of an Iridium 90 ms frame and far less than a Thuraya 640 ms multiframe or an Inmarsat-C
8.64 s frame** — so single-block capability is mostly *detect + classify + partial decode of one burst*, not full
system-information decode. This is a research reference; everything here is validated on **synthetic** signals (and the
GPL/AGPL oracles noted) and must be labelled `*-like*` until a real capture verifies it.

Band summary (all MSS L-band, India relevant):

| System | Downlink (sat→MES) | Uplink (MES→sat, **handset**) | Access | Notes |
|---|---|---|---|---|
| **Iridium** | 1616–1626.5 MHz | 1616–1626.5 MHz (**same band**, TDMA/FDMA) | TDMA+FDMA, DQPSK 25 ksym/s | duplex 1616–~1626.0; simplex 1626.0–1626.5 |
| **Thuraya / GMR-1** | 1525–1559 MHz | 1626.5–1660.5 MHz | TDMA, π/4-CQPSK 23.4 ksym/s | 31.25 kHz raster; UL–DL offset 101.5 MHz |
| **Inmarsat (GMR-2+/C/Aero/BGAN)** | 1525–1559 MHz | 1626.5–1660.5 MHz | varies | IsatPhone, Inmarsat-C 1200 sym/s, BGAN |
| **Globalstar CDMA** | 2483.5–2500 MHz (S-band) | 1610–1618.725 MHz | CDMA (IS-95-like), 1.2288 Mcps | UL overlaps Iridium 1616–1618.725 |
| **Globalstar simplex / SPOT** | 6875–7055 MHz (C-band) | 1611.25 / 1613.75 / 1616.25 / 1618.78 MHz | DSSS BPSK | trackers/PLBs, not phones |

> **Pluto caveat (from CLAUDE.md / FINDINGS):** no preselector. Strong 1.8–2.6 GHz downlinks fold into 330–830 MHz via
> 3rd/5th LO harmonics; the ADC clips in 2122–2185 and 2280–2405 MHz. 1610–1626.5 MHz is **not** in the listed clip
> bands, so L-band sat is capturable directly, but confirm any 1616–1618.7 MHz detection isn't a fold of a strong
> downlink. +7.4–10 ppm Pluto TCXO ≈ +12–16 kHz at 1.62 GHz → CFO search must cover ±20 kHz plus Iridium Doppler
> (±36 kHz LEO).

---

## 1. IRIDIUM

LEO constellation (66+ sats, 780 km). The one L-band service where **uplink and downlink share the band**, so a bare
"burst at 1620 MHz, 25 kBd, DQPSK" cannot by itself be called handset vs satellite — the **unique word** decides.

### 1.1 PHY constants (confirmed: sigidwiki, gr-iridium `iridium.h`, patents US5596315/US5603079)

| Parameter | Value |
|---|---|
| Band | 1616.0–1626.5 MHz |
| Channel raster | **41.667 kHz** (= 10 MHz / 240). Centre of channel n (1..240): `1616 + 0.0208333·(2n−1)` MHz |
| Sub-band / FA grid (toolkit) | `base=1616 MHz`, `channel_width = 1e7/(30·8) = 41 666.67 Hz`; 30 sub-bands × 8 freq-accesses (duplex) + simplex sub-band "S" with 12 FAs |
| Simplex channels | 1626.0–1626.5 MHz, 12 accesses. **Ring alert (ch 7) = 1626.270833 MHz**; messaging ch 11 = 1626.4375, ch 10 = 1626.395833, ch 4 = 1626.145833, ch 3 = 1626.104167 MHz |
| Symbol rate | **25 000 sym/s** (DQPSK → 50 kbps) |
| Modulation | differentially-encoded QPSK (π/2? no — plain DQPSK), RRC roll-off α = 0.4, 51-tap (gr-iridium) |
| Occupied BW | ~31.5 kHz |
| Frame period | **90 ms** L-band frame = 2250 symbols @ 25 kBd. Structure: 1 ringing/paging slot + 4 uplink slots + 4 downlink slots; simplex channels live in the guard between DL and UL blocks and occupy **20.48 ms** of the 90 ms (US5603079) |
| Burst | ~**8.28 ms** (~9.2 % duty) |
| Duplex vs simplex split (toolkit) | `f_simplex` lower bound ≈ 1626.104 MHz, `f_duplex` upper bound ≈ 1625.979 MHz; packets above → simplex frame types, below → duplex |

**34-ms-block reality:** one block holds ≈ 0.38 of a 90 ms frame → expect **0–several complete bursts** (each 8.28 ms).
Each burst is self-contained (preamble+UW+payload), so per-burst decode is feasible; multi-frame reassembly (voice,
SBD) is not.

### 1.2 Burst structure, preamble and Unique Words (confirmed: gr-iridium, iridium-toolkit `bitsparser.py`)

Symbol constellation used for sync generation: `s0 = (+1,+1)`, `s1 = (−1,−1)` (i.e. BPSK on the I=Q diagonal).
UW is expressed as **quaternary symbols 0/1/2/3** = DQPSK phase states.

- **Preamble:** short = **16 symbols**; long = **64 symbols** (simplex & broadcast frames).
  - Downlink preamble: `s0` repeated (constant phase).
  - Uplink preamble: `s1, s0` alternating (a tone at ±symbol-rate/2).
- **Unique words (12 symbols each):**
  - `UW_DOWNLINK = [0,2,2,2,2,0,0,0,2,0,0,2]`
  - `UW_UPLINK   = [2,2,0,0,0,2,0,0,2,0,2,2]`
  - Iridium-NEXT variants: `NXT_UW_DOWNLINK=[2,2,0,2,2,0,2,0,2,0,2,2]`, `NXT_UW_UPLINK=[0,2,0,0,0,0,0,0,2,0,2,0]`
- **Access codes (24 bits = 12 symbols, BPSK bit form used by the toolkit after demod):**
  - DL `iridium_access = 001100000011000011110011` (0x789h in BPSK)
  - UL `uplink_access  = 110011000011110011111100` (0xc4b)
  - NEXT DL `110011110011111111111100` (0xdab); NEXT UL `001111000000000011111111` (0x40a)
- **Direction decision (the handset discriminator):** demod to 12 quaternary symbols, compare to UW_DL vs UW_UL with
  symbol-distance (a 90° diff counts as 1); accept if min distance ≤ 2. UL match ⇒ **handset/MES uplink burst**.

Frame length gates (gr-iridium): duplex frames 131–191 symbols (IBC is exactly 131); simplex 80–444 symbols
(single-page IRA ≥ 80).

### 1.3 Demod / bit recovery (confirmed vs gr-iridium `iridium_qpsk_demod_impl.cc`)

1. Burst detect (energy/FFT), estimate CFO via FFT over preamble (square the signal: BPSK preamble → tone).
2. RRC matched filter (α 0.4), timing recovery, 1 sample/symbol.
3. First-order PLL to remove residual CFO; QPSK quadrant hard-decision:
   `(I≥0,Q≥0)→0, (I≥0,Q<0)→3, (I<0,Q<0)→2, (I<0,Q≥0)→1`.
4. **Differential decode:** `bits = (sym − prev_sym) mod 4`, then map `{0→0,1→2,2→3,3→1}`; each symbol → 2 bits
   (`b1 = sym&2`, `b0 = sym&1`).
5. Strip access code; branch on frame type.

### 1.4 FEC — BCH and interleaving (confirmed: iridium-toolkit `bch.py`, decoder reproduces oracle)

Generator polynomials (integer → binary, MSB = highest degree):

| Use | Poly (dec) | Binary | Degree |
|---|---|---|---|
| IBC header (6→?) | 29 | `11101` | 4 |
| LCW block 1 (BCH(7,3)) | 29 | `11101` | 4 |
| LCW block 3 (BCH(26,21)) | 41 | `101001` | 5 |
| LCW block 2 (BCH ~14/13) | 465 | `111010001` | 8 |
| Ring-alert / broadcast (RA, BC) | **1207** | `10010110111` | 10 |
| Messaging (MS) | **1897** | `11101101001` | 10 |
| ACCH | 3545 | `110111011001` | 11 |

**Standard protected block = 32 bits:** `[21 data][10 BCH][1 even-parity]`. Decode: take bits 0..30, BCH-correct
(≤2 errors via syndrome table), then check overall even parity across `data+bch+parity`. Uncorrectable → cut packet.

**Interleavers** (toolkit, bit-exact — I reproduced them in the encoder):
- 2-way `de_interleave(group)`: `symbols[k]=group[2k+1]+group[2k]`; `even = symbols[-2::-2]`, `odd = symbols[-1::-2]`;
  returns `(odd, even)`. Used for MS/BC 64-bit groups.
- 3-way `de_interleave3` (RA first 96 bits): same symbol swap, strided by 3 → `(first, second, third)`.
- LCW `de_interleave_lcw`: explicit 46-entry permutation table (see `bitsparser.py` `de_interleave_lcw`), splits into
  7+13+26.

Reed–Solomon (RS(15,9)/GF(2⁴)) is used only for BACH/voice-vocoder frames (`rs.py`, `rs6.py`) — out of scope for
detection/paging decode.

### 1.5 Frame types (freq-gated, confirmed vs toolkit)

| Type | Band | Payload |
|---|---|---|
| **IRA** (ring alert) | simplex (>1626.1) | beam position + TMSI paging — *best single-burst target* |
| **IMS** (messaging) | simplex | pager/SBD messaging |
| **ITL** (time/location) | simplex | 96-bit UW `11` + 94 zeros header |
| **IBC** (broadcast) | duplex (<1626.0) | sat/beam ID, system time, TMSI expiry, channel assignments |
| **LCW / IDA / VOC / VOD** | duplex | link-control, voice, SBD traffic |
| **IAQ** (acquisition) | duplex uplink | **handset access request** (BPSK, CRC-16) |
| **ISY** (sync) | either | all-sync burst |

### 1.6 IRA field layout (exact, confirmed — decoder output matches oracle)

After BCH, `bitstream_bch` (MSB-first):

| Field | Bits | Meaning |
|---|---|---|
| `ra_sat` | 0:7 | satellite ID (0–127, Iridium-internal) |
| `ra_cell` | 7:13 | spot-beam ID |
| `ra_pos_x` | 13:25 | signed 12-bit ECEF x (bit 13 = sign, ×(1<<11)) |
| `ra_pos_y` | 25:37 | signed 12-bit ECEF y |
| `ra_pos_z` | 37:49 | signed 12-bit ECEF z |
| `ra_int` | 49:56 | ring-alert interval (always 48 live) |
| `ra_ts` | 56:57 | broadcast slot (1 or 4) |
| `ra_eip` | 57:58 | 1 bit (EPI?) |
| `ra_bc_sb` | 58:63 | broadcast-channel sub-band |

Geometry: `lat = atan2(z, √(x²+y²))·180/π` (geocentric, ~0.2° low), `lon = atan2(y, x)·180/π`,
`alt_km = √(x²+y²+z²)·4` (display as `alt − 6378 + 23`). **Paging pages** (42 bits each, repeat until all-1s END):
`tmsi 0:32`, `zero1 32:34`, `msc_id 34:39`, `zero2 39:42`.

### 1.7 IBC field layout (exact, confirmed)

Header = 6 bits, BCH(poly 29) → `bc_type`. Then 4 × 42-bit blocks:
- **Block 1:** `sv_id 0:7`, `beam 7:13`, `unk 13:14`, `slot 14:15`, `sv_blocking 15:16`, `acqu_classes 16:32`,
  `acqu_subband 32:37`, `acqu_channels 37:40`, `unk 40:42`.
- **Block 2** (`type = bits 0:6`): type 1 → `iri_time = bits 10:42` (L-Band Frame Counter);
  type 2 → `tmsi_expiry 10:42`; type 0 → `max_uplink_pwr 36:42`.
- **Assignment blocks** (`type 0:3`): classic (type 0) = `random_id 3:11`, `timeslot 11:13 (+1)`, `ul_sb 13:18`,
  `dl_sb 18:23`, `access 23:26 (+1)`, `dtoa 26:34 (signed)`, `dfoa 34:40`.

Time conversion (ERA3 epoch, 2025-01-14 onward): `unix = iri_time·90/1000 + 1739556857`.

### 1.8 IAQ — handset acquisition (confirmed)

Uplink, duplex band, BPSK symbols only. Descrambled symbols map `{00→'0', 11→'1', 01/10→error}`. Fields:
`rid` (random id, bit-reordered from syms 4–11), 16-bit value, CRC: `iaq_crc16 = crcmod(poly=0x15101, init=0,
non-reflected, xorout=0)`, compared `>>2`. CRC-OK IAQ ⇒ **confirmed handset trying to access** — strong SIGINT signal.

### 1.9 Worked test vectors (generated by my independent encoder, decoded clean by iridium-toolkit oracle)

Encoder written from the spec notes above (no GPL code copied); oracle = iridium-toolkit (BSD-2, test only). These are
the *internal bit* payloads (before access code + DQPSK re-encode) — useful as C++/Python decoder unit vectors.

```
IRA  internal bits (224), freq 1626270833 Hz:
0000101011001010001011000010000110011111110000010111010010010000
0111110000011011000011000101100011101110100101011101100100110101
0010001100110000111011001010100011111111111111111111111111111111
11111111111111111111111111111111
→ oracle: IRA ... DL sat:016 beam:37 xyz=(+1396,+0711,+0862)
          pos=(+28.82/+026.99) alt=797 RAI:48 ?10 bc_sb:07
          P01: PAGE(tmsi:113dbaf8 msc_id:02) {OK}

IBC  internal bits (262), freq 1621698688 Hz:
0000001100110100101000001000110011101111000000110011111100110111
1111000000000100101110110010110010101011001100000010000000011011
0000000000100101101001000100110001000100100110000100110001001110
0001000010001100010000001000000000000000000000000000000000000000
010011
→ oracle: IBC ... DL bc:0 sat:028 cell:24 ... aq_sb:19 aq_ch:2
          time:2025-03-11T06:29:13.40Z [0 Rid:119 ts:1 ul_sb:22
          dl_sb:22 access:3 dtoa:+001 dfoa:00 00] []

IAQ  internal bits (52), uplink 1618000000 Hz:
0000000000111111001111110000001111001111001111001100
→ oracle: IAQ ... UL 0000 Rid:095 CRC:OK
```

The DQPSK-delta symbol strings (quarter-turns) for each are archived in the scratchpad script `iri_enc.py`; for a
synthetic IQ capture, map each internal bit-pair → QPSK symbol, differentially encode, prepend the DL/UL access code
(and 16/64-sym preamble), RRC-shape (α 0.4) at 25 kBd, place at the channel centre, add CFO/Doppler and noise.

### 1.10 Implementation recipe — Iridium

```
detector : FFT/energy burst finder in 1616–1626.5 MHz; narrow (~31.5 kHz) burst, 8.28 ms-ish.
           cheap trigger = x² carrier line (BPSK preamble) OR 25 kHz symbol line (reuse block FFT).
           gate on channel raster: snap centre to 1616 + 41.667·k kHz.
sync     : square→FFT CFO est; RRC+PLL→symbols; correlate 12-sym UW (DL & UL); dir = better match (dist ≤ 2).
decode   : strip access code → freq-gate frame type → BCH(1207/1897/29/41/465) per 32-bit block + parity →
           IRA/IBC/IAQ field extraction above.
fields   : IRA → sat, beam, ECEF→lat/lon/alt, TMSI pages; IBC → sat, beam, system time, assignments;
           IAQ → rid + CRC-16 OK.
validity : BCH all-blocks correctable + parity even; IAQ CRC-16 match. CRC/BCH-verified → drop "-like".
label    : "Iridium downlink (IRA/IBC)" or **"Iridium handset uplink (IAQ/UL burst)"** — the licensable event.
           Below sync threshold: "Iridium-like TDMA burst" (already partially in README).
```

---

## 2. THURAYA (GMR-1)

GEO (geostationary) MSS, GSM-derived. **Primary spec: ETSI GMR-1 family TS 101 376-x** (openly downloadable, see
sources; the spec text is public — implement from it, not from osmo-gmr code). Confirmed values below are quoted from
**TS 101 376-5-1 V1.2.1 (GMR-1 05.001)**.

### 2.1 PHY constants (confirmed, TS 101 376-5-1)

| Parameter | Value |
|---|---|
| Downlink band | 1525–1559 MHz |
| Uplink band (**handset**) | 1626.5–1660.5 MHz; **UL = DL + 101.5 MHz** (paired) |
| Carrier raster | **31.25 kHz**; carriers centred on integer multiples of 31.25 kHz; 1087 paired carriers; sub-band = 5 carriers |
| Symbol rate | **23.4 ksym/s** for *all* burst types (T = 1/23.4 ms) |
| Modulation | **π/4-CQPSK** (traffic, CCCH, BCCH); **π/4-CBPSK** (FACCH3, SDCCH); π/4-DBPSK (DKAB keep-alive); **6-PSK** (BACH); all RRC α = 0.35. **FCCH = real dual-chirp** spanning 3 slots |
| TDMA frame | **40 ms = 24 timeslots** |
| Multiframe | 640 ms = 16 frames; Superframe = 2.56 s = 64 frames; Hyperframe = 3 h 28 min 53.76 s |
| FEC | convolutional K=5 (tail 4 zero bits), some K=7 tail-biting; RS(15,9)/GF(2⁴) for BACH |

### 2.2 Burst types (confirmed, Table 6.1)

| Channel | Mod | Repeat (ms) | Info+CRC | Code rate | Coded bits | Interleave | Burst |
|---|---|---|---|---|---|---|---|
| FCCH | dual chirp | 320 | — | — | — | — | FCCH (3-slot chirp) |
| BCCH | π/4-CQPSK | 320 | 192+16 | 1/2 | 424 | 53 | BCCH (424 bits) |
| GBCH | π/4-CQPSK | 80 | 108+16 | 1/2 | 264 | 32 | DC2 (132) |
| PCH/AGCH/CBCH | π/4-CQPSK | 640/var | 192+16 (184+16 CBCH) | 1/2 | 424 | 53 | DC6 (432) |
| RACH | π/4-CQPSK | — | class I 16+8 / II 123+12 | 1/4 | 494 | 14/33 | RACH (494) |
| TCH3/6/9 | π/4-CQPSK | 40 | 80 / 2.4–9.6 kbps | var | 104/420/648 | 24/52/81 | NT3/6/9 |
| SDCCH | π/4-CBPSK | up 320 | 84+16 | 1/4 | 416 | 52 | SDCCH (208) |

### 2.3 What is detectable / decodable in 34 ms

- **FCCH (dual real chirp over 3 slots = 3·(40/24) ms ≈ 5 ms):** a chirp is the *cheapest, most distinctive* GMR-1
  signature — a sweeping tone, not a steady carrier. A matched-filter / dechirp (multiply by conj-chirp → tone) gives
  a clean detect inside one block and is unlike anything else in L-band. **This is the recommended Thuraya detector.**
  (Exact chirp rate/bandwidth not in 05.001 §7 — in TS 101 376-5-2/-5-4; marked *unconfirmed* until read.)
- **BCCH burst (424 coded bits ≈ 424/2/23.4k ≈ 9 ms):** one BCCH burst fits in a block. After CFO/timing it is
  π/4-CQPSK 23.4 kBd; demod + deinterleave(53) + K=5 rate-1/2 Viterbi + CRC-16 → **BCCH system information** (spot-beam
  ID, satellite position, RACH offset, LAI). This yields operator/beam identity. Full SI needs multiframe (several
  blocks) — single-block = one burst, partial.
- **Uplink RACH (494 bits) / TCH:** π/4-CQPSK bursts in 1626.5–1660.5 MHz ⇒ **handset transmitting**. A RACH or a
  steady TCH3/6/9 burst train on the 31.25 kHz raster in the uplink band is the licensable event.

### 2.4 Implementation recipe — Thuraya/GMR-1

```
detector : (a) FCCH dechirp in 1525–1559 (DL) — strongest tell; (b) π/4-CQPSK 23.4 kBd burst on 31.25 kHz raster.
           Uplink band 1626.5–1660.5 MHz burst ⇒ "GMR-1 handset uplink-like".
sync     : FCCH → coarse freq/time; then find BCCH on an associated carrier.
decode   : π/4-CQPSK demod (RRC 0.35, 23.4 kBd) → deinterleave(depth per channel) → Viterbi K=5 r=1/2 →
           CRC-16 → BCCH/RACH fields.
fields   : BCCH → spot beam, satellite ephemeris/position, LAI, RACH/Doppler offsets.
validity : CRC-16 on 192+16 BCCH block. Verified → "Thuraya/GMR-1 BCCH"; else "GMR-1-like π/4-CQPSK burst".
label    : "Thuraya/GMR-1 downlink" / **"GMR-1 handset uplink-like"**.
```

> **Licence note:** `osmo-gmr` (osmocom) implements GMR-1 L1 (FCCH sync, π/4-C(B)QPSK demod, conv-coding, BCCH
> analyzer) but is **GPL** — read as a learning oracle only, never copy into `c++/`. The ETSI specs are the free,
> licence-clean source; implement from them.

---

## 3. INMARSAT

GEO MSS. Several distinct air interfaces share the band (DL 1525–1559, UL 1626.5–1660.5 MHz). Practical decodability
varies widely.

| Service | Air interface | Decodable in 34 ms? |
|---|---|---|
| **IsatPhone Pro / 2** (handset) | GMR-2+ family (Inmarsat's variant; *not* identical to Thuraya GMR-1) — **unconfirmed** which sub-release | Detect π/4-shifted PSK bursts in UL band; decode = research-hard. Treat as **detect/classify only**, "IsatPhone-like" |
| **Inmarsat-C** (data/telex, EGC SafetyNET) | 1200 sym/s (D)BPSK, TDM NCS channel | **Yes, partial** — see §3.1; most tractable Inmarsat target |
| **Inmarsat Aero** | Classic Aero 600/1200/10500 bps (AMBE voice + ACARS-like) | detect/classify; full decode out of scope |
| **BGAN / IsatData Pro** | TDM/TDMA, QPSK/16-APSK bearer-adaptive, OFDM-like framing | detect/classify only |

### 3.1 Inmarsat-C TDM (NCS common channel) — confirmed enough to attempt

(Confirmed: sigidwiki Inmarsat-C_TDM, Wavecom SAT-C-TDM decoder notes, patent US5774462.)

| Parameter | Value |
|---|---|
| Modulation | **1200 sym/s (D)BPSK** (sources differ BPSK vs DBPSK — try both) |
| NCS common channel (NCSC) downlink | per ocean region: AOR-W 1537.70, IOR 1537.10 MHz (TDM band 1537.0–1545.0; sigidwiki quotes 1537.70–1539.685) |
| Terminal uplink | 1626.5–1646.5 MHz |
| **TDM frame** | **8.64 s = 64 slots × 162 symbols**; each slot = **2 interleaved unique-word symbols + 160 data symbols** (patent US5774462) |
| FEC | **rate-1/2 convolutional** + interleaving + scrambling (UW spread across the whole frame — must receive a full frame to sync) |

**34-ms reality:** 34 ms ≈ 41 symbols at 1200 sym/s — **far less than one 162-symbol slot, let alone an 8.64 s frame.**
So in one Pluto block Inmarsat-C is **detect/classify only** (a steady ~1.2 kBd BPSK carrier at a known NCSC
frequency). Decoding needs a long dwell (several seconds) via `StreamAnalyzer`, not a single block. EGC SafetyNET
content rides this TDM; decoding it requires the full-frame receiver.

### 3.2 Implementation recipe — Inmarsat

```
detector : narrowband BPSK carrier at NCSC freqs (1537.0–1539.7 MHz) ~1.2 kBd → "Inmarsat-C TDM-like".
           UL band (1626.5–1660.5) π/4-PSK burst ⇒ "Inmarsat IsatPhone handset uplink-like".
sync     : (C only, long dwell) acquire 1200-baud BPSK, find interleaved UW over a full 8.64 s frame.
decode   : (C, long dwell) deinterleave → Viterbi r=1/2 → descramble → EGC/message frames.
fields   : Inmarsat-C: message type, EGC SafetyNET addressing (if full-frame).
validity : frame CRC / UW consistency. Single block = classify only.
label    : "Inmarsat-C TDM-like" / **"IsatPhone handset uplink-like"** / "BGAN-like" (detect only).
```

---

## 4. GLOBALSTAR (CDMA phones) and 5. GLOBALSTAR SIMPLEX / SPOT

### 4. Globalstar CDMA (voice/data phones — GSP-1600/1700)

(Confirmed: FCC orders DOC-204216A1 / FCC-08-248A1; air interface is IS-95-derived.)

| Parameter | Value |
|---|---|
| Uplink (**handset**) | 1610–1618.725 MHz (1610–1617.775 exclusive; 1617.775–1618.725 **shared with Iridium**) |
| Downlink | 2483.5–2500 MHz (S-band) |
| Channel | 1.23 MHz wide, **QPSK CDMA, 1.2288 Mcps** (IS-95-like) |
| Discriminator vs Iridium | Globalstar UL = **wideband (~1.23 MHz) spread CDMA**; Iridium = **narrow (~31.5 kHz) TDMA burst**. In the 1616–1618.7 MHz overlap, bandwidth alone separates them. |

**34-ms reality:** a 1.23 MHz CDMA channel fits in-band; without the PN/Walsh codes you cannot despread/decode (not
public). So: **detect a ~1.23 MHz noise-like spread signal in 1610–1618.725 MHz, measure chip rate ≈ 1.2288 Mcps
(cyclostationary feature / 1.23 MHz occupied BW) → "Globalstar CDMA uplink-like"**. That alone = a handset in use.
GNSS at 1575.42 (GPS L1) / 1602 (GLONASS) and the 1610–1626.5 window are distinct; RNSS 1610.6–1613.8 (radio-astronomy)
is quiet.

```
detector : wide (~1.23 MHz) spread/noise-like energy in 1610–1618.725 MHz; chip-rate feature at 1.2288 Mcps.
validity : can't decode (codes proprietary). label "-like" always.
label    : **"Globalstar CDMA handset uplink-like"**.
```

### 5. Globalstar simplex / SPOT / SmartOne (trackers & PLBs) — FULLY decodable

(Confirmed & **verified by me**: Colby Moore, Black Hat USA 2015 white paper + PoC‖GTFO 0x09 article. I re-derived the
PN from the LFSR params, reproduced the published PN, and verified the CRC-24 and GPS decode on the published sample
packet.)

| Parameter | Value (confirmed) |
|---|---|
| Uplink channels | A 1611.25, B 1613.75, C 1616.25, **D 1618.78** MHz (note: source says **1618.78**, not 1618.75; mark D *unconfirmed*). Live devices mostly channel A. |
| Downlink | 6875–7055 MHz (C-band, bent-pipe) — not L-band |
| Modulation | **BPSK, non-differential**, 100.04 bit/s data |
| Spreading | **DSSS, 255-chip m-sequence, 1.25 Mcps**; PN repeats **49× per data bit**, never crosses a bit boundary (so despread ≡ decode as fast BPSK) |
| PN LFSR | m-sequence, degree 8, **Galois LFSR, mask 166 (dec), seed 59** — *verified: reproduces the published 255-chip PN exactly*. Equivalent Fibonacci recurrence tap mask `0b1100101`. |
| Polarisation | LHCP |
| Packet | **144 bits**, duration ≈ **1.44 s** (144 / 100.04); sent as multiple transmissions for reliability |

**Packet layout (144 bits, confirmed & CRC-verified):**

| Field | Bits | Notes |
|---|---|---|
| Preamble | 0:10 | `0000001011` |
| ESN | 10:36 | 3-bit manufacturer ID + 23-bit unit ID |
| Message # | 36:40 | mod-16, NV-stored |
| Packet # | 40:44 | # packets in message |
| Packet seq # | 44:48 | sequence within message |
| User data | 48:120 | 72 bits (9 bytes, MSB first) |
| CRC-24 | 120:144 | poly **0o114377431 = 0x131FF19**, init 0xFFFFFF, bit-reflected input per the published algorithm, final `~crc & 0xFFFFFF` |

CRC algorithm (verified against sample → `0xc8293`): feed 14 bytes starting at the ESN (skip first 2 preamble bits of
byte 0 via `&0x3f`), MSB-first, XOR each byte `<<16`, 8 shifts, reduce when bit 24 set.

**GPS decode (verified):** `lat = int(user_data[8:32],2)·90/2²³`; `lon = 360 − int(user_data[32:56],2)·180/2²³`
(i.e. 0–360° west-referenced).

**34-ms reality:** the 144-bit packet lasts ~1.44 s ⇒ **a single Pluto block holds only ~2.4 % of one packet** (~34 ms
= ~3.4 data bits spread). So in one block SPOT is **detect only** (despread-correlate the known 255-chip PN at
1.25 Mcps → correlation peak). Full packet decode (ESN, GPS) needs a ≥1.5 s dwell via `StreamAnalyzer`. The PN
correlation is a *very* robust detector even at negative SNR (DSSS processing gain ≈ 24 dB).

```
detector : correlate the known 255-chip PN @1.25 Mcps over the input at 1611.25/1613.75/1616.25/1618.78 MHz.
           sharp autocorrelation peak (off-peak = −1) → "Globalstar-simplex/SPOT-like". cheap & SNR-robust.
decode   : (≥1.5 s dwell) track PN phase, despread → BPSK 100.04 bps → 144-bit packet → CRC-24 verify → ESN + GPS.
fields   : ESN (mfr+unit), message/packet/seq #, lat/lon (if GPS-format user data).
validity : CRC-24 match → drop "-like"; expose ESN as device identity.
label    : "Globalstar simplex (SPOT/SmartOne) tracker" + ESN; single block → "-like" (detect only).
```

**SPOT test vector (verified by me):** the 144-bit packet from the PoC‖GTFO article decodes to
ESN mfr=1/unit=2541472 (0xA6C7A0), CRC `0xc8293` (matches field), lat 0.709°, lon −0.022° — use as the Python/C++
decoder unit vector. Full bit string in scratchpad `spot.py`.

---

## 6. LEGAL STATUS IN INDIA (why uplinks matter)

- **General rule:** operating *any* satellite phone in India requires a **prior licence from the Department of
  Telecommunications (DoT)**; unlicensed use is an offence under **Section 6 of the Indian Wireless Telegraphy Act 1933**
  and **Section 20 of the Indian Telegraph Act 1885**. Customs clear handsets only against a DoT licence.
- **Thuraya: banned** for general use (consistently stated across Indian maritime circulars and the Indian embassy
  notices).
- **Iridium: restricted.** Historically treated as banned for general use; since a **2023 DGS order**, Iridium-based
  equipment is permitted **only for GMDSS distress & safety** at sea. General/handset use still needs DoT permission.
- **Inmarsat: the permitted route** — DoT permissions are largely restricted to Inmarsat; other sat equipment is
  barred within 12 nm of the Indian coast.
- **Globalstar / SPOT:** not separately called out, but falls under the same DoT licensing regime.

**SIGINT implication:** a **handset *uplink*** (Iridium IAQ/UL burst, GMR-1 RACH/TCH in 1626.5–1660.5, Globalstar CDMA
1610–1618.7, SPOT simplex) is the on-air evidence that an unlicensed sat device is *transmitting* from the monitored
area — the reportable event. Downlinks only show a satellite is serving/paging the area.

Sources: [1][18][19][20].

---

## 7. Cross-system detector / disambiguation cheat-sheet

| Observation | Likely system |
|---|---|
| Narrow (~31.5 kHz) TDMA burst, 25 kBd DQPSK, 1616–1626.5 MHz, on 41.667 kHz raster | **Iridium** (UW decides UL vs DL) |
| Real chirp sweep in 1525–1559 MHz | **Thuraya/GMR-1 FCCH** |
| π/4-CQPSK 23.4 kBd burst on 31.25 kHz raster, 1626.5–1660.5 MHz | **GMR-1 handset uplink** (Thuraya) |
| Steady ~1.2 kBd BPSK at 1537.0–1539.7 MHz | **Inmarsat-C TDM** |
| Wide ~1.23 MHz spread/noise-like, 1610–1618.7 MHz | **Globalstar CDMA uplink** |
| 255-chip PN correlation @1.25 Mcps, 1611.25/1613.75/1616.25 MHz | **Globalstar simplex / SPOT** |
| CW/BPSK at 1575.42 / 1602 / 1176 / 1610.6–1613.8 | GNSS / RNSS — **not** a sat phone (see 02_gnss.md) |

Overlap resolver (1616–1618.725 MHz, Iridium ∩ Globalstar): **bandwidth** — Iridium ≈ 31.5 kHz TDMA burst vs
Globalstar ≈ 1.23 MHz CDMA vs SPOT DSSS PN-correlatable. Reuse the block FFT / Channelizer for the occupied-BW measure
(latency goal), gate PN/UW correlation behind the cheap energy+BW trigger.

---

## 8. Synthetic-signal generation notes (for tests)

- **Iridium:** map internal bit-pairs → QPSK syms → differential encode → prepend 16/64-sym preamble + DL/UL access
  code → RRC α0.4 @25 kBd → shift to channel centre (1616+41.667·k kHz) → add CFO (±20 kHz ppm + ±36 kHz Doppler) +
  AWGN. Use the §1.9 vectors; assert decoder output == oracle lines.
- **GMR-1:** π/4-CQPSK/CBPSK @23.4 kBd, RRC α0.35; conv-encode K=5 r=1/2, interleave (depth per channel), CRC-16. FCCH:
  synthesise a dual real chirp over 3 slots (chirp params *unconfirmed* — read TS 101 376-5-2/-5-4 first). Verify
  against osmo-gmr (GPL, oracle only) once chirp params are confirmed.
- **Inmarsat-C:** 1200-baud BPSK, 162-sym slots (2 UW + 160 data), r=1/2 conv + interleave + scramble; need a
  multi-second synthetic stream to exercise the full-frame receiver.
- **SPOT:** BPSK 100.04 bps → spread with the verified 255-chip PN @1.25 Mcps (49 PN/bit) → carrier → add CFO + AWGN;
  assert PN correlation peak and (full-length) CRC-24/ESN/GPS == §5 vector.

---

## 9. What is implementable now vs blocked

**Implementable from confirmed data (synthetic-validated):**
- Iridium burst detect + UW UL/DL discrimination + IRA/IBC/IAQ decode (BCH + CRC verified against oracle). Full
  recipe, exact bit layouts, test vectors in hand.
- Globalstar simplex / SPOT: PN detector (single block) + full packet decode (ESN, GPS, CRC-24) on a dwell — all
  constants verified by me.
- GMR-1/Thuraya: π/4-CQPSK burst + FCCH-chirp detector + BCCH framing from the ETSI spec (FEC/interleave depths
  confirmed).
- Inmarsat-C detect/classify; Globalstar CDMA detect/classify.

**Blocked / unconfirmed:**
- GMR-1 FCCH exact **chirp rate & bandwidth** (need TS 101 376-5-2/-5-4). BCCH bit-level SI element layout (need
  TS 101 376-4-x RR/MM).
- IsatPhone exact air-interface sub-release (GMR-2+ vs GMR-1 3G) — mark "-like", detect only.
- Globalstar CDMA & BGAN payloads: proprietary codes → never decodable, "-like" only.
- SPOT channel D = **1618.78** (source) vs 1618.75 (brief) — verify on capture.
- All labels stay `*-like*` until a real India capture (promotion to `proven` per batch rules).

---

## Sources

1. iridium-toolkit — muccc. https://github.com/muccc/iridium-toolkit — **Licence: BSD 2-Clause** (README: "(c) Sec &
   schneider … 2-Clause BSD"). *Permissive — usable as reference & test oracle; I used it only as a decode oracle for
   the §1.9 vectors.* Files read: `bitsparser.py`, `bch.py`, `util.py`, `FORMAT.md`, `iridium.txt`, frame-structure SVG,
   `tests/test_parser.py`.
2. gr-iridium — muccc. https://github.com/muccc/gr-iridium — **Licence: GPL-3** (source headers "GNU General Public
   License … version 3"). *GPL — read-only oracle; do NOT port into c++/.* Files read: `lib/iridium.h`,
   `burst_downmix_impl.cc`, `iridium_qpsk_demod_impl.cc`, `iridium_frame_printer_impl.cc`, `docs/README.iridium`.
3. Signal Identification Wiki — Iridium. https://sigidwiki.com/wiki/Iridium (raster, simplex channel freqs, 25 kBd,
   8.28/90 ms). Community-maintained.
4. US Patent 5,596,315 (Motorola) — Iridium message unit / frame. https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/5596315 — public patent.
5. US Patent 5,603,079 — simplex messaging in guard band (90 ms frame, 20.48 ms simplex). https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/5603079 — public patent.
6. ETSI TS 101 376-5-1 V1.2.1 (GMR-1 05.001, physical-layer general description). https://www.etsi.org/deliver/etsi_ts/101300_101399/1013760501/01.02.01_60/ts_1013760501v010201p.pdf — **ETSI, freely downloadable; spec text implementable (not GPL).** (Verified band, 31.25 kHz raster, 23.4 kBd, 40 ms frame, burst table.)
7. ETSI GMR-1 family index (TS 101 376-5-2 burst/mux, -5-4 modulation, -4-x RR/MM — *to read for FCCH chirp + SI
   layouts*). https://www.etsi.org/standards (search "GMR-1").
8. osmo-gmr — Osmocom. https://osmocom.org/projects/osmo-gmr / https://gitea.osmocom.org/sdr/osmo-gmr — **Licence: GPL**
   (Osmocom). *Learning/oracle only; implement from ETSI specs.* (28C3 talk "Introducing osmo-gmr" for FCCH/π/4-CQPSK.)
9. Colby Moore, "Spread Spectrum Satcom Hacking: Attacking the GlobalStar Simplex Data Service," Black Hat USA 2015 —
   white paper. https://www.blackhat.com/docs/us-15/materials/us-15-Moore-Spread-Spectrum-Satcom-Hacking-Attacking-The-GlobalStar-Simplex-Data-Service-wp.pdf (PN seq, LFSR mask 166/seed 59, 144-bit packet).
10. Same, slides. https://www.blackhat.com/docs/us-15/materials/us-15-Moore-Spread-Spectrum-Satcom-Hacking-Attacking-The-GlobalStar-Simplex-Data-Service.pdf
11. Colby Moore, "Breaking Globalstar Satellite Communications," PoC‖GTFO 0x09 §3. https://mcfp.felk.cvut.cz/publicDatasets/pocorgtfo/contents/articles/09-03.pdf (channels, 100.04 bps, 49 PN/bit, packet layout, CRC-24 poly 0o114377431, GPS decode). PoC code: https://github.com/synack/globalstar — *repo now offline; article text used.* Licence of article/code unspecified → treated as reference; my encoder/CRC/PN re-derived independently and verified.
12. sigidwiki — Inmarsat-C TDM. https://sigidwiki.com/wiki/Inmarsat-C_TDM (1200 sym/s, NCSC freqs).
13. Wavecom decoder help — SAT-C-TDM / SAT-C-TDMA. https://www.wavecom.ch/content/ext/DecoderOnlineHelp/worddocuments/satctdmsatctdma1.htm (1200 sym/s, 8.64 s frame, scrambled/encoded/interleaved).
14. US Patent 5,774,462 — Inmarsat unique-word recognition (64 slots × 162 sym, 2 UW + 160 data, 8.64 s). https://patents.google.com/patent/US5774462
15. FCC Order DOC-204216A1 & FCC-08-248A1 — Globalstar/Iridium L-band allocation (1610–1618.725 shared, 2483.5–2500,
    1.23 Mc/s QPSK CDMA). https://docs.fcc.gov/public/attachments/DOC-204216A1.pdf , https://docs.fcc.gov/public/attachments/FCC-08-248A1.pdf
16. "Signal design and system operation of Globalstar versus IS-95 CDMA," Wireless Networks 6 (2000). https://doi.org/10.1023/A:1019164909332 (IS-95-derived air interface) — *abstract only read; mark detailed channelisation unconfirmed.*
17. Globalstar STX-3 simplex module (datasheet/manual) — product pages. https://www.outfittersatellite.com/products/globalstar-stx-3-simplex-modem-1 (burst/repeat scheduling).
18. Embassy of India, Tokyo — "Restrictions on use of satellite phones in India." https://www.indembassy-tokyo.gov.in/eoityo_listview/Mzg3 (DoT licence required; Thuraya/Iridium banned).
19. Indian maritime circular (Adani/Kattupalli) — "Usage of Thuraya, Iridium … in Indian Waters" (Sec 6 Wireless Act / Sec 20 Telegraph Act). https://www.adaniports.com/-/media/Project/Ports/PortsAndTerminals/Kattupalli-Port-Documents/Public-Notices/22-04-14-TC-03-Usage-of-Thuraya-Iridium-and-other-satellite-communication-in-Indian-Waters.pdf
20. DGS Order (2012/2023) on Thuraya/Iridium & GMDSS distress use. https://betadgs.dgshipping.gov.in/ (2023 order permits Iridium for GMDSS distress only) — *exact 2023 order text unconfirmed; verify against primary DGS/DoT notification.*
