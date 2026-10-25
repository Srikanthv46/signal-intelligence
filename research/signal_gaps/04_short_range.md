# 04 — Short-range: 315/433 MHz OOK/FSK (EV1527, PT2262, KeeLoq, TPMS), Bluetooth BR/EDR, BLE trackers

Batch S, items **S6** (315/433 MHz OOK/FSK devices) and **S7** (Bluetooth BR/EDR LAP + BLE trackers).
Engineering reference built from datasheets, the Bluetooth Core Specification, vendor/IETF/academic documents, and
GPL/AGPL decoders used **only as fact sources / oracles** (never copied; C++ stays GPL-free). Every source is numbered
in §Sources with its licence. Items not confirmed against a datasheet or a real capture are marked **[UNCONFIRMED]**.

Scope reminder: input is an SDR IQ block of 2^21 complex samples. At 61.44 MSps that is 34 ms (too short for most of
these sub-GHz bursts — see §0.3); the practical path for 315/433 MHz is a **narrowband dwell at 2–4 MSps**, where one
block is 0.5–1.0 s and holds the full multi-frame remote/TPMS burst. Bluetooth is 2.4 GHz, so a 34 ms / 61.44 MSps
block (56 MHz usable) is the right input there.

---

## 0. Pulse-slicer front end (shared by all sub-GHz OOK/ASK/FSK protocols)

All the 315/433 protocols in §1 are **pulse-width / pulse-position coded**. One generic slicer feeds every decoder.

### 0.1 Pipeline
1. **Band/trigger.** These devices sit in the Indian SRD band **433.05–434.79 MHz** (ISM, license-free, max 10 mW ERP
   typical) and on **315 MHz** for some imported cars/fobs (315 MHz is *not* a legal SRD band in India, but imported
   vehicle fobs/TPMS transmit there). Gate on a narrowband carrier/energy burst in those ranges.
2. **Envelope (OOK/ASK).** `env = |x|` (or `|x|^2`), low-pass at ~2× the fastest expected symbol rate, then an
   adaptive threshold (e.g. midway between a sliding noise floor and peak, or Otsu on the burst). Output = run-length
   list of (mark duration, space duration) in microseconds.
3. **FSK path (TPMS, KeeLoq FSK variants).** For 2-FSK, take the instantaneous frequency `d/dt angle(x)` and slice its
   sign → a bit stream at the known data rate; then line-decode (Manchester / differential Manchester / PCM). Do NOT
   use the amplitude envelope for FSK.
4. **Line decode / PWM classify.** Convert runs to symbols by comparing against the protocol's short/long reference:
   - **PWM (OOK):** each bit = one mark+space pair; `short-mark/long-space` = one logic level, `long-mark/short-space`
     = the other (EV1527, PT2262, KeeLoq).
   - **PPM:** fixed mark, data in the gap length.
   - **Manchester:** 01→0 or 10→1 (polarity per protocol); half-bit = the short reference.
   - **Differential Manchester:** transition-at-start = one bit, no-transition = the other (Toyota/Pacific TPMS).
5. **Framing.** Find preamble/sync, slice the payload to the fixed bit count, run the integrity check (CRC/checksum/
   parity), map fields, assign label + confidence.

### 0.2 Reference timing table (units: microseconds unless noted)

| Protocol | Mod | "short" | "long" | sync/gap | bits | line code | integrity | Src |
|---|---|---|---|---|---|---|---|---|
| EV1527/HS1527/RT1527 | OOK/ASK PWM | 1 Te mark | 3 Te mark | sync = 1 Te high + 31 Te low | 24 (20 ID + 4 data) | PWM (bit0=HLLL, bit1=HHHL in Te) | none (learning code) | [1][2][14] |
| PT2262/SC2262 | OOK/ASK PWM | 4α | 12α (1 bit = 32α) | sync = 128α | 12 tri-state AD | PWM, 3 levels (0/1/F) | none (fixed) | [3][4][14] |
| PT2260 | OOK/ASK PWM | (1 bit = 1024α; see note) | — | sync 4096α | 12 tri-state | PWM | none | [15] |
| KeeLoq HCS200/300/301 | OOK PWM (HCS301 also FSK opt.) | 1 Te (≈400 µs typ) | bit = 3 Te (TBP≈1200) | preamble 23 Te, header 10 Te low, guard 39 Te | 66 (LSB first) | PWM (logic0=1Te h/2Te l; logic1=2Te h/1Te l) | none in clear (32-bit part encrypted) | [5][6][7][14] |
| TPMS (generic) | FSK or ASK, Manchester/diff-Manchester, PCM | ~50 µs/bit (FSK PCM) | — | per-sensor preamble | 64–160 | varies | CRC-8 / CRC-16 / sum | §1.3 |

Te = KeeLoq basic pulse element; α = PT2262 oscillator unit. KeeLoq Te typ 400 µs (range 280–620 µs); so one 66-bit
HCS301 word ≈ preamble 23 Te + header 10 Te + 66×3 Te ≈ (23+10+198)×400 µs ≈ **92 ms**, repeated ≥2× with a 17.5 ms
inter-packet gap → fits comfortably inside one 2–4 MSps dwell block.

### 0.3 What fits in one block / repetition

- **61.44 MSps block = 34 ms:** too short — an EV1527 frame alone (24 PWM bits at ~400 µs Te, bit=4 Te) is ~38 ms, a
  KeeLoq word ~92 ms, a TPMS burst 4–11 frames over tens of ms to a few seconds. **Flag `SHORT_BLOCK` and do not try to
  decode sub-GHz remotes/TPMS from a 34 ms block.** Detect-only (narrowband OOK/FSK burst present) is acceptable there.
- **2–4 MSps dwell block = 0.5–1.0 s:** holds the full multi-frame burst. This is the development/decode input for S6.
  Sample rate needs only ≥ ~4× the fastest line rate (≥ ~200 kSps for 50 µs TPMS bits; ≥ ~20 kSps for 400 µs remotes),
  so 2 MSps is ample.
- **Repetition:** remotes repeat the identical frame **4–20×** (EV1527 fobs often ~4–8; KeeLoq sends ≥2 with a 17.5 ms
  gap while the button is held). TPMS sends **bursts of 4–11 near-identical frames** every ~60 s (or on rapid pressure
  change / above ~40 km/h). Use repetition as a confidence booster: require ≥2 identical CRC-valid frames before a
  non-"-like" label; a lone frame → "-like".

---

## 1. 315 / 433.92 MHz devices

### 1.1 EV1527 / HS1527 / RT1527 (learning-code remotes) — IMPLEMENTABLE

Common in Indian gate/shutter remotes, wireless doorbells, cheap alarm PIRs, RF sockets. One-time-programmable, 24-bit
frame, **no checksum** (so validity rests on timing consistency + repetition + the fixed sync).

**Frame (sent MSB-first):** `[sync] C19 C18 … C0  D3 D2 D1 D0`
- **Sync bit:** 1 pulse of 1 Te high followed by **31 Te** low (long preamble gap). This long low is the frame
  delimiter — detect it to align.
- **20-bit code (C19..C0):** factory-random unique ID (2^20 ≈ 1.05M codes).
- **4-bit data (D3..D0):** the key/button state (mapped from input pins K0–K3).
- **Bit encoding (PWM, Te = oscillator period, ~1/80 kHz ≈ 12.5 µs internal but the *output* Te is set by Rosc; field
  Te commonly 300–500 µs):** logic **0 = 1 Te high + 3 Te low**; logic **1 = 3 Te high + 1 Te low**. Each bit = 4 Te.
  [done.land [14] quotes 275 µs high + 275 µs low for '0' and 275 high + 1225 low for '1' — same 1:3 / 3:1 ratio,
  Te≈275 µs.]

**Detector → slice → decode → validity → fields → label**
- Detector: OOK burst in 433.05–434.79 (or 315) MHz; find the 1-high/31-low sync gap; measure Te from the first few
  marks (median short mark).
- Slice: for each of 24 bits, ratio mark/(mark+space) ≈ 0.25 → '0', ≈ 0.75 → '1'. Reject bits whose ratio is ambiguous.
- Validity: (a) exactly 24 data bits between two syncs; (b) Te stable within ±30 %; (c) **≥2 identical repeats** of the
  same 24-bit word. No CRC → without ≥2 repeats, emit "-like" at ≤0.5 confidence (`unknown` per the never-guess rule).
- Fields: `id` = C19..C0 (hex, 5 nibbles), `data` = D3..D0, `repeats`, `te_us`, `rate`.
- Label: `EV1527/HS1527-like (OOK remote)`; promote to confirmed only on a real capture with a known fob.

**Synthetic generation for tests:** pick Te (e.g. 350 µs), random 20-bit ID + 4-bit data, render sync (1 Te / 31 Te)
then 24 PWM bits, repeat 6×, OOK-modulate onto a 433.92 MHz complex carrier at 2 MSps, add AWGN. Unit test: slicer
recovers the ID/data for SNR ≥ ~6 dB.

### 1.2 PT2262 / PT2260 / SC2262 tri-state fixed codes — IMPLEMENTABLE

Older fixed-code remotes/alarms; 12 **tri-state** address/data bits. Each of the 12 "AD" positions is 0, 1, or F
(floating). rtl_433 `generic_remote` [1] handles PT2260/PT2262/SC2260/SC2262/EV1527 as one family.

**PT2262 bit timing (α = oscillator unit, set by Rosc):** `1/8 bit = 4α`, `1 bit = 32α`, `4 bits = 128α`, **sync bit =
128α** (a short high + long low that delimits each 12-bit code word) [3][4].
**Tri-state symbol (two PWM half-bits per tri-state bit):**
- `0` = two "short" elements: `(4α high, 12α low)(4α high, 12α low)` → bits `00`
- `1` = two "long" elements: `(12α high, 4α low)(12α high, 4α low)` → bits `11`
- `F` (float) = one short + one long → bits `01`
(rtl_433 reads 25 PWM bits, inverts, then maps each 2-bit group: `00→'0'`, `11→'1'`, `01→'Z/open'`, `10→'X'` invalid
for SC226x but valid in EV1527 [1].)

**PT2260** uses much larger multiples (`1/8 bit = 128α`, `1 bit = 1024α`, `sync = 4096α`) [15] — do not mix with PT2262.

**Detector → … → label**
- Slice 12 tri-state symbols (24 PWM half-bits) + sync; map to the `0/1/F` string.
- Fields: `tristate` (12-char string over {0,1,F}), `address` (first 8), `data` (last 4), `code` (int), `alpha_us`.
- Validity: timing consistency + last bit always 1 (rtl_433 rejects `(b[3]&0x80)==0`) + reject all-zero ID/cmd to cut
  false positives [1]; require ≥2 repeats. No checksum → "-like".
- Label: `PT2262/SC2262-like (tri-state fixed-code remote)`.

**Synthetic:** encode a random 12-symbol tri-state word at α≈100–300 µs, 4–8 repeats, OOK at 433.92 MHz.

### 1.3 Microchip KeeLoq HCS200/HCS300/HCS301 — PARTIAL (clear fields only), IMPLEMENTABLE

Code-hopping garage/gate/car fobs. **Decode only the clear parts (serial, buttons, status); never attempt key recovery
or decrypt the 32-bit hopping code** (per task and CLAUDE "never guess"; the encrypted block is cryptographically
opaque without the manufacturer key).

**Modulation:** PWM OOK (HCS301 can also be FSK). rtl_433 `hcs200` uses `OOK_PULSE_PWM, s=370, l=772, gap=1500,
reset=9000, t=152` (µs) [6], i.e. Te≈386 µs, bit≈1.16 ms. Warm-up = **12 short pulses then a 4400 µs gap**; then the
66-bit word; **two packets with a 17500 µs gap** between them [6]. Datasheet timing: TE typ 400 µs (280–620), preamble
23 TE (≈9.2 ms typ), header 10 TE (≈4.0 ms low), PWM bit = 3 TE (TBP≈1200 µs), guard 39 TE [5][7].

**Bit encoding (PWM):** logic `0` = `1 TE high + 2 TE low`; logic `1` = `2 TE high + 1 TE low` (each bit 3 TE)
[5][7][14]. **Transmitted LSB-first.**

**66-bit code word layout (bit index, LSB first)** [6][5]:

| Bits | Field | Notes |
|---|---|---|
| 0–31 | Encrypted hopping code | **Do not decrypt.** Report as opaque 32-bit hex only. |
| 32–59 | Serial number (28-bit) | device identity (clear) |
| 60–63 | Button status | order S3,S0,S1,S2 (first→last). rtl_433: `btn_num = S3|S0<<2|S1|S2>>2`; `0x0F` = learn [6] |
| 64 | Battery low | 1 = low |
| 65 | Repeat / fixed | HCS200 row uses bit as repeat; datasheet calls last bit "fixed 1"/VLOW/RPT per part |

**Detector → slice → decode → validity → fields → label**
- Detector: OOK burst; find the 12-pulse warm-up + 4.4 ms gap, or the preamble/header; measure TE.
- Slice 66 PWM bits (LSB first). Preamble sanity: rtl_433 checks first row = `0xfff` (12 ones) [6].
- Validity: **length == 66**, preamble present, not all-0xFF, and ideally two packets (17.5 ms gap) matching serial +
  button. No CRC on the clear part → require the two-packet consistency for a confirmed-vs-"-like" call.
- Fields: `serial` (28-bit hex), `button`, `battery_ok`, `repeat`, `encrypted_hop` (32-bit hex, opaque), `te_us`.
- Label: `KeeLoq HCS200/300/301-like (code-hopping remote)`; note "hopping code encrypted, not decoded".

**Synthetic:** random 28-bit serial + button + random 32-bit "encrypted" block (we never verify it), render warm-up +
preamble + 66 PWM bits, two packets, OOK at 433.92 MHz. Test: serial/button/battery recovered for SNR ≥ ~6 dB.

### 1.4 TPMS (Tire Pressure Monitoring) — IMPLEMENTABLE (4 well-documented protocols below)

India mandated TPMS later than EU/US and it is largely **aftermarket or OEM on mid/upper trims**; both exist on Indian
roads. **315 MHz** sensors appear on imported/older US/JDM cars; **433.92 MHz** dominates EU-spec and Indian-market
sensors. OEM supplier landscape on India-sold cars (not exhaustive, **[UNCONFIRMED]** part-level): **Continental/VDO**
(large TPMS plant in Bangalore; Hyundai-Kia, BMW, FCA, PSA, Mazda, Mitsubishi families share the VDO protocol [16][11]),
**Schrader/Sensata** (GM, Nissan/Infiniti/Renault, BMW, Subaru families [8]), **Pacific Industries** (Toyota [10][12]),
**Huf**, **Bosch**; plus abundant **aftermarket Chinese kits** (Jansite/Steelmate/CUB and generic) [9]. Hyundai Creta
/ Kia Seltos use 433 MHz sensors (OE e.g. 52940-BV100 / 52940-L1100, **[UNCONFIRMED]** protocol = VDO TG1C family)
[13]. Most TPMS are **FSK + Manchester/differential-Manchester, ~50 µs/bit**, 8–20 payload bytes, with a CRC-8 or
CRC-16. **rtl_433 is used only as a fact source/oracle — all facts below are to be re-implemented from them, not
copied; rtl_433 is GPL-2 [1].**

Give the engine these four (each covers a large vehicle population and is fully specified):

#### (a) Toyota / Pacific PMV-C210 — 433.92 MHz, FSK, differential Manchester, CRC-8 [1 (tpms_toyota.c)]
- Mod: FSK PCM, ~52 samples/bit ("12–13 @250k"), differential Manchester, reset gap 150 µs.
- Preamble `0x55 0x3c` (`0101 0101 0011 11`), search `{0xa9,0xe0}` (12 bits); 14-bit sync + 72 Manchester data bits.
- Frame (9 bytes after diff-Manchester): `IIII` 32-bit ID · b4 bit7 status · **pressure1** = `(b4&0x7f)<<1 | b5>>7` ·
  **temp** = `(b5&0x7f)<<1 | b6>>7` · status = `(b4&0x80)|(b6&0x7f)` · b7 = inverted pressure (`b7^0xff`) · b8 = CRC-8.
- Formulas: **PSI = pressure1 × 0.25 − 7.0** (28 raw = 0 PSI); **°C = temp − 40**.
- Integrity: `crc8(b,8, poly=0x07, init=0x80)` == b8 **AND** pressure1 == (b7^0xff).

#### (b) Continental / VDO (Hyundai-Kia/BMW/FCA/PSA/Mazda…) — 433.92 MHz, FSK, Manchester, CRC-8 [1 (tpms_hyundai_vdo.c)]
- Mod: FSK (±38 kHz), 52 µs/bit (FCC test uses 42 µs), Manchester → 80 bits / 10 bytes, reset 150 µs.
- Preamble `55 55 55 56`; decoder inverts and searches `aa aa aa a9` (32 bits).
- Frame `UU IIIIIIII FR PP TT BB CC`: b0 state (0x20–0x23 seen) · b1–b4 32-bit ID · b5 hi nibble flags / lo nibble
  repeat · **b6 pressure** · **b7 temp** · b8 maybe-battery · b9 CRC-8.
- Formulas: **kPa = P × 1.375** (PSI = P/5); **°C = T − 50**.
- Integrity: `crc8(b,9, poly=0x07, init=0xaa)` over **bytes 0–8** (includes UU) == b9. ~10–11 packets/burst, ~every 64 s.

#### (c) Schrader (Nissan/Infiniti/Renault/Subaru SMD3MA4 family) — OOK, Manchester, parity [1 (schraeder.c)]
- Mod: OOK PCM/Manchester-zerobit, short 120 µs, reset 480 µs. (Generic Schrader variant is ASK, 68 bits.)
- SMD3MA4/NIS315G3: preamble `0xF5555555E` (search last 16 bits `55 5e`), Manchester `01→0,10→1`, 38 data bits:
  `1 fixed + 3 flags + 24 ID + 8 pressure + 2 checksum`. Pressure: SMD3MA4 **0.2 PSI/count**, NIS315G3 **0.25 PSI/count**;
  **no temperature**. Integrity: 2-bit-group addition checksum ≡ 1 (labelled PARITY). Flags: 0=learn,3=pressure-change,
  5=wake,7=driving. Example `{37}698e08eb48`.
- Generic Schrader (ASK, 68 bits, `S PF FI II II II PP TT CC`): 28-bit ID, pressure `b5×25` mbar (×0.1→kPa), temp
  `b6−50 °C`, **CRC-8 poly 0x07 init 0xF0** over 7 bytes == b7. Example `7 f6 70 3a 38 b2 00 49 49`.

#### (d) Aftermarket Chinese — Jansite TY02S (common in India kits) — 433.92 MHz, FSK, Manchester [1 (tpms_jansite.c)]
- Mod: FSK PCM, 52 µs/bit, Manchester → 56 bits / 7 bytes, reset 150 µs. Preamble (after invert) `aa aa a9` (24 bits).
- Frame `II II II IS PP TT CC`: 28-bit ID · 4-bit status · **P** pressure · **T** temp · CC checksum (algorithm not
  published — rtl_433 leaves it as TODO). Formulas (best-guess in rtl_433): **kPa ≈ P × 1.7**, **°C = T − 50**.
- Because the checksum is unverified, label **"-like"** and lean on repetition (ID stable across the 4-frame burst).

**[Oracle]** rtl_433 also has TRW (tpms_trw.c), Citroen/PSA (tpms_citroen.c: XOR checksum of b1..b9==0, kPa=b6×1.364,
°C=b7−50), Ford, Renault, GM, Kia (CRC-8 poly 0x07 init 0x76), Pacific PMV-107J @315 MHz (diff-Manchester, CRC-8 poly
0x13 init 0) — use as further oracles if a capture needs them.

**TPMS detector → … → label (shared):** FSK/OOK burst in 433.05–434.79 (or 315) MHz → instantaneous-frequency slice at
~52 µs/bit → try each protocol's preamble search → Manchester/diff-Manchester decode → CRC/checksum → on pass, emit
`<protocol> TPMS` with `id, pressure_kPa, temperature_C, flags, battery`; on preamble-match-but-CRC-fail, "-like".
**Validity gates (so noise is not labelled):** exact bit count, preamble match, CRC/checksum pass, plausibility
(pressure 0–450 kPa, temp −40..85 °C — Schrader-BMW uses exactly these bounds [8]). Require CRC pass OR ≥2 identical
frames for a confirmed label.

**Synthetic tests:** generate each frame from spec fields, Manchester/diff-Manchester encode, FSK/OOK-modulate at the
sensor's offset onto a 433.92/315 MHz carrier at 2 MSps, 4–6 repeats, AWGN. Assert CRC passes and fields round-trip.

### 1.5 Weather sensors / doorbells (well-documented, optional) — IMPLEMENTABLE if time permits
rtl_433 covers hundreds (Acurite, LaCrosse, Oregon Scientific, Nexus/temperature-humidity on 433.92 MHz; Ring/generic
doorbells reuse EV1527/PT2262) — all GPL-2 oracles [1]. For S6 the EV1527/PT2262 slicer already catches most Indian
doorbells. Treat named weather protocols as a later add; for now a generic 433 OOK sensor that does not match
EV1527/PT2262/KeeLoq/TPMS → `433 MHz OOK device (unknown)` with its measured timing as evidence.

---

## 2. Bluetooth BR/EDR (Classic) — access code / LAP recovery — IMPLEMENTABLE (verified below)

**PHY:** 79 RF channels, `f = 2402 + k` MHz, k=0..78; 1600 hops/s; 625 µs slots; **GFSK, 1 Mb/s, BT=0.5, modulation
index h = 0.28–0.35**. EDR packets keep the GFSK access code + header, then switch the payload to **π/4-DQPSK (2 Mb/s,
EDR2)** or **8DPSK (3 Mb/s, EDR3)**; a guard+sync precedes the DPSK section [17].

**Packet start = Access Code:** `[4-bit preamble] [64-bit sync word] [4-bit trailer]` = 72 bits (or 68 if no trailer).
Preamble `1010`/`0101` (chosen so the first sync bit continues the alternation); trailer likewise follows the last
sync bit [17].

### 2.1 Sync-word construction (Core Spec Vol 2 Part B §6.3.3) — verified against libbtbb [18]
The 64-bit sync word is a **(64,30) expurgated BCH code** over the 24-bit LAP (lower address part of BD_ADDR),
whitened by a fixed PN sequence and a Barker-code tail.

Let the LAP be `A = a23..a0`. Build 30 information bits `X = x0..x29`:
- `x0..x23 = LAP` (a0..a23),
- `x24..x29 = Barker-derived 6-bit tail`: **`001101`** if `a23 = 0`, else its complement **`110010`** (this is the
  6-bit barker extension that makes the trailer alternate) [17][18].

Then:
1. **PN:** `PN = 0x83848D96BBCC54FC` (64-bit). Its upper 34 bits `p34..p63` whiten the information part; the lower 30
   bits whiten the parity. (Equivalently, XOR the whole 64-bit codeword with PN at the end — see code below.) [18]
2. **Generator** (degree-34) `g(D)` = octal **`0260534236651`** = `0x585713DA9` [17][18].
3. Compute the (64,30) codeword: info bits in the high 30 positions, 34 parity bits = `(X·2^34) mod g`.
4. **Sync word = codeword XOR PN.**

**Verified Python (matches libbtbb's `DEFAULT_CODEWORD` for LAP 0 and its `btbb_gen_syncword` exactly):**
```python
PN = 0x83848D96BBCC54FC
G  = int('260534236651', 8)              # 0x585713DA9, degree 34
def polymod(a, g):
    dg = g.bit_length() - 1
    while a.bit_length() - 1 >= dg:
        a ^= g << (a.bit_length() - 1 - dg)
    return a
def syncword(lap):                       # returns 64-bit int; bit k = s_k, s_0 sent first
    barker = 0b101100 if ((lap >> 23) & 1) == 0 else 0b010011   # '001101'/'110010' as bits 0..5
    x  = lap | (barker << 24)            # 30 info bits x0..x29
    xt = x ^ (PN >> 34)                  # whiten info part
    c  = polymod(xt << 34, G)            # 34 parity bits
    return ((xt << 34) | c) ^ PN
# syncword(0x000000) == 0xB0000002C7820E7E  (== libbtbb DEFAULT_CODEWORD)  ✓ verified
# syncword(0x9E8B33) == 0x4E7A2CCE331A3AE2  (GIAC)                          ✓
```
The code's **minimum Hamming distance is 14** (sampled over 2×10^5 random LAP pairs — verified), so up to ~6 bit errors
are correctable in principle; practical detectors allow far fewer.

### 2.2 Detecting a sync word & extracting the LAP (what Ubertooth/libbtbb do — math, not their code) [18]
Slide a 64-bit window over the GFSK-sliced bit stream. For each window `w`:
1. **Barker/trailer check** on the top bits (`MAX_BARKER_ERRORS = 1` in libbtbb).
2. Candidate `lap = (w >> 34) & 0xFFFFFF` (the LAP sits in the information bits).
3. Recompute `syncword(lap)` and count `errs = popcount(w ^ syncword(lap))`. Accept if `errs ≤ max_ac_errors`
   (libbtbb `AC_ERROR_LIMIT = 5`, typically 2 in sniffers). Zero syndrome ⇒ exact codeword.
This both **validates** that a 64-bit window is a real access code and **gives the 24-bit LAP** in one step. For the
general-inquiry access code (GIAC) use the fixed reserved **LAP 0x9E8B33** [17].

### 2.3 Header, HEC and UAP recovery [17][18]
After the access code comes an **18-bit header**, sent with **rate-1/3 repetition** (each bit 3×, majority vote → 18
bits): `LT_ADDR(3) · TYPE(4) · FLOW(1) · ARQN(1) · SEQN(1) · HEC(8)`.
- **Whitening:** header (and payload) are XORed with a data-whitening sequence seeded by **CLK1 (the 6 LSBs of the
  master clock)**. libbtbb uses a `WHITENING_DATA` table indexed by `INDICES[clock & 0x3f]` (+skip for the payload) [18].
  The clock is unknown to a passive sniffer, so you try candidate CLK values.
- **HEC:** 8-bit CRC over the first 10 header bits, **initialised with the UAP** (upper address part, 8 bits) [17].
  `uap_from_hec` reverses the HEC LFSR to back out the UAP: process the data bits, `if (hec & 0x80) hec ^= 0x65;
  hec = (hec<<1) | (((hec>>7) ^ bit) & 1)`, final `reverse(hec)` → UAP [18].
- **UAP recovery:** for each CLK1 candidate (0–63), de-whiten the header, de-repeat, run `uap_from_hec`, and keep CLK
  values whose resulting UAP is self-consistent across **multiple packets of the same piconet** (same LAP). A single
  packet is **ambiguous** (many CLK/UAP pairs pass); UAP converges only after several packets and/or using the payload
  CRC [18]. The payload HDLC/CRC-16 (`reg = (reverse(UAP)<<8)&0xff00` init [18]) further confirms UAP+CLK.

### 2.4 What one 34 ms / 61.44 MSps block gives; EDR; vendor limitation
- **One block:** 56 MHz usable spans ~56 of the 79 BT channels at once. In 34 ms (≈54 slots) you may catch **several
  packets on different hop channels**. For each clean GFSK burst: detect the access code → **LAP** (+ errors), header
  type, and attempt **UAP**. Full UAP/CLKN usually needs packets accumulated across blocks.
- **EDR detection:** after the GFSK header, a constant-amplitude **phase-modulated** section (π/4-DQPSK or 8DPSK) with
  a different symbol constellation ⇒ flag `edr` (2 or 3 Mb/s) by detecting the modulation switch; the access code/header
  decode is identical to BR [17].
- **BD_ADDR → vendor limitation:** only the **LAP (24 bit)** and, with work, **UAP (8 bit)** are on air. The **NAP (16
  bit, upper)** is **never transmitted**, so the OUI (NAP+UAP top) cannot be fully recovered ⇒ **vendor cannot be
  resolved from BR/EDR sniffing** (unlike BLE, where the full public address, if used, gives the OUI). Report
  `lap`, `uap` (if resolved, else "ambiguous"), `giac` flag, `edr` flag; **do not claim a vendor**.

**Synthetic tests:** build a packet = preamble + `syncword(lap)` + trailer + rate-1/3 header (with a correct
HEC seeded by a chosen UAP, whitened by a chosen CLK1), GFSK-modulate (h≈0.32, BT=0.5) at 1 Mb/s on channel k, AWGN.
Assert: sync detector finds the window, recovers `lap`, and UAP brute-force returns the chosen UAP for the right CLK1.
Vectors: `LAP 0x9E8B33 → 0x4E7A2CCE331A3AE2`, `LAP 0x000000 → 0xB0000002C7820E7E`.

---

## 3. BLE trackers & notable devices (from advertising) — IMPLEMENTABLE

The library already CRC-decodes BLE adverts (access address, de-whiten, CRC-24, PDU type, AdvA, AD structures incl.
local name, company ID, 16-bit UUIDs, TX power) — see `receivers/rc_link/plugin.py`, `ble_coded.py`,
`gfsk_coherent.py`, `core/remote_id.py` [ADDR: schema row "Bluetooth LE advertising"]. This section adds **content
recognizers** on top of decoded AD structures. Each has strict validity checks so noise is never "-like"-labelled.

### 3.0 AD-structure primitives
A decoded advert is a list of AD structures `[len][type][data…]`. The recognizers key on:
- **Manufacturer Specific Data** (AD type `0xFF`): first 2 bytes = company ID (LE). Apple = **0x004C**, Microsoft =
  **0x0006**, Samsung = **0x0075**.
- **Service Data – 16-bit UUID** (AD type `0x16`): first 2 bytes = UUID (LE).
- **16-bit Service UUIDs list** (AD types `0x02/0x03`).

### 3.1 Apple Find My / AirTag offline-finding (manufacturer 0x004C, type 0x12) — IMPLEMENTABLE [19][20][21]
Inside Apple's `0xFF` manufacturer data, a TLV of **type `0x12` (Find My)**:

**Separated / "offline finding" variant (length byte `0x19` = 25):**
| Offset (within 0x12 payload) | Field | Notes |
|---|---|---|
| 0 | type = `0x12` | Find My |
| 1 | length = `0x19` (25) | 25 bytes follow |
| 2 | **status** | bit2 "maintained" = owner connected within last 15 min; bits6–7 = battery (full/med/low/crit) when maintained [20] |
| 3–24 | public key bytes 6–27 | 22 bytes of the EC P-224 public-key x-coordinate [20] |
| 25 | key bits | bits0–1 = bits 6–7 of key byte 0 (rest reserved) [20] |
| 26 | hint | byte 5 of the primary key's BT_ADDR (usually 0x00) [20] |

Key bytes 0–5 are carried in the **BLE advertising address** itself (the top 6 bytes of the key); the address is a
random/static type with the two MSBs forced to `11`. **Near-owner / "nearby" variant:** AirGuard filters Apple mfg
data `12 19 10` (so **byte after status often `0x10`** for AirTag) vs AirPods `12 19 18`; connection state is read as
**OVERMATURE_OFFLINE when mfg[1]==0x19** [19]. A short `0x12 0x02 …` form (length 2) is the "near owner/registered"
status-only beacon [UNCONFIRMED exact bytes].

**Device-type nibble** (AirGuard, mfg byte 2 status): `deviceType = (status & 0x30) >> 4` → **0 = AppleDevice/Find My
generic, 1 = AirTag, 3 = AirPods**; battery `= (status >> 6) & 3` [19].

**Validity:** company==0x004C, TLV type==0x12, length ∈ {0x19, 0x02}, total AD length consistent, advertising address
type == random-static (top bits 11). On all passing → `Apple Find My / AirTag-like`, report status/battery/deviceType;
**no serial/owner** is recoverable (key is rotating & anonymous).

### 3.2 Other Apple Continuity types (manufacturer 0x004C) — IMPLEMENTABLE (fingerprint-level) [21][22]
Same `0xFF`/0x004C container, one or more TLVs `[type][len][data]`. Types (furiousMAC/PETS, repo GPL-2 → **facts/oracle
only** [21]): `0x02` iBeacon, `0x05` AirDrop, `0x07` Proximity Pairing (AirPods), `0x09` AirPlay Source, `0x0A` AirPlay
Target [UNCONFIRMED], `0x0B` Magic Switch, `0x0C` Handoff, `0x0D` Tethering Target, `0x0E` Tethering Source, `0x0F`
Nearby Action, `0x10` Nearby Info, `0x12` Find My.
- **0x07 Proximity Pairing (AirPods):** `07 LL 01 [model:2] [status] [batt] [charge+case] [lid-counter] [color] 00
  [16B encrypted]`. Model is `uint16` (e.g. **AirPods Pro = 0x0E20** [21]); battery nibbles R/L; color 0x00=white.
  **Model-ID table is only partially published** → report `model_id` hex + a small known-model lookup, else "-like".
- **0x10 Nearby Info / 0x0C Handoff:** activity/state flags; good for "Apple device present" but **do not reveal
  identity** (ephemeral). Use for a device-class hint only.

**Validity:** company==0x004C, each TLV's length consistent to the end of the field; `0x02` iBeacon must be exactly
`02 15 [16B UUID][major:2][minor:2][txpower:1]` (len byte 0x15) — this is the strongest Apple check.

### 3.3 Samsung SmartTag offline-finding (Service Data UUID 0xFD5A) — IMPLEMENTABLE (verified layout) [23][19]
Registered SmartTag broadcasts **service data under 16-bit UUID 0xFD5A** (non-registered: 0xFD59). 20-byte OF payload
(Yu et al. USENIX'23, Fig. 3 [23]):

| Byte(s) | Field | Detail |
|---|---|---|
| 0 | info | bits0–3 = version, bit4 = adv type, **bits5–7 = tag state** |
| 1–3 | **aging counter** (24-bit) | `agingCounter = (tagTime − 1593648000)/900`; updates every 15 min |
| 4–11 | **privacy ID** (8 bytes) | rotating anonymous ID (pool size 1000) |
| 12 | flags | bits0–3 region ID, bit4 E2E-encryption, **bit5 UWB**, bits6–7 battery level |
| 13–15 | reserved | |
| 16–19 | **signature** (4 bytes) | truncated AES over first 16 bytes |

**Tag state (bits5–7 of byte0)** [23]: `1=premature-lost, 2=lost, 3=overmature-lost, 4=paired-one, 5=connected-one,
6=connected-two`. AirGuard reads the bit order swapped: `!b5&b6&b7→overmature, !b5&b6&!b7→offline, !b5&!b6&b7→premature`
[19]; battery from byte12 bits6–7 `11=full,10=med,01=low,00=very-low` [19]. UWB present ⇒ SmartTag+ vs SmartTag.
**Example (Yu et al.):** `15 6b fa 00 c8 40 62 b2 8f 00 e2 60 c3 00 00 00 ad 01 8b 47` (byte0=0x15 → state 0 in their
capture; aging `6bfa00`; privacy ID `c84062b28f00e260`; flags `c3`; sig `ad018b47`).
**Validity:** service UUID==0xFD5A, payload length ≥ 20, tag-state ∈ {1..6}; aging counter monotone across a session.
`Samsung SmartTag-like`, report state/battery/UWB; privacy ID is the rotating token (not a stable identity).

### 3.4 Tile (Service Data UUID 0xFEED; also 0xFEEC/0xFD84) — IMPLEMENTABLE (fingerprint) [19][24]
AirGuard filter: service UUID **0xFEED**, service data starting `02 00` (mask `FF FF`) [19]. (Tile also uses **0xFEEC**
for its "activate/registration" service and **0xFD84** registered to Tile with the SIG [24] **[UNCONFIRMED which is on
which firmware]**.) Tile IDs are not in the clear here; treat as presence + `Tile-like`. **Validity:** UUID∈{0xFEED,
0xFEEC,0xFD84}, service-data length consistent.

### 3.5 Google Find My Device Network / FMDN (Service Data UUID 0xFEAA, frame 0x40/0x41) — IMPLEMENTABLE [25][19]
(Note: 0xFEAA is shared with **Eddystone** — disambiguate by the frame-type byte.) Service data under 0xFEAA:
| Octet (service data) | Value | Meaning |
|---|---|---|
| 0 | `0x40` or `0x41` | **frame type** (0x41 = Unwanted-Tracking-Protection on) |
| 1..20 | 20 bytes (160-bit curve) | **Ephemeral Identifier** (EIK-derived) |
| 21 (opt.) | hashed flags | bits5–6 battery (00 unsupported,01 normal,10 low,11 critical), bit7 UTP; XOR'd with LSB of SHA256(r) |
(256-bit-curve variant: 32-byte EID, flags at octet 33.) AirGuard subtype from an advertising-flags value: **0x02 =
phone, 0x06 = tag** [19]. Advertise ≥ every 2 s; address rotates on average every 1024 s (near-owner), fixed in UTP.
**Validity:** service UUID==0xFEAA, **frame byte ∈ {0x40,0x41}**, service-data length ∈ {24,25} (160-bit) or {36,37}
(256-bit) → this length+frame test is what separates FMDN from Eddystone. `Google Find My (FMDN)-like`.

### 3.6 Eddystone (Service Data UUID 0xFEAA, frame 0x00/0x10/0x20/0x30) — IMPLEMENTABLE [26]
Same 0xFEAA, different frame byte: `0x00` UID (10-byte namespace + 6-byte instance + TX power), `0x10` URL
(TX power + scheme byte + encoded URL), `0x20` TLM (telemetry: battery mV, temp, adv count, uptime), `0x30` EID.
**Validity:** frame byte in that set and the per-frame fixed lengths (UID = 0x00 + txpwr + 10 + 6 + 2 RFU).
Disambiguation from FMDN: frame byte (FMDN uses 0x40/0x41). `Eddystone-<UID/URL/TLM/EID>`.

### 3.7 iBeacon (manufacturer 0x004C, type 0x02 len 0x15) — IMPLEMENTABLE [27]
`4C00 02 15 [UUID:16] [major:2] [minor:2] [tx_power:1]`. **Validity:** company 0x004C, bytes `02 15`, total 25 bytes.
Report proximity UUID / major / minor / measured power. Strongest, zero-false-positive check in the Apple family.

### 3.8 Microsoft Swift Pair (manufacturer 0x0006) — IMPLEMENTABLE [28]
Inside `0xFF`/0x0006 vendor section: `[beaconID] [subScenario] [reserved RSSI=0x80] [optional display name…]`.
**Sub-scenario:** `0x00` = LE-only, `0x01` = BR/EDR-only (then a **little-endian BD_ADDR** in the main packet),
`0x02` = LE+BR/EDR Secure Connections. Fast advertising (~30 ms) while pairing. **Validity:** company==0x0006,
reserved byte==0x80, sub-scenario ∈ {0,1,2}. Report `swift_pair` + display name (ASCII). `Microsoft Swift Pair-like`.

### 3.9 Google Fast Pair (Service Data UUID 0xFE2C) — IMPLEMENTABLE [29]
- **Discoverable (pairing):** service data = **24-bit model ID** (3 bytes) → look up the (public, partial) model-ID
  registry; report `model_id` hex. Interval ≤ 100 ms.
- **Not-discoverable:** `[version/flags=0x00] [len|type nibble] [account-key Bloom filter] [0x11 salt type][salt:2]`.
  Bloom filter size `s = floor(1.2·n)+3`; bits set from `SHA256(key‖salt)` split into eight 32-bit ints mod (s·8).
  Interval ≤ 250 ms. **Validity:** UUID==0xFE2C; discoverable form length==3; non-discoverable: flags byte 0x00 and a
  self-consistent TLV chain ending in the 2-byte salt. `Google Fast Pair-like` (+ model ID when discoverable).

### 3.10 Chipolo / Pebblebee / others — IMPLEMENTABLE (fingerprint) [19]
- **Chipolo:** service UUID **0xFE33**; firmware ≤ D8 exposes an offline bit in service-data byte1 bit0 [19].
- **Pebblebee:** its own offline-finding service UUID (AirGuard `PebbleBee.offlineFindingServiceUUID`) [19] [UNCONFIRMED
  exact 16-bit value — read from a capture].
- **Samsung Find-My-Mobile (phones/other):** service UUID **0xFD69** [19].

### 3.11 IETF DULT (Detecting Unwanted Location Trackers), accessory protocol — IMPLEMENTABLE (cross-vendor) [30]
`draft-ietf-dult-accessory-protocol-00` (2024-11-04, expired 2025-05; protocol v1.0.0) defines a **common** advertising
payload so any platform can detect trackers:
| Field | Bytes | Notes |
|---|---|---|
| MAC address | 0–5 | |
| Flags TLV | 6–8 | optional |
| Service Data TLV | 9–12 | type 0x16, **UUID 0xFCB2** |
| **Network ID** | 13 | per-manufacturer registered value (registry removed from -00; **values UNCONFIRMED**) |
| **Near-owner bit** | 14 | **bit0 = 1 ⇒ near-owner mode; 0 ⇒ separated** (DULT "near owner"); bits1–7 reserved |
| company payload | 15–36 | variable |
Address rotation: 15 min near-owner, 24 h separated. Max adv interval 4 s (≤2 s recommended). Non-owner control GATT
service `15190001-12F4-C226-88ED-2AC5579F2A85` (shared with Google FMDN sound). **Validity:** service UUID==0xFCB2,
byte14 bits1–7 == 0. `DULT tracker-like (near-owner / separated)` — this is the most future-proof recognizer; prefer
it as the generic fallback for any 0xFCB2 advert whose vendor is unknown.

---

## 4. Integration notes for S6/S7

- **Receiver placement:** sub-GHz OOK/FSK belongs in a new `rc_link`-adjacent receiver (or a `short_range` receiver)
  that imports only numpy/scipy + `spectrum_intel.core`; BLE-content recognizers extend the existing `rc_link` BLE path
  (operate on already-decoded AD structures from `core/remote_id.ble_ad` / `plugin.py`). BR/EDR LAP is a new 2.4 GHz
  GFSK path alongside the existing BLE GFSK detector.
- **Schema:** add `technology`/`family` codes (e.g. `tpms`, `ook_remote`, `bt_classic`, `ble_tracker`) in
  `schema/technologies.json`, `descriptions.json`, and `capabilities.json` (tag `research`, `sources:
  research/signal_gaps/04_short_range.md`). Response additions (ids, pressure_kPa, temperature_C, lap, uap, model_id,
  tracker_vendor) are **additive** (minor schema bump, Python + C++ together).
- **Confidence / labelling:** CRC/checksum pass OR ≥2 identical repeats ⇒ confirmed label; else "-like" ≤0.5 ⇒
  `unknown` with evidence (measured timing, matched preamble, partial fields). Never claim a BR/EDR vendor (no NAP).
- **Latency:** the pulse slicer runs once per burst off the envelope/inst-freq already needed for detection; BR/EDR LAP
  check is a 64-bit sliding popcount — both cheap. Gate TPMS protocol trials behind a preamble match so only one
  Manchester decode runs per burst.
- **C++ GPL cleanliness:** all §1 facts came from datasheets [5][7][14][15] and GPL rtl_433 **as an oracle** [1] — the
  CRC polynomials/inits, bit layouts and formulas are protocol facts, re-implemented from scratch. §2 came from the
  Core Spec [17] and were **verified** against GPL libbtbb [18] (oracle only). §3 came from vendor/IETF/CC-BY docs
  [25][28][29][30], the GPL-2 furiousMAC repo [21] and Apache-2.0 AirGuard [19] (facts/oracle). Nothing GPL is copied
  into C++.

---

## Sources

1. rtl_433 device decoders (schraeder.c, tpms_toyota.c, tpms_hyundai_vdo.c, tpms_jansite.c, tpms_citroen.c, tpms_kia.c,
   tpms_pmv107j.c, tpms_elantra2012.c, generic_remote.c, hcs200.c) — https://github.com/merbanan/rtl_433 — **GPL-2.0-or-later. Fact source / oracle only; not copied.**
2. EV1527 reference (adom.inc wiki) — https://wiki.adom.inc/ray/ev1527 — docs, licence unstated.
3. PT2262 / HX2262 datasheet (Princeton / Holtek) — https://www.mikrocontroller.net/attachment/599366/HX2262.pdf — vendor datasheet.
4. PT2262 datasheet copy — https://static.chipdip.ru/lib/159/DOC013159025.pdf — vendor datasheet.
5. Microchip HCS301 datasheet (DS21143) — https://www.alldatasheet.net/datasheet-pdf/pdf/74925/MICROCHIP/HCS301.html — vendor datasheet.
6. rtl_433 hcs200.c (HCS200/300 timing & 66-bit layout) — https://github.com/merbanan/rtl_433/blob/master/src/devices/hcs200.c — **GPL-2. Oracle only.** Cites Microchip HCS200 (DS40138C) & HCS300 (DS21137G) datasheets.
7. Microchip HCS200 datasheet DS40138 — http://ww1.microchip.com/downloads/en/devicedoc/40138c.pdf — vendor datasheet.
8. rtl_433 schraeder.c (Schrader generic/EG53MA4/SMD3MA4/NIS315G3/MRXBC5A4 BMW) — https://github.com/merbanan/rtl_433/blob/master/src/devices/schraeder.c — **GPL-2. Oracle only.**
9. rtl_433 tpms_jansite.c (aftermarket Jansite TY02S) — https://github.com/merbanan/rtl_433/blob/master/src/devices/tpms_jansite.c — **GPL-2. Oracle only.**
10. rtl_433 tpms_toyota.c (Pacific PMV-C210) — https://github.com/merbanan/rtl_433/blob/master/src/devices/tpms_toyota.c — **GPL-2. Oracle only.**
11. rtl_433 tpms_hyundai_vdo.c (Continental/VDO TG1C family) — https://github.com/merbanan/rtl_433/blob/master/src/devices/tpms_hyundai_vdo.c — **GPL-2. Oracle only.**
12. rtl_433 tpms_pmv107j.c (Pacific PMV-107J, 315 MHz, Toyota) — https://github.com/merbanan/rtl_433/blob/master/src/devices/tpms_pmv107j.c — **GPL-2. Oracle only.**
13. TPMS OE part-number listings (Creta/Seltos/Mahindra) — ebay/kiapartsnow third-party listings — **[UNCONFIRMED]** part/protocol mapping.
14. EV1527 timing (done.land) — https://done.land/components/data/datatransmission/wireless/shortrangedevice/am/ask/ev1527/ — community docs.
15. PT2260 datasheet — https://www.mikrocontroller.net/attachment/33295/PT2260.pdf — vendor datasheet.
16. "Continental expands TPMS production capacity in India" — S&P Global AutoTechInsight — trade news (supplier presence).
17. Bluetooth Core Specification, Vol 2 (BR/EDR Controller), Part B §6.3 (access code, sync word, PN, generator, Barker), header/HEC/whitening — https://www.bluetooth.com/specifications/specs/core-specification/ — Bluetooth SIG (normative spec).
18. libbtbb (bluetooth_packet.c) — https://github.com/greatscottgadgets/libbtbb — **GPL-2. Oracle only (used to VERIFY the §2 math; not copied).** Confirms PN=0x83848D96BBCC54FC, generator octal 0260534236651, DEFAULT_CODEWORD=0xB0000002C7820E7E, uap_from_hec (0x65), error limits.
19. AirGuard (seemoo-lab) device-type detectors (AirTag.kt, AppleFindMy.kt, AirPods.kt, SamsungTracker.kt, GoogleFindMyNetwork.kt, Tile.kt, Chipolo.kt, PebbleBee.kt, DeviceManager.kt) — https://github.com/seemoo-lab/AirGuard — **Apache-2.0 (safe to use; facts/oracle).**
20. furiousMAC continuity — Find My message (findmy.md) — https://github.com/furiousMAC/continuity/blob/master/messages/findmy.md — **GPL-2. Fact source only.**
21. furiousMAC continuity — message-type list & proximity_pairing.md (AirPods 0x07, model IDs) — https://github.com/furiousMAC/continuity — **GPL-2. Fact source only.** Based on Martin et al. (PETS 2019) & Celosia/Cunche (PETS 2020).
22. Martin et al., "Handoff All Your Privacy", PoPETs 2019; Celosia & Cunche, "Discontinued Privacy", PoPETs 2020 — academic papers (Apple Continuity).
23. Yu, Henderson, Tiu, Haines, "Privacy Analysis of Samsung's Crowd-Sourced Bluetooth Location Tracking System", USENIX Security 2023 — https://arxiv.org/abs/2210.14702 — academic paper (FD5A layout, Fig. 3).
24. Bluetooth SIG 16-bit UUID assignments (0xFEED/0xFEEC Tile, 0xFD84 Tile, 0xFD5A/0xFD59 Samsung, 0xFD69 Samsung, 0xFE33 Chipolo, 0xFE2C Google Fast Pair, 0xFEAA Google/Eddystone, 0xFCB2 DULT) — https://www.bluetooth.com/specifications/assigned-numbers/ — Bluetooth SIG.
25. Google Find My Device / FHN (FMDN) specification — https://developers.google.com/nearby/fast-pair/specifications/extensions/fmdn — **CC BY 4.0 (code Apache-2.0).**
26. Eddystone protocol specification — https://github.com/google/eddystone/blob/master/protocol-specification.md — **Apache-2.0.**
27. Apple iBeacon (Getting Started with iBeacon / Proximity Beacon spec) — https://developer.apple.com/ibeacon/ — Apple docs.
28. Microsoft Swift Pair (component guidelines) — https://learn.microsoft.com/en-us/windows-hardware/design/component-guidelines/bluetooth-swift-pair — Microsoft docs.
29. Google Fast Pair — Provider advertising signal — https://developers.google.com/nearby/fast-pair/specifications/service/provider — **CC BY 4.0 (code Apache-2.0).**
30. IETF draft-ietf-dult-accessory-protocol-00 (DULT) — https://www.ietf.org/archive/id/draft-ietf-dult-accessory-protocol-00.txt — IETF Internet-Draft (IETF Trust / BSD-like).
