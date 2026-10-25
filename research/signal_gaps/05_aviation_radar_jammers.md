# 05: aviation, radars, jammers, AM voice (batch S, items S8–S13)

Research reference for PLAN items S8 (radio altimeter), S9 (ILS glide slope), S10 (named radars and pulse
deinterleaving), S11 (Mode A/C and IFF replies), S12 (noise, spot, cellular and GNSS jammers) and S13 (AM voice).
Written 2026-10-25. No code was changed.

**Conventions**
- **[V]** means verified in this session against a primary or authoritative text (the extracted ITU-R PDF, the
  Mitch et al. paper, or an executed check).
- **[M]** means from memory of the standard and consistent with secondary sources, but the primary text was not
  re-read. Check it before relying on it.
- **[U]** means unconfirmed or inferred. Treat it as an assumption to test.
- "Block" means 2^21 complex samples, 34.13 ms at 61.44 MSps, with at most 56 MHz analog bandwidth (one Pluto
  channel, no direction finding).
- Labels follow the library rules:
  - "-like" unless a parity, CRC or FEC check (or an exact standard structure, for example a Gillham-valid reply
    with framing) verifies the label;
  - below 0.5 confidence the label is `unknown`, with evidence;
  - every item is tagged `research` until a real capture promotes it.

## 0. What already exists (read before implementing)

| Piece | Where | Relevant behaviour / limits |
|---|---|---|
| Pulse groups | `features/measure.py::pulse_train` | 1 µs envelope at ≥1 MSps on the signal's own band. Threshold 10 dB over the median. Gaps < 20 µs join one group (so DME pairs and radar short+long pulses are one group). PRI = the shortest interval cluster (±15 %, ≥ half of intervals, ≥100 µs). LFM fit on the widest pulse ≥10 µs (R² ≥ 0.8). **There is no per-pulse PDW list, no deinterleaving and no frequency per pulse.** |
| Radar / DME labels | `classify/rules.py::r_pulsed` | Runs only in `bands_india.PULSED_BANDS`. Labels: "DME/TACAN pulse pairs-like", "Pulsed radar (S-band ATC)-like", "Pulsed radar-like". Nothing within 1.5 MHz of 1030/1090 is DME. `r_adsb` gives "SSR interrogation (1030 MHz)-like". |
| Swept-carrier jammer | `receivers/jammer/` | STFT with 64 bins (0.96 MHz × 1.04 µs rows). Cyclic line search between max(300 Hz, 8 sweeps per block) and 120 kHz. Fold needs ≥40 phases per period (so periods ≥ ~40 µs). In 4200–4400 MHz it labels "Swept carrier (FMCW-like)". **Gaps:** (a) radio altimeters at 12–1623 Hz give fewer than 8 sweeps per block and fall below the trigger; (b) GNSS PPDs at 8–27 µs periods give fewer than 40 phases per period at 1.04 µs rows. Both need the dedicated detectors below. |
| Mode S | `receivers/adsb/` | CRC-verified DF frames. DF4 (AC13 altitude) and DF5 (ID13 squawk) cross-check Mode A/C. Verified Mode S frames must mask the time span they occupy before Mode A/C bracket search (S11). |
| FM voice | `receivers/voice_fm/` | Narrow-peak candidates in 1 kHz bins (≤40 kHz wide). Envelope constancy cv ≤ 0.35. Rejects an envelope that varies "more than the noise explains" (cv > max(0.12, 1.6/√(2·SNR))). **That rejected class is the AM voice candidate class**: S13 should take candidates that voice_fm rejects for envelope variation, not duplicate the candidate search. |
| Bands | `core/bands_india.py` | `radalt` 4200–4400, `radar_s_atc` 2700–2900, `radar_s_marine` 2900–3100, `radar_s_def` 3100–3300, `radar_5400`, `unii_2e` (weather 5600–5650), `arns_l`, `rnss_l2_radar`, `arns_radar`, `fixed_radar_1350`. **Missing:** `ils_gs` 328.6–335.4 MHz and `uhf_mil_air` 225–400 MHz. Add them as rows; the current table starts at 406 MHz for the government bands. |

---

## 1. Radio altimeters, 4200–4400 MHz (S8)

### 1.1 Characteristics: ITU-R M.2059-0 (02/2014), Tables 1 and 2 [V]

All types are centred nominally on 4300 MHz.

| Type | Modulation | Chirp BW (MHz) | Waveform repetition | Peak power | Max drift | Offset between units | 3 / 20 / 40 dB emission BW (MHz) | Antenna |
|---|---|---|---|---|---|---|---|---|
| A1 analogue | FMCW | 104 | 49–51 Hz | 0.6 W | ±15 MHz | 5 MHz | 110 / 120 / 180 | 10 dBi, 40–60° |
| A2 | FMCW | 132.8 | 150 Hz | 1 W | ±15 MHz | 5 MHz | 162.8 / 170 / 180 | 10 dBi, 55° |
| A3 | FMCW | 133 | **12–1623 Hz (variable with altitude)** | 0.1–0.25 W | ±20 MHz | 0 | 171 / 181 / 191 | 10 dBi |
| A4 | Pulsed | – | 10 000 pps, PW 130 ns | 100 W | – | – | 8 / 44 / 130 | 13 dBi |
| A5 | Pulsed | – | 20 000 pps, PW 200 ns | 5 W | – | – | 7 / 29 / 108 | 11 dBi |
| A6 | Pulsed | – | 6 000 pps, PW 75 ns | 40 W | – | – | 15 / 51 / 131 | 11 dBi |
| D1 digital | FMCW | 150 | 143 Hz fixed | 0.4 W | ±0.22 MHz | span offset by SDI −2.5/0/+2.5 MHz; timing adjusted on interference | 150 / 153 / 180 | 11 dBi |
| D2 | FMCW | 176.8 | 1000 Hz fixed | 0.1 W | ±0.129 MHz | −5/0/+5 MHz by system number, **plus a per-system linear frequency-hop pattern** | 177 / 180 / 190 | 10 dBi |
| D3 | FMCW | 133 | 100–4700 Hz | 0.1–1 W | ±0.22 MHz | – | 175 / 185 / 196 | 8–11 dBi |
| D4 | Pulsed | – | 25 000 pps, PW 30 or 225 ns | 5 W | – | – | 5 or 31 / 26 or 105 / 106 or 195 | 13 dBi |

Further facts from the same text:
- **Waveform.** All FMCW units use a **triangle**. M.2059 §2.1 and Annex 2 §1.1 say "fixed or variable duration of the
  triangle FMCW waveform". Variable-period units (A3, D3) hold the beat frequency constant: the period grows with
  altitude.
- **Repetition frequency.** The table's "waveform repetition frequency" is read here as the full triangle frequency
  (up + down), so one ramp lasts 1/(2·f_rep) [U: the table does not define it; Table 2 D2's hop pattern and D3's
  range make a sawtooth possible for some digital units].
- **Analogue units:** open-loop VCO, typically up to ±25 MHz over temperature.
- **Installation:** up to 3 units per aircraft, all transmitting at once.
- **Pulsed altimeters** use "spread spectrum techniques" and may occupy the whole 200 MHz.
- **Sweep limits:** f1 = 4300 − Bs/2 − drift and f2 = 4300 + Bs/2 + drift (M.2059 footnote 6).

### 1.2 What a 56 MHz receiver sees

Assuming a triangle, the sweep rate is r = 2·B·f_rep. For a centred capture window W = 56 MHz (in practice the
usable ~52 MHz after the anti-alias roll-off), one crossing lasts t_x = W/r.

| Type | r (MHz/ms) | Crossing of 56 MHz | Ramp length | Crossings per block |
|---|---|---|---|---|
| A1 (50 Hz) | 10.4 | 5.4 ms | 10 ms | 3–4 (1–2 up, 1–2 down) |
| A2 (150 Hz) | 39.8 | 1.41 ms | 3.33 ms | ~10 |
| A3 (12–1623 Hz) | 3.2–432 | 17.5 ms–0.13 ms | 41.7 ms–0.31 ms | 0–110 |
| D1 (143 Hz) | 42.9 | 1.31 ms | 3.5 ms | ~10 |
| D2 (1000 Hz) | 353.6 (triangle) or 176.8 (sawtooth) [U] | 0.16 / 0.32 ms | 0.5 / 1 ms | ~68 / ~34 |
| D3 (100–4700 Hz) | 26.6–1250 | 2.1 ms–45 µs | 5 ms–106 µs | 7–320 |

Consequences:
- **The sweep always extends beyond the capture.** B is 104–177 MHz, larger than 56 MHz, so the far turn-around
  points are not seen when the window sits mid-band. With fc at 4300 MHz the window covers 4272–4328 MHz, and both
  turns of an A1 sweep (4248–4352 MHz) lie outside.
- **The chirp is a pure CW tone sweeping linearly.** A dechirp test gives an almost ideal line: over 5.4 ms the
  time-bandwidth product is 56e6 × 5.4e-3 ≈ 3e5, so the processing gain is enormous.
- **Optimal STFT window per crossing:** T_w ≈ 1/√r.
  - A1: 1/√(1.04e10) ≈ 9.8 µs;
  - A2: 5 µs;
  - D3 at its fastest: 0.9 µs.

  A 64-point STFT (1.04 µs) suits fast types and smears slow ones only mildly: the peak stays within one bin over
  ~10 rows at A1.
- **Received level.** An aircraft at h km gives roughly −29 dBm EIRP-normalised at 1 km on the boresight, so on the
  ground it is only visible when the aircraft is near (approach path, a few km, side lobes) [U, to be measured on
  the bench]. Also expect ground test sets (ramp testers) at airports.

### 1.3 Measuring rate, period and shape from crossings

**Detector (ridge tracker, independent of the jammer fold):**
1. STFT (reuse the jammer's 64-bin STFT, or 256-point at 61.44 MSps for slow types). Per row, take the peak bin,
   its power over the bin's median floor, and a sub-bin parabolic interpolation.
2. Seed a track on rows where the peak is ≥ 8 dB over the floor.
3. Grow the track while the peak moves in one direction with a step consistent with the current slope (Kalman or
   a running LS fit, gate ±2 bins).
4. Keep a track that spans ≥ 0.5·W (≥ 28 MHz) with LS-fit R² ≥ 0.995.
5. **Refine the rate by dechirping.** For the track's central 80 %, y(t)·exp(−iπ r t²) for r in r0·(1 ± 2 %) on a
   fine grid. Maximise the FFT peak (golden section). Expect accuracy ~0.01 % on A1.

**Per crossing:**
- `t_c`: time at the window centre;
- sign s = ±1 (up / down);
- rate r̂;
- the in-capture extent f_min…f_max.

**Combining crossings** (sorted by t_c, with W the window and f_w the window centre frequency):
- **Sawtooth.** All s are equal. Period T = median(Δt_c) between same-sign crossings. Flyback: a dark interval, and
  no down crossing.
- **Triangle.** Signs alternate.
  - Same-sign spacing gives T.
  - An up crossing followed by a down crossing gives the top turn-around:
    f_top = f_w + r̂·(t_d − t_u)/2.
  - A down crossing followed by an up crossing gives the bottom:
    f_bot = f_w − r̂·(t_u − t_d)/2.
  - Then B = f_top − f_bot and T = 2B/r̂. With only one (up, down) pair, report f_top or f_bot alone, and
    T ≥ (t_d − t_u) + W/r̂.
  - Check: 2B·f_rep = r̂ should match within 2 %.
- **Too slow** (fewer than 2 crossings): report r̂ only, and bound T > block length − t_x.
- **Variable period** (A3, D3): T changes between blocks (altitude change), while r̂ stays proportional to 1/T at
  constant B. Report it across `StreamAnalyzer` blocks when available.
- **Several units on one aircraft:** parallel tracks with the same |r̂| (±1 %), offset in frequency by
  ~2.5 / 5 MHz (D1, D2, A1, A2) or offset in time. Report as one emitter group with `units = n`.

**Pulsed altimeters** (A4–A6, D4) go through pulse_train / PDW extraction (§4.4):
- PW 30–225 ns, i.e. 2–14 samples at 61.44 MSps; the 1 µs envelope of `pulse_train` is **too coarse**, so use the
  raw-rate envelope (§4.4);
- PRI 40–167 µs, so 200–850 pulses per block;
- 3 dB BW 5–31 MHz;
- band 4200–4400.

  Label "Pulsed radio altimeter-like". Distinguish it from a C-band radar: PW below 0.3 µs and PRF 6–25 kpps in
  this band have no other allocated user.

### 1.4 Altimeter or swept-carrier jammer

| Feature | Radio altimeter | Swept jammer |
|---|---|---|
| Band | 4200–4400, sweep inside 4300 ± (B/2 + drift) | anywhere: 2.4 / 5.8 GHz, cellular, GNSS |
| Rate r | 3–1300 MHz/ms (table above) | Typically 1–100 kHz sweep frequency over 20–100 MHz, so 40–20 000 MHz/ms (bench: 10 kHz over ~80 MHz ≈ 1600 MHz/ms) |
| Period | 0.2–83 ms | 10–1000 µs |
| Span vs window | > 100 MHz, always beyond 56 MHz | often inside the band it targets |
| Linearity | Very high: digital units are DDS-like, analogue are linearised VCOs | VCO with a sawtooth or triangle tuning voltage; nonlinear (curved trajectory) and jittered (~2 µs sweep-to-sweep on the bench unit) |
| Envelope | Constant CW | Constant, but often noise-modulated or with spurs |
| Presence | Intermittent (aircraft passing), level varies over seconds | Continuous while switched on |

**Rule:**
- In `radalt`, a track with r in [2, 1500] MHz/ms, residual from a line < 0.1 % of the span, and either the period
  in [0.2, 85] ms or too slow to see a second crossing → "Radio altimeter (FMCW)-like", confidence 0.7. Raise it to
  0.8 when ≥2 crossings give B̂ in 95–190 MHz, or when the period matches a table value.
- Otherwise in `radalt` keep the existing "Swept carrier (FMCW-like)".
- Outside `radalt` an FMCW-like sweep is the jammer receiver's job, but add a note: "FMCW radar outside the
  altimeter band would be called a jammer" (README limitation, unchanged).

### 1.5 Synthetic tests

Generator:
```
f(t) = fc0 + B/2 · tri(2π f_rep t + φ0)
x = A·exp(i·2π∫f) + noise,  f_rep from the table, φ0 random
```
- Window centre offsets 0, ±40 and ±80 MHz from 4300 MHz.
- SNR per MHz from 0 to 20 dB.
- Two units offset 5 MHz.
- An A3 with 12 Hz (one partial crossing).
- Pulsed: Gaussian-shaped 75–225 ns pulses at 6–25 kpps, with a BPSK chip code of 13 chips (Barker) [U, because
  "spread spectrum" is unspecified].

Negatives:
- the existing jammer synthetics (`receivers/jammer/synth.py`);
- an LFM radar pulse of 50 µs (track too short);
- noise.

---

## 2. ILS glide slope, 329.15–335.0 MHz (S9)

### 2.1 Signal (ICAO Annex 10 Vol I §3.1.3 localizer, §3.1.5 glide path) [M unless marked]

| Item | Localizer | Glide path (GP / GS) |
|---|---|---|
| Carrier | 108.10–111.95 MHz, odd tenths only (x.10, x.15, x.30, x.35 …), 40 channels | 329.15–335.00 MHz, 150 kHz grid, paired 1:1 (table §2.2) [V: band edges, NTIA compendium] |
| Polarisation | horizontal | horizontal |
| Modulation | A8W: carrier AM by 90 Hz and 150 Hz | same |
| Depth of each tone on course | **20 %** | **40 %** (tolerance 37.5–42.5 %) |
| Sum (SDM) | 40 % (≤ 60 % with ident / voice) | 80 % (the user's "40 % each" is the GP value) |
| Which tone dominates | 90 Hz on the **left** of course as seen from an approaching aircraft, 150 Hz on the right [V: Wikipedia localizer] | **90 Hz above the path, 150 Hz below** |
| DDM sensitivity | 0.155 at the course-sector half-width (sector set to ~210 m / 700 ft wide at the threshold) | 0.0875 at angular displacement 0.07θ–0.14θ (nominal 0.12θ) above / below the path, θ = path angle ~3° |
| Tone tolerance | 90 / 150 Hz ±2.5 % (Cat I), ±1.5 % (Cat II), ±1 % (Cat III) | same |
| Tone phase | 90 and 150 Hz harmonically locked to a 30 Hz reference, so the 150/90 ratio is exactly 5/3 | same |
| Ident | Morse at 1020 ± 50 Hz, depth 5–15 %, ~7 words/min, ≥ 6 times a minute; optional voice (ATIS) 300–3000 Hz | **none** (the localizer identifies the ILS) |
| Two-frequency (capture effect) | course and clearance carriers symmetric about the assigned frequency, separation 5–14 kHz | separation **4–32 kHz** |
| Frequency tolerance | ±0.005 % single; ±0.002 % per carrier for two-frequency [U] | ±0.005 % [U] |
| Coverage | 25 NM ±10°, 17 NM ±35° | 10 NM, ±8° azimuth, 0.45θ–1.75θ elevation; ≥ 400 µV/m |

**DDM definition (Annex 10 Vol I ch. 1):** "the percentage modulation depth of the larger signal minus the
percentage modulation depth of the smaller signal, divided by 100". It is unsigned in ICAO; the sign is the side.

**Library convention (signed):** `ddm = m90 − m150` (fractions).
- Glide path: ddm > 0 means the receiver is above the path ("fly down"); ddm < 0 means below ("fly up").
- Localizer: ddm > 0 means left of course as seen from the approach ("fly right").
- Also report `sdm = m90 + m150`.

For a ground receiver near an airport, the GP usually reads strongly 150-dominant (below the path), and the DDM is a
geometric fact about where the receiver sits; it says nothing about the beacon's health.

### 2.2 Localizer ↔ glide-path pairing (40 ILS channels)

Every x.x5 channel pairs with its x.x0 neighbour's GP frequency − 150 kHz.

| LOC | GS | LOC | GS | LOC | GS | LOC | GS |
|---|---|---|---|---|---|---|---|
| 108.10 | 334.70 | 109.10 | 331.40 | 110.10 | 334.40 | 111.10 | 331.70 |
| 108.15 | 334.55 | 109.15 | 331.25 | 110.15 | 334.25 | 111.15 | 331.55 |
| 108.30 | 334.10 | 109.30 | 332.00 | 110.30 | 335.00 | 111.30 | 332.30 |
| 108.35 | 333.95 | 109.35 | 331.85 | 110.35 | 334.85 | 111.35 | 332.15 |
| 108.50 | 329.90 | 109.50 | 332.60 | 110.50 | 329.60 | 111.50 | 332.90 |
| 108.55 | 329.75 | 109.55 | 332.45 | 110.55 | 329.45 | 111.55 | 332.75 |
| 108.70 | 330.50 | 109.70 | 333.20 | 110.70 | 330.20 | 111.70 | 333.50 |
| 108.75 | 330.35 | 109.75 | 333.05 | 110.75 | 330.05 | 111.75 | 333.35 |
| 108.90 | 329.30 | 109.90 | 333.80 | 110.90 | 330.80 | 111.90 | 331.10 |
| 108.95 | 329.15 | 109.95 | 333.65 | 110.95 | 330.65 | 111.95 | 330.95 |

Provenance:
- [V] for the x.x0 pairs (two independent web tables agree) and for 38/40 of the x.x5 pairs.
- 108.15 and 108.75 differ in one source (335.55 / 330.75), which are typos: 335.55 lies outside the GP band, and both
  break the −150 kHz rule. 334.55 and 330.35 are used here [U: confirm against the Annex 10 Vol I §3.1.6 table].

### 2.3 Tone estimator: 90 and 150 Hz are orthogonal over 1/30 s [V, arithmetic]

- Periods: 11.11 ms (90 Hz) and 6.67 ms (150 Hz). GCD(90, 150) = 30 Hz, so the common period is LCM = 33.33 ms. The
  window holds 3 and 5 whole cycles, so the DFT bins k = 3 and k = 5 of a 1/30 s window are exactly orthogonal to
  each other and to DC.
- A block (34.13 ms) holds one full 33.33 ms window plus 0.8 ms.
- Leakage when the tones are off nominal (±2.5 % Cat I tolerance): 150 Hz at +2.5 % lies 0.125 bin off bin 5. Its
  leakage into bin 3 at distance 2.125 bins is ≈ sin(π·0.125)/(π·2.125) ≈ 0.057 of its amplitude, about 2.3 % DDM
  bias at m150 = 0.4. Fixed DFT bins are therefore not enough for Cat I tolerances. Use a parametric fit.

**Recommended estimator (envelope, locked-ratio LS):**
```
y   = channel(x, f_carrier ± 1 kHz)          # decimate to fs_e = 6 kHz (N = 200 per 1/30 s)
a   = |y|                                    # envelope is CFO-invariant; ppm only moves the carrier
use the last N_w = floor(len(a) / (fs_e/30)) * (fs_e/30) samples (one full 1/30 s window in a block)
for s in linspace(0.97, 1.03, 61):           # common scale: tones are 90 s and 150 s (ratio locked to 5/3)
    H = [1, cos(2π90s t), sin(2π90s t), cos(2π150s t), sin(2π150s t)]
    θ = lstsq(H, a); keep s with min residual   (refine with a parabola on the 3 best)
a0 = θ0; m90 = hypot(θ1, θ2)/a0; m150 = hypot(θ3, θ4)/a0
ddm = m90 − m150; sdm = m90 + m150
σ_m ≈ sqrt(2/(N_w · SNR_env))   # report as uncertainty; SNR_env = a0² / var(residual)
```
- Use a squared-envelope variant (|y|²) only at high SNR. |y| is less biased by noise when SNR ≥ 10 dB; at lower SNR
  apply the Rician bias correction a0² ← a0² − σ_n².
- **Two-frequency systems:** the carriers are 4–32 kHz apart. Detect up to 3 carrier lines within ±20 kHz of the
  channel (the FFT of the block has 29 Hz resolution) and estimate each carrier separately with its own ±1.5 kHz
  filter (the tones need only ±150 Hz). The course carrier is usually the stronger near the path; label the carriers
  `course` and `clearance` only by relative level [U].
- **Tone frequencies:** report s·90 and s·150. |s − 1| > 0.03 means "not ILS".

### 2.4 Recipe

1. **Detector.** A carrier line in 328.6–335.4 MHz within ±(5 kHz + ppm·f) of a 150 kHz-grid channel from the
   table. Narrow (99 % power within ±2 kHz), with no FM voice structure (voice_fm rejects it, because its envelope
   varies).
2. **Measure:** carrier frequency and offset, the number of carriers and their spacing, m90, m150, DDM, SDM, tone
   scale s, and SNR.
3. **"Decode" (fields):**
   - channel (GS frequency → paired LOC frequency and the ILS channel number);
   - DDM and side ("above / below path");
   - SDM;
   - two-frequency yes/no.
4. **Validity:**
   - both tones present, each ≥ 5 σ;
   - SDM within 0.6–1.0 (GP) or 0.3–0.6 (an echo of the localizer; impossible at 330 MHz, so ignore);
   - s within ±3 %;
   - a 30 Hz-locked phase check: the 90 and 150 Hz phases fitted on the same t must satisfy the ICAO "zero
     crossings in the same direction" lock, φ150 − (5/3)·φ90 ≡ 0 (mod 2π/3) [U: verify the exact ICAO phase
     clause before using it as a hard gate].
5. **Label:**
   - "ILS glide slope-like" (confidence 0.8 with tones, the grid and the band);
   - add "(two-frequency)";
   - status `partial` (parameters measured, no payload);
   - a carrier on the grid without tones → `unknown`, with evidence.
6. **Synthetic:**
   ```
   x = A·(1 + m90·cos(2π90t+φ) + m150·cos(2π150t+(5/3)φ+ψ))·e^{i2πΔf t} + n
   ```
   with m90, m150 ∈ {(0.4, 0.4), (0.35, 0.45), (0.48, 0.32)}, tone scale 0.975–1.025, carrier offsets for
   ±10 ppm, two-frequency ±4 and ±16 kHz, SNR 0–30 dB.
   - Expected: DDM error < 0.005 at 20 dB.
   - Negatives: CW, NBFM voice, VOR-like 30 Hz AM, a 2-tone 100/200 Hz AM.

The pipeline input range starts at 300 MHz, so the GP is in range; the localizer (108–112 MHz) belongs to
`06_vhf_below_300.md` / S21. Pair the two when both are seen.

---

## 3. Mode A/C replies (1090 MHz) and interrogations (1030 MHz), military IFF (S11)

Main source: ICAO Annex 10 Vol IV §3.1.1 (SSR Mode A/C) [M]. Most numbers are also in RTCA DO-181 / EUROCAE ED-73
(Mode S transponder MOPS).

### 3.1 Reply format [M, widely reproduced]

| Pulse | Time after F1 (µs) | | Pulse | Time after F1 (µs) |
|---|---|---|---|---|
| F1 | 0.00 | | B1 | 11.60 |
| C1 | 1.45 | | D1 | 13.05 |
| A1 | 2.90 | | B2 | 14.50 |
| C2 | 4.35 | | D2 | 15.95 |
| A2 | 5.80 | | B4 | 17.40 |
| C4 | 7.25 | | D4 | 18.85 |
| A4 | 8.70 | | F2 | 20.30 |
| X | 10.15 (not transmitted in civil SSR; must be 0) | | SPI | 24.65 (= F2 + 4.35) |

- Pulse width 0.45 ± 0.10 µs; rise 0.05–0.1 µs; decay 0.05–0.2 µs.
- Position tolerance ±0.10 µs relative to F1 for F2 and the code pulses, ±0.15 µs between any two code pulses [M].
- Pulse-to-pulse amplitude within 1 dB of one another.
- Reply delay 3.0 ± 0.5 µs after the leading edge of P3.
- Carrier **1090 ± 3 MHz** for ATCRBS (Mode A/C-only) transponders; Mode S transponders ±1 MHz.
- Reply rate: SSR interrogations (A/C interlace) at PRF ~100–450 Hz per radar, so a transponder in a beam answers
  ~15–30 times per dwell. Uncorrelated replies (FRUIT) from other radars are common.

**Mode A code:** four octal digits ABCD with A = 4·A4 + 2·A2 + A1 (likewise B, C and D). There are 4096 codes.
Special codes: 7500 (unlawful interference), 7600 (radio failure), 7700 (emergency), 2000 (ICAO, no code assigned
on entering SSR coverage) [M], and 7000 (conspicuity in ICAO EUR/others) [U for India's VFR conspicuity code; 1200 is
US].

### 3.2 Mode C: Gillham altitude decoding [V, executed in the scratchpad]

- **Bits.** D2 D4 A1 A2 A4 B1 B2 B4 (D2 the MSB) are an 8-bit reflected Gray code n500 of 500 ft steps. C1 C2 C4
  form a 5-state code for 100 ft steps. D1 is 0 in all defined codes; D1 = 1 is invalid, and would only be needed
  above 126,750 ft.
- **100 ft part.** Read c = (C1, C2, C4) as a 3-bit pattern. In order these patterns map to n100 = 1…5:
  001 → 1, 011 → 2, 010 → 3, 110 → 4, 100 → 5. The patterns 000, 101 and 111 are **invalid**.
- **Reflection.** If n500 is odd, n100 ← 6 − n100.
- **Altitude:** alt_ft = 500·n500 + 100·n100 − 1300.
```
def gillham_alt(p):                          # p: dict pulse -> 0/1
    if p.D1 or (p.C1|p.C2|p.C4) == 0: return None
    g = bits(D2,D4,A1,A2,A4,B1,B2,B4)        # MSB first
    n500 = gray_to_binary(g)                 # b_i = XOR of g_j for j >= i
    n100 = {0b001:1, 0b011:2, 0b010:3, 0b110:4, 0b100:5}.get(bits(C1,C2,C4))
    if n100 is None: return None
    if n500 & 1: n100 = 6 - n100
    return 500*n500 + 100*n100 - 1300
```
**Check run:**
- All 4096 codes give **1280 valid codes**, mapping one-to-one onto **−1200 ft … +126,700 ft in contiguous 100 ft
  steps**.
- Known points agree with the published table:

  | Altitude | Code (D1D2D4 A1A2A4 B1B2B4 C1C2C4) |
  |---|---|
  | −1200 ft | 000 000 000 001 |
  | −1000 ft | 000 000 000 010 |
  | 0 ft | 000 000 011 010 |
  | 1000 ft | 000 000 110 010 |
  | 126,700 ft | 010 000 000 001 |
- The ICAO table begins at −1000 ft. Report −1200 and −1100 as `valid=False, reason="below ICAO table"` [U whether
  transponders ever emit them].
- The algorithm is the standard one (also in dump1090, GPL, which was not used). The pseudo-code above is written
  from the encoding definition, so it is free for the C++ port.
- Resolution 100 ft; the encoder's own rounding is ±50 ft.

### 3.3 Mode A or Mode C without seeing the interrogation

Both replies have the same format, so a single reply is ambiguous. Evidence, in order:
1. **C1 = C2 = C4 = 0** → Mode A only: no altitude has c = 000. Likewise D1 = 1 or c ∈ {101, 111}.
2. **Plausibility.** Gillham-valid and in −1000…60,000 ft is about 620 of 4096 codes (15 %). A random Mode A code
   passes about 15 % of the time, so a valid altitude alone is weak evidence.
3. **Pairing per transponder within the block:**
   - Cluster the replies by carrier offset. The ATCRBS ±3 MHz tolerance makes transponders separable: estimate the
     frequency of every pulse from its phase slope (±20 kHz at 0.45 µs and 20 dB), then the mean per reply. Then
     cluster by amplitude (±1.5 dB within a dwell).
   - In an interlaced A/C dwell one transponder gives **two** stable codes. The one that is Gillham-valid and
     consistent is C; the other is A. If both are valid, report both interpretations.
4. **Across blocks** (StreamAnalyzer or a sweep repeat): the altitude changes ≤ ~6000 ft/min, while the Mode A code
   is constant.
5. **Mode S cross-check.** If the adsb receiver decoded DF5 (ID13) or DF4/DF20 (AC13 with M=0, Q=0) in the same
   block, an equal squawk or altitude resolves it. DF17 TC28 gives the squawk too.
6. **1030 MHz in the block.** Mode A interrogations are P1–P3 = 8 µs and Mode C 21 µs. Pair them by a constant lag:
   reply time − interrogation time = 3 µs + (d_IT + d_TR − d_IR)/c is constant (±0.2 µs) over a dwell, for one
   interrogator and transponder.

   The two frequencies are 60 MHz apart, so they only fit at fc = 1060 MHz with 61.44 MSps, and both lie at ±30 MHz,
   **outside the 56 MHz analog bandwidth** (attenuated, but strong interrogators remain visible) [U: measure the
   Pluto roll-off at ±30 MHz].

**Output per reply:**
- `code_octal`, `ident` (SPI), `x_pulse`;
- `mode`: "A", "C", or "A_or_C";
- `altitude_ft` when Gillham-valid;
- `emergency` (7500 / 7600 / 7700, only when the mode is A or A_or_C);
- `f_offset_hz`, `amplitude_db`, `t_s`;
- flags `garbled` / `phantom_risk`.

**Labels:**
- "SSR Mode A/C reply" (confidence 0.85, status `success`) when the framing passes all tolerance checks: F1/F2 20.3
  ± 0.1 µs, pulse widths, the 1.45 µs grid, X = 0, and amplitudes within 2 dB.
- Per transponder: "Mode A code 1234" and "Mode C 10 300 ft" when resolved.

### 3.4 Garbling and phantoms

- **Garble.** Two replies whose brackets overlap (start times < 20.3 µs + 0.45 µs apart). Code pulses from one fall
  on the other's grid when Δt is a multiple of 1.45 µs (synchronous garble, which cannot be separated by timing).
- **Phantom.** A false F1…F2 pair at 20.3 µs formed by F2 of one reply and a code pulse of another, or by two code
  pulses of interleaved replies; or a bracket inside a Mode S reply (PPM chips at 0.5 µs multiples) or inside a DME
  pair train.

**Degarbling recipe (fs ≥ 10 MSps around 1090 ± 4 MHz):**
1. Extract pulses (PDWs): leading edge at the 50 % (−6 dB) crossing, width, amplitude, and the frequency from the
   phase slope inside the pulse.
2. Mask the time spans of CRC-verified Mode S frames, and pulses of width > 0.8 µs (DME 3.5 µs, Mode S
   interrogations).
3. Candidate brackets: pairs (i, j) with t_j − t_i = 20.3 ± 0.1 µs and |A_i − A_j| ≤ 2 dB.
4. Fill the code slots at 1.45·k ± 0.1 µs **using only pulses within ±2 dB of the bracket amplitude and within
   ±200 kHz of its frequency**. A pulse that fits two brackets is shared and flagged `garbled`.
5. **Reject when:**
   - a pulse lies between grid positions at 0.725 µs offsets (for example 1.45k + 0.725) with the bracket's
     amplitude and frequency;
   - X is set;
   - more than one candidate bracket shares F1 or F2 and both fit (the phantom), unless the frequency separates
     them.
6. **SPI** is a pulse at F2 + 4.35 µs matching amplitude and frequency.
7. Confidence: lower to 0.6 if any slot was shared.

The per-pulse carrier frequency is the IQ advantage over classic video-only degarblers. Spread across transponders is
up to ±3 MHz; the measurement precision is σ_f ≈ √(6/(SNR·T³))/(2π·√fs) per pulse. At 20 dB, 0.45 µs and 61.44 MSps
this gives tens of kHz [U, measure on synthetic data].

### 3.5 Interrogations at 1030 MHz [M]

| Mode | P1–P3 spacing | Use |
|---|---|---|
| 1 | 3 ± 0.1 µs | military: mission / aircraft type |
| 2 | 5 ± 0.2 µs | military: unit / tail |
| 3/A | 8 ± 0.2 µs | identity (civil Mode A) |
| B | 17 µs | rarely used |
| C | 21 ± 0.2 µs | altitude |
| D | 25 µs | unassigned |

(The task brief swapped Modes 1 and 2: Mode 1 is 3 µs and Mode 2 is 5 µs.)
- P1 and P3 are 0.8 ± 0.1 µs. The SLS P2 is 2.0 µs after P1 and at control (omni) level.
- Mode A/C/S all-call: P4 1.6 µs wide, 2 µs after P3. Mode A/C-only all-call: P4 0.8 µs.
- A Mode S interrogation is P1 and P2 (2 µs) followed by the P6 DPSK block.
- Interrogator PRF is ~100–450 Hz.
- **Recipe.** In the 1030 band, pulses of 0.7–0.9 µs; pairs at the spacings above (±0.2 µs); a P2 present at
  2.0 µs. Label "SSR interrogation Mode A/C/1/2-like" with the mode counts and the PRI (this extends the existing
  `r_adsb` 1030 rule). Interlace pattern: for example A, C, A, C or A, A, C.

### 3.6 Military IFF replies [M / U]

- **Mode 3/A** is civil Mode A (same 4096 codes).
- **Mode 2:** the same 12-bit reply format, 4096 codes.
- **Mode 1:** the same framing. Only 32 codes (2 octal digits) are used [U: which pulses carry it; commonly the A and
  B pulses, first digit 0–7 and second 0–3].
- Without the interrogation, Modes 1, 2 and 3/A replies are **indistinguishable from Mode A/C**. A 3 or 5 µs
  interrogation at 1030 in the same block is the only cue. Report the military mode only when an interrogation
  pairing (§3.3 item 6) exists; otherwise A_or_C.
- **Modes 4 and 5** are cryptographic. Mode 5 uses Mode S-like waveforms with spread / MSK elements; the public
  details are insufficient. Do not attempt decoding. An unknown, non-ATCRBS, non-Mode S pulse structure at 1030/1090
  gets `unknown`, with evidence [U].

### 3.7 Synthetic tests

- Generate replies with 0.45 µs trapezoid pulses (0.07 µs edges) on the grid, a random code, carrier offsets
  U(−3, 3) MHz, and 10–30 dB.
- Garble two replies offset by k·1.45 µs (k = 1…13) and by non-grid offsets.
- Phantom: replies 20.3 µs apart such that F2 of one coincides with F1 of the next.
- Mix with Mode S DF17 frames (reuse the adsb synth if present) and DME pairs at 1090 ± 30 MHz.
- Expected: no reply decoded from DF17 interiors; correct codes on non-grid garble; flags on grid garble.
- Gillham: all 1280 valid codes round-trip; the 2816 others are rejected.

---

## 4. Radars in India, 300 MHz–6 GHz (S10)

### 4.1 Allocation and emitter classes

| Band (MHz) | `bands_india` key | Typical emitters (India) | PW | PRF / PRI | Scan | Notes / source |
|---|---|---|---|---|---|---|
| 420–450 | `govt_mobile_435` | radiolocation (military), UHF wind profilers in 440–450 (WRC-97 Res. 217) [U for India] | 1–4 µs, phase-coded | 5–25 kHz | fixed beams | ITU-R M.1227 [M]. "404 MHz" profilers are the US NOAA network, not India [U]. |
| 1215–1400 | `rnss_l2_radar`, `arns_radar`, `fixed_radar_1350` | en-route ATC (AAI long-range ARSR), air-defence surveillance, 1280–1290 MHz boundary-layer wind profilers (NARL Gadanki LAWP 1280 MHz [U]) | ATC: 1–2 µs short + 50–150 µs LFM/NLFM long; profiler 0.5–3 µs coded | ATC 250–1000 Hz (staggered); profiler 5–25 kHz | ATC 5–6 rpm (10–12 s) | ITU-R M.1463 (not re-read) [M] |
| 2700–2900 | `radar_s_atc` | **ATC ASR** (AAI terminal radars), **IMD S-band DWR** (BEL/ISRO DWR Mk II, 9 m dish, dual-pol) | ASR: 0.4–1.03 µs (tube) or 1 µs + 20–89 µs LFM/NLFM (solid state, chirp 1.3–3 MHz); DWR 0.5–4.7 µs (IMD tender 0.5/1/2/4 µs) | ASR 722–1172 pps, staggered (for example 10/8 batches at two PRFs every 1.4°); DWR 250–1300 pps, dual / staggered / interleaved PRF | ASR 12–15 rpm (4–5 s); DWR 1–6 rpm, volume scan 6–15 min | M.1464-2 Table 1 Radars A–F2 [V]; M.1849-2 Table 6 [V]; IMD tender [V: PW only] |
| 2900–3100 | `radar_s_marine` | shipborne / coastal S-band navigation radar (IMO), coast guard VTS | magnetron 0.05–1.2 µs; solid-state up to ~100 µs with chirp | 500–3000 Hz (short pulses have the higher PRF) | 20–48 rpm (1.2–3 s) | ITU-R M.1313 / M.1460 [M] |
| 2700–3400 | `radar_s_def` | naval / air-defence / multifunction (M.1464 Table 2: Radars I–M) | 0.1–1000 µs | 300 Hz–20 kHz | 30–360°/s or electronic | M.1464-2 Table 2 [V]: compression up to 20 000, chirp up to 100 MHz, frequency agile |
| 4200–4400 | `radalt` | radio altimeters | §1 | §1 | – | M.2059 [V] |
| 5250–5850 | `radar_5400`, `unii_2e` | **IMD C-band DWR** (5600–5650), military 5.4 GHz radars, airport / TDWR-type [U for India] | met: 0.5–3.3 µs (operational); compression ~40 µs; TDWR-like 1.1 µs | met 250–1200 (up to 2400); TDWR-like 2000 | 1–6 rpm (met), 21–24°/s TDWR | M.1849-2 Table 7 [V]; DFS in Wi-Fi [M] |

From M.1849-2 [V]:
- "Most ground meteorological radars operate within 5600–5650 MHz", generally limited to 5430–5725.
- European C-band survey: PW 0.5–2.5 µs, PRF 250–1200 Hz (up to 2400), rotation 1–6 rpm.
- Schemes: fixed, staggered, double- and triple-interleaved PRF.

**What one block can and cannot show:**
- **The scan period is not measurable in a single 34 ms block.** It needs ≥2 main-beam passes over ≥ 2–15 s, from
  StreamAnalyzer or sweep repeats.
- **Main-beam illumination is rare.** Beam dwell = beamwidth / scan rate:
  - ASR: 1.4° at 72–90°/s is ~16–19 ms;
  - DWR: 1° at 36°/s is ~28 ms.

  So the chance that a given block catches the main beam is ~(dwell + 34 ms)/scan period ≈ 1–2 % (ASR) and lower
  for DWR volume scans.
- **Usually the block sees side and back lobes** (−25 to −40 dB) of nearby radars, continuously. A 34 ms block then
  holds:
  - ASR: ~25–40 pulses;
  - DWR: ~8–40;
  - marine: ~17–100;
  - profilers, MFR and pulsed altimeters: hundreds.

**Label rules** (extend `r_pulsed`; every label "-like"; `technology=pulsed_radar`, `fields.radar_class`):

| Condition (band ∧ measurements) | Label | Conf. |
|---|---|---|
| 2700–2900 ∧ (PW ≤ 1.2 µs paired with a 20–100 µs LFM/NLFM pulse in the same PRI) ∧ PRI 0.83–1.4 ms (± stagger) | "Pulsed radar (S-band ATC, solid-state)-like" | 0.7 |
| 2700–2900 ∧ PW 0.4–1.2 µs single ∧ PRI 0.83–1.4 ms ∧ stagger | "Pulsed radar (S-band ATC)-like" (existing) | 0.65 |
| 2700–3000 ∧ PW 0.5–4.7 µs ∧ PRI 0.77–4 ms (250–1300 Hz) ∧ (fixed, dual or staggered PRF ratio 2:3, 3:4 or 4:5) | "Pulsed radar (weather, S-band)-like" | 0.6 |
| 5600–5650 ∧ PW 0.5–3.3 µs ∧ PRI 0.42–4 ms | "Pulsed radar (weather, C-band)-like" | 0.65 |
| 5250–5850 other ∧ PRI | "Pulsed radar (C-band)-like" | 0.55 |
| 2900–3100 ∧ PW 0.05–1.2 µs ∧ PRI 0.33–2 ms | "Pulsed radar (marine S-band)-like" | 0.6 |
| 1215–1400 ∧ PW ≥ 1 µs (+ long LFM) ∧ PRI 1–4 ms | "Pulsed radar (L-band surveillance)-like" | 0.6 |
| 1270–1295 / 440–450 ∧ PRI 40–200 µs ∧ PW 0.5–4 µs ∧ phase-coded | "Wind profiler radar-like" | 0.55 |
| 4200–4400 ∧ PW < 0.3 µs ∧ PRF 6–25 kpps | "Pulsed radio altimeter-like" | 0.65 |
| else in a radar band, ≥ 3 pulses with a PRI | "Pulsed radar-like" (existing) | 0.5 |

The weather and ATC ranges overlap in 2700–2900. Prefer:
- "weather" when PW > 1.5 µs or when a dual-PRF ratio is found;
- "ATC" when a short pulse plus an LFM long pulse appear, or the stagger has more than 2 levels.

Otherwise give the generic label and list both candidates in the evidence.

### 4.2 Deinterleaving: algorithms

Input: the PDW list of a block, (TOA, PW, amplitude, frequency, intra-pulse slope), N up to ~2000.

**Step 0: pre-cluster before any PRI analysis** (literature consensus; Wiley *ELINT*):
- Frequency: separate by the detected signal box. Inside a box, cluster per-pulse frequencies with tolerance
  max(0.5 MHz, 1/PW).
- PW: log-ratio within 30 %, or ±0.2 µs for PW < 1 µs.
- Intra-pulse modulation: CW vs LFM (sign of the slope).
- Amplitude: use only loosely (±10 dB). A scanning beam modulates amplitude by tens of dB within a dwell, so
  amplitude separates emitters only within lobes.

**CDIF** (cumulative difference histogram; Mardia 1989, IEE Proc. F 136(4)) [M]:
- For difference levels c = 1, 2, …, histogram t_{n+c} − t_n and **accumulate** across levels.
- After each level, test the bins in ascending τ against T(τ) = x·(E − c)·e^{−τ/(kN)}, where E is the number of
  pulses, N the number of bins, and x, k empirical constants (x ≈ 0.3–0.5, k ≈ 0.2–0.4 [U: tune on synthetic
  data]). Also require the bin at 2τ to pass (subharmonic check).
- On a hit, run a **sequence search**: start at a pulse and step by τ ± tolerance, demanding ≥ 3–5 consecutive hits.
  Remove the extracted pulses and restart at c = 1.

**SDIF** (sequential difference histogram; Milojević & Popović 1992, IEE Proc. F 139(1)) [M]:
- Same as CDIF but **without accumulation**: one histogram per level c.
- Threshold T(τ) = x·(E − c)·e^{−τ/(kN)}.
- If only one bin passes at a level, take it as a PRI candidate and sequence-search immediately. That is faster and
  fewer harmonics are picked.

**PRI transform** (Nishiguchi & Kobayashi 2000, IEEE TAES 36(2)) [M]:
- D(τ_k) = Σ_{m<n, t_n−t_m ∈ bin k} exp(i·2π·t_n/(t_n − t_m)).
- |D| suppresses subharmonics (true-PRI phases add, harmonic phases cancel). The modified version (bin-dependent
  phase shift and overlapping bins) tolerates jitter up to ~±15 %.
- Cost O(N²) over pairs with a difference below τ_max: fine for N ≤ 1000 in C++ (≈5e5 pairs).

**Staggered PRI:**
- A K-level stagger shows K + 1 histogram peaks: the K sub-intervals and their sum, the frame period P = Σ Δ_i.
- Sequence-search on P first (the most periodic), then read the sub-interval pattern within each frame: the ordered
  Δ sequence, its period and the stagger ratio. M.1464 ATC example [V]: 10 pulses at one PRF and 8 at a lower one
  per 1.4° (block stagger). In 34 ms that is ~1 block change, visible as two PRI clusters with a step change, not an
  interleave.
- **Dual-PRF weather radars** switch per radial or per elevation. Within 34 ms, expect either one PRF, or two
  interleaved PRIs with a ratio of 2:3 / 3:4 / 4:5.

**Jitter:** a broad histogram peak (> 5 % width). Report `pri_jitter = σ/mean`. The PRI transform's jitter mode finds
it; CDIF at a wide tolerance also works.

**Achievable in one block:**
- **Constant or staggered PRI with ≥ 8 pulses per emitter:** reliable for 2–3 interleaved emitters that are
  separable by frequency or PW.
- **Two emitters on the same frequency and PW:** SDIF with ≥ 10 pulses each.
- **Jitter > 10 % with fewer than 20 pulses:** unreliable. Report a PRI range, not an emitter.
- **Fewer than 4 pulses** (main-beam edge): report the pulses only.

### 4.3 Recommended simple algorithm (for both Python and C++, deterministic)

```
pdws = extract_pdws(block)                  # §4.4
groups = cluster(pdws, key=(signal_box, freq±max(0.5MHz,1/PW), log PW ±30 %, CW/LFM))
for g in groups (sorted by size desc):
    r = g.toa (sorted); emitters = []
    while len(r) >= 4:
        # SDIF levels c = 1..4
        for c in 1..min(4, len(r)-1):
            h = histogram(r[c:] - r[:-c], bins = 0.25 µs from 20 µs to min(20 ms, span/3))
            thr = 0.4 * (len(r) - c) * exp(-τ / (0.3 * nbins))       # x=0.4, k=0.3 [tune]
            cand = [τ for τ in peaks(h) if h[τ] >= max(thr(τ), 3)]
            if cand: break
        if not cand: break
        τ = min(cand) refined by median of member differences
        seq = sequence_search(r, τ, tol=max(1 µs, 0.02 τ), min_run=4, allow_missing=2)
        if not seq:
            # stagger test: does the frame period P (sum of 2..4 consecutive Δ) recur?
            seq, pattern = stagger_search(r, cand, max_levels=8)
        if not seq: break
        emitters.append(describe(seq))      # PRI, stagger pattern, jitter, PW, f, LFM, n
        r = r \ seq
    leftover -> "pulses without PRI"
```
**Thresholds:**
- PDW detection: 10 dB over the median (as now), plus ≥ 2 samples above threshold.
- TOA from the −6 dB leading edge, interpolated (σ ≈ rise/√SNR, ~20–50 ns).
- PRI tolerance max(1 µs, 2 %); sequence acceptance ≥ 4 consecutive pulses or ≥ 60 % of the expected count across
  the span.
- Pulses per emitter for a confident label: ≥ 6.

**Integrate with `pulse_train`:**
- keep its group features for the existing rules;
- add `pdws` / `emitters` as new fields (an additive schema change);
- run the deinterleaver only where `pulsed_band` is true and N ≥ 6;
- latency is ~N log N plus the histograms, < 2 ms for N = 1000 in C++.

### 4.4 PDW extraction at full rate

The 1 µs envelope misses PW < 1 µs (altimeters at 30–225 ns, marine radars at 50 ns) and merges a DME pair.

Per signal box:
1. Channelize to ≥ 2·BW (minimum 10 MSps). Envelope |y|², smoothed by a moving average of max(1 sample, 0.1 µs).
2. Hysteresis threshold: on at +10 dB over the median, off at +6 dB.
3. Per pulse:
   - TOA (−6 dB edge, linear interpolation);
   - PW (−6 dB edges);
   - peak amplitude;
   - frequency (weighted phase slope over the central 80 %);
   - LFM slope and R² (existing fit, applied per pulse ≥ 5 µs);
   - optional phase-code flag: phase jumps > 60° at regular intervals (Barker-13 / complementary codes on profilers
     and pulsed altimeters).
4. Merge pulses closer than 0.2 µs (ringing).

### 4.5 Synthetic tests

Build pulse trains from the table rows:
- an ASR with 1 µs + 55 µs NLFM per PRI and 5-level stagger 970–1170 pps;
- a DWR at 1 µs with dual PRF 900/1200 interleaved, or 4.7 µs at 318 pps;
- a marine radar at 0.08 µs, 2000 pps;
- a C-band DWR at 2 µs, 1200 pps;
- a profiler: 1 µs Barker-13, PRI 100 µs.

Conditions:
- amplitude envelope of a scanning sinc² beam plus a −30 dB sidelobe floor;
- missing pulses at 10 %;
- TOA jitter 50 ns;
- 2–3 emitters interleaved in one block.

Expected:
- PRI within 0.5 %;
- the stagger pattern recovered (up to rotation);
- no emitter from Wi-Fi bursts (unii_2e) or random bursts (Poisson TOAs at the same rate). The false-emitter rate
  on Poisson must be < 1 % per block; set x and k from this.

---

## 5. Jammers beyond swept carriers (S12)

### 5.1 Classes and how they look

| Class | Typical implementation | IQ statistics | Spectral shape | Distinguishing features |
|---|---|---|---|---|
| **Swept CW** (existing) | VCO + sawtooth/triangle | constant envelope | fold shows a moving line | existing `jammer` receiver |
| **Fast swept CW** (GNSS PPD, some cellular) | VCO, sawtooth, period 8–27 µs | constant envelope | flat-topped "comb" with lines at k/T_sweep | instantaneous frequency is a sawtooth; needs a sub-µs time resolution |
| **Noise-FM (noise-modulated VCO)** | VCO driven by filtered noise (± a slow sweep) | **constant envelope** (kurtosis ≈ 1) | smooth, roughly flat or Gaussian-shaped over the deviation | |x|² nearly constant while the bandwidth is wide: no legitimate wideband service is constant-envelope **and** random-IF (FSK has discrete IF levels) |
| **Barrage noise (DDS / AWGN upconverted)** | baseband noise source + PA | complex Gaussian: \|x\|² exponential, E\|x\|⁴/E²\|x\|² = 2 | flat plateau with sharp edges matching a target band | Gaussian, **no** CP correlation, **no** cyclic features, **no** spectral lines; stationary over the block |
| **Spot** | CW or narrow noise on one channel | CW: constant; narrow noise: Gaussian | a single line or a narrow plateau | occupies exactly one channel of a service; continuous; no modulation structure |
| **Pulsed / follower** | reactive; transmits on detection | bursty | follows the victim | out of scope (no reference) [U] |

### 5.2 GNSS personal privacy devices (PPD)

From Mitch et al. 2011 [V], 18 commercial jammers:
- **All** used a single swept tone (linear FM), mostly sawtooth.
- L1 sweep period: **8–12 µs for 15 of 18** (9 µs is typical); two were 26–27 µs; one 1 µs.
- Sweep range about L1: typically +3.6…+31 MHz above and −4.4…−31 MHz below, i.e. under 20 MHz total for most;
  some were not centred on L1 at all (for example +18.6/−4.6: they do not even cross 1575.42).
- Six also jammed L2 (sweep period 1–9 µs); none jammed L5.
- 7 of 18 showed aberrant behaviour: sweep ripples, frequency jumps, centre wander > 2.5 MHz in 41 µs.
- Range up to 6–9 km (isotropic estimate).

Kraus et al. 2011 (in-car jammer survey) found the same chirp class and modelled it for cancellation [M].

**Detector** (the jammer receiver's STFT is too coarse: 9 µs ≈ 9 rows):
1. In 1559–1610 MHz (and 1215–1250 for L2), the signal box has bandwidth ≥ 2 MHz and is continuous over the block.
2. Channelize the box to ≥ 2·span (≤ 61.44 MSps), take the instantaneous frequency IF(t) = arg(y[n]·y*[n−1])·fs/2π,
   and smooth it over 0.1 µs.
3. The autocorrelation of IF(t) over lags 1–100 µs peaks at T_sweep (expect 8–30 µs).
4. Fold IF on T, giving the trajectory. Take the slope(s): sawtooth (one ramp + flyback) or triangle.
5. Validity:
   - fold coherence ≥ 0.8 (fraction of IF variance explained by the fold);
   - span ≥ 1 MHz;
   - T in 0.5–100 µs.
6. Label "GNSS jammer (swept, PPD-like)", confidence 0.85. Report T, span, centre offset from 1575.42, covers L1
   yes/no, and shape.

**Any signal above the noise floor in 1575.42 ± 12 MHz is interference.** GNSS is about 20–30 dB below the thermal
noise in this bandwidth, so even without a sweep, report "GNSS-band interference" with the measured class (CW, noise,
swept).

Synthetic: a sawtooth sweep of span 10–40 MHz, period 8–27 µs, ±0.5 µs period ripple, centre offset ±10 MHz, and
SNR per MHz from −5 to +20 dB.

### 5.3 Cellular / prison / exam jammers (India)

- **Targets:** the downlinks (handset receive):
  - B28 DL 758–803;
  - B5 DL 869–894 (in practice 869–890);
  - B8 DL 925–960;
  - B3 DL 1805–1880;
  - B1 DL 2110–2170;
  - B40 2300–2400 (TDD);
  - B41 / n41 2496–2690 (TDD);
  - n78 3300–3670;
  - often also 2.4 / 5.8 GHz Wi-Fi.
  [V: band edges from `bands_india.py`; U: which bands real Indian units cover.]
- **Implementations** (vendor literature, not peer-reviewed [U]): one module per band, each either
  - a VCO swept by a sawtooth or triangle at 1–100 kHz across the whole DL band, often with added noise
    (noise-modulated sweep), or
  - a band-limited noise generator.

  Output 1–10 W per band.
- **Legal [U]:** in India only government agencies may procure and operate jammers (DoT / Cabinet Secretariat
  policy; prisons, security, some exam centres with permission). Privately, they are illegal. The label must stay
  descriptive.
- **Detection:**
  - A swept module triggers the existing jammer receiver. Add a band-match field: `covers_band = "B3 DL"` when the
    sweep extremes are within 2 MHz of the band edges.
  - A noise module goes to the barrage detector (§5.4).
  - Evidence of jamming rather than traffic: an LTE/NR cell's PSS, CRS or SSB disappears or is buried; a DL band
    plateau with no PSS/SSB detected although the band is normally occupied; |x|² Gaussian with no 66.7 µs (LTE) CP
    correlation.

### 5.4 Barrage noise detector (recipe)

Per signal box of bandwidth ≥ 1 MHz that is continuous (present in ≥ 95 % of 100 µs frames):
1. **Plateau.** PSD from the block FFT (reuse the channelizer):
   - flatness inside (90 % of bins within ±3 dB of the median);
   - edges falling ≥ 15 dB within 5 % of the width.
2. **Gaussianity.** On the box channelized at 1.25·BW:
   - k = E|x|⁴/(E|x|²)² should be in 2 ± 0.1 for Gaussian (constant envelope gives 1, OFDM gives ~1.9–2.0, so this
     **alone does not separate OFDM**);
   - also the I/Q marginal kurtosis 3 ± 0.2 and I/Q balance.
3. **No structure** (this separates OFDM, DVB-T, LTE, Wi-Fi and spread spectrum):
   - **CP correlation:** max over Tu ∈ {3.2 µs Wi-Fi, 66.67 µs LTE, 1/SCS for NR 15/30 kHz, DVB-T2 2k–32k at 8 MHz,
     and a blind search from 1 to 1000 µs} of the normalised |R(Tu)| < 3/√N_eff. OFDM gives ≥ 0.05–0.2.
   - **Cyclic spectrum:** no line in FFT(|x|²) or FFT(x²) above 10 dB over its median between 1 kHz and BW
     (symbol rates, chip rates, frame rates such as 1 kHz LTE subframes and 100 Hz SSB bursts).
   - **No narrow spectral lines** (pilots, carriers): max/median of the PSD < 6 dB inside the plateau.
   - **Stationarity:** power per 1 ms within ±1 dB across the block (TDD bands: no UL/DL gaps at 5/10 ms).
4. **Not self-made:**
   - ADC clipping fraction < 1e-4 (clipping spreads a strong signal into a wide pedestal; the existing normalised
     clipping check runs first);
   - not the Pluto's capture-locked USB self-noise (+4–5 dB floor from ~12 ms into each buffer, with 8 kHz
     microframe modulation; FINDINGS): reject if the plateau starts at a fixed buffer offset or carries an 8 kHz
     line in |x|²;
   - plateau ≥ 6 dB above the receiver floor of an empty neighbouring step (when a sweep is available).
5. **Label:**
   - "Jammer (barrage noise)-like", confidence 0.6. Raise to 0.75 when the plateau edges match an allocation's
     edges within 2 MHz (`covers_band`) or when a known service is visibly suppressed.
   - Noise-FM (k ≈ 1, wide, random IF with a Gaussian-like IF histogram): "Jammer (noise-FM)-like", 0.7.
   - Otherwise "Wideband noise (unidentified)" → `unknown`, with evidence.

**False-positive risks:**
- OFDM at low SNR (CP correlation weak), especially LTE TDD and DVB-T2 → step 3, plus the existing LTE/NR detectors
  run first and own their band.
- Spread spectrum (UMTS, GNSS-like, LoRa) → chip-rate cyclic lines in x², and LoRa chirps in IF.
- Microwave ovens (2.45 GHz, 50 Hz-gated noise-like) → the 50/100 Hz on/off gating: not stationary.
- Wideband noise from saturation → step 4.
- Analog FM video → line rates.
- Radar sidelobe clutter → pulsed, not stationary.

Synthetic:
- filtered complex AWGN plateaus of 5–60 MHz with edges at band limits;
- noise-FM: IF = filtered Gaussian, deviation 10–40 MHz;
- negatives: LTE 20 MHz at 0–10 dB SNR (CP present), Wi-Fi bursts, a DVB-T2 8 MHz synthetic, clipped CW (ADC
  rails), and a noise-only block with the USB self-noise step.

### 5.5 Spot jammers

A CW line or a narrow (< 1 MHz) Gaussian plateau, continuous over the block, exactly on a service's channel:
- the 2.4 GHz FHSS RC channel set;
- GNSS L1 ± 1 MHz;
- a cellular PSS channel;
- 433 / 868 / 915 MHz control links.

The existing CW rule reports the CW. Raise it to "Jammer (spot)-like" (confidence 0.55) only when one holds:
- (a) it sits in a GNSS band (always interference); or
- (b) it is ≥ 20 dB above the band's legitimate signals and continuous, on a channel of a protocol that never
  transmits CW.

Otherwise keep "CW carrier". This is a high false-positive risk: receiver clock lines and LO leakage are CW, so
always run the existing receiver-line check (`ctx.cache["receiver_lines"]`).

### 5.6 Drone jammers (counter-UAS)

- **Targets:** 2400–2483.5 and 5725–5850 (control and video), GNSS L1 1560–1610 (sometimes L2 / L5 and GLONASS),
  sometimes 433 / 868 / 915 MHz.
- **Implementations:** swept VCO or noise per band, handheld 10–30 W [U: no peer-reviewed waveform measurement
  found].
- **Detection:** the existing swept detector (bench `LIVE_JAMMER_SWEEP_2460MHz`: 10 kHz, up-down, ~40 MHz seen)
  plus §5.4. **Cross-band evidence:** a simultaneous swept or noise signal in GNSS L1 and in 2.4 / 5.8 GHz within
  one sweep step (needs the sweep context, not one block) → `fields.multi_band = true`, label "Jammer (counter-UAS
  multi-band)-like".

---

## 6. AM voice: UHF military air and VHF civil air (S13)

### 6.1 Signal facts

| Item | VHF civil air | UHF military air |
|---|---|---|
| Band | 118.000–136.975 MHz (below the library's 300 MHz input: S14 / file 06) | **225–400 MHz** (NATO "UHF air"); only 300–400 is inside the current input range |
| Emission | A3E DSB-AM, full carrier | A3E DSB-AM, full carrier; some Have Quick FH (AM per hop, ~25 kHz channels, hops many per second) [M] |
| Channel spacing | 25 kHz; 8.33 kHz in Europe (named "channels", for example 118.005 = 118.000 MHz at 8.33). India: 25 kHz [U] | 25 kHz [M] |
| Audio band | 300–2500 Hz at 25 kHz (350–2500 Hz at 8.33) [M] | similar [U] |
| Modulation depth | ground ≥ 85 % peak, aircraft ≥ 70 % peak (ICAO Annex 10 Vol III Pt II §2.2.1.4 / 2.3.1.4) [M] | – |
| Offset carrier (climax) | 25 kHz channels with 2–5 ground transmitters: 2-carrier ±5 kHz; 3-carrier 0 / ±7.3 kHz; 4-carrier ±2.5 / ±7.5 kHz; 5-carrier 0 / ±4 / ±8 kHz (Annex 10 Vol V) [M: verify the exact offsets] | – |
| Guard | 121.5 MHz | **243.0 MHz** (military guard, = 2 × 121.5) |
| Other users in 225–400 | – | military SATCOM UHF (240–270, 290–320 MHz), tactical FM / digital, Have Quick, **ILS GS 329–335**, **COSPAS 406** above |

**India [U]:** the IAF and Navy use UHF AM for air-to-ground ATC at military airfields (public AIP entries list UHF
frequencies for joint civil/military aerodromes). Details beyond the AIP lists are not public. The AIP (aim-india.aai.aero)
lists the VHF civil frequencies per aerodrome; use it for channel validation, not as a hard requirement.

### 6.2 Detecting AM voice (vs FM, digital, CW)

On a candidate at f_c (reuse voice_fm's candidates, or narrow peaks ≤ 25 kHz wide in the bands above), channelized
to ±10 kHz at 24–48 kSps:
1. **Carrier line.** For DSB-AM with m_rms ≤ 0.5 the carrier holds ≥ 1/(1 + m²/2) ≥ 67 % of the power. The block
   FFT's 29 Hz bins show a line ≥ 15 dB above the sidebands' per-bin level.
2. **Carrier phase.** Estimate the carrier phase φ(t) with a 2nd-order PLL, or the mean phase per 5 ms with a
   frequency fit. Then z = y·e^{−iφ̂}.
3. **The I/Q test** (the core discriminator). In DSB-AM the information lies only on I, so after phase alignment:
   - Q_ratio = var(Im z) / var(Re z − mean). AM gives ≈ noise level (Q_ratio → σ²/(σ² + A²m²)). FM and PM put
     energy in Q.
   - Equivalent symmetric-sideband check: |Y(f_c + f)| = |Y(f_c − f)| **and** the sidebands are conjugate-coherent.
     The cross-spectrum of the upper and the mirrored lower sideband has coherence > 0.8 in 300–3000 Hz.
4. **Envelope beyond noise.** With a(t) = |y| and SNR measured from the carrier vs the out-of-channel floor:
   cv² = var(a)/mean(a)²; the noise-only expectation is 1/(2·SNR).
   - Voice AM: cv ≥ max(0.12, 1.6/√(2·SNR)) (the same threshold voice_fm uses to *reject*);
   - the envelope spectrum carries power in 300–3000 Hz, ≥ 6 dB over 4–8 kHz.
5. **Speech.** Envelope activity per 20 ms frame (syllabic 2–8 Hz modulation of the 300–3000 Hz energy) and kurtosis
   (speech is leptokurtic, > 4).
6. **Not digital:**
   - no symbol-rate line in the squared derivative of the envelope or of z (VDL2 D8PSK 10.5 ksym/s, ACARS MSK 2400
     audio tones on AM — ACARS **is** AM with 1200/2400 Hz tones, so check a 2400 Hz line in the envelope plus the
     ACARS preamble; file 06 handles that);
   - not CW (cv ≈ noise and no envelope band content).

**Classifier:**

| AM voice | FM voice | CW | Digital |
|---|---|---|---|
| carrier line; Q_ratio < 0.2; cv > noise; envelope 300–3000 Hz | cv ≤ noise; discriminator in 300–3400 Hz; Q carries the signal | cv ≈ noise, no band content | symbol lines |

**Label:** "AM voice (air band, UHF military)-like" in 225–400, "AM voice (VHF air)-like" in 118–137, otherwise
"AM voice-like". Confidence 0.7, raised to 0.8 with speech activity in ≥ 30 % of frames. Report the offset-carrier
group when 2–5 carriers sit within ±8 kHz of one channel centre with a common audio, coherence > 0.5 between their
envelopes ("offset-carrier (climax) n=3").

### 6.3 Demodulation, squelch, quality

**Demod:**
- Coherent demodulation is preferred: s(t) = Re(z) − DC (6 dB better than the envelope at low SNR, and no
  envelope-detector threshold).
- When several offset carriers are in the channel, use the envelope instead (coherent demodulation on one carrier
  makes the others beat).
- Then a 300–2500 Hz (25 kHz) or 350–2500 Hz (8.33 kHz) band-pass, AGC, and decimation to 8 kSps (16-bit) for the
  GUI.
- Expose `AMDemodulator` alongside `VoiceDemodulator`, with the same interface for signal-vault.

**Squelch:** carrier power ≥ 6 dB over the floor within ±5 kHz. Aircraft transmit the carrier only while keyed, so
carrier presence equals the keyed time. Bridge gaps < 100 ms. Report `keyed_s` and the segments.

**Quality metric (no reference):**
- audio SNR ≈ (power 300–2500 Hz during speech frames) / (power 300–2500 Hz in non-speech frames of the keyed
  segment) or, without pauses, / (power in 2.8–3.5 kHz × bandwidth ratio);
- report `audio_snr_db`, `mod_depth_rms`, `mod_depth_peak` (99th percentile of |s|/A), and over-modulation (fraction
  of a(t) < 0.02·mean, i.e. carrier pinch-off).
- ITU-T P.563 is not used (licence and complexity).

### 6.4 Synthetic

- x = A·(1 + m·s(t))·e^{i(2πΔf t + φ_n(t))} + n, with s = band-limited speech (use the voice_fm synth's speech
  surrogate), m_peak 0.7–0.95, phase noise 1° rms, SNR 5–30 dB.
- Offset-carrier with 2–5 carriers at the Annex 10 offsets and ±50 Hz.
- Negatives: NBFM speech (voice_fm synth), CW, ACARS-like AM-MSK tones, 2-tone beats (the voice_fm EMI case), the
  ILS GS carrier (90/150 Hz only, no 300–3000 Hz content).
- Expected: AM 12/12 at ≥ 10 dB; FM never AM; GS never AM voice.

---

## 7. Integration notes

**Schema (additive, minor bump), technologies, family and status:**

| Signal | Technology | Family | Status |
|---|---|---|---|
| radio altimeter | `radio_altimeter` | `aviation` | `partial` |
| ILS glide slope | `ils_glide_slope` | `aviation` | `partial` |
| SSR Mode A/C | `ssr_mode_ac` | `aviation` | `success` when the framing and Gillham (or Mode A) checks pass |
| deinterleaved radar | `pulsed_radar` with `fields.radar_class` | | |
| barrage jammer | `noise_jammer` | `jamming` | |
| GNSS PPD | `gnss_jammer` | `jamming` | |
| AM voice | `am_voice` | `voice` | |

Each one is `research` in `capabilities.json` with `sources: research/signal_gaps/05_aviation_radar_jammers.md`.

**Receivers:**
- extend `jammer` with the PPD and barrage detectors;
- new `aviation` receiver for the altimeter, GS and Mode A/C, or Mode A/C inside `adsb` (it shares the 1090
  channelizer and the Mode S masks; preferred);
- the radar deinterleaver lives in `features/` (pipeline-wide), the labels in `classify/rules.py`;
- AM voice: a new `voice_am` receiver, or inside `voice_fm` (it shares candidates; the plugin decides by the I/Q
  test).

**Bands:** add `ils_gs` (328.6–335.4) and `uhf_mil_air` (300–400 within range). Make `bands_india` list
`govt_mobile_*` overlapping with 225–400 when S14 extends the range.

**Latency targets (C++):**
- altimeter ridge tracker: reuse the STFT, < 2 ms;
- GS: one channelizer slice + a 61-point LS, < 1 ms;
- Mode A/C: per 1090 burst, < 1 ms;
- deinterleave: < 2 ms;
- barrage: CP correlation over ~20 lags plus 2 FFTs per box, ~5 ms per wideband box (gate on plateau flatness
  first);
- AM: per candidate, < 1 ms.

**Licences:** everything is implementable from the standards and papers below. No third-party code is needed beyond
pocketfft. dump1090 (GPL-2+) and pyModeS (GPL-3) are only oracles for Mode S / Gillham cross-checks.

## 8. Open items (to confirm)

1. The ITU-R M.2059 "waveform repetition frequency" means a full triangle period; whether D2 (with its hop pattern)
   is a sawtooth.
2. The Annex 10 Vol I §3.1.6 pairing for 108.15 and 108.75 (334.55 / 330.35 assumed). The GP phase-lock clause and
   the two-frequency tolerances.
3. The Annex 10 Vol V offset-carrier offsets. The AM modulation-depth minima (85 % ground, 70 % air). India's 8.33 kHz
   status.
4. The Mode 1 reply pulse assignment (32 codes). India's VFR conspicuity code.
5. The CDIF / SDIF constants x and k: tune on Poisson negatives.
6. Indian radar specifics (AAI ASR/ARSR models and PRFs; IMD DWR frequencies per site; any TDWR at Indian airports;
   Indian wind profiler frequencies).
7. Real-world waveforms of Indian cellular and drone jammers: no measurement literature found; a bench capture is
   needed.
8. Pluto roll-off at ±30 MHz (for 1030 + 1090 in one block at fc = 1060 MHz).

## Sources

| # | Source | URL | Licence / access |
|---|---|---|---|
| 1 | ITU-R M.2059-0 (02/2014), *Operational and technical characteristics and protection criteria of radio altimeters utilizing the band 4 200–4 400 MHz*. Tables 1–2 and footnote 6 extracted and read. | https://www.itu.int/rec/R-REC-M.2059-0-201402-I | © ITU, free download, facts only |
| 2 | ITU-R M.1464-2 (02/2015), *Characteristics of radiolocation radars, and characteristics and protection criteria … aeronautical radionavigation and meteorological radars … 2 700–2 900 MHz*. Tables 1–2 read. | https://www.itu.int/rec/R-REC-M.1464-2-201502-I | © ITU, free |
| 3 | ITU-R M.1849-2 (01/2019), *Technical and operational aspects of ground-based meteorological radars*. §4.1 and Tables 6–7 read. | https://www.itu.int/rec/R-REC-M.1849-2-201901-S | © ITU, free |
| 4 | ITU-R M.1638-1 (01/2015), radars in 5 250–5 850 MHz (downloaded, not tabulated here) | https://www.itu.int/rec/R-REC-M.1638 | © ITU, free |
| 5 | ITU-R M.1463 (1 215–1 400 MHz radars), M.1460 / M.1313 (2 900–3 100 MHz and maritime radars), M.1227 (wind profilers ~1 GHz). Cited, not re-read. | https://www.itu.int/rec/R-REC-M.1463 · https://www.itu.int/rec/R-REC-M.1460 · https://www.itu.int/rec/R-REC-M.1313 · https://www.itu.int/rec/R-REC-M.1227 | © ITU, free |
| 6 | ICAO Annex 10 Vol I (radio navigation aids: ILS §3.1), Vol III (communication systems: VHF AM), Vol IV (surveillance: SSR §3.1.1), Vol V (frequency utilisation: offset carriers) | https://store.icao.int (Annex 10) | © ICAO, paid; facts only |
| 7 | RTCA DO-181 / EUROCAE ED-73 (Mode S / ATCRBS transponder MOPS); RTCA DO-155 / EUROCAE ED-30 (radio altimeter MOPS) | https://www.rtca.org · https://www.eurocae.net | paid; facts only |
| 8 | NTIA Federal Spectrum Use compendium 328.6–335.4 MHz (glide slope band) | https://www.ntia.doc.gov/files/ntia/publications/compendium/0328.60-0335.40_01DEC15.pdf | US Gov, public domain |
| 9 | ILS LOC/GS pairing lists (cross-checked; two typos noted) | https://hfunderground.com/wiki/index.php/ILS · https://www.flightsimbooks.com/flightsimhandbook/CHAPTER_11_10_Navigation_Radios.php | web pages; facts only |
| 10 | Wikipedia, *Instrument landing system localizer* (90 Hz left / 150 Hz right); *Gillham code* (bit layout, example codes −1200…126 700 ft) | https://en.wikipedia.org/wiki/Instrument_landing_system_localizer · https://en.wikipedia.org/wiki/Gillham_code | CC BY-SA (facts only) |
| 11 | Mitch, Dougherty, Psiaki, Powell, O'Hanlon, Bhatti, Humphreys, *Signal Characteristics of Civil GPS Jammers*, ION GNSS 2011. Table 1 read. | https://radionavlab.ae.utexas.edu/images/stories/files/papers/Signal%20Characteristics%20Civil%20GPS%20Jammers_Mitch.pdf | © ION/authors; author copy |
| 12 | Kraus, Bauernfeind, Eissfeller, *Survey of In-Car Jammers – Analysis and Modeling of the RF Signals and IF Samples*, ION GNSS 2011 | https://www.ion.org/publications/abstract.cfm?articleID=9744 [U: URL id] | © ION, paid |
| 13 | Borio, Dovis, Kuusniemi, Lo Presti, *Impact and Detection of GNSS Jammers on Consumer Grade Satellite Navigation Receivers*, Proc. IEEE 104(6), 2016 | https://doi.org/10.1109/JPROC.2016.2543266 | © IEEE |
| 14 | Mardia, *New techniques for the deinterleaving of repetitive sequences*, IEE Proc. F 136(4):149–154, 1989 (CDIF) | https://doi.org/10.1049/ip-f-2.1989.0025 | © IET |
| 15 | Milojević & Popović, *Improved algorithm for the deinterleaving of radar pulses*, IEE Proc. F 139(1):98–104, 1992 (SDIF) | https://doi.org/10.1049/ip-f-2.1992.0012 | © IET |
| 16 | Nishiguchi & Kobayashi, *Improved algorithm for estimating pulse repetition intervals*, IEEE TAES 36(2):407–421, 2000 (PRI transform) | https://doi.org/10.1109/7.845217 | © IEEE |
| 17 | R. G. Wiley, *ELINT: The Interception and Analysis of Radar Signals*, Artech House, 2006 (PDW clustering, stagger / jitter analysis) | ISBN 978-1-58053-925-2 | © Artech, book |
| 18 | IMD tender pre-bid minutes, S-band polarimetric DWR (pulse widths 0.5 / 1 / 2 / 4 µs selectable) | https://internal.imd.gov.in/tenders/20250124_prebid_309.pdf | Govt of India, public |
| 19 | BEL, S-band DWR Mk II (ISRO technology, 9 m dish, dual-pol, 500 km) | https://bel-india.in/?p=2968 | vendor page |
| 20 | ITU/WMO training: E. Allaix, *Weather radar* (S 2700–2900, C 5350–5725; PRF 250–1200, PW 0.5–2.5 µs) | https://itu.int/en/ITU-R/seminars/Global-ITU-WMO/Documents/Training-Workshop-singapore/Presentations/Day-2/Eric%20Allaix_Weather%20Radar.pdf | © ITU/WMO, free |
| 21 | dump1090 (Mode A/C and Gillham reference; not used, not copied) | https://github.com/antirez/dump1090 · https://github.com/flightaware/dump1090 | BSD-3 (antirez) / GPL-2+ (FlightAware fork): oracle only |
| 22 | pyModeS (Mode S altitude / identity cross-check) | https://github.com/junzis/pyModeS | GPL-3: test oracle only |
| 23 | AAI AIP India (aerodrome VHF/UHF frequencies, for channel validation) | https://aim-india.aai.aero | Govt of India, public |
