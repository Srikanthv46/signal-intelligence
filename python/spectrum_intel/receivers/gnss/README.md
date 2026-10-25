# gnss — GPS L1 C/A satellites, GNSS jamming and spoofing indicators (research)

**Validation: research.** Built from IS-GPS-200 and the acquisition / indicator design in
`research/signal_gaps/02_gnss.md`. Checked against the ICD self-test octals (all 32 PRNs) and synthetic skies only.

**What it does.** When a block covers 1575.42 MHz (±1.1 MHz, fs ≥ 2.1 MSps, ≥ 16 ms), the receiver:
1. Crops the block FFT to ±1.024 MHz around L1 / (1 + receiver ppm), giving 2.048 MSps.
2. **Gate**, at 1 sample per chip, over all 32 PRNs:
   - Doppler bins of 125 Hz over ±8.2 kHz when the receiver ppm is known, ±25 kHz without it;
   - 4 ms coherent folds × 4 non-coherent, with one FFT per fold shared by every PRN;
   - a χ²(2M) threshold at P_fa 10⁻³ per PRN per block.
3. **Refine** each hit at 2 samples per chip over the whole block:
   - code phase (chips), Doppler (Hz), C/N0 (dB-Hz, with parabolic straddle correction plus a 1.8 dB implementation
     loss measured on synthetic skies);
   - the second correlation peak.
   - Hits ≥ 18 dB weaker than a PRN at its Doppler (mod 1 kHz) are dropped as Gold-code cross-correlation.
4. **Indicators:**
   - **Interference:** in-band J/N, i.e. L1 ±1.023 MHz against the capture's floor outside 1559–1610 MHz (needs a
     wide capture). Also narrow lines near L1, excluding the receiver's own lines.
   - **Spoofing:**
     - S1: C/N0 above 55 dB-Hz;
     - S2: ≥ 4 PRNs above 45 dB-Hz within 1.5 dB of each other;
     - S3: J/N ≥ 1.5 dB with ≥ 4 strong PRNs;
     - S4: a second correlation peak (≥ −15 dB, above threshold, outside ±1.5 chips / ±250 Hz) in ≥ 2 PRNs.
   - **Swept jammer ("personal privacy device")**, on captures ≥ 30 MSps: a 32-point FFT ridge over the first 400 µs
     around L1. The sweep period comes from the ridge's autocorrelation (0.5–50 µs, shortest strong peak), and the fold
     on that period gives span, centre, shape (sawtooth up / down, triangle) and whether the sweep crosses L1. Valid
     when the fold explains ≥ 60 % of the ridge variance and spans ≥ 1 MHz. Synthetic: 8–27 µs sweeps of 10–30 MHz
     found, no false alarm in 20 noise blocks.
   - **Jamming:** a swept jammer, or J/N ≥ 6 dB / narrow lines in band without ≥ 4 strong satellites.

**Labels** (status `partial`; the navigation message is not decoded):

| Label | When |
|---|---|
| `GNSS: GPS L1 C/A (n satellites)` | ≥ 2 PRNs, or 1 PRN at ≥ 40 dB-Hz |
| `GNSS spoofing-like (GPS L1 C/A, n satellites)` | S3, or ≥ 2 of S1 / S2 / S4 with ≥ 4 PRNs |
| `GNSS jamming-like (GPS L1, swept PPD)` | a swept jammer (confidence 0.85) |
| `GNSS jamming-like (GPS L1)` | other jamming indicators |

No satellites and no indicators produces no signal. The receiver never guesses from silence, for example an antenna
indoors.

**Fields:**
- `satellites[]`: `prn`, `cn0_dbhz`, `doppler_hz`, `code_phase_chips`, `peak_over_noise`, `threshold`,
  `second_peak_db`, `multipeak`.
- `receiver_offset_hz_estimate`: present when no receiver ppm was given. The Dopplers are then apparent.
- `search`, `interference`, `indicators`.

**Accuracy (synthetic):**
- Code phase within 0.5 chip, Doppler within 60 Hz, C/N0 within 1.5 dB at 40–48 dB-Hz.
- Detects about 38 dB-Hz and above in one 34 ms block.
- No false satellite on noise.

**Limitations:**
- GPS L1 C/A only. BeiDou B1I, GLONASS L1OF and NavIC L5 / S are in the research file (codes verified) but not yet
  built, so the single-constellation spoofing check (S7) is missing.
- Swept jammers: sweeps of ~1 µs (one device in the published survey) and periods that ripple strongly within
  400 µs are not folded; captures below 30 MSps are not searched.
- The Pluto needs an external bias-tee to power an active GNSS antenna.
