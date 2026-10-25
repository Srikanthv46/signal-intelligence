# navaids — ILS glide slope, ILS localizer, VOR (research)

**Validation: research.** Built from ICAO Annex 10 Vol I (via `research/signal_gaps/05_aviation_radar_jammers.md`
section 2 and `06_vhf_below_300.md` section 9) and checked on synthetic carriers only.

**Detects** the ground stations that guide aircraft:
- **ILS glide slope:** 40 channels, 329.15–335.00 MHz.
- **ILS localizer:** 108.10–111.95 MHz, odd tenths, ICAO-paired with the glide path.
- **VOR:** 108.00–117.95 MHz on a 50 kHz grid; below 112 MHz only the even tenths.

A channel counts when a carrier line stands ≥ 15 dB over the median of ±35 kHz around it, with the receiver ppm
applied. A second line 4–32 kHz away marks a two-frequency (capture-effect) installation.

**Measures:**

| Kind | Fields |
|---|---|
| ILS (`ils_gs`, `ils_loc`) | `m90`, `m150` (tone depths), `ddm` = m90 − m150 and `side` (glide path: above / below the path; localizer: left / right of course as seen on the approach), `sdm`, `tone_scale` (the fitted 90 / 150 Hz scale), `depth_sigma`, `paired_mhz`, `two_frequency_hz` |
| VOR | `radial_deg` (bearing from the station), `am_depth` (30 Hz), `deviation_hz` (9960 Hz subcarrier FM), `subcarrier_db` |

**How:**
- **ILS:** least-squares fit of the envelope over whole 1/30 s windows. 90 and 150 Hz are both harmonics of 30 Hz,
  so they are orthogonal over that window. The fit uses a common tone scale searched over 0.97–1.03, which removes
  the leakage that fixed DFT bins have at the ±2.5 % tone tolerance.
- **VOR:**
  - the 30 Hz AM is fitted on the envelope;
  - the 9960 Hz subcarrier is cut from the envelope's FFT, and its instantaneous frequency is fitted at 30 Hz;
  - the radial is the AM phase minus the FM phase. A Doppler VOR gives the same result by design.

**Validity:**
- **ILS:** both tones ≥ 5 σ, and an SDM that fits the kind: 0.6–1.0 for the glide path, 0.25–0.65 for the localizer.
- **VOR:** subcarrier ≥ 10 dB over its surroundings, AM depth 0.15–0.45, deviation 300–700 Hz.
- Anything short of this is labelled "-like" (confidence 0.55). A plain carrier on the grid, or a channel without
  the 9960 Hz subcarrier, is not reported.

**Status** is `partial`: these stations carry no payload, and the Morse ident (1020 Hz) needs a dwell of several seconds.

**Accuracy (synthetic):**
- **ILS DDM:** within 0.005 from one 1/30 s window at 30 dB.
- **VOR radial:** within 2° from one 34 ms block at 25 dB.

**Limitations:**
- No real ILS or VOR capture has been tested yet.
- The localizer and VOR lie below 300 MHz, so they need a receiver that tunes there (Pluto with the AD9364 range).
- The DDM gives only the receiver's side of the course or path, not a position.
