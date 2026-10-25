# jammer — swept-carrier jammers

**Detects** barrage jammers of the common kind: one carrier (a VCO) swept back and forth across a band by a
periodic tuning voltage, transmitting continuously. Anywhere in the capture, with fs ≥ 20 MHz.
**Measures** (there is no payload to decode):
- the sweep rate and period (about ±10 ppm on fast sweeps, ±100 ppm when the block holds only ~40 sweeps);
- the frequencies the carrier reaches inside the capture, and whether it continues below / above it;
- the shape (up-down = triangle or sine, sawtooth, or undetermined when the carrier leaves the band);
- the time in band, its dominance and the mean trajectory (time in the period → frequency).

## How
1. **Occupancy.** Hann STFT power with fs/64 bins and hop 32, summed over frame pairs: 0.96 MHz × 1.04 µs rows
   at 61.44 MSps. The outer 1/16 of the band each side (anti-alias roll-off) is dropped.
2. **Trigger.** A swept carrier lights every bin it crosses once or twice per sweep, so each bin's occupancy is
   periodic. The cyclic spectrum of the per-bin *log* occupancy, whitened per bin (bursty bins do not dominate)
   and summed over bins, is searched for lines ≥ 10 dB over the local median, between max(300 Hz, 8 sweeps per
   block) and the row rate / 8 (120 kHz). Up to 3 lines (not harmonics of each other) are refined on their
   harmonics and tried in turn.
3. **Period.** The true period is a multiple d ≤ 8 of a line's period. A fold on a fraction of the true period
   averages shifted copies of the pattern (smoother); a fold on a multiple repeats it. The period is the shortest
   whose phase structure (variance of fold / per-bin mean − 1, linear so a k-times diluted pattern scores 1/k) is
   ≥ 80 % of the best. It is refused when a fold on 2–5 × that period shows clearly more structure (a pulsed
   signal folded on a fraction of its pulse period). The rate is corrected from the drift of 2–4 block segments'
   folds against the whole-block fold.
4. **Confirmation** on the fold (≤ 128 phases), per phase the peak bin over that bin's floor (10th percentile over
   phases):
   - the carrier is present (peak ≥ max(3 dB, the strong half's median − 10 dB)) for ≥ 30 % of the period;
   - it spans ≥ 10 MHz (FSK and FM video are narrower);
   - it moves continuously: on the trajectory after a 5-phase circular median (none below 40 phases per period, so
     fast sweeps keep their turns), at most one flyback and at most two band exits per period (two when clipping
     images split the out-of-band gap). A jump of more than max(4 bins, 0.3 × span) between consecutive phases is a
     flyback only when some bin between the two peaks is dark in both phases (more than 10 dB under the weaker
     peak): a real sweep's fast transit, smeared in the fold by its ~2 µs sweep-to-sweep jitter, keeps those bins
     within 1–6 dB of the peaks (bench jammer), while a hop or a sawtooth's retrace leaves them dark (a strong
     hopper's switching splatter stays ~30 dB down);
   - it dominates the span it sweeps: the 75th percentile over phases of the peak over the span's median ≥ 4 dB
     (ADC clipping spreads a strong carrier's power across the band, so the median phase can be lower);
   - analog video line rates (15625 / 15734.26 Hz and harmonics) are refused.
5. **Band.** From the smoothed trajectory outward over bins the fold shows lit at some phase (each sweep's extremes
   wander ~2 MHz around the mean), at most 3 bins (~3 MHz): ADC clipping spreads a saturating jammer far wider, and
   the band must not take in a neighbour. It is extended to the capture edge on the side the carrier leaves through. The signal owns
   that band: the energy detector's pieces of the sweep inside it are absorbed.

About 0.1–0.3 s per 2M-sample block in Python (one 64-point STFT, one cyclic FFT, ~15 sparse folds).

## Labels
- **"Jammer (swept carrier)"**, status `success`, confidence 0.9 (0.75 when the sweep is under 20 MHz or the
  dominance under 6 dB). Technology `swept_jammer`, family `jamming`.
- **"Swept carrier (FMCW-like)"**, confidence 0.6, in the radio-altimeter band (4200–4400 MHz), where a continuous
  sweep is the band's own FMCW radar. Technology `swept_carrier`.

## Validation
- `tests/receivers/test_det_jammer.py` (pytest): triangle, sawtooth, sine and noise-FM sweeps of 1.2–50 kHz over
  20–50 MHz, jittered and clipped sweeps leaving the band, a sweep over Wi-Fi traffic, down to +4 dB per MHz;
  nothing on noise, CW, pulsed and pulsed-LFM radar, TDD OFDM (1 and 0.1 ms), periodic and random hoppers, Wi-Fi
  traffic, FM video or a 4 MHz sweep; nothing (rather than a wrong sweep) at +2 dB per MHz.
- Bench: `testdata/live/LIVE_JAMMER_SWEEP_2460MHz` (a jammer, confirmed by the user, 2026-10-20): 10000.25 Hz
  (long-baseline reference 10000.09–10000.21 Hz), seen 2435.3–2474.6 MHz and continuing below the capture,
  up-down, in band 88 % of the period; it replaces nine "unknown" energy boxes.
- The whole signal-vault session (10 blocks 0.32 s apart, `captures/2026-10-19/pluto/2460MHz_JAMMING_...`; tested
  when present): 10/10 found, 10000.13–10000.27 Hz, top 2469–2478 MHz, up-down, no flyback. One block
  (`LIVE_JAMMER_SWEEP_2460MHz_jitter`) had been refused as a hopper (3 "jumps" between dwell frequencies) before
  jumps were tested for bridging, and another had a one-phase outlier at 2484 MHz before the trajectory was smoothed.

## Limits
- One sweep per block (the strongest that confirms).
- Sensitivity about +4 dB per MHz over the noise; a sweep much weaker than the random traffic sharing its bins
  (e.g. 12 dB under 50 %-duty Wi-Fi) is not seen in the mean fold.
- Noise jammers are handled by the barrage decoder below. Spot jammers (fixed carrier) and very fast sweeps
  (> 120 kHz, period < 8 µs; GNSS PPDs have their own detector in `gnss`) are not detected here;
  frequency-hopping jammers are refused like any hopper.
- The jammer identity is inferred from the waveform (no payload exists); FMCW radars outside the altimeter band
  would also be called jammers.

## Noise jammers: barrage and noise-FM (`barrage.py`, plugin `jammer_barrage`; research, signal gaps S12)

**Validation: research** (`research/signal_gaps/05_aviation_radar_jammers.md` section 5.4; synthetic only).

A per-signal decoder on the energy detector's boxes that no receiver explained: ≥ 1 MHz wide, continuous over the
block, ≥ 6 dB. It tests what a noise jammer is and no service is:
1. **Plateau** (block FFT, inner 90 % in 64 sub-bands): ≥ 90 % within ±3 dB of the median, the strips beyond each
   edge ≥ 6 dB lower, no narrow line.
2. **Statistics** at 1.25 × the width: envelope k = E|x|⁴ / E²|x|² (Gaussian 2, constant envelope 1), I / Q kurtosis.
3. **No structure:** the autocorrelation of the Hann-tapered spectrum stays < 0.025 over max(1 µs, 8 / BW) to
   1000 µs (an OFDM cyclic prefix gives its CP fraction: LTE at 0 dB ≈ 0.035, Wi-Fi-like ≈ 0.09); no cyclic line in
   |x|², x² or (constant envelope) the instantaneous frequency (each bin against its 512-bin block, P_fa ≈ 1e-4);
   power per 1 ms within ±1 dB (no TDD gaps).
4. **Labels:** "Jammer (barrage noise)-like" (Gaussian, flat, no structure, ≥ 10 MHz or edges on an allocation),
   confidence 0.6, 0.75 with `covers_band` (an allocation in `core/bands_india` whose edges are within 2 MHz);
   "Jammer (noise-FM)-like" (constant envelope, instantaneous-frequency spread ≥ 15 % of the width with a
   Gaussian-like kurtosis 2.2–4.5), 0.7. Technology `noise_jammer`, status `partial`.

**Synthetic results:** barrage plateaus of 20–30 MHz at 8–15 dB and a noise-FM jammer are labelled; LTE (0–5 dB),
NR 30 kHz, Wi-Fi-like OFDM (0 dB), TDD bursts, a swept carrier with FM noise, analog FM video and receiver noise are
not. About 0.1–0.2 s per candidate box in Python.

**Limits:** a 5 MHz Gaussian plateau is not called a jammer unless its edges match an allocation (one WCDMA
carrier looks the same); the energy detector may split a noise-FM spectrum into pieces, and only the widest piece is
labelled; ADC clipping and the Pluto's USB self-noise are excluded only through the stationarity and cyclic tests;
no real noise jammer has been captured.
