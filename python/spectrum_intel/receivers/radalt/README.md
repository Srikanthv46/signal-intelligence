# radalt — FMCW radio altimeters, 4200–4400 MHz (research)

**Validation: research.** Built from ITU-R M.2059 (via `research/signal_gaps/05_aviation_radar_jammers.md`
section 1) and checked on synthetic sweeps only.

**Detects** the height-finding radar under aircraft and helicopters. An FMCW altimeter sweeps 104–177 MHz, wider
than any capture, so a block sees the sweep as straight frequency ramps ("crossings") through the band. Gate: the
capture overlaps 4200–4400 MHz and fs ≥ 20 MSps.

**Measures** (`features` / decoded fields):

| Field | Meaning |
|---|---|
| `n_crossings`, `r2_min` | straight crossings found, worst linear-fit R² |
| `rate_mhz_per_ms` | sweep rate (median of the crossing slopes) |
| `shape` | `triangle`, `sawtooth up` / `sawtooth down`, or `single crossing` (too slow for a second one in the block) |
| `f_top_mhz`, `f_bottom_mhz`, `sweep_bw_mhz` | turnarounds from up/down crossing pairs, f_c ± r Δt / 2, and their difference |
| `period_ms` | 2 B / r for a triangle with both turnarounds; otherwise the spacing of same-direction crossings, in whole periods |
| `seen_mhz` | the part of the sweep inside the capture |

**How:**
1. Ridge: 128-point Hann FFT every 64 samples; per frame the peak bin (parabolic) when it stands ≥ 10 dB over the
   frame median.
2. Tracks: a ridge frame joins the first open track it continues (≤ 4 bins per elapsed frame, gaps ≤ 8 frames), so
   fades and noise peaks do not cut a sweep. A track needs ≥ 5 frames, on in ≥ half the frames it spans (chance
   noise chains are sparse), and a linear fit with R² ≥ 0.99. A track that fails the fit is split once at its
   frequency extreme: that is a turnaround inside the capture. Same-direction pieces whose centre crossings lie
   within 20 µs merge. A crossing covers ≥ 15 MHz.
3. Sweep parameters as in the table; the rate must be 2–1500 MHz/ms.

**Confidence:** 0.7; 0.8 when the bandwidth is 95–190 MHz or the period matches an M.2059 type (20, 6.67, 7 or
1 ms, ± 5 %). The label is always "Radio altimeter (FMCW)-like (...)": the band and the sweep are evidence, not a
decoded identity. **Status** `partial`: an altimeter carries no payload. The receiver owns the band, so the sweep's
pieces do not become `unknown` boxes. The jammer receiver looks for µs-period sweeps and, in this band, calls one
"Swept carrier (FMCW-like)", never a jammer.

**Accuracy (synthetic, 61.44 MSps, 34 ms):**
- triangle 104 MHz / 50 Hz, 133 MHz / 300 Hz, 150 MHz / 143 Hz: bandwidth within 0.5 MHz, period within 10 µs, at
  amplitudes down to −10 dB per-sample SNR (a turnaround inside the capture included);
- sawtooth 176.8 MHz / 1 kHz: period 1.000 ms;
- no detection on 60 noise blocks or on a 20 MHz wideband block with a CW carrier.

**Limitations:**
- No real altimeter capture has been tested yet.
- Below about −10 dB per-sample SNR the 128-point ridge breaks up; a chirp-matched search would go further.
- A pulsed (non-FMCW) altimeter is not recognised here.
- Several altimeters in one block (multiple aircraft) are reported as one signal with the median rate.
