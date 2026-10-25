# analog_vtx — analog FPV video (FM PAL / NTSC)

**Detects** FM-modulated composite video transmitters (FPV cameras) anywhere in a wideband capture, far
below the FM threshold (Pd >= 0.9 from about -12.5 dB in-band SNR on synthetic VTXs).
**Decodes** the video to a picture; there is no CRC for analog video, so line and field timing locks are the
verification.

## How
1. Time-frequency power at about 1 µs × 1 MHz. Each frequency column is folded on the PAL (15625 Hz) and NTSC
   (15734.26 Hz) line periods over a ±500 ppm grid. The sync tip parks the carrier at one frequency for 4.7 µs
   every line, which gives a line-locked z-score. A score of 10 or more is an emitter.
2. Line rate is refined (fold sharpness + drift fit), along with the sync-tip frequency, the occupied band, the
   field timing (lines per field from the vertical sync) and a match to the standard FPV channel tables
   (A/B/E/F/R/L/U/O/H, 1.2 and 2.4 GHz).
3. An FM discriminator on the occupied band gives deviation, polarity and the audio subcarrier.
4. `video_frame()`: the discriminator output is sampled on the measured line clock and the active lines are
   stacked into a picture.

- **Colour standard** (`detector.chroma_standard`; **research**: `research/drone_transmitters/01_analog_video.md`,
  synthetic validation only, tagged `research_fields` in reports):
  - **Gate:** per line, the burst gate (5.7–7.5 µs after the sync leading edge), with a straight line removed (the
    porch level and the tail of the sync edge: in a 1.8 µs gate their leakage swamps every frequency; the real PAL
    VTX `LIVE_ANALOG_VTX_PAL_5740MHz` read −0.7 dB / colour false before, 6.7 dB / PAL after), is demodulated at each
    subcarrier of the measured line family, a fixed multiple of the measured line rate:
    - 525 lines: NTSC 227.5·fH, PAL-M 227.25·fH;
    - 625 lines: PAL 283.75·fH + 25 Hz, PAL-N 229.25·fH + 25 Hz.
  - **Coherence:** the line-to-line phasors, constant (NTSC) or ±45° alternating (PAL family), give a coherence
    (FFT over lines, ±0.06 cycles/line for the line-rate error) and a burst swing.
  - **PAL vs PAL-N:** decided by the burst energy at 4.43 vs 3.58 MHz. The phasor test only sees fsc mod fH.
  - **Energy reference:** the median over a 0.25 MHz grid from 1.5 MHz to 0.45 × the slice rate, ≥ 1 MHz from every
    candidate subcarrier (the old fixed set included 5.4 MHz, above Nyquist at ~10 MSps, and 4.0 MHz, inside the
    burst's main lobe). Gate: coherence ≥ 0.35, energy ≥ 6 dB, swing (PAL family > 25°, NTSC < 20°).
  - **Edges:** both H-sync edges are tried (the sync time can lock on the trailing edge at low SNR).
  - **Outputs:** `color` (true / false; null below 10 dB SNR: not measurable), `color_standard`, `chroma_fsc_hz`
    (the camera crystal: 227.25 × the line rate, ppm and all), `burst_swing_deg`, `chroma_coherence`,
    `chroma_candidates`; label `pal_m` / `pal_n` ("Analog FPV video (PAL-M)").
  - **Synthetic results** (`python/tests/test_analog_color_standard.py`): all 4 standards and monochrome NTSC / PAL correct;
    colour standard correct down to 7 dB in-band SNR (PAL-N 4 dB); 2.3 ms per gate in Python.
  - **Real transmitters:** NTSC (Analog_20260804, coherence 0.958, 9.0 dB, swing 0.3°) and PAL (user-confirmed colour
    camera, LIVE_ANALOG_VTX_PAL_5740MHz: coherence 0.875, 6.7 dB, swing 45.5°, fsc 4,433,716 Hz).
- **SECAM** (`detector.secam_test`; **research** (ITU-R BT.470), synthetic only, 2026-10-21): no burst. The back porch
  carries the FM subcarrier at its rest frequency, Db 4.250 MHz and Dr 4.40625 MHz on alternate lines.
  - Run only on 625-line video where no PAL-family burst passes.
  - Per line, the energies at the two rest frequencies over a 4.2 µs gate (5.8–10 µs) give
    r = (e_Db − e_Dr)/(e_Db + e_Dr). SECAM alternates it line by line (|mean((−1)^n r)| ~0.6–0.7, t ≥ 6); PAL and
    noise give ~0.
  - The rest frequencies are then measured per line parity on a 5 kHz grid (`secam_db_hz`, `secam_dr_hz`: the
    camera's subcarrier generators).
  - Label `secam` ("Analog FPV video (SECAM)").
  - Synthetic results: detected at 20 and 10 dB in-band SNR, rest frequencies within 5 kHz; PAL / PAL-N never read
    as SECAM. ~10 ms in Python.
- **Audio** (`detector.audio_features`; **research**, synthetic only, 2026-10-21, PLAN A5): when the audio FM subcarrier
  (5.5 / 6.0 / 6.5 / 7.0 MHz) stands ≥ 15 dB above its surroundings, it is sliced ±100 kHz and FM-discriminated.
  Outputs:
  - `audio_subcarrier_hz`: the exact subcarrier;
  - `audio_dev_rms_khz`: the audio level;
  - `audio_peak_hz` / `audio_peak_db`: the dominant tone between 50 Hz and 8 kHz, to ~2 Hz from a zero-padded Hann
    periodogram. On a drone microphone this is typically the propeller blade-pass frequency.
  - `audio_flatness`: tone ~0, noise ~1.

  Synthetic results: tones of 180, 440 and 1250 Hz within 5 Hz at 20 and 12 dB.
- **Channel table:** + `LB` (Lowband, 5333–5613 MHz, 40 MHz step; research), lowest priority; `channel_note` when matched.

## API (depends on `spectrum_intel.core` only)
```python
from spectrum_intel.core.emission import load_capture
from spectrum_intel.receivers.analog_vtx import detect, video_frame
x, fs, fc, meta = load_capture("cap.sigmf-meta")
for e in detect(x, fs, fc):           # Emission(family="analog_vtx", label="pal"|"ntsc"|"analog_fm_video")
    print(e.label, e.center_mhz, e.features["line_rate_hz"], e.features.get("channel"))
    img = video_frame(x, fs, fc, e)   # uint8 [lines, 720]
```
Streaming video (`video.py`), one field per block. Each field is locked on the vertical sync, with black/white
levels taken from the back porch. Built for weak VTXs (below the FM threshold, where FM "clicks" turn into
black/white dots):
- **discriminator:** the lag-1 product y[n]·y*[n−1] is low-passed (7-tap, 6 MHz) before the angle, so each phase
  step counts by its power and the low-amplitude samples where noise spins the phase (the clicks) count little;
- **line clock (an H-PLL over the block):** the sync phase from the fold at the detector's rate, then every
  line's sync leading edge, placed by the 4.7 µs sync box (the sync is the lowest level, so a bright picture
  cannot pull it) and timed by a short step filter within ±0.6 µs; a robust straight-line fit over all ~500
  lines gives the line period and phase (meta `sync_jitter_us`). A weak block is weighted against the
  detector's rate (±2 ppm prior) unless they disagree by more than 4σ. The old per-block fold-sharpness refit
  (±40 ppm) was driven by the picture: 14–39 ppm off even at 14 dB SNR (a sheared picture), −31.9 ppm for
  +2.1 on the weak 1.2 GHz VTX;
- **noise:** the robust spread of the flat sync tip over white − black (meta `noise`); above 0.09 (about 5 dB
  in-band SNR) FM clicks dominate and the field is 3×3-median filtered (meta `denoise`), white is taken after it;
  above 0.25 colour is killed (`color_killer: "noise"`): on synthetic scenes grey is then closer to the truth.

When the VTX sends a colour burst, the field also comes in colour (`rgb`):
- the discriminator output, low-passed to 6 MHz, is sampled on 1024 points per line and mixed with the
  line-locked subcarrier (NTSC 227.5 × fl, PAL 283.7516 × fl); a 95-tap 1 MHz filter gives the chroma;
- the phase reference is each line's burst (the mean over the burst gate) averaged over the ±3 nearest lines of
  the same parity, like a TV's burst PLL; a weak VTX's burst is noisy on every single line. PAL: the reference is
  the mean of two lines' (the ±45° swing cancels), the swing gives the V switch, and a two-line delay line
  averages U/V;
- the burst gate is tried at 0, +0.2 and +0.4 µs after the standard's (cameras put the burst up to ~0.4 µs late
  relative to the 50 % sync edge; all three real VTXs pick +0.4) and the strongest field burst is kept
  (meta `burst_gate_us`); the off-subcarrier reference is floored at its noise expectation, so a reference small by
  chance cannot make a gate look like a burst;
- the colour picture comes from the plain discriminator (the weighted one's short filter attenuates the subcarrier);
- levels come from the sync depth, as in a TV (NTSC 100 IRE white, 7.5 IRE setup; PAL 700 / 300 mV). Saturation
  comes from the burst (ACC: 20 IRE / 150 mV). A picture whiter than 1.3× the white the sync depth allows means a
  compressed sync (the weak 1.2 GHz VTX: 3.9× the depth instead of 2.5×; a washed-out picture): the depth is then
  taken from the picture's white (meta `sync_compressed`);
- colour killer: `rgb` is None unless the field's burst, averaged coherently over its lines (lag-2 products, so
  slow phase drift and the PAL swing do not cancel it), is ≥ 4× the same measure 1 MHz below the subcarrier
  (the audio subcarrier is above). Without a burst it stays ≤ 1.64 (60 synthetic noise draws); the user's live VTX
  gives 28 (17 dB) and 15 (weak, 10.7 dB).
  A per-line test (the first version) turned that VTX grey. `color=False` turns it off; `image` (luma) is unchanged either way.
```python
from spectrum_intel.receivers.analog_vtx import FieldDemodulator
demod = FieldDemodulator.from_fields(report_signal["decoded"]["fields"])   # or .from_emission(e)
frame = demod.process(iq_block, fs, fc)       # VideoFrame(image uint8 [240|288, 720], meta, rgb) or None
frame.rgb                                     # uint8 [240|288, 720, 3], or None (no colour burst)
```
signal-vault's Signal Intelligence tab uses this for its live **▶ View feed** window.

`synth.py` generates synthetic VTXs (standard, deviation, clock error, audio) for testing.

## Pipeline plugin (`plugin.py`)
- Runs when the capture overlaps 1.06–1.38, 2.33–2.53 or 4.9–6.0 GHz and has at least 20 MHz of band
  (min(sample rate, analog bandwidth); `schema/requirements.json`): a narrower block cannot hold the FM video (the
  response carries SAMPLE_RATE_TOO_LOW), and a long narrow dwell would cost minutes in the line fold.
- Rejects a capture-chain artifact (`detector._capture_locked`): a "line clock" within 2 ppm of fs / N, with its
  "sync tip" within 100 kHz of the LO and no measured field sync — LO leakage modulated by a periodic process of the
  capture (live Pluto at 8 MSps: exactly fs / 512 = 15625 Hz, sync tip 31–43 kHz from the LO, reported as a PAL
  picture at −13 to −18 dB). A real VTX meets all three only if its crystal error cancels the fs / N offset.
- Produces one signal per VTX with status `success`, and saves the picture as `<name>_vtx<MHz>_video.png`.
- Owns its band, and via its `arbitrate` hook it absorbs weak receiver detections inside its FM spectrum. Those
  are the VTX's own sidebands; the margin is max(6 MHz, 0.2 × bandwidth).

## Validation
`tests/receivers/test_det_analog_vtx.py`:
- **Pd vs SNR** on noise, Wi-Fi and DJI backgrounds; PAL/NTSC classification 152/152.
- **Accuracy:** line rate 0.9 ppm rms; centre frequency 0.28 MHz rms.
- **False alarms:** none on 6 real captures with no VTX.
- **Real capture:** the NTSC VTX at 5015.68 MHz in `testdata/Analog_20260804-161815`.
- **Live null:** `testdata/live/LIVE_ADSB_1090MHz` (the capture-locked line above): no VTX.
- **Colour** (`tests/receivers/test_analog_color.py`): synthetic 75 % colour bars, NTSC and PAL, −60 / +80 ppm,
  clean and at 12 dB in-band SNR, with an audio subcarrier (PAL: 5.5 MHz, 1.07 MHz above the colour subcarrier):
  every bar within 24/255 per RGB channel (typically ≤ 15). No burst → grey (colour killer). The real NTSC VTX:
  burst 23× the reference, a natural-colour picture (a hand: skin tones, a red-orange part, blue stripes), luma
  unchanged. The user's NTSC VTX live at 5945 MHz (`testdata/live/LIVE_ANALOG_VTX_5945MHz`, 17 dB, and `_weak`,
  10.7 dB with FM clicks): both in colour. Cost: about +90 ms per field in Python (the luma field takes about 300 ms).
- **Low SNR** (`tests/receivers/test_analog_video.py`): synthetic NTSC VTX, the demodulator 25 ppm off the true
  rate: line rate within 0.5 ppm clean / 14 dB, 1 ppm at 6 dB, 3 ppm at 3 dB; noise estimate monotonic; median on at
  3 dB, off at 14 dB; colour killed at 0 dB. Field PSNR against the clean picture (3 scenes, shift-tolerant), old →
  new: 0 dB 13.8 → 14.3, 2 dB 14.1 → 16.4, 4 dB 15.6 → 17.6, 6 dB 18.7 → 20.0, 10 dB 19.2 → 25.4, 14 dB
  20.9 → 28.1 dB. Real weak VTX `testdata/live/LIVE_ANALOG_VTX_1200MHz_weak` (1199.84 MHz, sync-tip noise 0.29):
  line rate +4.2 ppm (detector +2.1; old −31.9), straight picture edges, median filtered, grey.
- **Speed:** a field in 86 ms (strong live VTX) / 136 ms (weak 1.2 GHz) with the pipeline's channelizer, was 175 /
  207 ms (the ±40 ppm fold search is gone).
- **Known misses:** peak deviation error up to 3.7 MHz (limit 1.5); audio subcarrier on a busy Wi-Fi background
  9/15. Below about 2 dB in-band SNR one field is still grainy (FM clicks); a static scene would gain from averaging
  fields across blocks (not done yet).

## 3.3 GHz analog band (2026-10-23, PLAN G9; research)
The plugin searches `core.bands_india.ANALOG_VTX_RANGES` (adds 3150–3490 MHz) and the channel table names the two
8-channel 3.3 GHz groups (3G3A 3330–3470, 3G3B 3170–3310 MHz, 20 MHz steps; one vendor's VRX specification). The wide
4.9–6.0 GHz 80 / 96-channel groups are only published as images: those VTXs are detected but get no channel name.
