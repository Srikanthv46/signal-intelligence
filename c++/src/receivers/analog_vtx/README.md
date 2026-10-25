# analog_vtx (C++) — analog FPV video (FM PAL / NTSC)

Port of `python/spectrum_intel/receivers/analog_vtx` (detector.py, video.py, plugin.py). Same method; see the
Python README for the physics.

## What it does
- **Detects** FM composite-video transmitters anywhere in the capture, far below the FM threshold: a ~1 µs × ~1 MHz
  STFT is folded on the PAL / NTSC line period over ±500 ppm; the sync tip gives a line-locked z-score (≥ 10 →
  emitter).
- **Characterises** each VTX: line rate (fold sharpness + drift regression), sync-tip frequency, line-locked
  spectrum (centre, band, SNR, peak white → deviation), field timing (lines per field from the vertical sync),
  FM-discriminator features (levels, colour burst, audio subcarrier), FPV channel table.
- **Decodes** the video to a picture (`Signal::video_frame`, uint8 values 0–255 stored as `Image` floats).
- Plugin: runs on 1.06–1.38 / 2.33–2.53 / 4.9–6.0 GHz, owns its band, and `arbitrate` absorbs weak rc_link /
  digital_vtx detections inside the VTX's FM spectrum (margin max(6 MHz, 0.2 × bandwidth)).

## C++ API (`namespace si::rx::analog_vtx`, depends on core/ and io/ only)
```cpp
#include "receivers/analog_vtx/detector.hpp"
#include "receivers/analog_vtx/video.hpp"
#include "receivers/analog_vtx/png.hpp"
std::vector<si::Emission> em = detect(iq, fs, fc);            // label pal | ntsc | analog_fm_video
std::shared_ptr<si::Image> img = video_frame(iq, fs, fc, em[0]);   // [lines, 720], or nullptr
write_png("vtx.png", *img);                                   // 8-bit grayscale PNG

// streaming, one field per block (GUI / edge viewer; python FieldDemodulator): luma, plus colour
// when the VTX sends a colour burst (fr->rgb, [rows, cols, 3], empty otherwise; set_color(false) turns it off)
FieldDemodulator demod = FieldDemodulator::from_fields(report_signal["decoded"]["fields"]);  // or from_emission(e)
if (auto fr = demod.process(block, fs, fc)) show(fr->image, fr->rows, fr->cols);   // fr->meta as Python
```
Building blocks (`tf_power`, `fold`, `zscore`, `Scanner`, `refine_timing`, `field_timing`,
`line_locked_spectrum`, `spectrum_features`, `fm_discriminate`, `demod_features`, `match_channel`) are public for
tests.

## Implementation notes
- dtypes as numpy: STFT power float32 (np.fft on complex64 → `fft::np_fft`), folds float32-accumulated in time
  order (scipy sparse product) then float64; NEP 50 promotions reproduced (e.g. float32 arrays divided by numpy
  int64/float64 scalars become float64).
- Round-trip of numpy vs Python floats in the features: `round(np.float64, n)` → `np::round`, Python floats →
  `py_round`.
- numpy temporary elision: `y[1:] * np.conj(y[:-1])` on ≥ 256 KiB is computed by numpy as `conj(prev) * cur`;
  the fma complex multiply is not symmetric, so the discriminators mirror the operand order.
- Speed: every (standard, ppm) scan hypothesis is folded once in parallel and its z-score map kept, so the
  repeated scans with growing column masks only redo the argmax; grid line-rate search, drift blocks, the
  FieldDemodulator rate fit and rasters run on `si::pmap` / `parallel_for` (ordered results). `np.mod(x, 1.0)`
  is computed as `x - floor(x)` (bit-identical, much faster than fmod).
- Pictures: `png.hpp` writes stored-deflate PNG (no zlib). The bytes differ from PIL's compressed file; the
  decoded pixels are identical (checked: PIL reads it back equal to Python's `video_frame`).
- FieldDemodulator: the 31-tap low-pass is a direct float64 convolution (scipy: oaconvolve float32 × float64);
  the power-weighted discriminator's 7-tap FIR runs in complex64 as Python's; the line clock (per-line sync edge,
  robust `np::polyfit` line fit, detector-rate prior), sync-tip noise, 3×3 median (`ndi::median_filter`, nearest),
  the three burst gates with the floored reference, compressed-sync levels and the noise colour killer follow
  `video.py`. Fields agree with Python within 1 grey level on < 0.1 % of pixels, identical field timing and meta.

## Colour standard (research, 2026-10-21)
`chroma_standard()` ports python `chroma_standard` (NTSC / PAL-M / PAL / PAL-N from line-locked burst phasors,
both H-sync edges tried, colour null below 10 dB; 2026-10-22: per-line detrend of the burst gate and the in-band
reference grid, as Python; vectors `gen_chroma.py` and `gen_analog_vtx.py` regenerated, all identical). The `LB` Lowband table and the plugin fields (color,
color_standard, chroma_*, channel_note) and the PAL-M / PAL-N labels are ported too.
`test_rx_analog_vtx.cpp` "colour standard on synthetic composite" checks it against
`tests/vectors/analog_vtx/gen_chroma.py` (synthetic only, no capture): 5 cases, 60 assertions identical.
The capture-based analog vectors (gen_analog_vtx.py on testdata) were **not** regenerated: the user asked not to use
local captures. Regenerate them (and run parity) once allowed.

## Parity / tests
- `tools/parity.py --receivers analog_vtx,voice_fm`: 16/16 inputs identical; the NTSC VTX in
  `Analog_20260804-161815` matches in every field; no false VTX on the D_* drone captures or sweep snapshots.
- `tests/unit/test_rx_analog_vtx.cpp` (vectors: `tests/vectors/analog_vtx/gen_analog_vtx.py`): STFT power,
  normalised power and best fold bit-identical; scan profile (all 216 hypotheses) to 1e-9; every Emission
  feature (incl. ones not in the report: burst_db, vbi_us, levels) within the last printed decimal; the picture
  pixel-identical; FieldDemodulator on two blocks.
- Timing (receiver_s, Analog capture, loaded machine): C++ 1.3–1.6 s vs Python 2.8 s; no-VTX capture 0.2–0.3 s
  vs 0.5 s. What remains is pocketfft on the slice lengths numpy uses (some have large prime factors), which the
  bit-identity requires.
