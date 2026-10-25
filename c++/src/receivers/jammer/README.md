# jammer (C++) — swept-carrier jammers

Port of `python/spectrum_intel/receivers/jammer` (detector.py, plugin.py). Same method; see the Python README for
the reasoning and the validation.

## What it does
- **Detects** one carrier swept periodically across ≥ 10 MHz, continuously (in the captured band ≥ 30 % of the
  period), anywhere in a capture with fs ≥ 20 MHz: a cyclic line in the per-bin STFT occupancy gives the rate, the
  STFT folded on the true period (the shortest multiple of the line period with ≥ 80 % of the best phase
  structure) shows one carrier moving continuously.
- Jumps between fold phases count as flybacks only when the bins between the two peaks are dark (> 10 dB under the
  weaker peak); the trajectory is a 5-phase circular median (none below 40 phases); the band extends ≤ 3 bins beyond it.
- **Measures** sweep rate / period (drift-corrected over 2–4 block segments), the frequencies reached, the side it
  leaves the band through, shape (up-down / sawtooth / undetermined), time in band, dominance and the trajectory.
- Plugin: "Jammer (swept carrier)" (technology `swept_jammer`, family `jamming`), or "Swept carrier (FMCW-like)"
  in the radio-altimeter band; owns its band, so the energy detector's pieces of the sweep are absorbed.

## C++ API (`namespace si::rx::jammer`, depends on core/ only)
```cpp
#include "receivers/jammer/detector.hpp"
std::vector<si::Emission> em = detect(iq, fs, fc);   // family "jammer", label "swept_carrier", or empty
Sweep s = measure(iq, fs, fc);                        // internals: period, fold trajectory, diagnostics (s.tried)
```

## Implementation notes
- dtypes as numpy: occupancy (complex64 STFT, `fft::fft_rows`), its log, the per-column cyclic `rfft` and the
  whitened sum (pairwise) are float32; folds accumulate float32 in frame order (scipy sparse product) and continue in
  float64; scalars are Python floats (`py_round` for every reported value).
- `std::log(float)` may differ from numpy's float32 log by an ulp; it only enters the trigger statistics.
- About 20–40 ms per 2M-sample block (one 64-point STFT, 56 cyclic FFTs, ~15 folds).

## Noise jammers (`barrage.hpp`, plugin `jammer_barrage`; research)
Port of `python/spectrum_intel/receivers/jammer/barrage.py`. A per-signal decoder on unexplained energy boxes ≥ 1 MHz,
continuous, ≥ 6 dB: flat plateau, Gaussian envelope (k = E|x|⁴ / E²|x|² ≈ 2) or constant envelope with a Gaussian-like
instantaneous frequency, no cyclic prefix (Hann-tapered autocorrelation < 0.025), no cyclic line, stationary per 1 ms →
"Jammer (barrage noise)-like" (0.6, 0.75 with `covers_band`) or "Jammer (noise-FM)-like" (0.7).
```cpp
#include "receivers/jammer/barrage.hpp"
std::optional<Json> r = barrage::analyse(channelizer, f0, f1);   // (*r)["class"]: "barrage" | "noise_fm" | null
```
- numpy dtypes and order (complex128 analysis of the complex64 slice, `std::pow(x, 4.0)` as numpy's power, the
  numpy temporary-elision operand order in the discriminator); FFT tests on the first 2^17 samples: ~60 ms per
  candidate box. Vectors: `tests/vectors/signal_gaps/barrage*`.
