# radalt (C++) — FMCW radio altimeters, 4200–4400 MHz (research)

Port of `python/spectrum_intel/receivers/radalt` (detector.py, plugin.py). Same method; see the Python README for
the reasoning and the validation (synthetic sweeps only).

## What it does
- **Detects** the sweep's straight crossings through the capture (fs ≥ 20 MSps): 128-point Hann ridge every 64
  samples, gap-tolerant tracks (≤ 4 bins per elapsed frame, gaps ≤ 8 frames, ≥ half the spanned frames on),
  R² ≥ 0.99 linear fits, a V split at its extreme, ≥ 15 MHz covered.
- **Measures** rate, shape (triangle / sawtooth / single crossing), turnarounds, sweep bandwidth and period.
- Plugin: "Radio altimeter (FMCW)-like (...)", technology `radio_altimeter`, owns its band.

## C++ API (`namespace si::rx::radalt`, depends on core/ only)
```cpp
#include "receivers/radalt/detector.hpp"
std::vector<si::Emission> em = detect(iq, fs, fc);   // at most one emission
```

## Implementation notes
- The ridge rows run in blocks of 1024 frames on `parallel_for` (no per-frame allocation): ~45 ms per 2^21 block
  at 61.44 MSps on 2 threads; numpy order of operations elsewhere (pairwise sums, `py_round`).
- vectors in `tests/vectors/signal_gaps/radalt*` (`test_signal_gaps.cpp`).
