# gnss (C++) — GPS L1 C/A acquisition, GNSS jamming / spoofing indicators (research)

Port of `python/spectrum_intel/receivers/gnss` (codes.py, acquire.py, plugin.py). Same method; see the Python
README for the reasoning and the validation (synthetic skies only).

## What it does
- **Acquires** GPS L1 C/A satellites in one block: 32 PRNs (IS-GPS-200 Gold codes), 4 ms coherent folds × 4
  non-coherent, χ²(2M) threshold at P_fa 1e-3 per PRN; each hit refined at 2 samples / chip → code phase, Doppler,
  C/N0 (with the measured 1.8 dB implementation loss); a cross-correlation guard drops ghosts of stronger PRNs.
- **Indicators:** in-band J/N, narrow lines near L1, spoofing S1–S4, and a swept jammer (personal privacy device,
  fs ≥ 30 MSps: STFT ridge over 400 µs, period from its autocorrelation, fold coherence ≥ 0.6).
- Plugin (`capture_level`, after `lte` for the receiver ppm): "GNSS: GPS L1 C/A (n satellites)", "GNSS jamming-like
  (GPS L1[, swept PPD])" or "GNSS spoofing-like (...)".

## C++ API (`namespace si::rx::gnss`, depends on core/ only)
```cpp
#include "receivers/gnss/acquire.hpp"
std::vector<uint8_t> c = ca_code(prn);                         // 1023 chips, 0/1
std::optional<Json> r = analyse(channelizer, ppm, exclude_hz);  // {"satellites", "interference", "indicators"}
std::optional<Json> sw = swept_ppd(channelizer, f_l1);
```

## Implementation notes
- Same dtypes as numpy; one FFT per fold shared by all PRNs; the reference vectors are in
  `tests/vectors/signal_gaps/gnss*` (`test_signal_gaps.cpp`).
