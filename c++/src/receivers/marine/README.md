# marine (C++) — AIS (research)

Port of `python/spectrum_intel/receivers/marine` (ais.py, detector.py, plugin.py). Same method; see the Python
README for the reasoning and the validation (synthetic channels and a gpsd golden sentence).

## What it does
- **Detects and decodes** AIS on 161.975 / 162.025 MHz: GMSK 9600 bit/s via a ±10 kHz discriminator with one-bit
  integration, NRZI, HDLC flags and bit stuffing, CRC-16/X.25 (residue 0xF0B8) → message types 1/2/3, 4, 5, 18,
  21, 24 (MMSI, MID country (419 = India), position, SOG / COG / heading, name, call sign, ship type, dimensions,
  destination).
- Plugin: one signal for all AIS messages in the block, "AIS (n stations)", status success (every message
  CRC-verified).

## C++ API (`namespace si::rx::marine`, depends on core/ only)
```cpp
#include "receivers/marine/ais.hpp"
#include "receivers/marine/detector.hpp"
auto fr = ais::frames(bits); Json f = ais::fields(octets);
std::vector<si::Emission> em = detect(channelizer, ppm);
```

## Implementation notes
- vectors in `tests/vectors/signal_gaps/ais.json` (`test_signal_gaps.cpp`).
