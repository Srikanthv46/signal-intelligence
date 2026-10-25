# broadcast_fm (C++) — FM broadcast, stereo pilot and RDS (research)

Port of `python/spectrum_intel/receivers/broadcast_fm` (rds.py, detector.py, plugin.py). Same method; see the
Python README for the reasoning and the validation (synthetic stations and the IEC 62106 code tables).

## What it does
- **Detects** carriers on the 100 kHz raster 87.6–107.9 MHz (≥ 10 dB over the band floor, local maxima), measures
  the deviation (15–110 kHz) and the 19 kHz stereo pilot on the MPX.
- **Decodes RDS** (≥ 45 ms): 57 kHz BPSK (3 × the measured pilot), biphase matched filter at 1187.5 bit/s,
  differential decoding, (26,16) block sync on the offset words, error-free blocks only → PI, PTY, TP / TA, PS name,
  RadioText (2A), clock time (4A).
- Plugin: "FM broadcast (<channel> MHz[, stereo][, RDS '<PS>' or RDS PI <pi>])", status success with RDS, partial
  without.

## C++ API (`namespace si::rx::broadcast_fm`, depends on core/ only)
```cpp
#include "receivers/broadcast_fm/detector.hpp"
#include "receivers/broadcast_fm/rds.hpp"
std::vector<si::Emission> em = detect(channelizer);
auto blks = rds::blocks(bits); Json info = rds::decode(rds::groups(blks));
```

## Implementation notes
- float32 discriminator as numpy (`np::cmul` with the temporary-elision operand order); vectors in
  `tests/vectors/signal_gaps/broadcast_fm.json`.
