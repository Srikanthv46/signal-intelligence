# navaids (C++) — ILS glide slope / localizer and VOR (research)

Port of `python/spectrum_intel/receivers/navaids` (detector.py, plugin.py). Same method; see the Python README for
the reasoning and the validation (synthetic carriers only).

## What it does
- **Channels:** ILS glide path 329.15–335.00 MHz and localizer 108.10–111.95 MHz (ICAO pairing), VOR 108–117.95 MHz.
- **ILS:** 90 / 150 Hz least squares over whole 1/30 s windows with a common tone scale (0.97–1.03) → m90, m150,
  DDM, SDM, side, two-frequency installations.
- **VOR:** 30 Hz AM and the 9960 Hz subcarrier's 30 Hz FM → radial, AM depth, deviation.
- Plugin: "ILS glide slope (...)", "ILS localizer (...)", "VOR beacon (...)" ("-like" short of the validity rules);
  owns its band.

## C++ API (`namespace si::rx::navaids`, depends on core/ only)
```cpp
#include "receivers/navaids/detector.hpp"
std::vector<Channel> ch = channels();
std::vector<si::Emission> em = detect(channelizer, ppm);
```

## Implementation notes
- vectors in `tests/vectors/signal_gaps/navaids.json` (`test_signal_gaps.cpp`).
