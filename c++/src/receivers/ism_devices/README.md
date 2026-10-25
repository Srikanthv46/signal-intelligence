# ism_devices (C++) — 315 / 433 MHz remotes and TPMS (research)

Port of `python/spectrum_intel/receivers/ism_devices` (tpms.py, ook.py, detector.py, plugin.py). Same method; see
the Python README for the reasoning and the validation (synthetic devices only).

## What it does
- **Candidates:** bursty carriers in 314–316 and 433.05–434.79 MHz from 256-point frames normalised per bin.
- **OOK** (envelope pulses, debounced): EV1527 / PT2262 fixed codes (repeats as evidence, tri-state when it fits),
  KeeLoq HCS clear serial / buttons / battery (the encrypted hopping code is not reported).
- **FSK** (two-level discriminator, 52 µs symbols): TPMS Continental / VDO and Toyota PMV-C210, CRC-8 verified.
- Plugin: "TPMS (..., ID ...)", "OOK fixed-code remote (EV1527 / PT2262-like, ID ...)", "KeeLoq code-hopping
  remote-like (serial ...)", "OOK / FSK device-like (unrecognised protocol)"; `arbitrate` attributes another
  receiver's unverified guess on the same burst (e.g. rc_link "SiK-like") to the decoded device.

## C++ API (`namespace si::rx::ism_devices`, depends on core/ only)
```cpp
#include "receivers/ism_devices/detector.hpp"
#include "receivers/ism_devices/ook.hpp"
#include "receivers/ism_devices/tpms.hpp"
std::vector<si::Emission> em = detect(channelizer);
std::optional<Json> t = tpms::decode(symbols);          // either polarity
std::optional<Json> r = ook::fixed_code(pulses);        // ook::keeloq(pulses)
```

## Implementation notes
- float32 envelope / discriminator angle as numpy (`np::cabs`, `np::cmul` with the temporary-elision operand
  order), float64 moving averages by cumulative sum; vectors in `tests/vectors/signal_gaps/ism_*`.
