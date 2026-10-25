# adsb (C++) — Mode S / ADS-B 1090 MHz

Port of `python/spectrum_intel/receivers/adsb` (modes.py, messages.py, tables.py, tracks.py, plugin.py). Same
method and the same report; the Python README describes the protocol coverage.

## What it does
- **Demodulation** (`modes.cpp`, python `modes.demod`): the capture is channelised to 1090 MHz ± 2 MHz at
  8 MSps; on the magnitude envelope an 8 µs preamble template (pulses 0 / 1 / 3.5 / 4.5 µs) is correlated, local
  maxima (±3 samples, wrapping like `np.roll`) above 6× the median noise start a frame; PPM bits by comparing
  the two half-µs sums; CRC-24 (0xFFF409): DF17/18 residual 0, DF11 residual < 128 (interrogator code), DF0/4/5/
  16/20/21 kept only when the address/parity matches an address confirmed by DF11/17/18 (dump1090 rule).
- **Fields of every frame** (`messages.cpp`, python `messages.fields`): ADS-B identification / category,
  surface movement, baro and GNSS altitude + CPR fields, airborne velocity (all subtypes), emergency / ACAS RA
  broadcast, target state (TC29), operational status (TC31); DF11 interrogator, ACAS DF0/16 replies, DF4/5/20/21
  flight status, altitude, identity; Comm-B register inference and decoding (BDS 1,0 1,7 2,0 3,0 4,0 4,4 4,5
  5,0 6,0).
- **Per aircraft** (`tracks.cpp`, python `tracks.tracks`): CPR position (global from an even/odd pair ≤ 10 s,
  then local against the previous fix ≤ 30 s or the receiver position), > 1000 kt jumps flagged, kinematics,
  intent, integrity, capabilities, BDS50/60 disambiguation with the aircraft's own ADS-B state, wind and static
  air temperature, interrogation environment (ACAS, radar dwells and revisit, DF11 interrogators), flags.
- Plugin `adsb` (package `adsb`): capture level, one receiver signal (1089–1091 MHz) with every frame as a burst.

## Standards used (no pyModeS code)
The Python package calls pyModeS (GPL-3) for the encodings; the C++ library implements them from the public
standards in `codec.cpp`:
- ICAO Annex 10 Vol IV: CRC-24 parity (3.1.2.3.3), 13-bit altitude code incl. Gillham (3.1.2.6.5.4), identity
  code (3.1.2.6.7.1), 6-bit character set.
- DO-260B / ICAO Doc 9871: extended-squitter ME formats (typecode, identification, surface movement
  quantisation, airborne position altitude, airborne velocity subtypes 1–4, GNSS–baro difference, emergency
  Mode A code), CPR decoding (A.1.7: NL function, global even/odd and local decoding against a reference),
  Comm-B MB register layouts (BDS 1,7 capability bits, 2,0, 4,0, 4,4, 4,5, 5,0, 6,0) and their status bits.
- ICAO Doc 7488 standard atmosphere and the compressible pitot relations (TAS ↔ CAS ↔ Mach) for the BDS 6,0
  IAS/Mach consistency checks.

The conventions that are not in the standards but show in the report were matched to the Python output
black-box (pyModeS used only as a test oracle, through the Python receiver): value types (int vs float),
'_' for a space in callsigns (unassigned codes dropped in ADS-B identification, kept as '#' in BDS 2,0), the
register-inference plausibility limits (roll ≤ 50°, GS/TAS ≤ 600 kt, |TAS−GS| ≤ 200 kt, IAS ≤ 500 kt,
Mach ≤ 1, |vertical rate| ≤ 6000 ft/min, IAS within 20 kt of the Mach-derived CAS, wind ≤ 250 kt, temperature
−80…60 °C, reserved bits zero, BDS 4,4 requires a wind report, BDS 1,7 requires BDS 2,0 support), the
BDS50/60 choice by distance of the (TAS-converted) velocity vectors to the ADS-B ground vector, floored
division in local CPR, unrounded CPR positions, and the identity code rejected (`field_error`) when all 13
bits are equal.

## C++ API (`namespace si::rx::adsb`, depends on core/ only)
```cpp
#include "receivers/adsb/modes.hpp"      // demod(env, thr) -> frames (Json, Python's dict order); crc24, parity_residual
#include "receivers/adsb/messages.hpp"   // fields(frames): decoded fields in place
#include "receivers/adsb/tracks.hpp"     // tracks(frames_with_t, ref) -> aircraft (Json array)
#include "receivers/adsb/codec.hpp"      // the standard encodings (altitude, squawk, CPR, Comm-B, ...)
```
Frames and aircraft are `Json` objects built in Python's key order, since they go into the report verbatim.

## Implementation notes (exactness)
- The envelope is `float32` (`np.abs` of the complex64 slice, `np::cabs`), the threshold `float32`; the
  correlation is float64 (`np.correlate` of float32 with the float64 template; products are exact, so the
  sum order does not matter); bit sums and the SNR follow numpy's float32 reductions (`np::pairwise_sum`,
  `np::mean`) and NEP 50 promotion (`on / noise` in float32).
- CPR, haversine and the wind / temperature formulas keep Python's evaluation order (`math.radians` = x·(π/180),
  `x ** 2` = `pow`, floored `%` and `//`), `round()` via `py_round` / `py_round_int`.
- Python raises `KeyError: 'BDS50,BDS60'` in `tracks._resolve_commb` when an ambiguous BDS50/60 reply cannot be
  resolved (no ground track/speed in the BDS 5,0 part, no heading or airspeed in the 6,0 part, or no usable
  distance) while the aircraft has an ADS-B velocity and altitude within 10 s; the whole plugin then fails
  (`errors.adsb`). The C++ throws the same message, so the report is the same. This is a Python bug worth
  fixing there first (it hits ~25 % of the synthetic multi-aircraft scenarios); C++ follows Python.
- The report's `method` and `steps` strings still say "pyModeS", as Python's do (report parity).

## Parity and tests
- `c++/tools/parity.py --receivers adsb`: 16/16 inputs identical (adsb 1/1: real_1090_commb, 6 frames
  DF17/20/21, BDS 4,0 / 5,0 / ambiguous 5,0-6,0; 142 energy signals).
- `tests/unit/test_rx_adsb.cpp` (vectors: `tests/vectors/adsb/gen_adsb.py`): CRC and Annex 10 codes; demod on
  the real_1090_commb envelope (full and windows); `fields()` on 2427 frames (random and structured frames of
  every DF, typecode, velocity subtype, TC28/29/31 subtype and Comm-B register, incl. a `field_error`);
  `tracks()` on 110 synthetic scenarios (CPR global/local incl. surface, poles, ±180°, NL edges, implausible
  jumps, BDS50/60 resolution, met, interrogation, flags; 30 raising the Python KeyError); 800 local + 800
  global CPR decodes and 656 NL values (zone edges ± 1 ulp, 0°, ±87°); 1040 BDS50/60 resolutions — all
  identical to Python's JSON. The generator also runs larger fuzz sets (`ADSB_SEED=7 ADSB_SCALE=8`: 19k frames,
  1280 scenarios, 12k+12k CPR, 6.6k BDS50/60 — all identical; not committed).
- Timing (real_1090_commb, arm64 Mac): receiver 0.02 s in C++ vs 0.07 s in Python (mostly the channeliser slice).

## UAT 978 (2026-10-23, PLAN G14; research)
`uat.{hpp,cpp}` and the capture-level plugin `uat` (python `UatDecoder`): CPFSK phase-change slicer at 4 samples / bit,
sync within 4 bit errors, `core/rs` Reed-Solomon (own implementation), DO-282B fields. == python on
`tests/vectors/g10_g14` (RS cases, synthetic demod, dump978 sample messages).

## SSR Mode A/C replies (research, signal gaps S11)
`modeac.{hpp,cpp}` (python `modeac.py`): framing pulses 20.3 µs apart, 13 positions on the 1.45 µs grid, SPI, Gillham
altitude; `mode_ac` fields on the Mode S signal, or an "SSR Mode A/C replies[-like]" signal without Mode S. Vectors:
`tests/vectors/signal_gaps/modeac*` (Gillham table on all 2048 words, decode / summarise, a capture through the plugin).
