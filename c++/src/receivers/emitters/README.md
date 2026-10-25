# emitters (C++) — drone emitter association

Port of `python/spectrum_intel/receivers/emitters` (association.py, bursts.py, families.py, plugin.py). Same
method; see the Python README for the physics.

## What it does
Runs after every receiver. It takes the `Emission` records of the drone families (DJI OcuSync/DroneID, Herelink,
analog/digital VTX, RC FHSS, telemetry) and answers two questions: which bursts come from one radio (a hopping link
seen on several channels), and which radios form one drone system.
- **Fingerprints** (`measure_fingerprints`), for 15 kHz-SCS OFDM bursts:
  - CP-train symbol grid and CP-phase fine CFO;
  - integer subcarrier from a matched filter or the DC null;
  - re-timing on learned reference symbols (symbols that repeat burst to burst), to about 10 ns and about 20 Hz.
- **Pairwise score** (`same_emitter_score`): a sum of log-likelihood ratios from these cues:
  - frame-clock lattice timing;
  - oscillator ppm;
  - sample-clock ppm;
  - simultaneity veto;
  - power;
  - shape mismatch;
  - decoded identity, ZC root, sync word.
- **Grouping** (`associate`):
  - average-linkage clustering on the LLR matrix, never across a vetoed pair;
  - per emitter: hop pattern, frame fit (frame-clock ppm), oscillator fingerprint, confidence and evidence;
  - separation from each other emitter of the same type;
  - links to other-type emitters that share the oscillator.
- **Plugin output:** `report["emitters"]`, and `decoded.fields.emitter_ids` on each member signal.

## C++ API (`namespace si::rx::emitters`, depends on core/ only)
```cpp
#include "receivers/emitters/association.hpp"
#include "receivers/emitters/bursts.hpp"
std::vector<si::Emission> em = extract_bursts(iq, fs, fc);   // standalone 15 kHz-OFDM burst source
measure_fingerprints(iq, fs, fc, em);                        // refines em[i].features in place
std::vector<Emitter> emitters = associate(em, &assoc_models());
auto [llr, evidence] = same_emitter_score(em[0], em[1]);
```
The Python test-scene helpers (`synth_scene`, `inject_copy`) and `report()` are not ported.

## Implementation notes
- The operation order follows numpy:
  - `np::cmul` / `np::cdiv`, where numpy's complex divide multiplies by `1/b`;
  - pairwise sums;
  - the complex `exp` uses `__sincos` on Apple;
  - numpy temporary elision on ≥ 256 KiB arrays (operands swap in `y[:-N] * conj(y[N:])`).
- `math.hypot`, LAPACK lstsq/inv, BLAS convolve and polyfit agree to rounding, which is far below the printed
  precision.

## Parity / tests
- `tools/parity.py` (all receivers): 16/16 inputs identical. The emitters appear on these inputs:
  - DJI uplink: 5 bursts, reference-symbol timing;
  - OcuSync video;
  - Herelink control and video, with the frame clock;
  - analog VTX.
- `tests/unit/test_rx_emitters.cpp` checks three things against Python vectors:
  - the real-capture path on D_20260916-171343;
  - a synthetic two-radio capture;
  - a hand-built multi-family emission list, which exercises every cue, the veto, singletons and links.
