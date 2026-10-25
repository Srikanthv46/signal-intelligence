# dji (C++) — DJI OcuSync video, RC uplink, DroneID

Port of `python/spectrum_intel/receivers/dji` (detector.py, plugin.py). Same method; the Python README and
the detector.py docstring describe the physics.

## What it does
- **Video downlink** (OcuSync 2 / O3, 15 kHz SCS; O4, 30 kHz SCS): the CP metric of every (profile, centre) is
  folded on the nominal symbol period over the whole capture; local maxima over centre (z ≥ 8) are segmented
  into bursts on the known grid; the carrier comes from the CP phase + band position + DJI's empty DC bin. DJI
  test: 15 kHz, a length-601/1201 Zadoff-Chu slot symbol in DJI's mapping (DD coherence → root candidates →
  exact-sequence correlation over ±24 lags and ±2 bins); 30 kHz, a reference symbol repeated 4 symbols later.
  Features: bandwidth, frame period/phase, clock ppm (least-squares fit), CFO per burst, ZC root.
  The carrier estimate re-centres its slice until it stops moving; a link that failed the DJI test hides only
  the same profile at the same refined carrier (±0.5 MHz). When the 15 kHz ZC test fails on the capture-level
  grid, the link is re-timed burst by burst (`local_bursts`, `retimed_link`) and the exact test runs on the
  symbols whose DD coherence stands > 8 sigma out (`sparse_zc_evidence`, evidence `zc_any_symbol`, no clock
  ppm). Burst extent is robust to overlapping signals (`segment_multi` over the full band + 4 sub-bands,
  `bridge` across masked gaps ≤ 0.6 ms, band position from the per-bin median over symbols, in-band RC uplink
  hops — detected first — excluded); the DJI evidence and clock fit use only the full-band symbols
  (`ev_syms`). An exact ZC match > 0.3 above its reference is accepted without the slot statistic.
  O4: positional repeated-reference test (best rotation of slot positions {r, r+8}); low-SNR bursts completed
  from the TDD frame structure (`frame_complete`); links with evidence z 5–8 reported with
  `evidence_tier = like` and the plugin label `…-like`. A reported link that re-times mid-capture keeps its
  bursts on further grids (`regrid`: further peaks of its own CP fold, segmented at the same carrier, runs not
  overlapping the link); the O4 spectrum-edge walk closes dips < 0.6 MHz (18 MHz links were labelled 10 MHz).
  Per link (all its grids): the band of a link without a verified ZC length from `core/bursts::fit_band` (in-band
  power on vs off), the bursts extended / joined by `core/bursts::energy_bursts` (feature `band_from`).
- **RC uplink**: 1.4 MHz 15 kHz-OFDM hops from `core/ofdm15k` (scan, refine, merge, hop tracker); DJI if two
  length-73 ZC symbols of one root sit 6 symbols apart. IQ-imbalance images dropped; frame period over the hops.
- **DroneID**: 15.36 MSps slices every 3 MHz, full 9-symbol CP schedule score, ZC roots 600 / 147 in symbols 4
  and 6, carrier from the DC null + CP phase.
- Plugin `dji` (package `dji`): capture-level `EmissionReceiver` on ism_24 / ism_58 / unii_12 / unii_2e.

## C++ API (`namespace si::rx::dji`, depends on core/ only)
```cpp
#include "receivers/dji/detector.hpp"
std::vector<si::Emission> em = detect(iq, fs, fc);   // labels video_10MHz | video_20MHz | video_o4_* | uplink | droneid
cdvec z = zc_freq(root, L);                          // DJI ZC mapping
auto [corr, root] = zc_root_search(X, L);            // root -1 = None
```
`detail::` exposes the building blocks used by the tests (`FastSlicer`, `cp_metric`, `fold`, `robust_z`,
`segment`, `zc_dd`, `root_from_dd`, `droneid_cp_score`, `frame_period`).

## Implementation notes (numpy exactness)
- The channelizer is Python's `_FastChannelizer` (bins by index arithmetic, `np.fft.ifft` complex128) over
  `SpectralChannelizer::X()` with the spur notches; stateless, so it is the thread-safe `ofdm15k::Slicer`.
- numpy elides the temporary in `a * np.conj(b)` when it is ≥ 256 KiB and then computes `conj(b) * a`; the
  fma complex multiply is not symmetric, so `cp_metric`, the DroneID score and the frequency shift `y * exp(..)`
  switch operand order at 16384 complex128 elements, as numpy does.
- `np.exp(1j*θ)` is libm `cexp`; on Apple that equals `__sincos`, not `sin()` (1 ulp apart on ~0.01 % of
  angles), so `expj` uses `__sincos` there. Phases are built in numpy's order: `((-2π·f)·n)·(1/ofs)` (complex
  division by a real scalar multiplies by the reciprocal).
- `round(np.float64, n)` → `np::round`, Python floats → `py_round`, tracked per feature (e.g. `frame_period_ms`
  is a Python float when it comes from `_cycle`, an np.float64 from the 0.5 ms grid search).
- Speed: the scan folds with one cached fold index per profile (Python `_FOLD_IDX`) and shares |y|² between the
  15/30 kHz profiles of one slice; scan offsets, ZC-root hypotheses, O4 symbol correlations, uplink bursts and
  DroneID offsets run on `si::pmap` / `parallel_for` (ordered results).
- BLAS-based numpy calls (`np.vdot` in the O4 test, `np.linalg.lstsq` / `A @ sol` in the clock fit,
  `np.convolve`) agree to rounding only; their outputs are thresholded or rounded to 1–2 decimals.

## Deviation: DroneID payload
Python hands each DroneID burst to the external `droneid-decoder-reference/src/droneid_api` (turbo decode +
CRC). That repository has no licence file and no burst on testdata decodes, so it is not ported (it could not be
verified bit-exactly). C++ reports `decode_status = "import_failed: No module named 'droneid_api'"`, exactly
what Python reports when the decoder is absent, with `decoded = null` and confidence 0.85 (Python: 0.99 with a
payload, `serial` / `device_type` features). On a capture with a decodable DroneID burst the C++ report lacks
the payload.

## Parity and tests
- 2026-10-13 (further grids, O4 bandwidth): 28/28 identical with all receivers; tmpdata 8/8 identical (D_20260609-151042 ch0–3 now `video_o4_20MHz` over both grids).
- `c++/tools/parity.py --receivers dji`: 28/28 inputs identical (5 dji signals: Mini 2 20 MHz video, Mini 4
  Pro O4 video, 3 RC uplink hop groups; 275 energy signals). Also identical on the two 2026-10-10 captures of
  the missed links (`tmpdata/`: D_20260515-192148 and the 2484 MHz test_signal capture).
- `tests/unit/test_rx_dji.cpp` (vectors: `tests/vectors/dji/gen_dji.py`): ZC helpers, CP metric / fold /
  robust z / segmentation, frame period, a synthetic DroneID burst (score and emission) and `detect()` on the
  three DJI testdata captures — every emission and feature bit-identical to Python.
- Timing (3 DJI captures, arm64 Mac, detect only incl. capture load): Python 18.5 s; C++ 10.9 s (1 thread), 6.5 s (2), 4.9 s (4).
