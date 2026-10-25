# digital_vtx (C++) — Herelink, HDZero and any CP-OFDM video/data link

Port of `python/spectrum_intel/receivers/digital_vtx` (detector.py, plugin.py). Same method; the Python
README and the detector.py docstring describe the physics.

## What it does
- **Herelink** (CubePilot, LTE Rel-8 downlink PHY in the 2.4 GHz ISM band): the capture FFT (spur bins
  notched) is sliced at every 10 kHz carrier hypothesis to 1.92 MSps; the energy-normalised PSS matched filter
  (roots 25/29/34) is clipped, folded on the 5 ms half-frame, dilated ±1 sample and normalised by
  max(row median, global median). Candidates (z > 7, one per ±400 kHz) are refined: carrier/root search over
  ±28 subcarriers × 2.5 kHz then 250 Hz, PSS instances on a 5 ms grid (least-squares frame period and clock
  ppm), SSS over all half-frames → PCI and 10 ms frame phase, a PSS-purity test (rejects DJI's 5 ms ZC slot),
  20/10 MHz mode from the PSD shape, per-subframe power → control (sf 0/5), video (DL runs) and uplink
  (activity in the silent subframes) bursts, CP-phase CFO (link and per burst).
- **Any CP-OFDM link (Tu >= 10 us), blind** (`blind.cpp`, python `blind.py`): time-AGC'd lag-Tu ACF scan over
  10/20 MHz slices, symbol-rate line / per-window fold confirmation, sub-sample Tu with integer-lag slice rates,
  CP cross-spectrum band, bursts on local symbol grids, TDD frame completion, then band and bursts from in-band
  power (`core/bursts`: on-vs-off rectangle fit, energy runs seeded by the CP runs), duplicate / IQ-image removal;
  numerology table with roles (drone / infrastructure / wifi).
- **Hopping narrowband CP-OFDM links (RC uplinks), blind** (`blind_hops`): `core/ofdm15k` window scan at every blind
  link's numerology and at 15 / 30 kHz; per hit envelope edges, on-vs-off band (rectangle fit, re-centred once),
  CP-fold verification (z > 6); hop cadence tracking; one numerology per hop; dwell / width consistency; hop
  channels. Emissions `ofdm_hop_<SCS>kHz`, one report signal per channel.
- **Transmit antennas** (`streams.cpp`, python `streams.py`): on blind links with SNR ≥ 20 dB, the slice kept by
  `blind_links` (`Link::slice`, `keep_slice_snr`). Per burst it folds the CP grid and takes the CP-phase CFO, then
  computes per-subcarrier kurtosis and the Alamouti pair power correlation (time / per-block frequency) on runs of
  ≥ 20 full-band symbols. Result in `features.tx_antennas`: `2tx_stbc_like` / `2tx_sfbc_like` / `1_stream_like` /
  `unknown`, with a label suffix `2-TX STBC-like`. Grid rounding is half to even (`nearbyint`), as in numpy.
- `hdzero.cpp` (python `hdzero.py`, 2026-10-22): `slice_plan`, `hdzero_links` (CP fold on the fractional grid,
  fractional CFO, per-symbol FFT with the half-sample start correction, continual-pilot coherence over the integer
  offsets, preamble blocks → frame period) and `detail::hdzero_detect` (family `hdzero`, label `hdzero_video`). The
  plugin labels it `HDZero digital FPV video link` with `catalog::hdzero_model`. Unit test against
  `vectors/digital_vtx` (`hdzero.cf32`: a synthetic link at R1 + 20 kHz, ~10 dB): slice plans, search log, links and
  detect() features identical to Python. `parity.py` 2026-10-22: 48 / 48 inputs identical, incl. LIVE_HDZERO_5658MHz.
- Plugin `digital_vtx` (package `digital_vtx`): capture-level `EmissionReceiver`, in range of the 2.4 and
  5 GHz FPV ranges; keeps Herelink above 2400 MHz (confidence ≥ 0.5) and unknown OFDM with confidence ≥ 0.6,
  blind links outside cellular bands whose role is not wifi; groups blind links by (family, label, carrier);
  `arbitrate` attributes a blind link to a DJI video link / Herelink of the same SCS covering > half its band, a
  blind hop channel to a DJI uplink hop it overlaps (or to a DJI video link whose bursts and band hold most of its
  hops). The analog_vtx plugin absorbs weak digital_vtx detections inside a VTX's FM
  spectrum (done there, as in Python).

## C++ API (`namespace si::rx::digital_vtx`, depends on core/ only)
```cpp
#include "receivers/digital_vtx/detector.hpp"
std::vector<si::Emission> em = detect(iq, fs, fc, /*report_wifi=*/true);  // herelink | digital_vtx | wifi_ofdm
auto links = blind::blind_links(slicer, fs, fc, usable, noise);   // blind.hpp building blocks
auto [name, hits] = classify_numerology(tu_us, bw_mhz);
```
Files: `detector.cpp` (entry, slice, PSS template), `herelink.cpp` (fold search, refinement, bursts),
`generic.cpp` (numerology table, blind links -> emissions), `blind.cpp` (blind CP-OFDM detector), `common.hpp` (internal helpers), `plugin.cpp`.

## Exactness notes
- Fold search: float32 throughout as numpy/scipy (complex64 slices, scipy single-precision IFFT, float64
  cumsum of the float32 energy, float32 mean/clip/fold sum); rows go through pocketfft's multi-row transform
  (`fft::fft_rows`) in batches of 16, per-row identical to scipy's batched `ifft(axis=1)`.
- `np.exp(1j*θ)` via `__sincos` (cexp); θ as numpy's complex multiplies/divisions leave it
  (`((2π·df)·n)·(1/ofs)` in the carrier search); `round()` follows the numpy-vs-Python float type of each value.
- BLAS-backed dot products (`np.correlate`, `np.vdot`, `np.convolve`) are summed sequentially (Python sums
  them in OpenBLAS order): last-bit differences only, no decision on testdata depends on them.

## Parity
`parity.py` (all receivers): 28/28 inputs identical, 2026-10-13 (digital_vtx 2/2: Herelink PCI 411, the
Phantom 4 Lightbridge-like link); tmpdata (`dji`, `digital_vtx`): 8/8 identical, incl. the PMP 450-like links.
Unit tests: `tests/unit/test_rx_digital_vtx.cpp` (vectors from `tests/vectors/digital_vtx/gen_digital_vtx.py`):
PSS template, slice, fold search, carrier search, blind scan / confirmation / segmentation / frame period and
completion / band, numerology table and detect() on a synthetic capture; tx_streams on two synthetic link slices
(1 TX 64-QAM → `1_stream_like`, STBC 16-QAM → `2tx_stbc_like`; values to 1e-3, notes identical).
2026-10-21 (tx_antennas): `parity.py --inputs LIVE_OFDM_VTX_5805MHz`: identical (digital_vtx 2/2, incl.
`tx_antennas` = `2tx_stbc_like`).

Timing (Herelink capture, 2 M samples at 61.44 MSps): receiver 5.4 s with 2 threads (3.8 s with 4) vs
Python detect 9.8 s (4 threads, scipy workers = all cores). The PSS fold (16 k IFFTs of 62 500 points) is
~90 % of it.
