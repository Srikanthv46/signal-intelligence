# lte (C++) — LTE downlink cell search

Port of `python/spectrum_intel/receivers/lte` (plugin `lte`, capture level, owns its band).

**What:** every LTE cell on the 100 kHz EARFCN raster inside an always-on energy signal of an LTE band
(B28/B5/B8/B3/B1/B40/B41): PSS N_ID2, SSS N_ID1 (FDD/TDD x normal/extended CP) -> PCI, duplex, CP,
subframe-0 time, carrier frequency error, channel bandwidth from energy. MIB/SIB are not decoded (as Python).

**How (`cellsearch.hpp`):** `Slice` = a 3.84 MSps channel slice (core `Channelizer::slice`) with its FFT; a
frequency offset is an integer-bin rotation; the PSS correlation (`|ifft(Y*conj(T))|^2 / energy`) is folded
over 5 ms and scored as a robust z (median/MAD) and isolation. `search_cells` probes the rasters in parallel
(`si::pmap`, ordered), NMS per (carrier, N_ID2), refines CFO +-3 kHz, `sss_decode` accumulates the 168 x
{SF0, SF5} SSS bank over half-frames; ZC-ambiguity suppression as Python. `wide_ppm_search` (+-500 kHz,
offsets in parallel) measures the receiver ppm when the caller gives none; `bandwidth_from_energy` reads the
coarse spectrogram.

**Cache contract:** publishes `ctx.cache["ppm"]` (Json `{"ppm", "z", "carrier_hz", "raster_hz", "source"}`, or the
caller's `sweep_ppm`, or null) for `nr` and `gsm` (their `after()` names `lte`).

**C++ API:** `si::rx::lte::{Slice, search_cells, sss_decode, sss_decisive, instances, wide_ppm_search,
bandwidth_from_energy, template_bank}`.

**Parity:** `c++/tools/parity.py --receivers lte,nr,gsm`: 16/16 inputs identical. No testdata input holds an LTE
cell, so the decode path is checked by `tests/unit/test_rx_lte.cpp` (synthetic FDD cell from
`tests/vectors/lte/gen_lte.py`: cells, SSS z, instances, wide ppm search, bandwidth) and by a synthetic
15.36 MSps B8 capture compared end to end with and without a caller ppm (identical reports, PCI 172).

**Deviations:** sliding-window energies use a running float64 sum (re-anchored every 1024 samples) instead of
numpy's dot-product convolution, and the SSS bank product is summed sequentially (numpy einsum order): both
agree to ~1e-15 relative; no decision moved on the vectors or testdata.

## 2026-10-10 (ported with Python)
- Candidates: LTE-band energy >= 3 dB, continuous or bursty; continuous >= 6 dB searched as before, the others
  (widest first, rasters already searched skipped) report only SSS-resolved cells; cells only on an LTE band raster.
- SSS: maximum-ratio combining with the PSS channel estimate smoothed over 15 subcarriers per side of DC.
- PSS offsets 0 Hz when the receiver error is known and 2.5 ppm of the carrier is < 3 kHz, else +-5 kHz.
- Slice FFT length `next_fast_len(L + N)` in complex128, window energy by cumulative sum (as numpy).
- `subframe0_time_s` is null when no subframe-0 PSS instance is in the block (Python raised and lost every cell).
