# nr (C++) — 5G NR SS/PBCH block receiver

Port of `python/spectrum_intel/receivers/nr` (plugin `nr_ssb`, capture level, `after = lte`, owns its band).

**What:** NR cells on the GSCN raster of the NR bands in the capture (n78, n28, n40, n41, n3, n1, n5, n8):
PSS N_ID2, differential SSS N_ID1 -> PCI (co-channel cells listed), PBCH DMRS i_SSB, PBCH LLRs (DMRS channel
estimate, MRC, combined per half-frame), BCH decode -> MIB (SFN, SCS common, k_SSB, CORESET#0 / search space
#0, cellBarred, ...).

**How:** `ssb.hpp`: `SSBSlice` (256 x SCS slice + FFT); `best()` scores the three PSS templates per CFO offset
with one batched inverse FFT (`fft::fft_rows`, as scipy's `ifft(axis=1)`), robust z with the exact radix-select
median of `np::median_inplace`. The GSCN probes run in parallel (`si::pmap`, ordered); one template bank per SCS is
shared by every probe. `decode_ssb_multi`: PSS instances, 4x240 SSB grids, joint (N_ID2, N_ID1), `pbch_llr`,
descrambling per i_SSB and L_max, BCH decode. Receiver ppm from `ctx.cache["ppm"]` (lte) or `["sweep_ppm"]`.

`bch.hpp` implements BCH decoding from the standards, no py3gpp code: TS 38.212 5.4.1 rate recovery
(repetition combining, sub-block de-interleaving), 5.3.1 polar code (reliability sequence Table 5.3.1.2-1,
frozen set for K = 56 incl. CRC24C, N = 512, n_PC = 0, input interleaver Table 5.3.1.1-1), 5.1 CRC-24C,
7.1.2 / 7.1.1 PBCH payload descrambling and de-interleaving. PSS/SSS m-sequences from TS 38.211 7.4.2,
DMRS from 7.4.1.4.1 (core `seq3gpp::gold/qpsk`).

**Polar decoder:** the Python reference calls `nrPolarDecode(..., L=8)`, but py3gpp's implementation ignores L
and runs a single-path successive-cancellation decoder with the min-sum f function. The C++ does exactly that
(SC, list size 1, min-sum), so decoded bits and CRC outcome are bit-identical to Python
(`tests/unit/test_rx_nr.cpp`, 12 noisy/noise PBCH blocks: 6 pass, 6 fail CRC in both). A CA-SCL decoder would
decode more weak MIBs but would break parity; it belongs in Python first.

**Report text:** steps and method strings say "py3gpp" exactly as the Python report does (parity); the C++
contains no py3gpp code.

**C++ API:** `si::rx::nr::{SSBSlice, decode_ssb_multi, pss_seq, sss_seq, pbch_dmrs, sss_scores, pbch_llr, peaks,
gscn, bch_decode, polar_decode, rate_recover_polar, crc24c_ok, polar_info_positions, polar_input_interleaver}`.

**Parity:** 16/16 inputs identical (cellular: n28 PCI 75 MIB at 770.88 MHz, n78 PCI 82 CRC failed at 3380 MHz;
no NR on the other inputs). Unit vectors: `tests/vectors/nr/gen_nr.py` (Python receiver + py3gpp as oracle):
Gold/QPSK, PSS/SSS/DMRS, reliability table, info set, interleaver, rate recovery, polar decode, BCH, SSS
scores, PBCH LLRs, peaks, GSCN, PSS templates.

**Deviations:** window energies (running float64 sum) and the SSS/PBCH dot products (sequential sums vs BLAS)
agree with numpy to ~1e-15 relative; ties in argsort of float64 metrics are resolved as a reversed stable sort
(numpy's quicksort order is not reproduced for exact ties, which do not occur on real data).
