# gsm (C++) — GSM BCCH carrier receiver

Port of `python/spectrum_intel/receivers/gsm` (plugin `gsm`, capture level, `after = lte`, owns its band).

**What:** every 200 kHz ARFCN inside GSM-900 / DCS-1800 downlink energy regions (100 kHz - 6 MHz wide):
FCCH (pure +67.708 kHz tone) -> GSM BCCH carrier confirmed, residual CFO; SCH one TDMA frame later
(differential GMSK detection against the 64-bit training sequence, rate-1/2 K=5 Viterbi, CRC-10) -> BSIC
(NCC, BCC), T1/T2/T3', frame number.

**How (`sch.hpp`):** `find_fcch` (sliding tone purity, plateau start), `demod_sch` (+-20 symbol search, both
polarities), `viterbi`, `crc10_ok`, `sch_fields`. The plugin channelises each ARFCN at 4 samples/symbol with
the receiver ppm from `ctx.cache["ppm"]` (lte) or `["sweep_ppm"]`; ARFCNs are decoded in parallel
(`si::pmap`, ordered) after the same in-order de-duplication as Python.

**C++ API:** `si::rx::gsm::{find_fcch, demod_sch, viterbi, crc10_ok, sch_fields, arfcn_freq}`.

**Parity:** 16/16 inputs identical (cellular 928 MHz: 8 BCCH carriers, FCCH found, SCH not decoded, as Python).
The SCH success path is checked by `tests/unit/test_rx_gsm.cpp` on synthetic FCCH + SCH bursts
(`tests/vectors/gsm/gen_gsm.py`: FCCH start/CFO/purity, soft bits, Viterbi bits, CRC, BSIC/FN fields).

**Deviations:** the FCCH sliding sums use running sums re-anchored every 64 samples instead of numpy's
convolution (~1e-13 relative); no decision moved on the vectors or testdata.
