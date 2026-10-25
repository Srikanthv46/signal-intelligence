# pmr_digital (C++) — digital land-mobile radio (TETRA, DMR, NXDN)

Port of `python/spectrum_intel/receivers/pmr_digital` (tetra.py, fsk4.py, detector.py, plugin.py).

## What it does
- **TETRA downlink**: 18 ksym/s line of the matched-filtered envelope; pi/4-DQPSK (RRC 0.35, 4th-power carrier
  estimate, the quarter-turn ambiguity resolved by the training sequences); synchronisation (38-bit) and normal
  (22-bit) training sequences on the 510-bit slot grid; **BSCH** of every synchronisation burst (descramble e=0,
  de-interleave K=120/a=11, RCPC 2/3 -> rate-1/4 K=5 soft Viterbi, CRC-CCITT) -> SYNC PDU: **MCC, MNC, colour
  code**, system code, TN/FN/MN.
- **DMR**: FM discriminator at 48 kSps; 9 SYNC patterns; per burst RRC 0.2 timing/levels from the sync (correlation
  >= 0.9, level 0.5–1.5 x 648 Hz); slot type Golay(20,8) -> **colour code, data type**; BPTC(196,96) with
  Hamming(13,9)/(15,11) single-error correction; CRC-CCITT with the data-type mask; CSBK opcode / FID when valid.
- **NXDN**: frame sync word 0xCDF59 at 4800 / 2400 sym/s on the 192-symbol frame grid.
- The plugin drops voice_fm's analog reading of a carrier verified here (`arbitrate`).

## C++ API (`namespace si::rx::pmr_digital`, depends on core/ only)
```cpp
#include "receivers/pmr_digital/detector.hpp"
std::vector<si::Emission> em = detect(channelizer, band /* optional (lo, hi) */, exclude_hz);
// lower level: tetra::demod / tetra::find_training / tetra::decode_bsch, fsk4::discriminator / dmr / nxdn
```

## Implementation notes
- Same dtypes as numpy: complex64 channels from `core/dsp` Channelizer; complex128 matched filter and float64
  soft bits for TETRA; float32 discriminator and float64 filters for 4FSK; complex products / divisions / `exp`
  as numpy (`np::cmul`, `np::cdiv`, libm `cexp`; numpy's temporary-elision operand order on >= 256 KiB arrays).
- Python's ties are explicit (TETRA identity: the most frequent, first seen; DMR colour code: the most frequent,
  lowest), so both implementations choose the same.

## Parity
Checked with `c++/tools/parity.py --receivers pmr_digital` on all testdata (TETRA, DMR and NXDN live captures in
`testdata/live`; null on the rest).
