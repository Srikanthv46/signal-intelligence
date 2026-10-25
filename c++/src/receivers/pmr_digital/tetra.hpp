// TETRA downlink (ETSI EN 300 392-2) (python: receivers/pmr_digital/tetra.py): pi/4-DQPSK demodulation,
// synchronisation bursts, BSCH decoding (SYNC PDU: MCC, MNC, colour code, TN/FN/MN). Depends on core/ only.
#pragma once

#include <optional>
#include <utility>
#include <vector>

#include "core/types.hpp"

namespace si::rx::pmr_digital::tetra {

constexpr double kSr = 18000.0;     // symbol rate
constexpr int kSps = 8;
constexpr int64_t kSlotBits = 510;  // 255 symbols per timeslot
constexpr int kBschBits = 120;      // SB block 1 (type-5 bits), right before the sync training sequence

// python symbol_line(y, fs): 18 ksym/s line of the matched-filtered envelope power over its neighbourhood
double symbol_line(const cfvec& y, double fs);

struct Demod {
  std::vector<int8_t> bits;  // hard bits
  dvec soft;                 // + = 1
  double err = 0;            // rms phase error (rad)
  double foff = 0;           // carrier offset of y (Hz)
};
// python demod(y, quarter)
Demod demod(const cfvec& y, int quarter = 0, int block_syms = 2000);

struct Training {
  std::vector<int64_t> sync, normal;  // bit positions (sync: 38-bit y <= 4 errors; normal: n, p <= 2 errors)
};
// python find_training(bits)
Training find_training(const std::vector<int8_t>& bits);
// python slot_consistent(pos): largest number of hits sharing one phase modulo the 510-bit slot
int64_t slot_consistent(const std::vector<int64_t>& pos);

// python decode_bsch(soft120): the SYNC PDU fields (ordered as python) when the CRC passes
std::optional<Json> decode_bsch(const double* soft120);

}  // namespace si::rx::pmr_digital::tetra
