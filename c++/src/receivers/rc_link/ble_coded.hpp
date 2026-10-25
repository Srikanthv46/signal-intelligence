// BLE LE Coded PHY (Long Range, S=8 / S=2) for Remote ID over Bluetooth 5 (python receivers/rc_link/ble_coded.py;
// research, synthetic validation only, PLAN T2b).
#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include "core/types.hpp"

namespace si::rx::rc_link::coded {

std::vector<uint8_t> conv_encode(const std::vector<uint8_t>& bits);
std::vector<uint8_t> known_symbols(const std::vector<uint8_t>& aa_bits, int ci = 0);
// soft Viterbi, K = 4: llr is n x 2 (row-major), positive -> 1
std::vector<uint8_t> viterbi(const dvec& llr);
// python find_coded(p, aa_bits, top=2, min_ncorr=0.4)
std::vector<std::pair<int64_t, double>> find_coded(const cdvec& p, const std::vector<uint8_t>& aa_bits);

struct Decoded {
  int64_t aa_bit_errors = 0;
  int ci = -1;
  bool has_bits = false;
  std::vector<uint8_t> bits;  // block 2 (PDU + CRC, whitened)
  double h = 0.5;
};
// python decode(y, s0, aa_bits)
std::optional<Decoded> decode(const cfvec& y, int64_t s0, const std::vector<uint8_t>& aa_bits, int max_pdu_octets = 257);
// python parse_extended(by): the dict with adv_data as raw bytes in `data`
Json parse_extended(const std::vector<uint8_t>& by, std::vector<uint8_t>& data);

}  // namespace si::rx::rc_link::coded
