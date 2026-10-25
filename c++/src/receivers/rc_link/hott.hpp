// Graupner HoTT 2.4 GHz RC link (CC2500 MSK + FEC; python receivers/rc_link/hott.py; research, synthetic only,
// PLAN G12).
#pragma once

#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

#include "core/types.hpp"

namespace si::rx::rc_link::hott {

constexpr double kRate = 249.939e3, kDev = kRate / 4;
constexpr double kF0 = 2404.0264e6, kStep = 1.014038e6;
constexpr int kLenTx = 50, kLenTlm = 22;
// python coded_len(n_data): coded bits on air
int coded_len(int n_data);
// python fec_encode(data): interleaved coded bits (terminator added)
std::vector<uint8_t> fec_encode(const std::vector<uint8_t>& data);
// python fec_decode(bits, n_bytes): (decoded bytes, path metric)
std::pair<std::vector<uint8_t>, int> fec_decode(const std::vector<uint8_t>& bits, int n_bytes);
// python try_decode(y, ofs, f_hz=None, rate=RATE): the dict or null (CRC-verified only)
Json try_decode(const cfvec& y, double ofs, std::optional<double> f_hz = std::nullopt, double rate = kRate);

}  // namespace si::rx::rc_link::hott
