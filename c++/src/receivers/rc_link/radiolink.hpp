// RadioLink 2.4 GHz RC links (CC2500; python receivers/rc_link/radiolink.py; research, synthetic only, PLAN G4).
#pragma once

#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

#include "core/types.hpp"

namespace si::rx::rc_link::radiolink {

constexpr double kRate = 37.987e3, kDev = 38.09e3;
constexpr double kF0 = 2401.99994e6, kStep = 399.902e3;
constexpr int kMaxLen = 0x21;
// python pn9(n): CC2500 PN9 whitening bytes
std::vector<uint8_t> pn9(size_t n);
// python hop_table(tx_id): the 16 hop channels and the replacement channel
std::pair<std::vector<int>, int> hop_table(const uint8_t* tx_id);
// python try_decode(y, ofs, f_hz=None, rate=RATE): the dict or null (CRC-verified only)
Json try_decode(const cfvec& y, double ofs, std::optional<double> f_hz = std::nullopt, double rate = kRate);

}  // namespace si::rx::rc_link::radiolink
