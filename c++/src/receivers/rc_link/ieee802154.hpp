// IEEE 802.15.4 2.4 GHz O-QPSK frames (python receivers/rc_link/ieee802154.py; research, synthetic only, PLAN T4).
#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include "core/types.hpp"

namespace si::rx::rc_link::ieee802154 {

constexpr double kChipRate = 2e6;
std::optional<int> channel(double f_hz, double tol_hz = 150e3);
uint16_t crc16_kermit(const uint8_t* p, size_t n);
Json parse_mac(const std::vector<uint8_t>& m);
// python demod(y, ofs): the dict or null (ofs must be 8 MHz: 4 samples per chip)
Json demod(const cfvec& y, double ofs);

}  // namespace si::rx::rc_link::ieee802154
