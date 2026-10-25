// FlySky AFHDS 2A (A7105) passive identification (python receivers/rc_link/afhds2a.py; research, synthetic only,
// PLAN R7).
#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include "core/types.hpp"

namespace si::rx::rc_link::afhds2a {

constexpr double kRate = 500e3, kDev = 186e3;
constexpr int kPayload = 38;
uint16_t crc16(const uint8_t* p, size_t n, unsigned poly, unsigned init);
// python afhds1_hops(id3): the 16 AFHDS 1 hop channels of TX ID byte 0 (as sent)
std::vector<int> afhds1_hops(unsigned id3);
constexpr double kA7105F0 = 2400.000977e6, kA7105Step = 0.5e6;
// python try_identify(y, ofs, rate=RATE, f_hz=None): the dict or null (AFHDS 1 first, then AFHDS 2A)
Json try_identify(const cfvec& y, double ofs, double rate = kRate, std::optional<double> f_hz = std::nullopt);

}  // namespace si::rx::rc_link::afhds2a
