// XN297 toy-drone links and the Bayang protocol (python receivers/rc_link/xn297.py; research, synthetic only, PLAN R8).
#pragma once

#include <cstdint>
#include <vector>

#include "core/types.hpp"

namespace si::rx::rc_link::xn297 {

constexpr double kRate = 1e6;
const std::vector<uint8_t>& scramble();
uint16_t crc16(const uint8_t* p, size_t n, unsigned init = 0xB5D2);
// python try_decode(y, ofs, rate=RATE): the dict or null
Json try_decode(const cfvec& y, double ofs, double rate = kRate);

}  // namespace si::rx::rc_link::xn297
