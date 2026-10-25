// Hubsan toy / hobby drones (A7105; python receivers/rc_link/hubsan.py; research, synthetic only, PLAN G5).
#pragma once

#include <optional>

#include "core/types.hpp"

namespace si::rx::rc_link::hubsan {

constexpr double kRate = 100e3, kDev = 186e3;
constexpr double kF0 = 2400.000977e6, kStep = 0.5e6;
constexpr int kLen = 16;
// python try_decode(y, ofs, f_hz=None, rate=RATE): the dict or null
Json try_decode(const cfvec& y, double ofs, std::optional<double> f_hz = std::nullopt, double rate = kRate);

}  // namespace si::rx::rc_link::hubsan
