// FLARM beacons: detection and the clear header only (python receivers/rc_link/flarm.py; research, PLAN T5).
#pragma once

#include "core/types.hpp"

namespace si::rx::rc_link::flarm {

constexpr double kChipRate = 100e3, kDev = 50e3;
bool in_band(double f_hz);
// python try_identify(y, ofs, chip_rate, max_err=2): the dict or null
Json try_identify(const cfvec& y, double ofs, double chip_rate = kChipRate, int max_err = 2);

}  // namespace si::rx::rc_link::flarm
