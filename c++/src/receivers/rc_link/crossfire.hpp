// TBS Crossfire 150 Hz FSK air packets (python receivers/rc_link/crossfire.py; research, synthetic only, PLAN R9).
#pragma once

#include "core/types.hpp"

namespace si::rx::rc_link::crossfire {

constexpr double kRate = 85.1e3, kDev = 42.48e3;
// python grid_channel(f_hz): {band, channel, link} or null
Json grid_channel(double f_hz, double tol_hz = 20e3);
// python try_decode(y, ofs, f_hz): the dict or null
Json try_decode(const cfvec& y, double ofs, double f_hz);

}  // namespace si::rx::rc_link::crossfire
