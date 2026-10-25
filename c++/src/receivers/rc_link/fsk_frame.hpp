// Blind frame structure of any (G)FSK packet (python receivers/rc_link/fsk_frame.py; research, PLAN F1).
#pragma once

#include "core/types.hpp"

namespace si::rx::rc_link::fsk_frame {

constexpr int kMinPre = 12, kSyncBits = 32;
// python frame(y, ofs, rate): the dict or null
Json frame(const cfvec& y, double ofs, double rate);

}  // namespace si::rx::rc_link::fsk_frame
