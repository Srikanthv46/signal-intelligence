// Bluetooth BR/EDR access codes (research, signal gaps S7; python: receivers/rc_link/bredr.py): (64,30) BCH sync
// word over the LAP -> a verified access code and the 24-bit LAP.
#pragma once

#include <cstdint>
#include <optional>
#include <tuple>
#include <vector>

#include "core/types.hpp"

namespace si::rx::rc_link::bredr {

constexpr uint32_t kGiac = 0x9E8B33, kLiac = 0x9E8B00;

// python syncword(lap): bit k = the k-th bit sent
uint64_t syncword(uint32_t lap);
// python find(bits): (position, LAP, bit errors) of the first access code
std::optional<std::tuple<int64_t, uint32_t, int>> find(const std::vector<int64_t>& bits);
// python try_decode(y, ofs): null when no access code
Json try_decode(const cfvec& y, double ofs);

}  // namespace si::rx::rc_link::bredr
