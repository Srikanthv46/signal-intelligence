// mLRS frames: passive identification (python receivers/rc_link/mlrs.py; research, synthetic only, PLAN R3 / T6).
#pragma once

#include <cstdint>
#include <vector>

#include "core/types.hpp"

namespace si::rx::rc_link::mlrs {

constexpr int kFrameLen = 91;
struct Mode { const char* name; const char* band; int sf; double bw; int cr_den, preamble, frame; int period_us; bool decodable; };
const std::vector<Mode>& modes();   // python MLRS_LORA (dict order)
uint16_t crc16_mcrf4xx(const uint8_t* p, size_t n, unsigned init = 0xFFFF);
// python parse(frame): the dict or null
Json parse(const std::vector<uint8_t>& frame);

}  // namespace si::rx::rc_link::mlrs
