// FrSky D8 / D16 (CC2500, ACCST) passive identification (python receivers/rc_link/cc2500.py; research, synthetic only,
// PLAN R6).
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "core/types.hpp"

namespace si::rx::rc_link::cc2500 {

// python RATES (nominal bit rates, in dict order)
struct Rate { const char* variant; double bps; };
constexpr Rate kRates[] = {{"D8", 31.04e3}, {"D16 v1", 70.0e3}, {"D16 v2", 76.97e3}, {"D16 LBT", 99.98e3}};
// python RATES.get(variant.replace(" LBT v1", " LBT"), RATES["D16 v1"])
double nominal_rate(const std::string& variant);
uint16_t crc16_cc2500(const uint8_t* p, size_t n);
// python _find_sync(bits, sync_bytes): positions after the sync (<= 2 bit errors, either polarity) after >= 6 of 7
// alternating preamble transitions
struct SyncHit { int64_t pos; bool inv; int err; };
std::vector<SyncHit> find_sync(const std::vector<uint8_t>& bits, const uint8_t* sync, int nsync);
uint16_t crc16_d16(const uint8_t* p, size_t n);
// python try_identify(y, ofs, rate): the dict or null
Json try_identify(const cfvec& y, double ofs, double rate);

}  // namespace si::rx::rc_link::cc2500
