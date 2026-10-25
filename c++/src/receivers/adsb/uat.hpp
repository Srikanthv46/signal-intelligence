// UAT 978 MHz ADS-B downlink (DO-282B): demodulation, Reed-Solomon correction and message decode (python
// receivers/adsb/uat.py; research, synthetic + dump978 sample messages only, PLAN G14).
#pragma once

#include <cstdint>
#include <vector>

#include "core/types.hpp"

namespace si::rx::adsb::uat {

constexpr double kFUat = 978e6, kRate = 1.041667e6, kDev = 312.5e3;
constexpr int kSpb = 4;
constexpr double kFs = kSpb * kRate;
// python parse(m): corrected payload (18 / 34 bytes) -> fields
Json parse(const std::vector<uint8_t>& m);
// python demod(y): RS-corrected downlink messages at kFs
std::vector<Json> demod(const cfvec& y, int max_err = 4);

}  // namespace si::rx::adsb::uat
