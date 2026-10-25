// Hitec Optima / Minima (CC2500 MSK) and Futaba S-FHSS (CC2500 2-FSK) 2.4 GHz RC links (python
// receivers/rc_link/hitec.py; research, synthetic only, PLAN G12).
#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include "core/types.hpp"

namespace si::rx::rc_link::hitec {

constexpr double kRate = 249.939e3, kDev = kRate / 4;
constexpr double kF0 = 2405.59985e6, kStep = 0.299927e6;
constexpr double kSfRate = 128.143e3, kSfDev = 38.086e3;
constexpr double kSfF0 = 2403.9987e6, kSfStep = 1.49963e6;
// python crc16_cc2400(data, init)
uint16_t crc16_cc2400(const uint8_t* p, size_t n, uint16_t init);
// python hitec_try_decode(y, ofs, f_hz=None, rate=RATE): the dict or null
Json hitec_try_decode(const cfvec& y, double ofs, std::optional<double> f_hz = std::nullopt, double rate = kRate);
// python sfhss_try_decode(y, ofs, f_hz=None, rate=SF_RATE): the dict or null
Json sfhss_try_decode(const cfvec& y, double ofs, std::optional<double> f_hz = std::nullopt, double rate = kSfRate);

}  // namespace si::rx::rc_link::hitec
