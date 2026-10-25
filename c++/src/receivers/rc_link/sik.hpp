// SiK telemetry radio packets (python receivers/rc_link/sik.py; research, synthetic validation only).
#pragma once

#include <cstdint>
#include <vector>

#include "core/types.hpp"

namespace si::rx::rc_link::sik {

extern const std::vector<double> kRates;  // SiK air rates (b/s)
// one burst at ~8 samples/bit centred on the FSK carrier -> decoded dict, PHY-only dict, or null
Json sik_try_decode(const cfvec& y, double ofs, double rate);
// helpers exposed for tests
Json decode_packet(const std::vector<uint8_t>& bits_after_sync, bool inverted);
std::vector<std::pair<int64_t, bool>> find_sync(const std::vector<uint8_t>& bits);
// python demod_bits(y, ofs, rate): hard GFSK decisions at every symbol phase [(phase, bits)]
std::vector<std::pair<double, std::vector<uint8_t>>> demod_bits(const cfvec& y, double ofs, double rate);
// python two_level_threshold(v): 1-D 2-means decision level (independent of the ones / zeros ratio)
double two_level_threshold(const std::vector<double>& v, int iters = 6);
// python demod_bits_mid(y, ofs, rate): as demod_bits with the two-level threshold
std::vector<std::pair<double, std::vector<uint8_t>>> demod_bits_mid(const cfvec& y, double ofs, double rate);
// python WHITEN["pn9_msb_first"] (260 bytes)
const std::vector<uint8_t>& whiten_pn9_msb();

}  // namespace si::rx::rc_link::sik
