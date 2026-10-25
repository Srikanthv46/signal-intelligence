// Mode S / ADS-B 1090 MHz: preamble correlation, PPM bits, CRC-24 (python: receivers/adsb/modes.py).
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "core/types.hpp"

namespace si::rx::adsb {

constexpr double kFs = 8e6;        // envelope rate
constexpr int kSpu = 8;            // samples per us
constexpr double kFAdsb = 1090e6;

// Mode S CRC-24 (generator 0xFFF409, Annex 10 Vol IV 3.1.2.3.3) of the data bits (0/1 values)
uint32_t crc24(const uint8_t* bits, size_t n);
// CRC of the data bits xor the transmitted parity field (the last 24 bits)
uint32_t parity_residual(const uint8_t* bits, size_t n);

// python demod(env, thr): frames that pass parity, as the Python dicts (t_s, df, bits, hex,
// parity_residual, snr_db, aa, confirmed). env = |slice| at kFs (float32), thr float32.
std::vector<Json> demod(const fvec& env, float thr);

}  // namespace si::rx::adsb
