// Sub-audible squelch signalling on an FM discriminator output: CTCSS tones and DCS codes
// (python: receivers/voice_fm/subaudio.py).
//
// DCS: a 23-bit Golay(23,12) codeword sent continuously as NRZ FSK at 134.4 bit/s; bits 0-8
// the octal code (LSB first), bits 9-11 the marker 100, bits 12-22 the parity of generator
// 0xC75. A 1 bit is taken as a positive frequency shift ('N').
#pragma once

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "core/types.hpp"

namespace si::rx::voice_fm {

constexpr double kDcsRate = 134.4;
constexpr uint32_t kGolayG = 0xC75;
constexpr uint32_t kMask23 = 0x7FFFFF;
constexpr double kSubFs = 2000.0;

uint32_t golay_syndrome(uint32_t w, uint32_t g = kGolayG);
uint32_t dcs_encode(const std::string& code_octal);
// every (code, polarity) whose transmitted bit stream is the same as word's, sorted
std::vector<std::pair<std::string, std::string>> dcs_aliases(uint32_t word);

dvec to_subaudio(const fvec& disc, double fs);                 // 0-300 Hz band at kSubFs
std::optional<Json> dcs(const dvec& a, int min_words = 2);
std::optional<Json> ctcss(const dvec& a, double min_s = 0.3);
std::optional<Json> subaudible(const fvec& disc, double fs);   // DCS first, then CTCSS

}  // namespace si::rx::voice_fm
