// Bell 202 AFSK packet radio (APRS / AX.25) on an FM voice channel (research, signal gaps S19;
// python: receivers/voice_fm/afsk.py).
#pragma once

#include <vector>

#include "core/types.hpp"

namespace si::rx::voice_fm::afsk {

// python decode(d, fs): d = FM discriminator output (Hz, float32) at fs
std::vector<Json> decode(const fvec& d, double fs);

}  // namespace si::rx::voice_fm::afsk
