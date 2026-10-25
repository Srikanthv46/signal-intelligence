// FM broadcast carriers, stereo pilot and RDS (research, signal gaps S20; python: receivers/broadcast_fm/detector.py).
#pragma once

#include <utility>
#include <vector>

#include "core/dsp.hpp"
#include "core/model.hpp"
#include "core/types.hpp"

namespace si::rx::broadcast_fm {

// python candidates(ch): (channel MHz, dB over the band floor)
std::vector<std::pair<double, double>> candidates(const Channelizer& ch);
// python rds_bits(mpx, fs, f57): differentially decoded RDS bits
std::vector<int> rds_bits(const fvec& mpx, double fs, double f57);
// python detect(ch)
std::vector<Emission> detect(const Channelizer& ch);

}  // namespace si::rx::broadcast_fm
