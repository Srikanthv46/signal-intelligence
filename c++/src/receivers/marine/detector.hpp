// AIS receiver (research, signal gaps S15; python: receivers/marine/detector.py).
#pragma once

#include <vector>

#include "core/dsp.hpp"
#include "core/model.hpp"
#include "core/types.hpp"

namespace si::rx::marine {

// python channel_frames(ch, f_meas)
std::vector<Json> channel_frames(const Channelizer& ch, double f_meas);
// python detect(ch, ppm)
std::vector<Emission> detect(const Channelizer& ch, double ppm = 0.0);

}  // namespace si::rx::marine
