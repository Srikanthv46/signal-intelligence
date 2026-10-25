// ILS glide slope / localizer and VOR (research, signal gaps S9 / S21; python: receivers/navaids/detector.py).
#pragma once

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "core/dsp.hpp"
#include "core/model.hpp"
#include "core/types.hpp"

namespace si::rx::navaids {

struct Channel {
  std::string kind;              // ils_gs | ils_loc | vor
  double mhz;
  std::optional<double> paired;  // the other half of an ILS pair
};
// python channels()
std::vector<Channel> channels();

// python detect(ch, ppm)
std::vector<Emission> detect(const Channelizer& ch, double ppm = 0.0);

}  // namespace si::rx::navaids
