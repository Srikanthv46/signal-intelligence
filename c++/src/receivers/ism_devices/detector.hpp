// 315 / 433 MHz short-range devices (research, signal gaps S6; python: receivers/ism_devices/detector.py).
#pragma once

#include <tuple>
#include <vector>

#include "core/dsp.hpp"
#include "core/model.hpp"
#include "core/types.hpp"

namespace si::rx::ism_devices {

// python BANDS
extern const std::vector<std::pair<double, double>> kBands;

struct Candidate {
  double f;
  int64_t frames;
  double peak_db;
};
// python candidates(ch)
std::vector<Candidate> candidates(const Channelizer& ch);

// python detect(ch): families "ism_devices", labels tpms | ook_fixed_code | keeloq | ook_device | fsk_device
std::vector<Emission> detect(const Channelizer& ch);

}  // namespace si::rx::ism_devices
