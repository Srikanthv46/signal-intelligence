// Barrage-noise and noise-FM jammers (research, signal gaps S12; python: receivers/jammer/barrage.py).
#pragma once

#include <optional>
#include <string>
#include <vector>

#include "core/dsp.hpp"
#include "core/types.hpp"

namespace si::rx::jammer::barrage {

// python covers_band(f0, f1): allocations whose edges are within 2 MHz of the box's
std::vector<std::string> covers_band(double f0, double f1);
// python analyse(ch, f0, f1): measurements + "class" ("barrage" | "noise_fm" | null); nullopt below 1 MHz
std::optional<Json> analyse(const Channelizer& ch, double f0, double f1);

}  // namespace si::rx::jammer::barrage
