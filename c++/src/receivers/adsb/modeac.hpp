// SSR Mode A/C replies at 1090 MHz (research, signal gaps S11; python: receivers/adsb/modeac.py).
#pragma once

#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "core/types.hpp"

namespace si::rx::adsb::modeac {

constexpr double kFrameUs = 20.3;
constexpr double kStepUs = 1.45;
constexpr double kSpiUs = 24.65;

// pulse bits by position name (C1 A1 C2 A2 C4 A4 X B1 D1 B2 D2 B4 D4)
using Bits = std::map<std::string, int>;

// python gillham_ft(p): Mode C altitude in ft, or nullopt when not a valid Gillham code
std::optional<int> gillham_ft(const Bits& p);
// python squawk(p): Mode A identity ABCD
std::string squawk(const Bits& p);

// python decode(env, fs, busy, noise): replies in a magnitude envelope (noise < 0: median of env)
std::vector<Json> decode(const fvec& env, double fs, const std::vector<std::pair<double, double>>& busy,
                         double noise = -1.0);
// python summarise(replies)
Json summarise(const std::vector<Json>& replies);

}  // namespace si::rx::adsb::modeac
