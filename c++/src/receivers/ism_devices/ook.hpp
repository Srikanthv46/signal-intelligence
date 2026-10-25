// OOK remotes on a pulse list (research, signal gaps S6; python: receivers/ism_devices/ook.py).
#pragma once

#include <optional>
#include <vector>

#include "core/types.hpp"

namespace si::rx::ism_devices::ook {

struct Pulse {
  double t_us, mark_us, space_us;   // space_us = inf after the last mark
};

// python fixed_code(p): EV1527 / PT2262
std::optional<Json> fixed_code(const std::vector<Pulse>& p);
// python keeloq(p): HCS200 / 300 / 301 clear fields (the encrypted part is not reported)
std::optional<Json> keeloq(const std::vector<Pulse>& p);

}  // namespace si::rx::ism_devices::ook
