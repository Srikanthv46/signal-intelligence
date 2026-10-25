// Rules over the device tracks (python: receivers/devices/threats.py):
//   findings(devices, frames, links, drone_links) -> findings, most serious first
// Severity: critical | high | medium | low | info.
#pragma once

#include <vector>

#include "core/types.hpp"
#include "receivers/devices/tracker.hpp"

namespace si::rx::devices {

// drone_links: the signals' wifi drone_links (PLAN W1); the behaviour-based ones become findings
Json findings(const Json& devices, const std::vector<Frame>& frames, const Json& links = Json::array(),
              const Json& drone_links = Json::array());

}  // namespace si::rx::devices
