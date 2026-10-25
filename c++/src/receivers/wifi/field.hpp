// Drone Wi-Fi links among the 802.11 traffic of one signal (python receivers/wifi/field.py; research, PLAN W1).
#pragma once

#include <vector>

#include "core/types.hpp"

namespace si::rx::wifi::field {

// frames: python field.py frames {t0, t1, dur_us, bytes (or null), ta (or null), macs [], ctrl, video}
// -> one profile per BSS (AP, clients, airtime / PPDUs / bytes per direction, drone_link verdict), heaviest first
Json links(const std::vector<Json>& frames, double block_s);

}  // namespace si::rx::wifi::field
