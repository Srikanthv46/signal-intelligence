// Item trackers and notable devices from a decoded BLE advertisement (research, signal gaps S7;
// python: core/ble_trackers.py).
#pragma once

#include <optional>

#include "core/types.hpp"

namespace si::ble_trackers {

// python recognize(ad): the strongest tracker / device match of one parsed advert (remote_id::parse_ble_ad output)
std::optional<Json> recognize(const Json& ad);

}  // namespace si::ble_trackers
