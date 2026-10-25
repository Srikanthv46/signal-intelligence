// ASTM F3411 / ASD-STAN prEN 4709-002 Open Drone ID (broadcast Remote ID) and Bluetooth LE advertising data
// (python core/remote_id.py). Json objects with the Python dicts' keys in Python's order; nullopt = None.
#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "core/types.hpp"

namespace si::remote_id {

// China GB 46750-2025 packet (from the 0xFF type byte; research, PLAN G14), nullopt when not well formed
std::optional<Json> parse_gb46750(const std::vector<uint8_t>& b, const std::string& transport = "Wi-Fi Beacon vendor element");
// [message counter][message pack | message]; transport names where it was carried
std::optional<Json> parse_odid(const std::vector<uint8_t>& body, const std::string& transport = "Wi-Fi Beacon vendor element");
// advertising data after AdvA -> AD fields (flags, local_name, company, service data, ..., remote_id)
Json parse_ble_ad(const std::vector<uint8_t>& adv_data);

}  // namespace si::remote_id
