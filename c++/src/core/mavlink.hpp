// MAVLink v1 / v2 frame scanning and message decoding (python core/mavlink.py; research:
// research/drone_transmitters/05_sik_mavlink_details.md). Json with the Python dicts' keys in Python's order.
#pragma once

#include <cstdint>
#include <vector>

#include "core/types.hpp"

namespace si::mavlink {

uint16_t x25(const uint8_t* p, size_t n, uint16_t crc = 0xFFFF);
// every frame in data (python scan)
Json scan(const std::vector<uint8_t>& data);
// latest identity / position / state over CRC-valid frames (python summary)
Json summary(const Json& frames);

}  // namespace si::mavlink
