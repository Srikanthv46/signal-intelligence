// TPMS frames on a PCM symbol stream (research, signal gaps S6; python: receivers/ism_devices/tpms.py).
#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include "core/types.hpp"

namespace si::rx::ism_devices::tpms {

uint8_t crc8(const std::vector<uint8_t>& data, uint8_t poly, uint8_t init);
std::vector<int> manchester(const std::vector<int>& s, size_t start, size_t max_bits);
std::optional<Json> vdo(const std::vector<int>& s);
std::optional<Json> toyota(const std::vector<int>& s);
// python decode(s): the first CRC-valid frame, either polarity
std::optional<Json> decode(const std::vector<int>& s);

}  // namespace si::rx::ism_devices::tpms
