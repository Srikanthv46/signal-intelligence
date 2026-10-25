// AIS codec (research, signal gaps S15; python: receivers/marine/ais.py).
#pragma once

#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

#include "core/types.hpp"

namespace si::rx::marine::ais {

uint16_t crc_x25(const std::vector<uint8_t>& octets);
bool residue_ok(const std::vector<uint8_t>& octets);
std::vector<int> unnrzi(const std::vector<int>& levels);
// python frames(bits): (opening flag index, data octets without the FCS) of the CRC-valid frames
std::vector<std::pair<int64_t, std::vector<uint8_t>>> frames(const std::vector<int>& bits);
// python fields(octs, nbits)
Json fields(const std::vector<uint8_t>& octs, std::optional<size_t> nbits = std::nullopt);

}  // namespace si::rx::marine::ais
