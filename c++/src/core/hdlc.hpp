// HDLC framing, AX.25 UI frames and APRS fields (research, signal gaps S19; python: core/hdlc.py).
#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "core/types.hpp"

namespace si::hdlc {

uint16_t crc_x25(const std::vector<uint8_t>& octets);
std::vector<int> unnrzi(const std::vector<int>& levels);
// python deframe(bits, min_octets, max_bits)
std::vector<std::pair<int64_t, std::vector<uint8_t>>> deframe(const std::vector<int>& bits, size_t min_octets = 18,
                                                               size_t max_bits = 3000);

struct Ax25 {
  std::string dest, src;
  std::vector<std::string> path;
  std::string info;   // latin-1 bytes
};
std::optional<Ax25> ax25(const std::vector<uint8_t>& octs);
// python aprs({"dest", "info"}): APRS fields
Json aprs(const std::string& dest, const std::string& info);
// latin-1 bytes -> UTF-8 (python bytes.decode("latin-1") as a JSON string)
std::string latin1_utf8(const std::string& s);

}  // namespace si::hdlc
