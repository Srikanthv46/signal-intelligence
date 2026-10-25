// Drone identification broadcast inside 802.11 frames (python receivers/wifi/droneid.py).
//
//   parse_odid(body)   ASTM F3411 / ASD-STAN prEN 4709-002 Open Drone ID (OUI fa:0b:bc, type 0x0d)
//   parse_dji(body)    DJI's 91-byte DroneID payload, located by its CRC-16/Kermit
//   from_vendor_ie     vendor-specific element -> Remote ID
//   signature          platform guess from the IEEE vendor name and/or the SSID ("-like" upstream)
//   summary            one operator line for a Remote ID dict
//
// Results are Json objects with the Python dicts' keys, in Python's order; std::nullopt
// where Python returns None.
#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "core/types.hpp"

namespace si::rx::wifi::droneid {

constexpr uint16_t kDjiCrcInit = 0x3692;
constexpr size_t kDjiLen = 91;

uint16_t crc16_kermit(const uint8_t* data, size_t n, uint16_t crc = kDjiCrcInit);

std::optional<Json> parse_odid(const std::vector<uint8_t>& body);
std::optional<Json> parse_dji(const std::vector<uint8_t>& body);
// `body` is everything after OUI + vendor type
std::optional<Json> from_vendor_ie(const std::string& oui, int vendor_type, const std::vector<uint8_t>& body);
// python from_nan_action(b): Remote ID in a Wi-Fi NAN service discovery frame (action body from the category byte)
std::optional<Json> from_nan_action(const std::vector<uint8_t>& b);
// walk raw elements (tag, len, value) for a Remote ID vendor element
std::optional<Json> scan_vendor_ies(const std::vector<uint8_t>& tlv);
// vendor_name / ssid: nullopt = None
std::optional<Json> signature(const std::optional<std::string>& vendor_name, const std::optional<std::string>& ssid);
std::string summary(const Json& rid);

}  // namespace si::rx::wifi::droneid
