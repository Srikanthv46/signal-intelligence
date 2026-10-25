// IEEE 802.11 MAC frame dissector (python receivers/wifi/mac80211.py), shared by the DSSS and
// OFDM Wi-Fi receivers.
//
//   parse(mpdu_without_fcs) -> Json      flat fields + nested element dicts + a one-line "summary",
//                                       keys and values exactly as the Python dict
//   parse_elements(tlv)                  tagged elements (SSID, rates, RSN, HT/VHT/HE, WPS, P2P,
//                                       vendor elements incl. Open Drone ID / DJI DroneID)
//   vendor(mac) / oui_name(oui)          IEEE registry names (compiled-in table, oui_table.cpp)
//   security(info) / phy_generation(info) / summary(info)
//   crc32 / fcs_ok                       802.11 FCS (CRC-32, zlib.crc32)
#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "core/types.hpp"

namespace si::rx::wifi::mac80211 {

// zlib.crc32(data, crc)
uint32_t crc32(const uint8_t* data, size_t n, uint32_t crc = 0);
// python dsss.fcs_ok: len > 4 and crc32(mpdu[:-4]) == little-endian mpdu[-4:]
bool fcs_ok(const std::vector<uint8_t>& mpdu_with_fcs);

// Dissect one MPDU (without FCS).
Json parse(const std::vector<uint8_t>& mpdu);
Json parse(const uint8_t* mpdu, size_t n);
// Tagged elements -> dict (unknown IDs listed in element_ids).
Json parse_elements(const uint8_t* tlv, size_t n);

// Registry vendor of a unicast, universally administered address ("aa:bb:cc:dd:ee:ff");
// "broadcast" / "multicast" / "locally-administered (randomized/private)" / "unknown OUI".
std::optional<std::string> vendor(const std::optional<std::string>& mac);
// Name of a 3-byte organisation identifier ("00:50:f2") of a vendor-specific element.
std::string oui_name(const std::string& oui);

std::optional<std::string> security(const Json& info);
std::optional<std::string> phy_generation(const Json& info);
std::string summary(const Json& info);

}  // namespace si::rx::wifi::mac80211
