// Drone data links found inside other receivers' frames (PLAN G1 / G2 / G7) against python
// (tests/vectors/drone_links/gen_drone_links.py): MAVLink over Wi-Fi (UDP / TCP) and over transparent LoRa.
#include <doctest/doctest.h>

#include <fstream>
#include <string>
#include <vector>

#include "receivers/rc_link/lora_phy.hpp"
#include "receivers/wifi/mac80211.hpp"

using namespace si;

namespace {
const Json& ref() {
  static const Json r = [] {
    std::ifstream f(std::string(SI_TESTVEC) + "/drone_links/ref.json");
    REQUIRE(f.good());
    return Json::parse(f);
  }();
  return r;
}
std::vector<uint8_t> unhex(const std::string& h) {
  std::vector<uint8_t> v;
  for (size_t i = 0; i + 1 < h.size(); i += 2) v.push_back(static_cast<uint8_t>(std::stoi(h.substr(i, 2), nullptr, 16)));
  return v;
}
}  // namespace

TEST_CASE("drone links: MAVLink over Wi-Fi == python") {
  for (const auto& [k, c] : ref()["mpdu"].items()) {
    CAPTURE(k);
    const Json got = rx::wifi::mac80211::parse(unhex(c["hex"].get<std::string>()));
    for (const char* key : {"ipv4", "drone", "summary", "wfb_ng", "esp_now", "elrs_backpack", "action"}) {
      CAPTURE(key);
      CHECK((got.contains(key) ? got[key] : Json(nullptr)).dump() ==
            (c["out"].contains(key) ? c["out"][key] : Json(nullptr)).dump());
    }
  }
}

TEST_CASE("drone links: MAVLink over LoRa == python") {
  for (const auto& [k, c] : ref()["lora"].items()) {
    CAPTURE(k);
    CHECK(rx::rc_link::lora_phy::mavlink_in(unhex(c["hex"].get<std::string>())).dump() == c["out"].dump());
  }
}
