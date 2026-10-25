// Remote ID over Wi-Fi NAN service discovery (wifi droneid::from_nan_action via mac80211::parse; research, PLAN T3)
// vs python (tests/vectors/wifi_nan/gen_wifi_nan.py).
#include <doctest/doctest.h>

#include <fstream>
#include <string>

#include "receivers/wifi/mac80211.hpp"

using namespace si;

TEST_CASE("wifi: Remote ID over NAN service discovery == python") {
  std::ifstream jf(std::string(SI_TESTVEC) + "/wifi_nan/nan.json");
  REQUIRE(jf);
  const Json ref = Json::parse(jf);
  const std::string h = ref["hex"];
  std::vector<uint8_t> b;
  for (size_t i = 0; i + 1 < h.size(); i += 2) b.push_back(static_cast<uint8_t>(std::stoi(h.substr(i, 2), nullptr, 16)));
  const Json r = rx::wifi::mac80211::parse(b);
  REQUIRE(r.contains("remote_id"));
  CHECK(r["remote_id"].dump() == ref["remote_id"].dump());
  CHECK(r["drone"].dump() == ref["drone"].dump());
}
