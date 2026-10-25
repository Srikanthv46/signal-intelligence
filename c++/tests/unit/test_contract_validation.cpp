// schema/capabilities.json: proven / research validation tag (python api/contract.py validation()).
#include <doctest/doctest.h>

#include "api/contract.hpp"

using namespace si;

TEST_CASE("contract: validation tag (proven / research)") {
  const Json a = contract::validation("analog_fpv_video", "analog_vtx", Json{{"standard", "pal"}});
  CHECK(a["level"] == "proven");
  CHECK(a.contains("proof"));
  CHECK(!a.contains("research_fields"));
  const Json d = contract::validation("ofdm_data_link", "digital_vtx", Json{{"tx_antennas", {{"verdict", "2tx_stbc_like"}}}});
  CHECK(d["level"] == "proven");
  CHECK(d["research_fields"] == Json::array({"tx_antennas"}));
  CHECK(d["research_sources"] == Json::array({"bench/tx_antennas/FINDINGS.md"}));
  // a field capability of another receiver does not apply
  CHECK(!contract::validation("ofdm_data_link", "dji", Json{{"tx_antennas", 2}}).contains("research_fields"));
  // a technology missing from the registry is research-level
  CHECK(contract::validation("newthing", "", Json::object())["level"] == "research");
  // every technology code of technologies.json has a registry entry
  for (const auto& r : contract::technologies().at("rules"))
    CHECK(contract::capabilities().at("technologies").contains(r.at("technology").get<std::string>()));
}
