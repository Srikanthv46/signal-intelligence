// rc_link/models (research, PLAN M5) vs python (tests/vectors/rc_models/gen_rc_models.py).
#include <doctest/doctest.h>

#include <fstream>
#include <string>

#include "receivers/rc_link/models.hpp"

using namespace si;

TEST_CASE("rc_link: model-level rules == python") {
  std::ifstream jf(std::string(SI_TESTVEC) + "/rc_models/cases.json");
  REQUIRE(jf);
  for (const auto& c : Json::parse(jf)) {
    std::vector<Emission> em;
    for (const auto& x : c["em"]) {
      Emission e;
      e.family = "rc_fhss";
      e.label = x["label"].get<std::string>();
      e.center_mhz = x["center_mhz"].get<double>();
      e.decoded = x["decoded"];
      e.features = x["features"];
      em.push_back(e);
    }
    CHECK(rx::rc_link::models::model(em).dump() == c["out"].dump());
  }
}
