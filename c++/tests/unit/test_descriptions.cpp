// Plain-language signal descriptions (schema/descriptions.json) against python api.contract.description
// (tests/vectors/descriptions/gen_descriptions.py: every label of the latest regression + one per variant).
#include <doctest/doctest.h>

#include <fstream>
#include <string>

#include "api/contract.hpp"

using namespace si;

TEST_CASE("descriptions: every regression label and variant == python") {
  std::ifstream f(std::string(SI_TESTVEC) + "/descriptions/cases.json");
  REQUIRE(f.good());
  const Json cases = Json::parse(f);
  REQUIRE(cases.size() > 50);
  for (const auto& c : cases) {
    CAPTURE(c["label"].get<std::string>());
    const std::string dec = c["decoder"].is_string() ? c["decoder"].get<std::string>() : "";
    CHECK(contract::description(c["label"].get<std::string>(), dec, c["flags"]) == c["text"].get<std::string>());
  }
}
