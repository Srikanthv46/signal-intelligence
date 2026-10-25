// core/drone_models (research, PLAN M0-M2) vs python (tests/vectors/drone_models/gen_drone_models.py).
#include <doctest/doctest.h>

#include <fstream>
#include <string>

#include "core/drone_models.hpp"

using namespace si;

TEST_CASE("drone_models: product types, serials, model blocks == python") {
  std::ifstream jf(std::string(SI_TESTVEC) + "/drone_models/cases.json");
  REQUIRE(jf);
  const Json ref = Json::parse(jf);
  for (const auto& c : ref["dji"]) {
    const std::string n = drone_models::dji_product(c[0].get<int64_t>());
    CHECK((c[1].is_null() ? n.empty() : n == c[1].get<std::string>()));
    CHECK(drone_models::dji_product_source(c[0].get<int64_t>()) == c[2].get<std::string>());
  }
  for (const auto& c : ref["serial_models"]) CHECK(drone_models::serial_model(c[0].get<std::string>()).dump() == c[1].dump());
  for (const auto& c : ref["serials"]) CHECK(drone_models::serial_manufacturer(c[0].get<std::string>()).dump() == c[1].dump());
  for (const auto& c : ref["blocks"]) {
    const Json& a = c[0];
    const Json got = drone_models::model_block(a[0].get<std::string>(), a[1].get<std::vector<std::string>>(),
                                               a[2].get<std::vector<std::string>>(), a[3].get<std::string>(),
                                               a[4].is_null() ? "" : a[4].get<std::string>());
    CHECK(got.dump() == c[1].dump());
  }
  for (const auto& c : ref["rids"]) CHECK(drone_models::remote_id_model(c[0]).dump() == c[1].dump());
}
