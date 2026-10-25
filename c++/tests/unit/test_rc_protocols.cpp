// RC / telemetry protocol research catalog (rc_link/protocols) == python on fixed measurement tuples (no capture).
#include <doctest/doctest.h>

#include <fstream>
#include <string>

#include "receivers/rc_link/protocols.hpp"

using namespace si;

TEST_CASE("rc_link: protocol research catalog == python") {
  std::ifstream f(std::string(SI_TESTVEC) + "/rc_protocols/ref.json");
  REQUIRE(f);
  const Json ref = Json::parse(f);
  auto opt = [](const Json& v) { return v.is_null() ? std::optional<double>() : std::optional<double>(v.get<double>()); };
  for (const auto& c : ref) {
    const Json& a = c["args"];
    const Json got = rx::rc_link::protocols::candidates(a[0].get<double>(), opt(a[1]), opt(a[2]), opt(a[3]));
    CHECK(got.dump() == c["out"].dump());
  }
}
