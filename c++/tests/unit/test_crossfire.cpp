// TBS Crossfire 150 Hz (rc_link/crossfire): research capability vs python on a synthetic capture
// (tests/vectors/crossfire/gen_crossfire.py).
#include <doctest/doctest.h>

#include <fstream>
#include <string>

#include "receivers/rc_link/crossfire.hpp"
#include "receivers/rc_link/detector.hpp"

using namespace si;

TEST_CASE("crossfire: grid channels and the rc_link pipeline == python") {
  std::ifstream jf(std::string(SI_TESTVEC) + "/crossfire/cases.json");
  REQUIRE(jf);
  const Json ref = Json::parse(jf);
  for (const auto& g : ref["grid"]) CHECK(rx::rc_link::crossfire::grid_channel(g[0].get<double>()).dump() == g[1].dump());
  std::ifstream f(std::string(SI_TESTVEC) + "/crossfire/capture.cf32", std::ios::binary | std::ios::ate);
  REQUIRE(f);
  const size_t n = static_cast<size_t>(f.tellg()) / sizeof(cf);
  cfvec x(n);
  f.seekg(0);
  f.read(reinterpret_cast<char*>(x.data()), static_cast<std::streamsize>(n * sizeof(cf)));
  Json got = Json::array();
  for (const auto& e : rx::rc_link::detect(x, ref["fs"], ref["fc"]))
    if (e.label == "crossfire_sync") got.push_back(e.decoded);
  CHECK(got.dump() == ref["decoded"].dump());
}
