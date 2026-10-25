// FLARM (rc_link/flarm): research capability vs python (tests/vectors/flarm/gen_flarm.py).
#include <doctest/doctest.h>

#include <fstream>
#include <string>

#include "receivers/rc_link/detector.hpp"
#include "receivers/rc_link/flarm.hpp"

using namespace si;

TEST_CASE("flarm: clear header and the rc_link pipeline == python") {
  std::ifstream jf(std::string(SI_TESTVEC) + "/flarm/cases.json");
  REQUIRE(jf);
  const Json ref = Json::parse(jf);
  for (const auto& c : ref["standalone"]) {
    cfvec y;
    for (size_t i = 0; i < c["iq"][0].size(); ++i) y.push_back(cf(c["iq"][0][i].get<float>(), c["iq"][1][i].get<float>()));
    CHECK(rx::rc_link::flarm::try_identify(y, 8e5).dump() == c["out"].dump());
  }
  std::ifstream f(std::string(SI_TESTVEC) + "/flarm/capture.cf32", std::ios::binary | std::ios::ate);
  REQUIRE(f);
  const size_t n = static_cast<size_t>(f.tellg()) / sizeof(cf);
  cfvec x(n);
  f.seekg(0);
  f.read(reinterpret_cast<char*>(x.data()), static_cast<std::streamsize>(n * sizeof(cf)));
  Json labels = Json::array(), got = Json::array();
  for (const auto& e : rx::rc_link::detect(x, ref["fs"], ref["fc"]))
    if (e.label == "flarm_like") labels.push_back(e.label), got.push_back(e.decoded);
  CHECK(labels.dump() == ref["labels"].dump());
  CHECK(got.dump() == ref["decoded"].dump());
}
