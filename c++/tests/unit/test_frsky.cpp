// FrSky D8 / D16 (rc_link/cc2500): research capability vs python on a synthetic capture
// (tests/vectors/frsky/gen_frsky.py).
#include <doctest/doctest.h>

#include <fstream>
#include <string>

#include "receivers/rc_link/cc2500.hpp"
#include "receivers/rc_link/detector.hpp"

using namespace si;

TEST_CASE("frsky: CRCs and the rc_link pipeline == python") {
  std::ifstream jf(std::string(SI_TESTVEC) + "/frsky/cases.json");
  REQUIRE(jf);
  const Json ref = Json::parse(jf);
  for (const auto& c : ref["crc"]) {
    const std::string h = c[0].get<std::string>();
    std::vector<uint8_t> d;
    for (size_t i = 0; i + 1 < h.size(); i += 2) d.push_back(static_cast<uint8_t>(std::stoul(h.substr(i, 2), nullptr, 16)));
    CHECK(rx::rc_link::cc2500::crc16_cc2500(d.data(), d.size()) == c[1].get<int>());
    CHECK(rx::rc_link::cc2500::crc16_d16(d.data(), d.size()) == c[2].get<int>());
  }
  std::ifstream f(std::string(SI_TESTVEC) + "/frsky/capture.cf32", std::ios::binary | std::ios::ate);
  REQUIRE(f);
  const size_t n = static_cast<size_t>(f.tellg()) / sizeof(cf);
  cfvec x(n);
  f.seekg(0);
  f.read(reinterpret_cast<char*>(x.data()), static_cast<std::streamsize>(n * sizeof(cf)));
  Json labels = Json::array(), got = Json::array();
  for (const auto& e : rx::rc_link::detect(x, ref["fs"], ref["fc"]))
    if (e.label.rfind("frsky_", 0) == 0) labels.push_back(e.label), got.push_back(e.decoded);
  CHECK(labels.dump() == ref["labels"].dump());
  CHECK(got.dump() == ref["decoded"].dump());
}
