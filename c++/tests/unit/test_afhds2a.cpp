// FlySky AFHDS 2A (rc_link/afhds2a): research capability vs python on a synthetic capture
// (tests/vectors/afhds2a/gen_afhds2a.py).
#include <doctest/doctest.h>

#include <fstream>
#include <string>

#include "receivers/rc_link/afhds2a.hpp"
#include "receivers/rc_link/detector.hpp"

using namespace si;

TEST_CASE("afhds2a: CRCs and the rc_link pipeline == python") {
  std::ifstream jf(std::string(SI_TESTVEC) + "/afhds2a/cases.json");
  REQUIRE(jf);
  const Json ref = Json::parse(jf);
  for (const auto& c : ref["crc"]) {
    const std::string h = c[0].get<std::string>();
    std::vector<uint8_t> d;
    for (size_t i = 0; i + 1 < h.size(); i += 2) d.push_back(static_cast<uint8_t>(std::stoul(h.substr(i, 2), nullptr, 16)));
    CHECK(rx::rc_link::afhds2a::crc16(d.data(), d.size(), c[1].get<unsigned>(), c[2].get<unsigned>()) == c[3].get<int>());
  }
  std::ifstream f(std::string(SI_TESTVEC) + "/afhds2a/capture.cf32", std::ios::binary | std::ios::ate);
  REQUIRE(f);
  const size_t n = static_cast<size_t>(f.tellg()) / sizeof(cf);
  cfvec x(n);
  f.seekg(0);
  f.read(reinterpret_cast<char*>(x.data()), static_cast<std::streamsize>(n * sizeof(cf)));
  Json labels = Json::array(), got = Json::array();
  for (const auto& e : rx::rc_link::detect(x, ref["fs"], ref["fc"]))
    if (e.label.rfind("afhds2a", 0) == 0) labels.push_back(e.label), got.push_back(e.decoded);
  CHECK(labels.dump() == ref["labels"].dump());
  CHECK(got.dump() == ref["decoded"].dump());
}
