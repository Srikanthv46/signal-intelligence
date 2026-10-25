// XN297 / Bayang (rc_link/xn297): research capability vs python on a synthetic capture
// (tests/vectors/xn297/gen_xn297.py).
#include <doctest/doctest.h>

#include <fstream>
#include <string>

#include "receivers/rc_link/detector.hpp"
#include "receivers/rc_link/xn297.hpp"

using namespace si;

TEST_CASE("xn297: scramble table, CRC and the rc_link pipeline == python") {
  std::ifstream jf(std::string(SI_TESTVEC) + "/xn297/cases.json");
  REQUIRE(jf);
  const Json ref = Json::parse(jf);
  CHECK(Json(rx::rc_link::xn297::scramble()).dump() == ref["scramble"].dump());
  for (const auto& c : ref["crc"]) {
    const std::string h = c[0].get<std::string>();
    std::vector<uint8_t> d;
    for (size_t i = 0; i + 1 < h.size(); i += 2) d.push_back(static_cast<uint8_t>(std::stoul(h.substr(i, 2), nullptr, 16)));
    CHECK(rx::rc_link::xn297::crc16(d.data(), d.size()) == c[1].get<int>());
  }
  std::ifstream f(std::string(SI_TESTVEC) + "/xn297/capture.cf32", std::ios::binary | std::ios::ate);
  REQUIRE(f);
  const size_t n = static_cast<size_t>(f.tellg()) / sizeof(cf);
  cfvec x(n);
  f.seekg(0);
  f.read(reinterpret_cast<char*>(x.data()), static_cast<std::streamsize>(n * sizeof(cf)));
  Json labels = Json::array(), got = Json::array();
  for (const auto& e : rx::rc_link::detect(x, ref["fs"], ref["fc"]))
    if (e.label == "bayang" || e.label == "xn297_like") labels.push_back(e.label), got.push_back(e.decoded);
  for (const auto& c : ref["standalone"]) {
    cfvec y;
    for (size_t i = 0; i < c["iq"][0].size(); ++i) y.push_back(cf(c["iq"][0][i].get<float>(), c["iq"][1][i].get<float>()));
    CHECK(rx::rc_link::xn297::try_decode(y, 8e6).dump() == c["out"].dump());
  }
  CHECK(labels.dump() == ref["labels"].dump());
  CHECK(got.dump() == ref["decoded"].dump());
}
