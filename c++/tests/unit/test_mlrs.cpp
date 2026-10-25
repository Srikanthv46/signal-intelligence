// mLRS (rc_link/mlrs): research capability vs python (tests/vectors/mlrs/gen_mlrs.py).
#include <doctest/doctest.h>

#include <fstream>
#include <string>

#include "receivers/rc_link/detector.hpp"
#include "receivers/rc_link/mlrs.hpp"

using namespace si;

namespace {
std::vector<uint8_t> unhex(const std::string& h) {
  std::vector<uint8_t> d;
  for (size_t i = 0; i + 1 < h.size(); i += 2) d.push_back(static_cast<uint8_t>(std::stoul(h.substr(i, 2), nullptr, 16)));
  return d;
}
}  // namespace

TEST_CASE("mlrs: CRC, frame parse and the rc_link pipeline == python") {
  std::ifstream jf(std::string(SI_TESTVEC) + "/mlrs/cases.json");
  REQUIRE(jf);
  const Json ref = Json::parse(jf);
  for (const auto& c : ref["crc"]) {
    const auto d = unhex(c[0].get<std::string>());
    CHECK(rx::rc_link::mlrs::crc16_mcrf4xx(d.data(), d.size()) == c[1].get<int>());
  }
  for (const auto& c : ref["parse"]) CHECK(rx::rc_link::mlrs::parse(unhex(c[0].get<std::string>())).dump() == c[1].dump());
  for (const auto& cap : ref["captures"]) {
    std::ifstream f(std::string(SI_TESTVEC) + "/mlrs/" + cap["file"].get<std::string>(), std::ios::binary | std::ios::ate);
    REQUIRE(f);
    const size_t n = static_cast<size_t>(f.tellg()) / sizeof(cf);
    cfvec x(n);
    f.seekg(0);
    f.read(reinterpret_cast<char*>(x.data()), static_cast<std::streamsize>(n * sizeof(cf)));
    Json labels = Json::array(), got = Json::array(), modes = Json::array();
    for (const auto& e : rx::rc_link::detect(x, cap["fs"], cap["fc"]))
      if (e.label.rfind("mlrs", 0) == 0) labels.push_back(e.label), got.push_back(e.decoded), modes.push_back(e.features["mlrs_modes"]);
    CHECK(labels.dump() == cap["labels"].dump());
    CHECK(got.dump() == cap["decoded"].dump());
    CHECK(modes.dump() == cap["modes"].dump());
  }
}
