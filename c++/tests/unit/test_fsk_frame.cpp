// Blind FSK frame structure (rc_link/fsk_frame): research capability vs python on a synthetic capture
// (tests/vectors/fsk_frame/gen_fsk_frame.py).
#include <doctest/doctest.h>

#include <fstream>
#include <string>

#include "receivers/rc_link/detector.hpp"

using namespace si;

TEST_CASE("fsk_frame: blind frame structure through the rc_link pipeline == python") {
  std::ifstream jf(std::string(SI_TESTVEC) + "/fsk_frame/cases.json");
  REQUIRE(jf);
  const Json ref = Json::parse(jf);
  std::ifstream f(std::string(SI_TESTVEC) + "/fsk_frame/capture.cf32", std::ios::binary | std::ios::ate);
  REQUIRE(f);
  const size_t n = static_cast<size_t>(f.tellg()) / sizeof(cf);
  cfvec x(n);
  f.seekg(0);
  f.read(reinterpret_cast<char*>(x.data()), static_cast<std::streamsize>(n * sizeof(cf)));
  Json labels = Json::array(), frames = Json::array();
  for (const auto& e : rx::rc_link::detect(x, ref["fs"], ref["fc"]))
    if (e.features.contains("frame")) labels.push_back(e.label), frames.push_back(e.features["frame"]);
  CHECK(labels.dump() == ref["labels"].dump());
  CHECK(frames.dump() == ref["frames"].dump());
}
