// MAVLink (core/mavlink) and SiK telemetry (rc_link/sik): research capabilities vs python on synthetic vectors.
#include <doctest/doctest.h>

#include <fstream>
#include <string>

#include "core/mavlink.hpp"
#include "receivers/rc_link/sik.hpp"

using namespace si;

namespace {
Json load(const std::string& n) {
  std::ifstream f(std::string(SI_TESTVEC) + "/sik/" + n);
  REQUIRE(f);
  return Json::parse(f);
}
}  // namespace

TEST_CASE("mavlink: scan + summary == python") {
  const Json ref = load("mavlink.json");
  const std::string h = ref["hex"].get<std::string>();
  std::vector<uint8_t> b;
  for (size_t i = 0; i + 1 < h.size(); i += 2) b.push_back(static_cast<uint8_t>(std::stoi(h.substr(i, 2), nullptr, 16)));
  const Json fr = mavlink::scan(b);
  CHECK(fr.dump() == ref["scan"].dump());
  CHECK(mavlink::summary(fr).dump() == ref["summary"].dump());
}

TEST_CASE("sik: synthetic packet decoded == python") {
  const Json ref = load("sik.json");
  std::ifstream f(std::string(SI_TESTVEC) + "/sik/sik.cf32", std::ios::binary | std::ios::ate);
  REQUIRE(f);
  const size_t n = static_cast<size_t>(f.tellg()) / sizeof(cf);
  cfvec y(n);
  f.seekg(0);
  f.read(reinterpret_cast<char*>(y.data()), static_cast<std::streamsize>(n * sizeof(cf)));
  const Json got = rx::rc_link::sik::sik_try_decode(y, ref["fs"], ref["rate"]);
  REQUIRE(got.is_object());
  const Json& want = ref["out"];
  for (const char* k : {"netid", "length", "whitening", "hw_crc_variant", "mavlink_crc_ok", "data_hex", "polarity_inverted"})
    CHECK(got[k] == want[k]);
  CHECK(got["mavlink_frames"].dump() == want["mavlink_frames"].dump());
}
