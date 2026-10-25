// ExpressLRS SYNC packets (rc_link/elrs, lora_phy::decode_implicit): research capability vs python on synthetic
// vectors (tests/vectors/elrs_sync/gen_elrs_sync.py), incl. the rc_link pipeline on a 915 MHz capture.
#include <doctest/doctest.h>

#include <fstream>
#include <string>

#include "receivers/rc_link/detector.hpp"
#include "receivers/rc_link/elrs.hpp"
#include "receivers/rc_link/lora_phy.hpp"

using namespace si;

namespace {
std::vector<uint8_t> unhex(const std::string& h) {
  std::vector<uint8_t> b;
  for (size_t i = 0; i + 1 < h.size(); i += 2) b.push_back(static_cast<uint8_t>(std::stoi(h.substr(i, 2), nullptr, 16)));
  return b;
}
}  // namespace

TEST_CASE("elrs: SYNC parse, implicit-header decode and the rc_link pipeline == python") {
  std::ifstream jf(std::string(SI_TESTVEC) + "/elrs_sync/cases.json");
  REQUIRE(jf);
  const Json ref = Json::parse(jf);
  for (const auto& c : ref["parse"]) CHECK(rx::rc_link::elrs::parse(unhex(c["hex"])).dump() == c["out"].dump());
  for (const auto& c : ref["implicit"]) {
    const auto d = rx::rc_link::lora_phy::decode_implicit(c["syms"].get<std::vector<int64_t>>(), c["sf"], c["cr"], 8);
    REQUIRE(d);
    CHECK(unhex(c["payload"]) == d->payload);
    CHECK(d->fec_corrected == c["corr"].get<int>());
    CHECK(d->fec_detected_errors == c["bad"].get<int>());
    CHECK(d->symbols_used == c["used"].get<int64_t>());
  }
  std::ifstream f(std::string(SI_TESTVEC) + "/elrs_sync/capture.cf32", std::ios::binary | std::ios::ate);
  REQUIRE(f);
  const size_t n = static_cast<size_t>(f.tellg()) / sizeof(cf);
  cfvec x(n);
  f.seekg(0);
  f.read(reinterpret_cast<char*>(x.data()), static_cast<std::streamsize>(n * sizeof(cf)));
  Json got = Json::array();
  for (const auto& e : rx::rc_link::detect(x, ref["capture"]["fs"], ref["capture"]["fc"]))
    if (e.decoded.is_object() && e.decoded.value("protocol", std::string()) == "ExpressLRS") got.push_back(e.decoded);
  CHECK(got.dump() == ref["capture"]["decoded"].dump());
}
