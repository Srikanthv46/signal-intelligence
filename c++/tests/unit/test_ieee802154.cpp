// IEEE 802.15.4 O-QPSK (rc_link/ieee802154): research capability vs python on synthetic frames
// (tests/vectors/ieee802154/gen_ieee802154.py).
#include <doctest/doctest.h>

#include <fstream>
#include <string>

#include "receivers/rc_link/detector.hpp"
#include "receivers/rc_link/ieee802154.hpp"

using namespace si;
namespace Z = si::rx::rc_link::ieee802154;

namespace {
std::vector<uint8_t> unhex(const std::string& h) {
  std::vector<uint8_t> d;
  for (size_t i = 0; i + 1 < h.size(); i += 2) d.push_back(static_cast<uint8_t>(std::stoul(h.substr(i, 2), nullptr, 16)));
  return d;
}
}  // namespace

TEST_CASE("ieee802154: channels, CRC, MAC parse and demod == python") {
  std::ifstream jf(std::string(SI_TESTVEC) + "/ieee802154/cases.json");
  REQUIRE(jf);
  const Json ref = Json::parse(jf);
  for (const auto& c : ref["channels"]) {
    const auto k = Z::channel(c[0].get<double>());
    CHECK((c[1].is_null() ? !k.has_value() : (k && *k == c[1].get<int>())));
  }
  for (const auto& c : ref["crc"]) {
    const auto d = unhex(c[0].get<std::string>());
    CHECK(Z::crc16_kermit(d.data(), d.size()) == c[1].get<int>());
  }
  for (const auto& c : ref["mac"]) CHECK(Z::parse_mac(unhex(c[0].get<std::string>())).dump() == c[1].dump());
  std::ifstream f(std::string(SI_TESTVEC) + "/ieee802154/frames.cf32", std::ios::binary);
  REQUIRE(f);
  for (size_t i = 0; i < ref["lengths"].size(); ++i) {
    const size_t n = ref["lengths"][i].get<size_t>();
    std::vector<std::complex<float>> buf(n);
    f.read(reinterpret_cast<char*>(buf.data()), static_cast<std::streamsize>(n * sizeof(buf[0])));
    cfvec y(buf.begin(), buf.end());
    CHECK(Z::demod(y, 8e6).dump() == ref["out"][i].dump());
  }
}

TEST_CASE("ieee802154: the rc_link pipeline == python") {
  std::ifstream jf(std::string(SI_TESTVEC) + "/ieee802154/cases.json");
  REQUIRE(jf);
  const Json ref = Json::parse(jf);
  std::ifstream f(std::string(SI_TESTVEC) + "/ieee802154/capture.cf32", std::ios::binary | std::ios::ate);
  REQUIRE(f);
  const size_t n = static_cast<size_t>(f.tellg()) / sizeof(cf);
  cfvec x(n);
  f.seekg(0);
  f.read(reinterpret_cast<char*>(x.data()), static_cast<std::streamsize>(n * sizeof(cf)));
  Json labels = Json::array(), got = Json::array();
  for (const auto& e : rx::rc_link::detect(x, ref["fs"], ref["fc"]))
    if (e.label.rfind("ieee802154", 0) == 0) labels.push_back(e.label), got.push_back(e.decoded);
  CHECK(labels.dump() == ref["labels"].dump());
  CHECK(got.dump() == ref["decoded"].dump());
}
