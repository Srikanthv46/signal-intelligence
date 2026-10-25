// COFDM / DVB-T receiver (research capability) vs python on synthetic DVB-T vectors (tests/vectors/cofdm; no capture).
#include <doctest/doctest.h>

#include <cmath>
#include <fstream>
#include <string>

#include "receivers/cofdm/detector.hpp"

using namespace si;

namespace {
const Json& ref() {
  static const Json j = [] {
    std::ifstream f(std::string(SI_TESTVEC) + "/cofdm/ref.json");
    REQUIRE(f);
    return Json::parse(f);
  }();
  return j;
}
void same(const Json& a, const Json& b, const std::string& path) {
  INFO(path);
  if (a.is_number() && b.is_number()) {
    CHECK(std::fabs(a.get<double>() - b.get<double>()) <= 1e-6 * std::fabs(b.get<double>()) + 0.011);
  } else if (a.is_object() && b.is_object()) {
    CHECK(a.size() == b.size());
    for (auto it = b.begin(); it != b.end(); ++it) {
      REQUIRE(a.contains(it.key()));
      same(a[it.key()], it.value(), path + "." + it.key());
    }
  } else {
    CHECK(a == b);
  }
}
}  // namespace

TEST_CASE("cofdm: synthetic DVB-T detected and TPS decoded == python") {
  for (const auto& c : ref()["cases"]) {
    std::ifstream f(std::string(SI_TESTVEC) + "/cofdm/" + c["file"].get<std::string>(), std::ios::binary | std::ios::ate);
    REQUIRE(f);
    const size_t n = static_cast<size_t>(f.tellg()) / sizeof(cf);
    cfvec x(n);
    f.seekg(0);
    f.read(reinterpret_cast<char*>(x.data()), static_cast<std::streamsize>(n * sizeof(cf)));
    const auto em = rx::cofdm::detect(x, ref()["fs"], ref()["fc"]);
    const Json& want = c["emissions"];
    REQUIRE(em.size() == want.size());
    for (size_t i = 0; i < em.size(); ++i) {
      CHECK(em[i].label == want[i]["label"].get<std::string>());
      CHECK(em[i].f_lo_mhz == doctest::Approx(want[i]["f_lo_mhz"].get<double>()).epsilon(1e-9));
      CHECK(em[i].f_hi_mhz == doctest::Approx(want[i]["f_hi_mhz"].get<double>()).epsilon(1e-9));
      same(em[i].features, want[i]["features"], "features");
    }
  }
}

TEST_CASE("cofdm: TPS bit decoder == python") {
  std::vector<int> bits;
  for (const auto& b : ref()["tps_bits"]) bits.push_back(b.get<int>());
  CHECK(rx::cofdm::tps_decode(bits).dump() == ref()["tps_out"].dump());
}
