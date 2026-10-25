// DVB-T2 P1 preamble (cofdm/dvbt2_p1): research capability vs python on synthetic vectors
// (tests/vectors/dvbt2_p1/gen_dvbt2_p1.py), incl. cofdm detect on a 15.36 MSps capture.
#include <doctest/doctest.h>

#include <cmath>
#include <fstream>
#include <string>

#include "receivers/cofdm/detector.hpp"
#include "receivers/cofdm/dvbt2_p1.hpp"

using namespace si;

namespace {
cfvec load(const std::string& n) {
  std::ifstream f(std::string(SI_TESTVEC) + "/dvbt2_p1/" + n, std::ios::binary | std::ios::ate);
  REQUIRE(f);
  const size_t m = static_cast<size_t>(f.tellg()) / sizeof(cf);
  cfvec y(m);
  f.seekg(0);
  f.read(reinterpret_cast<char*>(y.data()), static_cast<std::streamsize>(m * sizeof(cf)));
  return y;
}
}  // namespace

TEST_CASE("dvbt2: P1 find / decode and cofdm detect == python") {
  std::ifstream jf(std::string(SI_TESTVEC) + "/dvbt2_p1/cases.json");
  REQUIRE(jf);
  const Json ref = Json::parse(jf);
  const cfvec all = load("p1.cf32");
  for (const auto& c : ref["cases"]) {
    const size_t o = c["offset"], n = c["n"];
    const cdvec y(all.begin() + static_cast<ptrdiff_t>(o), all.begin() + static_cast<ptrdiff_t>(o + n));
    const auto f = rx::cofdm::p1::find(y, 0.35);
    REQUIRE(f.size() == c["found"].size());
    for (size_t i = 0; i < f.size(); ++i) {
      CHECK(f[i].start == c["found"][i][0].get<int64_t>());
      CHECK(f[i].metric == doctest::Approx(c["found"][i][1].get<double>()).epsilon(1e-9));
      CHECK(f[i].decoded.dump() == c["found"][i][2].dump());
    }
  }
  const cfvec x = load("capture.cf32");
  const auto em = rx::cofdm::detect(x, ref["capture"]["fs"], ref["capture"]["fc"]);
  const Json& want = ref["capture"]["emissions"];
  REQUIRE(em.size() == want.size());
  for (size_t i = 0; i < em.size(); ++i) {
    CHECK(em[i].label == want[i]["label"].get<std::string>());
    CHECK(em[i].center_mhz == doctest::Approx(want[i]["center"].get<double>()));
    CHECK(em[i].features.dump() == want[i]["features"].dump());
  }
}
