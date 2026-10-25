// receivers/gsm against the Python receiver on synthetic FCCH + SCH bursts (tests/vectors/gsm/gen_gsm.py).
#include <doctest/doctest.h>

#include <cmath>
#include <fstream>

#include "receivers/gsm/sch.hpp"

using namespace si;

namespace {

const Json& vec() {
  static const Json j = [] {
    std::ifstream f(std::string(SI_TESTVEC) + "/gsm/gsm.json");
    REQUIRE_MESSAGE(f.good(), "missing tests/vectors/gsm/gsm.json (run gen_gsm.py)");
    return Json::parse(f);
  }();
  return j;
}

template <class C>
std::vector<C> cvec(const Json& a) {
  std::vector<C> v;
  for (const auto& p : a) v.emplace_back(p[0].get<double>(), p[1].get<double>());
  return v;
}

}  // namespace

TEST_CASE("gsm: FCCH search, SCH demodulation, Viterbi, CRC-10 and fields equal Python") {
  int decoded = 0;
  for (const auto& c : vec()["cases"]) {
    const cfvec y = cvec<cf>(c["y"]);
    const auto fc = rx::gsm::find_fcch(y);
    REQUIRE(fc.has_value() == c.contains("fcch"));
    if (!fc) continue;
    CHECK(fc->start == c["fcch"][0].get<int64_t>());
    CHECK(std::fabs(fc->cfo - c["fcch"][1].get<double>()) < 1e-6);
    CHECK(std::fabs(fc->purity - c["fcch"][2].get<double>()) < 1e-9);
    const cdvec y2 = cvec<cd>(c["y2"]);
    const auto d = rx::gsm::demod_sch(y2, fc->start + 5000);
    REQUIRE(d.has_value() == c.contains("soft"));
    if (!d) continue;
    CHECK(d->ts_match == c["ts_match"].get<int>());
    CHECK(d->sample == c["sample"].get<int64_t>());
    const dvec ref = c["soft"].get<dvec>();
    REQUIRE(d->soft.size() == ref.size());
    for (size_t i = 0; i < ref.size(); ++i) CHECK(std::fabs(d->soft[i] - ref[i]) < 1e-12);
    const auto bits = rx::gsm::viterbi(d->soft);
    CHECK(bits == c["bits"].get<std::vector<int>>());
    const std::vector<int> b25(bits.begin(), bits.begin() + 25), p10(bits.begin() + 25, bits.begin() + 35);
    const bool ok = rx::gsm::crc10_ok(b25, p10);
    CHECK(ok == c["crc_ok"].get<bool>());
    if (ok) {
      CHECK(rx::gsm::sch_fields(b25) == c["fields"]);
      ++decoded;
    }
  }
  CHECK(decoded > 0);
  CHECK(rx::gsm::viterbi(vec()["vit_in"].get<dvec>()) == vec()["vit_out"].get<std::vector<int>>());
  CHECK(rx::gsm::sch_fields(vec()["fields_in"].get<std::vector<int>>()) == vec()["fields"]);
}
