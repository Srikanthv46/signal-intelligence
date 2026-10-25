// ExpressLRS FHSS (rc_link/elrs; research, PLAN R2d) vs python (tests/vectors/elrs_fhss/gen_elrs_fhss.py).
#include <doctest/doctest.h>

#include <fstream>
#include <string>

#include "receivers/rc_link/detector.hpp"
#include "receivers/rc_link/elrs.hpp"

using namespace si;
namespace E = si::rx::rc_link::elrs;

TEST_CASE("elrs: FHSS sequence, UID2/UID3 search and hop recovery == python") {
  std::ifstream jf(std::string(SI_TESTVEC) + "/elrs_fhss/cases.json");
  REQUIRE(jf);
  const Json ref = Json::parse(jf);
  for (const auto& c : ref["seqs"])
    CHECK(Json(E::fhss_sequence(c[0].get<uint32_t>(), c[1].get<int>())).dump() == c[2].dump());
  for (const auto& c : ref["searches"]) {
    const Json& a = c[0];
    std::vector<std::pair<int64_t, int>> obs;
    for (const auto& o : a[3]) obs.push_back({o[0].get<int64_t>(), o[1].get<int>()});
    const auto got = E::uid23_candidates(a[0].get<int>(), a[1].get<int>(), a[2].get<int>(), obs, a[4].get<int>(),
                                         a[5].is_null() ? -1 : a[5].get<int>());
    Json gj = Json::array();
    for (const auto& [u2, u3] : got) gj.push_back(Json::array({u2, u3}));
    CHECK(gj.dump() == c[1].dump());
  }
  for (const auto& c : ref["recover"]) {
    const Json& a = c[0];
    std::vector<std::pair<double, int>> pk;
    for (const auto& p : a[6]) pk.push_back({p[0].get<double>(), p[1].get<int>()});
    CHECK(E::recover_hops(a[0], a[1].get<double>(), a[2].get<int>(), a[3].get<double>(), a[4].get<int>(), a[5].get<int>(), pk)
              .dump() == c[1].dump());
  }
}

TEST_CASE("elrs: hop recovery through the rc_link pipeline == python") {
  std::ifstream jf(std::string(SI_TESTVEC) + "/elrs_fhss/cases.json");
  REQUIRE(jf);
  const Json ref = Json::parse(jf)["pipeline"];
  std::ifstream f(std::string(SI_TESTVEC) + "/elrs_fhss/capture.cf32", std::ios::binary | std::ios::ate);
  REQUIRE(f);
  const size_t n = static_cast<size_t>(f.tellg()) / sizeof(cf);
  cfvec x(n);
  f.seekg(0);
  f.read(reinterpret_cast<char*>(x.data()), static_cast<std::streamsize>(n * sizeof(cf)));
  Json got = Json::array();
  for (const auto& e : rx::rc_link::detect(x, ref["fs"], ref["fc"]))
    if (e.decoded.is_object() && e.decoded.value("packet_type", Json(nullptr)) == "SYNC") got.push_back(e.decoded);
  CHECK(got.dump() == ref["sync"].dump());
}
