// receivers/emitters against python receivers.emitters:
//  * tests/vectors/emitters/gen_assoc.py: extract_bursts -> measure_fingerprints -> associate on
//    the DJI RC-uplink capture testdata/D_20260916-171343 (learned reference-symbol timing);
//  * tests/vectors/emitters/gen_emitters.py: a synthetic two-radio capture (iq_ci16.bin) and a
//    hand-built multi-family emission list (every cue, veto, singleton, linked emitters).
#include <doctest/doctest.h>

#include <cmath>
#include <fstream>
#include <string>

#include "io/sigmf.hpp"
#include "receivers/emitters/association.hpp"
#include "receivers/emitters/bursts.hpp"
#include "receivers/emitters/families.hpp"

using namespace si;
using namespace si::rx::emitters;

namespace {

Json load(const std::string& rel) {
  std::ifstream f(std::string(SI_TESTVEC) + "/emitters/" + rel);
  REQUIRE(f.good());
  return Json::parse(f);
}

Emission from_json(const Json& j) {
  Emission e;
  e.family = j["family"];
  e.label = j["label"];
  e.start_us = j["start_us"];
  e.end_us = j["end_us"];
  e.f_lo_mhz = j.value("f_lo_mhz", 0.0);
  e.f_hi_mhz = j.value("f_hi_mhz", 0.0);
  e.center_mhz = j["center_mhz"];
  if (!j["snr_db"].is_null()) e.snr_db = j["snr_db"].get<double>();
  e.confidence = j.value("confidence", 0.5);
  if (j.contains("decoded")) e.decoded = j["decoded"];
  e.features = j["features"].is_null() ? Json::object() : j["features"];
  return e;
}

// python values vs C++: numbers to 1e-6 relative (1e-9 absolute), everything else exact
void same(const Json& py, const Json& cc, const std::string& path) {
  CAPTURE(path);
  if (py.is_number() && cc.is_number()) {
    const double a = py.get<double>(), b = cc.get<double>();
    CHECK(std::fabs(a - b) <= std::max(1e-6 * std::fabs(a), 1e-9));
  } else if (py.is_object() && cc.is_object()) {
    CHECK(py.size() == cc.size());
    for (const auto& [k, v] : py.items()) {
      const std::string key = path + "." + k;
      REQUIRE_MESSAGE(cc.contains(k), key);
      same(v, cc[k], path + "." + k);
    }
  } else if (py.is_array() && cc.is_array()) {
    REQUIRE(py.size() == cc.size());
    for (size_t i = 0; i < py.size(); ++i) same(py[i], cc[i], path + "[" + std::to_string(i) + "]");
  } else {
    CHECK(py == cc);
  }
}

void same_emitters(const Json& py, const std::vector<Emitter>& cc, const std::string& path) {
  REQUIRE(py.size() == cc.size());
  for (size_t i = 0; i < cc.size(); ++i) same(py[i], cc[i].to_dict(), path + "[" + std::to_string(i) + "]");
}

}  // namespace

TEST_CASE("emitters: bursts, fingerprints, association on D_20260916-171343") {
  const Json ref = load("assoc_D_20260916-171343.json");
  const Capture cap = load_sigmf(std::string(SI_TESTDATA) + "/D_20260916-171343/D_20260916-171343.sigmf-meta");
  auto em = extract_bursts(cap.iq(), cap.fs, cap.fc);
  REQUIRE(em.size() == ref["emissions"].size());
  measure_fingerprints(cap.iq(), cap.fs, cap.fc, em);
  for (size_t i = 0; i < em.size(); ++i) {
    const Json& r = ref["emissions"][i];
    CHECK(em[i].label == r["label"].get<std::string>());
    same(r["start_us"], em[i].start_us, "start_us");
    same(r["center_mhz"], em[i].center_mhz, "center_mhz");
    same(r["snr_db"], *em[i].snr_db, "snr_db");
    same(r["features"], em[i].features, "emission[" + std::to_string(i) + "].features");
  }
  same_emitters(ref["emitters"], associate(em), "emitters");
  const auto [llr, ev] = same_emitter_score(em[0], em[1]);
  same(ref["pair01"]["llr"], llr, "pair01.llr");
  Json evj = Json::array();
  for (const auto& c : ev) evj.push_back(c.to_json());
  same(ref["pair01"]["evidence"], evj, "pair01.evidence");
}

TEST_CASE("emitters: synthetic two-radio capture") {
  const Json ref = load("ref.json")["capture"];
  std::ifstream f(std::string(SI_TESTVEC) + "/emitters/iq_ci16.bin", std::ios::binary | std::ios::ate);
  REQUIRE(f.good());
  std::vector<int16_t> q(static_cast<size_t>(f.tellg()) / sizeof(int16_t));
  f.seekg(0);
  f.read(reinterpret_cast<char*>(q.data()), static_cast<std::streamsize>(q.size() * sizeof(int16_t)));
  cfvec x(q.size() / 2);
  for (size_t i = 0; i < x.size(); ++i) x[i] = cf(q[2 * i] / 2048.0f, q[2 * i + 1] / 2048.0f);
  const double fs = ref["fs"], fc = ref["fc"];
  std::vector<Emission> em;
  for (const auto& j : ref["emissions_in"]) em.push_back(from_json(j));
  measure_fingerprints(x, fs, fc, em);
  for (size_t i = 0; i < em.size(); ++i) {
    same(ref["emissions_out"][i]["features"], em[i].features, "emissions_out[" + std::to_string(i) + "].features");
    same(ref["emissions_out"][i]["snr_db"], *em[i].snr_db, "snr_db");
  }
  same_emitters(ref["emitters"], associate(em, &assoc_models()), "capture.emitters");
  const auto b = extract_bursts(x, fs, fc);
  REQUIRE(b.size() == ref["bursts"].size());
  for (size_t i = 0; i < b.size(); ++i) {
    CHECK(b[i].label == ref["bursts"][i]["label"].get<std::string>());
    same(ref["bursts"][i]["start_us"], b[i].start_us, "bursts.start_us");
    same(ref["bursts"][i]["f_lo_mhz"], b[i].f_lo_mhz, "bursts.f_lo_mhz");
    same(ref["bursts"][i]["features"], b[i].features, "bursts.features");
  }
}

TEST_CASE("emitters: association on a hand-built multi-family emission list") {
  const Json ref = load("ref.json")["assoc"];
  std::vector<Emission> S;
  for (const auto& j : ref["emissions"]) S.push_back(from_json(j));
  const auto A = associate_full(S, &assoc_models());
  for (size_t i = 0; i < S.size(); ++i)
    for (size_t j = 0; j < S.size(); ++j) same(ref["M"][i][j], A.M(i, j), "M");
  same_emitters(ref["emitters"], A.emitters, "assoc.emitters");
  for (const auto& p : ref["pairs"]) {
    const auto [llr, ev] = same_emitter_score(S[p[0].get<size_t>()], S[p[1].get<size_t>()], &assoc_models());
    Json evj = Json::array();
    for (const auto& c : ev) evj.push_back(c.to_json());
    same(p[2], evj, "pair " + p[0].dump() + "," + p[1].dump());
  }
  same_emitters(ref["emitters_default_models"], associate(S), "emitters_default_models");
}
