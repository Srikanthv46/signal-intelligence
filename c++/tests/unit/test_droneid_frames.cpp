// DJI DroneID 9-symbol vs 8-symbol (legacy) frames (dji detect_droneid; research, PLAN D4) vs python on synthetic
// 15.36 MSps captures (tests/vectors/droneid_frames/gen_droneid_frames.py).
#include <doctest/doctest.h>

#include <cmath>
#include <fstream>
#include <string>

#include "receivers/dji/detector.hpp"

using namespace si;

TEST_CASE("dji: DroneID 9-symbol / 8-symbol frames == python") {
  std::ifstream jf(std::string(SI_TESTVEC) + "/droneid_frames/cases.json");
  REQUIRE(jf);
  const Json ref = Json::parse(jf);
  for (const Json& c : ref["cases"]) {
    std::ifstream f(std::string(SI_TESTVEC) + "/droneid_frames/" + c["file"].get<std::string>(), std::ios::binary | std::ios::ate);
    REQUIRE(f);
    const size_t n = static_cast<size_t>(f.tellg()) / sizeof(cf);
    cfvec x(n);
    f.seekg(0);
    f.read(reinterpret_cast<char*>(x.data()), static_cast<std::streamsize>(n * sizeof(cf)));
    std::vector<Emission> em;
    for (auto& e : rx::dji::detect(x, ref["fs"], ref["fc"], false))
      if (e.family == "dji_droneid") em.push_back(e);
    const Json& want = c["emissions"];
    REQUIRE(em.size() == want.size());
    for (size_t i = 0; i < em.size(); ++i) {
      CHECK(em[i].start_us == doctest::Approx(want[i]["start_us"].get<double>()));
      CHECK(em[i].end_us == doctest::Approx(want[i]["end_us"].get<double>()));
      CHECK(std::fabs(em[i].center_mhz - want[i]["center_mhz"].get<double>()) <= 1e-4);
      const Json& w = want[i]["features"];
      for (const char* k : {"n_symbols", "droneid_frame", "burst_us", "zc_corr", "symbol0_cp", "cp_score"})
        CHECK(em[i].features[k].dump() == w[k].dump());
    }
  }
}
