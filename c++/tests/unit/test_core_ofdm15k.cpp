// core/ofdm15k against python core.ofdm15k on a DJI RC-uplink capture
// (tests/vectors/core/gen_ofdm15k.py): same hits, bursts, tracked hops and labels.
#include <doctest/doctest.h>

#include <cmath>
#include <fstream>
#include <string>

#include "core/ofdm15k.hpp"
#include "core/spectral.hpp"
#include "io/sigmf.hpp"

using namespace si;

namespace {

bool near(double a, double b, double rel = 1e-9) { return std::fabs(a - b) <= rel * std::max(1.0, std::fabs(b)); }

void check_burst(const ofdm15k::Burst& b, const Json& j) {
  CHECK(b.profile == j["profile"].get<std::string>());
  CHECK(near(b.start_s, j["start_s"].get<double>()));
  CHECK(near(b.end_s, j["end_s"].get<double>()));
  CHECK(near(b.f_off, j["f_off"].get<double>(), 1e-6));
  CHECK(near(b.bw_hz, j["bw_hz"].get<double>()));
  CHECK(b.nsym == j["nsym"].get<int64_t>());
  CHECK(b.cp_peaks == j["cp_peaks"].get<int64_t>());
  CHECK(b.period_us == j["period_us"].get<double>());
  CHECK(b.cp_scatter_us == j["cp_scatter_us"].get<double>());
  CHECK(b.tracked == j.value("tracked", false));
}

}  // namespace

TEST_CASE("ofdm15k: scan, refine, merge, track on D_20260916-171343") {
  std::ifstream f(std::string(SI_TESTVEC) + "/core/ofdm15k_D_20260916-171343.json");
  REQUIRE(f.good());
  const Json ref = Json::parse(f);
  const Capture cap = load_sigmf(std::string(SI_TESTDATA) + "/D_20260916-171343/D_20260916-171343.sigmf-meta");
  const double usable_half = std::min(cap.fs, 56e6) / 2;
  const auto spurs = find_spurs(cap.iq(), cap.fs, cap.fc);
  REQUIRE(spurs.size() == ref["spurs"].size());
  std::vector<std::pair<double, double>> notch;
  for (const auto& s : spurs) notch.emplace_back(s.offset_hz, 30e3);
  cdvec x(cap.n());
  for (size_t i = 0; i < x.size(); ++i) x[i] = cd(cap.iq()[i].real(), cap.iq()[i].imag());
  const SpectralChannelizer ch(x, cap.fs, notch);
  const ofdm15k::Slicer slicer = [&](double off, double half, double ofs) { return ch.slice(off, half, ofs); };

  const auto scan = ofdm15k::ofdm15k_scan(slicer, ofdm15k::profiles(), 0.5e-3, 7.0, usable_half);
  REQUIRE(scan.hits.size() == ref["hits"].size());
  for (size_t i = 0; i < scan.hits.size(); ++i) {
    CAPTURE(i);
    const Json& h = ref["hits"][i];
    CHECK(scan.hits[i].profile == h["profile"].get<std::string>());
    CHECK(scan.hits[i].f_off == h["f_off"].get<double>());
    CHECK(scan.hits[i].t0_s == h["t0_s"].get<double>());
    CHECK(near(scan.hits[i].score, h["score"].get<double>(), 1e-9));
  }
  for (const auto& [name, st] : ref["stats"].items()) {
    CHECK(near(scan.stats.at(name).median, st["median"].get<double>(), 1e-9));
    CHECK(near(scan.stats.at(name).sigma, st["sigma"].get<double>(), 1e-9));
  }
  std::vector<ofdm15k::Burst> refined;
  for (size_t i = 0; i < scan.hits.size(); ++i) {
    CAPTURE(i);
    const auto b = ofdm15k::refine_ofdm15k(slicer, scan.hits[i]);
    REQUIRE(b.has_value() == !ref["refined"][i].is_null());
    if (b) {
      check_burst(*b, ref["refined"][i]);
      refined.push_back(*b);
    }
  }
  const auto ob = ofdm15k::merge_ofdm15k(refined);
  REQUIRE(ob.size() == ref["merged"].size());
  const auto tracked = ofdm15k::track_hops(slicer, ob, scan.stats, cap.duration(), 4.0, usable_half);
  REQUIRE(tracked.size() == ref["tracked"].size());
  for (size_t i = 0; i < tracked.size(); ++i) check_burst(tracked[i], ref["tracked"][i]);
  auto all = ob;
  all.insert(all.end(), tracked.begin(), tracked.end());
  const auto fin = ofdm15k::merge_ofdm15k(all);
  REQUIRE(fin.size() == ref["final"].size());
  for (size_t i = 0; i < fin.size(); ++i) {
    CAPTURE(i);
    check_burst(fin[i], ref["final"][i]);
    CHECK(ofdm15k::label_ofdm15k(fin[i]).first == ref["labels"][i][0].get<std::string>());
  }
}
