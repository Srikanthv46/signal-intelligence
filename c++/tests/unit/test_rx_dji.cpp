// receivers/dji against python receivers.dji (tests/vectors/dji/gen_dji.py): ZC helpers, CP
// metric / fold / segmentation, frame period, a synthetic DroneID burst and detect() on the
// three DJI testdata captures (OcuSync 2 video, O4 video, RC uplink).
#include <doctest/doctest.h>

#include <cmath>
#include <fstream>
#include <string>

#include "core/fft.hpp"
#include "io/sigmf.hpp"
#include "receivers/dji/detector.hpp"

using namespace si;
using namespace si::rx::dji;

namespace {

const Json& ref() {
  static const Json r = [] {
    std::ifstream f(std::string(SI_TESTVEC) + "/dji/ref.json");
    REQUIRE(f.good());
    return Json::parse(f);
  }();
  return r;
}

template <class T>
std::vector<T> read_bin(const std::string& name) {
  std::ifstream f(std::string(SI_TESTVEC) + "/dji/" + name, std::ios::binary | std::ios::ate);
  REQUIRE(f.good());
  const size_t n = static_cast<size_t>(f.tellg()) / sizeof(T);
  std::vector<T> v(n);
  f.seekg(0);
  f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(n * sizeof(T)));
  return v;
}

cdvec cplx(const Json& j) {
  cdvec v;
  for (const auto& z : j) v.emplace_back(z[0].get<double>(), z[1].get<double>());
  return v;
}

bool near(double a, double b, double rel = 1e-12) { return std::fabs(a - b) <= rel * std::max(1.0, std::fabs(b)); }

void check_emissions(const std::vector<Emission>& em, const Json& ref_em) {
  REQUIRE(em.size() == ref_em.size());
  for (size_t i = 0; i < em.size(); ++i) {
    CAPTURE(i);
    const Json& r = ref_em[i];
    CHECK(em[i].family == r["family"].get<std::string>());
    CHECK(em[i].label == r["label"].get<std::string>());
    CHECK(em[i].start_us == r["start_us"].get<double>());
    CHECK(em[i].end_us == r["end_us"].get<double>());
    CHECK(em[i].f_lo_mhz == r["f_lo_mhz"].get<double>());
    CHECK(em[i].f_hi_mhz == r["f_hi_mhz"].get<double>());
    CHECK(em[i].center_mhz == r["center_mhz"].get<double>());
    CHECK(em[i].snr_db.has_value() == !r["snr_db"].is_null());
    if (em[i].snr_db && !r["snr_db"].is_null()) CHECK(*em[i].snr_db == r["snr_db"].get<double>());
    CHECK(near(em[i].confidence, r["confidence"].get<double>()));
    CHECK(em[i].features.dump() == r["features"].dump());
  }
}

}  // namespace

TEST_CASE("dji: Zadoff-Chu helpers") {
  const Json& R = ref();
  const cdvec z = zc_freq(137, 601);
  cdvec pick(z.begin(), z.begin() + 4);
  pick.insert(pick.end(), z.begin() + 298, z.begin() + 303);
  pick.insert(pick.end(), z.end() - 3, z.end());
  const cdvec zr = cplx(R["zc_freq_137_601"]);
  REQUIRE(pick.size() == zr.size());
  for (size_t i = 0; i < pick.size(); ++i) {
    CAPTURE(i);
    CHECK(pick[i] == zr[i]);
  }
  const Json& s = R["zc_search_73"];
  const cdvec X = cplx(s["X"]);
  const auto [cc, r] = zc_root_search(X, 73);
  CHECK(r == s["root"].get<int64_t>());
  CHECK(near(cc, s["cc"].get<double>(), 1e-12));
  const auto [cc1, r1] = zc_root_search(X, 73, {50});
  CHECK(r1 == R["zc_search_73_fixed"]["root"].get<int64_t>());
  CHECK(near(cc1, R["zc_search_73_fixed"]["cc"].get<double>(), 1e-12));

  const Json& d = R["zc_dd"];
  const auto [acc, mag] = detail::zc_dd(cplx(d["X"]), 601);
  CHECK(acc.real() == d["acc"][0].get<double>());
  CHECK(acc.imag() == d["acc"][1].get<double>());
  CHECK(mag == d["mag"].get<double>());
  CHECK(detail::root_from_dd(acc, 601) == d["root"].get<int64_t>());
}

TEST_CASE("dji: CP metric, fold, robust z, segmentation") {
  const Json& c = ref()["cp"];
  const cdvec y = read_bin<cd>("cp_y.c128");
  REQUIRE(y.size() == c["n"].get<size_t>());
  const auto co = detail::cp_metric(y, 1024, 72);
  for (size_t i = 0; i < 5; ++i) CHECK(co.m[i] == c["m_head"][i].get<double>());
  CHECK(co.cs.back().real() == c["cs_last"][0].get<double>());
  CHECK(co.cs.back().imag() == c["cs_last"][1].get<double>());
  double s = 0;  // np.sum is pairwise; compare to rounding
  for (double v : co.m) s += v;
  CHECK(near(s, c["m_sum"].get<double>(), 1e-12));
  const double P = 0.5e-3 * 15.36e6 / 7;
  const dvec fz = detail::robust_z(detail::fold(co.m, P));
  const size_t k = static_cast<size_t>(std::max_element(fz.begin(), fz.end()) - fz.begin());
  CHECK(k == c["fold_argmax"].get<size_t>());
  CHECK(fz[k] == c["fold_max"].get<double>());
  const auto sg = detail::segment(co.m, static_cast<int64_t>(k), P);
  REQUIRE(sg.runs.size() == c["runs"].size());
  for (size_t i = 0; i < sg.runs.size(); ++i) {
    CHECK(sg.runs[i].first == c["runs"][i][0].get<size_t>());
    CHECK(sg.runs[i].second == c["runs"][i][1].get<size_t>());
  }
  for (size_t i = 0; i < 5; ++i) CHECK(sg.t[i] == c["t"][i].get<double>());
  REQUIRE(sg.on.size() == c["on"].size());
  for (size_t i = 0; i < sg.on.size(); ++i) {
    CHECK(sg.on[i] == c["on"][i].get<double>());
    CHECK(sg.off[i] == c["off"][i].get<double>());
  }
}

TEST_CASE("dji: frame period") {
  for (const auto& fp : ref()["frame_period"]) {
    const auto T = detail::frame_period(fp["s"].get<dvec>());
    CHECK(T.has_value() == !fp["T"].is_null());
    if (T && !fp["T"].is_null()) CHECK(T->T == fp["T"].get<double>());
  }
}

TEST_CASE("dji: synthetic DroneID burst") {
  const Json& d = ref()["droneid"];
  const cfvec x = read_bin<cf>("droneid_cf32.bin");
  const double fs = d["fs"].get<double>(), fc = d["fc"].get<double>();
  cdvec xd(x.begin(), x.end());
  const cdvec X = fft::fft(xd);
  const detail::FastSlicer ch(X, fs);
  const dvec sc = detail::droneid_cp_score(ch.slice(-0.18e6, 5.0e6, fs));
  const size_t k = static_cast<size_t>(std::max_element(sc.begin(), sc.end()) - sc.begin());
  CHECK(k == d["score_argmax"].get<size_t>());
  CHECK(near(sc[k], d["score_max"].get<double>(), 1e-12));
  check_emissions(detect(x, fs, fc, false), d["emissions"]);
}

TEST_CASE("dji: detect() on the DJI testdata captures") {
  for (const auto& [rel, em] : ref()["captures"].items()) {
    CAPTURE(rel);
    const Capture cap = load_sigmf(std::string(SI_TESTDATA) + "/" + rel + ".sigmf-meta");
    check_emissions(detect(cap.iq(), cap.fs, cap.fc, false), em);
  }
}
