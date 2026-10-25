// receivers/digital_vtx against python receivers.digital_vtx (tests/vectors/digital_vtx/
// gen_digital_vtx.py): PSS template, FFT slice, PSS fold search, carrier search, the blind
// CP-OFDM link building blocks (blind.hpp), numerology table and detect() on a synthetic capture.
#include <doctest/doctest.h>

#include <cmath>
#include <fstream>
#include <string>
#include <tuple>
#include <algorithm>
#include <optional>

#include "core/bands_india.hpp"
#include "core/fft.hpp"
#include "receivers/digital_vtx/blind.hpp"
#include "receivers/digital_vtx/detector.hpp"
#include "receivers/digital_vtx/hdzero.hpp"
#include "receivers/digital_vtx/catalog.hpp"
#include "receivers/digital_vtx/streams.hpp"

using namespace si;
using namespace si::rx::digital_vtx;

namespace {

const Json& ref() {
  static const Json r = [] {
    std::ifstream f(std::string(SI_TESTVEC) + "/digital_vtx/ref.json");
    REQUIRE(f.good());
    return Json::parse(f);
  }();
  return r;
}

const cfvec& synth() {
  static const cfvec x = [] {
    std::ifstream f(std::string(SI_TESTVEC) + "/digital_vtx/synth.cf32", std::ios::binary | std::ios::ate);
    REQUIRE(f.good());
    cfvec v(static_cast<size_t>(f.tellg()) / sizeof(cf));
    f.seekg(0);
    f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(cf)));
    return v;
  }();
  return x;
}

const cdvec& spectrum() {
  static const cdvec X = [] {
    cdvec xd(synth().begin(), synth().end());
    return fft::fft(xd);
  }();
  return X;
}

void check_cplx(const cdvec& got, const Json& want) {
  REQUIRE(got.size() >= want.size());
  for (size_t i = 0; i < want.size(); ++i) {
    CHECK(got[i].real() == want[i][0].get<double>());
    CHECK(got[i].imag() == want[i][1].get<double>());
  }
}

}  // namespace

TEST_CASE("digital_vtx: PSS template and FFT slice") {
  check_cplx(detail::pss_time(29, 1.92e6), ref()["pss_time_29_192"]);
  check_cplx(detail::pss_time(34, 3.84e6), ref()["pss_time_34_384"]);
  const Json& s = ref()["slice"];
  double ofs = 0;
  const cdvec y = detail::slice(spectrum(), ref()["fs"], 1.2345e6, 0.5e6, 0.96e6, &ofs);
  CHECK(y.size() == s["n"].get<size_t>());
  CHECK(ofs == s["ofs"].get<double>());
  check_cplx(y, s["y"]);
}

TEST_CASE("digital_vtx: PSS fold search and carrier search") {
  const Json& w = ref()["fold"];
  const auto sr = detail::pss_fold_search(spectrum(), ref()["fs"], 3.0e6);
  CHECK(sr.ofs == w["ofs"].get<double>());
  CHECK(sr.M == w["M"].get<int64_t>());
  CHECK(sr.P == w["P"].get<double>());
  CHECK(sr.nt == w["nt"].get<size_t>());
  CHECK(sr.global_median == w["global_median"].get<double>());
  REQUIRE(sr.offs.size() == w["offs"].size());
  for (size_t i = 0; i < sr.offs.size(); ++i) CHECK(sr.offs[i] == w["offs"][i].get<double>());
  for (int u : {25, 29, 34}) {
    const Json& z = w["Z"][std::to_string(u)];
    const Json& t = w["T"][std::to_string(u)];
    size_t bad = 0;
    for (size_t i = 0; i < z.size(); ++i)
      bad += sr.Z.at(u)[i] != z[i].get<double>() || sr.T.at(u)[i] != t[i].get<int64_t>();
    CHECK_MESSAGE(bad == 0, "root " << u);
  }
  const Json& c = ref()["carrier"];
  const auto hit = detail::carrier_search(spectrum(), ref()["fs"], c["f0"], c["t_inst"].get<dvec>());
  CHECK(hit.f_c == c["f_c"].get<double>());
  CHECK(hit.root == c["root"].get<int>());
  CHECK(hit.score == doctest::Approx(c["score"].get<double>()).epsilon(1e-12));
}

namespace {
blind::Slicer slicer() {
  return [](double f, double h, double r, double* o) { return detail::slice(spectrum(), ref()["fs"], f, h, r, o); };
}
void check_runs(const std::vector<blind::Run>& got, const Json& want) {
  REQUIRE(got.size() == want.size());
  for (size_t i = 0; i < got.size(); ++i) {
    CHECK(got[i].t0 == doctest::Approx(want[i][0].get<double>()).epsilon(1e-12));
    CHECK(got[i].t1 == doctest::Approx(want[i][1].get<double>()).epsilon(1e-12));
    CHECK(got[i].n == want[i][2].get<int64_t>());
  }
}
}  // namespace

TEST_CASE("digital_vtx blind: rate, scan, confirmation") {
  const double t[][2] = {{66.667e-6, 18e6}, {33.3335e-6, 9e6}, {16.6667e-6, 1e6}};
  for (size_t i = 0; i < 3; ++i) {
    const auto [rate, k] = blind::rate_for(t[i][0], t[i][1]);
    CHECK(rate == ref()["rate_for"][i][0].get<double>());
    CHECK(k == ref()["rate_for"][i][1].get<int64_t>());
  }
  const double usable = ref()["fs"].get<double>() / 2;
  auto hits = blind::scan(slicer(), usable);
  REQUIRE(hits.size() == ref()["n_scan"].get<size_t>());
  std::stable_sort(hits.begin(), hits.end(), [](const blind::Hit& a, const blind::Hit& b) {
    return std::make_tuple(-a.z, a.W, a.off, a.N) < std::make_tuple(-b.z, b.W, b.off, b.N);
  });
  for (size_t i = 0; i < ref()["scan"].size(); ++i) {
    const Json& w = ref()["scan"][i];
    CHECK(hits[i].N == w["N"].get<int64_t>());
    CHECK(hits[i].off == w["off"].get<double>());
    CHECK(hits[i].ofs == w["ofs"].get<double>());
    CHECK(hits[i].z == doctest::Approx(w["z"].get<double>()).epsilon(1e-9));
  }
  const Json& c = ref()["confirm"];
  const auto cf = blind::confirm(slicer(), c["off"], c["W"], c["N"], c["ofs"]);
  REQUIRE(cf.has_value() == !c["result"].is_null());
  CHECK(cf->tu == doctest::Approx(c["result"][0].get<double>()).epsilon(1e-12));
  CHECK(cf->ts == doctest::Approx(c["result"][1].get<double>()).epsilon(1e-12));
  CHECK(cf->line_z == doctest::Approx(c["result"][2].get<double>()).epsilon(1e-9));
  CHECK(cf->nwin == c["result"][3].get<int64_t>());
}

TEST_CASE("digital_vtx blind: bursts, frame, band") {
  const Json& c = ref()["confirm"];
  const double tu = c["result"][0], ts = c["result"][1];
  const auto sg = blind::segment(slicer(), c["off"], c["W"], tu, ts);
  const Json& s = ref()["segment"];
  REQUIRE(sg.tk.size() == s["n_tk"].get<size_t>());
  for (size_t i = 0; i < s["tk"].size(); ++i) {
    CHECK(sg.tk[i] == doctest::Approx(s["tk"][i].get<double>()).epsilon(1e-12));
    CHECK(sg.zr[i] == doctest::Approx(s["zr"][i].get<double>()).epsilon(1e-9));
  }
  CHECK(sg.ts == doctest::Approx(s["ts"].get<double>()).epsilon(1e-12));
  CHECK(sg.ngrid == s["ngrid"].get<int64_t>());
  CHECK(sg.frac == doctest::Approx(s["frac"].get<double>()).epsilon(1e-9));
  check_runs(sg.runs, s["runs"]);
  const auto [T, con] = blind::frame_period(sg.tk, sg.zr, sg.ts);
  REQUIRE(T.has_value() == !ref()["frame_period"][0].is_null());
  CHECK(*T == doctest::Approx(ref()["frame_period"][0].get<double>()).epsilon(1e-12));
  CHECK(con == doctest::Approx(ref()["frame_period"][1].get<double>()).epsilon(1e-9));
  check_runs(blind::frame_complete(sg.runs, sg.tk, sg.zr, sg.ts, *T), ref()["frame_complete"]);
  dvec tsym;
  for (size_t i = 0; i < sg.tk.size(); ++i)
    if (sg.zr[i] > 2) tsym.push_back(sg.tk[i]);
  const auto [fr, cc] = blind::band_cp(slicer(), c["off"], c["W"], tu, sg.ts, tsym);
  const Json& b = ref()["band_cp"];
  REQUIRE(cc.size() == b["n"].get<size_t>());
  for (size_t i = 0; i < b["c"].size(); ++i) {
    CHECK(fr[8 * i] == doctest::Approx(b["fr"][i].get<double>()).epsilon(1e-12));
    CHECK(cc[8 * i] == doctest::Approx(b["c"][i].get<double>()).epsilon(1e-9));
  }
  const double usable = ref()["fs"].get<double>() / 2;
  const auto mb = blind::measure_band(slicer(), c["off"], c["W"], tu, sg.ts, tsym, usable);
  REQUIRE(mb.has_value());
  CHECK(mb->first == ref()["measure_band"][0].get<double>());
  CHECK(mb->second == ref()["measure_band"][1].get<double>());
}

TEST_CASE("digital_vtx: numerology table") {
  for (const auto& [k, v] : ref()["classify"].items()) {
    const auto bar = k.find('|');
    const double t = std::stod(k.substr(0, bar)), bw = std::stod(k.substr(bar + 1));
    const auto got = classify_numerology(t, bw > 0 ? std::optional<double>(bw) : std::nullopt);
    CHECK(got.first == v[0].get<std::string>());
    CHECK(Json(got.second) == v[1]);
  }
}

TEST_CASE("digital_vtx: detect() on the synthetic capture") {
  const auto em = detect(synth(), ref()["fs"], ref()["fc"], false);
  const Json& w = ref()["detect"];
  REQUIRE(em.size() == w.size());
  for (size_t i = 0; i < em.size(); ++i) {
    const Json& e = w[i];
    CHECK(em[i].family == e["family"].get<std::string>());
    CHECK(em[i].label == e["label"].get<std::string>());
    CHECK(em[i].start_us == e["start_us"].get<double>());
    CHECK(em[i].end_us == e["end_us"].get<double>());
    CHECK(em[i].f_lo_mhz == e["f_lo_mhz"].get<double>());
    CHECK(em[i].f_hi_mhz == e["f_hi_mhz"].get<double>());
    CHECK(em[i].center_mhz == e["center_mhz"].get<double>());
    CHECK(em[i].confidence == e["confidence"].get<double>());
    CHECK((em[i].snr_db ? Json(*em[i].snr_db) : Json(nullptr)) == e["snr_db"]);
    CHECK(em[i].features == e["features"]);
  }
}

TEST_CASE("digital_vtx: transmit streams (1 TX vs Alamouti STBC)") {
  for (const auto& [kind, r] : ref()["tx_streams"].items()) {
    std::ifstream f(std::string(SI_TESTVEC) + "/digital_vtx/streams_" + kind + ".cf32", std::ios::binary | std::ios::ate);
    REQUIRE(f);
    const auto n = static_cast<size_t>(f.tellg()) / sizeof(cf);
    cfvec v(n);
    f.seekg(0);
    f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(n * sizeof(cf)));
    const cdvec y(v.begin(), v.end());
    std::vector<blind::Run> runs;
    for (const auto& q : r["runs"]) runs.push_back({q[0].get<double>(), q[1].get<double>(), q[2].get<int64_t>()});
    const Json got = streams::tx_streams(y, r["ofs"], r["tu"], r["ts"], r["bw"], runs, 30.0, r["dur"]);
    const Json& want = r["out"];
    CHECK(got["verdict"] == want["verdict"]);
    CHECK(got["tx_antennas"] == want["tx_antennas"]);
    for (const char* k : {"segments", "symbols", "subcarriers"}) CHECK(got[k] == want[k]);
    for (const char* k : {"kurtosis", "time_pair_delta", "freq_pair_delta"})
      CHECK(std::fabs(got[k].get<double>() - want[k].get<double>()) <= 1e-3 + 1e-9);
    CHECK(got["note"] == want["note"]);
  }
}

TEST_CASE("digital_vtx: research catalog candidates == python") {
  for (const auto& c : ref()["catalog"]) {
    const Json& a = c["args"];
    auto opt = [](const Json& v) { return v.is_null() ? std::optional<double>() : std::optional<double>(v.get<double>()); };
    const std::optional<int> na = a[4].is_null() ? std::optional<int>() : std::optional<int>(a[4].get<int>());
    const Json got = catalog::candidates(a[0].get<double>(), opt(a[1]), opt(a[2]), opt(a[3]), na, opt(a[5]));
    CHECK(got.dump() == c["out"].dump());
    CHECK(catalog::model_candidates(got, opt(a[3]), opt(a[2]), na).dump() == c["model"].dump());
  }
}

TEST_CASE("digital_vtx: DJI uplink hop-timing model == python") {
  for (const auto& c : ref()["hop_model"]) {
    const Json& a = c[0];
    auto opt = [](const Json& v) { return v.is_null() ? std::optional<double>() : std::optional<double>(v.get<double>()); };
    CHECK(catalog::hop_model(opt(a[0]), opt(a[1]), opt(a[2])).dump() == c[1].dump());
  }
}

TEST_CASE("digital_vtx: HDZero slice plan, channel search and detect() == python") {
  const Json& r = ref()["hdzero"];
  for (const auto& p : r["plans"]) {
    const auto got = hdzero::slice_plan(p[0].get<int64_t>(), p[1].get<double>());
    CHECK(got.first == p[2].get<int64_t>());
    CHECK(got.second == p[3].get<int64_t>());
  }
  std::ifstream f(std::string(SI_TESTVEC) + "/digital_vtx/hdzero.cf32", std::ios::binary | std::ios::ate);
  REQUIRE(f);
  cfvec v(static_cast<size_t>(f.tellg()) / sizeof(cf));
  f.seekg(0);
  f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(cf)));
  const double fs = r["fs"], fc = r["fc"];
  const cdvec X = fft::fft(cdvec(v.begin(), v.end()));
  std::vector<hdzero::LogEntry> log;
  const hdzero::Slicer sl = [&](double ff, double h, double rr, double* o) { return detail::slice(X, fs, ff, h, rr, o); };
  const auto links = hdzero::hdzero_links(sl, fs, fc, fs / 2, static_cast<int64_t>(v.size()), &log);
  REQUIRE(log.size() == r["log"].size());
  for (size_t i = 0; i < log.size(); ++i) {
    const Json& w = r["log"][i];
    CHECK(log[i].channel == w["channel"].get<std::string>());
    CHECK(log[i].shift == w["shift"].get<int64_t>());
    CHECK(log[i].accepted == w["accepted"].get<bool>());
    CHECK(py_round(log[i].pilot_coherence, 3) == w["pilot_coherence"].get<double>());
    CHECK(py_round(log[i].data_coherence, 3) == w["data_coherence"].get<double>());
    CHECK(py_round(log[i].fold_contrast, 2) == w["fold_contrast"].get<double>());
  }
  REQUIRE(links.size() == r["links"].size());
  for (size_t i = 0; i < links.size(); ++i) {
    const Json& w = r["links"][i];
    CHECK(links[i].channel == w["channel"].get<std::string>());
    CHECK(links[i].cfo_hz == doctest::Approx(w["cfo_hz"].get<double>()).epsilon(1e-9));
    CHECK(links[i].p_band == doctest::Approx(w["p_band"].get<double>()).epsilon(1e-9));
    CHECK(links[i].n_symbols == w["n_symbols"].get<int64_t>());
    CHECK(Json(links[i].preamble_ms) == w["preamble_ms"]);
    CHECK(Json(links[i].preamble_symbols) == w["preamble_symbols"]);
    CHECK((links[i].frame_symbols ? Json(*links[i].frame_symbols) : Json(nullptr)) == w["frame_symbols"]);
  }
  const auto em = detect(v, fs, fc, false);
  const Json& we = r["detect"];
  REQUIRE(em.size() == we.size());
  for (size_t i = 0; i < em.size(); ++i) {
    CHECK(em[i].family == we[i]["family"].get<std::string>());
    CHECK(em[i].label == we[i]["label"].get<std::string>());
    CHECK(em[i].f_lo_mhz == we[i]["f_lo_mhz"].get<double>());
    CHECK(em[i].f_hi_mhz == we[i]["f_hi_mhz"].get<double>());
    CHECK((em[i].snr_db ? Json(*em[i].snr_db) : Json(nullptr)) == we[i]["snr_db"]);
    CHECK(em[i].features.dump() == we[i]["features"].dump());
  }
  CHECK(catalog::hdzero_model("27 MHz (normal)").dump() == r["model"].dump());
}

TEST_CASE("digital_vtx: wider search bands and band context (PLAN G8)") {
  CHECK(bands::uav_link_band(1437e6).value_or("").rfind("1.4 GHz UAV data-link band", 0) == 0);
  CHECK(bands::uav_link_band(1380e6).value_or("").rfind("L-band", 0) == 0);
  CHECK(bands::uav_link_band(2250e6).value_or("").rfind("1.6-2.33 GHz", 0) == 0);
  CHECK(bands::uav_link_band(4600e6).value_or("").rfind("4.4-4.9 GHz", 0) == 0);
  CHECK(!bands::uav_link_band(5800e6));
  CHECK(bands::cellular_at(1840e6));   // B3 downlink
  CHECK(bands::cellular_at(2350e6));   // B40
  CHECK(!bands::cellular_at(2250e6));
  CHECK(!bands::cellular_at(1437e6));
}
