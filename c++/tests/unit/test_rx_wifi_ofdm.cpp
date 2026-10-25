// receivers/wifi 802.11 OFDM (ofdm_rx, payload, ldpc, ofdm_scaled, ofdm_report) against vectors
// written by the Python receiver (tests/vectors/wifi/gen_wifi_ofdm.py).
#include <doctest/doctest.h>

#include <cmath>
#include <fstream>
#include <random>
#include <sstream>
#include <string>

#include "core/signal.hpp"
#include "receivers/wifi/ldpc.hpp"
#include "receivers/wifi/ofdm_report.hpp"
#include "receivers/wifi/ofdm_rx.hpp"
#include "receivers/wifi/ofdm_scaled.hpp"
#include "receivers/wifi/payload.hpp"

using namespace si;
using namespace si::rx::wifi;

namespace {

std::string vec_path(const std::string& name) { return std::string(SI_TESTVEC) + "/wifi/" + name; }

Json load_json(const std::string& name) {
  std::ifstream f(vec_path(name));
  REQUIRE_MESSAGE(f.good(), "missing vector " << name << " (run tests/vectors/wifi/gen_wifi_ofdm.py)");
  std::stringstream ss;
  ss << f.rdbuf();
  return Json::parse(ss.str());
}

cfvec load_cf32(const std::string& name) {
  std::ifstream f(vec_path(name), std::ios::binary | std::ios::ate);
  REQUIRE_MESSAGE(f.good(), "missing vector " << name);
  const auto bytes = static_cast<size_t>(f.tellg());
  cfvec v(bytes / sizeof(cf));
  f.seekg(0);
  f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(cf)));
  return v;
}

template <class T>
std::vector<T> vec(const Json& j) {
  return j.get<std::vector<T>>();
}
cdvec cvec(const Json& j) {
  cdvec v;
  for (const auto& p : j) v.push_back({p[0].get<double>(), p[1].get<double>()});
  return v;
}

// Python value == C++ value: numbers within rel (0 = exact), everything else equal
void same(const Json& py, const Json& cc, const std::string& path, double rel = 0) {
  if (py.is_number() && cc.is_number()) {
    const double a = py.get<double>(), b = cc.get<double>();
    CHECK_MESSAGE(std::fabs(a - b) <= rel * std::fabs(a), path << ": " << a << " vs " << b);
    return;
  }
  if (py.is_object() && cc.is_object()) {
    CHECK_MESSAGE(py.size() == cc.size(), path << ": keys " << py.dump().substr(0, 200) << " vs " << cc.dump().substr(0, 200));
    for (const auto& [k, v] : py.items()) {
      if (!cc.contains(k)) {
        CHECK_MESSAGE(false, path << "." << k << " missing in C++");
        continue;
      }
      same(v, cc[k], path + "." + k, rel);
    }
    return;
  }
  if (py.is_array() && cc.is_array()) {
    REQUIRE_MESSAGE(py.size() == cc.size(), path << ": length " << py.size() << " vs " << cc.size());
    for (size_t i = 0; i < py.size(); ++i) same(py[i], cc[i], path + "[" + std::to_string(i) + "]", rel);
    return;
  }
  CHECK_MESSAGE(py == cc, path << ": " << py.dump() << " vs " << cc.dump());
}

}  // namespace

TEST_CASE("wifi ofdm: Viterbi, convolutional encoder, CRC-8") {
  const Json j = load_json("ofdm_blocks.json")["viterbi"];
  const dvec soft = vec<double>(j["soft"]);
  CHECK(ofdm::viterbi(soft) == vec<uint8_t>(j["free"]));
  CHECK(ofdm::viterbi(soft, true) == vec<uint8_t>(j["terminated"]));
  CHECK(ofdm::conv_encode(vec<uint8_t>(j["bits"])) == vec<uint8_t>(j["encoded"]));
  const auto bits = vec<uint8_t>(j["bits"]);
  CHECK(ofdm::crc8_sig(bits.data(), 34) == vec<uint8_t>(j["crc8"]));
}

TEST_CASE("wifi ofdm: limit_denominator, channel plan, resampler") {
  const Json b = load_json("ofdm_blocks.json");
  for (const auto& c : b["limit_denominator"]) {
    const auto r = ofdm::limit_denominator(c[0].get<double>(), 4096);
    CHECK(r.first == c[1].get<int64_t>());
    CHECK(r.second == c[2].get<int64_t>());
  }
  for (const auto& c : b["channel_plan"]) {
    std::vector<std::string> names;
    for (const auto& ch : ofdm::channel_plan(c[0].get<double>(), c[1].get<double>(), 60e6)) names.push_back(ch.name);
    CHECK(names == vec<std::string>(c[2]));
  }
  // the receiver's resampler is bit-identical to sig::resample_poly
  std::mt19937 g(5);
  std::normal_distribution<float> nd;
  cfvec x(7001);
  for (auto& v : x) v = cf(nd(g), nd(g));
  for (auto [u, d] : {std::pair<int, int>{125, 384}, {125, 192}, {275, 1536}, {3, 1}, {1, 3}}) {
    const cfvec a = ofdm::detail::resample_cf(x, u, d), r = sig::resample_poly(x, u, d);
    REQUIRE(a.size() == r.size());
    size_t bad = 0;
    for (size_t i = 0; i < a.size(); ++i) bad += !(a[i] == r[i]);
    CHECK_MESSAGE(bad == 0, u << "/" << d);
  }
}

TEST_CASE("wifi ofdm: payload helpers") {
  const Json b = load_json("ofdm_blocks.json");
  for (const auto& c : b["deinterleave"])
    CHECK(payload::deinterleave_indices(c[0], c[1], c[2]).second == vec<int>(c[3]));
  const dvec s = vec<double>(b["depuncture_in"]);
  for (const auto& c : b["depuncture"])
    CHECK(payload::depuncture(s, {c[0][0].get<int>(), c[0][1].get<int>()}) == vec<double>(c[1]));
  const auto [d, seed] = payload::descramble(vec<uint8_t>(b["descramble"]["in"]));
  CHECK(d == vec<uint8_t>(b["descramble"]["out"]));
  CHECK(seed == vec<uint8_t>(b["descramble"]["seed"]));
  const Json& l = b["llr"];
  for (int nb : {1, 2, 4, 6}) {
    const dvec got = payload::llr_exact(cvec(l["sym"]), nb, vec<double>(l["var"]));
    same(l[std::to_string(nb)], Json(got), "llr" + std::to_string(nb), 1e-13);
  }
  const Json& sm = b["smooth"];
  for (int w : {1, 3}) {
    const cdvec got = payload::smooth_channel(cvec(sm["h"]), vec<int>(sm["sc"]), w);
    Json gj = Json::array();
    for (const auto& z : got) gj.push_back(Json::array({z.real(), z.imag()}));
    same(sm[std::to_string(w)], gj, "smooth" + std::to_string(w), 1e-12);
  }
  const auto parts = payload::deaggregate(vec<uint8_t>(b["deaggregate"]["in"]));
  const Json& po = b["deaggregate"]["out"];
  REQUIRE(parts.size() == po.size());
  for (size_t i = 0; i < parts.size(); ++i) {
    CHECK(parts[i].first == vec<uint8_t>(po[i][0]));
    CHECK(parts[i].second == po[i][1].get<bool>());
  }
}

TEST_CASE("wifi ofdm: LDPC parameters and min-sum decoding") {
  const Json b = load_json("ofdm_blocks.json");
  for (const auto& c : b["ldpc_params"]) {
    const auto P = ldpc::params(c[0], c[1], {c[2][0].get<int>(), c[2][1].get<int>()}, c[3]);
    const Json& e = c[4];
    CHECK(P.npld == e["npld"].get<int64_t>());
    CHECK(P.navbits == e["navbits"].get<int64_t>());
    CHECK(P.nsym == e["nsym"].get<int64_t>());
    CHECK(P.ncw == e["ncw"].get<int64_t>());
    CHECK(P.L == e["L"].get<int64_t>());
    CHECK(P.nshrt == e["nshrt"].get<int64_t>());
    CHECK(P.npunc == e["npunc"].get<int64_t>());
    CHECK(P.nrep == e["nrep"].get<int64_t>());
  }
  for (const auto& c : b["ldpc_decode"]) {
    const auto r = ldpc::decode(vec<double>(c["llr"]), c["L"], {c["rate"][0].get<int>(), c["rate"][1].get<int>()});
    CHECK(r.hard == vec<uint8_t>(c["hard"]));
    CHECK(r.ok == c["ok"].get<bool>());
    CHECK(r.used == c["used"].get<int>());
    if (r.ok) CHECK(r.hard == vec<uint8_t>(c["cw"]));
  }
  const Json& st = b["ldpc_stream"];
  const ldpc::Rate rate{5, 6};
  const auto P = ldpc::params(st["length"], st["ncbps"], rate, st["m_stbc"]);
  CHECK(ldpc::offsets(P) == vec<int64_t>(st["offsets"]));
  const auto r = ldpc::decode_stream(vec<double>(st["llr"]), P, rate);
  CHECK(r.info == vec<uint8_t>(st["info"]));
  CHECK(r.ok == vec<bool>(st["ok"]));
  CHECK(r.used == vec<int>(st["used"]));
}

TEST_CASE("wifi ofdm: synthetic scene, detector and records") {
  const Json j = load_json("ofdm_scene.json");
  const cfvec x = load_cf32("ofdm_scene.cf32");
  const double fs = j["fs"], fc = j["fc"];
  ofdm::WifiDetector det(fs, fc);
  const auto pk = det.run(x, 2);
  CHECK(det.n_rejects() == j["n_rejects"].get<size_t>());
  REQUIRE(pk.size() == j["packets"].size());
  for (size_t i = 0; i < pk.size(); ++i) {
    const auto& p = pk[i];
    Json c = Json::object();
    c["channel"] = p.channel;
    c["center_mhz"] = p.center_mhz;
    c["start_sample"] = p.start_sample;
    c["end_sample"] = p.end_sample;
    c["start_us"] = p.start_us;
    c["end_us"] = p.end_us;
    c["fmt"] = p.fmt;
    c["lsig_rate_mbps"] = p.lsig_rate_mbps;
    c["lsig_length"] = p.lsig_length;
    c["f_lo_mhz"] = p.f_lo_mhz;
    c["f_hi_mhz"] = p.f_hi_mhz;
    c["snr_db"] = p.snr_db;
    c["cfo_khz"] = p.cfo_khz;
    c["stf_metric"] = p.stf_metric;
    c["mcs"] = p.mcs ? Json(*p.mcs) : Json(nullptr);
    c["sig_crc_ok"] = p.sig_crc_ok ? Json(*p.sig_crc_ok) : Json(nullptr);
    c["energy_end_us"] = p.energy_end_us ? Json(*p.energy_end_us) : Json(nullptr);
    c["lsig_evm_db"] = p.lsig_evm_db;
    c["lsig_bit_agreement"] = p.lsig_bit_agreement;
    same(j["packets"][i], c, "packet" + std::to_string(i));
  }
  auto recs = ofdm_report::records(pk, fs, 20e6, false, 2);
  ofdm_report::attribute_vht(recs);
  same(j["records"], Json(recs), "records");
}

TEST_CASE("wifi ofdm: time-scaled 802.11 (10 MHz clock)") {
  const Json j = load_json("ofdm_scaled.json");
  const cfvec x = load_cf32("ofdm_scaled.cf32");
  const double fs = j["fs"], fc = j["fc"];
  const auto [bands, unions] = scaled::find_bands(x, fs, fc);
  REQUIRE(bands.size() == j["bands"].size());
  REQUIRE(unions.size() == j["unions"].size());
  for (size_t i = 0; i < bands.size(); ++i) {
    CHECK(bands[i].f_lo == j["bands"][i]["f_lo"].get<double>());
    CHECK(bands[i].f_hi == j["bands"][i]["f_hi"].get<double>());
    CHECK(bands[i].n_frames == j["bands"][i]["n_frames"].get<int64_t>());
  }
  std::vector<const scaled::Band*> cand;
  for (size_t i = 0; i < std::min<size_t>(12, bands.size()); ++i) cand.push_back(&bands[i]);
  for (size_t i = 0; i < std::min<size_t>(4, unions.size()); ++i) cand.push_back(&unions[i]);
  REQUIRE(cand.size() == j["numerology"].size());
  for (size_t i = 0; i < cand.size(); ++i) {
    const auto nm = scaled::numerology(x, fs, fc, *cand[i]);
    const Json& e = j["numerology"][i];
    REQUIRE(nm.has_value() == !e.is_null());
    if (!nm) continue;
    Json got = nm->json();
    if (nm->reject) got["reject"] = *nm->reject;
    same(e, got, "numerology" + std::to_string(i), 1e-9);
  }
  const auto links = scaled::scan(x, fs, fc);
  REQUIRE(links.size() == j["links"].size());
  for (size_t i = 0; i < links.size(); ++i) {
    const Json& e = j["links"][i];
    CHECK(links[i].clock_hz == e["clock_hz"].get<double>());
    CHECK(links[i].inverted == e["inverted"].get<bool>());
    CHECK(links[i].n_rejects == e["n_rejects"].get<size_t>());
    CHECK(links[i].center_hz == doctest::Approx(e["center_hz"].get<double>()).epsilon(1e-12));
    same(e["records"], Json(ofdm_report::records(links[i].packets, fs, links[i].clock_hz, true)), "scaled records");
  }
}
