// receivers/wifi (MAC dissector, DroneID, link attribution, 802.11b DSSS/CCK) against vectors
// written by the Python receiver (tests/vectors/wifi/gen_wifi.py).
#include <doctest/doctest.h>

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

#include "core/dsp.hpp"
#include "core/signal.hpp"
#include "receivers/wifi/droneid.hpp"
#include "receivers/wifi/dsss.hpp"
#include "receivers/wifi/links.hpp"
#include "receivers/wifi/mac80211.hpp"
#include "receivers/wifi/oui_table.hpp"

using namespace si;
using namespace si::rx::wifi;

namespace {

std::string vec_path(const std::string& name) { return std::string(SI_TESTVEC) + "/wifi/" + name; }

Json load_json(const std::string& name) {
  std::ifstream f(vec_path(name));
  REQUIRE_MESSAGE(f.good(), "missing vector " << name << " (run tests/vectors/wifi/gen_wifi.py)");
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

std::vector<uint8_t> unhex(const std::string& h) {
  std::vector<uint8_t> out;
  for (size_t i = 0; i + 1 < h.size(); i += 2) out.push_back(static_cast<uint8_t>(std::stoi(h.substr(i, 2), nullptr, 16)));
  return out;
}

template <class T>
std::vector<T> as_vec(const Json& j) {
  return j.get<std::vector<T>>();
}

// Structural comparison: same keys in the same order, same strings / bools / nulls, numbers
// equal within rel_tol (0: exact). Returns the first difference ("" when equal).
std::string diff(const Json& a, const Json& b, const std::string& path = "$", double rel_tol = 0.0) {
  if (a.is_object() && b.is_object()) {
    if (a.size() != b.size()) {
      std::string ka, kb;
      for (auto it = a.begin(); it != a.end(); ++it) ka += it.key() + ",";
      for (auto it = b.begin(); it != b.end(); ++it) kb += it.key() + ",";
      return path + ": keys {" + ka + "} vs {" + kb + "}";
    }
    auto ia = a.begin();
    auto ib = b.begin();
    for (; ia != a.end(); ++ia, ++ib) {
      if (ia.key() != ib.key()) return path + ": key " + ia.key() + " vs " + ib.key();
      auto d = diff(ia.value(), ib.value(), path + "." + ia.key(), rel_tol);
      if (!d.empty()) return d;
    }
    return "";
  }
  if (a.is_array() && b.is_array()) {
    if (a.size() != b.size()) return path + ": " + std::to_string(a.size()) + " vs " + std::to_string(b.size()) + " items";
    for (size_t i = 0; i < a.size(); ++i) {
      auto d = diff(a[i], b[i], path + "[" + std::to_string(i) + "]", rel_tol);
      if (!d.empty()) return d;
    }
    return "";
  }
  if (a.is_number() && b.is_number() && !a.is_boolean() && !b.is_boolean()) {
    const double x = a.get<double>(), y = b.get<double>();
    if (x == y || std::fabs(x - y) <= rel_tol * std::max(std::fabs(x), std::fabs(y))) return "";
    return path + ": " + a.dump() + " vs " + b.dump();
  }
  if (a != b) return path + ": " + a.dump().substr(0, 200) + " vs " + b.dump().substr(0, 200);
  return "";
}

}  // namespace

TEST_CASE("wifi mac80211: parse() equals python on real and synthetic MPDUs") {
  const Json v = load_json("mac_frames.json");
  CHECK(detail::oui_registry_size() == v["oui_db_size"].get<size_t>());
  size_t n = 0;
  for (const auto& c : v["frames"]) {
    CAPTURE(c["name"].get<std::string>());
    CAPTURE(c["hex"].get<std::string>().substr(0, 120));
    const Json got = mac80211::parse(unhex(c["hex"].get<std::string>()));
    CHECK_MESSAGE(diff(c["parsed"], got) == "", diff(c["parsed"], got));
    ++n;
  }
  CHECK(n >= 80);
  for (const auto& c : v["vendor"]) {
    auto got = mac80211::vendor(c["mac"].get<std::string>());
    CHECK(Json(got ? Json(*got) : Json(nullptr)) == c["vendor"]);
  }
  for (const auto& c : v["oui_name"]) CHECK(mac80211::oui_name(c["oui"].get<std::string>()) == c["name"]);
  for (const auto& c : v["crc32"]) {
    auto b = unhex(c["hex"].get<std::string>());
    CHECK(mac80211::crc32(b.data(), b.size()) == c["crc"].get<uint32_t>());
  }
}

TEST_CASE("wifi droneid: Open Drone ID, DJI DroneID, signatures") {
  const Json v = load_json("droneid.json");
  for (auto it = v["cases"].begin(); it != v["cases"].end(); ++it) {
    const std::string name = it.key();
    const Json& c = it.value();
    CAPTURE(name);
    const auto b = unhex(c["hex"].get<std::string>());
    const auto got = name.rfind("odid", 0) == 0 ? droneid::parse_odid(b) : droneid::parse_dji(b);
    const Json gj = got ? *got : Json(nullptr);
    CHECK_MESSAGE(diff(c["out"], gj) == "", diff(c["out"], gj));
    if (got) CHECK(droneid::summary(*got) == c["summary"].get<std::string>());
  }
  auto opt = [](const Json& j) { return j.is_null() ? std::optional<std::string>() : j.get<std::string>(); };
  for (const auto& c : v["signature"]) {
    auto got = droneid::signature(opt(c["vendor"]), opt(c["ssid"]));
    const Json gj = got ? *got : Json(nullptr);
    CHECK_MESSAGE(diff(c["out"], gj) == "", c.dump());
  }
  // every registry name: exactly the ones Python's regexes match
  size_t hits = 0;
  for (const auto& c : v["vendor_hits"]) {
    auto got = droneid::signature(c["vendor"].get<std::string>(), std::nullopt);
    REQUIRE(got);
    CHECK(diff(c["out"], *got) == "");
    ++hits;
  }
  CHECK(hits == v["vendor_hits"].size());
}

TEST_CASE("wifi links: attribute() equals python") {
  const Json v = load_json("links.json");
  std::vector<Json> recs = v["input"].get<std::vector<Json>>();
  const Json rows = links::attribute(recs, 16.0);
  CHECK_MESSAGE(diff(v["rows"], rows, "$", 1e-12) == "", diff(v["rows"], rows, "$", 1e-12));
  REQUIRE(recs.size() == v["records"].size());
  for (size_t i = 0; i < recs.size(); ++i) CHECK_MESSAGE(diff(v["records"][i], recs[i]) == "", diff(v["records"][i], recs[i]));
}

TEST_CASE("wifi dsss: descrambler, CRC-16, bytes, CCK codebooks, DBPSK/DQPSK, SFD search") {
  const Json v = load_json("dsss_blocks.json");
  using dsss::Bits;
  CHECK(dsss::descramble(as_vec<uint8_t>(v["descramble_in"])) == as_vec<uint8_t>(v["descramble_out"]));
  const auto hdr = as_vec<uint8_t>(v["crc_hdr"]), bad = as_vec<uint8_t>(v["crc_bad"]);
  CHECK(dsss::crc16_ok(hdr.data(), hdr.data() + 32) == v["crc_ok"].get<bool>());
  CHECK(dsss::crc16_ok(bad.data(), bad.data() + 32) == v["crc_bad_ok"].get<bool>());
  CHECK(v["crc_ok"].get<bool>());
  CHECK(dsss::bits_to_bytes(as_vec<uint8_t>(v["bytes_in"])) == as_vec<uint8_t>(v["bytes_out"]));
  for (const auto& p : v["psdu_bytes"])
    CHECK(dsss::psdu_bytes(p[0].get<double>(), p[1].get<int64_t>(), p[2].get<int64_t>()) == p[3].get<int>());
  for (const auto& [book, key] : {std::pair<const dsss::CckBook*, std::string>{&dsss::cck11(), "cck11"}, {&dsss::cck55(), "cck55"}}) {
    const auto re = v[key + "_re"], im = v[key + "_im"], bits = v[key + "_bits"];
    REQUIRE(re.size() == book->ncode());
    for (size_t m = 0; m < book->ncode(); ++m) {
      CHECK(as_vec<uint8_t>(bits[m]) == book->bits[m]);
      for (size_t j = 0; j < 8; ++j) {
        CHECK(book->C[m * 8 + j].real() == re[m][j].get<double>());
        CHECK(book->C[m * 8 + j].imag() == im[m][j].get<double>());
      }
    }
  }
  auto cvec = [&](const char* r, const char* i) {
    const auto a = as_vec<double>(v[r]), b = as_vec<double>(v[i]);
    cdvec z(a.size());
    for (size_t k = 0; k < a.size(); ++k) z[k] = cd(a[k], b[k]);
    return z;
  };
  const cdvec syms = cvec("syms_re", "syms_im");
  CHECK(dsss::dbpsk_bits(syms) == as_vec<uint8_t>(v["dbpsk"]));
  CHECK(dsss::diff_bits(syms, 0.013, 1) == as_vec<uint8_t>(v["diff1"]));
  CHECK(dsss::diff_bits(syms, -0.2, 2) == as_vec<uint8_t>(v["diff2"]));
  const cdvec sfd = cvec("sfd_re", "sfd_im");
  const Bits raw = dsss::dbpsk_bits(sfd);
  CHECK(raw == as_vec<uint8_t>(v["sfd_raw"]));
  Json sfds = Json::array();
  for (const auto& [e, p] : dsss::find_sfds(dsss::descramble(raw))) sfds.push_back(Json::array({e, p}));
  CHECK(sfds == v["sfds"]);
  CHECK(!sfds.empty());
  Bits desc(128, 1);
  for (int i = 0; i < 16; ++i) desc.push_back((0xF3A0 >> i) & 1);
  const auto rebuilt = dsss::rebuild_scrambled(Bits(raw.begin(), raw.begin() + 144), desc);
  CHECK(rebuilt.first == as_vec<uint8_t>(v["rebuilt"]));
  CHECK(rebuilt.second == v["agree"].get<int64_t>());
}

TEST_CASE("wifi dsss: synthetic captures (every rate, correction, excision) equal python scan()") {
  for (const char* name : {"rates", "weak"}) {
    CAPTURE(name);
    const Json v = load_json(std::string("scene_") + name + ".json");
    const cfvec x = load_cf32(std::string("scene_") + name + ".cf32");
    REQUIRE(x.size() == v["n"].get<size_t>());
    const double fs = v["fs"].get<double>(), fc = v["fc"].get<double>();
    const Channelizer chz(x, fs, fc);
    // detection on every channel
    for (auto it = v["candidates"].begin(); it != v["candidates"].end(); ++it) {
      const std::string ch = it.key();
      const Json& cands = it.value();
      CAPTURE(ch);
      const double f = 2407e6 + 5e6 * std::stoi(ch);
      const auto got = dsss::candidates(chz.slice(f, 11e6, dsss::kFs).y);
      REQUIRE(got.size() == cands.size());
      for (size_t i = 0; i < got.size(); ++i) {
        CHECK(got[i].a == cands[i][0].get<int64_t>());
        CHECK(got[i].b == cands[i][1].get<int64_t>());
        CHECK(got[i].rho == doctest::Approx(cands[i][2].get<double>()).epsilon(1e-12));
      }
    }
    // Barker and chip matched filters, tracking (channel 6)
    auto sl = chz.slice(2437e6, 11e6, dsss::kFsDec);
    const dsss::Chan chan(sl.y, sl.fs, 2437e6);
    const auto idx = as_vec<size_t>(v["mf_idx"]);
    for (size_t k = 0; k < idx.size(); ++k) {
      CHECK(chan.c[idx[k]].real() == doctest::Approx(v["c_re"][k].get<double>()).epsilon(1e-12));
      CHECK(chan.c[idx[k]].imag() == doctest::Approx(v["c_im"][k].get<double>()).epsilon(1e-12));
      CHECK(chan.ym[idx[k]].real() == doctest::Approx(v["ym_re"][k].get<double>()).epsilon(1e-12));
      CHECK(chan.ym[idx[k]].imag() == doctest::Approx(v["ym_im"][k].get<double>()).epsilon(1e-12));
    }
    const auto tr = chan.track(v["track_tau"].get<double>(), 60, std::nullopt, v["track_A"].get<double>());
    REQUIRE(tr.s.size() == v["track_s_re"].size());
    for (size_t k = 0; k < tr.s.size(); ++k) {
      CHECK(tr.taus[k] == doctest::Approx(v["track_taus"][k].get<double>()).epsilon(1e-12));
      CHECK(tr.s[k].real() == doctest::Approx(v["track_s_re"][k].get<double>()).epsilon(1e-9));
    }
    CHECK(tr.T == doctest::Approx(v["track_T"].get<double>()).epsilon(1e-12));
    // the whole receiver
    const auto res = dsss::scan(chz, fc, fs);
    CHECK(res.n_candidates == v["n_candidates"].get<size_t>());
    const Json got = res.ppdus, exp = v["ppdus"];
    CHECK_MESSAGE(diff(exp, got, "$", 1e-9) == "", diff(exp, got, "$", 1e-9));
    const Json ugot = res.unexplained;
    CHECK_MESSAGE(diff(v["unexplained"], ugot, "$", 1e-12) == "", diff(v["unexplained"], ugot, "$", 1e-12));
    size_t ok = 0;
    for (const auto& r : res.ppdus) ok += r.value("fcs_ok", false);
    CHECK(ok >= 5);
  }
}

TEST_CASE("wifi dsss: local oaconvolve is bit-identical to sig::oaconvolve") {
  cfvec x(50000);
  for (size_t i = 0; i < x.size(); ++i) x[i] = cf(std::sin(0.37f * i) + 0.01f * (i % 7), std::cos(0.11f * i * i));
  for (size_t taps : {4, 22, 44}) {
    dvec k(taps);
    for (size_t i = 0; i < taps; ++i) k[i] = dsss::kBarker[(10 - i * 11 / taps) % 11] / (1 + i % 3);
    const cdvec ref = sig::oaconvolve(x, k, np::ConvMode::Full);
    const cdvec got = dsss::detail::oaconvolve_cf_real(x, k, 2);
    REQUIRE(got.size() == ref.size());
    size_t bad = 0;
    for (size_t i = 0; i < ref.size(); ++i) bad += !(got[i] == ref[i]);
    CHECK(bad == 0);
  }
}
