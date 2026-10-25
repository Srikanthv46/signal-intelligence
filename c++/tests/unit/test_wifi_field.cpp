// PLAN batch W (Wi-Fi field gaps) against python (tests/vectors/wifi_field/gen_wifi_field.py): HE-SIG-A (W2),
// drone / DIY hotspot SSIDs (W3), Tello SDK and video streams on open links (W4), link profiles + drone-link
// verdicts (W1).
#include <doctest/doctest.h>

#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "core/signal.hpp"
#include "receivers/wifi/field.hpp"
#include "receivers/wifi/mac80211.hpp"
#include "receivers/wifi/ofdm_report.hpp"
#include "receivers/wifi/ofdm_rx.hpp"

using namespace si;
using namespace si::rx::wifi;

namespace {
const Json& ref() {
  static const Json r = [] {
    std::ifstream f(std::string(SI_TESTVEC) + "/wifi_field/ref.json");
    REQUIRE(f.good());
    return Json::parse(f);
  }();
  return r;
}
cfvec load_cf32(const std::string& name) {
  std::ifstream f(std::string(SI_TESTVEC) + "/wifi_field/" + name, std::ios::binary | std::ios::ate);
  REQUIRE(f.good());
  cfvec v(static_cast<size_t>(f.tellg()) / sizeof(cf));
  f.seekg(0);
  f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(cf)));
  return v;
}
std::vector<uint8_t> unhex(const std::string& h) {
  std::vector<uint8_t> v;
  for (size_t i = 0; i + 1 < h.size(); i += 2) v.push_back(static_cast<uint8_t>(std::stoi(h.substr(i, 2), nullptr, 16)));
  return v;
}
Json opt(const std::optional<int>& v) { return v ? Json(*v) : Json(nullptr); }

// Python value == C++ value: numbers within rel, everything else equal (object keys compared both ways)
void same(const Json& py, const Json& cc, const std::string& path, double rel = 1e-9) {
  if (py.is_number() && cc.is_number()) {
    const double a = py.get<double>(), b = cc.get<double>();
    CHECK_MESSAGE(std::fabs(a - b) <= rel * std::fabs(a) + 1e-12, path << ": " << a << " vs " << b);
  } else if (py.is_object() && cc.is_object()) {
    CHECK_MESSAGE(py.size() == cc.size(), path << ": keys " << py.dump() << " vs " << cc.dump());
    for (const auto& [k, v] : py.items()) {
      CHECK_MESSAGE(cc.contains(k), path << "." << k << " missing");
      if (cc.contains(k)) same(v, cc[k], path + "." + k, rel);
    }
  } else if (py.is_array() && cc.is_array()) {
    REQUIRE_MESSAGE(py.size() == cc.size(), path << ": size " << py.size() << " vs " << cc.size());
    for (size_t i = 0; i < py.size(); ++i) same(py[i], cc[i], path + "[" + std::to_string(i) + "]", rel);
  } else {
    CHECK_MESSAGE(py == cc, path << ": " << py.dump() << " vs " << cc.dump());
  }
}
}  // namespace

TEST_CASE("wifi field W2: HE-SIG-A parse == python") {
  for (const auto& c : ref()["he_parse"]) {
    const auto h = ofdm::parse_he_sig_a(c["bits"].get<std::vector<uint8_t>>(), c["lsig_length"].get<int>());
    const Json& o = c["out"];
    REQUIRE(h.has_value() == !o.is_null());
    if (!h) continue;
    CHECK(h->format == o["he_format"].get<std::string>());
    CHECK(h->uplink == o["uplink"].get<bool>());
    CHECK(h->bss_color == o["bss_color"].get<int>());
    CHECK(h->bw == o["bw"].get<int>());
    CHECK(opt(h->txop_us) == o["txop_us"]);
    CHECK(opt(h->mcs) == (o.contains("mcs") ? o["mcs"] : Json(nullptr)));
    CHECK(opt(h->nsts) == (o.contains("nsts") ? o["nsts"] : Json(nullptr)));
    if (o.contains("ldpc")) CHECK(h->ldpc == o["ldpc"].get<int>());
    if (o.contains("dcm")) CHECK(h->dcm == o["dcm"].get<int>());
  }
}

TEST_CASE("wifi field W2: HE scene, detector, records and summary == python") {
  const Json& j = ref()["he_scene"];
  const cfvec x = load_cf32("he_scene.cf32");
  const double fs = j["fs"], fc = j["fc"];
  ofdm::WifiDetector det(fs, fc);
  const auto pk = det.run(x, 2);
  REQUIRE(pk.size() == j["packets"].size());
  for (size_t i = 0; i < pk.size(); ++i) {
    const auto& p = pk[i];
    Json c = {{"fmt", p.fmt}, {"start_sample", p.start_sample}, {"lsig_length", p.lsig_length},
              {"mcs", opt(p.mcs)}, {"sig_crc_ok", p.sig_crc_ok ? Json(*p.sig_crc_ok) : Json(nullptr)},
              {"bandwidth_mhz", p.bandwidth_mhz}, {"f_lo_mhz", p.f_lo_mhz}, {"f_hi_mhz", p.f_hi_mhz}};
    same(j["packets"][i], c, "packet" + std::to_string(i));
  }
  auto recs = ofdm_report::records(pk, fs, 20e6, false, 2);
  ofdm_report::attribute_vht(recs);
  same(j["records"], Json(recs), "records", 1e-6);
  std::vector<Json*> mine;
  for (auto& r : recs) mine.push_back(&r);
  const auto s = ofdm_report::summarise_ofdm(mine, recs.size(), 0, &recs);
  Json he = Json::object();
  for (const auto& [k, v] : s.fields.items())
    if (k.rfind("he_", 0) == 0) he[k] = v;
  same(j["summary"], he, "summary");
}

TEST_CASE("wifi field W3 / W4: SSID signatures, Tello SDK, video streams == python") {
  for (const auto& [k, c] : ref()["mpdu"].items()) {
    CAPTURE(k);
    const Json got = rx::wifi::mac80211::parse(unhex(c["hex"].get<std::string>()));
    for (const char* key : {"ipv4", "drone", "summary", "ssid"}) {
      CAPTURE(key);
      CHECK((got.contains(key) ? got[key] : Json(nullptr)).dump() == c["out"][key].dump());
    }
  }
}

TEST_CASE("wifi field W1: link profiles and drone-link verdicts == python") {
  for (const auto& [k, c] : ref()["links"].items()) {
    CAPTURE(k);
    const auto frames = c["frames"].get<std::vector<Json>>();
    same(c["out"], rx::wifi::field::links(frames, c["block_s"].get<double>()), k);
  }
}

TEST_CASE("wifi field W1: drone-style link end to end (records -> summary) == python") {
  const Json& j = ref()["link_scene"];
  const cfvec x = load_cf32("link_scene.cf32");
  const double fs = j["fs"], fc = j["fc"];
  ofdm::WifiDetector det(fs, fc);
  auto recs = ofdm_report::records(det.run(x, 2), fs, 20e6, false, 2);
  ofdm_report::attribute_vht(recs);
  const std::vector<Json> all = recs;
  std::vector<Json*> mine;
  for (auto& r : recs) mine.push_back(&r);
  const auto s = ofdm_report::summarise_ofdm(mine, recs.size(), 0, &all, j["block_s"].get<double>());
  CHECK(*s.label == j["label"].get<std::string>());
  CHECK(*s.evidence == j["evidence"].get<std::string>());
  same(j["wifi_links"], s.fields["wifi_links"], "wifi_links", 1e-6);
  same(j["drone_links"], s.fields["drone_links"], "drone_links", 1e-6);
}
