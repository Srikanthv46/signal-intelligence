// PLAN G10 / G12 / G14 (research) vs python (tests/vectors/g10_g14/gen_g10_g14.py): HoTT / Hitec / S-FHSS and
// ADS-L / OGN through rc_link detect + plugin fields + model blocks; FANET through the LoRa identify(); Reed-Solomon;
// UAT demodulation and message fields; China GB 46750 Remote ID.
#include <doctest/doctest.h>

#include <fstream>
#include <string>

#include "core/remote_id.hpp"
#include "core/rs.hpp"
#include "decode/emission_adapter.hpp"
#include "decode/registry.hpp"
#include "receivers/adsb/uat.hpp"
#include "receivers/rc_link/detector.hpp"
#include "receivers/rc_link/lora_phy.hpp"

using namespace si;

namespace {
Json ref() {
  std::ifstream jf(std::string(SI_TESTVEC) + "/g10_g14/cases.json");
  REQUIRE(jf);
  return Json::parse(jf);
}

cfvec load(const std::string& name) {
  std::ifstream f(std::string(SI_TESTVEC) + "/g10_g14/" + name, std::ios::binary | std::ios::ate);
  REQUIRE(f);
  const size_t n = static_cast<size_t>(f.tellg()) / sizeof(cf);
  cfvec x(n);
  f.seekg(0);
  f.read(reinterpret_cast<char*>(x.data()), static_cast<std::streamsize>(n * sizeof(cf)));
  return x;
}

std::vector<uint8_t> unhex(const std::string& h) {
  std::vector<uint8_t> b;
  for (size_t i = 0; i + 1 < h.size(); i += 2) b.push_back(static_cast<uint8_t>(std::stoul(h.substr(i, 2), nullptr, 16)));
  return b;
}

void check_capture(const Json& c, const std::string& file, const std::vector<std::string>& labs) {
  const cfvec x = load(file);
  Json labels = Json::array(), got = Json::array();
  std::vector<Emission> em;
  for (const auto& e : rx::rc_link::detect(x, c["fs"], c["fc"]))
    for (const auto& l : labs)
      if (e.label == l) {
        labels.push_back(e.label), got.push_back(e.decoded), em.push_back(e);
        break;
      }
  CHECK(labels.dump() == c["labels"].dump());
  CHECK(got.dump() == c["decoded"].dump());
  const PluginList pl = plugins({"rc_link"});
  const EmissionReceiver* rc = nullptr;
  for (const auto& p : pl)
    if (p->name() == "rc_link") rc = dynamic_cast<const EmissionReceiver*>(p.get());
  REQUIRE(rc);
  for (const auto& [lab, want] : c["fields"].items()) {
    std::vector<Emission> g;
    for (const auto& e : em)
      if (e.label == lab) g.push_back(e);
    const auto r = rc->summarise(g, nullptr);
    for (const auto& [k, v] : want.items()) {
      CAPTURE(lab);
      CAPTURE(k);
      CHECK(r.fields[k].dump() == v.dump());
    }
  }
}
}  // namespace

TEST_CASE("g12: HoTT / Hitec / S-FHSS through rc_link == python") { check_capture(ref()["c24"], "c24.cf32", {"hott", "hitec", "sfhss"}); }

TEST_CASE("g10: ADS-L / OGN through rc_link == python") { check_capture(ref()["c868"], "c868.cf32", {"adsl", "ogn_like"}); }

TEST_CASE("g10: FANET through the LoRa identify() == python") {
  for (const auto& c : ref()["fanet"]) {
    const auto d = rx::rc_link::lora_phy::identify(c[0].get<int>(), unhex(c[1].get<std::string>()));
    CAPTURE(c[1].get<std::string>());
    CHECK((d ? *d : Json(nullptr)).dump() == c[2].dump());
  }
}

TEST_CASE("g14: Reed-Solomon encode / decode == python") {
  for (const auto& c : ref()["rs"]) {
    const rs::RS code(0x187, 120, 1, c[0].get<int>());
    const auto data = unhex(c[1].get<std::string>()), cw = unhex(c[2].get<std::string>());
    auto full = data;
    const auto par = code.encode(data);
    full.insert(full.end(), par.begin(), par.end());
    CHECK(std::vector<uint8_t>(cw.size() == full.size() ? full : cw) == full);
    const auto d = code.decode(cw);
    if (c[3].is_null()) {
      CHECK(!d);
    } else {
      REQUIRE(d);
      CHECK(d->first == unhex(c[3].get<std::string>()));
      CHECK(d->second == c[4].get<int>());
    }
  }
}

TEST_CASE("g14: UAT demodulation and fields == python") {
  const Json r = ref()["uat"];
  CHECK(Json(rx::adsb::uat::demod(load("uat.cf32"))).dump() == r["demod"].dump());
  for (const auto& c : r["parse"]) CHECK(rx::adsb::uat::parse(unhex(c[0].get<std::string>())).dump() == c[1].dump());
}

TEST_CASE("g14: China GB 46750 Remote ID == python") {
  const Json r = ref()["gb46750"];
  auto j = [](const std::optional<Json>& d) { return d ? *d : Json(nullptr); };
  CHECK(j(remote_id::parse_odid(unhex(r[0][0].get<std::string>()))).dump() == r[0][1].dump());
  CHECK(j(remote_id::parse_gb46750(unhex(r[1][0].get<std::string>()))).dump() == r[1][1].dump());
  CHECK(j(remote_id::parse_gb46750(unhex(r[2][0].get<std::string>()))).dump() == r[2][1].dump());
}
