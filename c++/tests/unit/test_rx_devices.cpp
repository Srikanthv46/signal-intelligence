// receivers/devices against python receivers.devices (tests/vectors/devices/gen_devices.py):
// synthetic decoded 802.11 signals covering every tracker path and threat rule -> report extras
// (devices, findings, device_summary, device_groups), per-signal device_macs, finalize_report().
#include <doctest/doctest.h>

#include <fstream>
#include <string>

#include "decode/registry.hpp"
#include "receivers/devices/devices.hpp"

using namespace si;

namespace {

const Json& cases() {
  static const Json r = [] {
    std::ifstream f(std::string(SI_TESTVEC) + "/devices/cases.json");
    REQUIRE(f.good());
    return Json::parse(f);
  }();
  return r;
}

std::vector<Signal> signals_of(const Json& c) {
  std::vector<Signal> v;
  for (const auto& s : c["signals"]) {
    Signal x;
    x.decoded.method = s["method"].get<std::string>();
    x.decoded.fields = s["fields"];
    v.push_back(x);
  }
  return v;
}

// Json equality plus identical serialisation (catches 24 vs 24.0 in keys and strings)
void same(const Json& a, const Json& b, const std::string& what) {
  INFO(what);
  CHECK(a == b);
  if (a != b) {
    const std::string sa = a.dump(), sb = b.dump();
    size_t i = 0;
    while (i < sa.size() && i < sb.size() && sa[i] == sb[i]) ++i;
    MESSAGE("C++    ..." << sa.substr(i > 80 ? i - 80 : 0, 240));
    MESSAGE("python ..." << sb.substr(i > 80 ? i - 80 : 0, 240));
  }
}

}  // namespace

TEST_CASE("devices: tracker, links and findings match python") {
  for (auto it = cases().begin(); it != cases().end(); ++it) {
    const Json& c = it.value();
    std::vector<Signal> sigs = signals_of(c);
    const Json extras = rx::devices::device_report(sigs);
    same(extras, c["extras"], it.key() + " extras");
    for (size_t i = 0; i < sigs.size(); ++i) {
      const Json& want = c["device_macs"][i];
      const Json got = sigs[i].decoded.fields.contains("device_macs") ? sigs[i].decoded.fields["device_macs"]
                                                                       : Json(nullptr);
      same(got, want, it.key() + " device_macs " + std::to_string(i));
    }
    // finalize_report through the registered plugin
    Json rep = Json::object();
    Json ss = Json::array();
    for (size_t i = 0; i < sigs.size(); ++i) {
      Json f = Json::object();
      if (sigs[i].decoded.fields.contains("device_macs")) f["device_macs"] = sigs[i].decoded.fields["device_macs"];
      ss.push_back({{"id", c["signals"][i]["id"]}, {"decoded", {{"fields", f}}}});
    }
    rep["signals"] = ss;
    for (auto e = extras.begin(); e != extras.end(); ++e) rep[e.key()] = e.value();
    auto plugs = plugins({"devices"});
    REQUIRE(plugs.size() == 1);
    CHECK(plugs[0]->name() == "devices");
    CHECK(plugs[0]->capture_level());
    plugs[0]->finalize_report(rep);
    same(rep, c["finalized"], it.key() + " finalized");
  }
}

TEST_CASE("devices: a signal without a decoder method fails like python (IndexError)") {
  std::vector<Signal> sigs(1);
  CHECK_THROWS_WITH(rx::devices::device_report(sigs), "IndexError: list index out of range");
}
