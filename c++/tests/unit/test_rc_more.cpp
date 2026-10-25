// PLAN G3-G6 RC decoders (rc_link radiolink / AFHDS 1 / hubsan / toys / ELRS GFSK; research) vs python on a synthetic
// capture (tests/vectors/rc_more/gen_rc_more.py): hop tables, detect's emissions, plugin fields and model blocks.
#include <doctest/doctest.h>

#include <fstream>
#include <string>

#include "decode/emission_adapter.hpp"
#include "decode/registry.hpp"
#include "receivers/rc_link/afhds2a.hpp"
#include "receivers/rc_link/detector.hpp"
#include "receivers/rc_link/radiolink.hpp"

using namespace si;

TEST_CASE("rc_more: hop tables and the rc_link pipeline == python") {
  std::ifstream jf(std::string(SI_TESTVEC) + "/rc_more/cases.json");
  REQUIRE(jf);
  const Json ref = Json::parse(jf);
  for (const auto& c : ref["radiolink_hops"]) {
    const std::string h = c[0].get<std::string>();
    uint8_t id[4];
    for (int i = 0; i < 4; ++i) id[i] = static_cast<uint8_t>(std::stoul(h.substr(2 * i, 2), nullptr, 16));
    const auto [hf, rep] = rx::rc_link::radiolink::hop_table(id);
    CHECK(Json(hf) == c[1]);
    CHECK(rep == c[2].get<int>());
  }
  for (const auto& c : ref["afhds1_hops"]) CHECK(Json(rx::rc_link::afhds2a::afhds1_hops(c[0].get<unsigned>())) == c[1]);
  std::ifstream f(std::string(SI_TESTVEC) + "/rc_more/capture.cf32", std::ios::binary | std::ios::ate);
  REQUIRE(f);
  const size_t n = static_cast<size_t>(f.tellg()) / sizeof(cf);
  cfvec x(n);
  f.seekg(0);
  f.read(reinterpret_cast<char*>(x.data()), static_cast<std::streamsize>(n * sizeof(cf)));
  Json labels = Json::array(), got = Json::array();
  std::vector<Emission> em;
  for (const auto& e : rx::rc_link::detect(x, ref["fs"], ref["fc"]))
    for (const char* l : {"radiolink", "afhds1", "hubsan", "toy_drone", "elrs_gfsk"})
      if (e.label.rfind(l, 0) == 0) {
        labels.push_back(e.label), got.push_back(e.decoded), em.push_back(e);
        break;
      }
  CHECK(labels.dump() == ref["labels"].dump());
  CHECK(got.dump() == ref["decoded"].dump());
  const PluginList pl = plugins({"rc_link"});
  const EmissionReceiver* rc = nullptr;
  for (const auto& p : pl)
    if (p->name() == "rc_link") rc = dynamic_cast<const EmissionReceiver*>(p.get());
  REQUIRE(rc);
  for (const auto& [lab, want] : ref["fields"].items()) {
    std::vector<Emission> g;
    for (const auto& e : em)
      if (e.label.rfind(lab, 0) == 0) g.push_back(e);
    const auto r = rc->summarise(g, nullptr);
    for (const auto& [k, v] : want.items()) {
      CAPTURE(lab);
      CAPTURE(k);
      CHECK(r.fields[k].dump() == v.dump());
    }
  }
}
