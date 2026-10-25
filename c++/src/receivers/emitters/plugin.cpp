// Pipeline plugin: drone emitter association (python: receivers/emitters/plugin.py). Which
// receiver-defined signals are one radio and which radios form one drone system. Adds
// report["emitters"] and signal.decoded.fields.emitter_ids; reads only Signal::emissions.
#include <algorithm>
#include <memory>

#include "decode/plugin.hpp"
#include "receivers/emitters/association.hpp"
#include "receivers/emitters/families.hpp"

namespace si::rx::emitters {
namespace {

class EmitterAssociation : public Decoder {
 public:
  std::string name() const override { return "emitters"; }
  std::string package() const override { return "emitters"; }
  std::string protocol() const override { return "drone emitter association"; }

  Json postprocess(DecodeContext& ctx, std::vector<Signal>& signals) const override {
    std::vector<std::pair<size_t, size_t>> ref;   // (signal, emission) of every drone-family emission
    for (size_t s = 0; s < signals.size(); ++s)
      for (size_t j = 0; j < signals[s].emissions.size(); ++j)
        if (DRONE_FAMILIES().count(signals[s].emissions[j].family)) ref.emplace_back(s, j);
    if (ref.empty()) return Json{{"emitters", Json::array()}};
    auto em_of = [&](const std::pair<size_t, size_t>& r) -> Emission& { return signals[r.first].emissions[r.second]; };
    std::stable_sort(ref.begin(), ref.end(), [&](const auto& a, const auto& b) {
      const Emission &x = em_of(a), &y = em_of(b);
      return x.start_us < y.start_us || (x.start_us == y.start_us && x.center_mhz < y.center_mhz);
    });
    std::vector<Emission> em;
    for (const auto& r : ref) em.push_back(em_of(r));
    measure_fingerprints(ctx.cap.iq(), ctx.cap.fs, ctx.cap.fc, em);
    for (size_t k = 0; k < ref.size(); ++k) em_of(ref[k]) = em[k];   // python refines them in place
    const auto emitters = associate(em, &assoc_models());
    Json out = Json::array();
    for (size_t k = 0; k < emitters.size(); ++k) {
      const Emitter& m = emitters[k];
      std::vector<size_t> sigs;
      for (const int64_t i : m.members)
        if (std::find(sigs.begin(), sigs.end(), ref[i].first) == sigs.end()) sigs.push_back(ref[i].first);
      const std::string eid = "E" + std::to_string(k);
      for (const size_t s : sigs) {   // one signal can hold several radios' bursts (Herelink control + video)
        Json& fields = signals[s].decoded.fields;
        if (!fields.contains("emitter_ids")) fields["emitter_ids"] = Json::array();
        Json& ids = fields["emitter_ids"];
        if (std::find(ids.begin(), ids.end(), Json(eid)) == ids.end()) ids.push_back(eid);
      }
      Json d = m.to_dict();
      d["emitter_id"] = eid;
      const auto it = DRONE_FAMILIES().find(m.family);
      d["operator_meaning"] = it == DRONE_FAMILIES().end() ? m.family : it->second;
      d["n_bursts"] = d["members"].size();
      d.erase("members");
      out.push_back(d);
    }
    return Json{{"emitters", out}};
  }

  void finalize_report(Json& rep) const override {
    if (!rep.contains("emitters")) return;
    for (Json& em : rep["emitters"]) {
      Json ids = Json::array();
      for (const Json& s : rep["signals"]) {
        const Json& f = s["decoded"]["fields"];
        if (!f.contains("emitter_ids")) continue;
        const Json& e = f["emitter_ids"];
        if (std::find(e.begin(), e.end(), em["emitter_id"]) != e.end()) ids.push_back(s["id"]);
      }
      em["signal_ids"] = ids;
    }
  }
};

}  // namespace

void register_plugins(PluginList& out) { out.push_back(std::make_unique<EmitterAssociation>()); }

}  // namespace si::rx::emitters
