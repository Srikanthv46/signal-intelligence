// Pipeline plugin: DJI OcuSync video / RC uplink / DroneID (python: receivers/dji/plugin.py).
#include <map>
#include <memory>

#include "core/drone_models.hpp"
#include "decode/emission_adapter.hpp"
#include "receivers/dji/detector.hpp"

namespace si::rx::dji {
namespace {

class DjiDecoder : public EmissionReceiver {
 public:
  std::string name() const override { return "dji"; }
  std::string package() const override { return "dji"; }
  std::string protocol() const override { return "DJI OcuSync/DroneID"; }
  std::string key() const override { return "dji"; }
  std::set<std::string> bands() const override { return {"ism_24", "ism_58", "unii_12", "unii_2e"}; }
  std::string method() const override { return "receivers.dji (research receiver)"; }
  std::vector<std::string> steps() const override {
    return {"15/30 kHz CP-OFDM fold on the 0.5 ms slot + DJI Zadoff-Chu root search (video downlink)",
            "RC uplink: 72-subcarrier ZC hops", "DroneID: 9-symbol 15.36 MSps burst, ZC sync, turbo decode + CRC"};
  }

  std::vector<Emission> run_detector(DecodeContext& ctx) const override {
    return detect(ctx.cap.iq(), ctx.cap.fs, ctx.cap.fc, true);
  }

  std::string pretty(const std::string& top, const std::vector<Emission>& em) const override {
    // python: "-like" when every emission of the group comes from a "-like" (inferred) link
    bool like = true;
    for (const auto& e : em) like = like && e.features.value("evidence_tier", std::string()) == "like";
    std::string s = top == "droneid" ? "DJI DroneID" : top == "uplink" ? "DJI OcuSync RC uplink" : "DJI OcuSync " + top;
    return like ? s + "-like" : s;
  }

  DecodeResult summarise(const std::vector<Emission>& em, const Signal* sig) const override {
    DecodeResult r = EmissionReceiver::summarise(em, sig);
    // research (PLAN D4): 9-symbol OcuSync 2/3 vs 8-symbol legacy DroneID frame (python: most common, then max)
    std::map<std::string, int> cnt;
    for (const auto& e : em)
      if (e.features.contains("droneid_frame") && e.features["droneid_frame"].is_string())
        ++cnt[e.features["droneid_frame"].get<std::string>()];
    if (!cnt.empty()) {
      auto best = cnt.begin();
      for (auto it = cnt.begin(); it != cnt.end(); ++it)
        if (it->second > best->second || (it->second == best->second && it->first > best->first)) best = it;
      r.fields["droneid_frame"] = best->first;
    }
    // model level (research, PLAN M1): the decoded product-type byte wins; the frame layout gives the generation
    const Json* dt = nullptr;
    for (const auto& e : em)
      if (e.features.contains("device_type") && !e.features["device_type"].is_null()) { dt = &e.features["device_type"]; break; }
    const bool name_given = dt && dt->is_string() &&
                            dt->get<std::string>().find_first_not_of(" 0123456789") != std::string::npos;
    if (name_given) {
      // python: the external decoder already maps the byte to a name; "Unknown" = unmapped
      const std::string nm = dt->get<std::string>();
      const bool known = nm != "Unknown";
      r.fields["model"] = drone_models::model_block(
          "DJI", known ? std::vector<std::string>{nm} : std::vector<std::string>{},
          {"DroneID product type: " + nm + " (the decoder's own table)"}, "decoded",
          known ? "" : "product-type code not in the decoder's table");
    } else if (dt) {
      const std::string code = dt->is_string() ? dt->get<std::string>() : dt->dump();
      const std::string name = dt->is_number_integer() ? drone_models::dji_product(dt->get<int64_t>()) : "";
      const std::string ev = "DroneID product type " + code +
                             (name.empty() ? std::string(" (code not in the public table)")
                                           : " = " + name + " (" + drone_models::dji_product_source(dt->get<int64_t>()) + ")");
      r.fields["model"] = drone_models::model_block("DJI", name.empty() ? std::vector<std::string>{} : std::vector<std::string>{name},
                                                    {ev}, "decoded", name.empty() ? "product-type code not in the public table" : "");
    } else if (!cnt.empty()) {
      // file 13 §4 (NDSS 2023, RUB legacy flag): the 8-symbol frame is the Mavic Pro / Mavic 2 family
      const std::string fr = r.fields["droneid_frame"].get<std::string>();
      const bool legacy = fr.rfind("8-symbol", 0) == 0;
      const std::vector<std::string> cands =
          legacy ? std::vector<std::string>{"Mavic Pro", "Mavic Pro Platinum", "Mavic 2 Pro", "Mavic 2 Zoom",
                                            "Mavic 2 Enterprise", "Mavic 2 Enterprise Advanced", "Mini 4 Pro"}
                 : std::vector<std::string>{};
      r.fields["model"] = drone_models::model_block(
          "DJI", cands, {"DroneID " + fr}, "inferred",
          legacy ? "the 8-symbol frame is sent by the Mavic Pro / Mavic 2 family and by O4 aircraft (Mini 4 Pro measured; "
                   "other O4 models not measured)"
                 : "the DroneID frame layout gives the generation only (payload not decoded)");
    }
    return r;
  }

  // python: f"OFDM{f' {scs} kHz SCS' if scs else ''} (ZC sync)"
  std::string modulation(const std::vector<Emission>& em) const override {
    const Json& f = em[0].features;
    std::string s;
    if (f.contains("scs_khz") && f["scs_khz"].is_number() && f["scs_khz"].get<double>() != 0)
      s = " " + std::to_string(f["scs_khz"].get<int64_t>()) + " kHz SCS";
    return "OFDM" + s + " (ZC sync)";
  }
};

}  // namespace

void register_plugins(si::PluginList& out) { out.push_back(std::make_unique<DjiDecoder>()); }

}  // namespace si::rx::dji
