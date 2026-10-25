// Pipeline plugin: AIS vessel / station reports (research; python: receivers/marine/plugin.py).
#include <algorithm>
#include <map>
#include <memory>

#include "core/pyfmt.hpp"
#include "decode/emission_adapter.hpp"
#include "receivers/marine/detector.hpp"

namespace si::rx::marine {

namespace {

const char* const kMerge[] = {"flag", "mid", "kind", "nav_status", "sog_kn", "lat", "lon", "cog_deg",
                              "heading_deg", "name", "callsign", "ship_type", "destination", "imo",
                              "length_m", "beam_m", "draught_m", "aid_type", "utc"};

class MarineDecoder : public EmissionReceiver {
 public:
  std::string name() const override { return "marine"; }
  std::string package() const override { return "marine"; }
  std::string protocol() const override { return "AIS"; }
  std::string key() const override { return "marine"; }
  std::set<std::string> bands() const override { return {"marine_vhf"}; }
  std::string method() const override { return "receivers.marine: AIS GMSK / HDLC / CRC-16 X.25 (research)"; }
  bool owns_band() const override { return true; }
  std::vector<std::string> steps() const override {
    return {"AIS 1 / 2 (161.975 / 162.025 MHz, receiver ppm applied): +-10 kHz at 48 kSps, FM discriminator, "
            "offset removed (32-bit moving average), integrated over one bit",
            "per sample phase: NRZ levels -> NRZI decoding -> training + flag 0x7E, de-stuffing to the closing "
            "flag, octets LSB first, CRC-16/X.25 residue 0xF0B8",
            "message fields (ITU-R M.1371-5): types 1/2/3 position, 4 base station, 5 static / voyage, 18 class B, "
            "21 aid to navigation, 24 class B static; MID 419 = India"};
  }
  std::vector<Emission> run_detector(DecodeContext& ctx) const override {
    const Json* cal = ctx.cache.get<Json>("ppm");
    if (!cal || cal->is_null() || cal->empty()) cal = ctx.cache.get<Json>("sweep_ppm");
    double ppm = 0.0;
    if (cal && cal->is_object() && cal->contains("ppm") && (*cal)["ppm"].is_number()) ppm = (*cal)["ppm"].get<double>();
    return detect(ctx.ch, ppm);
  }
  std::string group_key(const Emission&) const override { return "ais"; }

  DecodeResult summarise(const std::vector<Emission>& em0, const Signal*) const override {
    std::vector<Emission> em = em0;
    std::stable_sort(em.begin(), em.end(), [](const Emission& a, const Emission& b) { return a.start_us < b.start_us; });
    std::vector<std::string> order;
    std::map<std::string, Json> vessels;
    Json msgs = Json::array();
    for (const auto& e : em) {
      const Json& d = e.decoded;
      Json m = Json::object();
      m["t_s"] = py_round(e.start_us * 1e-6, 5);
      m["channel"] = e.features["channel"];
      for (auto it = d.begin(); it != d.end(); ++it)
        if (it.key() != "crc_ok") m[it.key()] = it.value();
      msgs.push_back(m);
      const std::string mm = d["mmsi"].get<std::string>();
      if (!vessels.count(mm)) {
        vessels[mm] = Json{{"mmsi", mm}, {"messages", 0}, {"msg_types", Json::array()}};
        order.push_back(mm);
      }
      Json& v = vessels[mm];
      v["messages"] = v["messages"].get<int64_t>() + 1;
      if (std::find(v["msg_types"].begin(), v["msg_types"].end(), d["msg_type"]) == v["msg_types"].end())
        v["msg_types"].push_back(d["msg_type"]);
      for (const char* k : kMerge)
        if (d.contains(k) && !d[k].is_null()) v[k] = d[k];
    }
    Json vs = Json::array(), names = Json::array();
    for (const auto& k : order) {
      vs.push_back(vessels[k]);
      names.push_back(vessels[k].contains("name") && truthy_str(vessels[k]["name"]) ? vessels[k]["name"] : vessels[k]["mmsi"]);
    }
    Json fields = Json::object();
    fields["n_messages"] = static_cast<int64_t>(msgs.size());
    fields["vessels"] = vs;
    Json top = Json::array();
    for (size_t i = 0; i < msgs.size() && i < 100; ++i) top.push_back(msgs[i]);
    fields["messages"] = top;
    double lo = em[0].f_lo_mhz, hi = em[0].f_hi_mhz;
    for (const auto& x : em) {
      lo = std::min(lo, x.f_lo_mhz);
      hi = std::max(hi, x.f_hi_mhz);
    }
    const size_t n = vs.size();
    return result("success", method(), steps(), fields, "every message verified by its CRC-16/X.25",
                  fmt("AIS (%zu station%s)", n, n != 1 ? "s" : ""), 0.97, std::string("GMSK 9600 bit/s (AIS)"),
                  fmt("%zu CRC-valid AIS message(s) from ", msgs.size()) + py_repr(names),
                  std::make_pair(lo * 1e6, hi * 1e6));
  }

 private:
  static bool truthy_str(const Json& j) { return j.is_string() && !j.get<std::string>().empty(); }
};

}  // namespace

void register_plugins(si::PluginList& out) { out.push_back(std::make_unique<MarineDecoder>()); }

}  // namespace si::rx::marine
