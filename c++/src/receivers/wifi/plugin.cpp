// Pipeline plugins for Wi-Fi (python receivers/wifi/plugin.py: PLUGINS = _OFDM + _DSSS):
// wifi_ofdm, wifi_ofdm_scaled (plugin_ofdm.py) and wifi_dsss (plugin_dsss.py), in that order.
#include <algorithm>
#include <memory>
#include <set>

#include "core/bands_india.hpp"
#include "core/npcompat.hpp"
#include "decode/plugin.hpp"
#include "receivers/wifi/ofdm_report.hpp"
#include "receivers/wifi/ofdm_rx.hpp"
#include "receivers/wifi/ofdm_scaled.hpp"
#include "receivers/wifi/dsss.hpp"
#include "receivers/wifi/dsss_report.hpp"

namespace si::rx::wifi {
namespace {

DecodeResult to_result(const dsss::Summary& s) {
  return result(s.status, s.method, s.steps, s.fields, s.reason, s.label, s.confidence, s.modulation, s.evidence,
                s.span);
}

// ------------------------------------------------------------------ 802.11 OFDM
const std::set<std::string>& wifi_bands() {
  static const std::set<std::string> b = {"ism_24", "unii_12", "unii_2e", "ism_58", "its_59", "wifi6e"};
  return b;
}
bool in_wifi_bands(double lo, double hi) {
  for (const auto& k : bands::band_keys(lo, hi))
    if (wifi_bands().count(k)) return true;
  return false;
}

// python plugin_ofdm._scan / _scan_scaled cache entries (records are mutated by links.attribute,
// as the Python dicts are)
struct OfdmScan {
  std::vector<Json> recs;
  size_t n_rejects = 0;
};
struct ScaledScan {
  std::vector<std::pair<scaled::Link, std::vector<Json>>> links;
};

double block_s(const DecodeContext& ctx) { return static_cast<double>(ctx.cap.iq().size()) / ctx.cap.fs; }

double median_snr(const std::vector<Json*>& g) {
  dvec v;
  for (const Json* p : g) v.push_back((*p)["snr_db"].get<double>());
  return np::median(v);
}

std::vector<std::pair<double, double>> bursts_of(const std::vector<Json*>& g) {
  std::vector<std::pair<double, double>> b;
  for (const Json* p : g) b.push_back({(*p)["t0"].get<double>(), (*p)["t1"].get<double>()});
  return b;
}

// python plugin_ofdm.WifiOfdmDecoder
class WifiOfdmDecoder : public Decoder {
 public:
  std::string name() const override { return "wifi_ofdm"; }
  std::string package() const override { return "wifi"; }
  std::string protocol() const override { return "IEEE 802.11 OFDM"; }
  bool capture_level() const override { return true; }

  bool can_decode(const Signal& sig, DecodeContext&) const override {
    return in_wifi_bands(sig.f0, sig.f1) && sig.bw() >= 2e6;
  }

  std::vector<Signal> capture_signals(DecodeContext& ctx) const override {
    OfdmScan& sc = scan(ctx);
    std::vector<std::pair<std::pair<std::string, std::string>, std::vector<Json*>>> groups;
    for (auto& p : sc.recs) {
      const auto key = std::make_pair(p["channel"].get<std::string>(), p["bw_mhz"].dump());
      auto it = std::find_if(groups.begin(), groups.end(), [&](const auto& g) { return g.first == key; });
      if (it == groups.end()) groups.push_back({key, {&p}});
      else it->second.push_back(&p);
    }
    std::vector<Signal> out;
    for (const auto& [key, g] : groups) {
      const auto res = ofdm_report::summarise_ofdm(g, sc.recs.size(), 0, &sc.recs, block_s(ctx));
      out.push_back(receiver_signal(bursts_of(g), res.span->first, res.span->second, to_result(res), median_snr(g),
                                    name()));
    }
    return out;
  }

  DecodeResult decode(const Signal& sig, DecodeContext& ctx) const override {
    OfdmScan& sc = scan(ctx);
    std::vector<Json*> mine;
    for (auto& p : sc.recs)
      if (overlaps(sig, p["t0"].get<double>(), p["t1"].get<double>(), p["f0"].get<double>(), p["f1"].get<double>()))
        mine.push_back(&p);
    return to_result(ofdm_report::summarise_ofdm(mine, sc.recs.size(), sc.n_rejects, &sc.recs, block_s(ctx)));
  }

 private:
  static OfdmScan& scan(DecodeContext& ctx) {
    using Ptr = std::shared_ptr<OfdmScan>;
    if (Ptr* p = ctx.cache.get_mut<Ptr>("wifi_ofdm")) return **p;
    ofdm::WifiDetector det(ctx.cap.fs, ctx.cap.fc);
    auto res = std::make_shared<OfdmScan>();
    {
      std::vector<ofdm::Packet> pk;
      if (in_wifi_bands(ctx.cap.fc - ctx.cap.fs / 2, ctx.cap.fc + ctx.cap.fs / 2)) pk = det.run(ctx.cap.iq());
      res->n_rejects = det.n_rejects();
      // 802.11 off the standard grid / outside the Wi-Fi bands: blind 1 MHz centres only where a 3.2 us CP
      // correlation remains after the known PPDUs are blanked (python plugin_ofdm._scan)
      auto extra = ofdm::offgrid_plan(ctx.cap.iq(), ctx.cap.fs, ctx.cap.fc, pk, 60e6, &ctx.ch.X());
      if (!extra.empty()) {
        ofdm::WifiDetector d2(ctx.cap.fs, ctx.cap.fc, {}, extra);
        pk = ofdm::merge_packets(std::move(pk), d2.run(ctx.cap.iq()));
        res->n_rejects += d2.n_rejects();
      }
      res->recs = ofdm_report::records(pk, ctx.cap.fs, 20e6, false);
    }
    ofdm_report::attribute_vht(res->recs);
    ctx.cache.set<Ptr>("wifi_ofdm", res);
    return *res;
  }
};

// python plugin_ofdm.WifiOfdmScaledDecoder
class WifiOfdmScaledDecoder : public Decoder {
 public:
  std::string name() const override { return "wifi_ofdm_scaled"; }
  std::string package() const override { return "wifi"; }
  std::string protocol() const override { return "IEEE 802.11 OFDM, non-standard numerology"; }
  bool capture_level() const override { return true; }

  std::vector<Signal> capture_signals(DecodeContext& ctx) const override {
    ScaledScan& sc = scan(ctx);
    std::vector<Signal> out;
    for (auto& [ln, recs] : sc.links) {
      if (recs.empty()) continue;
      std::vector<Json*> g;
      for (auto& r : recs) g.push_back(&r);
      const auto res = ofdm_report::summarise_scaled(ln, g, block_s(ctx));
      out.push_back(receiver_signal(bursts_of(g), res.span->first, res.span->second, to_result(res), median_snr(g),
                                    name()));
    }
    return out;
  }

 private:
  static ScaledScan& scan(DecodeContext& ctx) {
    using Ptr = std::shared_ptr<ScaledScan>;
    if (Ptr* p = ctx.cache.get_mut<Ptr>("wifi_ofdm_scaled")) return **p;
    auto res = std::make_shared<ScaledScan>();
    for (auto& ln : scaled::scan(ctx.cap.iq(), ctx.cap.fs, ctx.cap.fc)) {
      auto recs = ofdm_report::records(ln.packets, ctx.cap.fs, ln.clock_hz, true);
      ln.packets.clear();  // the records hold what is needed; free the IQ state
      res->links.push_back({std::move(ln), std::move(recs)});
    }
    ctx.cache.set<Ptr>("wifi_ofdm_scaled", res);
    return *res;
  }
};

// ------------------------------------------------------------------ 802.11b DSSS
// python plugin_dsss.WifiDsssDecoder
class WifiDsssDecoder : public Decoder {
 public:
  std::string name() const override { return "wifi_dsss"; }
  std::string package() const override { return "wifi"; }
  std::string protocol() const override { return "IEEE 802.11b DSSS/CCK"; }
  bool capture_level() const override { return true; }

  bool can_decode(const Signal& sig, DecodeContext&) const override {
    return 2401e6 <= sig.fc() && sig.fc() <= 2495e6 && sig.bw() >= 3e6;
  }

  std::vector<Signal> capture_signals(DecodeContext& ctx) const override {
    if (!(ctx.cap.fc - ctx.cap.fs / 2 < 2495e6 && ctx.cap.fc + ctx.cap.fs / 2 > 2401e6)) return {};
    const auto& sc = scan(ctx);
    std::vector<const Json*> pk;
    for (const auto& p : sc.ppdus)
      if (p.contains("plcp_crc_ok") && p["plcp_crc_ok"].get<bool>()) pk.push_back(&p);
    std::vector<std::pair<int, std::vector<const Json*>>> groups;  // by channel, first-seen order
    for (const Json* p : pk) {
      const int ch = (*p)["channel"].get<int>();
      auto it = std::find_if(groups.begin(), groups.end(), [&](const auto& g) { return g.first == ch; });
      if (it == groups.end()) groups.push_back({ch, {p}});
      else it->second.push_back(p);
    }
    std::vector<Signal> out;
    for (const auto& [ch, g] : groups) {
      const dsss::Summary res = dsss::summarise(g, pk.size(), block_s(ctx));
      std::vector<std::pair<double, double>> bursts;
      for (const Json* p : g) bursts.push_back({(*p)["start_s"].get<double>(), (*p)["end_s"].get<double>()});
      out.push_back(receiver_signal(bursts, res.span->first, res.span->second, to_result(res), std::nullopt, name()));
    }
    return out;
  }

  DecodeResult decode(const Signal& sig, DecodeContext& ctx) const override {
    const auto& sc = scan(ctx);
    std::vector<const Json*> mine;
    for (const auto& p : sc.ppdus) {
      if (!(p.contains("plcp_crc_ok") && p["plcp_crc_ok"].get<bool>())) continue;
      const double f = p["ch_freq_mhz"].get<double>() * 1e6;
      if (overlaps(sig, p["start_s"].get<double>(), p["end_s"].get<double>(), f - 11e6, f + 11e6)) mine.push_back(&p);
    }
    if (!mine.empty()) return to_result(dsss::summarise(mine, sc.ppdus.size(), block_s(ctx)));
    std::vector<const Json*> regs;
    for (const auto& r : sc.unexplained) {
      const double f = r["ch_freq_mhz"].get<double>() * 1e6;
      if (overlaps(sig, r["t0"].get<double>(), r["t1"].get<double>(), f - 11e6, f + 11e6)) regs.push_back(&r);
    }
    if (!regs.empty()) return to_result(dsss::barker_only(regs));
    return to_result(dsss::no_dsss());
  }

 private:
  // python dsss.scan(ctx) with its ctx.cache entries
  static const dsss::ScanResult& scan(DecodeContext& ctx) {
    using Ptr = std::shared_ptr<const dsss::ScanResult>;
    if (const Ptr* p = ctx.cache.get<Ptr>("wifi_dsss")) return **p;
    auto res = std::make_shared<const dsss::ScanResult>(dsss::scan(ctx.ch, ctx.cap.fc, ctx.cap.fs));
    ctx.cache.set<Ptr>("wifi_dsss", res);
    ctx.cache.set<size_t>("wifi_dsss_candidates", res->n_candidates);
    return *res;
  }
};

}  // namespace

void register_plugins(PluginList& out) {
  out.push_back(std::make_unique<WifiOfdmDecoder>());
  out.push_back(std::make_unique<WifiOfdmScaledDecoder>());
  out.push_back(std::make_unique<WifiDsssDecoder>());
}

}  // namespace si::rx::wifi
