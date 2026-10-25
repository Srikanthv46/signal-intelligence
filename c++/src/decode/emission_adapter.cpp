#include "decode/emission_adapter.hpp"

#include <algorithm>
#include <map>

#include "core/bands_india.hpp"
#include "core/npcompat.hpp"
#include "core/pyfmt.hpp"

namespace si {
namespace {

bool intersects(const std::set<std::string>& a, const std::vector<std::string>& b) {
  for (const auto& x : b)
    if (a.count(x)) return true;
  return false;
}

}  // namespace

bool in_ranges(double lo, double hi, const std::vector<std::pair<double, double>>& ranges) {
  for (const auto& [a, b] : ranges)
    if (lo < b && hi > a) return true;
  return false;
}

Json scalar_features(const Json& features) {
  Json out = Json::object();
  for (auto it = features.begin(); it != features.end(); ++it)
    if (!(it.value().is_array() || it.value().is_object()) || it.key() == "sync_syms") out[it.key()] = it.value();
  return out;
}

bool EmissionReceiver::in_range(double lo, double hi) const { return intersects(bands(), bands::band_keys(lo, hi)); }

bool EmissionReceiver::can_decode(const Signal& sig, DecodeContext& ctx) const {
  const double lo = ctx.cap.fc - ctx.cap.fs / 2, hi = ctx.cap.fc + ctx.cap.fs / 2;
  return intersects(bands(), bands::band_keys(lo, hi)) && intersects(bands(), bands::band_keys(sig.f0, sig.f1));
}

const std::vector<Emission>& EmissionReceiver::emissions(DecodeContext& ctx) const {
  static const std::vector<Emission> none;
  const double lo = ctx.cap.fc - ctx.cap.fs / 2, hi = ctx.cap.fc + ctx.cap.fs / 2;
  if (!in_range(lo, hi)) return none;
  if (!ctx.cache.contains(key())) ctx.cache.set(key(), run_detector(ctx));
  return *ctx.cache.get<std::vector<Emission>>(key());
}

bool EmissionReceiver::keep(const Emission& e) const {
  const bool crc = e.decoded.is_object() && e.decoded.contains("crc_ok") && e.decoded["crc_ok"].is_boolean() &&
                   e.decoded["crc_ok"].get<bool>();
  return e.confidence >= 0.5 || crc;
}

std::string EmissionReceiver::group_key(const Emission& e) const {
  // python: (e.label.split("_sf")[0] if "lora" in e.label else e.label, round(e.center_mhz, 1))
  std::string lab = e.label;
  if (lab.find("lora") != std::string::npos) lab = lab.substr(0, lab.find("_sf"));
  return lab + "|" + py_float_repr(py_round(e.center_mhz, 1));
}

std::vector<Signal> EmissionReceiver::capture_signals(DecodeContext& ctx) const {
  std::vector<std::string> order;
  std::map<std::string, std::vector<Emission>> groups;
  for (const auto& e : emissions(ctx)) {
    if (!keep(e)) continue;  // weak energy-only candidates stay with the energy detector
    const std::string k = group_key(e);
    if (!groups.count(k)) order.push_back(k);
    groups[k].push_back(e);
  }
  std::vector<Signal> out;
  for (const auto& k : order) {
    const auto& g = groups[k];
    DecodeResult res = summarise(g, nullptr);
    const auto [f0, f1] = *res.span;
    dvec snr;
    for (const auto& e : g)
      if (e.snr_db) snr.push_back(*e.snr_db);
    std::vector<std::pair<double, double>> bursts;
    for (const auto& e : g) bursts.push_back({e.start_us * 1e-6, e.end_us * 1e-6});
    Signal s = receiver_signal(bursts, f0, f1, res, snr.empty() ? std::nullopt : std::optional<double>(np::median(snr)),
                               name());
    s.emissions = g;  // for emitter association (not serialised)
    if (owns_band()) s.features["owns_band"] = true;
    finish(s, g, ctx);
    out.push_back(std::move(s));
  }
  return out;
}

DecodeResult EmissionReceiver::decode(const Signal& sig, DecodeContext& ctx) const {
  std::vector<Emission> em;
  for (const auto& e : emissions(ctx))
    if (overlaps(sig, e.start_us * 1e-6, e.end_us * 1e-6, e.f_lo_mhz * 1e6, e.f_hi_mhz * 1e6)) em.push_back(e);
  if (!em.empty()) return summarise(em, &sig);
  return result("not_possible", method(), steps(), Json::object(),
                protocol() + " receiver found no matching burst in this signal's span");
}

DecodeResult EmissionReceiver::summarise(const std::vector<Emission>& em, const Signal*) const {
  // Counter(e.label): insertion order, most_common ties -> first seen
  Json labels = Json::object();
  for (const auto& e : em) labels[e.label] = labels.value(e.label, 0) + 1;
  std::string top;
  int best = -1;
  for (auto it = labels.begin(); it != labels.end(); ++it)
    if (it.value().get<int>() > best) {
      best = it.value().get<int>();
      top = it.key();
    }
  size_t ok = 0;
  for (const auto& e : em) {
    if (!e.decoded.is_object() || e.decoded.empty()) continue;  // python: `if e.decoded` (non-empty dict)
    auto truthy = [&](const char* k) {
      if (!e.decoded.contains(k)) return false;
      const Json& v = e.decoded[k];
      if (v.is_boolean()) return v.get<bool>();
      if (v.is_number()) return v.get<double>() != 0;
      if (v.is_string()) return !v.get<std::string>().empty();
      return !v.is_null() && !v.empty();
    };
    if (truthy("crc_ok") || truthy("crc") || truthy("serial")) ++ok;
  }
  Json bursts = Json::array();
  for (size_t i = 0; i < em.size() && i < 100; ++i) {
    const auto& e = em[i];
    bursts.push_back({{"t_start_s", e.start_us * 1e-6},
                      {"t_end_s", e.end_us * 1e-6},
                      {"f_center_hz", e.center_mhz * 1e6},
                      {"label", e.label},
                      {"snr_db", e.snr_db ? Json(*e.snr_db) : Json(nullptr)},
                      {"features", scalar_features(e.features)},
                      {"decoded", e.decoded}});
  }
  Json fields = {{"n_bursts", em.size()}, {"labels", labels}, {"bursts", bursts}};
  double conf = 0;
  double lo = em[0].f_lo_mhz, hi = em[0].f_hi_mhz;
  for (const auto& e : em) {
    conf = std::max(conf, e.confidence);
    lo = std::min(lo, e.f_lo_mhz);
    hi = std::max(hi, e.f_hi_mhz);
  }
  std::string ev = std::to_string(em.size()) + " burst(s) matched by " + protocol() + " receiver: " + py_repr(labels);
  if (ok) ev += "; " + std::to_string(ok) + " CRC-verified payload(s)";
  return result(ok ? "success" : "partial", method(), steps(), fields,
                ok ? "" : "physical layer identified/measured; no CRC-verified payload", pretty(top, em),
                std::min(0.97, conf + (ok ? 0.1 : 0.0)), modulation(em), ev,
                std::make_pair(lo * 1e6, hi * 1e6));
}

std::string EmissionReceiver::modulation(const std::vector<Emission>& em) const {
  const Json& f = em[0].features;
  return f.contains("modulation") && f["modulation"].is_string() ? f["modulation"].get<std::string>() : "";
}

}  // namespace si
