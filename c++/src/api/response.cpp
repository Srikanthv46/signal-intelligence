#include "api/response.hpp"

#include <algorithm>
#include <cstdio>
#include <map>
#include <set>

#include "api/contract.hpp"

namespace si::api {
namespace {
bool truthy_str(const Json& v) { return v.is_string() && !v.get<std::string>().empty(); }

const std::set<std::string> kTimingFeatures = {"bursts", "duty_cycle", "burst_period_s", "n_bursts", "continuous",
                                               "decoder"};

Json utc(const std::optional<UtcTime>& t0, double s) {
  if (!t0) return nullptr;
  return t0->plus(s).iso_z();
}

double num(const Json& j, double dflt = 0.0) { return j.is_number() ? j.get<double>() : dflt; }

Json frequency(const Json& lo, const Json& hi, double fc) {
  const double l = num(lo), h = num(hi);
  const int64_t c = py_round_int(0.5 * (l + h));
  return {{"center_hz", c},
          {"low_hz", lo},
          {"high_hz", hi},
          {"bandwidth_hz", py_round_int(h - l)},
          {"offset_hz", py_round_int(static_cast<double>(c) - fc)}};
}

Json time_block(const Json& t0s, const Json& t1s, const std::optional<UtcTime>& start) {
  const double a = num(t0s), b = num(t1s);
  return {{"start_s", t0s},
          {"end_s", t1s},
          {"duration_s", py_round(b - a, 7)},
          {"start_utc", utc(start, a)},
          {"end_utc", utc(start, b)}};
}

bool empty_value(const Json& v) {
  return v.is_null() || (v.is_array() && v.empty()) || (v.is_object() && v.empty()) ||
         (v.is_string() && v.get<std::string>().empty());
}

Json signal(const Json& s, double fc, const std::optional<UtcTime>& start, bool compact) {
  const Json& f = s.at("features");
  const std::string label = s.at("label").get<std::string>();
  const std::string decoder = f.contains("decoder") && f["decoder"].is_string() ? f["decoder"].get<std::string>() : "";
  auto [tech, fam] = contract::technology(label, decoder);
  const Json& flags = s.at("flags");
  auto has = [&](const char* x) {
    return std::any_of(flags.begin(), flags.end(), [&](const Json& v) { return v == x; });
  };
  bool artifact = fam == "artifact";
  for (const char* x : contract::kArtifactFlags) artifact = artifact || has(x);
  if (artifact && fam != "artifact") {
    for (const char* x : contract::kArtifactFlags)
      if (has(x)) {
        tech = x;
        break;
      }
    fam = "artifact";
  }
  const Json& dec = s.at("decoded");
  const Json fields = dec.contains("fields") && dec["fields"].is_object() ? dec["fields"] : Json::object();
  Json t = time_block(s.at("t_start_s"), s.at("t_end_s"), start);
  t["n_bursts"] = f.contains("n_bursts") ? f["n_bursts"] : Json(f.value("bursts", Json::array()).size());
  t["duty_cycle"] = f.value("duty_cycle", Json(nullptr));
  t["burst_period_s"] = f.value("burst_period_s", Json(nullptr));
  t["continuous"] = f.contains("continuous") && f["continuous"].is_boolean() && f["continuous"].get<bool>();
  Json ids = Json::object();
  for (const auto& k : contract::identifier_keys())
    if (fields.contains(k) && !empty_value(fields[k])) ids[k] = fields[k];
  const std::string status = dec.at("status").get<std::string>();
  Json out = {
      {"id", s.at("id")},
      {"category", artifact ? "artifact" : "emission"},
      {"technology", tech},
      {"family", fam},
      {"label", label},
      {"description", s.contains("description") && truthy_str(s["description"]) ? s["description"] : Json(contract::description(label, decoder, flags))},
      {"confidence", s.at("confidence")},
      {"inferred", status == "not_possible" && label != "unknown"},
      {"frequency", frequency(s.at("f_low_hz"), s.at("f_high_hz"), fc)},
      {"time", t},
      {"snr_db", s.at("snr_db")},
      {"modulation", s.at("modulation")},
      {"decode",
       {{"status", status},
        {"reason", dec.value("reason", Json("")).is_string() ? dec.value("reason", Json("")) : Json("")},
        {"method", dec.value("method", Json(nullptr))},
        {"fields", fields}}},
      {"identifiers", ids},
      {"validation", contract::validation(tech, f.contains("decoder") && f["decoder"].is_string() ? f["decoder"].get<std::string>() : std::string(), fields)},
      {"evidence", s.at("evidence")},
      {"band_allocation", s.at("band_allocation_india")},
      {"flags", flags},
  };
  if (!compact) {
    out["decode"]["steps"] = dec.value("steps", Json::array());
    Json bursts = Json::array();
    for (const auto& b : f.value("bursts", Json::array()))
      bursts.push_back({{"start_s", b.at("t_start_s")},
                        {"end_s", b.at("t_end_s")},
                        {"low_hz", b.at("f_low_hz")},
                        {"high_hz", b.at("f_high_hz")},
                        {"snr_db", b.at("snr_db")}});
    out["time"]["bursts"] = bursts;
    out["detection"] = s.at("detection");
    Json feats = Json::object();
    for (auto it = f.begin(); it != f.end(); ++it)
      if (!kTimingFeatures.count(it.key())) feats[it.key()] = it.value();
    out["features"] = feats;
  }
  return out;
}

Json weak(const Json& w, double fc, const std::optional<UtcTime>& start) {
  return {{"frequency", frequency(w.at("f_low_hz"), w.at("f_high_hz"), fc)},
          {"time", time_block(w.at("t_start_s"), w.at("t_end_s"), start)},
          {"snr_db", w.at("snr_db")},
          {"rule_label", w.at("rule_label")},
          {"flags", w.at("flags")}};
}

std::string fmt_g(double v) {
  char b[64];
  std::snprintf(b, sizeof b, "%g", v);
  return b;
}

std::string fmt_3g(double v) {
  char b[64];
  std::snprintf(b, sizeof b, "%.3g", v);
  return b;
}

void coverage(const Packages& packages, double fs, std::optional<double> bw, std::pair<double, double> band,
              double dur, const Json& timings, Json& cov, Json& warn) {
  const Json& req = contract::requirements().at("receivers");
  const double rate = bw ? std::min(fs, *bw) : fs;
  for (const auto& [pkg, names] : packages) {
    Json r = req.contains(pkg) ? req[pkg]
                               : Json{{"min_sample_rate_hz", 0}, {"min_duration_s", 0},
                                      {"recommended_duration_s", 0}, {"bands_hz", Json::array()}, {"why", ""}};
    const Json bands = r.value("bands_hz", Json::array());
    bool in_band = bands.empty();
    for (const auto& b : bands)
      if (b[0].get<double>() < band.second && b[1].get<double>() > band.first) in_band = true;
    const double min_rate = r.at("min_sample_rate_hz").get<double>();
    const double min_dur = r.at("min_duration_s").get<double>();
    const double rec_dur = r.at("recommended_duration_s").get<double>();
    const bool rate_ok = rate >= min_rate, dur_ok = dur >= min_dur, full = dur >= rec_dur;
    bool err = false;
    double secs = 0;
    bool any_secs = false;
    for (const auto& n : names) {
      if (timings.contains(n + "_error")) err = true;
      if (timings.contains("receiver_s") && timings["receiver_s"].contains(n) && timings["receiver_s"][n].is_number()) {
        secs += timings["receiver_s"][n].get<double>();
        any_secs = true;
      }
    }
    const std::string why = r.value("why", std::string());
    cov.push_back({{"receiver", pkg},
                   {"status", err ? "error" : "ran"},
                   {"in_band", in_band},
                   {"sample_rate_ok", rate_ok},
                   {"duration_ok", dur_ok},
                   {"full_duration", full},
                   {"min_sample_rate_hz", r.at("min_sample_rate_hz")},
                   {"min_duration_s", r.at("min_duration_s")},
                   {"recommended_duration_s", r.at("recommended_duration_s")},
                   {"time_s", any_secs ? Json(py_round(secs, 3)) : Json(nullptr)},
                   {"note", (rate_ok && full) ? std::string() : why}});
    if (in_band && !rate_ok)
      warn.push_back({{"code", "SAMPLE_RATE_TOO_LOW"},
                      {"receiver", pkg},
                      {"message", pkg + " needs >= " + fmt_g(min_rate / 1e6) + " MHz of bandwidth; this block has " +
                                      fmt_g(rate / 1e6) + " MHz"}});
    if (in_band && !dur_ok)
      warn.push_back({{"code", "BLOCK_TOO_SHORT_FOR_RECEIVER"},
                      {"receiver", pkg},
                      {"message", pkg + " needs >= " + fmt_g(min_dur * 1e3) + " ms; this block is " +
                                      fmt_3g(dur * 1e3) + " ms (" + why + ")"}});
  }
}

Json engine_input(const BlockMeta& b, double fc_used, double fs, size_t n, double dur,
                  std::pair<double, double> band) {
  return {{"source_id", b.source_id ? Json(*b.source_id) : Json(nullptr)},
          {"center_frequency_hz", b.center_frequency_hz},
          {"center_frequency_used_hz", fc_used},
          {"sample_rate_hz", fs},
          {"bandwidth_hz", b.bandwidth_hz ? Json(*b.bandwidth_hz) : Json(nullptr)},
          {"n_samples", n},
          {"duration_s", dur},
          {"start_time_utc", utc(b.start, 0.0)},
          {"gain_db", b.gain_db ? Json(*b.gain_db) : Json(nullptr)},
          {"receiver_ppm", b.receiver_ppm ? Json(*b.receiver_ppm) : Json(nullptr)},
          {"analysed_band_hz", Json::array({band.first, band.second})}};
}

}  // namespace

Json build_response(const Json& rep, const Analysis& a, const BlockMeta& block, const Json& engine,
                    const Packages& packages, double total_s, bool compact) {
  const double fc_used = a.cap.fc;
  const Json& q = rep.at("capture_quality");
  const std::pair<double, double> band{q.at("usable_band_hz")[0].get<double>(), q.at("usable_band_hz")[1].get<double>()};
  const double dur = rep.at("duration_s").get<double>();
  const Json& T = rep.at("processing").at("timings_s");

  Json signals = Json::array();
  for (const auto& s : rep.at("signals")) signals.push_back(signal(s, fc_used, block.start, compact));
  Json weaks = Json::array();
  for (const auto& w : rep.value("weak_detections", Json::array())) weaks.push_back(weak(w, fc_used, block.start));

  Json cov = Json::array(), warnings = Json::array();
  coverage(packages, a.cap.fs, block.bandwidth_hz, band, dur, T, cov, warnings);
  const Json& preq = contract::requirements().at("pipeline");
  const double rec_n = preq.at("recommended_samples").get<double>();
  const double rec_d = preq.at("recommended_duration_s").get<double>();
  if (static_cast<double>(a.cap.n()) < rec_n && dur < rec_d) {
    warnings.insert(warnings.begin(),
                    Json{{"code", "SHORT_BLOCK"},
                         {"receiver", nullptr},
                         {"message", std::to_string(a.cap.n()) + " samples (" + fmt_3g(dur * 1e3) +
                                         " ms): below the recommended " + std::to_string(static_cast<int64_t>(rec_n)) +
                                         " samples / " + fmt_g(rec_d * 1e3) +
                                         " ms; expect fewer decodes and more spurious energy detections"}});
  }
  if (q.value("clipping", false))
    warnings.push_back({{"code", "ADC_CLIPPING"},
                        {"receiver", nullptr},
                        {"message", fmt_3g(num(q.value("clipped_fraction", Json(0.0))) * 100) +
                                        " % of samples at ADC full scale: intermodulation products may appear as "
                                        "signals; lower the gain"}});
  const Json lc = rep.value("frequency_label_check", Json(nullptr));
  if (lc.is_object() && lc.value("status", std::string()).rfind("corrected", 0) == 0) {
    char b[64];
    std::snprintf(b, sizeof b, "%.6f", fc_used / 1e6);
    warnings.push_back({{"code", "CENTER_FREQUENCY_CORRECTED"},
                        {"receiver", lc.value("plugin", Json(nullptr))},
                        {"message", "decoded content proves the block's centre frequency wrong by " +
                                        lc.value("offset_mhz", Json(nullptr)).dump() + " MHz; analysed at " + b +
                                        " MHz"}});
  }
  Json errors = Json::array();
  for (auto it = T.begin(); it != T.end(); ++it) {
    const std::string& k = it.key();
    if (k.size() > 6 && k.compare(k.size() - 6, 6, "_error") == 0 && it.value().is_string())
      errors.push_back({{"code", "RECEIVER_ERROR"}, {"receiver", k.substr(0, k.size() - 6)}, {"message", it.value()}});
  }

  std::map<std::string, int> by_status_m, by_family_m;
  std::vector<std::string> status_order, family_order;
  std::set<std::string> techs;
  int n_art = 0, n_unk = 0;
  for (const auto& s : signals) {
    const std::string st = s["decode"]["status"], fam = s["family"];
    if (!by_status_m.count(st)) status_order.push_back(st);
    if (!by_family_m.count(fam)) family_order.push_back(fam);
    ++by_status_m[st];
    ++by_family_m[fam];
    techs.insert(s["technology"].get<std::string>());
    n_art += s["category"] == "artifact";
    n_unk += s["technology"] == "unknown";
  }
  Json by_status = Json::object(), by_family = Json::object();
  for (const auto& k : status_order) by_status[k] = by_status_m[k];
  for (const auto& k : family_order) by_family[k] = by_family_m[k];

  static const std::set<std::string> consumed = {"file", "source_path", "center_freq_hz", "sample_rate_hz",
                                                 "capture_time_utc", "frequency_label_check", "duration_s", "gain_db",
                                                 "capture_quality", "processing", "signals", "weak_detections",
                                                 "weak_detections_note", "emitters", "devices", "findings"};
  Json ext = Json::object();
  for (auto it = rep.begin(); it != rep.end(); ++it)
    if (!consumed.count(it.key())) ext[it.key()] = it.value();
  Json stages = Json::object();
  for (auto it = T.begin(); it != T.end(); ++it) {
    const std::string& k = it.key();
    if (k.size() > 2 && k.compare(k.size() - 2, 2, "_s") == 0 && it.value().is_number())
      stages[k.substr(0, k.size() - 2)] = py_round(it.value().get<double>() * 1e3, 1);
  }
  Json rx = Json::object();
  if (T.contains("receiver_s"))
    for (auto it = T["receiver_s"].begin(); it != T["receiver_s"].end(); ++it)
      rx[it.key()] = py_round(it.value().get<double>() * 1e3, 1);

  Json noise_med = a.det.info.contains("noise_floor_dbfs_per_bin_median") ? a.det.info["noise_floor_dbfs_per_bin_median"]
                                                                          : Json(nullptr);
  return {
      {"schema_version", contract::kSchemaVersion},
      {"kind", "spectrum_intel.analysis"},
      {"request_id", block.request_id ? Json(*block.request_id) : Json(nullptr)},
      {"status", "ok"},
      {"engine", engine},
      {"input", engine_input(block, fc_used, a.cap.fs, a.cap.n(), dur, band)},
      {"quality",
       {{"adc_clipping", q.value("clipping", false)},
        {"clipped_fraction", q.value("clipped_fraction", Json(nullptr))},
        {"peak_abs", q.value("peak_abs", Json(nullptr))},
        {"full_scale", q.value("full_scale", Json(nullptr))},
        {"rms_dbfs", q.value("rms_dbfs", Json(nullptr))},
        {"noise_floor_dbfs_per_bin_median", noise_med},
        {"frequency_label_check", (lc.is_object() && !lc.empty()) ? lc : Json(nullptr)}}},
      {"coverage", cov},
      {"summary",
       {{"n_signals", signals.size()},
        {"n_weak_detections", weaks.size()},
        {"n_artifacts", n_art},
        {"n_unknown", n_unk},
        {"by_status", by_status},
        {"by_family", by_family},
        {"technologies", Json(std::vector<std::string>(techs.begin(), techs.end()))}}},
      {"signals", signals},
      {"weak_detections", weaks},
      {"emitters", rep.value("emitters", Json::array())},
      {"devices", rep.value("devices", Json::array())},
      {"findings", rep.value("findings", Json::array())},
      {"extensions", ext},
      {"warnings", warnings},
      {"errors", errors},
      {"timing_ms", {{"total", py_round(total_s * 1e3, 1)}, {"stages", stages}, {"receivers", rx}}},
  };
}

Json failure_response(const BlockMeta& b, const Json& engine, const std::string& status, const std::string& code,
                      const std::string& message, double total_s) {
  const double fs = b.sample_rate_hz > 0 ? b.sample_rate_hz : 1.0;
  const double fc = b.center_frequency_hz;
  return {
      {"schema_version", contract::kSchemaVersion},
      {"kind", "spectrum_intel.analysis"},
      {"request_id", b.request_id ? Json(*b.request_id) : Json(nullptr)},
      {"status", status},
      {"engine", engine},
      {"input",
       {{"source_id", b.source_id ? Json(*b.source_id) : Json(nullptr)},
        {"center_frequency_hz", fc},
        {"sample_rate_hz", fs},
        {"bandwidth_hz", b.bandwidth_hz ? Json(*b.bandwidth_hz) : Json(nullptr)},
        {"n_samples", b.n_samples},
        {"duration_s", static_cast<double>(b.n_samples) / fs},
        {"start_time_utc", utc(b.start, 0.0)},
        {"gain_db", b.gain_db ? Json(*b.gain_db) : Json(nullptr)},
        {"receiver_ppm", b.receiver_ppm ? Json(*b.receiver_ppm) : Json(nullptr)},
        {"analysed_band_hz", Json::array({fc - fs / 2, fc + fs / 2})}}},
      {"quality", Json::object()},
      {"coverage", Json::array()},
      {"summary",
       {{"n_signals", 0},
        {"n_weak_detections", 0},
        {"n_artifacts", 0},
        {"n_unknown", 0},
        {"by_status", Json::object()},
        {"by_family", Json::object()},
        {"technologies", Json::array()}}},
      {"signals", Json::array()},
      {"weak_detections", Json::array()},
      {"emitters", Json::array()},
      {"devices", Json::array()},
      {"findings", Json::array()},
      {"extensions", Json::object()},
      {"warnings", Json::array()},
      {"errors", Json::array({Json{{"code", code}, {"receiver", nullptr}, {"message", message}}})},
      {"timing_ms", {{"total", py_round(total_s * 1e3, 1)}}},
  };
}

}  // namespace si::api
