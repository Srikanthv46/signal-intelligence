#include "report/writer.hpp"

#include <algorithm>
#include <filesystem>
#include <set>

#include "api/contract.hpp"
#include "core/bands_india.hpp"

namespace si {
namespace {

const double kWeakSnrDb = 6.0;
const char* const kWeakNote =
    "energy-only burst detections below 6 dB in-band SNR with no decoder or rule match: at this level the CFAR "
    "false-alarm rate is set by noise-floor estimate error (~1 dB), so they are listed for completeness but not "
    "reported as signals";

const std::vector<std::string> kPipelineSteps = {
    "detect: coarse+fine STFT, noise floor (sweep template / opening), per-frame gain g(t), continuous + burst CFAR, "
    "merge, valley split, group bursts by frequency span, LO-leak / IQ-image flags",
    "features: timing, envelope, cyclic lines, OFDM CP autocorrelation, FSK",
    "receivers (plugins, see processing.plugins): capture-wide receiver runs; receiver-defined signals absorb "
    "matching energy detections; plugin arbitration, post-processing (e.g. emitter association)",
    "per-signal decoders on remaining signals, then India band rules"};

namespace {
// proven / research tag (python writer._validation): technology from the label and the defining receiver
Json report_validation(const std::string& label, const Json& feats, const Json& fields) {
  const std::string dec = feats.contains("decoder") && feats["decoder"].is_string() ? feats["decoder"].get<std::string>() : "";
  return contract::validation(contract::technology(label, dec).first, dec, fields);
}
}  // namespace

Json signal_dict(size_t i, const Signal& s) {
  const Json& f = s.features;
  std::vector<Box> bursts = s.bursts;
  std::stable_sort(bursts.begin(), bursts.end(), [](const Box& a, const Box& b) { return a.t0 < b.t0; });
  Json feats = {{"symbol_rate", f.value("symbol_rate", Json(nullptr))},
                {"burst_period_s", f.value("burst_period_s", Json(nullptr))},
                {"duty_cycle", f.value("duty_cycle", Json(nullptr))}};
  for (auto it = f.begin(); it != f.end(); ++it) {
    const std::string& k = it.key();
    if (feats.contains(k) || k == "classifier_candidates" || (!k.empty() && k[0] == '_')) continue;
    feats[k] = it.value();
  }
  Json bl = Json::array();
  for (size_t j = 0; j < bursts.size() && j < 300; ++j) {
    const Box& b = bursts[j];
    bl.push_back({{"t_start_s", py_round(b.t0, 7)},
                  {"t_end_s", py_round(b.t1, 7)},
                  {"f_low_hz", py_round_int(b.f0)},
                  {"f_high_hz", py_round_int(b.f1)},
                  {"snr_db", b.mean_snr_db}});
  }
  feats["bursts"] = bl;
  Json det = s.detection;
  if (f.contains("classifier_candidates") && !f["classifier_candidates"].is_null())
    det["classifier_candidates"] = f["classifier_candidates"];
  const std::string dec = feats.contains("decoder") && feats["decoder"].is_string() ? feats["decoder"].get<std::string>() : "";
  return {{"id", i},
          {"label", s.label},
          {"description", contract::description(s.label, dec, Json(s.flags))},
          {"confidence", s.confidence},
          {"t_start_s", py_round(s.t0(), 7)},
          {"t_end_s", py_round(s.t1(), 7)},
          {"f_low_hz", py_round_int(s.f0)},
          {"f_high_hz", py_round_int(s.f1)},
          {"f_center_hz", py_round_int(s.fc())},
          {"bandwidth_hz", py_round_int(s.bw())},
          {"snr_db", s.snr_db},
          {"modulation", s.modulation},
          {"features", feats},
          {"evidence", s.evidence},
          {"band_allocation_india", bands::describe(s.fc())},
          {"flags", s.flags},
          {"detection", det},
          {"decoded",
           {{"status", s.decoded.status},
            {"reason", s.decoded.reason},
            {"method", s.decoded.method ? Json(*s.decoded.method) : Json(nullptr)},
            {"steps", s.decoded.steps},
            {"fields", s.decoded.fields},
            {"validation", report_validation(s.label, feats, s.decoded.fields)}}}};
}

bool is_weak(const Json& s) {
  if (s["decoded"]["status"] != "not_possible") return false;
  const Json kinds = s["detection"].value("kinds", Json::array());
  if (std::find(kinds.begin(), kinds.end(), Json("continuous")) != kinds.end()) return false;
  for (const auto& fl : s["flags"])
    if (fl == "lo_leakage" || fl == "receiver_spur") return false;
  const double snr = s["snr_db"].is_number() ? s["snr_db"].get<double>() : 0.0;
  return snr < kWeakSnrDb;
}

void split_weak(Json& rep) {
  Json keep = Json::array();
  Json weak = rep.value("weak_detections", Json::array());
  for (auto& s : rep["signals"]) {
    if (is_weak(s)) {
      weak.push_back({{"t_start_s", s["t_start_s"]},
                      {"t_end_s", s["t_end_s"]},
                      {"f_low_hz", s["f_low_hz"]},
                      {"f_high_hz", s["f_high_hz"]},
                      {"snr_db", s["snr_db"]},
                      {"rule_label", s["label"]},
                      {"flags", s["flags"]}});
    } else {
      keep.push_back(std::move(s));
    }
  }
  for (size_t i = 0; i < keep.size(); ++i) keep[i]["id"] = i;
  rep["signals"] = keep;
  rep["weak_detections"] = weak;
  rep["weak_detections_note"] = kWeakNote;
}

}  // namespace

Json build_report(const Analysis& a) {
  const Capture& cap = a.cap;
  Json quality = a.quality;
  quality["noise_floor"] = a.det.info.contains("noise") ? a.det.info["noise"] : Json(nullptr);
  quality["usable_band_hz"] = Json::array({py_round_int(a.det.usable_hz.first), py_round_int(a.det.usable_hz.second)});
  Json label_check = a.label_check ? *a.label_check
                                   : Json{{"status", "not_checked"},
                                          {"note", "no plugin evidence (e.g. CRC-valid BLE adverts)"}};
  Json source = nullptr;
  if (!cap.meta_path.empty() && std::filesystem::exists(cap.meta_path))
    source = std::filesystem::absolute(cap.meta_path).lexically_normal().string();
  Json signals = Json::array();
  for (size_t i = 0; i < a.signals.size(); ++i) signals.push_back(signal_dict(i, a.signals[i]));
  Json plugins = Json::array();
  for (const auto* p : a.plugins) plugins.push_back(p->name());
  Json rep = {{"file", cap.name},
              {"source_path", source},
              {"center_freq_hz", cap.fc},
              {"sample_rate_hz", cap.fs},
              {"capture_time_utc", cap.t0 ? cap.t0->isoformat() : std::string()},
              {"frequency_label_check", label_check},
              {"duration_s", py_round(cap.duration(), 6)},
              {"gain_db", cap.gain_db ? Json(*cap.gain_db) : Json(nullptr)},
              {"capture_quality", quality},
              {"processing", {{"timings_s", a.timings}, {"pipeline", kPipelineSteps}}},
              {"signals", signals}};
  for (auto it = a.extras.begin(); it != a.extras.end(); ++it) rep[it.key()] = it.value();
  rep["processing"]["plugins"] = plugins;
  split_weak(rep);
  for (const auto* p : a.plugins) p->finalize_report(rep);
  return rep;
}

}  // namespace si
