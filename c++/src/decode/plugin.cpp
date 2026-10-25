#include "decode/plugin.hpp"

#include <algorithm>

namespace si {

DecodeResult result(std::string status, std::string method, std::vector<std::string> steps, Json fields,
                    std::string reason, std::optional<std::string> label, std::optional<double> confidence,
                    std::optional<std::string> modulation, std::optional<std::string> evidence,
                    std::optional<std::pair<double, double>> span) {
  DecodeResult r;
  r.status = std::move(status);
  r.method = std::move(method);
  r.steps = std::move(steps);
  r.fields = fields.is_null() ? Json::object() : std::move(fields);
  r.reason = std::move(reason);
  r.label = std::move(label);
  r.confidence = confidence;
  r.modulation = std::move(modulation);
  r.evidence = std::move(evidence);
  r.span = span;
  return r;
}

Signal receiver_signal(const std::vector<std::pair<double, double>>& bursts, double f0, double f1,
                       const DecodeResult& res, std::optional<double> snr_db, const std::string& decoder) {
  auto sorted = bursts;
  std::sort(sorted.begin(), sorted.end());
  Signal s;
  const double snr = snr_db.value_or(0.0);
  for (const auto& [t0, t1] : sorted) {
    Box b;
    b.t0 = t0;
    b.t1 = t1;
    b.f0 = f0;
    b.f1 = f1;
    b.peak_snr_db = snr;
    b.mean_snr_db = snr;
    b.res = "receiver";
    s.bursts.push_back(b);
  }
  s.f0 = f0;
  s.f1 = f1;
  s.snr_db = snr_db ? py_round(*snr_db, 1) : 0.0;
  s.decoded.status = res.status;
  s.decoded.reason = res.reason;
  s.decoded.method = res.method;
  s.decoded.steps = res.steps;
  s.decoded.fields = res.fields;
  s.label = res.label.value_or("unknown");
  s.confidence = res.confidence.value_or(0.0);
  s.modulation = res.modulation.value_or("");
  s.evidence = res.evidence.value_or("");
  const size_t n = s.bursts.size();
  s.detection = {{"method", "protocol receiver (" + decoder + ")"},
                 {"description", "defined by " + std::to_string(n) + " burst(s) the " + decoder +
                                     " receiver synchronised to and decoded (preamble/sync correlation + header "
                                     "checks), not by energy alone"},
                 {"n_bursts", n}};
  s.features["decoder"] = decoder;
  return s;
}

bool overlaps(const Signal& sig, double t0, double t1, double f0, double f1, double f_frac) {
  const double fov = std::max(0.0, std::min(sig.f1, f1) - std::max(sig.f0, f0));
  if (fov < f_frac * std::min(sig.bw(), f1 - f0)) return false;
  for (const auto& b : sig.bursts)
    if (std::min(b.t1, t1) > std::max(b.t0, t0) - 50e-6) return true;
  return false;
}

}  // namespace si
