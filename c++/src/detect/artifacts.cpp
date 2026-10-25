#include "detect/artifacts.hpp"

#include <algorithm>
#include <cmath>
#include <map>

#include "core/npcompat.hpp"
#include "core/pyfmt.hpp"

namespace si::detect {

const double kAdcFullScale = 2048.0;

Json clipping(const cfvec& x, double full_scale) {
  const size_t n = x.size();
  fvec a(n), p2(n);
  float peak = 0.0f;
  for (size_t i = 0; i < n; ++i) {
    a[i] = std::max(std::fabs(x[i].real()), std::fabs(x[i].imag()));
    peak = std::max(peak, a[i]);
    const float m = np::cabs(x[i]);
    p2[i] = m * m;
  }
  // raw counts if the peak is in the ADC range; else assume normalised [-1, 1]
  const double fs_ = peak > 4 ? full_scale : 1.0;
  // raw counts reach full_scale - 1 (2047); normalised samples reach ~1.0 (python float vs float32 array)
  const float thr = static_cast<float>(0.995 * (fs_ > 1 ? fs_ - 1 : fs_));
  size_t clipped = 0;
  for (size_t i = 0; i < n; ++i) clipped += a[i] >= thr;
  const double frac = n ? static_cast<double>(clipped) / static_cast<double>(n) : 0.0;
  const float ms = n ? np::sum(p2.data(), n) / static_cast<float>(n) : 0.0f;  // float32 mean
  const double rms = static_cast<double>(std::sqrt(ms));
  return {{"peak_abs", py_round(static_cast<double>(peak), 1)},
          {"full_scale", fs_},
          {"clipped_fraction", frac},
          {"clipping", frac > 1e-6},
          {"rms_dbfs", py_round(20.0 * std::log10(rms / fs_ + 1e-12), 1)}};
}

void flag_receiver_spurs(std::vector<Signal>& sigs, const Json& lines, double fs, std::optional<double> ppm,
                         double max_bw) {
  for (auto& s : sigs) {
    if (!s.flags.empty()) continue;
    const Json* best = nullptr;
    Json on = Json::array();
    for (const auto& ln : lines) {
      const double f = ln["freq_hz"].get<double>();
      if (s.f0 - 1e3 <= f && f <= s.f1 + 1e3) {
        on.push_back(f);
        if (!best || ln["peak_db"].get<double>() > (*best)["peak_db"].get<double>()) best = &ln;
      }
    }
    if (!best) continue;
    if (s.bw() <= max_bw) {
      s.flags.push_back("receiver_spur");
      s.features["spur_kind"] = (*best)["kind"];
      s.features["spur_freq_hz"] = (*best)["freq_hz"];
      s.features["spur_sample_rate_hz"] = fs;
      if (!(*best)["harmonic"].is_null()) s.features["spur_harmonic"] = (*best)["harmonic"];
      if ((*best)["kind"] == "reference_locked") s.features["spur_receiver_ppm"] = py_round(*ppm, 2);
    } else {
      s.features["contains_receiver_spur_hz"] = on;
    }
  }
}

void flag_lo_leakage(std::vector<Signal>& sigs, double fc, double df_fine) {
  for (auto& s : sigs)
    if (s.bw() <= 6 * df_fine + 1 && std::fabs(s.fc() - fc) <= 2 * df_fine) s.flags.push_back("lo_leakage");
}

void flag_iq_images(std::vector<Signal>& sigs, double fc, double min_rej_db) {
  // An image is the mirror (2fc - f) of a stronger signal with the same burst timing, at least
  // min_rej_db weaker (AD936x image rejection is ~40-50 dB: only very strong emitters leave images).
  for (size_t i = 0; i < sigs.size(); ++i) {
    Signal& s = sigs[i];
    if (s.has_flag("lo_leakage")) continue;
    const double m0 = 2 * fc - s.f1, m1 = 2 * fc - s.f0;
    double on = 0.0;
    for (const auto& b : s.bursts) on += b.t1 - b.t0;
    for (size_t j = 0; j < sigs.size(); ++j) {
      const Signal& o = sigs[j];
      if (j == i || o.snr_db - s.snr_db < min_rej_db) continue;
      const double fov = std::max(0.0, std::min(m1, o.f1) - std::max(m0, o.f0));
      if (fov < 0.6 * std::max(s.bw(), 1.0)) continue;
      double tov = 0.0;
      for (const auto& a : s.bursts)
        for (const auto& b : o.bursts) tov += std::max(0.0, std::min(a.t1, b.t1) - std::max(a.t0, b.t0));
      if (tov >= 0.7 * on) {
        s.flags.push_back("iq_image");
        s.features["image_of_f_center_hz"] = py_round_int(o.fc());
        break;
      }
    }
  }
}

std::vector<Signal> group_images(std::vector<Signal> sigs) {
  // one artifact entry per mirrored source instead of one per fragment
  std::vector<Signal> keep;
  std::vector<int64_t> order;
  std::map<int64_t, std::vector<Signal>> groups;
  for (auto& s : sigs) {
    if (s.has_flag("iq_image")) {
      const int64_t src = s.features["image_of_f_center_hz"].get<int64_t>();
      if (!groups.count(src)) order.push_back(src);
      groups[src].push_back(std::move(s));
    } else {
      keep.push_back(std::move(s));
    }
  }
  for (int64_t src : order) {
    auto& g = groups[src];
    if (g.size() == 1) {
      keep.push_back(std::move(g[0]));
      continue;
    }
    Signal m;
    m.f0 = g[0].f0;
    m.f1 = g[0].f1;
    m.snr_db = g[0].snr_db;
    for (const auto& x : g) {
      m.bursts.insert(m.bursts.end(), x.bursts.begin(), x.bursts.end());
      m.f0 = std::min(m.f0, x.f0);
      m.f1 = std::max(m.f1, x.f1);
      m.snr_db = std::max(m.snr_db, x.snr_db);
    }
    m.flags = {"iq_image"};
    m.features["image_of_f_center_hz"] = src;
    m.features["n_fragments"] = g.size();
    m.detection = g[0].detection;
    m.detection["n_bursts"] = m.bursts.size();
    m.detection["description"] = std::to_string(g.size()) + " energy fragments mirrored about the LO from the signal at " +
                                 fmt("%.3f", static_cast<double>(src) / 1e6) + " MHz, merged into one artifact entry";
    keep.push_back(std::move(m));
  }
  return keep;
}

}  // namespace si::detect
