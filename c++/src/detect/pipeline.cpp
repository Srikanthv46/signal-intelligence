#include "detect/pipeline.hpp"

#include <algorithm>
#include <cmath>
#include <set>

#include "core/npcompat.hpp"
#include "core/pyfmt.hpp"
#include "detect/artifacts.hpp"
#include "detect/cfar.hpp"
#include "detect/merge.hpp"
#include "detect/noise.hpp"
#include "detect/stft.hpp"

namespace si::detect {

const size_t kCoarseNfft = 1024;
const size_t kFineNfft = 16384;
const size_t kFineHop = 16384;

namespace {

std::string passes_repr() {
  std::string s = "[";
  for (size_t i = 0; i < kBurstPasses.size(); ++i) {
    const auto& p = kBurstPasses[i];
    if (i) s += ", ";
    s += "(" + std::to_string(p.nt) + ", " + std::to_string(p.nf) + ", " + py_float_repr(p.seed_db) + ", " +
         py_float_repr(p.grow_db) + ")";
  }
  return s + "]";
}

}  // namespace

DetectResult detect(const Capture& cap, double usable_frac) {
  const cfvec& x = cap.iq();
  const double fs = cap.fs, fc = cap.fc;
  DetectResult R;
  Spectrogram C = stft(x, fs, fc, kCoarseNfft);
  const dvec prof = profile_db(C.P);
  auto [floor_c, ninfo] = estimate_floor(prof, fc, C.df());
  Spectrogram F = stft(x, fs, fc, kFineNfft, kFineHop);
  // floor_f = 10 log10(interp(F.f, C.f, 10^(floor_c/10)))
  dvec lin_c(floor_c.size());
  for (size_t i = 0; i < floor_c.size(); ++i) lin_c[i] = std::pow(10.0, floor_c[i] / 10.0);
  dvec floor_f = np::interp(F.f, C.f, lin_c);
  for (auto& v : floor_f) v = 10.0 * std::log10(v);
  const dvec prof_f = profile_db(F.P);

  // strong switched emitters lift the whole floor while on: per-frame gain g(t)
  const fvec g_c = time_gain(C.P, floor_c, prof);
  dvec tc(C.t.size()), tf(F.t.size());
  for (size_t i = 0; i < tc.size(); ++i) tc[i] = C.t[i] + C.frame_len() / 2;
  for (size_t i = 0; i < tf.size(); ++i) tf[i] = F.t[i] + F.frame_len() / 2;
  dvec g_c64(g_c.begin(), g_c.end());
  const dvec g_f64 = np::interp(tf, tc, g_c64);
  fvec g_f(g_f64.size());
  for (size_t i = 0; i < g_f.size(); ++i) g_f[i] = static_cast<float>(g_f64[i]);
  for (size_t r = 0; r < C.P.rows; ++r)
    for (size_t c = 0; c < C.P.cols; ++c) C.P(r, c) /= g_c[r];
  for (size_t r = 0; r < F.P.rows; ++r)
    for (size_t c = 0; c < F.P.cols; ++c) F.P(r, c) /= g_f[r];

  // always-on emitters first, then bursts above max(floor, always-on level)
  auto cc = detect_continuous(C, prof, floor_c, "coarse");
  auto cf_ = detect_continuous(F, prof_f, floor_f, "fine", 1, 3.0);
  dvec base_c(prof.size()), base_f(prof_f.size());
  for (size_t i = 0; i < base_c.size(); ++i) base_c[i] = std::max(floor_c[i], prof[i]);
  for (size_t i = 0; i < base_f.size(); ++i) base_f[i] = std::max(floor_f[i], prof_f[i]);
  auto bc = detect_bursts(C, floor_c, base_c, "coarse");
  auto bf = detect_bursts(F, floor_f, base_f, "fine", {kBurstPasses[0]}, 6, 15.0, 3.0);
  std::vector<Box> boxes = fuse_resolutions(cc, cf_);
  auto fb = fuse_resolutions(bc, bf);
  boxes.insert(boxes.end(), fb.begin(), fb.end());
  boxes = merge_fragments(boxes);
  std::vector<Box> split;
  for (const auto& b : boxes) {
    if (b.res == "coarse") {
      auto parts = split_valleys(b, C, floor_c);
      split.insert(split.end(), parts.begin(), parts.end());
    } else {
      split.push_back(b);
    }
  }
  std::vector<Signal> sigs = absorb_into_continuous(group_signals(split));

  double half = usable_frac * fs / 2;
  if (cap.bandwidth_hz && *cap.bandwidth_hz > 0) half = std::min(half, *cap.bandwidth_hz / 2);  // the analog filter
  const double lo = fc - half, hi = fc + half;
  const std::string thresholds =
      "continuous: 20th-pct-over-time profile > " + py_float_repr(kContSeedDb) + " dB (grow " +
      py_float_repr(kContGrowDb) +
      " dB) above floor; bursts: hysteresis CFAR on P/max(floor, continuous level), passes (t x f cells, seed/grow dB) " +
      passes_repr() + "; fine-res bursts need 3 dB flank contrast";
  for (auto& s : sigs) {
    if (s.f0 < lo || s.f1 > hi) s.flags.push_back("truncated_at_capture_edge");
    std::set<std::string> res, kinds;
    double peak = -1e300;
    for (const auto& b : s.bursts) {
      res.insert(b.res);
      kinds.insert(b.kind);
      peak = std::max(peak, b.peak_snr_db);
    }
    std::string spec = res.count("coarse") ? "coarse 1024-pt Hann, 60 kHz x 16.7 us" : "";
    if (res.size() == 2) spec += " / ";
    if (res.count("fine")) spec += "fine 16384-pt Hann, 3.75 kHz x 267 us";
    s.detection = {{"method", "STFT energy detection, hysteresis CFAR"},
                   {"spectrogram", spec},
                   {"noise_floor", ninfo["method"].get<std::string>() + "; per-frame gain g(t) from noise-only bins"},
                   {"thresholds", thresholds},
                   {"kinds", std::vector<std::string>(kinds.begin(), kinds.end())},
                   {"n_bursts", s.bursts.size()},
                   {"description", std::to_string(s.bursts.size()) +
                                       " time-frequency region(s) rose above the estimated noise floor (peak " +
                                       py_float_repr(peak) +
                                       " dB) and were connected/merged into one emitter by shared frequency span."}};
  }
  flag_lo_leakage(sigs, fc, F.df());
  flag_iq_images(sigs, fc);
  sigs = group_images(std::move(sigs));
  std::stable_sort(sigs.begin(), sigs.end(),
                   [](const Signal& a, const Signal& b) { return a.f0 != b.f0 ? a.f0 < b.f0 : a.t0() < b.t0(); });

  // time_gain_db_p95, noise_floor_dbfs_per_bin_median
  fvec gcopy = g_c;
  const float p95 = np::percentile(gcopy.data(), gcopy.size(), 95.0);
  ninfo["time_gain_db_p95"] = py_round(static_cast<double>(10.0f * std::log10(p95)), 2);
  R.info = {{"noise", ninfo},
            {"n_continuous", Json::array({cc.size(), cf_.size()})},
            {"n_bursts", Json::array({bc.size(), bf.size()})},
            {"noise_floor_dbfs_per_bin_median", py_round(np::median(floor_c.data(), floor_c.size()), 1)}};
  R.signals = std::move(sigs);
  R.coarse = std::move(C);
  R.floor_coarse_db = std::move(floor_c);
  R.fine = std::move(F);
  R.floor_fine_db = std::move(floor_f);
  R.usable_hz = {lo, hi};
  return R;
}

}  // namespace si::detect
