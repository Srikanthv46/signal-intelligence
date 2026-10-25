#include "detect/merge.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "core/ndimage.hpp"
#include "core/npcompat.hpp"

namespace si::detect {

const double kNarrowHz = 300e3;

namespace {

double overlap(double a0, double a1, double b0, double b1) { return std::max(0.0, std::min(a1, b1) - std::max(a0, b0)); }

double median_of(std::vector<double> v) { return np::median(v.data(), v.size()); }

}  // namespace

std::vector<Box> fuse_resolutions(const std::vector<Box>& coarse, const std::vector<Box>& fine) {
  std::vector<Box> wide, fine_n, out;
  for (const auto& b : coarse)
    if (b.bw() >= kNarrowHz) wide.push_back(b);
  out = wide;
  for (const auto& b : fine)
    if (b.bw() < kNarrowHz) fine_n.push_back(b);
  auto inside_wide = [&](const Box& b) {
    for (const auto& w : wide)
      if (overlap(b.f0, b.f1, w.f0, w.f1) > 0.5 * b.bw() && overlap(b.t0, b.t1, w.t0, w.t1) > 0.5 * b.dur())
        return true;
    return false;
  };
  for (const auto& b : fine_n)
    if (!inside_wide(b)) out.push_back(b);
  for (const auto& b : coarse) {
    if (!(b.bw() < kNarrowHz) || inside_wide(b)) continue;
    bool seen_fine = false;
    for (const auto& f : fine_n)
      if (overlap(b.f0 - 10e3, b.f1 + 10e3, f.f0, f.f1) > 0 && overlap(b.t0, b.t1, f.t0 - 300e-6, f.t1 + 300e-6) > 0) {
        seen_fine = true;
        break;
      }
    if (!seen_fine) out.push_back(b);
  }
  return out;
}

std::vector<Box> merge_fragments(std::vector<Box> boxes, double f_gap, double t_gap) {
  bool changed = true;
  while (changed) {
    changed = false;
    std::stable_sort(boxes.begin(), boxes.end(),
                     [](const Box& a, const Box& b) { return a.f0 != b.f0 ? a.f0 < b.f0 : a.t0 < b.t0; });
    std::vector<Box> out;
    for (const auto& b : boxes) {
      bool merged = false;
      for (auto& a : out) {
        const double tov = overlap(a.t0, a.t1, b.t0, b.t1);
        const double fov = overlap(a.f0, a.f1, b.f0, b.f1);
        if (a.kind != b.kind) continue;
        const bool side = tov > 0.5 * std::min(a.dur(), b.dur()) && std::max(a.f0, b.f0) - std::min(a.f1, b.f1) <= f_gap &&
                          a.res == b.res;
        bool seq = fov > 0.8 * std::min(a.bw(), b.bw()) && std::max(a.t0, b.t0) - std::min(a.t1, b.t1) <= t_gap;
        if (seq) {
          const double r = a.bw() / b.bw();
          seq = 0.6 < r && r < 1.6;
        }
        if (side || seq || (tov > 0 && fov > 0.5 * std::min(a.bw(), b.bw()))) {
          Box m;
          m.t0 = std::min(a.t0, b.t0);
          m.t1 = std::max(a.t1, b.t1);
          m.f0 = std::min(a.f0, b.f0);
          m.f1 = std::max(a.f1, b.f1);
          m.peak_snr_db = std::max(a.peak_snr_db, b.peak_snr_db);
          m.mean_snr_db = a.n_cells >= b.n_cells ? a.mean_snr_db : b.mean_snr_db;
          m.res = a.bw() >= b.bw() ? a.res : b.res;
          m.n_cells = a.n_cells + b.n_cells;
          m.kind = a.kind;
          a = m;
          changed = true;
          merged = true;
          break;
        }
      }
      if (!merged) out.push_back(b);
    }
    boxes = std::move(out);
  }
  return boxes;
}

std::vector<Box> split_valleys(const Box& b, const Spectrogram& S, const dvec& floor_db, double depth_db,
                               double min_part_hz) {
  if (b.bw() < 2 * min_part_hz) return {b};
  const double df = S.df();
  const size_t i0 = static_cast<size_t>(S.f_index(b.f0 + df / 2));
  const size_t i1 = static_cast<size_t>(S.f_index(b.f1 - df / 2)) + 1;
  const size_t j0 = static_cast<size_t>(S.t_index(b.t0));
  const size_t j1 = static_cast<size_t>(S.t_index(b.t1 - S.frame_len())) + 1;
  if (i1 <= i0 || j1 <= j0) return {b};
  // P[j0:j1, i0:i1].mean(0): float32, rows added in order (reduction along a non-contiguous axis)
  const size_t n = i1 - i0;
  fvec colsum(n, 0.0f);
  for (size_t r = j0; r < j1; ++r) {
    const float* p = S.P.row(r) + i0;
    for (size_t c = 0; c < n; ++c) colsum[c] += p[c];
  }
  dvec raw(n);
  const float cnt = static_cast<float>(j1 - j0);
  for (size_t c = 0; c < n; ++c) {
    const float m = colsum[c] / cnt;
    raw[c] = static_cast<double>(10.0f * std::log10(m + 1e-30f)) - floor_db[i0 + c];
  }
  const dvec prof = ndi::uniform_filter1d(raw, 3, ndi::Mode::Reflect);
  const size_t kmin = std::max<size_t>(2, static_cast<size_t>(min_part_hz / df));
  auto max_of = [&](size_t a, size_t c) {
    double m = -std::numeric_limits<double>::infinity();
    for (size_t i = a; i < c; ++i) m = std::max(m, prof[i]);
    return m;
  };
  std::vector<size_t> cuts;
  for (size_t k = kmin; k + kmin < prof.size(); ++k) {
    const size_t a = k >= 2 ? k - 2 : 0;
    double mn = std::numeric_limits<double>::infinity();
    for (size_t i = a; i < std::min(prof.size(), k + 3); ++i) mn = std::min(mn, prof[i]);
    if (prof[k] != mn) continue;
    const double left = max_of(0, k), right = max_of(k + 1, prof.size());
    const double lo0 = max_of(k >= 2 * kmin ? k - 2 * kmin : 0, k);
    const double lo1 = max_of(k + 1, std::min(prof.size(), k + 1 + 2 * kmin));
    if (std::min(lo0, lo1) - prof[k] >= depth_db && std::min(left, right) - prof[k] >= depth_db)
      if (cuts.empty() || k - cuts.back() >= kmin) cuts.push_back(k);
  }
  if (cuts.empty()) return {b};
  std::vector<size_t> edges{0};
  edges.insert(edges.end(), cuts.begin(), cuts.end());
  edges.push_back(prof.size());
  std::vector<Box> out;
  for (size_t e = 0; e + 1 < edges.size(); ++e) {
    const size_t a = edges[e], c = edges[e + 1];
    Box x;
    x.t0 = b.t0;
    x.t1 = b.t1;
    x.f0 = S.f[i0 + a] - df / 2;
    x.f1 = S.f[i0 + c - 1] + df / 2;
    x.peak_snr_db = py_round(max_of(a, c), 1);
    x.mean_snr_db = py_round(median_of(dvec(prof.begin() + a, prof.begin() + c)), 1);
    x.res = b.res;
    x.n_cells = b.n_cells * static_cast<int64_t>(c - a) / static_cast<int64_t>(prof.size());
    x.kind = b.kind;
    out.push_back(x);
  }
  return out;
}

std::vector<Signal> group_signals(std::vector<Box> boxes, double iou_min) {
  std::stable_sort(boxes.begin(), boxes.end(), [](const Box& a, const Box& b) { return a.t0 < b.t0; });
  struct Group {
    std::vector<Box> boxes;
    double f0, f1;  // median spans, updated on append
  };
  std::vector<Group> groups;
  for (const auto& b : boxes) {
    Group* best = nullptr;
    double best_iou = iou_min;
    for (auto& g : groups) {
      const double inter = overlap(g.f0, g.f1, b.f0, b.f1);
      const double iou = inter / (std::max(g.f1, b.f1) - std::min(g.f0, b.f0));
      if (iou >= best_iou) {
        best = &g;
        best_iou = iou;
      }
    }
    if (!best) {
      groups.push_back({{b}, b.f0, b.f1});
    } else {
      best->boxes.push_back(b);
      dvec f0s, f1s;
      for (const auto& x : best->boxes) {
        f0s.push_back(x.f0);
        f1s.push_back(x.f1);
      }
      best->f0 = median_of(f0s);
      best->f1 = median_of(f1s);
    }
  }
  std::vector<Signal> sigs;
  for (auto& g : groups) {
    dvec w, v;
    for (const auto& x : g.boxes) {
      w.push_back(std::max(x.dur(), 1e-6));
      v.push_back(std::pow(10.0, x.mean_snr_db / 10.0));
    }
    dvec vw(v.size());
    for (size_t i = 0; i < v.size(); ++i) vw[i] = v[i] * w[i];
    const double avg = np::sum(vw.data(), vw.size()) / np::sum(w.data(), w.size());  // np.average(v, weights=w)
    Signal s;
    s.bursts = std::move(g.boxes);
    s.f0 = g.f0;
    s.f1 = g.f1;
    s.snr_db = py_round(10.0 * std::log10(avg), 1);
    sigs.push_back(std::move(s));
  }
  return sigs;
}

std::vector<Signal> absorb_into_continuous(std::vector<Signal> sigs, double frac) {
  auto is_cont = [](const Signal& s) {
    return std::any_of(s.bursts.begin(), s.bursts.end(), [](const Box& b) { return b.kind == "continuous"; });
  };
  std::vector<size_t> cont;
  for (size_t i = 0; i < sigs.size(); ++i)
    if (is_cont(sigs[i])) cont.push_back(i);
  std::vector<uint8_t> keep(sigs.size(), 1);
  for (size_t i = 0; i < sigs.size(); ++i) {
    if (is_cont(sigs[i])) continue;
    const Signal& s = sigs[i];
    for (size_t c : cont) {
      Signal& h = sigs[c];
      if (std::max(0.0, std::min(s.f1, h.f1) - std::max(s.f0, h.f0)) >= frac * s.bw()) {
        h.bursts.insert(h.bursts.end(), s.bursts.begin(), s.bursts.end());
        const int64_t prev = h.features.value("n_bursts_above_continuous_level", int64_t(0));
        h.features["n_bursts_above_continuous_level"] = prev + static_cast<int64_t>(s.bursts.size());
        keep[i] = 0;
        break;
      }
    }
  }
  std::vector<Signal> out;
  for (size_t i = 0; i < sigs.size(); ++i)
    if (keep[i]) out.push_back(std::move(sigs[i]));
  return out;
}

namespace {
constexpr double kPtPriTol = 0.01;        // pieces of one pulse train: same PRI within 1 %
constexpr double kPtTTolS = 5e-6;         // and the same pulse starts within 5 us
constexpr double kPtMatchFrac = 0.5;      // for at least half of the shorter train's pulses
constexpr double kPtFGapHz = 2e6;         // spans at most 2 MHz apart
constexpr double kPtSliverPadS = 20e-6;   // a sliver's every burst holds a pulse start (+-20 us) of the train

bool truthy_num(const Json& f, const char* k) {
  return f.contains(k) && f[k].is_number() && f[k].get<double>() != 0;
}
dvec pulse_times(const Signal& s) {
  dvec t;
  if (s.features.contains("_pulse_t_s") && s.features["_pulse_t_s"].is_array())
    for (const auto& v : s.features["_pulse_t_s"]) t.push_back(v.get<double>());
  return t;
}
double coincide(const dvec& ta, const dvec& tb) {
  const dvec& a = ta.size() <= tb.size() ? ta : tb;
  const dvec& b = ta.size() <= tb.size() ? tb : ta;
  if (a.empty()) return 0.0;
  size_t n = 0;
  for (double t : a)
    for (double u : b)
      if (std::fabs(u - t) <= kPtTTolS) {
        ++n;
        break;
      }
  return static_cast<double>(n) / static_cast<double>(a.size());
}
}  // namespace

std::pair<std::vector<Signal>, std::vector<size_t>> merge_pulse_trains(std::vector<Signal> sigs) {
  std::vector<size_t> hosts;
  for (size_t i = 0; i < sigs.size(); ++i)
    if (truthy_num(sigs[i].features, "pri_s") && !pulse_times(sigs[i]).empty()) hosts.push_back(i);
  std::vector<std::vector<size_t>> groups;
  for (size_t i : hosts) {
    const Signal& s = sigs[i];
    const double pri_s = s.features["pri_s"].get<double>();
    std::vector<size_t>* hit = nullptr;
    for (auto& g : groups) {
      for (size_t j : g) {
        const Signal& h = sigs[j];
        const double pri_h = h.features["pri_s"].get<double>();
        if (std::fabs(pri_s - pri_h) <= kPtPriTol * pri_h && std::max(s.f0, h.f0) - std::min(s.f1, h.f1) <= kPtFGapHz &&
            coincide(pulse_times(s), pulse_times(h)) >= kPtMatchFrac) {
          hit = &g;
          break;
        }
      }
      if (hit) break;
    }
    if (hit)
      hit->push_back(i);
    else
      groups.push_back({i});
  }
  std::vector<bool> used(sigs.size(), false);
  std::vector<Signal> merged;
  for (const auto& g : groups) {
    double f0 = std::numeric_limits<double>::infinity(), f1 = -f0;
    dvec t;
    for (size_t j : g) {
      f0 = std::min(f0, sigs[j].f0);
      f1 = std::max(f1, sigs[j].f1);
      for (double v : pulse_times(sigs[j])) t.push_back(v);
    }
    std::sort(t.begin(), t.end());
    std::vector<size_t> slivers;
    for (size_t i = 0; i < sigs.size(); ++i) {
      const Signal& s = sigs[i];
      if (std::find(g.begin(), g.end(), i) != g.end() || used[i] || truthy_num(s.features, "pri_s")) continue;
      if (std::any_of(s.bursts.begin(), s.bursts.end(), [](const Box& b) { return b.kind == "continuous"; })) continue;
      if (!(s.f0 >= f0 - kPtFGapHz && s.f1 <= f1 + kPtFGapHz)) continue;
      const bool on = std::all_of(s.bursts.begin(), s.bursts.end(), [&](const Box& b) {
        return std::any_of(t.begin(), t.end(),
                           [&](double v) { return v >= b.t0 - kPtSliverPadS && v <= b.t1 + kPtSliverPadS; });
      });
      if (on) slivers.push_back(i);
    }
    if (g.size() + slivers.size() < 2) continue;
    std::vector<size_t> parts = g;
    parts.insert(parts.end(), slivers.begin(), slivers.end());
    size_t top = g[0];
    for (size_t j : g)
      if (sigs[j].snr_db > sigs[top].snr_db) top = j;
    Signal m;
    for (size_t p : parts) m.bursts.insert(m.bursts.end(), sigs[p].bursts.begin(), sigs[p].bursts.end());
    m.f0 = std::numeric_limits<double>::infinity();
    m.f1 = -m.f0;
    for (size_t p : parts) {
      m.f0 = std::min(m.f0, sigs[p].f0);
      m.f1 = std::max(m.f1, sigs[p].f1);
    }
    m.snr_db = sigs[top].snr_db;
    m.detection = sigs[top].detection;
    const std::string d0 = m.detection.contains("description") ? m.detection["description"].get<std::string>() : "";
    m.detection["description"] = d0 + "; " + std::to_string(parts.size()) +
                                 " energy detections with one pulse train (same PRI and pulse times) merged into one "
                                 "emitter";
    for (size_t p : parts) used[p] = true;
    merged.push_back(std::move(m));
  }
  std::vector<Signal> out;
  for (size_t i = 0; i < sigs.size(); ++i)
    if (!used[i]) out.push_back(std::move(sigs[i]));
  std::vector<size_t> idx;
  for (auto& m : merged) {
    idx.push_back(out.size());
    out.push_back(std::move(m));
  }
  return {std::move(out), idx};
}

}  // namespace si::detect
