#include "detect/cfar.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "core/ndimage.hpp"
#include "core/npcompat.hpp"

namespace si::detect {

const std::vector<BurstPass> kBurstPasses = {{3, 3, 6.0, 3.0}, {15, 5, 3.0, 1.5}};
const double kContSeedDb = 3.0, kContGrowDb = 1.5;

namespace {

std::vector<std::pair<size_t, size_t>> runs(const std::vector<uint8_t>& mask) {
  std::vector<std::pair<size_t, size_t>> out;
  size_t i = 0;
  while (i < mask.size()) {
    if (!mask[i]) {
      ++i;
      continue;
    }
    size_t j = i;
    while (j < mask.size() && mask[j]) ++j;
    out.push_back({i, j});
    i = j;
  }
  return out;
}

double slice_mean(const dvec& v, size_t a, size_t b) { return np::mean(v.data() + a, b - a); }

}  // namespace

double mean2d(const dmat& m, size_t r0, size_t r1, size_t c0, size_t c1) {
  double acc = 0.0;
  for (size_t r = r0; r < r1; ++r) acc += np::sum(m.row(r) + c0, c1 - c0);
  return acc / static_cast<double>((r1 - r0) * (c1 - c0));
}

std::vector<Box> detect_continuous(const Spectrogram& S, const dvec& prof_db, const dvec& floor_db,
                                   const std::string& res, size_t min_bins, std::optional<double> flank_db) {
  dvec diff(prof_db.size());
  for (size_t i = 0; i < diff.size(); ++i) diff[i] = prof_db[i] - floor_db[i];
  const dvec excess = ndi::uniform_filter1d(diff, 3, ndi::Mode::Reflect);
  std::vector<uint8_t> seed(excess.size()), grow(excess.size());
  for (size_t i = 0; i < excess.size(); ++i) {
    seed[i] = excess[i] > kContSeedDb;
    grow[i] = excess[i] > kContGrowDb;
  }
  const double T = S.t.back() + S.frame_len();
  const double df = S.df();
  std::vector<Box> out;
  for (const auto& [a, b] : runs(grow)) {
    if (b - a < min_bins || !std::any_of(seed.begin() + a, seed.begin() + b, [](uint8_t v) { return v; })) continue;
    if (flank_db) {
      // narrow tone test: the run must stand above the profile beside it
      const size_t w = std::max<size_t>(b - a, 4);
      const double ninf = -std::numeric_limits<double>::infinity();
      const int64_t lo = static_cast<int64_t>(a) - 1 - static_cast<int64_t>(w);
      const double fl0 = lo >= 0 ? slice_mean(excess, static_cast<size_t>(lo), a - 1) : ninf;
      const double fl1 = b + 1 + w <= excess.size() ? slice_mean(excess, b + 1, b + 1 + w) : ninf;
      if (slice_mean(excess, a, b) - std::max(fl0, fl1) < *flank_db) continue;
    }
    dvec lin(b - a);
    double mx = -std::numeric_limits<double>::infinity();
    for (size_t i = a; i < b; ++i) {
      lin[i - a] = std::pow(10.0, excess[i] / 10.0);
      mx = std::max(mx, excess[i]);
    }
    Box box;
    box.t0 = 0.0;
    box.t1 = T;
    box.f0 = S.f[a] - df / 2;
    box.f1 = S.f[b - 1] + df / 2;
    box.peak_snr_db = py_round(mx, 1);
    box.mean_snr_db = py_round(10.0 * std::log10(std::max(np::mean(lin.data(), lin.size()) - 1.0, 1e-3)), 1);
    box.res = res;
    box.n_cells = static_cast<int64_t>(b - a);
    box.kind = "continuous";
    out.push_back(box);
  }
  return out;
}

namespace {

// Narrowband test: the box must stand min_db above the spectrum on both sides over the
// same time span (rejects fragments of a weak wideband signal).
bool flank_ok(const dmat& sm, size_t t0, size_t t1, size_t f0, size_t f1, double min_db) {
  const int64_t w = std::max<int64_t>(static_cast<int64_t>(f1 - f0), 4);
  const int64_t F = static_cast<int64_t>(sm.cols);
  const double inner = mean2d(sm, t0, t1, f0, f1);
  double best = -std::numeric_limits<double>::infinity();
  bool any = false;
  const int64_t spans[2][2] = {{static_cast<int64_t>(f0) - 1 - w, static_cast<int64_t>(f0) - 1},
                               {static_cast<int64_t>(f1) + 1, static_cast<int64_t>(f1) + 1 + w}};
  for (const auto& sp : spans) {
    const int64_t a = std::max<int64_t>(0, sp[0]), b = std::min<int64_t>(F, sp[1]);
    if (b - a >= 2) {
      best = std::max(best, mean2d(sm, t0, t1, static_cast<size_t>(a), static_cast<size_t>(b)));
      any = true;
    }
  }
  return !any || inner >= std::pow(10.0, min_db / 10.0) * best;
}

}  // namespace

std::vector<Box> detect_bursts(const Spectrogram& S, const dvec& floor_db, const dvec& baseline_db,
                               const std::string& res, const std::vector<BurstPass>& passes, size_t min_cells,
                               double strong_db, std::optional<double> flank_db) {
  const size_t T = S.P.rows, F = S.P.cols;
  dvec base_lin(F), floor_lin(F);
  for (size_t c = 0; c < F; ++c) {
    base_lin[c] = std::pow(10.0, baseline_db[c] / 10.0);
    floor_lin[c] = std::pow(10.0, floor_db[c] / 10.0);
  }
  dmat rel(T, F), snr(T, F);  // float32 P / float64 baseline -> float64
  for (size_t r = 0; r < T; ++r)
    for (size_t c = 0; c < F; ++c) {
      rel(r, c) = static_cast<double>(S.P(r, c)) / base_lin[c];
      snr(r, c) = static_cast<double>(S.P(r, c)) / floor_lin[c];
    }
  bmat seed(T, F, 0), grow(T, F, 0);
  dmat sm0;
  for (size_t p = 0; p < passes.size(); ++p) {
    dmat sm = ndi::uniform_filter(rel, passes[p].nt, passes[p].nf, ndi::Mode::Nearest);
    const double sd = std::pow(10.0, passes[p].seed_db / 10.0), gr = std::pow(10.0, passes[p].grow_db / 10.0);
    for (size_t i = 0; i < sm.size(); ++i) {
      seed.data[i] |= sm.data[i] > sd;
      grow.data[i] |= sm.data[i] > gr;
    }
    if (p == 0) sm0 = std::move(sm);
  }
  auto [lab, n] = ndi::label(grow, ndi::ones(3, 3));
  if (n == 0) return {};
  std::vector<uint8_t> has_seed(static_cast<size_t>(n) + 1, 0);
  for (size_t i = 0; i < lab.size(); ++i)
    if (seed.data[i]) has_seed[static_cast<size_t>(lab.data[i])] = 1;
  has_seed[0] = 0;
  const dmat snr_floor = ndi::uniform_filter(snr, 3, 3, ndi::Mode::Nearest);
  const auto objs = ndi::find_objects(lab, n);
  const double df = S.df(), flen = S.frame_len();
  std::vector<Box> boxes;
  for (int32_t k = 1; k <= n; ++k) {
    if (!has_seed[static_cast<size_t>(k)]) continue;
    const auto& o = objs[static_cast<size_t>(k - 1)];
    int64_t cells = 0;
    double pk = -std::numeric_limits<double>::infinity();
    dvec lin;
    for (size_t r = o.r0; r < o.r1; ++r)
      for (size_t c = o.c0; c < o.c1; ++c)
        if (lab(r, c) == k) {
          ++cells;
          pk = std::max(pk, sm0(r, c));
          lin.push_back(snr_floor(r, c));
        }
    const double pk_rel = 10.0 * std::log10(pk);
    if (cells < static_cast<int64_t>(min_cells) && pk_rel < strong_db) continue;
    if (flank_db && !flank_ok(sm0, o.r0, o.r1, o.c0, o.c1, *flank_db)) continue;
    const double lmax = *std::max_element(lin.begin(), lin.end());
    Box b;
    b.t0 = S.t[o.r0];
    b.t1 = S.t[o.r1 - 1] + flen;
    b.f0 = S.f[o.c0] - df / 2;
    b.f1 = S.f[o.c1 - 1] + df / 2;
    b.peak_snr_db = py_round(10.0 * std::log10(lmax), 1);
    b.mean_snr_db = py_round(10.0 * std::log10(std::max(np::mean(lin.data(), lin.size()) - 1.0, 1e-3)), 1);
    b.res = res;
    b.n_cells = cells;
    boxes.push_back(b);
  }
  return boxes;
}

}  // namespace si::detect
