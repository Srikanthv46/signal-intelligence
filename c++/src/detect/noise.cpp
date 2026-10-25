#include "detect/noise.hpp"

#include <algorithm>
#include <cmath>

#include "core/ndimage.hpp"
#include "core/npcompat.hpp"
#include "core/parallel.hpp"

namespace si::detect {

const double kQ = 20.0;
const double kPlateauSeedDb = 10.0;  // python PLATEAU_SEED_DB
const double kPlateauGrowDb = 3.0;   // python PLATEAU_GROW_DB
const double kQBiasDb = -10.0 * std::log10(-std::log(1.0 - kQ / 100.0));

dvec profile_db(const fmat& P) {
  // 10 * log10(percentile(P, Q, axis=0) + 1e-30) + Q_BIAS_DB: float32 up to the log, then float64
  dvec out(P.cols);
  parallel_for(P.cols, [&](size_t c) {
    fvec col(P.rows);
    for (size_t r = 0; r < P.rows; ++r) col[r] = P(r, c);
    const float p = np::percentile(col.data(), col.size(), kQ);
    const float l = 10.0f * std::log10(p + 1e-30f);
    out[c] = static_cast<double>(l) + kQBiasDb;
  });
  return out;
}

namespace {

// LO leakage sits on the DC bin of every block: bridge it with a straight line (np.linspace)
dvec patch_dc(const dvec& v, size_t half = 2) {
  dvec out = v;
  const size_t c = v.size() / 2;
  const size_t lo = c - half - 1, hi = c + half + 1;
  const size_t num = hi - lo + 1;
  const double step = (v[hi] - v[lo]) / static_cast<double>(num - 1);
  for (size_t i = 1; i + 1 < num; ++i) out[lo + i] = static_cast<double>(i) * step + v[lo];
  return out;
}

}  // namespace

dvec opening_floor(const dvec& prof_db, double df, double win_hz) {
  const size_t k = std::max<size_t>(3, static_cast<size_t>(static_cast<int64_t>(win_hz / df) | 1));
  dvec fl = ndi::grey_opening(patch_dc(prof_db), k, ndi::Mode::Nearest);
  return ndi::uniform_filter1d(fl, k, ndi::Mode::Nearest);
}

// python bridge_plateaus: where the opening sits >= 10 dB above the band's quiet level (10th percentile of the
// opening over the central 2/3) it stands on a wide always-on carrier; across each such run (grown down to 3 dB
// above the quiet level) the floor is interpolated from the quiet bins on either side (np.interp)
std::pair<dvec, size_t> bridge_plateaus(const dvec& op) {
  const size_t n = op.size();
  const double quiet = np::percentile(dvec(op.begin() + n / 6, op.begin() + (n - n / 6)), 10.0);
  bool any_seed = false;
  for (double v : op) any_seed = any_seed || v - quiet > kPlateauSeedDb;
  if (!any_seed) return {op, 0};
  dvec out = op;
  size_t nb = 0, i = 0;
  while (i < n) {
    if (!(op[i] - quiet > kPlateauGrowDb)) {
      ++i;
      continue;
    }
    size_t i1 = i;
    bool seed = false;
    for (; i1 < n && op[i1] - quiet > kPlateauGrowDb; ++i1)
      if (op[i1] - quiet > kPlateauSeedDb) seed = true;
    const size_t i0 = i;
    i = i1;
    if (!seed) continue;
    if (i0 == 0 && i1 == n) return {op, 0};  // no quiet bin to anchor on
    if (i0 == 0) {
      for (size_t k = i0; k < i1; ++k) out[k] = op[i1];
    } else if (i1 == n) {
      for (size_t k = i0; k < i1; ++k) out[k] = op[i0 - 1];
    } else {
      const double x0 = static_cast<double>(i0 - 1), x1 = static_cast<double>(i1);
      const double slope = (op[i1] - op[i0 - 1]) / (x1 - x0);
      for (size_t k = i0; k < i1; ++k) out[k] = slope * (static_cast<double>(k) - x0) + op[i0 - 1];
    }
    nb += i1 - i0;
  }
  return {out, nb};
}

std::pair<dvec, Json> estimate_floor(const dvec& prof_db, double fc, double df) {
  auto [op, nb] = bridge_plateaus(opening_floor(prof_db, df));
  Json info = {{"percentile_over_time", kQ},
               {"percentile_bias_db", py_round(kQBiasDb, 2)},
               {"opening_window_hz", 3e6},
               {"method", "per-snapshot: per-bin 20th percentile over time, morphological opening (3 MHz), "
                          "wide always-on carriers (>= 10 dB above the quiet level) bridged"},
               {"bridged_hz", py_round_int(static_cast<double>(nb) * df)}};
  return {op, info};
}

fvec time_gain(const fmat& P, const dvec& floor_db, const dvec& prof_db, double q, size_t smooth,
               double max_excess_db, size_t min_bins) {
  std::vector<size_t> nb;
  for (size_t i = 0; i < prof_db.size(); ++i)
    if (prof_db[i] - floor_db[i] < max_excess_db) nb.push_back(i);
  if (nb.size() < min_bins) return fvec(P.rows, 1.0f);
  dvec lin(nb.size());
  for (size_t j = 0; j < nb.size(); ++j) lin[j] = std::pow(10.0, floor_db[nb[j]] / 10.0);
  const double bias = -std::log(1.0 - q / 100.0);
  dvec g(P.rows);
  parallel_for(P.rows, [&](size_t r) {
    dvec row(nb.size());
    const float* p = P.row(r);
    for (size_t j = 0; j < nb.size(); ++j) row[j] = static_cast<double>(p[nb[j]]) / lin[j];
    g[r] = np::percentile(row.data(), row.size(), q) / bias;
  });
  g = ndi::uniform_filter1d(g, smooth, ndi::Mode::Nearest);
  fvec out(P.rows);
  for (size_t r = 0; r < P.rows; ++r) out[r] = static_cast<float>(std::max(g[r], 1.0));
  return out;
}

}  // namespace si::detect
