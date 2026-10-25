// python: core/bursts.py (TFPower, rect_fit, fit_band, energy_bursts).
#include "core/bursts.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "core/fft.hpp"
#include "core/ndimage.hpp"
#include "core/npcompat.hpp"

namespace si::bursts {

TFPower::TFPower(const cfvec& x, double fs_, double frame_s) : fs(fs_) {
  L = static_cast<size_t>(std::max<int64_t>(16, static_cast<int64_t>(std::pow(2.0, static_cast<double>(py_round_int(std::log2(fs * frame_s)))))));
  n = x.size() / L;
  const dvec w = np::hanning(L);
  P = dmat(n, L);
  cdvec v(L);
  const size_t h = L / 2;
  for (size_t k = 0; k < n; ++k) {
    for (size_t i = 0; i < L; ++i) {
      const cf& s = x[k * L + i];
      v[i] = cd(static_cast<double>(s.real()) * w[i], static_cast<double>(s.imag()) * w[i]);
    }
    const cdvec F = fft::fft(v);
    double* row = P.row(k);
    for (size_t j = 0; j < L; ++j) {
      const double a = np::cabs(F[(j + L - h) % L]);
      row[j] = a * a;
    }
  }
  dt = static_cast<double>(L) / fs;
  f.resize(L);
  for (size_t j = 0; j < L; ++j) f[j] = static_cast<double>(static_cast<int64_t>(j) - static_cast<int64_t>(h)) * (fs / static_cast<double>(L));
}

std::vector<size_t> TFPower::cols(double lo, double hi) const {
  std::vector<size_t> c;
  for (size_t j = 0; j < L; ++j)
    if (f[j] >= lo && f[j] <= hi) c.push_back(j);
  return c;
}

std::vector<uint8_t> TFPower::frames(const Runs& runs) const {
  std::vector<uint8_t> on(n, 0);
  for (const auto& [t0, t1] : runs) {
    const int64_t a = static_cast<int64_t>(std::ceil(t0 / dt - 1e-9));
    const int64_t b = static_cast<int64_t>(std::floor(t1 / dt + 1e-9));
    for (int64_t k = std::max<int64_t>(a, 0); k < std::min<int64_t>(b, static_cast<int64_t>(n)); ++k) on[k] = 1;
  }
  return on;
}

std::vector<uint8_t> TFPower::touched(const Runs& runs, int64_t pad) const {
  std::vector<uint8_t> m(n, 0);
  for (const auto& [t0, t1] : runs) {
    const int64_t a = static_cast<int64_t>(std::floor(t0 / dt)) - pad;
    const int64_t b = static_cast<int64_t>(std::ceil(t1 / dt)) + pad;
    for (int64_t k = std::max<int64_t>(a, 0); k < std::min<int64_t>(b, static_cast<int64_t>(n)); ++k) m[k] = 1;
  }
  return m;
}

RectFit rect_fit(const dvec& d) {
  const size_t n = d.size();
  const dvec sm = ndi::uniform_filter1d(d, 5);
  double tau = *std::max_element(sm.begin(), sm.end()) / 2;
  RectFit r{0, n, 0.0, 0.0};
  for (int it = 0; it < 4; ++it) {
    double best = -std::numeric_limits<double>::infinity(), s = 0.0;
    size_t a = 0;
    for (size_t i = 0; i < n; ++i) {
      if (s <= 0) s = 0.0, a = i;
      s += d[i] - tau;
      if (s > best) best = s, r.lo = a, r.hi = i + 1;
    }
    r.inside = np::mean(d.data() + r.lo, r.hi - r.lo);
    dvec rest(d.begin(), d.begin() + r.lo);
    rest.insert(rest.end(), d.begin() + r.hi, d.end());
    r.outside = rest.empty() ? 0.0 : np::median(rest);
    tau = (r.inside + r.outside) / 2;
  }
  return r;
}

namespace {

dmat select(const TFPower& tf, const std::vector<uint8_t>& rows, const std::vector<size_t>& c) {
  size_t nr = 0;
  for (auto v : rows) nr += v;
  dmat m(nr, c.size());
  size_t r = 0;
  for (size_t k = 0; k < tf.n; ++k) {
    if (!rows[k]) continue;
    const double* p = tf.P.row(k);
    for (size_t j = 0; j < c.size(); ++j) m(r, j) = p[c[j]];
    ++r;
  }
  return m;
}

double lower_sigma(const dvec& v) { return (np::median(v) - np::percentile(v, 15.87)) + 1e-12; }

dvec pick(const dvec& v, const std::vector<uint8_t>& m) {
  dvec o;
  for (size_t i = 0; i < v.size(); ++i)
    if (m[i]) o.push_back(v[i]);
  return o;
}

// np.percentile(v, 25, method="lower"): the element of rank floor(0.25 (n - 1))
double lower_q25(dvec& v) {
  const size_t i = static_cast<size_t>(std::floor(0.25 * static_cast<double>(v.size() - 1)));
  std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(i), v.end());
  return v[i];
}

// mean over the contiguous column range [c0, c0 + len) of each frame
dvec band_mean(const TFPower& tf, size_t c0, size_t len) {
  dvec g(tf.n);
  for (size_t k = 0; k < tf.n; ++k) g[k] = np::mean(tf.P.row(k) + c0, len);
  return g;
}

}  // namespace

std::optional<Band> fit_band(const TFPower& tf, const Runs& runs, double lo, double hi, double search_lo,
                             double search_hi, double smooth_hz) {
  const auto on = tf.frames(runs);
  auto idle = tf.touched(runs, 2);
  for (auto& v : idle) v = !v;
  size_t non = 0, nidle = 0;
  for (size_t k = 0; k < tf.n; ++k) non += on[k], nidle += idle[k];
  if (non < 3 || nidle < 3) return std::nullopt;
  const auto c = tf.cols(search_lo, search_hi);
  if (c.size() < 8) return std::nullopt;
  const dvec A = np::median_axis0(select(tf, on, c));
  const dvec I = np::median_axis0(select(tf, idle, c));
  const size_t k = static_cast<size_t>(std::max<int64_t>(1, py_round_int(smooth_hz / (tf.fs / static_cast<double>(tf.L)))));
  const dvec Af = ndi::uniform_filter1d(A, k), If = ndi::uniform_filter1d(I, k);
  dvec d(c.size());
  for (size_t i = 0; i < c.size(); ++i) d[i] = 10 * std::log10(std::max(Af[i], 1e-30) / std::max(If[i], 1e-30));
  const RectFit r = rect_fit(d);
  const double cen = (lo + hi) / 2;
  const double fa = tf.f[c[r.lo]], fb = tf.f[c[r.hi - 1]];
  if (!(fa <= cen && cen <= fb)) return std::nullopt;
  dvec rest(d.begin(), d.begin() + r.lo);
  rest.insert(rest.end(), d.begin() + r.hi, d.end());
  double sd = 1e-9;
  if (rest.size() >= 8) {
    const double m = np::median(rest);
    dvec ad(rest.size());
    for (size_t i = 0; i < rest.size(); ++i) ad[i] = std::fabs(rest[i] - m);
    sd = 1.4826 * np::median(ad) + 1e-9;
  }
  const double z = (r.inside - r.outside) / sd * std::sqrt(static_cast<double>(r.hi - r.lo) / static_cast<double>(k));
  if (r.inside - r.outside < 1.0 || z < 6.0) return std::nullopt;
  const double df = tf.fs / static_cast<double>(tf.L);
  return Band{fa - df / 2, fb + df / 2, py_round(r.inside - r.outside, 2), z};
}

std::pair<Runs, EnergyInfo> energy_bursts(const TFPower& tf, Runs seeds, double lo, double hi, double usable_half,
                                          double min_sym_s) {
  std::sort(seeds.begin(), seeds.end());
  EnergyInfo info;
  if (seeds.empty() || hi - lo <= 0) return {seeds, info};
  // sub-band powers (python _subband_ratio): raw and time-smoothed
  const auto c = tf.cols(lo, hi);
  const size_t nsb = static_cast<size_t>(std::clamp<int64_t>(py_round_int((hi - lo) / kSbHz), 2, kMaxSb));
  std::vector<std::pair<size_t, size_t>> groups;  // np.array_split: (start, len)
  {
    const size_t q = c.size() / nsb, rem = c.size() % nsb;
    size_t s = 0;
    for (size_t g = 0; g < nsb; ++g) {
      const size_t len = q + (g < rem ? 1 : 0);
      groups.emplace_back(s, len);
      s += len;
    }
  }
  dmat raw(tf.n, nsb);
  for (size_t g = 0; g < nsb; ++g) {
    const dvec v = groups[g].second ? band_mean(tf, c[groups[g].first], groups[g].second) : dvec(tf.n, std::nan(""));
    for (size_t k = 0; k < tf.n; ++k) raw(k, g) = v[k];
  }
  const dmat sm = ndi::uniform_filter1d(raw, kSmooth, 0);
  const auto seed_on = tf.frames(seeds);
  auto away = tf.touched(seeds, 2);
  for (auto& v : away) v = !v;
  size_t naway = 0, non = 0;
  for (size_t k = 0; k < tf.n; ++k) naway += away[k], non += seed_on[k];
  if (static_cast<double>(naway) < 0.2 * static_cast<double>(tf.n) || non < 3) return {seeds, info};
  dmat sa(naway, nsb);
  for (size_t k = 0, r = 0; k < tf.n; ++k)
    if (away[k]) {
      for (size_t g = 0; g < nsb; ++g) sa(r, g) = sm(k, g);
      ++r;
    }
  dvec floor = np::percentile_axis0(sa, 25.0);
  for (auto& v : floor) v = std::max(v, 1e-30);
  dmat q(tf.n, nsb), qr(tf.n, nsb);
  for (size_t k = 0; k < tf.n; ++k)
    for (size_t g = 0; g < nsb; ++g) q(k, g) = sm(k, g) / floor[g] - 1.0, qr(k, g) = raw(k, g) / floor[g] - 1.0;
  // lit guard bands
  std::vector<std::vector<uint8_t>> lit;
  const double gw = std::max(kGuardHz, 0.1 * (hi - lo));
  const std::pair<double, double> gb[2] = {{lo - kGuardGapHz - gw, lo - kGuardGapHz}, {hi + kGuardGapHz, hi + kGuardGapHz + gw}};
  for (const auto& [a, b] : gb) {
    const auto gc = tf.cols(a, b);
    if (gc.size() < 3 || std::fabs(b) > usable_half || std::fabs(a) > usable_half) {
      lit.emplace_back(tf.n, 0);
      continue;
    }
    const dvec gs = ndi::uniform_filter1d(band_mean(tf, gc[0], gc.size()), kSmooth);
    const double gf = std::max(np::percentile(pick(gs, away), 25.0), 1e-30);
    dvec gr(tf.n);
    for (size_t k = 0; k < tf.n; ++k) gr[k] = gs[k] / gf - 1.0;
    const double ls = 6 * lower_sigma(pick(gr, away));
    std::vector<uint8_t> l(tf.n);
    for (size_t k = 0; k < tf.n; ++k) l[k] = gr[k] > ls;
    lit.push_back(std::move(l));
  }
  const size_t e = (nsb + 3) / 4;
  dvec r(tf.n), rr(tf.n), med(tf.n);
  dvec buf;
  for (size_t k = 0; k < tf.n; ++k) {
    size_t kl = lit[0][k] ? e : 0, kh = lit[1][k] ? nsb - e : nsb;
    if (kh < kl + 1) kl = 0, kh = nsb;
    buf.assign(q.row(k) + kl, q.row(k) + kh);
    r[k] = lower_q25(buf);
    buf.assign(qr.row(k) + kl, qr.row(k) + kh);
    rr[k] = lower_q25(buf);
    buf.assign(q.row(k), q.row(k) + nsb);
    med[k] = np::median(buf);
  }
  const double level = np::median(pick(r, seed_on));
  const double sig = lower_sigma(pick(r, away));
  info.has_level = true;
  info.level_db = py_round(10 * std::log10(1 + std::max(level, 0.0)), 2);
  info.sigma = py_round(sig, 4);
  if (level < 6 * sig) return {seeds, info};
  const double thr = std::max(level / 2, 4 * sig);
  std::vector<uint8_t> unk(tf.n), on(tf.n);
  for (size_t k = 0; k < tf.n; ++k) unk[k] = med[k] > 4 * level + 4 * sig;
  // widened by the smoothing length (python binary_dilation(unk, ones(SMOOTH)))
  unk = ndi::binary_dilation(unk, ndi::ones1(kSmooth));
  for (size_t k = 0; k < tf.n; ++k) on[k] = r[k] > thr && !unk[k];
  const int64_t max_unk = static_cast<int64_t>(kMaxBridgeS / tf.dt);
  const int64_t max_dip = std::max<int64_t>(1, static_cast<int64_t>(2 * min_sym_s / tf.dt));
  std::vector<std::pair<int64_t, int64_t>> runs;
  for (int64_t k = 0; k < static_cast<int64_t>(tf.n); ++k) {
    if (!on[k]) continue;
    if (!runs.empty()) {
      const int64_t b = runs.back().second, len = k - b;
      int64_t known = 0;
      for (int64_t g = b; g < k; ++g) known += !unk[g];
      if (len <= max_unk && known <= max_dip) {
        runs.back().second = k + 1;
        continue;
      }
    }
    runs.emplace_back(k, k + 1);
  }
  const int64_t N = static_cast<int64_t>(tf.n);
  std::vector<std::pair<int64_t, int64_t>> edged;
  for (auto [a, b] : runs) {
    while (b - a > 1 && rr[a] < thr && !unk[a]) ++a;
    while (b - a > 1 && rr[b - 1] < thr && !unk[b - 1]) --b;
    while (a > 0 && rr[a - 1] > thr && !unk[a - 1] && !(!edged.empty() && a - 1 < edged.back().second)) --a;
    while (b < N && rr[b] > thr && !unk[b]) ++b;
    edged.emplace_back(a, b);
  }
  const auto touch = tf.touched(seeds);
  Runs all;
  size_t n_mine = 0;
  for (const auto& [a, b] : edged) {
    bool t = false;
    for (int64_t k = a; k < b; ++k) t = t || touch[k];
    if (t) all.emplace_back(static_cast<double>(a) * tf.dt, static_cast<double>(b) * tf.dt), ++n_mine;
  }
  for (auto [t0, t1] : seeds) {
    int64_t a = static_cast<int64_t>(std::floor(t0 / tf.dt));
    int64_t b = std::min<int64_t>(static_cast<int64_t>(std::ceil(t1 / tf.dt)), N);
    if (b > a) {
      bool any_unk = false, any_on = false;
      for (int64_t k = a; k < b; ++k) any_unk = any_unk || unk[k], any_on = any_on || on[k];
      if (!any_unk && !any_on && np::mean(rr.data() + a, static_cast<size_t>(b - a)) < thr / 2) continue;
    }
    auto off = [&](int64_t k) { return r[k] < thr / 2 && !unk[k]; };
    while (b - a > 1 && off(a)) ++a, t0 = static_cast<double>(a) * tf.dt;
    while (b - a > 1 && off(b - 1)) --b, t1 = static_cast<double>(b) * tf.dt;
    all.emplace_back(t0, t1);
  }
  std::sort(all.begin(), all.end());
  Runs out;
  for (const auto& [a, b] : all) {
    if (!out.empty() && a <= out.back().second + std::max(tf.dt, min_sym_s))
      out.back().second = std::max(out.back().second, b);
    else
      out.emplace_back(a, b);
  }
  info.used = true;
  info.thr = py_round(thr, 4);
  info.n_energy = n_mine;
  return {out, info};
}

}  // namespace si::bursts
