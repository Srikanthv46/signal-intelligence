// FMCW radio altimeters (research, signal gaps S8; python: receivers/radalt/detector.py): STFT ridge -> straight
// crossings -> sweep rate, shape, turnarounds, bandwidth, period.
#include "receivers/radalt/detector.hpp"

#include <algorithm>
#include <cmath>
#include <optional>

#include "core/fft.hpp"
#include "core/npcompat.hpp"
#include "core/parallel.hpp"

namespace si::rx::radalt {

namespace {

constexpr double kBandLo = 4200e6, kBandHi = 4400e6;
constexpr int kNfft = 128, kHop = 64;
constexpr double kRidgeDb = 10.0;
constexpr double kMaxStepBins = 4;
constexpr size_t kMinRun = 5;
constexpr size_t kMaxGap = 8;
constexpr double kMergeS = 20e-6;
constexpr double kMinFill = 0.5;
constexpr double kMinSpanHz = 15e6;
constexpr double kMinR2 = 0.99;
constexpr double kRateLo = 2e9, kRateHi = 1500e9;
constexpr double kTablePeriods[] = {20e-3, 6.667e-3, 7e-3, 1e-3};

struct Ridge {
  dvec t, f;
  std::vector<bool> on;
};

// python ridge(x, fs)
std::optional<Ridge> ridge(const cfvec& x, double fs) {
  if (x.size() < static_cast<size_t>(kNfft)) return std::nullopt;
  const size_t nf = (x.size() - kNfft) / kHop + 1;
  if (nf < kMinRun) return std::nullopt;
  const dvec win = np::hanning(kNfft);
  Ridge r{dvec(nf), dvec(nf), std::vector<bool>(nf)};
  std::vector<char> on(nf);
  const size_t nblk = (nf + 1023) / 1024;
  parallel_for(nblk, [&](size_t blk) {              // rows are independent: blocks of 1024 frames per thread
    cdvec fr(kNfft), F(kNfft);
    dvec P(kNfft), tmp(kNfft);
    for (size_t row = blk * 1024; row < std::min(nf, (blk + 1) * 1024); ++row) {
      for (int i = 0; i < kNfft; ++i) {
        const cf v = x[row * kHop + static_cast<size_t>(i)];
        fr[i] = cd(static_cast<double>(v.real()) * win[i], static_cast<double>(v.imag()) * win[i]);
      }
      fft::c2c(fr.data(), F.data(), kNfft, true);
      for (int i = 0; i < kNfft; ++i) {
        const double a = np::cabs(F[i]);
        P[i] = a * a;
      }
      const int k = static_cast<int>(np::argmax(P));
      const double pk = P[k];
      tmp = P;
      on[row] = 10 * std::log10(pk / (np::median_inplace(tmp.data(), kNfft) + 1e-30) + 1e-30) >= kRidgeDb;
      const double a = std::log(P[(k - 1 + kNfft) % kNfft] + 1e-30), b = std::log(pk + 1e-30),
                   c = std::log(P[(k + 1) % kNfft] + 1e-30);
      const double den = a - 2 * b + c;
      const double frac = den < 0 ? 0.5 * (a - c) / den : 0.0;
      const double kk = k < kNfft / 2 ? k : k - kNfft;
      r.t[row] = (static_cast<double>(row) * kHop + kNfft / 2.0) / fs;
      r.f[row] = (kk + frac) * (fs / kNfft);
    }
  });
  for (size_t i = 0; i < nf; ++i) r.on[i] = on[i] != 0;
  return r;
}

struct Crossing {
  double t_s, rate, span, r2, f_lo, f_hi;
  std::vector<size_t> frames;
};

// python _fit(t, f, r)
std::optional<Crossing> fit(const Ridge& rd, const std::vector<size_t>& r) {
  const size_t n = r.size();
  dvec tt(n), ff(n);
  for (size_t i = 0; i < n; ++i) {
    tt[i] = rd.t[r[i]];
    ff[i] = rd.f[r[i]];
  }
  const double tm = np::mean(tt), fm = np::mean(ff);
  dvec a(n), b(n);
  for (size_t i = 0; i < n; ++i) {
    const double d = tt[i] - tm;
    a[i] = d * d;
    b[i] = d * (ff[i] - fm);
  }
  const double sxx = np::sum(a);
  if (sxx <= 0) return std::nullopt;
  const double slope = np::sum(b) / sxx;
  if (slope == 0) return std::nullopt;
  for (size_t i = 0; i < n; ++i) {
    const double res = ff[i] - (fm + slope * (tt[i] - tm));
    a[i] = res * res;
    const double d = ff[i] - fm;
    b[i] = d * d;
  }
  const double syy = np::sum(b);
  const double r2 = syy > 0 ? 1 - np::sum(a) / syy : 0.0;
  return Crossing{tm - fm / slope, slope, std::fabs(slope) * (tt[n - 1] - tt[0]), r2,
                  *std::min_element(ff.begin(), ff.end()), *std::max_element(ff.begin(), ff.end()), r};
}

// python crossings(t, f, on, fs)
std::vector<Crossing> crossings(const Ridge& rd, double fs) {
  const double step = kMaxStepBins * fs / kNfft;
  std::vector<std::vector<size_t>> open, done;
  for (size_t i = 0; i < rd.on.size(); ++i) {
    if (!rd.on[i]) continue;
    std::vector<std::vector<size_t>> keep;
    for (auto& r : open) (i - r.back() <= kMaxGap ? keep : done).push_back(std::move(r));
    open = std::move(keep);
    bool joined = false;
    for (auto& r : open) {
      if (std::fabs(rd.f[i] - rd.f[r.back()]) <= step * static_cast<double>(i - r.back())) {
        r.push_back(i);
        joined = true;
        break;
      }
    }
    if (!joined) open.push_back({i});
  }
  for (auto& r : open) done.push_back(std::move(r));
  std::vector<std::vector<size_t>> pieces;
  for (const auto& r : done) {
    if (r.size() < kMinRun || static_cast<double>(r.size()) < kMinFill * static_cast<double>(r.back() - r[0] + 1))
      continue;
    const auto c = fit(rd, r);
    if (c && c->r2 >= kMinR2) {
      pieces.push_back(r);
      continue;
    }
    dvec fr(r.size());
    for (size_t i = 0; i < r.size(); ++i) fr[i] = rd.f[r[i]];
    for (size_t k : {np::argmax(fr), np::argmin(fr)}) {
      if (kMinRun <= k + 1 && r.size() - k >= kMinRun) {
        pieces.emplace_back(r.begin(), r.begin() + static_cast<ptrdiff_t>(k) + 1);
        pieces.emplace_back(r.begin() + static_cast<ptrdiff_t>(k), r.end());
        break;
      }
    }
  }
  std::stable_sort(pieces.begin(), pieces.end(), [](const auto& a, const auto& b) { return a[0] < b[0]; });
  std::vector<Crossing> out;
  for (const auto& r : pieces) {
    const auto c = fit(rd, r);
    if (!c || c->r2 < kMinR2) continue;
    if (!out.empty() && (out.back().rate > 0) == (c->rate > 0) && std::fabs(out.back().t_s - c->t_s) < kMergeS) {
      std::vector<size_t> u = out.back().frames;
      u.insert(u.end(), r.begin(), r.end());
      std::sort(u.begin(), u.end());
      const auto m = fit(rd, u);
      if (m && m->r2 >= kMinR2) {
        out.back() = *m;
        continue;
      }
    }
    out.push_back(*c);
  }
  std::vector<Crossing> kept;
  for (auto& c : out)
    if (c.span >= kMinSpanHz) kept.push_back(std::move(c));
  return kept;
}

}  // namespace

std::vector<Emission> detect(const cfvec& x, double fs, double fc) {
  if (!(fc - fs / 2 < kBandHi && fc + fs / 2 > kBandLo) || fs < 20e6) return {};
  const auto rd = ridge(x, fs);
  if (!rd) return {};
  const auto cr = crossings(*rd, fs);
  if (cr.empty()) return {};
  dvec rates;
  for (const auto& c : cr) rates.push_back(std::fabs(c.rate));
  const double rate = np::median(rates);
  if (!(kRateLo <= rate && rate <= kRateHi)) return {};
  std::vector<int> signs;
  for (const auto& c : cr) signs.push_back(c.rate > 0 ? 1 : -1);
  double r2min = cr[0].r2;
  for (const auto& c : cr) r2min = std::min(r2min, c.r2);
  Json feats = Json::object();
  feats["n_crossings"] = static_cast<int64_t>(cr.size());
  feats["rate_mhz_per_ms"] = py_round(rate / 1e9, 3);
  feats["r2_min"] = py_round(r2min, 4);
  std::optional<double> period;
  for (int sg : {1, -1}) {
    dvec ts;
    for (size_t i = 0; i < cr.size(); ++i)
      if (signs[i] == sg) ts.push_back(cr[i].t_s);
    if (ts.size() >= 2) {   // spacings are whole periods when a crossing is missed
      const dvec d = np::diff(ts);
      const double dmin = *std::min_element(d.begin(), d.end());
      dvec q(d.size());
      for (size_t i = 0; i < d.size(); ++i) q[i] = d[i] / std::max(std::nearbyint(d[i] / dmin), 1.0);
      period = np::median(q);
      break;
    }
  }
  dvec tops, bots;
  for (size_t i = 0; i + 1 < cr.size(); ++i) {
    const auto &a = cr[i], &b = cr[i + 1];
    if (a.rate > 0 && 0 > b.rate)
      tops.push_back(fc + rate * (b.t_s - a.t_s) / 2);
    else if (a.rate < 0 && 0 < b.rate)
      bots.push_back(fc - rate * (b.t_s - a.t_s) / 2);
  }
  const bool both = std::find(signs.begin(), signs.end(), 1) != signs.end() &&
                    std::find(signs.begin(), signs.end(), -1) != signs.end();
  if (both)
    feats["shape"] = "triangle";
  else if (cr.size() >= 2)
    feats["shape"] = signs[0] > 0 ? "sawtooth up" : "sawtooth down";
  else
    feats["shape"] = "single crossing";
  std::optional<double> bw;
  if (!tops.empty()) feats["f_top_mhz"] = py_round(np::median(tops) / 1e6, 2);
  if (!bots.empty()) feats["f_bottom_mhz"] = py_round(np::median(bots) / 1e6, 2);
  if (!tops.empty() && !bots.empty()) {
    bw = np::median(tops) - np::median(bots);
    feats["sweep_bw_mhz"] = py_round(*bw / 1e6, 1);
    if (both) period = 2 * *bw / rate;
  }
  if (period) feats["period_ms"] = py_round(*period * 1e3, 3);
  double conf = 0.7;
  bool table = false;
  if (period)
    for (double p : kTablePeriods) table = table || std::fabs(*period - p) <= 0.05 * p;
  if ((bw && 95e6 <= *bw && *bw <= 190e6) || table) conf = 0.8;
  double lo = cr[0].f_lo, hi = cr[0].f_hi;
  for (const auto& c : cr) {
    lo = std::min(lo, c.f_lo);
    hi = std::max(hi, c.f_hi);
  }
  const double seen_lo = lo + fc, seen_hi = hi + fc;
  feats["seen_mhz"] = Json::array({py_round(seen_lo / 1e6, 2), py_round(seen_hi / 1e6, 2)});
  Emission e;
  e.family = "radalt";
  e.label = "radio_altimeter";
  e.start_us = rd->t.front() * 1e6;
  e.end_us = rd->t.back() * 1e6;
  e.f_lo_mhz = seen_lo / 1e6;
  e.f_hi_mhz = seen_hi / 1e6;
  e.center_mhz = (seen_lo + seen_hi) / 2e6;
  e.confidence = conf;
  e.features = feats;
  return {e};
}

}  // namespace si::rx::radalt
