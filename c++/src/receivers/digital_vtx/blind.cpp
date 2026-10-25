// Vendor-independent CP-OFDM link detector (python: receivers/digital_vtx/blind.py). Every step
// follows the python line by line, with numpy's arithmetic (FMA complex multiply, Smith division,
// pairwise sums, sequential cumsum / bincount, half-even rounding).
#include "receivers/digital_vtx/blind.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <tuple>

#include "core/fft.hpp"
#include "core/ofdm15k.hpp"
#include "core/spectral.hpp"
#include "core/ndimage.hpp"
#include "core/npcompat.hpp"
#include "core/signal.hpp"
#include "receivers/digital_vtx/common.hpp"

namespace si::rx::digital_vtx::blind {

namespace {
using detail::expj;
using detail::kPi;
using detail::trunc_i;

constexpr double WIDTHS[] = {10e6, 20e6};
constexpr double AGC_S = 20e-6;
constexpr double TU_SCAN0 = 8e-6, TU_SCAN1 = 120e-6;
constexpr double K_ACF = 7.0;
constexpr double LINE_Z = 10.0;
constexpr double WIN_Z = 8.0;
constexpr int64_t WIN_N = 3;
constexpr size_t MAX_CLUSTERS = 24;
constexpr double CP_RATIOS[] = {1.0 / 32, 1.0 / 16, 1.0 / 14, 1.0 / 12, 1.0 / 8, 1.0 / 6, 1.0 / 4};
constexpr double SEG_WIN_S = 2e-3;
constexpr double SEG_ZW = 4.0;
constexpr int64_t MIN_SYMS = 40, MIN_LONGEST = 10;
constexpr size_t ELIDE_CPLX = 256 * 1024 / 16;  // numpy reuses a temporary this large (a * conj(b) -> conj(b) * a)

inline int64_t npround_int(double x) { return static_cast<int64_t>(std::nearbyint(x)); }
inline int64_t floordiv(int64_t a, int64_t b) {
  const int64_t q = a / b;
  return (a % b != 0 && ((a < 0) != (b < 0))) ? q - 1 : q;
}

dvec abs2(const cdvec& y) {
  dvec a(y.size());
  for (size_t i = 0; i < y.size(); ++i) a[i] = np::cabs2(y[i]);
  return a;
}

// python _cp_metric(y, N, w) -> (|movsum(p, w)| / movsum(e, w), movsum(p, w))
struct Cp {
  dvec m;
  cdvec cs;
};
Cp cp_metric(const cdvec& y, size_t N, size_t w) {
  Cp r;
  if (y.size() <= N) return r;
  const size_t n = y.size() - N;
  const bool swap = n >= ELIDE_CPLX;
  cdvec c(n + 1);
  dvec ce(n + 1);
  c[0] = cd(0.0, 0.0);
  ce[0] = 0.0;
  for (size_t i = 0; i < n; ++i) {
    c[i + 1] = swap ? np::cmul(std::conj(y[i + N]), y[i]) : np::cmul(y[i], std::conj(y[i + N]));
    ce[i + 1] = (np::cabs2(y[i]) + np::cabs2(y[i + N])) / 2;
  }
  for (size_t i = 0; i < n; ++i) {
    c[i + 1] = c[i] + c[i + 1];
    ce[i + 1] = ce[i] + ce[i + 1];
  }
  if (n + 1 <= w) return r;
  const size_t nv = n + 1 - w;
  r.m.resize(nv);
  r.cs.resize(nv);
  for (size_t j = 0; j < nv; ++j) {
    const cd s(c[j + w].real() - c[j].real(), c[j + w].imag() - c[j].imag());
    r.cs[j] = s;
    r.m[j] = np::cabs(s) / ((ce[j + w] - ce[j]) + 1e-20);
  }
  return r;
}

cdvec agc(const cdvec& y, double ofs) {
  const dvec env = ndi::uniform_filter1d(abs2(y), static_cast<size_t>(std::max<int64_t>(4, trunc_i(AGC_S * ofs))));
  cdvec out(y.size());
  for (size_t i = 0; i < y.size(); ++i) out[i] = np::cdiv(y[i], cd(std::sqrt(env[i] + 1e-30), 0.0));
  return out;
}

// python _acf(y, Lmax): energy-normalised sum_n y[n+L] conj(y[n]), L = 0..Lmax
cdvec acf(const cdvec& y, size_t Lmax) {
  const size_t n = y.size();
  const size_t nf = size_t(1) << trunc_i(std::ceil(std::log2(static_cast<double>(n + Lmax + 1))));
  cdvec F = fft::fft(y.data(), n, nf);
  for (auto& v : F) {
    const double a = np::cabs(v);
    v = cd(a * a, 0.0);
  }
  const cdvec r = fft::ifft(F);
  dvec cs(n + 1, 0.0);
  for (size_t i = 0; i < n; ++i) cs[i + 1] = cs[i] + np::cabs2(y[i]);
  cdvec out(Lmax + 1);
  for (size_t L = 0; L <= Lmax; ++L) {
    const double en = 0.5 * ((cs[n - L] - cs[0]) + (cs[n] - cs[L]));
    out[L] = np::cdiv(std::conj(r[L]), cd(en + 1e-30, 0.0));
  }
  return out;
}

// python _symrate(m, N, ofs) -> (symbol period in samples, line z, period is a numpy float)
std::tuple<double, double, bool> symrate(const dvec& m, int64_t N, double ofs, double seg_s = 4e-3, double lo = 1.02,
                                  double hi = 1.30) {
  const size_t L = std::min<size_t>(static_cast<size_t>(trunc_i(seg_s * ofs)), m.size());
  const size_t nf = size_t(1) << trunc_i(std::ceil(std::log2(static_cast<double>(L))));
  const dvec w = np::hanning(L);
  dvec acc(nf / 2 + 1, 0.0);
  dvec s(L);
  for (size_t i = 0; i + L <= m.size(); i += L) {
    const double mu = np::mean(m.data() + i, L);
    for (size_t j = 0; j < L; ++j) s[j] = (m[i + j] - mu) * w[j];
    const cdvec F = fft::rfft(s.data(), L, nf);
    for (size_t j = 0; j < acc.size(); ++j) {
      const double a = np::cabs(F[j]);
      acc[j] += a * a;
    }
  }
  const double Nd = static_cast<double>(N), nfd = static_cast<double>(nf);
  dvec sel(acc.size()), wide;
  for (size_t j = 0; j < acc.size(); ++j) {
    const double a = static_cast<double>(j) / nfd;
    sel[j] = (a >= 1 / (hi * Nd) && a <= 1 / (lo * Nd)) ? acc[j] : -1.0;
    if (a >= 1 / (2.5 * Nd) && a <= 1 / (0.8 * Nd)) wide.push_back(acc[j]);
  }
  const size_t i = np::argmax(sel);
  const double z = acc[i] / (np::median(wide) + 1e-30);
  double d = 0.0;
  bool np_ = false;
  if (i > 0 && i < acc.size() - 1) {
    const double l = acc[i - 1], c = acc[i], r = acc[i + 1];
    const double den = (l - 2 * c) + r;
    np_ = den != 0;
    d = np_ ? (0.5 * (l - r)) / den : 0.0;
  }
  return {nfd / (static_cast<double>(i) + d), z, np_};
}

// python _tu_refine(y, ofs, N)
double tu_refine(const cdvec& y, double ofs, int64_t N) {
  double r[5];
  const int64_t n = static_cast<int64_t>(y.size());
  for (int64_t j = 0; j < 5; ++j) {
    const int64_t L = N - 2 + j;
    cd s(0.0, 0.0);
    for (int64_t i = 0; i + L < n; ++i) s += np::cmul(std::conj(y[static_cast<size_t>(i)]), y[static_cast<size_t>(i + L)]);
    r[j] = np::cabs(s);
  }
  int64_t j = 1;
  for (int64_t k = 2; k <= 3; ++k)
    if (r[k] > r[j]) j = k;
  const double den = (r[j - 1] - 2 * r[j]) + r[j + 1];
  const double d = den < 0 ? (0.5 * (r[j - 1] - r[j + 1])) / den : 0.0;
  return (static_cast<double>(N - 2 + j) + std::clamp(d, -0.5, 0.5)) / ofs;
}

// python window_phases(m, P, Wn, H) -> (phases, z)
// np.mod(x, P) for x >= 0, P > 0, bit-identical to fmod (the true remainder x - qP is representable, so fma
// returns it exactly) but without fmod's slow long division (PLAN P1)
inline double mod_pos(double x, double P) {
  double q = std::floor(x / P);
  double r = std::fma(-q, P, x);
  if (r < 0) r = std::fma(-(q - 1), P, x);
  else if (r >= P) r = std::fma(-(q + 1), P, x);
  return r;
}

std::pair<std::vector<int64_t>, dvec> window_phases(const dvec& m, double P, int64_t Wn, int64_t H) {
  const int64_t nW = std::max<int64_t>(1, floordiv(static_cast<int64_t>(m.size()) - Wn, H) + 1);
  const size_t nb = static_cast<size_t>(std::ceil(P));
  std::vector<int64_t> ph(static_cast<size_t>(nW));
  dvec zw(static_cast<size_t>(nW));
  dvec s(nb), f(nb), d(nb);
  std::vector<int64_t> cnt(nb);
  for (int64_t w = 0; w < nW; ++w) {
    const size_t a = static_cast<size_t>(w * H);
    const size_t b = std::min(m.size(), static_cast<size_t>(w * H + Wn));
    std::fill(s.begin(), s.end(), 0.0);
    std::fill(cnt.begin(), cnt.end(), 0);
    for (size_t i = a; i < b; ++i) {
      const size_t k = static_cast<size_t>(std::floor(mod_pos(static_cast<double>(i), P)));
      s[k] += m[i];
      ++cnt[k];
    }
    dvec g(nb);
    for (size_t k = 0; k < nb; ++k) g[k] = s[k] / static_cast<double>(std::max<int64_t>(cnt[k], 1));
    for (size_t k = 0; k < nb; ++k) f[k] = (g[k] + g[(k + nb - 1) % nb]) + g[(k + 1) % nb];
    const double med = np::median(f);
    const size_t k = np::argmax(f);
    for (size_t i = 0; i < nb; ++i) d[i] = std::fabs(f[i] - med);
    ph[static_cast<size_t>(w)] = static_cast<int64_t>(k);
    zw[static_cast<size_t>(w)] = (f[k] - med) / (1.4826 * np::median(d) + 1e-12);
  }
  return {ph, zw};
}

// python refine_period(m, P, ofs)
double refine_period(const dvec& m, double P, double ofs, double win_s = 1e-3, int iters = 4) {
  const int64_t Wn = trunc_i(win_s * ofs);
  for (int it = 0; it < iters; ++it) {
    const auto [ph, zw] = window_phases(m, P, Wn, Wn);
    dvec d;
    for (size_t i = 0; i + 1 < ph.size(); ++i) {
      if (!(zw[i] > 6 && zw[i + 1] > 6)) continue;
      const double v = np::mod(static_cast<double>(ph[i + 1] - ph[i]) + P / 2, P) - P / 2;
      if (std::fabs(v) < P / 4) d.push_back(v);
    }
    if (d.size() < 3) break;
    const double st = np::median(d);
    if (std::fabs(st) < 0.5) break;
    P = P + st / (static_cast<double>(Wn) / P);
  }
  return P;
}

double sum_range(const dvec& v, size_t a, size_t b) { return np::pairwise_sum(v.data() + a, b - a); }

}  // namespace

std::pair<double, int64_t> rate_for(double tu, double B, double cp) {
  const int64_t k = static_cast<int64_t>(std::ceil(tu * std::max({2.5e6, 1.25 * B, cp > 0 ? 16 / cp : 0.0})));
  return {static_cast<double>(k) / tu, k};
}

std::vector<Hit> scan(const Slicer& sl, double usable, const SkipFn& skip) {
  std::vector<Hit> hits;
  dvec widths;
  for (double w : WIDTHS)
    if (w <= 2 * usable) widths.push_back(w);
  if (widths.empty()) widths.push_back(2 * usable);  // a narrow capture: one slice
  for (double W : widths) {
    for (double off : np::arange(-usable + W / 2, usable - W / 2 + 1, W / 2)) {
      if (skip) {
        bool all = true;
        for (double f : np::arange(off - W / 2, off + W / 2 + 1, 1e6))
          if (!skip(f)) {
            all = false;
            break;
          }
        if (all) continue;
      }
      double ofs = 0;
      const cdvec y = agc(sl(off, W / 2, 1.25 * W, &ofs), ofs);
      const int64_t Lmin = trunc_i(TU_SCAN0 * ofs), Lmax = trunc_i(TU_SCAN1 * ofs);
      if (Lmax >= static_cast<int64_t>(y.size() / 2)) continue;
      const cdvec r = acf(y, static_cast<size_t>(Lmax));
      dvec a(static_cast<size_t>(Lmax - Lmin + 1));
      for (size_t i = 0; i < a.size(); ++i) a[i] = np::cabs(r[static_cast<size_t>(Lmin) + i]);
      const dvec base = ndi::median_filter(a, static_cast<size_t>(std::max<int64_t>(7, trunc_i(1.5e-6 * ofs)) | 1));
      dvec d(a.size());
      for (size_t i = 0; i < a.size(); ++i) d[i] = a[i] - base[i];
      const double md = np::median(d);
      dvec ad(d.size());
      for (size_t i = 0; i < d.size(); ++i) ad[i] = std::fabs(d[i] - md);
      const double sd = 1.4826 * np::median(ad) + 1e-12;
      sig::PeakOpts po;
      po.height_min = K_ACF * sd;
      po.distance = static_cast<double>(std::max<int64_t>(3, trunc_i(0.5e-6 * ofs)));
      for (int64_t p : sig::find_peaks(d, po).peaks)
        hits.push_back({W, off, p + Lmin, ofs, static_cast<double>(p + Lmin) / ofs, d[static_cast<size_t>(p)] / sd});
    }
  }
  return hits;
}

std::optional<Confirmed> confirm(const Slicer& sl, double off, double W, int64_t N0, double ofs0) {
  cdvec yr = sl(off, W / 2, ofs0, &ofs0);
  const cdvec y = agc(yr, ofs0);
  double P0 = 0, zl = 0;
  int64_t NN = 0;
  bool have = false, p0_np = false;
  for (int64_t n : {N0 - 1, N0, N0 + 1}) {
    const Cp c = cp_metric(y, static_cast<size_t>(n), static_cast<size_t>(std::max<int64_t>(3, n / 16)));
    const auto [P, z, pn] = symrate(c.m, n, ofs0);
    if (!have || z > zl) P0 = P, zl = z, NN = n, p0_np = pn, have = true;
  }
  const double tu = tu_refine(yr, ofs0, NN);
  auto [ofs, N] = rate_for(tu, W);
  yr = sl(off, W / 2, ofs, &ofs);
  const Cp c = cp_metric(yr, static_cast<size_t>(N), static_cast<size_t>(std::max<int64_t>(3, N / 16)));
  const int64_t Wn = trunc_i(1e-3 * ofs);
  int64_t nwin = 0;
  double Pw = 0;
  bool pw_np = false;
  dvec cands{P0 / ofs0 * ofs};
  for (double r : CP_RATIOS) cands.push_back(static_cast<double>(N) * (1 + r));
  const double Nd = static_cast<double>(N);
  for (size_t ci = 0; ci < cands.size(); ++ci) {
    const double Pc = cands[ci];
    if (!(1.02 * Nd <= Pc && Pc <= 1.30 * Nd)) continue;
    const double Pr = refine_period(c.m, Pc, ofs);
    const dvec zw = window_phases(c.m, Pr, Wn, Wn).second;
    int64_t k = 0;
    for (double z : zw) k += z > WIN_Z;
    if (k > nwin) nwin = k, Pw = Pr, pw_np = ci == 0 && p0_np;
  }
  if (zl >= LINE_Z) return Confirmed{tu, P0 / ofs0, zl, nwin, p0_np};
  if (nwin >= WIN_N) return Confirmed{tu, Pw / ofs, zl, nwin, pw_np};
  return std::nullopt;
}

Segmented segment(const Slicer& sl, double off, double B, double tu, double ts) {
  Segmented out;
  auto [ofs, N] = rate_for(tu, B, ts - tu);
  const cdvec y = sl(off, B / 2, ofs, &ofs);
  double P = ts * ofs;
  const int64_t cp = std::max<int64_t>(3, npround_int(P - static_cast<double>(N)));
  const Cp c = cp_metric(y, static_cast<size_t>(N), static_cast<size_t>(std::max<int64_t>(3, trunc_i(0.8 * static_cast<double>(cp)))));
  const dvec& m = c.m;
  P = refine_period(m, P, ofs);
  const int64_t Wn = trunc_i(SEG_WIN_S * ofs);
  const int64_t H = Wn / 2;
  const auto [ph, zw] = window_phases(m, P, Wn, H);
  std::vector<int64_t> act;
  for (size_t i = 0; i < zw.size(); ++i)
    if (zw[i] > SEG_ZW) act.push_back(static_cast<int64_t>(i));
  out.ts = P / ofs;
  if (act.empty()) return out;
  dvec sp;
  for (size_t i = 0; i < zw.size(); ++i)
    if (zw[i] > 6) sp.push_back(static_cast<double>(ph[i]));
  out.ngrid = 1;
  if (sp.size() > 1)
    for (size_t i = 0; i + 1 < sp.size(); ++i)
      out.ngrid += std::fabs(np::mod((sp[i + 1] - sp[i]) + P / 2, P) - P / 2) > static_cast<double>(cp);
  const int64_t nm = static_cast<int64_t>(m.size());
  auto smp = [&](double a) {
    const int64_t r = npround_int(a);
    const int64_t lo = std::max<int64_t>(r - 2, 0), hi = std::min<int64_t>(r + 3, nm);
    double v = m[static_cast<size_t>(lo)];
    for (int64_t i = lo + 1; i < hi; ++i) v = std::max(v, m[static_cast<size_t>(i)]);
    return v;
  };
  dvec tk, on, of, q;
  for (int64_t bk = 0; bk < nm / H; ++bk) {
    std::vector<int64_t> cand;
    for (int64_t w : {bk - 1, bk})
      if (w >= 0 && w < static_cast<int64_t>(zw.size()) && zw[static_cast<size_t>(w)] > SEG_ZW) cand.push_back(w);
    if (cand.empty()) {
      int64_t j = act[0];
      for (int64_t a : act)
        if (std::llabs(a - bk) < std::llabs(j - bk)) j = a;
      if (std::llabs(j - bk) > 2) continue;
      cand.push_back(j);
    }
    int64_t w = cand[0];
    for (int64_t i : cand)
      if (zw[static_cast<size_t>(i)] > zw[static_cast<size_t>(w)] ||
          (zw[static_cast<size_t>(i)] == zw[static_cast<size_t>(w)] && i < w))
        w = i;
    const double p0 = static_cast<double>(ph[static_cast<size_t>(w)]);
    int64_t k = static_cast<int64_t>(std::ceil((static_cast<double>(bk * H) - p0) / P));
    while (p0 + static_cast<double>(k) * P < static_cast<double>((bk + 1) * H) &&
           p0 + static_cast<double>(k) * P + P + 3 < static_cast<double>(nm)) {
      const double a = p0 + static_cast<double>(k) * P;
      tk.push_back(a);
      on.push_back(smp(a));
      of.push_back(smp(a + P / 2));
      q.push_back(smp(a + P / 4));
      ++k;
    }
  }
  const size_t K = tk.size();
  dvec nref(K);
  for (size_t i = 0; i < K; ++i) nref[i] = of[i] - q[i];
  const double mr0 = np::median(nref);
  dvec ad(K);
  for (size_t i = 0; i < K; ++i) ad[i] = std::fabs(nref[i] - mr0);
  const double sg = 1.4826 * np::median(ad) + 1e-9;
  dvec zr(K);
  for (size_t i = 0; i < K; ++i) zr[i] = (on[i] - of[i]) / sg;
  const dvec zs = ndi::uniform_filter1d(zr, 4);
  std::vector<uint8_t> mask(K);
  for (size_t i = 0; i < K; ++i) mask[i] = zs[i] * 2 > 4.0;
  mask = ndi::binary_closing(mask, ndi::ones1(3));
  const auto lab = ndi::label(mask);
  std::vector<Run> runs;
  for (const auto& [a0, b0] : ndi::find_objects(lab.lab, lab.n)) {
    size_t a = a0, b = b0;
    while (a < b && zr[a] < 1.0) ++a;
    while (b > a && zr[b - 1] < 1.0) --b;
    if (b - a >= 3 && np::mean(zr.data() + a, b - a) * std::sqrt(static_cast<double>(b - a)) >= 6)
      runs.push_back({tk[a] / ofs, (tk[b - 1] + P) / ofs, static_cast<int64_t>(b - a)});
  }
  for (const auto& r : runs) {
    if (!out.runs.empty() && r.t0 - out.runs.back().t1 <= 3.5 * P / ofs)
      out.runs.back() = {out.runs.back().t0, r.t1, out.runs.back().n + r.n};
    else
      out.runs.push_back(r);
  }
  cdvec g;
  for (size_t i = 0; i < K; ++i) {
    const int64_t t = trunc_i(tk[i]);
    if (zr[i] > 3 && t < static_cast<int64_t>(c.cs.size())) g.push_back(c.cs[static_cast<size_t>(t)]);
  }
  if (!g.empty()) {
    const cd s = np::pairwise_sum(g.data(), g.size());
    out.frac = -std::atan2(s.imag(), s.real()) / (2 * kPi) * (ofs / static_cast<double>(N));
  }
  out.tk.resize(K);
  for (size_t i = 0; i < K; ++i) out.tk[i] = tk[i] / ofs;
  out.zr = std::move(zr);
  return out;
}

namespace {
std::pair<double, dvec> fold_contrast(const dvec& tk, const dvec& e, double T, double ts) {
  const int64_t nb = std::max<int64_t>(8, npround_int(T / ts));
  dvec s(static_cast<size_t>(nb), 0.0);
  std::vector<int64_t> cnt(static_cast<size_t>(nb), 0);
  for (size_t i = 0; i < tk.size(); ++i) {
    const int64_t b = std::min<int64_t>(trunc_i(np::mod(tk[i], T) / T * static_cast<double>(nb)), nb - 1);
    s[static_cast<size_t>(b)] += e[i];
    ++cnt[static_cast<size_t>(b)];
  }
  const double sc = std::sqrt(std::max(tk.back() / T, 1.0));
  dvec zt(static_cast<size_t>(nb));
  for (size_t b = 0; b < zt.size(); ++b) zt[b] = s[b] / static_cast<double>(std::max<int64_t>(cnt[b], 1)) * sc;
  dvec so = zt;
  std::sort(so.begin(), so.end());
  const size_t top = static_cast<size_t>(std::max<int64_t>(nb / 4, 1));
  return {np::mean(so.data() + so.size() - top, top), zt};
}

dvec clip_e(const dvec& zr) {
  dvec e(zr.size());
  for (size_t i = 0; i < zr.size(); ++i) e[i] = std::min(std::max(zr[i], -3.0), 6.0);
  return e;
}
}  // namespace

std::pair<std::optional<double>, double> frame_period(const dvec& tk, const dvec& zr, double ts, double lo,
                                                      double hi) {
  if (tk.size() < 50) return {std::nullopt, 0.0};
  const dvec e = clip_e(zr);
  const size_t n = static_cast<size_t>(trunc_i(std::ceil(tk.back() / ts))) + 1;
  dvec a(n, 0.0);
  for (size_t i = 0; i < tk.size(); ++i) a[static_cast<size_t>(npround_int(tk[i] / ts))] = e[i];
  const double mu = np::mean(a);
  for (auto& v : a) v = v - mu;
  cdvec F = fft::rfft(a.data(), n, 2 * n);
  for (auto& v : F) {
    const double m = np::cabs(v);
    v = cd(m * m, 0.0);
  }
  dvec A = fft::irfft(F.data(), F.size(), 2 * n);
  A.resize(n);
  const double a0 = A[0] + 1e-30;
  for (auto& v : A) v = v / a0;
  A = ndi::uniform_filter1d(A, static_cast<size_t>(std::max<int64_t>(3, trunc_i(0.2e-3 / ts))));
  dvec sel(n);
  std::vector<uint8_t> ok(n);
  for (size_t L = 0; L < n; ++L) {
    const double t = static_cast<double>(L) * ts;
    ok[L] = t >= lo && t <= hi && L < n / 3;
    sel[L] = ok[L] ? A[L] : -1.0;
  }
  std::vector<int64_t> pk;
  for (int64_t p : sig::find_peaks(sel).peaks)
    if (ok[static_cast<size_t>(p)] && A[static_cast<size_t>(p)] > 0.1) pk.push_back(p);
  if (pk.empty()) return {std::nullopt, 0.0};
  double best = A[static_cast<size_t>(pk[0])];
  for (int64_t p : pk) best = std::max(best, A[static_cast<size_t>(p)]);
  int64_t p0 = -1;
  for (int64_t p : pk)
    if (A[static_cast<size_t>(p)] >= 0.7 * best && (p0 < 0 || p < p0)) p0 = p;
  const double T0 = static_cast<double>(p0) * ts;
  const dvec Ts = np::arange(0.94 * T0, 1.06 * T0, ts / 4);
  // fine search on the variance the fold explains (R^2, fixed bin count): the fold contrast scales with
  // sqrt(frames in the block), so it always grew towards the search's short edge (python frame_period; PLAN D7)
  const int64_t nb = std::max<int64_t>(8, npround_int(T0 / ts));
  const double emu = np::mean(e);
  dvec ec(e.size()), ec2(e.size());
  for (size_t i = 0; i < e.size(); ++i) {
    ec[i] = e[i] - emu;
    ec2[i] = ec[i] * ec[i];
  }
  const double tot = np::pairwise_sum(ec2.data(), ec2.size()) + 1e-30;
  size_t k = 0;
  double br = 0;
  dvec sm(static_cast<size_t>(nb)), q(static_cast<size_t>(nb));
  std::vector<int64_t> cnt(static_cast<size_t>(nb));
  for (size_t i = 0; i < Ts.size(); ++i) {
    const double T = Ts[i];
    std::fill(sm.begin(), sm.end(), 0.0);
    std::fill(cnt.begin(), cnt.end(), 0);
    for (size_t j = 0; j < tk.size(); ++j) {
      const int64_t b = std::min<int64_t>(trunc_i(np::mod(tk[j], T) / T * static_cast<double>(nb)), nb - 1);
      sm[static_cast<size_t>(b)] += ec[j];
      ++cnt[static_cast<size_t>(b)];
    }
    for (size_t b = 0; b < q.size(); ++b) q[b] = sm[b] * sm[b] / static_cast<double>(std::max<int64_t>(cnt[b], 1));
    const double r2 = np::pairwise_sum(q.data(), q.size()) / tot;
    if (i == 0 || r2 > br) br = r2, k = i;
  }
  return {Ts[k], fold_contrast(tk, e, Ts[k], ts).first};
}

std::vector<Run> frame_complete(const std::vector<Run>& runs, const dvec& tk, const dvec& zr, double ts, double T,
                                double min_con) {
  const dvec e = clip_e(zr);
  const auto [con, zt] = fold_contrast(tk, e, T, ts);
  if (con < min_con) return runs;
  const size_t nb = zt.size();
  const dvec sm = ndi::uniform_filter1d(zt, 5, ndi::Mode::Wrap);
  const size_t sh = np::argmin(sm);
  std::vector<uint8_t> rm(nb);
  for (size_t i = 0; i < nb; ++i) rm[i] = zt[(i + sh) % nb] > 3;
  const auto lab = ndi::label(rm);
  std::vector<std::pair<double, double>> regs;
  for (const auto& [a, b] : ndi::find_objects(lab.lab, lab.n))
    if (b - a >= 3)
      regs.emplace_back(static_cast<double>(a + sh) * T / static_cast<double>(nb),
                        static_cast<double>(b + sh) * T / static_cast<double>(nb));
  std::vector<Run> add;
  std::vector<std::pair<double, double>> present;
  const int64_t nf = trunc_i(tk.back() / T) + 1;
  for (int64_t f = -1; f < nf; ++f) {
    for (const auto& [a, b] : regs) {
      std::vector<size_t> sel;
      const double fT = static_cast<double>(f) * T;
      for (size_t i = 0; i < tk.size(); ++i)
        if (tk[i] >= fT + a && tk[i] < fT + b) sel.push_back(i);
      if (sel.size() < 3) continue;
      dvec seg(sel.size());
      for (size_t i = 0; i < sel.size(); ++i) seg[i] = e[sel[i]];
      if (np::mean(seg) * std::sqrt(static_cast<double>(seg.size())) < 4) continue;
      present.emplace_back(tk[sel.front()], tk[sel.back()] + ts);
      const dvec u = ndi::uniform_filter1d(seg, 4);
      std::vector<uint8_t> on(seg.size());
      for (size_t i = 0; i < seg.size(); ++i) on[i] = u[i] >= 1;
      std::vector<size_t> idx;
      for (size_t i = 0; i < on.size(); ++i)
        if (on[i]) idx.push_back(i);
      for (size_t i = 0; i + 1 < idx.size(); ++i)
        if (idx[i + 1] - idx[i] > 1 && idx[i + 1] - idx[i] <= 5)
          for (size_t j = idx[i]; j < idx[i + 1]; ++j) on[j] = 1;
      const auto lab2 = ndi::label(on);
      bool have = false;
      std::tuple<double, size_t, size_t> best{0.0, 0, 0};
      for (const auto& [a2, b2] : ndi::find_objects(lab2.lab, lab2.n)) {
        if (b2 - a2 < 3) continue;
        const std::tuple<double, size_t, size_t> key{sum_range(seg, a2, b2), a2, b2};
        if (!have || key > best) best = key, have = true;
      }
      if (!have) continue;
      size_t a2 = std::get<1>(best), b2 = std::get<2>(best);
      while (a2 < b2 && seg[a2] < 1) ++a2;
      while (b2 > a2 && seg[b2 - 1] < 1) --b2;
      if (b2 - a2 >= 3) add.push_back({tk[sel[a2]], tk[sel[b2 - 1]] + ts, static_cast<int64_t>(b2 - a2)});
    }
  }
  std::vector<Run> all = runs;
  all.insert(all.end(), add.begin(), add.end());
  std::sort(all.begin(), all.end(), [](const Run& x, const Run& y) {
    return std::tie(x.t0, x.t1, x.n) < std::tie(y.t0, y.t1, y.n);
  });
  std::vector<Run> out;
  for (const auto& r : all) {
    bool join = false;
    if (!out.empty()) {
      const double pe = out.back().t1;
      join = r.t0 <= pe + 1.5 * ts;
      if (!join && r.t0 - pe <= 6.5 * ts)
        for (const auto& [p0, p1] : present) join = join || (p0 <= pe && r.t0 <= p1);
    }
    if (join)
      out.back() = {out.back().t0, std::max(out.back().t1, r.t1), out.back().n + r.n};
    else
      out.push_back(r);
  }
  return out;
}

std::pair<dvec, dvec> band_cp(const Slicer& sl, double off, double W, double tu, double ts, const dvec& tsym) {
  auto [ofs, N] = rate_for(tu, 2 * W);
  const cdvec y = sl(off, W, ofs, &ofs);
  const int64_t cp = npround_int((ts - tu) * ofs);
  const int64_t M = std::max<int64_t>(8, trunc_i(0.9 * static_cast<double>(cp)));
  const size_t nf = std::max<size_t>(64, size_t(1) << trunc_i(std::ceil(std::log2(4.0 * static_cast<double>(M)))));
  const dvec win = np::hanning(static_cast<size_t>(M));
  cdvec acc(nf, cd(0.0, 0.0)), s1(static_cast<size_t>(M)), s2(static_cast<size_t>(M));
  int64_t n = 0;
  const int64_t ny = static_cast<int64_t>(y.size());
  for (double t : tsym) {
    const int64_t a = npround_int(t * ofs) + std::max<int64_t>(0, floordiv(cp - M, 2));
    if (a < 0 || a + N + M > ny) continue;
    for (int64_t i = 0; i < M; ++i) {
      const size_t ii = static_cast<size_t>(i);
      s2[ii] = np::cmul(y[static_cast<size_t>(a + N + i)], cd(win[ii], 0.0));
      s1[ii] = np::cmul(y[static_cast<size_t>(a + i)], cd(win[ii], 0.0));
    }
    const cdvec F2 = fft::fft(s2.data(), s2.size(), nf), F1 = fft::fft(s1.data(), s1.size(), nf);
    for (size_t i = 0; i < nf; ++i) acc[i] += np::cmul(F2[i], std::conj(F1[i]));
    ++n;
  }
  cdvec sh(nf);
  for (size_t i = 0; i < nf; ++i) sh[i] = acc[(i + nf - nf / 2) % nf];
  const cd tot = np::pairwise_sum(sh.data(), sh.size());
  const cd rot = expj(-std::atan2(tot.imag(), tot.real()));
  dvec fr(nf), c(nf);
  const double den = static_cast<double>(std::max<int64_t>(n, 1));
  for (size_t i = 0; i < nf; ++i) {
    c[i] = np::cmul(sh[i], rot).real() / den;
    fr[i] = off + (static_cast<double>(i) - static_cast<double>(nf) / 2) * ofs / static_cast<double>(nf);
  }
  return {fr, c};
}

std::optional<std::pair<double, double>> measure_band(const Slicer& sl, double off, double W, double tu, double ts,
                                                      const dvec& tsym, double usable) {
  std::optional<std::pair<double, double>> r;
  for (double hw0 : {W, 2 * W, 4 * W}) {
    const double hw = std::min(hw0, usable);
    if (r) off = hw < usable ? std::clamp((r->first + r->second) / 2, -usable + hw, usable - hw) : 0.0;
    const auto [fr, c] = band_cp(sl, off, hw, tu, ts, tsym);
    const dvec sm = ndi::uniform_filter1d(c, 3);
    const double p95 = np::percentile(sm, 95.0);
    dvec top;
    for (double v : sm)
      if (v > 0.5 * p95) top.push_back(v);
    if (top.empty() || np::median(top) <= 0) return std::nullopt;
    const double df = fr[1] - fr[0];
    const double thr = 0.25 * np::median(top);
    std::vector<uint8_t> occ(sm.size());
    for (size_t i = 0; i < sm.size(); ++i) occ[i] = sm[i] > thr;
    occ = ndi::binary_closing(occ, ndi::ones1(static_cast<size_t>(std::max<int64_t>(3, trunc_i(1e6 / df)))));
    const auto lab = ndi::label(occ);
    if (lab.n == 0) return std::nullopt;
    int bi = 1;
    double bs = 0;
    for (int i = 1; i <= lab.n; ++i) {
      dvec v;
      for (size_t j = 0; j < sm.size(); ++j)
        if (lab.lab[j] == i) v.push_back(sm[j]);
      const double s = np::sum(v);
      if (i == 1 || s > bs) bs = s, bi = i;
    }
    size_t i0 = sm.size(), i1 = 0;
    for (size_t j = 0; j < sm.size(); ++j)
      if (lab.lab[j] == bi) i0 = std::min(i0, j), i1 = j;
    r = std::make_pair(fr[i0], fr[i1]);
    if ((r->first - (off - hw) > 0.04 * 2 * hw && (off + hw) - r->second > 0.04 * 2 * hw) || hw >= usable) break;
  }
  return r;
}

std::vector<Link> blind_links(const Slicer& sl, double fs, double fc, double usable, double noise_pps,
                              const bursts::TFPower* tf, double keep_slice_snr, const SkipFn& skip) {
  (void)fc;
  std::vector<Hit> hits = scan(sl, usable, skip);
  std::stable_sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) {
    return std::make_tuple(-a.z, a.W, a.off, a.N) < std::make_tuple(-b.z, b.W, b.off, b.N);
  });
  std::vector<Hit> cl;
  for (const auto& h : hits) {
    bool placed = false;
    for (const auto& c : cl)
      if (std::fabs(h.tu - c.tu) < 0.015 * c.tu + 0.1e-6 && std::fabs(h.off - c.off) < std::min(h.W, c.W) / 2) {
        placed = true;
        break;
      }
    if (!placed) cl.push_back(h);
  }
  std::vector<Link> out;
  for (size_t ci = 0; ci < cl.size() && ci < MAX_CLUSTERS; ++ci) {
    const Hit& c = cl[ci];
    if (c.tu < TU_MIN || (skip && skip(c.off))) continue;
    bool inside = false;
    for (const auto& o : out) inside = inside || (std::fabs(c.tu - o.tu) < 0.015 * o.tu && o.lo <= c.off && c.off <= o.hi);
    if (inside) continue;
    const auto cf = confirm(sl, c.off, c.W, c.N, c.ofs);
    if (!cf) continue;
    const double tu = cf->tu;
    Segmented sg = segment(sl, c.off, c.W, tu, cf->ts);
    if (sg.runs.empty()) continue;
    dvec tsym;
    for (size_t i = 0; i < sg.tk.size(); ++i)
      if (sg.zr[i] > 2) tsym.push_back(sg.tk[i]);
    const auto bd = measure_band(sl, c.off, c.W, tu, sg.ts, tsym, usable);
    if (!bd) continue;
    double lo = bd->first, hi = bd->second, B = hi - lo, cen = (lo + hi) / 2;
    sg = segment(sl, cen, 1.02 * B, tu, sg.ts);
    const double ts = sg.ts;
    std::pair<std::optional<double>, double> fp{std::nullopt, 0.0};
    if (!sg.tk.empty()) fp = frame_period(sg.tk, sg.zr, ts);
    std::vector<Run> runs = sg.runs;
    if (fp.first) runs = frame_complete(runs, sg.tk, sg.zr, ts, *fp.first);
    int64_t tot = 0, longest = 0;
    for (const auto& r : runs) {
      const int64_t ns = npround_int((r.t1 - r.t0) / ts);
      tot += ns;
      longest = std::max(longest, ns);
    }
    if (runs.empty() || tot < MIN_SYMS || longest < MIN_LONGEST) continue;
    if (tf) {
      // band and burst extent from in-band power, seeded by the CP runs (core/bursts)
      bursts::Runs seeds;
      for (const auto& r : runs) seeds.emplace_back(r.t0, r.t1);
      bursts::EnergyInfo einfo;
      for (int it = 0; it < 2; ++it) {
        if (const auto fb = bursts::fit_band(*tf, seeds, lo, hi, std::max(lo - 0.5 * B, -usable), std::min(hi + 0.5 * B, usable))) {
          lo = fb->lo, hi = fb->hi;
          B = hi - lo, cen = (lo + hi) / 2;
        }
        if (it == 0) {
          auto [s2, inf] = bursts::energy_bursts(*tf, seeds, lo, hi, usable, ts);
          seeds = std::move(s2);
          einfo = inf;
          if (!einfo.used) break;
        }
      }
      if (einfo.used) {
        runs.clear();
        for (const auto& [a, b] : seeds) runs.push_back({a, b, npround_int((b - a) / ts)});
      }
    }
    bool dup = false;
    for (const auto& o : out)
      dup = dup || (std::fabs(o.tu - tu) < 0.015 * tu && std::min(o.hi, hi) - std::max(o.lo, lo) > 0.5 * std::min(o.hi - o.lo, B));
    if (dup) continue;
    double ofs2 = rate_for(tu, B).first;
    const cdvec yb = sl(cen, B / 2, ofs2, &ofs2);
    std::vector<uint8_t> inb(yb.size(), 0);
    for (const auto& r : runs) {
      const int64_t a = std::max<int64_t>(trunc_i(r.t0 * ofs2), 0);
      const int64_t b = std::min<int64_t>(trunc_i(r.t1 * ofs2), static_cast<int64_t>(yb.size()));
      for (int64_t i = a; i < b; ++i) inb[static_cast<size_t>(i)] = 1;
    }
    dvec pw;
    for (size_t i = 0; i < yb.size(); ++i)
      if (inb[i]) pw.push_back(np::cabs2(yb[i]));
    const double p_on = pw.empty() ? 0.0 : np::mean(pw);
    const double snr = 10 * std::log10(std::max(p_on / (noise_pps * B / fs + 1e-30) - 1, 1e-3));
    const double dur = static_cast<double>(yb.size()) / ofs2;
    double on_s = 0;
    for (const auto& r : runs) on_s += r.t1 - r.t0;
    Link lk;
    lk.tu = tu;
    lk.ts = ts;
    lk.lo = lo;
    lk.hi = hi;
    lk.cen = cen;
    lk.bw = B;
    lk.line_z = cf->line_z;
    lk.win_z = cf->nwin;
    lk.acf_z = c.z;
    lk.runs = runs;
    lk.n_grids = sg.ngrid;
    lk.frame_s = fp.first;
    lk.frame_contrast = fp.second;
    lk.snr_db = snr;
    lk.p_on = p_on;
    lk.duty = on_s / dur;
    lk.cfo_frac_hz = sg.frac;
    lk.ts_np = cf->ts_np;
    if (snr >= keep_slice_snr) {
      lk.slice = yb;
      lk.slice_rate = ofs2;
    }
    out.push_back(std::move(lk));
  }
  // IQ image: a mirrored copy about fc (same Tu), >= 10 dB weaker
  std::vector<Link> kept;
  for (const auto& o : out) {
    bool img = false;
    for (const auto& q : out)
      img = img || (&q != &o && std::fabs(q.tu - o.tu) < 0.015 * o.tu && std::fabs(q.cen + o.cen) < 0.5e6 && q.p_on > 10 * o.p_on);
    if (!img) kept.push_back(o);
  }
  return kept;
}

// ============================================================== hopping narrowband links (RC uplinks)
namespace {
constexpr double HOP_FAMILIES[][2] = {{66.667e-6, 144.0 / 2048}, {33.333e-6, 72.0 / 1024}};
constexpr double HOP_WIN_S = 0.25e-3;
constexpr double HOP_K = 5.0;
constexpr double HOP_TRACK_K = 3.5;
constexpr double HOP_Z = 6.0;
constexpr double HOP_CON_DB = 3.0;
constexpr int64_t HOP_MIN_SYMS = 4;
constexpr double HOP_MAX_S = 1.0e-3;
constexpr double HOP_BW0 = 0.5e6, HOP_BW1 = 8e6;
constexpr double HOP_PAD_S = 1e-3;
constexpr double HOP_DUR_TOL = 0.2, HOP_BW_TOL = 0.3;

double overlap(double a0, double a1, double b0, double b1) { return std::min(a1, b1) - std::max(a0, b0); }

// python _fold_z(y, N, cp, P)
double fold_z(const cdvec& y, size_t N, size_t cp, double P) {
  const Cp c = cp_metric(y, N, cp);
  const int64_t K = static_cast<int64_t>(std::floor(static_cast<double>(c.m.size()) / P));
  if (K < 3) return 0.0;
  const size_t Pi = static_cast<size_t>(trunc_i(P));
  dmat M(static_cast<size_t>(K), Pi);
  for (int64_t k = 0; k < K; ++k)
    for (size_t j = 0; j < Pi; ++j)
      M(static_cast<size_t>(k), j) = c.m[static_cast<size_t>(static_cast<int64_t>(static_cast<double>(k) * P + static_cast<double>(j)))];
  const dvec f = np::mean_axis0(M);
  const double med = np::median(f);
  dvec ad(f.size());
  for (size_t i = 0; i < f.size(); ++i) ad[i] = std::fabs(f[i] - med);
  const double mad = 1.4826 * np::median(ad) + 1e-12;
  return (*std::max_element(f.begin(), f.end()) - med) / mad;
}

// |fftshift(fft(rows of y[i0:i0+n*nf] * hanning))|^2, nullopt when < 1 frame
std::optional<dmat> spec_frames(const cdvec& y, double ofs, size_t nf, const dvec& win, double u0, double u1) {
  const int64_t i0 = std::max<int64_t>(trunc_i(u0 * ofs), 0);
  const int64_t i1 = std::min<int64_t>(trunc_i(u1 * ofs), static_cast<int64_t>(y.size()));
  const int64_t n = (i1 - i0) >= 0 ? (i1 - i0) / static_cast<int64_t>(nf) : -1;
  if (n < 1) return std::nullopt;
  dmat F(static_cast<size_t>(n), nf);
  cdvec v(nf);
  const size_t h = nf / 2;
  for (int64_t r = 0; r < n; ++r) {
    for (size_t i = 0; i < nf; ++i) {
      const cd& s = y[static_cast<size_t>(i0 + r * static_cast<int64_t>(nf)) + i];
      v[i] = cd(s.real() * win[i], s.imag() * win[i]);
    }
    const cdvec X = fft::fft(v);
    for (size_t j = 0; j < nf; ++j) {
      const double a = np::cabs(X[(j + nf - h) % nf]);
      F(static_cast<size_t>(r), j) = a * a;
    }
  }
  return F;
}
}  // namespace

double hop_width(double tu) { return std::max(1.0e6, 75.0 / tu); }

std::optional<Hop> refine_hop(const Slicer& sl, double f_off, double t0_s, double t1_s, double tu, double cpf, double W) {
  const double ts = tu * (1 + cpf);
  double o2 = 0;
  const cdvec yn = sl(f_off, W / 2, 1.25 * W, &o2);
  const int64_t a2 = std::max<int64_t>(0, trunc_i((t0_s - HOP_PAD_S) * o2));
  const int64_t b2 = std::min<int64_t>(static_cast<int64_t>(yn.size()), trunc_i((t1_s + HOP_PAD_S) * o2));
  const int64_t k = std::max<int64_t>(4, trunc_i(10e-6 * o2));
  const int64_t nb = static_cast<int64_t>(std::floor(static_cast<double>(b2 - a2) / static_cast<double>(k)));
  if (nb < 5) return std::nullopt;
  dvec e(static_cast<size_t>(nb));
  dvec blk(static_cast<size_t>(k));
  for (int64_t i = 0; i < nb; ++i) {
    for (int64_t j = 0; j < k; ++j) {
      const double a = np::cabs(yn[static_cast<size_t>(a2 + i * k + j)]);
      blk[static_cast<size_t>(j)] = a * a;
    }
    e[static_cast<size_t>(i)] = np::mean(blk);
  }
  const double th = std::sqrt(np::percentile(e, 90.0) * np::percentile(e, 10.0));
  std::vector<uint8_t> on(e.size());
  for (size_t i = 0; i < e.size(); ++i) on[i] = e[i] > th;
  on = ndi::binary_closing(on, ndi::ones1(3));
  int64_t pk = static_cast<int64_t>(std::floor(((t0_s + t1_s) / 2 * o2 - static_cast<double>(a2)) / static_cast<double>(k)));
  pk = std::clamp<int64_t>(pk, 0, nb - 1);
  if (!on[static_cast<size_t>(pk)]) {
    int64_t best = -1;
    for (int64_t i = 0; i < nb; ++i)
      if (on[static_cast<size_t>(i)] && (best < 0 || std::llabs(i - pk) < std::llabs(best - pk))) best = i;
    if (best < 0) return std::nullopt;
    pk = best;
  }
  int64_t s0 = pk, s1 = pk;
  while (s0 > 0 && on[static_cast<size_t>(s0 - 1)]) --s0;
  while (s1 < nb - 1 && on[static_cast<size_t>(s1 + 1)]) ++s1;
  const double t0 = static_cast<double>(a2 + s0 * k) / o2, t1 = static_cast<double>(a2 + (s1 + 1) * k) / o2;
  if (t1 - t0 < static_cast<double>(HOP_MIN_SYMS) * ts || t1 - t0 > HOP_MAX_S) return std::nullopt;
  // band: on-burst vs off-burst spectrum, re-centred once when it reaches the window edge
  double c = f_off, hw = 1.5 * W, lo = 0, hi = 0, inn = 0, outv = 0;
  bool edge = true;
  for (int it = 0; it < 2; ++it) {
    double ofs = 0;
    const cdvec y = sl(c, hw, 2.5 * hw, &ofs);
    const size_t nf = size_t(1) << trunc_i(std::ceil(std::log2(ofs * 32e-6)));
    const dvec win = np::hanning(nf);
    const auto fon = spec_frames(y, ofs, nf, win, t0, t1);
    std::vector<dmat> foff;
    if (auto f1 = spec_frames(y, ofs, nf, win, t0 - HOP_PAD_S, t0 - static_cast<double>(nf) / ofs)) foff.push_back(std::move(*f1));
    if (auto f2 = spec_frames(y, ofs, nf, win, t1 + static_cast<double>(nf) / ofs, t1 + HOP_PAD_S)) foff.push_back(std::move(*f2));
    if (!fon || foff.empty()) return std::nullopt;
    dmat off_all(0, nf);
    for (const auto& m : foff) {
      off_all.data.insert(off_all.data.end(), m.data.begin(), m.data.end());
      off_all.rows += m.rows;
    }
    const dvec mon = np::median_axis0(*fon), moff = np::median_axis0(off_all);
    const size_t h = nf / 2;
    dvec fb(nf);
    for (size_t j = 0; j < nf; ++j) fb[j] = c + static_cast<double>(static_cast<int64_t>(j) - static_cast<int64_t>(h)) * ofs / static_cast<double>(nf);
    std::vector<size_t> inw;
    for (size_t j = 0; j < nf; ++j)
      if (std::fabs(fb[j] - c) <= hw) inw.push_back(j);
    dvec a_on(inw.size()), a_off(inw.size());
    for (size_t i = 0; i < inw.size(); ++i) a_on[i] = mon[inw[i]], a_off[i] = moff[inw[i]];
    const dvec S_on = ndi::uniform_filter1d(a_on, 3), S_off = ndi::uniform_filter1d(a_off, 3);
    dvec d(inw.size());
    for (size_t i = 0; i < d.size(); ++i) d[i] = 10 * std::log10(std::max(S_on[i], 1e-30) / std::max(S_off[i], 1e-30));
    const bursts::RectFit r = bursts::rect_fit(d);
    const double half_bin = ofs / static_cast<double>(nf) / 2;
    lo = fb[inw[r.lo]] - half_bin, hi = fb[inw[r.hi - 1]] + half_bin;
    inn = r.inside, outv = r.outside;
    edge = r.lo <= 1 || r.hi >= inw.size() - 1;
    if (!edge) break;
    c = (lo + hi) / 2, hw = std::max(hw, hi - lo);
  }
  if (edge) return std::nullopt;
  const double bw = hi - lo;
  if (!(HOP_BW0 <= bw && bw <= HOP_BW1)) return std::nullopt;
  auto [o3, N] = rate_for(tu, bw, tu * cpf);
  const size_t cp = static_cast<size_t>(std::max<int64_t>(1, py_round_int(static_cast<double>(N) * cpf)));
  const cdvec y3 = sl((lo + hi) / 2, bw / 2, o3, &o3);
  const int64_t i0 = std::max<int64_t>(0, trunc_i(t0 * o3)), i1 = std::min<int64_t>(trunc_i(t1 * o3), static_cast<int64_t>(y3.size()));
  const cdvec seg = i1 > i0 ? cdvec(y3.begin() + i0, y3.begin() + i1) : cdvec();
  Hop out;
  out.t0 = t0, out.t1 = t1, out.lo = lo, out.hi = hi;
  out.z = fold_z(seg, static_cast<size_t>(N), cp, static_cast<double>(N) * (1 + cpf));
  out.contrast = inn - outv;
  return out;
}

std::optional<double> hop_period(const dvec& starts, double lo, double hi, double tol) {
  dvec s = starts;
  std::sort(s.begin(), s.end());
  dvec d;
  for (size_t i = 0; i < s.size(); ++i)
    for (size_t j = i + 1; j < s.size(); ++j) d.push_back(s[j] - s[i]);
  if (d.empty()) return std::nullopt;
  int64_t best = 0;
  std::optional<double> Tb;
  for (double T : np::arange(lo, hi + 1e-12, 5e-6)) {
    int64_t n = 0;
    for (double v : d) {
      const double q = std::nearbyint(v / T);
      n += (q >= 1 && std::fabs(v - q * T) < tol);
    }
    if (n >= best) best = n, Tb = T;
  }
  return best >= 3 ? Tb : std::nullopt;
}

std::vector<HopLink> blind_hops(const Slicer& sl, double fs, double fc, double usable, const std::vector<Link>& links,
                                const SkipFn& skip) {
  if (skip) {  // the whole usable band skipped: nothing to scan (the scan threshold is a whole-grid statistic)
    bool all = true;
    for (double f : np::arange(-usable, usable + 1, 1e6))
      if (!skip(f)) {
        all = false;
        break;
      }
    if (all) return {};
  }
  (void)fs;
  (void)fc;
  struct Num {
    double tu, cpf;
    std::optional<double> T;
  };
  std::vector<Num> nums;
  auto add = [&](double tu, double cpf, std::optional<double> T) {
    if (tu < TU_MIN) return;
    for (const auto& n : nums)
      if (std::fabs(tu - n.tu) < 0.015 * n.tu) return;
    nums.push_back({tu, cpf, T});
  };
  for (const auto& lk : links) add(lk.tu, (lk.ts - lk.tu) / lk.tu, lk.frame_s);
  for (const auto& h : HOP_FAMILIES) add(h[0], h[1], std::nullopt);
  const ofdm15k::Slicer ch = [&](double f, double h, double r) { return sl(f, h, r, nullptr); };
  std::optional<double> dur;
  struct Found {
    double tu, cpf, ts;
    std::optional<double> T;
    std::vector<const Link*> own;
    std::vector<Hop> hops;
  };
  std::vector<Found> found;
  for (const auto& nm : nums) {
    const double tu = nm.tu, cpf = nm.cpf;
    const double W = hop_width(tu);
    const auto [ofs, N] = rate_for(tu, W, tu * cpf);
    const size_t cp = static_cast<size_t>(std::max<int64_t>(1, py_round_int(static_cast<double>(N) * cpf)));
    const double ts = tu * (1 + cpf);
    const ofdm15k::Scan sc = ofdm15k::ofdm15k_scan(ch, {{"hop", ofs, static_cast<size_t>(N), cp, W, W / 2}}, HOP_WIN_S,
                                                   HOP_K, usable);
    std::vector<const Link*> own;
    for (const auto& lk : links)
      if (std::fabs(lk.tu - tu) < 0.015 * tu) own.push_back(&lk);
    if (!dur) dur = static_cast<double>(ch(0.0, W / 2, 1.25 * W).size()) / (1.25 * W);
    auto inside_link = [&](const Hop& r) {
      const double c = (r.lo + r.hi) / 2;
      for (const Link* lk : own) {
        if (!(lk->lo <= c && c <= lk->hi)) continue;
        for (const auto& ru : lk->runs)
          if (overlap(r.t0, r.t1, ru.t0, ru.t1) > 0.5 * (r.t1 - r.t0)) return true;
      }
      return false;
    };
    auto dup = [&](const Hop& r, const std::vector<Hop>& hops) {
      for (const auto& h : hops)
        if (overlap(r.t0, r.t1, h.t0, h.t1) > 0.5 * std::min(r.t1 - r.t0, h.t1 - h.t0) && overlap(r.lo, r.hi, h.lo, h.hi) > 0)
          return true;
      return false;
    };
    std::vector<ofdm15k::Hit> hits = sc.hits;
    std::stable_sort(hits.begin(), hits.end(), [](const ofdm15k::Hit& a, const ofdm15k::Hit& b) {
      return std::make_tuple(-a.score, a.t0_s, a.f_off) < std::make_tuple(-b.score, b.t0_s, b.f_off);
    });
    std::vector<Hop> hops;
    for (const auto& h : hits) {
      bool covered = false;
      for (const auto& g : hops) covered = covered || (h.t0_s < g.t1 && h.t1_s > g.t0 && g.lo <= h.f_off && h.f_off <= g.hi);
      if (covered || (skip && skip(h.f_off))) continue;
      auto r = refine_hop(sl, h.f_off, h.t0_s, h.t1_s, tu, cpf, W);
      if (!r || r->z < HOP_Z || r->contrast < HOP_CON_DB || inside_link(*r) || dup(*r, hops)) continue;
      r->tracked = false;
      hops.push_back(*r);
    }
    // hops the scan missed: re-test the times the frame cadence predicts, at a lower threshold
    if (hops.size() >= 2) {
      dvec st;
      for (const auto& h : hops) st.push_back(h.t0);
      std::optional<double> Tc = hop_period(st);
      if (!Tc) Tc = nm.T;
      if (Tc) {
        dvec dd, bb;
        for (const auto& h : hops) dd.push_back(h.t1 - h.t0), bb.push_back(h.hi - h.lo);
        const double dmed = np::median(dd), bmed = np::median(bb);
        std::set<double> cs;
        const int64_t kk = trunc_i(*dur / *Tc);
        for (const auto& h : hops)
          for (int64_t k = -kk - 1; k < kk + 2; ++k) {
            const double t = h.t0 + static_cast<double>(k) * *Tc;
            if (0 <= t && t <= *dur - dmed) cs.insert(py_round(t, 6));
          }
        dvec times;
        for (double t : cs) {
          bool ov = false;
          for (const auto& g : hops) ov = ov || overlap(t, t + dmed, g.t0, g.t1) > 0;
          if (ov) continue;
          if (times.empty() || t - times.back() > 0.1e-3) times.push_back(t);
        }
        const auto& stt = sc.stats.at("hop");
        const double thr = stt.median + HOP_TRACK_K * stt.sigma;
        const int64_t Wn = trunc_i(dmed * ofs);
        const dvec offs = np::arange(-usable + W / 2, usable - W / 2 + 1, W / 2);
        std::vector<std::pair<double, cdvec>> slices;
        for (double o : offs) slices.emplace_back(o, ch(o, W / 2 * 1.05, ofs));
        const dvec dts = np::arange(-0.15e-3, 0.151e-3, 0.05e-3);
        for (double t : times) {
          std::vector<std::tuple<double, double, double>> scs;
          for (const auto& [o, y] : slices) {
            double bs = -1.0, bt = 0.0;
            for (double dt : dts) {
              const int64_t a = trunc_i((t + dt) * ofs);
              if (a < 0 || a + Wn + N + static_cast<int64_t>(cp) >= static_cast<int64_t>(y.size())) continue;
              const cdvec seg(y.begin() + a, y.begin() + a + Wn + N + static_cast<int64_t>(cp));
              const Cp m = cp_metric(seg, static_cast<size_t>(N), cp);
              const dvec f = folded_scores(m.m, static_cast<size_t>(N) + cp, static_cast<size_t>(Wn), static_cast<size_t>(Wn));
              if (!f.empty() && f[0] > bs) bs = f[0], bt = t + dt;
            }
            if (bs > thr) scs.emplace_back(bs, o, bt);
          }
          std::sort(scs.begin(), scs.end(), std::greater<>());
          for (size_t i = 0; i < scs.size() && i < 4; ++i) {
            const double o = std::get<1>(scs[i]), t0 = std::get<2>(scs[i]);
            auto r = refine_hop(sl, o, t0, t0 + dmed, tu, cpf, W);
            if (!r || r->z < HOP_Z || r->contrast < HOP_CON_DB || inside_link(*r) || dup(*r, hops) ||
                std::fabs((r->t1 - r->t0) - dmed) > 0.25 * dmed || std::fabs((r->hi - r->lo) - bmed) > 0.35 * bmed)
              continue;
            r->tracked = true;
            hops.push_back(*r);
            break;
          }
        }
      }
    }
    for (auto& h : hops) h.num = found.size();
    found.push_back({tu, cpf, ts, nm.T, own, std::move(hops)});
  }
  // one hop, one numerology: the stronger CP fold wins where hops of two numerologies coincide
  std::vector<Hop> all;
  for (const auto& f : found) all.insert(all.end(), f.hops.begin(), f.hops.end());
  std::stable_sort(all.begin(), all.end(), [](const Hop& a, const Hop& b) {
    return std::make_tuple(-a.z, a.t0, a.lo) < std::make_tuple(-b.z, b.t0, b.lo);
  });
  std::vector<Hop> kept;
  for (const auto& h : all) {
    bool ov = false;
    for (const auto& g : kept)
      ov = ov || (overlap(h.t0, h.t1, g.t0, g.t1) > 0.5 * std::min(h.t1 - h.t0, g.t1 - g.t0) && overlap(h.lo, h.hi, g.lo, g.hi) > 0);
    if (!ov) kept.push_back(h);
  }
  std::vector<HopLink> out;
  for (size_t i = 0; i < found.size(); ++i) {
    const Found& f = found[i];
    std::vector<Hop> hops;
    for (const auto& h : kept)
      if (h.num == i) hops.push_back(h);
    std::stable_sort(hops.begin(), hops.end(), [](const Hop& a, const Hop& b) { return a.t0 < b.t0; });
    for (int it = 0; it < 2 && hops.size() >= 2; ++it) {
      dvec dd, bb;
      for (const auto& h : hops) dd.push_back(h.t1 - h.t0), bb.push_back(h.hi - h.lo);
      const double dm = np::median(dd), bm = np::median(bb);
      std::vector<Hop> k2;
      for (const auto& h : hops)
        if (std::fabs(h.t1 - h.t0 - dm) <= HOP_DUR_TOL * dm && std::fabs(h.hi - h.lo - bm) <= HOP_BW_TOL * bm) k2.push_back(h);
      hops = std::move(k2);
    }
    if (hops.size() < 2) continue;
    // hop channels: hops whose bands overlap by more than half
    std::vector<size_t> ord(hops.size());
    std::iota(ord.begin(), ord.end(), 0);
    std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return (hops[a].lo + hops[a].hi) / 2 < (hops[b].lo + hops[b].hi) / 2; });
    struct Chan {
      double lo, hi;
      std::vector<size_t> members;
    };
    std::vector<Chan> chans;
    for (size_t j : ord) {
      const Hop& h = hops[j];
      bool placed = false;
      for (auto& c : chans)
        if (overlap(h.lo, h.hi, c.lo, c.hi) > 0.5 * std::min(h.hi - h.lo, c.hi - c.lo)) {
          c.members.push_back(j);
          placed = true;
          break;
        }
      if (!placed) chans.push_back({h.lo, h.hi, {j}});
    }
    for (const auto& c : chans) {
      dvec cs;
      for (size_t j : c.members) cs.push_back((hops[j].lo + hops[j].hi) / 2);
      const double cen = np::median(cs);
      for (size_t j : c.members) hops[j].channel = cen;
    }
    dvec st, dd, bb;
    for (const auto& h : hops) st.push_back(h.t0), dd.push_back(h.t1 - h.t0), bb.push_back(h.hi - h.lo);
    HopLink hl;
    hl.tu = f.tu, hl.ts = f.ts, hl.cpf = f.cpf;
    hl.n_channels = static_cast<int64_t>(chans.size());
    hl.frame_s = hop_period(st);
    if (!hl.frame_s) hl.frame_s = f.T;
    for (const Link* lk : f.own) hl.paired.push_back(lk->cen);
    hl.dur_s = np::median(dd);
    hl.bw = np::median(bb);
    hl.hops = std::move(hops);
    out.push_back(std::move(hl));
  }
  return out;
}

}  // namespace si::rx::digital_vtx::blind
