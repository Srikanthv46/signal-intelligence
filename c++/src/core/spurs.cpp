#include "core/spurs.hpp"

#include <algorithm>
#include <cmath>
#include <string>

#include "core/npcompat.hpp"

namespace si::spurs {

namespace {

bool strong(const bands::Band& b) {
  static const char* kStrong[] = {"lte_dl", "nr_dl", "gsm_dl", "umts_dl", "lte_tdd", "nr_tdd"};
  for (const auto& c : b.classes)
    for (const char* s : kStrong)
      if (c == s) return true;
  return false;
}

// Python floor division / modulo on int64
inline int64_t floordiv(int64_t a, int64_t b) { return a / b - ((a % b != 0) && ((a < 0) != (b < 0))); }
inline int64_t pymod(int64_t a, int64_t b) {
  const int64_t r = a % b;
  return r < 0 ? r + b : r;
}

}  // namespace

std::pair<double, double> refine_line(const cfvec& X, double fs, double fc, double f_abs, double half_hz) {
  const int64_t L = static_cast<int64_t>(X.size());
  const double df = fs / static_cast<double>(L);
  const int64_t k0 = py_round_int((f_abs - fc) / df);
  const int64_t h = std::max<int64_t>(2, static_cast<int64_t>(std::ceil(half_hz / df)));
  std::vector<int64_t> k;
  for (int64_t i = k0 - h; i <= k0 + h; ++i)
    if (i > floordiv(-L, 2) && i < floordiv(L, 2)) k.push_back(i);
  if (k.size() < 5) return {f_abs, 0.0};
  dvec p(k.size());
  for (size_t i = 0; i < k.size(); ++i) {
    const double a = static_cast<double>(np::cabs(X[static_cast<size_t>(pymod(k[i], L))]));
    p[i] = a * a;
  }
  const size_t i = np::argmax(p.data(), p.size());
  double off = 0.0;
  if (i > 0 && i + 1 < p.size()) {
    const double a = std::log(p[i - 1] + 1e-30), b = std::log(p[i] + 1e-30), c = std::log(p[i + 1] + 1e-30);
    const double den = a - 2 * b + c;
    off = den < 0 ? 0.5 * (a - c) / den : 0.0;
  }
  dvec q = p;
  const double med = np::median_inplace(q.data(), q.size());
  return {fc + (static_cast<double>(k[i]) + off) * df, 10 * std::log10(p[i] / (med + 1e-30))};
}

std::optional<int64_t> sample_clock_harmonic(double f_abs, double fs, double tol_hz) {
  const int64_t k = py_round_int(f_abs / fs);
  if (k >= 2 && std::fabs(f_abs - static_cast<double>(k) * fs) <= tol_hz) return k;
  return std::nullopt;
}

bool reference_locked(double f_abs, std::optional<double> ppm, double tol_hz) {
  if (!ppm || std::fabs(*ppm) * 1e-6 * f_abs < 8 * tol_hz) return false;
  return std::fabs(f_abs - static_cast<double>(py_round_int(f_abs / kRefGridHz)) * kRefGridHz) <= tol_hz;
}

Json receiver_lines(const cfvec& X, double fs, double fc, std::optional<double> ppm, std::pair<double, double> band,
                    double win_hz) {
  const double lo = band.first, hi = band.second;
  struct Nom {
    double f;
    bool clock;
    int64_t k;
  };
  std::vector<Nom> nominal;
  for (int64_t k = std::max<int64_t>(2, static_cast<int64_t>(std::ceil(lo / fs)));
       k <= static_cast<int64_t>(std::floor(hi / fs)); ++k)
    nominal.push_back({static_cast<double>(k) * fs, true, k});
  if (ppm)
    for (int64_t n = static_cast<int64_t>(std::ceil(lo / kRefGridHz)); n <= static_cast<int64_t>(std::floor(hi / kRefGridHz));
         ++n)
      nominal.push_back({static_cast<double>(n) * kRefGridHz, false, 0});
  Json out = Json::array();
  for (const Nom& nm : nominal) {
    if (!(lo <= nm.f && nm.f <= hi) || std::fabs(nm.f - fc) < 3 * win_hz) continue;  // the DC bin is LO leakage
    const auto [f, pk] = refine_line(X, fs, fc, nm.f, win_hz);
    if (pk < kLineMinDb) continue;
    bool ok;
    if (nm.clock) {
      const auto k = sample_clock_harmonic(f, fs);
      ok = k && *k == nm.k;
    } else {
      ok = reference_locked(f, ppm);
    }
    if (ok && std::fabs(f - nm.f) <= kLineTolHz)
      out.push_back({{"freq_hz", np::round(f, 1)},
                     {"kind", nm.clock ? "sample_clock_harmonic" : "reference_locked"},
                     {"harmonic", nm.clock ? Json(nm.k) : Json(nullptr)},
                     {"peak_db", py_round(pk, 1)}});
  }
  return out;
}

Json harmonic_image_sources(double f0, double f1, double lo) {
  Json out = Json::array();
  for (int n : {3, 5}) {
    const double d0 = f0 - lo, d1 = f1 - lo;
    const bool inv = n % 4 == 3;
    double a = inv ? n * lo - d0 : n * lo + d0, b = inv ? n * lo - d1 : n * lo + d1;
    if (b < a) std::swap(a, b);
    for (const auto& band : bands::bands()) {
      if (strong(band) && a >= band.lo && b <= band.hi) {
        out.push_back({{"n", n},
                       {"inverted", inv},
                       {"source_low_hz", py_round_int(a)},
                       {"source_high_hz", py_round_int(b)},
                       {"band", band.key},
                       {"band_description", band.desc}});
        break;
      }
    }
  }
  return out;
}

bool in_strong_band(double f0, double f1) {
  for (const auto& b : bands::bands())
    if (strong(b) && f0 < b.hi && f1 > b.lo) return true;
  return false;
}

}  // namespace si::spurs
