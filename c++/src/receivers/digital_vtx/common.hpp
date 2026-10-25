// Internal helpers shared by the Herelink and generic halves of the digital_vtx detector.
#pragma once

#include <cmath>
#include <cstdint>

#include "core/npcompat.hpp"
#include "core/types.hpp"

namespace si::rx::digital_vtx::detail {

constexpr double kPi = 3.141592653589793;

// np.exp(0 + 1j*im): libm cexp (on Apple __sincos, which may differ from sin() in the last bit)
inline cd expj(double im) {
#ifdef __APPLE__
  double s, c;
  __sincos(im, &s, &c);
  return {c, s};
#else
  return {std::cos(im), std::sin(im)};
#endif
}

inline int64_t pymod(int64_t a, int64_t b) {
  const int64_t r = a % b;
  return r < 0 ? r + b : r;
}

// Python float // float
inline double py_floordiv(double vx, double wx) {
  const double mod = std::fmod(vx, wx);
  double div = (vx - mod) / wx;
  if (mod != 0.0 && ((wx < 0) != (mod < 0))) div -= 1.0;
  if (div == 0.0) return std::copysign(0.0, vx / wx);
  double fl = std::floor(div);
  if (div - fl > 0.5) fl += 1.0;
  return fl;
}

// round(v, n) where v is a numpy float64 (np.round: rint(v*10^n)/10^n) or a Python float
inline double rnd(double v, int n, bool numpy) { return numpy ? np::round(v, n) : py_round(v, n); }

// python int(v) (truncation toward zero)
inline int64_t trunc_i(double v) { return static_cast<int64_t>(v); }

}  // namespace si::rx::digital_vtx::detail

#include <vector>

#include "core/model.hpp"

namespace si::rx::digital_vtx::detail {
// python herelink_detect(x, fs, fc, X, ctx, info) / generic_detect(..., exclude)
std::vector<Emission> herelink_detect(size_t n_samples, double fs, double fc, const cdvec& X);
std::vector<Emission> generic_detect(const cfvec& x, double fs, double fc, const cdvec& X,
                                     const std::vector<Emission>& exclude);
// sum_k a[k] * conj(b[k]) sequentially (numpy correlate's dot on the conjugated template)
inline cd dot_conj(const cd* a, const cd* b, size_t n) {
  cd s(0.0, 0.0);
  for (size_t i = 0; i < n; ++i) s += np::cmul(a[i], std::conj(b[i]));
  return s;
}
// np.convolve(v, np.ones(nt), "valid"): sequential window sums
inline dvec box_valid(const double* v, size_t n, size_t nt) {
  if (n < nt) return {};
  dvec e(n - nt + 1);
  for (size_t k = 0; k + nt <= n; ++k) {
    double s = 0.0;
    for (size_t i = 0; i < nt; ++i) s += v[k + i] * 1.0;
    e[k] = s;
  }
  return e;
}
}  // namespace si::rx::digital_vtx::detail
