// numpy functions with numpy's exact semantics (dtype rules, summation order, edge
// cases), for the ports of the Python receivers. Float32 inputs stay float32 where numpy
// keeps float32 (NEP 50: Python scalars are weak), float64 stays float64.
//
// Exactness notes (numpy 2.2 on aarch64 / x86-64 with FMA):
//  * np.sum / np.mean / np.std over a contiguous axis use numpy's pairwise summation
//    (pairwise_sum): bit-identical. Reductions over a non-contiguous axis (mean(axis=0)
//    of a C-order matrix) add rows sequentially: *_axis0 do the same.
//  * np.abs(complex) is numpy's SIMD formula  max * sqrt(fma(min/max, min/max, 1))
//    (cabs, not std::abs/hypot); complex multiply is fma-based (cmul). Both verified
//    bit-exact against numpy on this machine.
//  * dot-product based functions (convolve, correlate, vdot, norm, BLAS solves) cannot
//    reproduce OpenBLAS's SIMD summation order: they agree to rounding (~1e-15 rel. for
//    float64, ~1e-7 for float32).
//  * argsort: argsort_quick ports numpy's default introsort (same order of ties as numpy
//    for non-SIMD dtypes); argsort_stable = kind="stable".
#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "core/types.hpp"

namespace si::np {

template <class T> struct is_cplx : std::false_type {};
template <class T> struct is_cplx<std::complex<T>> : std::true_type {};
template <class T> struct real_of { using type = T; };
template <class T> struct real_of<std::complex<T>> { using type = T; };
template <class T> using real_t = typename real_of<T>::type;

// ------------------------------------------------------------------ complex scalars
// np.abs(z) for complex64/complex128 arrays (numpy's SIMD loop):
// larger * sqrt(fma(ratio, ratio, 1)), ratio = smaller / larger; inf/nan as numpy.
// Branch-free (selects only) so loops over it vectorise; the same value as the branchy form
// for every input, inf/nan included.
template <class R>
inline R cabs(const std::complex<R>& z) {
  const R re = std::fabs(z.real()), im = std::fabs(z.imag());
  const R inf = std::numeric_limits<R>::infinity();
  const R larger = std::max(re, im), smaller = std::min(re, im);
  const R ratio = smaller / (larger == R(0) ? R(1) : larger);  // 0 + 0j -> 0
  const R a = std::sqrt(std::fma(ratio, ratio, R(1))) * larger;
  const R a_nan = (std::isnan(re) || std::isnan(im)) ? std::numeric_limits<R>::quiet_NaN() : a;
  return (re == inf || im == inf) ? inf : a_nan;
}
// np.abs(z) ** 2 (the square of the rounded magnitude, as numpy computes it)
template <class R>
inline R cabs2(const std::complex<R>& z) { R a = cabs(z); return a * a; }

// numpy complex multiply (fma-based SIMD loop): re = fma(ar, br, -(ai*bi)),
// im = fma(ar, bi, ai*br).
template <class R>
inline std::complex<R> cmul(const std::complex<R>& a, const std::complex<R>& b) {
  return {std::fma(a.real(), b.real(), -(a.imag() * b.imag())),
          std::fma(a.real(), b.imag(), a.imag() * b.real())};
}
// out[i] = cmul(a[i], b[i]) (vectorises; out may be a or b)
template <class R>
inline void cmul_n(const std::complex<R>* a, const std::complex<R>* b, std::complex<R>* out, size_t n) {
  const R* pa = reinterpret_cast<const R*>(a);
  const R* pb = reinterpret_cast<const R*>(b);
  R* po = reinterpret_cast<R*>(out);
  for (size_t i = 0; i < n; ++i) {
    const R ar = pa[2 * i], ai = pa[2 * i + 1], br = pb[2 * i], bi = pb[2 * i + 1];
    po[2 * i] = std::fma(ar, br, -(ai * bi));
    po[2 * i + 1] = std::fma(ar, bi, ai * br);
  }
}
// numpy complex divide (Smith's algorithm, npy loops).
template <class R>
inline std::complex<R> cdiv(const std::complex<R>& a, const std::complex<R>& b) {
  const R ar = a.real(), ai = a.imag(), br = b.real(), bi = b.imag();
  const R abr = std::fabs(br), abi = std::fabs(bi);
  if (abr >= abi) {
    if (abr == 0 && abi == 0) return {ar / abr, ai / abr};
    const R rat = bi / br;
    const R scl = R(1) / (br + bi * rat);
    return {(ar + ai * rat) * scl, (ai - ar * rat) * scl};
  }
  const R rat = br / bi;
  const R scl = R(1) / (bi + br * rat);
  return {(ar * rat + ai) * scl, (ai * rat - ar) * scl};
}
// np.sqrt(complex): C99 csqrt (FreeBSD msun algorithm, as npymath's fallback).
template <class R>
std::complex<R> csqrt(const std::complex<R>& z);

// Python/numpy float remainder (np.mod / np.remainder / Python %): sign of the divisor.
template <class T>
inline T mod(T a, T b) {
  T m = std::fmod(a, b);
  if (b == T(0)) return m;
  if (m != T(0)) {
    if ((b < T(0)) != (m < T(0))) m += b;
  } else {
    m = std::copysign(T(0), b);
  }
  return m;
}

// np.round(x, decimals) / round(numpy_float64, n): rint(x * 10**n) / 10**n (not Python's
// correctly rounded decimal round(), which py_round reproduces for Python floats)
inline double round(double x, int decimals) {
  if (decimals >= 0) {
    const double f = std::pow(10.0, decimals);
    return std::nearbyint(x * f) / f;
  }
  const double f = std::pow(10.0, -decimals);
  return std::nearbyint(x / f) * f;
}

// np.sinc(x) = sin(pi*y)/(pi*y), y = 1e-20 where x == 0
inline double sinc(double x) {
  double y = kPi * (x == 0 ? 1.0e-20 : x);
  return std::sin(y) / y;
}

// ------------------------------------------------------------------ pairwise summation
// numpy's pairwise_sum (umath loops): blocks of 8 accumulators up to 128 elements,
// recursive halving above. a[i * stride].
// numpy's reduction iterator hands the loop at most NPY_BUFSIZE = 8192 elements at a time,
// so a longer sum is 0 + pw(block 0) + pw(block 1) + ... (verified for float32/float64/complex64).
constexpr size_t kNpyBufsize = 8192;
template <class T>
T pairwise_sum(const T* a, size_t n, ptrdiff_t stride = 1) {
  static_assert(std::is_floating_point_v<T>, "real only");
  if (n > kNpyBufsize) {
    T res = T(0);
    for (size_t k = 0; k < n; k += kNpyBufsize)
      res += pairwise_sum(a + static_cast<ptrdiff_t>(k) * stride, std::min(kNpyBufsize, n - k), stride);
    return res;
  }
  if (n < 8) {
    T res = T(0);  // numpy: res = 0.; res += a[i] (signed zero aside)
    for (size_t i = 0; i < n; ++i) res += a[i * stride];
    return res;
  }
  if (n <= 128) {
    T r[8];
    for (int j = 0; j < 8; ++j) r[j] = a[j * stride];
    size_t i;
    for (i = 8; i < n - (n % 8); i += 8)
      for (int j = 0; j < 8; ++j) r[j] += a[(i + j) * stride];
    T res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
    for (; i < n; ++i) res += a[i * stride];
    return res;
  }
  size_t n2 = n / 2;
  n2 -= n2 % 8;
  return pairwise_sum(a, n2, stride) + pairwise_sum(a + n2 * stride, n - n2, stride);
}

// numpy's complex pairwise sum (counts in real lanes: 4 complex accumulators).
template <class R>
std::complex<R> pairwise_sum(const std::complex<R>* a, size_t count) {
  // n = number of real lanes (2 * count), as in numpy's CDOUBLE_pairwise_sum
  struct Impl {
    static void run(const R* f, size_t n, R* rr, R* ri) {
      if (n < 8) {
        *rr = R(0); *ri = R(0);
        for (size_t i = 0; i < n; i += 2) { *rr += f[i]; *ri += f[i + 1]; }
        return;
      }
      if (n <= 128) {
        R r[8];
        for (int j = 0; j < 8; ++j) r[j] = f[j];
        size_t i;
        for (i = 8; i < n - (n % 8); i += 8)
          for (int j = 0; j < 8; ++j) r[j] += f[i + j];
        *rr = ((r[0] + r[2]) + (r[4] + r[6]));
        *ri = ((r[1] + r[3]) + (r[5] + r[7]));
        for (; i < n; i += 2) { *rr += f[i]; *ri += f[i + 1]; }
        return;
      }
      size_t n2 = n / 2;
      n2 -= n2 % 8;
      R rr1, ri1, rr2, ri2;
      run(f, n2, &rr1, &ri1);
      run(f + n2, n - n2, &rr2, &ri2);
      *rr = rr1 + rr2; *ri = ri1 + ri2;
    }
  };
  R rr = R(0), ri = R(0);
  for (size_t k = 0; k < count; k += kNpyBufsize) {   // NPY_BUFSIZE blocks, as the real sum
    R br, bi;
    Impl::run(reinterpret_cast<const R*>(a + k), 2 * std::min(kNpyBufsize, count - k), &br, &bi);
    if (k == 0) { rr = br; ri = bi; } else { rr += br; ri += bi; }
  }
  return {rr, ri};
}

// np.sum(x) of a 1-D array
template <class T> T sum(const T* x, size_t n) { return pairwise_sum(x, n); }
template <class T> T sum(const std::vector<T>& x) { return pairwise_sum(x.data(), x.size()); }

// np.mean(x): sum / n in the array's dtype (complex: numpy's complex division by n+0j,
// which multiplies by 1/n).
template <class T>
T mean(const T* x, size_t n) {
  if constexpr (is_cplx<T>::value) {
    using R = real_t<T>;
    return cdiv(pairwise_sum(x, n), T(R(n), R(0)));
  } else {
    return pairwise_sum(x, n) / static_cast<T>(n);
  }
}
template <class T> T mean(const std::vector<T>& x) { return mean(x.data(), x.size()); }

// np.var(x, ddof) / np.std(x, ddof): mean, deviations squared (|d|^2 = re^2 + im^2 for
// complex), pairwise sum, / (n - ddof). Result in the real dtype.
template <class T>
real_t<T> var(const T* x, size_t n, int ddof = 0) {
  using R = real_t<T>;
  const T m = mean(x, n);
  std::vector<R> d(n);
  for (size_t i = 0; i < n; ++i) {
    if constexpr (is_cplx<T>::value) {
      T v = x[i] - m;
      R a = v.real() * v.real(), b = v.imag() * v.imag();
      d[i] = a + b;
    } else {
      T v = x[i] - m;
      d[i] = v * v;
    }
  }
  const R s = pairwise_sum(d.data(), n);
  const long dd = static_cast<long>(n) - ddof;
  return s / static_cast<R>(dd > 0 ? dd : 0);
}
template <class T> real_t<T> var(const std::vector<T>& x, int ddof = 0) { return var(x.data(), x.size(), ddof); }
template <class T> real_t<T> stddev(const T* x, size_t n, int ddof = 0) { return std::sqrt(var(x, n, ddof)); }
template <class T> real_t<T> stddev(const std::vector<T>& x, int ddof = 0) { return std::sqrt(var(x.data(), x.size(), ddof)); }

// m.sum(axis=1) / m.mean(axis=1): pairwise per (contiguous) row.
template <class T>
std::vector<T> sum_axis1(const Mat<T>& m) {
  std::vector<T> out(m.rows);
  for (size_t r = 0; r < m.rows; ++r) out[r] = pairwise_sum(m.row(r), m.cols);
  return out;
}
template <class T>
std::vector<T> mean_axis1(const Mat<T>& m) {
  std::vector<T> out(m.rows);
  for (size_t r = 0; r < m.rows; ++r) out[r] = mean(m.row(r), m.cols);
  return out;
}
// m.sum(axis=0) / m.mean(axis=0): numpy adds the rows one after the other.
template <class T>
std::vector<T> sum_axis0(const Mat<T>& m) {
  std::vector<T> out(m.cols, T(0));
  for (size_t r = 0; r < m.rows; ++r) {
    const T* p = m.row(r);
    for (size_t c = 0; c < m.cols; ++c) out[c] += p[c];
  }
  return out;
}
template <class T>
std::vector<T> mean_axis0(const Mat<T>& m) {
  std::vector<T> out = sum_axis0(m);
  for (auto& v : out) {
    if constexpr (is_cplx<T>::value) v = cdiv(v, T(real_t<T>(m.rows), 0));
    else v = v / static_cast<T>(m.rows);
  }
  return out;
}

// np.cumsum (sequential, in the dtype)
template <class T>
std::vector<T> cumsum(const std::vector<T>& x) {
  std::vector<T> y(x.size());
  T acc = T(0);
  for (size_t i = 0; i < x.size(); ++i) { acc += x[i]; y[i] = acc; }
  if (!x.empty()) y[0] = x[0];
  return y;
}

// np.diff (n=1)
template <class T>
std::vector<T> diff(const std::vector<T>& x) {
  std::vector<T> y(x.size() > 0 ? x.size() - 1 : 0);
  for (size_t i = 0; i + 1 < x.size(); ++i) y[i] = x[i + 1] - x[i];
  return y;
}

// ------------------------------------------------------------------ order statistics
// numpy's ordering for floats: NaN sorts last.
template <class T>
inline bool np_less(T a, T b) { return a < b || (b != b && a == a); }

namespace detail {
// order-preserving unsigned key of a float (-0.0 sorts just below +0.0)
inline uint64_t order_key(double v) {
  uint64_t u;
  std::memcpy(&u, &v, sizeof u);
  return (u >> 63) ? ~u : (u | (uint64_t(1) << 63));
}
inline uint32_t order_key(float v) {
  uint32_t u;
  std::memcpy(&u, &v, sizeof u);
  return (u >> 31) ? ~u : (u | (uint32_t(1) << 31));
}
}  // namespace detail

// np.median(x) (copy): middle element, or mean of the two middle ones ((a+b)/2 in the
// dtype). NaN anywhere -> NaN. Large float arrays: a radix histogram (<= 4096 bins over the
// key range) finds the bin(s) holding the middle rank(s), then nth_element runs on those
// elements only; the value is the same as nth_element over the whole array.
template <class T>
T median_inplace(T* x, size_t n) {
  if (n == 0) return std::numeric_limits<T>::quiet_NaN();
  bool nan = false;
  for (size_t i = 0; i < n; ++i) nan |= x[i] != x[i];
  if (nan) return std::numeric_limits<T>::quiet_NaN();
  const size_t h = n / 2;
  if constexpr (std::is_floating_point_v<T>) {
    if (n >= 1024) {
      using U = decltype(detail::order_key(T()));
      U kmin = std::numeric_limits<U>::max(), kmax = 0;
      for (size_t i = 0; i < n; ++i) {
        const U k = detail::order_key(x[i]);
        kmin = std::min(kmin, k);
        kmax = std::max(kmax, k);
      }
      int shift = 0;
      while (((kmax - kmin) >> shift) >= 4096) ++shift;
      uint32_t hist[4096] = {};
      for (size_t i = 0; i < n; ++i) ++hist[(detail::order_key(x[i]) - kmin) >> shift];
      const size_t r_lo = n % 2 ? h : h - 1;
      size_t cum = 0, b_lo = 0, b_hi = 0, before_lo = 0;
      bool found_lo = false;
      for (size_t b = 0; b < 4096; ++b) {
        const size_t c = hist[b];
        if (!found_lo && cum + c > r_lo) {
          b_lo = b;
          before_lo = cum;
          found_lo = true;
        }
        if (cum + c > h) {
          b_hi = b;
          break;
        }
        cum += c;
      }
      size_t w = 0;  // compact the elements of bins [b_lo, b_hi] to the front
      for (size_t i = 0; i < n; ++i) {
        const size_t b = (detail::order_key(x[i]) - kmin) >> shift;
        if (b >= b_lo && b <= b_hi) x[w++] = x[i];
      }
      const size_t k = h - before_lo;
      std::nth_element(x, x + k, x + w);
      if (n % 2) return x[k];
      return (*std::max_element(x, x + k) + x[k]) / T(2);
    }
  }
  std::nth_element(x, x + h, x + n);
  if (n % 2) return x[h];
  T lo = *std::max_element(x, x + h);
  return (lo + x[h]) / T(2);
}
// np.median(x) of a const array, using buf as the working copy (no allocation once grown)
template <class T>
T median_buf(const T* x, size_t n, std::vector<T>& buf) {
  buf.assign(x, x + n);
  return median_inplace(buf.data(), n);
}
template <class T> T median(std::vector<T> x) { return median_inplace(x.data(), x.size()); }
template <class T> T median(const T* x, size_t n) { std::vector<T> v(x, x + n); return median_inplace(v.data(), n); }

// np.percentile(x, q, method="linear") on a scratch copy (modified). Arithmetic in the
// dtype T as numpy does for float32/float64 inputs (q/100 and (n-1)*q in T; _lerp).
template <class T>
T percentile_inplace(T* x, size_t n, double q_percent) {
  static_assert(std::is_floating_point_v<T>, "float input");
  if (n == 0) return std::numeric_limits<T>::quiet_NaN();
  for (size_t i = 0; i < n; ++i)
    if (std::isnan(x[i])) return std::numeric_limits<T>::quiet_NaN();
  const T q = static_cast<T>(q_percent) / static_cast<T>(100);
  const T virt = static_cast<T>(n - 1) * q;
  if (!(virt < static_cast<T>(n - 1))) return *std::max_element(x, x + n);
  T prev_f = std::floor(virt);
  if (virt < T(0)) prev_f = 0;
  const size_t prev = static_cast<size_t>(prev_f);
  std::nth_element(x, x + prev, x + n);
  const T a = x[prev];
  const T b = *std::min_element(x + prev + 1, x + n);
  const T g = virt - prev_f;
  const T d = b - a;
  if (g >= T(0.5)) return b - d * (T(1) - g);
  return a + d * g;
}
template <class T> T percentile(std::vector<T> x, double q) { return percentile_inplace(x.data(), x.size(), q); }
template <class T> T percentile(const T* x, size_t n, double q) { std::vector<T> v(x, x + n); return percentile_inplace(v.data(), n, q); }
// np.quantile(x, q) with q in [0, 1] given in the dtype (q is not divided by 100).
template <class T>
T quantile(std::vector<T> x, double q) { return percentile_inplace(x.data(), x.size(), q * 100.0); }

// np.percentile(M, q, axis=0) -> (cols,) and axis=1 -> (rows,)
template <class T>
std::vector<T> percentile_axis0(const Mat<T>& m, double q) {
  std::vector<T> out(m.cols), col(m.rows);
  for (size_t c = 0; c < m.cols; ++c) {
    for (size_t r = 0; r < m.rows; ++r) col[r] = m(r, c);
    out[c] = percentile_inplace(col.data(), m.rows, q);
  }
  return out;
}
template <class T>
std::vector<T> percentile_axis1(const Mat<T>& m, double q) {
  std::vector<T> out(m.rows), row(m.cols);
  for (size_t r = 0; r < m.rows; ++r) {
    std::copy(m.row(r), m.row(r) + m.cols, row.begin());
    out[r] = percentile_inplace(row.data(), m.cols, q);
  }
  return out;
}
template <class T>
std::vector<T> median_axis0(const Mat<T>& m) {
  std::vector<T> out(m.cols), col(m.rows);
  for (size_t c = 0; c < m.cols; ++c) {
    for (size_t r = 0; r < m.rows; ++r) col[r] = m(r, c);
    out[c] = median_inplace(col.data(), m.rows);
  }
  return out;
}
template <class T>
std::vector<T> median_axis1(const Mat<T>& m) {
  std::vector<T> out(m.rows), row(m.cols);
  for (size_t r = 0; r < m.rows; ++r) {
    std::copy(m.row(r), m.row(r) + m.cols, row.begin());
    out[r] = median_inplace(row.data(), m.cols);
  }
  return out;
}

// np.argmax / np.argmin: first occurrence; the first NaN wins (numpy semantics).
template <class T>
size_t argmax(const T* x, size_t n) {
  if (n == 0) throw std::invalid_argument("argmax of an empty sequence");
  size_t best = 0;
  if (std::isnan(x[0])) return 0;
  for (size_t i = 1; i < n; ++i) {
    if (std::isnan(x[i])) return i;
    if (x[i] > x[best]) best = i;
  }
  return best;
}
template <class T> size_t argmax(const std::vector<T>& x) { return argmax(x.data(), x.size()); }
template <class T>
size_t argmin(const T* x, size_t n) {
  if (n == 0) throw std::invalid_argument("argmin of an empty sequence");
  size_t best = 0;
  if (std::isnan(x[0])) return 0;
  for (size_t i = 1; i < n; ++i) {
    if (std::isnan(x[i])) return i;
    if (x[i] < x[best]) best = i;
  }
  return best;
}
template <class T> size_t argmin(const std::vector<T>& x) { return argmin(x.data(), x.size()); }

// np.argsort(x, kind="stable"): ties keep input order.
template <class T>
std::vector<int64_t> argsort_stable(const std::vector<T>& x) {
  std::vector<int64_t> idx(x.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(), [&](int64_t a, int64_t b) { return np_less(x[a], x[b]); });
  return idx;
}

// np.argsort(x) default kind (numpy's introsort `aquicksort`: median-of-3 quicksort,
// insertion sort below 16 elements, heapsort past the depth limit). Ties come out in
// the same order as numpy (for arrays of <= 16 elements that is input order).
template <class T>
std::vector<int64_t> argsort_quick(const std::vector<T>& x);

// np.searchsorted(a, v, side) on a sorted a
template <class T>
size_t searchsorted(const std::vector<T>& a, T v, bool right = false) {
  auto it = right ? std::upper_bound(a.begin(), a.end(), v, [](T u, T w) { return np_less(u, w); })
                  : std::lower_bound(a.begin(), a.end(), v, [](T u, T w) { return np_less(u, w); });
  return static_cast<size_t>(it - a.begin());
}

// ------------------------------------------------------------------ ranges
// np.arange(start, stop, step) in dtype T: length ceil((stop-start)/step) in T; values
// start + i*delta with delta = (start+step) - start (numpy's fill).
template <class T>
std::vector<T> arange(T start, T stop, T step = T(1)) {
  const T len_f = (stop - start) / step;
  const double c = std::ceil(static_cast<double>(len_f));
  const size_t n = (std::isfinite(c) && c > 0) ? static_cast<size_t>(c) : 0;
  std::vector<T> y(n);
  if (n == 0) return y;
  y[0] = start;
  if (n == 1) return y;
  y[1] = start + step;
  const T delta = y[1] - y[0];
  for (size_t i = 2; i < n; ++i) y[i] = start + static_cast<T>(i) * delta;
  return y;
}

// np.linspace(start, stop, num, endpoint) in dtype T: i*step + start, last = stop.
template <class T>
std::vector<T> linspace(T start, T stop, size_t num, bool endpoint = true) {
  std::vector<T> y(num);
  if (num == 0) return y;
  const size_t div = endpoint ? num - 1 : num;
  const T delta = stop - start;
  if (div == 0) { y[0] = start; return y; }
  const T step = delta / static_cast<T>(div);
  if (step == T(0)) {
    for (size_t i = 0; i < num; ++i) y[i] = (static_cast<T>(i) / static_cast<T>(div)) * delta + start;
  } else {
    for (size_t i = 0; i < num; ++i) y[i] = static_cast<T>(i) * step + start;
  }
  if (endpoint && num > 1) y[num - 1] = stop;
  return y;
}

// np.hanning(M) (symmetric): 0.5 + 0.5*cos(pi*n/(M-1)), n = 1-M, 3-M, ..., M-1
dvec hanning(size_t M);
// np.hamming(M) (symmetric)
dvec hamming(size_t M);

// ------------------------------------------------------------------ interpolation
// np.interp(x, xp, fp, left=fp[0], right=fp[-1]) (float64; xp increasing).
dvec interp(const dvec& x, const dvec& xp, const dvec& fp,
            std::optional<double> left = std::nullopt, std::optional<double> right = std::nullopt);
double interp(double x, const dvec& xp, const dvec& fp,
              std::optional<double> left = std::nullopt, std::optional<double> right = std::nullopt);

// ------------------------------------------------------------------ histogram / bincount
template <class T>
struct Histogram {
  std::vector<int64_t> counts;
  std::vector<T> edges;
};
// np.histogram(v, bins=nbins, range=(lo, hi)): numpy's uniform-bin fast path, in dtype T
// (float32 data with float32 range values stays float32; pass T=double for float64).
template <class T>
Histogram<T> histogram(const std::vector<T>& v, size_t nbins, T lo, T hi);
// np.histogram(v, bins=edges): bins [e_i, e_{i+1}), the last bin closed.
template <class T, class E>
Histogram<E> histogram(const std::vector<T>& v, const std::vector<E>& edges);

// np.bincount(idx, minlength) / np.bincount(idx, weights, minlength) (weights -> float64)
std::vector<int64_t> bincount(const std::vector<int64_t>& idx, size_t minlength = 0);
dvec bincount(const std::vector<int64_t>& idx, const dvec& weights, size_t minlength = 0);

// ------------------------------------------------------------------ polynomials
// np.polyfit(x, y, deg) (least squares on the column-scaled Vandermonde matrix, SVD with
// rcond = len(x)*eps as numpy). Highest power first. T = double (float32 inputs: numpy
// solves in float64 and rounds the result; use polyfit<float>).
template <class T>
std::vector<T> polyfit(const std::vector<T>& x, const std::vector<T>& y, int deg);
// np.polyval(p, x): Horner, y = y*x + p[i]
template <class T>
T polyval(const std::vector<T>& p, T x) {
  T y = T(0);
  for (const T& c : p) y = y * x + c;
  return y;
}
template <class T>
std::vector<T> polyval(const std::vector<T>& p, const std::vector<T>& x) {
  std::vector<T> y(x.size());
  for (size_t i = 0; i < x.size(); ++i) y[i] = polyval(p, x[i]);
  return y;
}

// np.unwrap(p) (period 2*pi, discont pi), in dtype T
template <class T>
std::vector<T> unwrap(const std::vector<T>& p);

// ------------------------------------------------------------------ convolution
enum class ConvMode { Full, Same, Valid };
ConvMode conv_mode(const std::string& s);

// np.correlate(a, v, mode) (default mode "valid"): c[k] = sum_n a[n+k] * conj(v[n]).
template <class T>
std::vector<T> correlate(const std::vector<T>& a, const std::vector<T>& v, ConvMode mode = ConvMode::Valid);
// np.convolve(a, v, mode) (default "full"); no conjugation.
template <class T>
std::vector<T> convolve(const std::vector<T>& a, const std::vector<T>& v, ConvMode mode = ConvMode::Full);

// ------------------------------------------------------------------ linear algebra
// np.vdot(a, b) = sum conj(a) * b ; np.dot for real
template <class T> T vdot(const std::vector<T>& a, const std::vector<T>& b);
// np.linalg.norm(x) (2-norm of a vector): sqrt(x.x) (complex: re.re + im.im)
template <class T> real_t<T> norm(const T* x, size_t n);
template <class T> real_t<T> norm(const std::vector<T>& x) { return norm(x.data(), x.size()); }

// np.linalg.solve(A, b) / inv(A) (LU with partial pivoting as LAPACK getrf/getrs);
// T = double or complex<double>. Throws std::runtime_error("Singular matrix").
template <class T> std::vector<T> solve(const Mat<T>& A, const std::vector<T>& b);
template <class T> Mat<T> solve(const Mat<T>& A, const Mat<T>& B);
template <class T> Mat<T> inv(const Mat<T>& A);

// np.linalg.lstsq(A, b, rcond) (real A): SVD solution; singular values <= rcond*s_max are
// treated as zero (rcond < 0 -> numpy's default eps*max(M, N)). residuals: sum of
// squared residuals when rank == N < M (else empty), like numpy.
struct Lstsq {
  dvec x;
  dvec residuals;
  int rank = 0;
  dvec s;  // singular values, descending
};
Lstsq lstsq(const dmat& A, const dvec& b, double rcond = -1.0);
// complex right-hand side with a real A (solved per real/imag part)
cdvec lstsq_x(const dmat& A, const cdvec& b, double rcond = -1.0);

// matrix product helpers for the small systems the receivers build
template <class T>
Mat<T> matmul(const Mat<T>& A, const Mat<T>& B) {
  Mat<T> C(A.rows, B.cols, T(0));
  for (size_t i = 0; i < A.rows; ++i)
    for (size_t k = 0; k < A.cols; ++k) {
      const T a = A(i, k);
      for (size_t j = 0; j < B.cols; ++j) C(i, j) += a * B(k, j);
    }
  return C;
}

}  // namespace si::np
