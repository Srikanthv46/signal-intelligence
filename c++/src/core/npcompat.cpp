#include "core/npcompat.hpp"

#include <cfloat>

namespace si::np {

// ------------------------------------------------------------------ complex sqrt
// C99 csqrt, FreeBSD msun s_csqrt.c (Algorithm 312, CACM vol 10, Oct 1967), with its
// overflow scaling; numpy's npy_csqrt calls the platform csqrt, which agrees to rounding.
template <class R>
std::complex<R> csqrt(const std::complex<R>& z) {
  const R a = z.real(), b = z.imag();
  const R inf = std::numeric_limits<R>::infinity();
  if (a == 0 && b == 0) return {R(0), b};
  if (std::isinf(b)) return {inf, b};
  if (std::isnan(a)) { R t = (b - b) / (b - b); return {a, t}; }
  if (std::isinf(a)) {
    if (std::signbit(a)) return {std::fabs(b - b), std::copysign(a, b)};
    return {a, std::copysign(b - b, b)};
  }
  const R THRESH = std::numeric_limits<R>::max() / R(2.414213562373095048801688724209698);
  R ra = a, rb = b;
  bool scale = false;
  if (std::fabs(ra) >= THRESH || std::fabs(rb) >= THRESH) { ra *= R(0.25); rb *= R(0.25); scale = true; }
  std::complex<R> res;
  if (ra >= 0) {
    R t = std::sqrt((ra + std::hypot(ra, rb)) * R(0.5));
    res = {t, rb / (2 * t)};
  } else {
    R t = std::sqrt((-ra + std::hypot(ra, rb)) * R(0.5));
    res = {std::fabs(rb) / (2 * t), std::copysign(t, rb)};
  }
  if (scale) res *= R(2);
  return res;
}
template std::complex<float> csqrt(const std::complex<float>&);
template std::complex<double> csqrt(const std::complex<double>&);

// ------------------------------------------------------------------ argsort (numpy introsort)
namespace {
template <class T>
void aheapsort(const T* v, int64_t* tosort, int64_t n) {
  int64_t* a = tosort - 1;  // 1-based heap indexing as numpy
  int64_t i, j, l, tmp;
  for (l = n >> 1; l > 0; --l) {
    tmp = a[l];
    for (i = l, j = l << 1; j <= n;) {
      if (j < n && np_less(v[a[j]], v[a[j + 1]])) j += 1;
      if (np_less(v[tmp], v[a[j]])) { a[i] = a[j]; i = j; j += j; }
      else break;
    }
    a[i] = tmp;
  }
  for (; n > 1;) {
    tmp = a[n];
    a[n] = a[1];
    n -= 1;
    for (i = 1, j = 2; j <= n;) {
      if (j < n && np_less(v[a[j]], v[a[j + 1]])) j++;
      if (np_less(v[tmp], v[a[j]])) { a[i] = a[j]; i = j; j += j; }
      else break;
    }
    a[i] = tmp;
  }
}
}  // namespace

template <class T>
std::vector<int64_t> argsort_quick(const std::vector<T>& x) {
  const int64_t num = static_cast<int64_t>(x.size());
  std::vector<int64_t> idx(x.size());
  std::iota(idx.begin(), idx.end(), 0);
  if (num < 2) return idx;
  const T* v = x.data();
  constexpr int SMALL_QUICKSORT = 16;
  constexpr int PYA_QS_STACK = 128;
  int64_t* pl = idx.data();
  int64_t* pr = idx.data() + num - 1;
  int64_t* stack[PYA_QS_STACK];
  int64_t** sptr = stack;
  int depth[PYA_QS_STACK];
  int* psdepth = depth;
  int cdepth = 0;
  for (uint64_t u = static_cast<uint64_t>(num); u >>= 1;) ++cdepth;
  cdepth *= 2;
  int64_t *pm, *pi, *pj, *pk, vi;
  T vp;
  for (;;) {
    if (cdepth < 0) {
      aheapsort(v, pl, pr - pl + 1);
      goto stack_pop;
    }
    while ((pr - pl) > SMALL_QUICKSORT) {
      pm = pl + ((pr - pl) >> 1);
      if (np_less(v[*pm], v[*pl])) std::swap(*pm, *pl);
      if (np_less(v[*pr], v[*pm])) std::swap(*pr, *pm);
      if (np_less(v[*pm], v[*pl])) std::swap(*pm, *pl);
      vp = v[*pm];
      pi = pl;
      pj = pr - 1;
      std::swap(*pm, *pj);
      for (;;) {
        do { ++pi; } while (np_less(v[*pi], vp));
        do { --pj; } while (np_less(vp, v[*pj]));
        if (pi >= pj) break;
        std::swap(*pi, *pj);
      }
      pk = pr - 1;
      std::swap(*pi, *pk);
      if (pi - pl < pr - pi) {
        *sptr++ = pi + 1;
        *sptr++ = pr;
        pr = pi - 1;
      } else {
        *sptr++ = pl;
        *sptr++ = pi - 1;
        pl = pi + 1;
      }
      *psdepth++ = --cdepth;
    }
    for (pi = pl + 1; pi <= pr; ++pi) {  // insertion sort
      vi = *pi;
      vp = v[vi];
      pj = pi;
      pk = pi - 1;
      while (pj > pl && np_less(vp, v[*pk])) *pj-- = *pk--;
      *pj = vi;
    }
  stack_pop:
    if (sptr == stack) break;
    pr = *(--sptr);
    pl = *(--sptr);
    cdepth = *(--psdepth);
  }
  return idx;
}
template std::vector<int64_t> argsort_quick(const std::vector<float>&);
template std::vector<int64_t> argsort_quick(const std::vector<double>&);
template std::vector<int64_t> argsort_quick(const std::vector<int64_t>&);
template std::vector<int64_t> argsort_quick(const std::vector<int32_t>&);

// ------------------------------------------------------------------ windows
dvec hanning(size_t M) {
  if (M == 0) return {};
  if (M == 1) return {1.0};
  dvec w(M);
  for (size_t i = 0; i < M; ++i) {
    const double n = static_cast<double>(static_cast<int64_t>(2 * i) + 1 - static_cast<int64_t>(M));
    w[i] = 0.5 + 0.5 * std::cos(kPi * n / static_cast<double>(M - 1));
  }
  return w;
}

dvec hamming(size_t M) {
  if (M == 0) return {};
  if (M == 1) return {1.0};
  dvec w(M);
  for (size_t i = 0; i < M; ++i) {
    const double n = static_cast<double>(static_cast<int64_t>(2 * i) + 1 - static_cast<int64_t>(M));
    w[i] = 0.54 + 0.46 * std::cos(kPi * n / static_cast<double>(M - 1));
  }
  return w;
}

// ------------------------------------------------------------------ interp
namespace {
// numpy compiled interp: j = -1 below, n above, else largest j with xp[j] <= x
inline int64_t interp_index(double x, const dvec& xp) {
  const int64_t n = static_cast<int64_t>(xp.size());
  if (x < xp[0]) return -1;
  if (x > xp[n - 1]) return n;
  int64_t imin = 0, imax = n;
  while (imin < imax) {
    const int64_t imid = imin + ((imax - imin) >> 1);
    if (x >= xp[imid]) imin = imid + 1;
    else imax = imid;
  }
  return imin - 1;
}
inline double interp_one(double xv, const dvec& xp, const dvec& fp, double lval, double rval) {
  if (std::isnan(xv)) return xv;
  const int64_t n = static_cast<int64_t>(xp.size());
  const int64_t j = interp_index(xv, xp);
  if (j == -1) return lval;
  if (j == n) return rval;
  if (j == n - 1) return fp[j];
  if (xp[j] == xv) return fp[j];  // avoid non-finite interpolation
  const double slope = (fp[j + 1] - fp[j]) / (xp[j + 1] - xp[j]);
  double r = slope * (xv - xp[j]) + fp[j];
  if (std::isnan(r)) {
    r = slope * (xv - xp[j + 1]) + fp[j + 1];
    if (std::isnan(r) && fp[j] == fp[j + 1]) r = fp[j];
  }
  return r;
}
}  // namespace

dvec interp(const dvec& x, const dvec& xp, const dvec& fp, std::optional<double> left,
            std::optional<double> right) {
  if (xp.empty() || xp.size() != fp.size()) throw std::invalid_argument("interp: xp and fp must be non-empty, same length");
  const double lval = left ? *left : fp.front(), rval = right ? *right : fp.back();
  dvec y(x.size());
  if (xp.size() == 1) {
    for (size_t i = 0; i < x.size(); ++i)
      y[i] = std::isnan(x[i]) ? x[i] : (x[i] < xp[0] ? lval : (x[i] > xp[0] ? rval : fp[0]));
    return y;
  }
  for (size_t i = 0; i < x.size(); ++i) y[i] = interp_one(x[i], xp, fp, lval, rval);
  return y;
}

double interp(double x, const dvec& xp, const dvec& fp, std::optional<double> left, std::optional<double> right) {
  return interp(dvec{x}, xp, fp, left, right)[0];
}

// ------------------------------------------------------------------ histogram
template <class T>
Histogram<T> histogram(const std::vector<T>& v, size_t nbins, T lo, T hi) {
  if (nbins == 0) throw std::invalid_argument("histogram: bins must be positive");
  if (lo > hi) throw std::invalid_argument("max must be larger than min in range parameter.");
  T first = lo, last = hi;
  if (first == last) { first = first - T(0.5); last = last + T(0.5); }
  Histogram<T> h;
  h.edges = linspace<T>(first, last, nbins + 1, true);
  h.counts.assign(nbins, 0);
  const T denom = last - first;
  const T numer = static_cast<T>(nbins);
  const int64_t nb = static_cast<int64_t>(nbins);
  for (const T& a : v) {
    if (!(a >= first) || !(a <= last)) continue;
    const T f = ((a - first) / denom) * numer;
    int64_t idx = static_cast<int64_t>(f);
    if (idx == nb) idx -= 1;
    if (a < h.edges[idx]) idx -= 1;
    if (a >= h.edges[idx + 1] && idx != nb - 1) idx += 1;
    h.counts[idx] += 1;
  }
  return h;
}
template Histogram<float> histogram(const std::vector<float>&, size_t, float, float);
template Histogram<double> histogram(const std::vector<double>&, size_t, double, double);

template <class T, class E>
Histogram<E> histogram(const std::vector<T>& v, const std::vector<E>& edges) {
  if (edges.size() < 2) throw std::invalid_argument("histogram: need at least 2 edges");
  for (size_t i = 0; i + 1 < edges.size(); ++i)
    if (edges[i] > edges[i + 1]) throw std::invalid_argument("bins must increase monotonically");
  dvec s(v.begin(), v.end());
  std::sort(s.begin(), s.end(), [](double a, double b) { return np_less(a, b); });
  const size_t ne = edges.size();
  std::vector<int64_t> cum(ne);
  for (size_t i = 0; i < ne; ++i) {
    const double e = static_cast<double>(edges[i]);
    auto it = (i + 1 < ne) ? std::lower_bound(s.begin(), s.end(), e, [](double a, double b) { return np_less(a, b); })
                           : std::upper_bound(s.begin(), s.end(), e, [](double a, double b) { return np_less(a, b); });
    cum[i] = it - s.begin();
  }
  Histogram<E> h;
  h.edges = edges;
  h.counts.resize(ne - 1);
  for (size_t i = 0; i + 1 < ne; ++i) h.counts[i] = cum[i + 1] - cum[i];
  return h;
}
template Histogram<double> histogram(const std::vector<double>&, const std::vector<double>&);
template Histogram<float> histogram(const std::vector<float>&, const std::vector<float>&);
template Histogram<double> histogram(const std::vector<float>&, const std::vector<double>&);

// ------------------------------------------------------------------ bincount
std::vector<int64_t> bincount(const std::vector<int64_t>& idx, size_t minlength) {
  int64_t mx = -1;
  for (int64_t i : idx) {
    if (i < 0) throw std::invalid_argument("'list' argument must have no negative elements");
    mx = std::max(mx, i);
  }
  std::vector<int64_t> out(std::max<size_t>(minlength, static_cast<size_t>(mx + 1)), 0);
  for (int64_t i : idx) out[i] += 1;
  return out;
}

dvec bincount(const std::vector<int64_t>& idx, const dvec& w, size_t minlength) {
  if (w.size() != idx.size()) throw std::invalid_argument("The weights and list don't have the same length.");
  int64_t mx = -1;
  for (int64_t i : idx) {
    if (i < 0) throw std::invalid_argument("'list' argument must have no negative elements");
    mx = std::max(mx, i);
  }
  dvec out(std::max<size_t>(minlength, static_cast<size_t>(mx + 1)), 0.0);
  for (size_t k = 0; k < idx.size(); ++k) out[idx[k]] += w[k];
  return out;
}

// ------------------------------------------------------------------ SVD / lstsq
namespace {
// One-sided Jacobi SVD of A (m x n, m >= n): A = U diag(s) V^T, s descending.
void svd_jacobi(const dmat& A, dmat& U, dvec& s, dmat& V) {
  const size_t m = A.rows, n = A.cols;
  U = A;
  V = dmat(n, n, 0.0);
  for (size_t i = 0; i < n; ++i) V(i, i) = 1.0;
  const double eps = std::numeric_limits<double>::epsilon();
  for (int sweep = 0; sweep < 80; ++sweep) {
    bool rotated = false;
    for (size_t p = 0; p + 1 < n; ++p)
      for (size_t q = p + 1; q < n; ++q) {
        double alpha = 0, beta = 0, gamma = 0;
        for (size_t i = 0; i < m; ++i) {
          const double up = U(i, p), uq = U(i, q);
          alpha += up * up; beta += uq * uq; gamma += up * uq;
        }
        if (gamma == 0.0 || std::fabs(gamma) <= eps * std::sqrt(alpha * beta)) continue;
        rotated = true;
        const double zeta = (beta - alpha) / (2.0 * gamma);
        const double t = std::copysign(1.0, zeta) / (std::fabs(zeta) + std::sqrt(1.0 + zeta * zeta));
        const double c = 1.0 / std::sqrt(1.0 + t * t), sn = c * t;
        for (size_t i = 0; i < m; ++i) {
          const double up = U(i, p), uq = U(i, q);
          U(i, p) = c * up - sn * uq;
          U(i, q) = sn * up + c * uq;
        }
        for (size_t i = 0; i < n; ++i) {
          const double vp = V(i, p), vq = V(i, q);
          V(i, p) = c * vp - sn * vq;
          V(i, q) = sn * vp + c * vq;
        }
      }
    if (!rotated) break;
  }
  s.assign(n, 0.0);
  for (size_t j = 0; j < n; ++j) {
    double nn = 0;
    for (size_t i = 0; i < m; ++i) nn += U(i, j) * U(i, j);
    s[j] = std::sqrt(nn);
    if (s[j] > 0)
      for (size_t i = 0; i < m; ++i) U(i, j) /= s[j];
  }
  // sort descending
  std::vector<size_t> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return s[a] > s[b]; });
  dmat U2(m, n), V2(n, n);
  dvec s2(n);
  for (size_t j = 0; j < n; ++j) {
    s2[j] = s[ord[j]];
    for (size_t i = 0; i < m; ++i) U2(i, j) = U(i, ord[j]);
    for (size_t i = 0; i < n; ++i) V2(i, j) = V(i, ord[j]);
  }
  U = std::move(U2); V = std::move(V2); s = std::move(s2);
}

dmat transpose(const dmat& A) {
  dmat T(A.cols, A.rows);
  for (size_t i = 0; i < A.rows; ++i)
    for (size_t j = 0; j < A.cols; ++j) T(j, i) = A(i, j);
  return T;
}
}  // namespace

Lstsq lstsq(const dmat& A, const dvec& b, double rcond) {
  const size_t m = A.rows, n = A.cols;
  if (b.size() != m) throw std::invalid_argument("lstsq: incompatible dimensions");
  Lstsq r;
  r.x.assign(n, 0.0);
  if (m == 0 || n == 0) return r;
  if (rcond < 0) rcond = std::numeric_limits<double>::epsilon() * static_cast<double>(std::max(m, n));
  dmat U, V;
  dvec s;
  if (m >= n) {
    svd_jacobi(A, U, s, V);
  } else {  // A^T = U' S V'^T  =>  A = V' S U'^T
    dmat Up, Vp;
    svd_jacobi(transpose(A), Up, s, Vp);
    U = Vp; V = Up;
  }
  const size_t k = std::min(m, n);
  const double thr = rcond * (s.empty() ? 0.0 : s[0]);
  for (size_t j = 0; j < k; ++j) {
    if (!(s[j] > thr)) continue;
    r.rank += 1;
    double ub = 0;
    for (size_t i = 0; i < m; ++i) ub += U(i, j) * b[i];
    const double c = ub / s[j];
    for (size_t i = 0; i < n; ++i) r.x[i] += c * V(i, j);
  }
  r.s = s;
  r.s.resize(k);
  if (r.rank == static_cast<int>(n) && m > n) {
    double ss = 0;
    for (size_t i = 0; i < m; ++i) {
      double ax = 0;
      for (size_t j = 0; j < n; ++j) ax += A(i, j) * r.x[j];
      const double d = b[i] - ax;
      ss += d * d;
    }
    r.residuals = {ss};
  }
  return r;
}

cdvec lstsq_x(const dmat& A, const cdvec& b, double rcond) {
  dvec br(b.size()), bi(b.size());
  for (size_t i = 0; i < b.size(); ++i) { br[i] = b[i].real(); bi[i] = b[i].imag(); }
  const dvec xr = lstsq(A, br, rcond).x, xi = lstsq(A, bi, rcond).x;
  cdvec x(xr.size());
  for (size_t i = 0; i < x.size(); ++i) x[i] = {xr[i], xi[i]};
  return x;
}

// ------------------------------------------------------------------ polyfit
template <class T>
std::vector<T> polyfit(const std::vector<T>& x, const std::vector<T>& y, int deg) {
  if (deg < 0) throw std::invalid_argument("expected deg >= 0");
  if (x.empty()) throw std::invalid_argument("expected non-empty vector for x");
  if (x.size() != y.size()) throw std::invalid_argument("expected x and y to have same length");
  const size_t n = x.size(), order = static_cast<size_t>(deg) + 1;
  const T rcond = static_cast<T>(n) * std::numeric_limits<T>::epsilon();
  Mat<T> lhs(n, order);
  for (size_t i = 0; i < n; ++i) {  // np.vander: powers by multiply.accumulate
    T pw = T(1);
    lhs(i, order - 1) = T(1);
    for (size_t k = 1; k < order; ++k) {
      pw = (k == 1) ? x[i] : pw * x[i];
      lhs(i, order - 1 - k) = pw;
    }
  }
  std::vector<T> scale(order, T(0));  // sqrt((lhs*lhs).sum(axis=0)), rows added in order
  for (size_t i = 0; i < n; ++i)
    for (size_t k = 0; k < order; ++k) scale[k] += lhs(i, k) * lhs(i, k);
  for (auto& s : scale) s = std::sqrt(s);
  dmat A(n, order);
  dvec b(n);
  for (size_t i = 0; i < n; ++i) {
    for (size_t k = 0; k < order; ++k) A(i, k) = static_cast<double>(lhs(i, k) / scale[k]);
    b[i] = static_cast<double>(y[i]);
  }
  const dvec c = lstsq(A, b, static_cast<double>(rcond)).x;
  std::vector<T> out(order);
  for (size_t k = 0; k < order; ++k) out[k] = static_cast<T>(c[k]) / scale[k];
  return out;
}
template std::vector<double> polyfit(const std::vector<double>&, const std::vector<double>&, int);
template std::vector<float> polyfit(const std::vector<float>&, const std::vector<float>&, int);

// ------------------------------------------------------------------ unwrap
template <class T>
std::vector<T> unwrap(const std::vector<T>& p) {
  std::vector<T> up(p);
  if (p.size() < 2) return up;
  const T period = static_cast<T>(2 * kPi);
  const T hi = static_cast<T>(kPi), lo = -hi;
  const T discont = static_cast<T>(kPi);
  T acc = T(0);
  for (size_t i = 1; i < p.size(); ++i) {
    const T dd = p[i] - p[i - 1];
    T ddmod = mod<T>(dd - lo, period) + lo;
    if (ddmod == lo && dd > T(0)) ddmod = hi;
    T corr = ddmod - dd;
    if (std::fabs(dd) < discont) corr = T(0);
    acc = (i == 1) ? corr : acc + corr;
    up[i] = p[i] + acc;
  }
  return up;
}
template std::vector<double> unwrap(const std::vector<double>&);
template std::vector<float> unwrap(const std::vector<float>&);

// ------------------------------------------------------------------ convolve / correlate
ConvMode conv_mode(const std::string& s) {
  if (s == "full") return ConvMode::Full;
  if (s == "same") return ConvMode::Same;
  if (s == "valid") return ConvMode::Valid;
  throw std::invalid_argument("mode must be one of 'full', 'same', 'valid'");
}

namespace {
template <class T>
inline T dot_n(const T* a, const T* b, size_t n) {
  T s = T(0);
  if constexpr (is_cplx<T>::value) {
    for (size_t i = 0; i < n; ++i) s += cmul(a[i], b[i]);
  } else {
    for (size_t i = 0; i < n; ++i) s += a[i] * b[i];
  }
  return s;
}

// numpy's _pyarray_correlate for n1 >= n2 (the caller swaps)
template <class T>
std::vector<T> pyarray_correlate(const T* ap1, size_t n1, const T* ap2, size_t n2, ConvMode mode) {
  size_t length = n1, n = n2, n_left = 0, n_right = 0;
  switch (mode) {
    case ConvMode::Valid: length = length - n + 1; break;
    case ConvMode::Same: n_left = n / 2; n_right = n - n_left - 1; break;
    case ConvMode::Full: n_right = n - 1; n_left = n - 1; length = length + n - 1; break;
  }
  std::vector<T> out(length);
  size_t op = 0;
  const T* ip1 = ap1;
  const T* ip2 = ap2 + n_left;
  n = n - n_left;
  for (size_t i = 0; i < n_left; ++i) {
    out[op++] = dot_n(ip1, ip2, n);
    ++n;
    --ip2;
  }
  for (size_t i = 0; i < n1 - n2 + 1; ++i) {
    out[op++] = dot_n(ip1, ip2, n);
    ++ip1;
  }
  for (size_t i = 0; i < n_right; ++i) {
    --n;
    out[op++] = dot_n(ip1, ip2, n);
    ++ip1;
  }
  return out;
}
}  // namespace

template <class T>
std::vector<T> correlate(const std::vector<T>& a, const std::vector<T>& v, ConvMode mode) {
  if (a.empty()) throw std::invalid_argument("a cannot be empty");
  if (v.empty()) throw std::invalid_argument("v cannot be empty");
  std::vector<T> vc(v);
  if constexpr (is_cplx<T>::value)
    for (auto& z : vc) z = std::conj(z);
  if (a.size() >= vc.size()) return pyarray_correlate(a.data(), a.size(), vc.data(), vc.size(), mode);
  std::vector<T> r = pyarray_correlate(vc.data(), vc.size(), a.data(), a.size(), mode);
  std::reverse(r.begin(), r.end());
  return r;
}

template <class T>
std::vector<T> convolve(const std::vector<T>& a0, const std::vector<T>& v0, ConvMode mode) {
  if (a0.empty()) throw std::invalid_argument("a cannot be empty");
  if (v0.empty()) throw std::invalid_argument("v cannot be empty");
  const std::vector<T>& a = (v0.size() > a0.size()) ? v0 : a0;
  const std::vector<T>& v = (v0.size() > a0.size()) ? a0 : v0;
  std::vector<T> vr(v.rbegin(), v.rend());
  return pyarray_correlate(a.data(), a.size(), vr.data(), vr.size(), mode);
}

template std::vector<float> correlate(const std::vector<float>&, const std::vector<float>&, ConvMode);
template std::vector<double> correlate(const std::vector<double>&, const std::vector<double>&, ConvMode);
template std::vector<cf> correlate(const std::vector<cf>&, const std::vector<cf>&, ConvMode);
template std::vector<cd> correlate(const std::vector<cd>&, const std::vector<cd>&, ConvMode);
template std::vector<float> convolve(const std::vector<float>&, const std::vector<float>&, ConvMode);
template std::vector<double> convolve(const std::vector<double>&, const std::vector<double>&, ConvMode);
template std::vector<cf> convolve(const std::vector<cf>&, const std::vector<cf>&, ConvMode);
template std::vector<cd> convolve(const std::vector<cd>&, const std::vector<cd>&, ConvMode);

// ------------------------------------------------------------------ vdot / norm
template <class T>
T vdot(const std::vector<T>& a, const std::vector<T>& b) {
  if (a.size() != b.size()) throw std::invalid_argument("vdot: length mismatch");
  T s = T(0);
  for (size_t i = 0; i < a.size(); ++i) {
    if constexpr (is_cplx<T>::value) s += cmul(std::conj(a[i]), b[i]);
    else s += a[i] * b[i];
  }
  return s;
}
template float vdot(const std::vector<float>&, const std::vector<float>&);
template double vdot(const std::vector<double>&, const std::vector<double>&);
template cf vdot(const std::vector<cf>&, const std::vector<cf>&);
template cd vdot(const std::vector<cd>&, const std::vector<cd>&);

template <class T>
real_t<T> norm(const T* x, size_t n) {
  using R = real_t<T>;
  if constexpr (is_cplx<T>::value) {
    R a = 0, b = 0;
    for (size_t i = 0; i < n; ++i) { a += x[i].real() * x[i].real(); b += x[i].imag() * x[i].imag(); }
    return std::sqrt(a + b);
  } else {
    R a = 0;
    for (size_t i = 0; i < n; ++i) a += x[i] * x[i];
    return std::sqrt(a);
  }
}
template float norm(const float*, size_t);
template double norm(const double*, size_t);
template float norm(const cf*, size_t);
template double norm(const cd*, size_t);

// ------------------------------------------------------------------ LU solve / inv
namespace {
template <class T>
inline double cabs1(const T& v) {
  if constexpr (is_cplx<T>::value) return std::fabs(v.real()) + std::fabs(v.imag());
  else return std::fabs(v);
}

// LAPACK getrf (partial pivoting: first index of max |a| (cabs1 for complex), column
// scaled by the reciprocal pivot), in place; piv[k] = row swapped with k.
template <class T>
void lu_factor(Mat<T>& A, std::vector<size_t>& piv) {
  const size_t n = A.rows;
  piv.resize(n);
  for (size_t k = 0; k < n; ++k) {
    size_t p = k;
    double best = cabs1(A(k, k));
    for (size_t i = k + 1; i < n; ++i) {
      const double v = cabs1(A(i, k));
      if (v > best) { best = v; p = i; }
    }
    piv[k] = p;
    if (A(p, k) == T(0)) throw std::runtime_error("Singular matrix");
    if (p != k)
      for (size_t j = 0; j < n; ++j) std::swap(A(k, j), A(p, j));
    const T rcp = T(1) / A(k, k);
    for (size_t i = k + 1; i < n; ++i) A(i, k) *= rcp;
    for (size_t i = k + 1; i < n; ++i) {
      const T l = A(i, k);
      for (size_t j = k + 1; j < n; ++j) A(i, j) -= l * A(k, j);
    }
  }
}

template <class T>
void lu_solve_inplace(const Mat<T>& LU, const std::vector<size_t>& piv, T* b, size_t stride) {
  const size_t n = LU.rows;
  for (size_t k = 0; k < n; ++k)
    if (piv[k] != k) std::swap(b[k * stride], b[piv[k] * stride]);
  for (size_t i = 0; i < n; ++i)
    for (size_t k = 0; k < i; ++k) b[i * stride] -= LU(i, k) * b[k * stride];
  for (size_t ii = n; ii-- > 0;) {
    for (size_t k = ii + 1; k < n; ++k) b[ii * stride] -= LU(ii, k) * b[k * stride];
    b[ii * stride] /= LU(ii, ii);
  }
}
}  // namespace

template <class T>
std::vector<T> solve(const Mat<T>& A, const std::vector<T>& b) {
  if (A.rows != A.cols || b.size() != A.rows) throw std::invalid_argument("solve: shape mismatch");
  Mat<T> LU = A;
  std::vector<size_t> piv;
  lu_factor(LU, piv);
  std::vector<T> x = b;
  lu_solve_inplace(LU, piv, x.data(), 1);
  return x;
}

template <class T>
Mat<T> solve(const Mat<T>& A, const Mat<T>& B) {
  if (A.rows != A.cols || B.rows != A.rows) throw std::invalid_argument("solve: shape mismatch");
  Mat<T> LU = A;
  std::vector<size_t> piv;
  lu_factor(LU, piv);
  Mat<T> X = B;
  for (size_t c = 0; c < B.cols; ++c) lu_solve_inplace(LU, piv, X.data.data() + c, B.cols);
  return X;
}

template <class T>
Mat<T> inv(const Mat<T>& A) {
  Mat<T> I(A.rows, A.rows, T(0));
  for (size_t i = 0; i < A.rows; ++i) I(i, i) = T(1);
  return solve(A, I);
}

template dvec solve(const dmat&, const dvec&);
template cdvec solve(const Mat<cd>&, const cdvec&);
template dmat solve(const dmat&, const dmat&);
template Mat<cd> solve(const Mat<cd>&, const Mat<cd>&);
template dmat inv(const dmat&);
template Mat<cd> inv(const Mat<cd>&);

}  // namespace si::np
