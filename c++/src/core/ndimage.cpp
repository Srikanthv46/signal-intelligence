#include "core/ndimage.hpp"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <deque>
#include <stdexcept>

#include "core/npcompat.hpp"

namespace si::ndi {

Mode mode_from(const std::string& s) {
  if (s == "reflect" || s == "grid-mirror") return Mode::Reflect;
  if (s == "nearest") return Mode::Nearest;
  if (s == "mirror") return Mode::Mirror;
  if (s == "wrap" || s == "grid-wrap") return Mode::Wrap;
  if (s == "constant" || s == "grid-constant") return Mode::Constant;
  throw std::invalid_argument("boundary mode not supported: " + s);
}

int64_t ext_index(int64_t i, int64_t n, Mode m) {
  if (i >= 0 && i < n) return i;
  switch (m) {
    case Mode::Constant: return -1;
    case Mode::Nearest: return i < 0 ? 0 : n - 1;
    case Mode::Wrap: {
      int64_t k = i % n;
      return k < 0 ? k + n : k;
    }
    case Mode::Reflect: {
      const int64_t p = 2 * n;
      int64_t k = i % p;
      if (k < 0) k += p;
      return k >= n ? p - 1 - k : k;
    }
    case Mode::Mirror: {
      if (n == 1) return 0;
      const int64_t p = 2 * n - 2;
      int64_t k = i % p;
      if (k < 0) k += p;
      return k >= n ? p - k : k;
    }
  }
  return -1;
}

namespace {

void check_origin(size_t size, int origin) {
  const int64_t c = static_cast<int64_t>(size / 2) + origin;
  if (size == 0) throw std::invalid_argument("filter size must be >= 1");
  if (c < 0 || c >= static_cast<int64_t>(size)) throw std::invalid_argument("invalid origin");
}

// x extended by `before` samples on the left and `after` on the right, as type B
template <class B, class T>
void extend_line(const T* x, size_t n, size_t before, size_t after, Mode mode, B cval, std::vector<B>& buf) {
  buf.resize(n + before + after);
  const int64_t N = static_cast<int64_t>(n);
  for (size_t j = 0; j < buf.size(); ++j) {
    const int64_t i = static_cast<int64_t>(j) - static_cast<int64_t>(before);
    const int64_t k = ext_index(i, N, mode);
    buf[j] = k < 0 ? cval : static_cast<B>(x[k]);
  }
}

// ---------------------------------------------------------------- 1-D cores (contiguous)
template <class T>
void uniform_core(const T* x, size_t n, size_t size, Mode mode, double cval, int origin, T* out,
                  std::vector<double>& buf) {
  check_origin(size, origin);
  const size_t size1 = size / 2, size2 = size - size1 - 1;
  const size_t before = static_cast<size_t>(static_cast<int64_t>(size1) + origin);
  const size_t after = static_cast<size_t>(static_cast<int64_t>(size2) - origin);
  extend_line(x, n, before, after, mode, cval, buf);
  const double fsz = static_cast<double>(size);
  double tmp = 0.0;
  for (size_t ll = 0; ll < size; ++ll) tmp += buf[ll];
  tmp /= fsz;
  if (n == 0) return;
  out[0] = static_cast<T>(tmp);
  for (size_t ll = 1; ll < n; ++ll) {
    tmp += (buf[ll + size - 1] - buf[ll - 1]) / fsz;
    out[ll] = static_cast<T>(tmp);
  }
}

template <class T>
void minmax_core(const T* x, size_t n, size_t size, Mode mode, double cval, int origin, bool is_max, T* out,
                 std::vector<T>& buf) {
  check_origin(size, origin);
  const size_t size1 = size / 2, size2 = size - size1 - 1;
  const size_t before = static_cast<size_t>(static_cast<int64_t>(size1) + origin);
  const size_t after = static_cast<size_t>(static_cast<int64_t>(size2) - origin);
  extend_line(x, n, before, after, mode, static_cast<T>(cval), buf);
  std::deque<size_t> dq;  // indices into buf, values monotone
  auto better = [is_max](T a, T b) { return is_max ? (a >= b) : (a <= b); };
  for (size_t j = 0; j < buf.size(); ++j) {
    while (!dq.empty() && better(buf[j], buf[dq.back()])) dq.pop_back();
    dq.push_back(j);
    if (j + 1 >= size) {
      const size_t i = j + 1 - size;  // output index
      while (dq.front() < i) dq.pop_front();
      out[i] = buf[dq.front()];
    }
  }
}

template <class T>
void median_core(const T* x, size_t n, size_t size, Mode mode, double cval, int origin, T* out,
                 std::vector<T>& buf, std::vector<T>& win) {
  check_origin(size, origin);
  const size_t size1 = size / 2, size2 = size - size1 - 1;
  const size_t before = static_cast<size_t>(static_cast<int64_t>(size1) + origin);
  const size_t after = static_cast<size_t>(static_cast<int64_t>(size2) - origin);
  extend_line(x, n, before, after, mode, static_cast<T>(cval), buf);
  if (n == 0) return;
  const size_t rank = size / 2;
  win.assign(buf.begin(), buf.begin() + size);
  std::sort(win.begin(), win.end());
  out[0] = win[rank];
  for (size_t i = 1; i < n; ++i) {
    const T old = buf[i - 1], nw = buf[i + size - 1];
    if (!(old == nw)) {
      auto it = std::lower_bound(win.begin(), win.end(), old);
      if (it == win.end() || !(*it == old)) it = std::find(win.begin(), win.end(), old);  // NaN fallback
      win.erase(it);
      win.insert(std::upper_bound(win.begin(), win.end(), nw), nw);
    }
    out[i] = win[rank];
  }
}

void correlate_core(const double* x, size_t n, const dvec& w, Mode mode, double cval, int origin, double* out,
                    std::vector<double>& buf) {
  const size_t fs = w.size();
  check_origin(fs, origin);
  const int64_t size1 = static_cast<int64_t>(fs / 2), size2 = static_cast<int64_t>(fs) - size1 - 1;
  int symmetric = 0;
  if (fs & 1) {
    symmetric = 1;
    for (int64_t ii = 1; ii <= size1; ++ii)
      if (std::fabs(w[size1 + ii] - w[size1 - ii]) > DBL_EPSILON) { symmetric = 0; break; }
    if (symmetric == 0) {
      symmetric = -1;
      for (int64_t ii = 1; ii <= size1; ++ii)
        if (std::fabs(w[size1 + ii] + w[size1 - ii]) > DBL_EPSILON) { symmetric = 0; break; }
    }
  }
  const size_t before = static_cast<size_t>(size1 + origin);
  const size_t after = static_cast<size_t>(size2 - origin);
  extend_line(x, n, before, after, mode, cval, buf);
  const double* fw = w.data() + size1;
  for (size_t ll = 0; ll < n; ++ll) {
    const double* il = buf.data() + ll + size1;
    double o;
    if (symmetric > 0) {
      o = il[0] * fw[0];
      for (int64_t jj = -size1; jj < 0; ++jj) o += (il[jj] + il[-jj]) * fw[jj];
    } else if (symmetric < 0) {
      o = il[0] * fw[0];
      for (int64_t jj = -size1; jj < 0; ++jj) o += (il[jj] - il[-jj]) * fw[jj];
    } else {
      o = il[size2] * fw[size2];
      for (int64_t jj = -size1; jj < size2; ++jj) o += il[jj] * fw[jj];
    }
    out[ll] = o;
  }
}

// apply a contiguous-line function along `axis` of a matrix
template <class T, class F>
Mat<T> along_axis(const Mat<T>& x, int axis, F f) {
  if (axis != 0 && axis != 1) throw std::invalid_argument("axis must be 0 or 1");
  Mat<T> out(x.rows, x.cols);
  if (axis == 1) {
    for (size_t r = 0; r < x.rows; ++r) f(x.row(r), x.cols, out.row(r));
  } else {
    std::vector<T> line(x.rows), res(x.rows);
    for (size_t c = 0; c < x.cols; ++c) {
      for (size_t r = 0; r < x.rows; ++r) line[r] = x(r, c);
      f(line.data(), x.rows, res.data());
      for (size_t r = 0; r < x.rows; ++r) out(r, c) = res[r];
    }
  }
  return out;
}

}  // namespace

// ------------------------------------------------------------------ uniform
template <class T>
std::vector<T> uniform_filter1d(const std::vector<T>& x, size_t size, Mode mode, double cval, int origin) {
  std::vector<T> y(x.size());
  std::vector<double> buf;
  uniform_core(x.data(), x.size(), size, mode, cval, origin, y.data(), buf);
  return y;
}

template <class T>
Mat<T> uniform_filter1d(const Mat<T>& x, size_t size, int axis, Mode mode, double cval, int origin) {
  if (axis == 0 && x.rows > 0 && x.cols > 0) {
    // all columns at once, row by row (vectorises across columns, no strided copies); per
    // column the same running sum in the same order as uniform_core
    check_origin(size, origin);
    const size_t n = x.rows, C = x.cols;
    const int64_t before = static_cast<int64_t>(size / 2) + origin;
    const int64_t N = static_cast<int64_t>(n);
    const double fsz = static_cast<double>(size);
    Mat<T> out(n, C);
    std::vector<double> tmp(C, 0.0), ra(C), rb(C);
    auto fetch = [&](int64_t j, std::vector<double>& r) {  // row j of the extended column, as double
      const int64_t k = ext_index(j - before, N, mode);
      if (k < 0) {
        std::fill(r.begin(), r.end(), cval);
      } else {
        const T* xr = x.row(static_cast<size_t>(k));
        for (size_t c = 0; c < C; ++c) r[c] = static_cast<double>(xr[c]);
      }
    };
    for (size_t ll = 0; ll < size; ++ll) {
      fetch(static_cast<int64_t>(ll), ra);
      for (size_t c = 0; c < C; ++c) tmp[c] += ra[c];
    }
    T* o = out.row(0);
    for (size_t c = 0; c < C; ++c) {
      tmp[c] /= fsz;
      o[c] = static_cast<T>(tmp[c]);
    }
    for (size_t ll = 1; ll < n; ++ll) {
      fetch(static_cast<int64_t>(ll + size - 1), ra);
      fetch(static_cast<int64_t>(ll - 1), rb);
      o = out.row(ll);
      for (size_t c = 0; c < C; ++c) {
        tmp[c] += (ra[c] - rb[c]) / fsz;
        o[c] = static_cast<T>(tmp[c]);
      }
    }
    return out;
  }
  std::vector<double> buf;
  return along_axis(x, axis, [&](const T* in, size_t n, T* out) {
    uniform_core(in, n, size, mode, cval, origin, out, buf);
  });
}

template <class T>
Mat<T> uniform_filter(const Mat<T>& x, size_t s0, size_t s1, Mode mode, double cval) {
  Mat<T> y = x;
  if (s0 > 1) y = uniform_filter1d(y, s0, 0, mode, cval, 0);
  if (s1 > 1) y = uniform_filter1d(y, s1, 1, mode, cval, 0);
  return y;
}

// ------------------------------------------------------------------ rank
template <class T>
std::vector<T> median_filter(const std::vector<T>& x, size_t size, Mode mode, double cval, int origin) {
  std::vector<T> y(x.size()), buf, win;
  median_core(x.data(), x.size(), size, mode, cval, origin, y.data(), buf, win);
  return y;
}

template <class T>
Mat<T> median_filter(const Mat<T>& x, size_t s0, size_t s1, Mode mode, double cval) {
  if (s0 == 0 || s1 == 0) throw std::invalid_argument("filter size must be >= 1");
  if (s0 == 1 && s1 == 1) return x;
  std::vector<T> buf, win;
  if (s0 == 1)
    return along_axis(x, 1, [&](const T* in, size_t n, T* out) { median_core(in, n, s1, mode, cval, 0, out, buf, win); });
  if (s1 == 1)
    return along_axis(x, 0, [&](const T* in, size_t n, T* out) { median_core(in, n, s0, mode, cval, 0, out, buf, win); });
  Mat<T> y(x.rows, x.cols);
  const size_t rank = (s0 * s1) / 2;
  const int64_t R = static_cast<int64_t>(x.rows), C = static_cast<int64_t>(x.cols);
  const int64_t h0 = static_cast<int64_t>(s0 / 2), h1 = static_cast<int64_t>(s1 / 2);
  win.resize(s0 * s1);
  for (int64_t r = 0; r < R; ++r)
    for (int64_t c = 0; c < C; ++c) {
      size_t k = 0;
      for (int64_t i = 0; i < static_cast<int64_t>(s0); ++i) {
        const int64_t rr = ext_index(r - h0 + i, R, mode);
        for (int64_t j = 0; j < static_cast<int64_t>(s1); ++j) {
          const int64_t cc = ext_index(c - h1 + j, C, mode);
          win[k++] = (rr < 0 || cc < 0) ? static_cast<T>(cval) : x(rr, cc);
        }
      }
      std::nth_element(win.begin(), win.begin() + rank, win.end());
      y(r, c) = win[rank];
    }
  return y;
}

template <class T>
std::vector<T> minimum_filter1d(const std::vector<T>& x, size_t size, Mode mode, double cval, int origin) {
  std::vector<T> y(x.size()), buf;
  minmax_core(x.data(), x.size(), size, mode, cval, origin, false, y.data(), buf);
  return y;
}
template <class T>
std::vector<T> maximum_filter1d(const std::vector<T>& x, size_t size, Mode mode, double cval, int origin) {
  std::vector<T> y(x.size()), buf;
  minmax_core(x.data(), x.size(), size, mode, cval, origin, true, y.data(), buf);
  return y;
}
template <class T>
Mat<T> minimum_filter1d(const Mat<T>& x, size_t size, int axis, Mode mode, double cval, int origin) {
  std::vector<T> buf;
  return along_axis(x, axis, [&](const T* in, size_t n, T* out) { minmax_core(in, n, size, mode, cval, origin, false, out, buf); });
}
template <class T>
Mat<T> maximum_filter1d(const Mat<T>& x, size_t size, int axis, Mode mode, double cval, int origin) {
  std::vector<T> buf;
  return along_axis(x, axis, [&](const T* in, size_t n, T* out) { minmax_core(in, n, size, mode, cval, origin, true, out, buf); });
}
template <class T>
Mat<T> minimum_filter(const Mat<T>& x, size_t s0, size_t s1, Mode mode, double cval) {
  Mat<T> y = x;
  if (s0 > 1) y = minimum_filter1d(y, s0, 0, mode, cval, 0);
  if (s1 > 1) y = minimum_filter1d(y, s1, 1, mode, cval, 0);
  return y;
}
template <class T>
Mat<T> maximum_filter(const Mat<T>& x, size_t s0, size_t s1, Mode mode, double cval) {
  Mat<T> y = x;
  if (s0 > 1) y = maximum_filter1d(y, s0, 0, mode, cval, 0);
  if (s1 > 1) y = maximum_filter1d(y, s1, 1, mode, cval, 0);
  return y;
}

template <class T>
std::vector<T> grey_erosion(const std::vector<T>& x, size_t size, Mode mode, double cval) {
  return minimum_filter1d(x, size, mode, cval, 0);
}
template <class T>
std::vector<T> grey_dilation(const std::vector<T>& x, size_t size, Mode mode, double cval) {
  return maximum_filter1d(x, size, mode, cval, (size & 1) ? 0 : -1);  // scipy flips the origin
}
template <class T>
std::vector<T> grey_opening(const std::vector<T>& x, size_t size, Mode mode, double cval) {
  return grey_dilation(grey_erosion(x, size, mode, cval), size, mode, cval);
}
template <class T>
std::vector<T> grey_closing(const std::vector<T>& x, size_t size, Mode mode, double cval) {
  return grey_erosion(grey_dilation(x, size, mode, cval), size, mode, cval);
}

// ------------------------------------------------------------------ correlate / convolve
namespace {
template <class T>
void correlate_any(const T* x, size_t n, const dvec& w, Mode mode, double cval, int origin, T* out,
                   std::vector<double>& in, std::vector<double>& res, std::vector<double>& buf) {
  in.resize(n);
  res.resize(n);
  if constexpr (np::is_cplx<T>::value) {
    using R = typename T::value_type;
    for (size_t i = 0; i < n; ++i) in[i] = x[i].real();
    correlate_core(in.data(), n, w, mode, cval, origin, res.data(), buf);
    std::vector<double> re(res);
    for (size_t i = 0; i < n; ++i) in[i] = x[i].imag();
    correlate_core(in.data(), n, w, mode, 0.0, origin, res.data(), buf);
    for (size_t i = 0; i < n; ++i) out[i] = T(static_cast<R>(re[i]), static_cast<R>(res[i]));
  } else {
    for (size_t i = 0; i < n; ++i) in[i] = static_cast<double>(x[i]);
    correlate_core(in.data(), n, w, mode, cval, origin, res.data(), buf);
    for (size_t i = 0; i < n; ++i) out[i] = static_cast<T>(res[i]);
  }
}

// convolve1d = correlate1d with reversed weights and origin -origin (-1 more for even)
inline void conv_params(const dvec& w, int origin, dvec& wr, int& o) {
  wr.assign(w.rbegin(), w.rend());
  o = -origin;
  if (!(w.size() & 1)) o -= 1;
}
}  // namespace

template <class T>
std::vector<T> correlate1d(const std::vector<T>& x, const dvec& w, Mode mode, double cval, int origin) {
  std::vector<T> y(x.size());
  std::vector<double> in, res, buf;
  correlate_any(x.data(), x.size(), w, mode, cval, origin, y.data(), in, res, buf);
  return y;
}
template <class T>
std::vector<T> convolve1d(const std::vector<T>& x, const dvec& w, Mode mode, double cval, int origin) {
  dvec wr;
  int o;
  conv_params(w, origin, wr, o);
  return correlate1d(x, wr, mode, cval, o);
}
template <class T>
Mat<T> correlate1d(const Mat<T>& x, const dvec& w, int axis, Mode mode, double cval, int origin) {
  std::vector<double> in, res, buf;
  return along_axis(x, axis, [&](const T* a, size_t n, T* out) { correlate_any(a, n, w, mode, cval, origin, out, in, res, buf); });
}
template <class T>
Mat<T> convolve1d(const Mat<T>& x, const dvec& w, int axis, Mode mode, double cval, int origin) {
  dvec wr;
  int o;
  conv_params(w, origin, wr, o);
  return correlate1d(x, wr, axis, mode, cval, o);
}

// ------------------------------------------------------------------ labelling
bmat cross3() {
  bmat s(3, 3, 0);
  s(0, 1) = s(1, 0) = s(1, 1) = s(1, 2) = s(2, 1) = 1;
  return s;
}
bmat ones(size_t r, size_t c) { return bmat(r, c, 1); }

Labels label(const bmat& mask, const bmat& st) {
  if (st.rows != 3 || st.cols != 3) throw std::invalid_argument("Structure dimensions must be equal to 3");
  std::vector<std::pair<int, int>> nb;
  for (int dr = -1; dr <= 1; ++dr)
    for (int dc = -1; dc <= 1; ++dc)
      if ((dr || dc) && st(dr + 1, dc + 1)) nb.emplace_back(dr, dc);
  Labels L;
  L.lab = imat(mask.rows, mask.cols, 0);
  const int64_t R = static_cast<int64_t>(mask.rows), C = static_cast<int64_t>(mask.cols);
  std::vector<std::pair<int64_t, int64_t>> stack;
  for (int64_t r = 0; r < R; ++r)
    for (int64_t c = 0; c < C; ++c) {
      if (!mask(r, c) || L.lab(r, c)) continue;
      const int id = ++L.n;
      L.lab(r, c) = id;
      stack.assign(1, {r, c});
      while (!stack.empty()) {
        auto [pr, pc] = stack.back();
        stack.pop_back();
        for (auto [dr, dc] : nb) {
          const int64_t rr = pr + dr, cc = pc + dc;
          if (rr < 0 || rr >= R || cc < 0 || cc >= C) continue;
          if (mask(rr, cc) && !L.lab(rr, cc)) {
            L.lab(rr, cc) = id;
            stack.emplace_back(rr, cc);
          }
        }
      }
    }
  return L;
}

Labels1d label(const std::vector<uint8_t>& mask) {
  Labels1d L;
  L.lab.assign(mask.size(), 0);
  for (size_t i = 0; i < mask.size(); ++i) {
    if (!mask[i]) continue;
    if (i == 0 || !mask[i - 1]) ++L.n;
    L.lab[i] = L.n;
  }
  return L;
}

std::vector<Box2> find_objects(const imat& lab, int n) {
  std::vector<Box2> out(static_cast<size_t>(std::max(n, 0)));
  std::vector<uint8_t> seen(out.size(), 0);
  for (size_t r = 0; r < lab.rows; ++r)
    for (size_t c = 0; c < lab.cols; ++c) {
      const int v = lab(r, c);
      if (v <= 0 || v > n) continue;
      Box2& b = out[v - 1];
      if (!seen[v - 1]) {
        b = {r, r + 1, c, c + 1};
        seen[v - 1] = 1;
      } else {
        b.r0 = std::min(b.r0, r); b.r1 = std::max(b.r1, r + 1);
        b.c0 = std::min(b.c0, c); b.c1 = std::max(b.c1, c + 1);
      }
    }
  return out;
}

std::vector<std::pair<size_t, size_t>> find_objects(const std::vector<int32_t>& lab, int n) {
  std::vector<std::pair<size_t, size_t>> out(static_cast<size_t>(std::max(n, 0)), {0, 0});
  std::vector<uint8_t> seen(out.size(), 0);
  for (size_t i = 0; i < lab.size(); ++i) {
    const int v = lab[i];
    if (v <= 0 || v > n) continue;
    auto& b = out[v - 1];
    if (!seen[v - 1]) { b = {i, i + 1}; seen[v - 1] = 1; }
    else { b.first = std::min(b.first, i); b.second = std::max(b.second, i + 1); }
  }
  return out;
}

// ------------------------------------------------------------------ binary morphology
namespace {
// scipy NI_BinaryErosion, iterations=1, border_value 0:
//   erosion  (invert=false): out = AND over structure offsets (outside -> 0)
//   dilation (invert=true, reflected structure / shifted origin): out = OR (outside ignored)
bmat erode_core(const bmat& x, const bmat& st, int64_t o0, int64_t o1, bool invert) {
  std::vector<std::pair<int64_t, int64_t>> off;
  const int64_t c0 = static_cast<int64_t>(st.rows / 2) + o0, c1 = static_cast<int64_t>(st.cols / 2) + o1;
  for (size_t i = 0; i < st.rows; ++i)
    for (size_t j = 0; j < st.cols; ++j)
      if (st(i, j)) off.emplace_back(static_cast<int64_t>(i) - c0, static_cast<int64_t>(j) - c1);
  bmat out(x.rows, x.cols, 0);
  const int64_t R = static_cast<int64_t>(x.rows), C = static_cast<int64_t>(x.cols);
  for (int64_t r = 0; r < R; ++r)
    for (int64_t c = 0; c < C; ++c) {
      uint8_t o;
      if (!invert) {
        o = 1;
        for (auto [dr, dc] : off) {
          const int64_t rr = r + dr, cc = c + dc;
          if (rr < 0 || rr >= R || cc < 0 || cc >= C || !x(rr, cc)) { o = 0; break; }
        }
      } else {
        o = 0;
        for (auto [dr, dc] : off) {
          const int64_t rr = r + dr, cc = c + dc;
          if (rr < 0 || rr >= R || cc < 0 || cc >= C) continue;
          if (x(rr, cc)) { o = 1; break; }
        }
      }
      out(r, c) = o;
    }
  return out;
}

bmat to_row(const std::vector<uint8_t>& v) {
  bmat m(1, v.size());
  for (size_t i = 0; i < v.size(); ++i) m.data[i] = v[i] ? 1 : 0;
  return m;
}
std::vector<uint8_t> from_row(const bmat& m) { return m.data; }
}  // namespace

bmat binary_erosion(const bmat& x, const bmat& st) {
  if (st.empty()) throw std::invalid_argument("structure must not be empty");
  return erode_core(x, st, 0, 0, false);
}

bmat binary_dilation(const bmat& x, const bmat& st) {
  if (st.empty()) throw std::invalid_argument("structure must not be empty");
  bmat rs(st.rows, st.cols);
  for (size_t i = 0; i < st.rows; ++i)
    for (size_t j = 0; j < st.cols; ++j) rs(i, j) = st(st.rows - 1 - i, st.cols - 1 - j);
  const int64_t o0 = (st.rows & 1) ? 0 : -1, o1 = (st.cols & 1) ? 0 : -1;
  return erode_core(x, rs, o0, o1, true);
}

bmat binary_opening(const bmat& x, const bmat& st) { return binary_dilation(binary_erosion(x, st), st); }
bmat binary_closing(const bmat& x, const bmat& st) { return binary_erosion(binary_dilation(x, st), st); }

std::vector<uint8_t> binary_erosion(const std::vector<uint8_t>& x, const std::vector<uint8_t>& s) {
  return from_row(binary_erosion(to_row(x), to_row(s)));
}
std::vector<uint8_t> binary_dilation(const std::vector<uint8_t>& x, const std::vector<uint8_t>& s) {
  return from_row(binary_dilation(to_row(x), to_row(s)));
}
std::vector<uint8_t> binary_opening(const std::vector<uint8_t>& x, const std::vector<uint8_t>& s) {
  return from_row(binary_opening(to_row(x), to_row(s)));
}
std::vector<uint8_t> binary_closing(const std::vector<uint8_t>& x, const std::vector<uint8_t>& s) {
  return from_row(binary_closing(to_row(x), to_row(s)));
}

// ------------------------------------------------------------------ instantiations
#define SI_NDI_REAL(T)                                                                                   \
  template std::vector<T> uniform_filter1d(const std::vector<T>&, size_t, Mode, double, int);           \
  template Mat<T> uniform_filter1d(const Mat<T>&, size_t, int, Mode, double, int);                     \
  template Mat<T> uniform_filter(const Mat<T>&, size_t, size_t, Mode, double);                          \
  template std::vector<T> median_filter(const std::vector<T>&, size_t, Mode, double, int);              \
  template Mat<T> median_filter(const Mat<T>&, size_t, size_t, Mode, double);                           \
  template std::vector<T> minimum_filter1d(const std::vector<T>&, size_t, Mode, double, int);           \
  template std::vector<T> maximum_filter1d(const std::vector<T>&, size_t, Mode, double, int);           \
  template Mat<T> minimum_filter1d(const Mat<T>&, size_t, int, Mode, double, int);                     \
  template Mat<T> maximum_filter1d(const Mat<T>&, size_t, int, Mode, double, int);                     \
  template Mat<T> minimum_filter(const Mat<T>&, size_t, size_t, Mode, double);                          \
  template Mat<T> maximum_filter(const Mat<T>&, size_t, size_t, Mode, double);                          \
  template std::vector<T> grey_erosion(const std::vector<T>&, size_t, Mode, double);                    \
  template std::vector<T> grey_dilation(const std::vector<T>&, size_t, Mode, double);                   \
  template std::vector<T> grey_opening(const std::vector<T>&, size_t, Mode, double);                    \
  template std::vector<T> grey_closing(const std::vector<T>&, size_t, Mode, double);
SI_NDI_REAL(float)
SI_NDI_REAL(double)
#undef SI_NDI_REAL

#define SI_NDI_CORR(T)                                                                         \
  template std::vector<T> correlate1d(const std::vector<T>&, const dvec&, Mode, double, int); \
  template std::vector<T> convolve1d(const std::vector<T>&, const dvec&, Mode, double, int);  \
  template Mat<T> correlate1d(const Mat<T>&, const dvec&, int, Mode, double, int);           \
  template Mat<T> convolve1d(const Mat<T>&, const dvec&, int, Mode, double, int);
SI_NDI_CORR(float)
SI_NDI_CORR(double)
SI_NDI_CORR(cf)
SI_NDI_CORR(cd)
#undef SI_NDI_CORR

}  // namespace si::ndi
