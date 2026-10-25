#include "core/fft.hpp"

#include <algorithm>
#include <thread>

#include <pocketfft_hdronly.h>

#include "core/parallel.hpp"

namespace si::fft {
namespace {

template <class T>
void c2c_impl(const std::complex<T>* in, std::complex<T>* out, size_t n, bool forward, T scale) {
  if (n == 0) return;
  pocketfft::shape_t shape{n};
  pocketfft::stride_t stride{static_cast<ptrdiff_t>(sizeof(std::complex<T>))};
  pocketfft::c2c(shape, stride, stride, {0}, forward, in, out, scale, 1);
}

template <class T>
std::vector<std::complex<T>> padded(const std::complex<T>* x, size_t len, size_t n) {
  std::vector<std::complex<T>> y(n);
  std::copy(x, x + std::min(len, n), y.begin());
  return y;
}

template <class T>
std::vector<std::complex<T>> fwd(const std::complex<T>* x, size_t len, size_t n, bool forward) {
  if (n == 0) n = len;
  auto y = padded(x, len, n);
  // inverse scale as numpy / scipy pass it: the double 1/n, then cast to T
  c2c_impl(y.data(), y.data(), n, forward, forward ? T(1) : static_cast<T>(1.0 / static_cast<double>(n)));
  return y;
}

template <class T>
std::vector<std::complex<T>> rfft_impl(const T* x, size_t len, size_t n) {
  if (n == 0) n = len;
  std::vector<T> in(n, T(0));
  std::copy(x, x + std::min(len, n), in.begin());
  std::vector<std::complex<T>> out(n / 2 + 1);
  if (n == 0) return out;
  pocketfft::shape_t shape{n};
  pocketfft::stride_t sin{static_cast<ptrdiff_t>(sizeof(T))};
  pocketfft::stride_t sout{static_cast<ptrdiff_t>(sizeof(std::complex<T>))};
  pocketfft::r2c(shape, sin, sout, 0, true, in.data(), out.data(), T(1), 1);
  return out;
}

template <class T>
std::vector<T> irfft_impl(const std::complex<T>* X, size_t len, size_t n, T fct) {
  if (n == 0) n = 2 * (len - 1);
  std::vector<T> out(n);
  if (n == 0) return out;
  std::vector<std::complex<T>> in(n / 2 + 1, std::complex<T>(0, 0));
  std::copy(X, X + std::min(len, in.size()), in.begin());
  pocketfft::shape_t shape{n};
  pocketfft::stride_t sin{static_cast<ptrdiff_t>(sizeof(std::complex<T>))};
  pocketfft::stride_t sout{static_cast<ptrdiff_t>(sizeof(T))};
  pocketfft::c2r(shape, sin, sout, 0, false, in.data(), out.data(), fct, 1);
  return out;
}

cfvec to_c64(const cdvec& y) {
  cfvec out(y.size());
  for (size_t i = 0; i < y.size(); ++i) out[i] = cf(static_cast<float>(y[i].real()), static_cast<float>(y[i].imag()));
  return out;
}

template <class T>
void rows_impl(Mat<std::complex<T>>& m, bool forward, size_t nthreads) {
  if (m.rows == 0 || m.cols == 0) return;
  if (nthreads == 0) nthreads = static_cast<size_t>(threads());
  pocketfft::shape_t shape{m.rows, m.cols};
  const ptrdiff_t es = sizeof(std::complex<T>);
  pocketfft::stride_t stride{static_cast<ptrdiff_t>(m.cols) * es, es};
  const T scale = forward ? T(1) : static_cast<T>(1.0 / static_cast<double>(m.cols));
  pocketfft::c2c(shape, stride, stride, {1}, forward, m.data.data(), m.data.data(), scale,
                 std::min(nthreads, m.rows));
}

}  // namespace

void c2c(const cf* in, cf* out, size_t n, bool forward, float scale) { c2c_impl(in, out, n, forward, scale); }
void c2c(const cd* in, cd* out, size_t n, bool forward, double scale) { c2c_impl(in, out, n, forward, scale); }

cfvec fft(const cf* x, size_t len, size_t n) { return fwd(x, len, n, true); }
cfvec ifft(const cf* x, size_t len, size_t n) { return fwd(x, len, n, false); }
cdvec fft(const cd* x, size_t len, size_t n) { return fwd(x, len, n, true); }
cdvec ifft(const cd* x, size_t len, size_t n) { return fwd(x, len, n, false); }
cfvec rfft(const float* x, size_t len, size_t n) { return rfft_impl(x, len, n); }
cdvec rfft(const double* x, size_t len, size_t n) { return rfft_impl(x, len, n); }
fvec irfft(const cf* X, size_t len, size_t n) {
  const size_t m = n ? n : 2 * (len - 1);
  return irfft_impl(X, len, n, static_cast<float>(1.0 / static_cast<double>(m)));
}
dvec irfft(const cd* X, size_t len, size_t n) {
  const size_t m = n ? n : 2 * (len - 1);
  return irfft_impl(X, len, n, 1.0 / static_cast<double>(m));
}

// numpy 2.x passes fct=1 (a Python int) to its forward ufunc, which resolves to the
// double loop: complex64 in, transformed in double, cast to complex64 on output.
cfvec np_fft(const cf* x, size_t len, size_t n) {
  cdvec y(len);
  for (size_t i = 0; i < len; ++i) y[i] = cd(x[i].real(), x[i].imag());
  return to_c64(fft(y.data(), len, n));
}
cfvec np_rfft(const float* x, size_t len, size_t n) {
  dvec y(x, x + len);
  return to_c64(rfft(y.data(), len, n));
}
// numpy's inverse passes fct = reciprocal(n, dtype=float32): single-precision loop.
cfvec np_ifft(const cf* x, size_t len, size_t n) {
  if (n == 0) n = len;
  auto y = padded(x, len, n);
  c2c_impl(y.data(), y.data(), n, false, 1.0f / static_cast<float>(n));
  return y;
}
fvec np_irfft(const cf* X, size_t len, size_t n) {
  const size_t m = n ? n : 2 * (len - 1);
  return irfft_impl(X, len, n, 1.0f / static_cast<float>(m));
}

void fft_rows(cfmat& m, bool forward, size_t nthreads) { rows_impl(m, forward, nthreads); }
void fft_rows(Mat<cd>& m, bool forward, size_t nthreads) { rows_impl(m, forward, nthreads); }

dvec fftfreq(size_t n, double d) {
  dvec f(n);
  const double val = 1.0 / (static_cast<double>(n) * d);
  const int64_t N = static_cast<int64_t>((n - 1) / 2 + 1);
  for (int64_t i = 0; i < N; ++i) f[i] = static_cast<double>(i) * val;
  for (int64_t i = N; i < static_cast<int64_t>(n); ++i)
    f[i] = static_cast<double>(-static_cast<int64_t>(n / 2) + (i - N)) * val;
  return f;
}

dvec rfftfreq(size_t n, double d) {
  dvec f(n / 2 + 1);
  const double val = 1.0 / (static_cast<double>(n) * d);
  for (size_t i = 0; i < f.size(); ++i) f[i] = static_cast<double>(i) * val;
  return f;
}

size_t next_fast_len(size_t n) {
  if (n <= 6) return n;
  return pocketfft::detail::util::good_size_cmplx(n);
}
size_t good_size_real(size_t n) { return n == 0 ? 0 : pocketfft::detail::util::good_size_real(n); }
size_t good_size_cmplx(size_t n) { return n == 0 ? 0 : pocketfft::detail::util::good_size_cmplx(n); }

}  // namespace si::fft
