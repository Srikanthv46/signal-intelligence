// FFTs with numpy.fft / scipy.fft conventions, on pocketfft (third_party/pocketfft is the
// exact commit numpy 2.x and scipy vendor, and fft.cpp is compiled with the compiler's
// default FP contraction like their builds, so results are bit-identical on the same
// platform). Forward unscaled, inverse scaled by 1/n; n zero-pads or truncates.
//
// Precision of float32 input differs between the two Python libraries (measured, numpy 2.2):
//   scipy.fft.*            native single precision                 -> fft, ifft, rfft, irfft
//   numpy.fft.fft / rfft   computed in double, cast to complex64   -> np_fft, np_rfft
//   numpy.fft.ifft / irfft single precision, scale float32(1)/n    -> np_ifft, np_irfft
// complex128 / float64 input is the same in both (use fft, ifft, rfft, irfft).
// Only fft.cpp includes pocketfft.
#pragma once

#include "core/types.hpp"

namespace si::fft {

// In-place / out-of-place 1-D transforms of length n (pointers may alias).
void c2c(const cf* in, cf* out, size_t n, bool forward, float scale = 1.0f);
void c2c(const cd* in, cd* out, size_t n, bool forward, double scale = 1.0);

cfvec fft(const cf* x, size_t len, size_t n = 0);
cfvec ifft(const cf* x, size_t len, size_t n = 0);
cdvec fft(const cd* x, size_t len, size_t n = 0);
cdvec ifft(const cd* x, size_t len, size_t n = 0);
inline cfvec fft(const cfvec& x, size_t n = 0) { return fft(x.data(), x.size(), n); }
inline cfvec ifft(const cfvec& x, size_t n = 0) { return ifft(x.data(), x.size(), n); }
inline cdvec fft(const cdvec& x, size_t n = 0) { return fft(x.data(), x.size(), n); }
inline cdvec ifft(const cdvec& x, size_t n = 0) { return ifft(x.data(), x.size(), n); }

// Real input: n/2+1 bins (scipy.fft.rfft; float64 also numpy.fft.rfft).
cfvec rfft(const float* x, size_t len, size_t n = 0);
cdvec rfft(const double* x, size_t len, size_t n = 0);
inline cfvec rfft(const fvec& x, size_t n = 0) { return rfft(x.data(), x.size(), n); }
inline cdvec rfft(const dvec& x, size_t n = 0) { return rfft(x.data(), x.size(), n); }
// Inverse of rfft: len half-spectrum bins in, n real samples out (n = 0: 2*(len-1)).
fvec irfft(const cf* X, size_t len, size_t n = 0);
dvec irfft(const cd* X, size_t len, size_t n = 0);

// numpy.fft on float32 / complex64 input (see the table above).
cfvec np_fft(const cf* x, size_t len, size_t n = 0);
cfvec np_rfft(const float* x, size_t len, size_t n = 0);
cfvec np_ifft(const cf* x, size_t len, size_t n = 0);
fvec np_irfft(const cf* X, size_t len, size_t n = 0);
inline cfvec np_fft(const cfvec& x, size_t n = 0) { return np_fft(x.data(), x.size(), n); }
inline cfvec np_rfft(const fvec& x, size_t n = 0) { return np_rfft(x.data(), x.size(), n); }
inline cfvec np_ifft(const cfvec& x, size_t n = 0) { return np_ifft(x.data(), x.size(), n); }

// Every row of a (rows x cols) matrix, in place; up to nthreads threads (0 = auto).
void fft_rows(cfmat& m, bool forward, size_t nthreads = 0);
void fft_rows(Mat<cd>& m, bool forward, size_t nthreads = 0);

dvec fftfreq(size_t n, double d = 1.0);   // numpy.fft.fftfreq
dvec rfftfreq(size_t n, double d = 1.0);  // numpy.fft.rfftfreq

template <class T>
std::vector<T> fftshift(const std::vector<T>& x) {
  const size_t n = x.size(), h = n / 2;  // numpy: element 0 moves to index n//2
  std::vector<T> y(n);
  for (size_t i = 0; i < n; ++i) y[(i + h) % n] = x[i];
  return y;
}

template <class T>
std::vector<T> ifftshift(const std::vector<T>& x) {
  const size_t n = x.size(), h = n - n / 2;  // inverse of fftshift for odd n too
  std::vector<T> y(n);
  for (size_t i = 0; i < n; ++i) y[(i + h) % n] = x[i];
  return y;
}

// The best FFT length >= n made of small primes (scipy.fft.next_fast_len(n, real)).
size_t next_fast_len(size_t n);
size_t good_size_real(size_t n);
size_t good_size_cmplx(size_t n);

}  // namespace si::fft
