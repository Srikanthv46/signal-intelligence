#include "core/dsp.hpp"

#include <algorithm>
#include <stdexcept>

#include "core/fft.hpp"
#include "core/spectral.hpp"

namespace si {

// python Channelizer.__init__: scipy.fft.fft of the complex64 capture (single precision)
Channelizer::Channelizer(const cfvec& x, double fs, double fc) : fs_(fs), fc_(fc), L_(x.size()) {
  X_ = fft::fft(x);
  f_ = fft::fftfreq(L_, 1.0 / fs);
}

// python Channelizer.slice: bins with |f - f_off| <= half_bw (same double expression,
// found by index arithmetic instead of a scan of all L bins), placed at (k - k0) mod M,
// complex64 IFFT of length M, times M/L.
Channelizer::Slice Channelizer::slice(double f_center_abs, double half_bw, double out_fs) const {
  const double f_off = f_center_abs - fc_;
  const int64_t M = std::max<int64_t>(16, py_round_int(static_cast<double>(L_) * out_fs / fs_));
  const std::vector<size_t> sel = select_bins(f_, f_off, half_bw);
  const int64_t k0 = py_round_int(f_off / fs_ * static_cast<double>(L_));
  const int64_t L = static_cast<int64_t>(L_);
  cfvec Y(static_cast<size_t>(M), cf(0.f, 0.f));
  for (size_t i : sel) {
    const int64_t k = static_cast<int64_t>(i) > L / 2 ? static_cast<int64_t>(i) - L : static_cast<int64_t>(i);
    int64_t j = (k - k0) % M;
    if (j < 0) j += M;
    Y[static_cast<size_t>(j)] = X_[i];
  }
  Slice s;
  s.y = fft::ifft(Y);
  const float sc = static_cast<float>(static_cast<double>(M) / static_cast<double>(L_));
  for (auto& v : s.y) v = cf(v.real() * sc, v.imag() * sc);
  s.fs = static_cast<double>(M) * fs_ / static_cast<double>(L_);
  return s;
}

// python core.dsp.burst_mask: m[max(0, int((t0-pad)*fs)) : min(n, int((t1+pad)*fs))] = True
// (Python slice semantics, so a negative stop counts from the end).
std::vector<uint8_t> burst_mask(size_t n, double fs, const std::vector<Box>& bursts, double pad_s) {
  std::vector<uint8_t> m(n, 0);
  const int64_t N = static_cast<int64_t>(n);
  for (const Box& b : bursts) {
    int64_t a = std::max<int64_t>(0, static_cast<int64_t>((b.t0 - pad_s) * fs));
    int64_t e = std::min<int64_t>(N, static_cast<int64_t>((b.t1 + pad_s) * fs));
    if (e < 0) e = std::max<int64_t>(0, e + N);
    a = std::min(a, N);
    for (int64_t i = a; i < e; ++i) m[static_cast<size_t>(i)] = 1;
  }
  return m;
}

}  // namespace si
