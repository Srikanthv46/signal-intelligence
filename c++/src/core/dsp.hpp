// Shared DSP helpers of the pipeline and receivers (python: core/dsp.py).
#pragma once

#include <vector>

#include "core/model.hpp"
#include "core/types.hpp"

namespace si {

// One FFT of the whole capture; any (centre, width, rate) slice is then a small IFFT
// (brick-wall filter: fine for detection and feature extraction). complex64, like the
// Python version (scipy.fft on complex64).
class Channelizer {
 public:
  Channelizer(const cfvec& x, double fs, double fc);

  struct Slice {
    cfvec y;      // complex baseband of [f_center +- half_bw]
    double fs;    // its actual sample rate (M * fs / L)
  };
  // python Channelizer.slice(f_center_abs, half_bw, out_fs)
  Slice slice(double f_center_abs, double half_bw, double out_fs) const;

  const cfvec& X() const { return X_; }   // fft(x), complex64
  const dvec& f() const { return f_; }    // fftfreq(L, 1/fs)
  size_t L() const { return L_; }
  double fs() const { return fs_; }
  double fc() const { return fc_; }

 private:
  double fs_, fc_;
  size_t L_;
  cfvec X_;
  dvec f_;
};

// Boolean mask (n) of the samples inside the bursts (python core.dsp.burst_mask).
std::vector<uint8_t> burst_mask(size_t n, double fs, const std::vector<Box>& bursts, double pad_s = 0.0);

}  // namespace si
