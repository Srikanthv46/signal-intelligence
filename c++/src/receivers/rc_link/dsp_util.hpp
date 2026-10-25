// Small DSP helpers of the rc_link detector (python: receivers/rc_link/detector.py,
// "small DSP helpers"), with numpy's dtypes: where a Python expression mixes a numpy
// float32 with a Python float (NEP 50: the Python float is weak) the arithmetic is float32.
#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include "core/npcompat.hpp"
#include "core/types.hpp"

namespace si::rx::rc_link::detail {

inline int64_t pymod(int64_t a, int64_t b) {
  int64_t m = a % b;
  return m < 0 ? m + b : m;
}

// numpy's temporary elision: in `a * f(b)` the temporary f(b) (when it owns >= 256 KiB and
// has a's shape and dtype) is reused as the output, i.e. numpy evaluates f(b) * a. For
// complex64 the fma-based multiply is not commutative in the last bit, so the operand
// order of np::cmul follows that rule.
inline bool elides(size_t n_elems, size_t itemsize) { return n_elems * itemsize >= 256 * 1024; }
// a * conj(b) as numpy evaluates `a_arr * np.conj(b_arr)` of n elements
inline cf mul_conj(const cf& a, const cf& b, bool elided) {
  return elided ? np::cmul(std::conj(b), a) : np::cmul(a, std::conj(b));
}

// an f_off value as Python holds it: a Python float (double) or a numpy float32 (when it
// came from float32 arithmetic); round(f_off / fs * L) is then done in that precision
struct FOff {
  double v = 0;
  bool f32 = false;
};

struct SliceOut {
  cfvec y;
  double ofs = 0;    // exact output rate M*fs/L
  double f_act = 0;  // exact offset k0*fs/L
};

// python _slice(X, fs, f_off, half_bw, out_fs, t_adv): FFT-domain slice with signed bin indexing
SliceOut slice(const cfvec& X, double fs, FOff f_off, double half_bw, double out_fs, double t_adv = 0.0);

// python _upchirp(N) (complex64)
cfvec upchirp(int64_t N);

// python _peak_interp(S, k) for a float32 |FFT|^2 (the result is float32 arithmetic)
float peak_interp(const float* S, size_t n, size_t k);

// python _wrap(v, period) in dtype T (period/2 and period cast to T as NEP 50 does)
template <class T>
inline T wrap(T v, double period) {
  const T h = static_cast<T>(period / 2), p = static_cast<T>(period);
  return np::mod<T>(v + h, p) - h;
}

// numpy round(x) (half to even) -> int
inline int64_t rint_i(double x) { return static_cast<int64_t>(std::nearbyint(x)); }
inline int64_t rint_i(float x) { return static_cast<int64_t>(std::nearbyintf(x)); }

// python _sym_spec(z, start, N, pad=8): spectrum, fine peak, peak/noise-bin, s2
struct SymSpec {
  fvec S;            // |FFT|^2, 8N bins
  float pk = 0;      // fine peak position (bins of N)
  float r = 0;       // peak / noise-bin
  float s2 = 0;      // noise per sample (float32) ...
  bool s2_floor = false;  // ... or the Python float 1e-30
  // x / (N * s2) as numpy evaluates it
  float per_noise(float x, int64_t N) const {
    return s2_floor ? x / static_cast<float>(static_cast<double>(N) * 1e-30) : x / (static_cast<float>(N) * s2);
  }
};
std::optional<SymSpec> sym_spec(const cfvec& z, int64_t start, int64_t N, int pad = 8);

// python _local_fft(x, fs, t0_us, t1_us, margin_us): fft of the segment and its start (us)
std::pair<cfvec, double> local_fft(const cfvec& x, double fs, double t0_us, double t1_us, double margin_us);

// python _kmeans2(v) on float32 values
std::pair<float, float> kmeans2(const fvec& v, int it = 20);

// float32 pairwise mean and median as numpy
inline float mean_f(const fvec& v) { return np::pairwise_sum(v.data(), v.size()) / static_cast<float>(v.size()); }

}  // namespace si::rx::rc_link::detail
