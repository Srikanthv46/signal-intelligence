#include "receivers/rc_link/dsp_util.hpp"

#include <algorithm>
#include <cmath>

#include "core/fft.hpp"

namespace si::rx::rc_link::detail {

SliceOut slice(const cfvec& X, double fs, FOff f_off, double half_bw, double out_fs, double t_adv) {
  const int64_t L = static_cast<int64_t>(X.size());
  const int64_t M = py_round_int(static_cast<double>(L) * out_fs / fs);
  int64_t k0;
  if (f_off.f32) {
    const float q = (static_cast<float>(f_off.v) / static_cast<float>(fs)) * static_cast<float>(L);
    k0 = rint_i(q);
  } else {
    k0 = py_round_int(f_off.v / fs * static_cast<double>(L));
  }
  const int64_t hb = static_cast<int64_t>(half_bw / fs * static_cast<double>(L));
  cfvec Y(static_cast<size_t>(M), cf(0.f, 0.f));
  const double two_pi = 2 * kPi, inv_m = 1.0 / static_cast<double>(M);
  for (int64_t k = -hb; k <= hb; ++k) {
    cf v = X[static_cast<size_t>(pymod(k0 + k, L))];
    if (t_adv != 0.0) {
      // np.exp(2j*np.pi*k/M*t_adv) (complex128), v promoted to complex128
      const double th = ((two_pi * static_cast<double>(k)) * inv_m) * t_adv;
      const cd e(std::cos(th), std::sin(th));
      const cd p = np::cmul(cd(v.real(), v.imag()), e);
      v = cf(static_cast<float>(p.real()), static_cast<float>(p.imag()));
    }
    Y[static_cast<size_t>(pymod(k, M))] = v;
  }
  SliceOut o;
  o.y = fft::ifft(Y);
  const float sc = static_cast<float>(static_cast<double>(M) / static_cast<double>(L));
  for (auto& v : o.y) v = cf(v.real() * sc, v.imag() * sc);
  o.ofs = static_cast<double>(M) * fs / static_cast<double>(L);
  o.f_act = static_cast<double>(k0) * fs / static_cast<double>(L);
  return o;
}

cfvec upchirp(int64_t N) {
  cfvec u(static_cast<size_t>(N));
  for (int64_t n = 0; n < N; ++n) {
    const double v = kPi * (static_cast<double>(n * n) / static_cast<double>(N) - static_cast<double>(n));
    u[static_cast<size_t>(n)] = cf(static_cast<float>(std::cos(v)), static_cast<float>(std::sin(v)));
  }
  return u;
}

float peak_interp(const float* S, size_t n, size_t k) {
  const float e = static_cast<float>(1e-30);
  const float a = std::log(S[(k + n - 1) % n] + e), b = std::log(S[k] + e), c = std::log(S[(k + 1) % n] + e);
  const float d = (a - 2.0f * b) + c;
  if (d < 0) return static_cast<float>(k) + (0.5f * (a - c)) / d;
  return static_cast<float>(k);
}

std::optional<SymSpec> sym_spec(const cfvec& z, int64_t start, int64_t N, int pad) {
  if (start < 0 || start + N > static_cast<int64_t>(z.size())) return std::nullopt;
  const cf* seg = z.data() + start;
  const size_t nf = static_cast<size_t>(N) * pad;
  cfvec F = fft::fft(seg, static_cast<size_t>(N), nf);
  SymSpec p;
  p.S.resize(nf);
  for (size_t i = 0; i < nf; ++i) p.S[i] = np::cabs2(F[i]);
  fvec a2(static_cast<size_t>(N));
  for (int64_t i = 0; i < N; ++i) a2[i] = np::cabs2(seg[i]);
  const float E = np::pairwise_sum(a2.data(), a2.size());
  const size_t k = np::argmax(p.S);
  const float pk = p.S[k];
  const float s2 = (E - pk / static_cast<float>(N)) / static_cast<float>(N - 1);
  if (1e-30 > static_cast<double>(s2)) p.s2_floor = true;  // python max(s2, 1e-30)
  else p.s2 = s2;
  p.pk = peak_interp(p.S.data(), nf, k) / static_cast<float>(pad);
  p.r = p.per_noise(pk, N);
  return p;
}

std::pair<cfvec, double> local_fft(const cfvec& x, double fs, double t0_us, double t1_us, double margin_us) {
  const int64_t n = static_cast<int64_t>(x.size());
  const int64_t i0 = std::max<int64_t>(0, static_cast<int64_t>((t0_us - margin_us) * 1e-6 * fs));
  const int64_t i1 = std::min<int64_t>(n, static_cast<int64_t>((t1_us + margin_us) * 1e-6 * fs) + 1);
  cfvec X = i1 > i0 ? fft::fft(x.data() + i0, static_cast<size_t>(i1 - i0)) : cfvec();
  return {std::move(X), static_cast<double>(i0) / fs * 1e6};
}

std::pair<float, float> kmeans2(const fvec& v, int it) {
  float a = np::percentile(v, 25), b = np::percentile(v, 75);
  fvec in, out;
  in.reserve(v.size());
  out.reserve(v.size());
  for (int i = 0; i < it; ++i) {
    in.clear();
    out.clear();
    for (float x : v) (std::fabs(x - a) < std::fabs(x - b) ? in : out).push_back(x);
    if (out.empty() || in.empty()) break;
    a = mean_f(in);
    b = mean_f(out);
  }
  return {std::min(a, b), std::max(a, b)};
}

}  // namespace si::rx::rc_link::detail
