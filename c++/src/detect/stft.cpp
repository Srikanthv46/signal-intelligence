#include "detect/stft.hpp"

#include <cmath>

#include "core/fft.hpp"
#include "core/npcompat.hpp"
#include "core/parallel.hpp"

namespace si::detect {

Spectrogram stft(const cfvec& x, double fs, double fc, size_t nfft, size_t hop) {
  if (hop == 0) hop = nfft;
  Spectrogram S;
  S.nfft = nfft;
  S.hop = hop;
  S.fs = fs;
  S.fc = fc;
  const size_t T = x.size() >= nfft ? (x.size() - nfft) / hop + 1 : 0;
  // np.hanning(nfft).astype(float32): computed in double, then rounded
  fvec w(nfft);
  for (size_t n = 0; n < nfft; ++n)
    w[n] = static_cast<float>(nfft > 1 ? 0.5 - 0.5 * std::cos(2.0 * kPi * static_cast<double>(n) / (nfft - 1.0)) : 1.0);
  fvec w2(nfft);
  for (size_t n = 0; n < nfft; ++n) w2[n] = w[n] * w[n];
  const float wsum = np::sum(w2.data(), w2.size());  // numpy float32 pairwise sum

  cfmat F(T, nfft);
  for (size_t r = 0; r < T; ++r) {
    const cf* src = x.data() + r * hop;
    cf* dst = F.row(r);
    for (size_t n = 0; n < nfft; ++n) dst[n] = src[n] * w[n];  // complex64 * float32
  }
  fft::fft_rows(F, true);
  S.P = fmat(T, nfft);
  const size_t h = nfft / 2;
  parallel_for(T, [&](size_t r) {
    const cf* src = F.row(r);
    float* dst = S.P.row(r);
    for (size_t k = 0; k < nfft; ++k) {
      const float a = np::cabs(src[k]);  // np.abs(complex64)
      dst[(k + h) % nfft] = a * a / wsum;                        // fftshift
    }
  });
  S.f.resize(nfft);  // fc + (arange(nfft) - nfft // 2) * fs / nfft, evaluated left to right
  for (size_t k = 0; k < nfft; ++k)
    S.f[k] = fc + (static_cast<double>(k) - static_cast<double>(nfft / 2)) * fs / static_cast<double>(nfft);
  S.t.resize(T);
  for (size_t r = 0; r < T; ++r) S.t[r] = static_cast<double>(r * hop) / fs;
  return S;
}

}  // namespace si::detect
