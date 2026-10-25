#include "core/spectral.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "core/fft.hpp"
#include "core/ndimage.hpp"
#include "core/npcompat.hpp"

namespace si {

// ------------------------------------------------------------------ bin selection
std::vector<size_t> select_bins(const dvec& f, double f_off, double half_bw) {
  std::vector<size_t> out;
  const int64_t L = static_cast<int64_t>(f.size());
  if (L == 0 || !(half_bw >= 0)) return out;
  const int64_t N = (L - 1) / 2 + 1;   // bins 0..N-1 positive, N..L-1 negative
  const int64_t kmin = -(L / 2), kmax = N - 1;  // signed bin range
  auto idx = [&](int64_t k) { return static_cast<size_t>(k < 0 ? k + L : k); };
  auto in = [&](int64_t k) { return std::fabs(f[idx(k)] - f_off) <= half_bw; };
  // f is monotone in the signed bin k (f = k * df), so the selection is one run of k
  const double df = L > 1 ? f[1] - f[0] : 1.0;
  double kc_d = std::nearbyint(f_off / df);
  if (!std::isfinite(kc_d)) return out;
  kc_d = std::clamp(kc_d, static_cast<double>(kmin), static_cast<double>(kmax));
  const int64_t kc = static_cast<int64_t>(kc_d);
  int64_t start = kmax + 1;
  for (int64_t d : {0, -1, 1, -2, 2}) {
    const int64_t k = kc + d;
    if (k >= kmin && k <= kmax && in(k)) { start = k; break; }
  }
  if (start > kmax) return out;
  int64_t lo = start, hi = start;
  while (lo - 1 >= kmin && in(lo - 1)) --lo;
  while (hi + 1 <= kmax && in(hi + 1)) ++hi;
  // ascending array index: non-negative bins first, then the negative ones
  for (int64_t k = std::max<int64_t>(lo, 0); k <= hi; ++k) out.push_back(idx(k));
  for (int64_t k = lo; k <= std::min<int64_t>(hi, -1); ++k) out.push_back(idx(k));
  return out;
}

// ------------------------------------------------------------------ channelizer
SpectralChannelizer::SpectralChannelizer(const cdvec& x, double fs, const std::vector<std::pair<double, double>>& notch_hz)
    : fs_(fs), L_(x.size()) {
  X_ = fft::fft(x);
  f_ = fft::fftfreq(L_, 1.0 / fs);
  for (const auto& [f0, half] : notch_hz)
    for (size_t i : select_bins(f_, f0, half)) X_[i] = cd(0.0, 0.0);
}

cdvec SpectralChannelizer::slice(double f_off, double half_bw, double out_fs) const {
  const int64_t M = py_round_int(static_cast<double>(L_) * out_fs / fs_);
  if (M <= 0) throw std::invalid_argument("SpectralChannelizer.slice: output length must be positive");
  const std::vector<size_t> sel = select_bins(f_, f_off, half_bw);
  const int64_t k0 = py_round_int(f_off / fs_ * static_cast<double>(L_));
  const int64_t L = static_cast<int64_t>(L_);
  cdvec Y(static_cast<size_t>(M), cd(0.0, 0.0));
  for (size_t i : sel) {
    const int64_t k = static_cast<int64_t>(i) > L / 2 ? static_cast<int64_t>(i) - L : static_cast<int64_t>(i);
    int64_t j = (k - k0) % M;
    if (j < 0) j += M;
    Y[static_cast<size_t>(j)] = X_[i];
  }
  cdvec y = fft::ifft(Y);
  const double s = static_cast<double>(M) / static_cast<double>(L_);
  for (auto& v : y) v = cd(v.real() * s, v.imag() * s);
  return y;
}

// ------------------------------------------------------------------ helpers
namespace {
// |fftshift(fft(x[r*hop : r*hop+nfft] * hanning(nfft)))|**2 per frame (complex128 math)
template <class C>
dmat frames_power(const std::vector<C>& x, size_t nfft, size_t hop, size_t nframes) {
  const dvec w = np::hanning(nfft);
  dmat P(nframes, nfft);
  cdvec seg(nfft);
  const size_t h = nfft / 2;
  for (size_t r = 0; r < nframes; ++r) {
    const C* p = x.data() + r * hop;
    for (size_t i = 0; i < nfft; ++i) seg[i] = cd(static_cast<double>(p[i].real()) * w[i], static_cast<double>(p[i].imag()) * w[i]);
    fft::c2c(seg.data(), seg.data(), nfft, true);
    double* row = P.row(r);
    for (size_t i = 0; i < nfft; ++i) row[(i + h) % nfft] = np::cabs2(seg[i]);
  }
  return P;
}

template <class C>
std::vector<Spur> find_spurs_t(const std::vector<C>& x, double fs, double fc, size_t nfft, double min_prom_db) {
  const size_t n = x.size() / nfft;
  std::vector<Spur> out;
  if (n == 0) return out;
  const dmat S = frames_power(x, nfft, nfft, n);
  dvec psd = np::median_axis0(S);
  for (auto& v : psd) v = 10 * std::log10(v + 1e-20);
  const dvec base = ndi::median_filter(psd, 63);
  dvec f(nfft), prom(nfft);
  for (size_t i = 0; i < nfft; ++i) {
    f[i] = static_cast<double>(static_cast<int64_t>(i) - static_cast<int64_t>(nfft / 2)) / static_cast<double>(nfft) * fs;
    prom[i] = psd[i] - base[i];
  }
  for (size_t k = 0; k < nfft; ++k) {
    if (!(prom[k] > min_prom_db)) continue;
    if (prom[k] < prom[k > 0 ? k - 1 : 0] || prom[k] < prom[std::min(k + 1, nfft - 1)]) continue;
    Spur s;
    s.offset_hz = f[k];
    s.freq_mhz = py_round((fc + f[k]) / 1e6, 4);
    s.prominence_db = py_round(prom[k], 1);
    s.note = std::fabs(f[k]) < 2 * fs / static_cast<double>(nfft) ? "LO leakage (DC)" : "CW tone";
    out.push_back(s);
  }
  return out;
}

template <class C>
NoiseFloor noise_floor_t(const std::vector<C>& x, size_t nfft) {
  const size_t n = x.size() / nfft;
  dmat S = frames_power(x, nfft, nfft, n);
  dvec w2 = np::hanning(nfft);
  for (auto& v : w2) v = v * v;
  const double sw = np::pairwise_sum(w2.data(), nfft);
  for (auto& v : S.data) v /= sw;
  NoiseFloor r;
  r.nfft = nfft;
  dvec fb = np::percentile_axis0(S, 25);
  for (auto& v : fb) v /= 0.2877;
  r.floor_bin = ndi::median_filter(fb, 31);
  r.floor_power_per_sample = np::percentile(r.floor_bin, 20);
  return r;
}

template <class C>
SpectralSpectrogram spectrogram_t(const std::vector<C>& x, double fs, double fc, size_t nfft, size_t hop) {
  if (hop == 0) hop = nfft;
  SpectralSpectrogram r;
  const int64_t n = x.size() >= nfft ? static_cast<int64_t>((x.size() - nfft) / hop + 1) : 0;
  r.P = frames_power(x, nfft, hop, static_cast<size_t>(n));
  r.t_us.resize(static_cast<size_t>(n));
  for (int64_t i = 0; i < n; ++i)
    r.t_us[i] = (static_cast<double>(i * static_cast<int64_t>(hop)) + static_cast<double>(nfft) / 2) / fs * 1e6;
  r.f_mhz.resize(nfft);
  for (size_t i = 0; i < nfft; ++i)
    r.f_mhz[i] = (static_cast<double>(static_cast<int64_t>(i) - static_cast<int64_t>(nfft / 2)) / static_cast<double>(nfft) * fs + fc) / 1e6;
  return r;
}

template <class C>
double measure_snr_t(const std::vector<C>& x, double fs, double fc, double t0_us, double t1_us, double f_lo_mhz,
                     double f_hi_mhz, size_t nfft) {
  SpectralSpectrogram S = spectrogram_t(x, fs, fc, nfft, 0);
  dvec w2 = np::hanning(nfft);
  for (auto& v : w2) v = v * v;
  const double sw = np::pairwise_sum(w2.data(), nfft);
  for (auto& v : S.P.data) v /= sw;
  const dvec fl = noise_floor_t(x, nfft).floor_bin;
  std::vector<size_t> ti, fi;
  for (size_t i = 0; i < S.t_us.size(); ++i)
    if (S.t_us[i] >= t0_us && S.t_us[i] <= t1_us) ti.push_back(i);
  for (size_t i = 0; i < S.f_mhz.size(); ++i)
    if (S.f_mhz[i] >= f_lo_mhz && S.f_mhz[i] <= f_hi_mhz) fi.push_back(i);
  dmat sub(ti.size(), fi.size());
  for (size_t a = 0; a < ti.size(); ++a)
    for (size_t b = 0; b < fi.size(); ++b) sub(a, b) = S.P(ti[a], fi[b]);
  dvec sig = np::mean_axis0(sub);
  dvec flf(fi.size());
  for (size_t b = 0; b < fi.size(); ++b) { sig[b] -= fl[fi[b]]; flf[b] = fl[fi[b]]; }
  const double sm = np::mean(sig);
  return 10 * std::log10(std::max(sm, 1e-12) / np::mean(flf));
}
}  // namespace

Json to_json(const Spur& s) {
  Json j = Json::object();
  j["offset_hz"] = s.offset_hz;
  j["freq_mhz"] = s.freq_mhz;
  j["prominence_db"] = s.prominence_db;
  j["note"] = s.note;
  return j;
}

std::vector<Spur> find_spurs(const cfvec& x, double fs, double fc, size_t nfft, double p) { return find_spurs_t(x, fs, fc, nfft, p); }
std::vector<Spur> find_spurs(const cdvec& x, double fs, double fc, size_t nfft, double p) { return find_spurs_t(x, fs, fc, nfft, p); }

CpMetric cp_metric(const cd* y, size_t len, size_t N, size_t cp) {
  CpMetric r;
  if (len <= N) return r;
  const size_t n = len - N;
  r.c.resize(n);
  dvec e(n);
  for (size_t i = 0; i < n; ++i) {
    r.c[i] = np::cmul(y[i], std::conj(y[i + N]));
    e[i] = (np::cabs2(y[i]) + np::cabs2(y[i + N])) / 2;
  }
  if (n < cp || cp == 0) return r;
  const size_t nv = n - cp + 1;
  r.m.resize(nv);
  for (size_t j = 0; j < nv; ++j) {  // np.convolve(.., ones(cp), "valid"): window sums in order
    cd sc(0.0, 0.0);
    double se = 0.0;
    for (size_t k = 0; k < cp; ++k) { sc += r.c[j + k]; se += e[j + k]; }
    r.m[j] = np::cabs(sc) / (se + 1e-20);
  }
  return r;
}
CpMetric cp_metric(const cdvec& y, size_t N, size_t cp) { return cp_metric(y.data(), y.size(), N, cp); }

dvec folded_scores(const dvec& m, size_t sym, size_t W, size_t hop) {
  dvec out;
  if (m.size() <= W || hop == 0 || sym == 0) return out;
  const size_t k = W / sym;
  dvec prof(sym);
  for (size_t a = 0; a < m.size() - W; a += hop) {
    std::fill(prof.begin(), prof.end(), 0.0);
    for (size_t r = 0; r < k; ++r)
      for (size_t c = 0; c < sym; ++c) prof[c] += m[a + r * sym + c];
    for (auto& v : prof) v = v / static_cast<double>(k);
    out.push_back(*std::max_element(prof.begin(), prof.end()) - np::median(prof));
  }
  return out;
}

NoiseFloor noise_floor_psd(const cfvec& x, size_t nfft) { return noise_floor_t(x, nfft); }
NoiseFloor noise_floor_psd(const cdvec& x, size_t nfft) { return noise_floor_t(x, nfft); }

SpectralSpectrogram spectrogram(const cfvec& x, double fs, double fc, size_t nfft, size_t hop) { return spectrogram_t(x, fs, fc, nfft, hop); }
SpectralSpectrogram spectrogram(const cdvec& x, double fs, double fc, size_t nfft, size_t hop) { return spectrogram_t(x, fs, fc, nfft, hop); }

std::pair<double, double> robust_threshold(const dvec& v, double k_sigma) {
  const double med = np::median(v);
  dvec d(v.size());
  for (size_t i = 0; i < v.size(); ++i) d[i] = std::fabs(v[i] - med);
  return {med + k_sigma * 1.4826 * np::median(d), med};
}

double measure_snr(const cfvec& x, double fs, double fc, double t0, double t1, double flo, double fhi, size_t nfft) {
  return measure_snr_t(x, fs, fc, t0, t1, flo, fhi, nfft);
}
double measure_snr(const cdvec& x, double fs, double fc, double t0, double t1, double flo, double fhi, size_t nfft) {
  return measure_snr_t(x, fs, fc, t0, t1, flo, fhi, nfft);
}

}  // namespace si
