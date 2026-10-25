// python: receivers/analog_vtx/detector.py, ported line by line (same dtypes: the STFT power
// is float32 as numpy keeps it, folds and z-scores float64). Where a numpy expression mixes a
// float32 array with a numpy float64 scalar (strong under NEP 50) the C++ promotes the same way.
#include "receivers/analog_vtx/detector.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <tuple>

#include "core/fft.hpp"
#include "core/ndimage.hpp"
#include "core/npcompat.hpp"
#include "core/parallel.hpp"
#include "core/pyfmt.hpp"
#include "core/signal.hpp"

namespace si::rx::analog_vtx {

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

// Python's x % m for ints (result has the sign of m)
inline int64_t imod(int64_t a, int64_t m) {
  int64_t r = a % m;
  return r < 0 ? r + m : r;
}

// np.mod(x, 1.0) for finite x: x - floor(x) is exact (the fractional part of a double is
// representable) and rounds like fmod(x, 1) + 1 for negative x, so it is bit-identical and
// several times faster than fmod in the fold loops.
inline double mod1(double x) { return x - std::floor(x); }

// numpy float32 reductions over axis 0 (rows added one after the other) of the selected rows
fvec mean_rows(const fmat& M, const std::vector<uint8_t>* sel) {
  fvec acc(M.cols, 0.f);
  size_t n = 0;
  for (size_t r = 0; r < M.rows; ++r) {
    if (sel && !(*sel)[r]) continue;
    const float* p = M.row(r);
    for (size_t c = 0; c < M.cols; ++c) acc[c] += p[c];
    ++n;
  }
  for (auto& v : acc) v = static_cast<float>(static_cast<double>(v) / static_cast<double>(n));
  return acc;
}
fvec var_rows(const fmat& M, const std::vector<uint8_t>* sel) {
  const fvec m = mean_rows(M, sel);
  fvec acc(M.cols, 0.f);
  size_t n = 0;
  for (size_t r = 0; r < M.rows; ++r) {
    if (sel && !(*sel)[r]) continue;
    const float* p = M.row(r);
    for (size_t c = 0; c < M.cols; ++c) {
      const float d = p[c] - m[c];
      acc[c] += d * d;
    }
    ++n;
  }
  for (auto& v : acc) v = static_cast<float>(static_cast<double>(v) / static_cast<double>(n));
  return acc;
}
// np.median(M[sel, c]) for a float32 column
float median_col(const fmat& M, size_t c, const std::vector<uint8_t>* sel) {
  fvec v;
  v.reserve(M.rows);
  for (size_t r = 0; r < M.rows; ++r)
    if (!sel || (*sel)[r]) v.push_back(M(r, c));
  return np::median_inplace(v.data(), v.size());
}

// np.searchsorted(a_float32, python_float): compared in float64
template <class T>
size_t searchsorted_d(const std::vector<T>& a, double v) {
  size_t lo = 0, hi = a.size();
  while (lo < hi) {
    size_t mid = (lo + hi) / 2;
    if (static_cast<double>(a[mid]) < v) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

template <class T>
std::vector<T> cumsum_seq(const std::vector<T>& x) {
  std::vector<T> y(x.size());
  T s = T(0);
  for (size_t i = 0; i < x.size(); ++i) y[i] = (s += x[i]);
  return y;
}

// np.bincount(b, w, nb) / np.maximum(np.bincount(b, minlength=nb), 1)
dvec fold_prof(const std::vector<int64_t>& b, const double* e, size_t n, int nbf) {
  dvec S(nbf, 0.0);
  std::vector<int64_t> N(nbf, 0);
  for (size_t i = 0; i < n; ++i) {
    S[b[i]] += e[i];
    ++N[b[i]];
  }
  for (int i = 0; i < nbf; ++i) S[i] = S[i] / static_cast<double>(std::max<int64_t>(N[i], 1));
  return S;
}

std::vector<int64_t> phase_bins(const dvec& tt, double rate, int nb) {
  std::vector<int64_t> b(tt.size());
  for (size_t i = 0; i < tt.size(); ++i) b[i] = static_cast<int64_t>(mod1(tt[i] * rate) * nb);
  return b;
}

dvec time_axis(size_t n, double fs) {
  dvec tt(n);
  for (size_t i = 0; i < n; ++i) tt[i] = static_cast<double>(i) / fs;
  return tt;
}

double parabolic(double y0, double y1, double y2) { return 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2 + 1e-30); }

// python _xcorr_shift(a, b)[0]
double xcorr_shift(const dvec& a, const dvec& b) {
  const size_t n = a.size();
  const double ma = np::mean(a), mb = np::mean(b);
  cdvec A(n), B(n);
  for (size_t i = 0; i < n; ++i) {
    A[i] = cd(a[i] - ma, 0.0);
    B[i] = cd(b[i] - mb, 0.0);
  }
  A = fft::fft(A);
  B = fft::fft(B);
  for (size_t i = 0; i < n; ++i) A[i] = np::cmul(A[i], std::conj(B[i]));
  const cdvec ci = fft::ifft(A);
  dvec c(n);
  for (size_t i = 0; i < n; ++i) c[i] = ci[i].real();
  const size_t k = np::argmax(c);
  const double y0 = c[(k + n - 1) % n], y1 = c[k], y2 = c[(k + 1) % n];
  const double s = static_cast<double>(k) + parabolic(y0, y1, y2);
  const double nd = static_cast<double>(n);
  return np::mod(s + nd / 2, nd) - nd / 2;
}

// python _lead_edge(p, smooth)
dvec lead_edge(dvec p, int smooth) {
  if (smooth > 1) p = ndi::uniform_filter1d(p, smooth, ndi::Mode::Wrap);
  const size_t n = p.size();
  dvec d(n);
  for (size_t i = 0; i < n; ++i) d[i] = p[i] - p[(i + n - 1) % n];
  d = ndi::uniform_filter1d(d, 3, ndi::Mode::Wrap);
  for (auto& v : d) v = std::max(v, 0.0);
  return d;
}

// python _sync_start(prof, rate, nbf)
double sync_start(const dvec& prof, double rate, int nbf) {
  const int64_t w = std::max<int64_t>(py_round_int(kSyncS * rate * nbf), 1);
  const dvec sm = ndi::convolve1d(prof, dvec(static_cast<size_t>(w), 1.0), ndi::Mode::Wrap);
  const size_t c = np::argmax(sm);
  return np::mod(static_cast<double>(c) - static_cast<double>(w) / 2, static_cast<double>(nbf)) / nbf / rate;
}

double npr(double x, int n) { return np::round(x, n); }   // round(numpy float64, n)
inline int64_t pymod(int64_t a, int64_t b) { const int64_t r = a % b; return (r != 0 && ((r < 0) != (b < 0))) ? r + b : r; }
double pyr(double x, int n) { return py_round(x, n); }    // round(python float, n)

}  // namespace

// ============================================================== constants / channels
const std::vector<Nominal>& nominal() {
  static const std::vector<Nominal> v = {{"pal", 15625.0, 625, 4.43361875e6}, {"ntsc", 15734.264, 525, 3.579545e6}};
  return v;
}
const Nominal& nominal(const std::string& s) { return s == "pal" ? nominal()[0] : nominal()[1]; }

namespace {
const std::map<std::string, std::vector<int>>& channels() {
  static const std::map<std::string, std::vector<int>> c = {
      {"A", {5865, 5845, 5825, 5805, 5785, 5765, 5745, 5725}},
      {"B", {5733, 5752, 5771, 5790, 5809, 5828, 5847, 5866}},
      {"E", {5705, 5685, 5665, 5645, 5885, 5905, 5925, 5945}},
      {"F", {5740, 5760, 5780, 5800, 5820, 5840, 5860, 5880}},
      {"R", {5658, 5695, 5732, 5769, 5806, 5843, 5880, 5917}},
      {"L", {5362, 5399, 5436, 5473, 5510, 5547, 5584, 5621}},
      {"LB", {5333, 5373, 5413, 5453, 5493, 5533, 5573, 5613}},  // Lowband 40 MHz step (research)
      {"U", {5325, 5348, 5366, 5384, 5402, 5420, 5438, 5456}},
      {"O", {5474, 5492, 5510, 5528, 5546, 5564, 5582, 5600}},
      {"H", {5653, 5693, 5733, 5773, 5813, 5853, 5893, 5933}},
      {"1G3", {1080, 1120, 1160, 1200, 1240, 1280, 1320, 1360, 1258}},
      {"2G4", {2414, 2432, 2450, 2468}},
      {"2G4X", {2410, 2430, 2450, 2470, 2490, 2510, 2390, 2370}},
      {"3G3A", {3330, 3350, 3370, 3390, 3410, 3430, 3450, 3470}},  // 3.3 GHz, two 8-channel groups (PLAN G9)
      {"3G3B", {3170, 3190, 3210, 3230, 3250, 3270, 3290, 3310}},
  };
  return c;
}
const std::vector<std::string> kBandPriority = {"R", "F", "A", "B", "E", "L", "1G3", "2G4", "2G4X", "H", "U", "O", "LB",
                                                 "3G3A", "3G3B"};
}  // namespace

std::vector<std::pair<std::string, int>> match_channel(double f_mhz, double tol_mhz) {
  std::vector<std::tuple<double, int, std::string, int>> out;
  for (size_t p = 0; p < kBandPriority.size(); ++p) {
    const std::string& band = kBandPriority[p];
    const auto& tab = channels().at(band);
    for (size_t i = 0; i < tab.size(); ++i) {
      const double d = f_mhz - tab[i];
      if (std::fabs(d) <= tol_mhz) {
        std::string name = (band != "1G3" || i < 8) ? band + std::to_string(i + 1) : "1G3-1258";
        out.emplace_back(std::fabs(d), static_cast<int>(p), name, tab[i]);
      }
    }
  }
  std::sort(out.begin(), out.end());
  std::vector<std::pair<std::string, int>> r;
  for (const auto& t : out) r.emplace_back(std::get<2>(t), std::get<3>(t));
  return r;
}

// ============================================================== time-frequency fold
TfPower tf_power(const cf* x, size_t len, double fs) {
  size_t nfft = static_cast<size_t>(std::pow(2.0, std::nearbyint(std::log2(fs * 1e-6))));
  nfft = std::max<size_t>(nfft, 8);
  const size_t hop = nfft / 2;
  const size_t n = (len - nfft) / hop + 1;
  const dvec wd = np::hanning(nfft);
  fvec win(nfft);
  for (size_t j = 0; j < nfft; ++j) win[j] = static_cast<float>(wd[j]);
  TfPower r;
  r.P = fmat(n, nfft);
  parallel_for((n + 4095) / 4096, [&](size_t blk) {
    cfvec seg(nfft);
    const size_t h = nfft / 2;
    for (size_t i = blk * 4096; i < std::min(n, (blk + 1) * 4096); ++i) {
      const cf* xi = x + i * hop;
      for (size_t j = 0; j < nfft; ++j) seg[j] = cf(xi[j].real() * win[j], xi[j].imag() * win[j]);
      const cfvec X = fft::np_fft(seg);
      float* row = r.P.row(i);
      for (size_t j = 0; j < nfft; ++j) {
        const float a = np::cabs(X[j]);
        row[(j + h) % nfft] = a * a;   // fftshift
      }
    }
  });
  r.t.resize(n);
  for (size_t i = 0; i < n; ++i)
    r.t[i] = (static_cast<double>(static_cast<int64_t>(i * hop)) + static_cast<double>(nfft) / 2) / fs;
  r.f.resize(nfft);
  for (size_t j = 0; j < nfft; ++j)
    r.f[j] = static_cast<double>(static_cast<int64_t>(j) - static_cast<int64_t>(nfft / 2)) / static_cast<double>(nfft) * fs;
  fvec med = np::median_axis0(r.P);
  for (auto& m : med) m = m + 1e-30f;
  r.q = fmat(n, nfft);
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < nfft; ++j) r.q(i, j) = std::min(r.P(i, j) / med[j], kClip);
  return r;
}

Fold fold(const fmat& q, const dvec& t, double rate, int nb) {
  const size_t nbins = q.cols;
  fmat acc(nb, nbins, 0.f);
  Fold r;
  r.cnt.assign(nb, 0.0);
  std::vector<int64_t> cnt(nb, 0);
  for (size_t i = 0; i < t.size(); ++i) {
    const int64_t b = static_cast<int64_t>(mod1(t[i] * rate) * nb);
    float* a = acc.row(b);
    const float* qi = q.row(i);
    for (size_t k = 0; k < nbins; ++k) a[k] += qi[k];
    ++cnt[b];
  }
  r.F = dmat(nb, nbins);
  bool empty = false;
  for (int b = 0; b < nb; ++b) {
    r.cnt[b] = static_cast<double>(cnt[b]);
    empty |= cnt[b] == 0;
    const double d = std::max(r.cnt[b], 1.0);
    for (size_t k = 0; k < nbins; ++k) r.F(b, k) = static_cast<double>(acc(b, k)) / d;
  }
  if (empty) {
    for (size_t k = 0; k < nbins; ++k) {
      dvec col;
      for (int b = 0; b < nb; ++b)
        if (cnt[b] > 0) col.push_back(r.F(b, k));
      const double m = np::median(col);
      for (int b = 0; b < nb; ++b)
        if (cnt[b] == 0) r.F(b, k) = m;
    }
  }
  return r;
}

std::pair<dmat, dmat> zscore(const dmat& F, const fvec& q_sd, double n_per_bin, int nb_box) {
  const size_t nb = F.rows, nk = F.cols;
  const dvec med = np::median_axis0(F);
  dmat dev(nb, nk);
  for (size_t b = 0; b < nb; ++b)
    for (size_t k = 0; k < nk; ++k) dev(b, k) = std::fabs(F(b, k) - med[k]);
  const dvec mm = np::median_axis0(dev);
  const double sq = std::sqrt(n_per_bin);
  dvec den(nk);
  for (size_t k = 0; k < nk; ++k) {
    const double mad = 1.4826 * mm[k];
    const double fl = static_cast<double>(0.7f * q_sd[k]) / sq + 1e-12;
    den[k] = std::isnan(mad) || std::isnan(fl) ? std::numeric_limits<double>::quiet_NaN() : std::max(mad, fl);
  }
  dmat z(nb, nk);
  for (size_t b = 0; b < nb; ++b)
    for (size_t k = 0; k < nk; ++k) z(b, k) = (F(b, k) - med[k]) / den[k];
  const dvec w(static_cast<size_t>(nb_box), 1.0 / std::sqrt(static_cast<double>(nb_box)));
  dmat zb = ndi::convolve1d(z, w, 0, ndi::Mode::Wrap);
  dmat zs = zb;
  const double r2 = std::sqrt(2.0);
  for (size_t b = 0; b < nb; ++b)
    for (size_t k = 0; k + 1 < nk; ++k) zs(b, k) = std::max(zb(b, k), (zb(b, k) + zb(b, k + 1)) / r2);
  return {zb, zs};
}

Scanner::Scanner(const TfPower& tf, double ppm_max) : tf_(tf) {
  q_sd_ = var_rows(tf.q, nullptr);
  for (auto& v : q_sd_) v = std::sqrt(v);
  const double D = tf.t.back() - tf.t.front();
  const double step = 0.6e-6 / D * 1e6;
  const dvec ppms = np::arange<double>(-ppm_max, ppm_max + step / 2, step);
  for (size_t s = 0; s < nominal().size(); ++s)
    for (double p : ppms) hyp_.push_back({static_cast<int>(s), p, nominal()[s].fl * (1 + p * 1e-6), dmat()});
  parallel_for(hyp_.size(), [&](size_t i) {
    Hyp& h = hyp_[i];
    const int nb_box = static_cast<int>(py_round_int(kSyncS * nominal()[h.s].fl * kNB));
    const Fold F = fold(tf_.q, tf_.t, h.rate);
    h.zs = zscore(F.F, q_sd_, np::mean(F.cnt), nb_box).second;
  });
}

ScanBest Scanner::best(const std::vector<uint8_t>* col_mask) const {
  ScanBest best;
  bool have = false;
  size_t bi = 0;
  for (size_t i = 0; i < hyp_.size(); ++i) {
    const dmat& zs = hyp_[i].zs;
    double s = -kInf;
    size_t arg = 0;
    bool first = true;
    for (size_t b = 0; b < zs.rows; ++b)
      for (size_t k = 0; k < zs.cols; ++k) {
        const double v = (col_mask && (*col_mask)[k]) ? -kInf : zs(b, k);
        if (first || v > s || (std::isnan(v) && !std::isnan(s))) {
          s = v;
          arg = b * zs.cols + k;
          first = false;
          if (std::isnan(v)) goto done;
        }
      }
  done:
    if (!have || s > best.score) {
      have = true;
      bi = i;
      best.score = s;
      best.b0 = static_cast<int>(arg / zs.cols);
      best.k0 = static_cast<int>(arg % zs.cols);
    }
  }
  const Hyp& h = hyp_[bi];
  best.std_name = nominal()[h.s].name;
  best.ppm = h.ppm;
  best.rate = h.rate;
  const int nb_box = static_cast<int>(py_round_int(kSyncS * nominal()[h.s].fl * kNB));
  const Fold F = fold(tf_.q, tf_.t, h.rate);
  best.zb = zscore(F.F, q_sd_, np::mean(F.cnt), nb_box).first;
  return best;
}

dvec Scanner::profile() const {
  dvec out;
  for (const Hyp& h : hyp_) out.push_back(*std::max_element(h.zs.data.begin(), h.zs.data.end()));
  return out;
}

std::pair<int, int> emitter_columns(const dmat& zb, int b0, int k0, int nbins_win, double thr) {
  const int nb = static_cast<int>(zb.rows), nk = static_cast<int>(zb.cols);
  dvec s(nk, -kInf);
  for (int i = -nbins_win; i < 2 * nbins_win; ++i) {
    const int ph = static_cast<int>(imod(i + b0, nb));
    for (int k = 0; k < nk; ++k) s[k] = std::max(s[k], std::fabs(zb(ph, k)));
  }
  int lo = k0, hi = k0, gap = 0;
  while (lo > 0 && gap < 3) {
    --lo;
    gap = s[lo] > thr ? 0 : gap + 1;
  }
  lo += gap;
  gap = 0;
  while (hi < nk - 1 && gap < 3) {
    ++hi;
    gap = s[hi] > thr ? 0 : gap + 1;
  }
  hi -= gap;
  return {lo, hi};
}

// ============================================================== refinement
Band slice_band(const SpectralChannelizer& ch, double f_off, double half_bw, double out_fs) {
  const size_t L = ch.L();
  const double fs = ch.fs();
  const int64_t M = py_round_int(static_cast<double>(L) * out_fs / fs);
  const std::vector<size_t> sel = select_bins(ch.f(), f_off, half_bw);
  const int64_t k0 = py_round_int(f_off / fs * static_cast<double>(L));
  cdvec Y(static_cast<size_t>(M), cd(0, 0));
  for (size_t i : sel) {
    const int64_t kk = static_cast<int64_t>(std::nearbyint(ch.f()[i] / fs * static_cast<double>(L)));
    Y[static_cast<size_t>(imod(kk - k0, M))] = ch.X()[i];
  }
  Band b;
  b.y = fft::ifft(Y);
  const double sc = static_cast<double>(M) / static_cast<double>(L);
  for (auto& v : b.y) v = cd(v.real() * sc, v.imag() * sc);
  b.fs = static_cast<double>(M) * fs / static_cast<double>(L);
  return b;
}

double grid_line_rate(const dvec& e, double efs, double rate, double span_ppm) {
  const dvec tt = time_axis(e.size(), efs);
  auto sharp = [&](double r, int nb) {
    const std::vector<int64_t> b = phase_bins(tt, r, nb);
    dvec p = fold_prof(b, e.data(), e.size(), nb);
    p = ndi::uniform_filter1d(p, std::max(nb / 256, 1), ndi::Mode::Wrap);
    const double m = np::mean(p);
    dvec d(p.size());
    for (size_t i = 0; i < p.size(); ++i) {
      const double v = p[i] - m;
      d[i] = v * v;
    }
    return np::sum(d);
  };
  const std::tuple<double, double, int> passes[2] = {{span_ppm, 4.0, 256}, {6.0, 0.25, 1024}};
  for (const auto& [span, step, nb] : passes) {
    const dvec ppms = np::arange<double>(-span, span + step / 2, step);
    dvec met = pmap(ppms, [&, nb = nb](double p) { return sharp(rate * (1 + p * 1e-6), nb); });
    if (step < 1) met = ndi::uniform_filter1d(met, 3);
    const size_t k = np::argmax(met);
    double d = 0.0;
    if (0 < k && k < met.size() - 1) d = parabolic(met[k - 1], met[k], met[k + 1]);
    rate = rate * (1 + (ppms[k] + d * step) * 1e-6);
  }
  return rate;
}

std::pair<double, double> drift_line_rate(const dvec& e, double efs, double rate, int nblk, int nbf, int iters,
                                          int smooth) {
  const size_t n = e.size();
  const dvec tt = time_axis(n, efs);
  const dvec ed = np::linspace<double>(0.0, static_cast<double>(n), static_cast<size_t>(nblk + 1));
  std::vector<size_t> edges(ed.size());
  for (size_t i = 0; i < ed.size(); ++i) edges[i] = static_cast<size_t>(ed[i]);
  double sig = kInf;
  for (int it = 0; it < iters; ++it) {
    const std::vector<int64_t> b = phase_bins(tt, rate, nbf);
    dvec S(nbf, 0.0);
    std::vector<int64_t> N(nbf, 0);
    for (size_t i = 0; i < n; ++i) {
      S[b[i]] += e[i];
      ++N[b[i]];
    }
    dvec sh(nblk), tc(nblk);
    std::vector<size_t> blocks(nblk);
    for (int j = 0; j < nblk; ++j) blocks[j] = j;
    const auto res = pmap(blocks, [&](size_t j) {
      const size_t a = edges[j], z = edges[j + 1];
      dvec Sb(nbf, 0.0);
      std::vector<int64_t> Nb(nbf, 0);
      for (size_t i = a; i < z; ++i) {
        Sb[b[i]] += e[i];
        ++Nb[b[i]];
      }
      dvec rest(nbf), own(nbf);
      for (int i = 0; i < nbf; ++i) {
        rest[i] = (S[i] - Sb[i]) / static_cast<double>(std::max<int64_t>(N[i] - Nb[i], 1));
        own[i] = Sb[i] / static_cast<double>(std::max<int64_t>(Nb[i], 1));
      }
      const double s_ = xcorr_shift(lead_edge(own, smooth), lead_edge(rest, smooth));
      return std::make_pair(s_ / nbf / rate, tt[(a + z) / 2]);
    });
    for (int j = 0; j < nblk; ++j) {
      sh[j] = res[j].first;
      tc[j] = res[j].second;
    }
    const double tm = np::mean(tc);
    dmat A(nblk, 2);
    dvec tcc(nblk);
    for (int j = 0; j < nblk; ++j) {
      tcc[j] = tc[j] - tm;
      A(j, 0) = tcc[j];
      A(j, 1) = 1.0;
    }
    const dvec coef = np::lstsq(A, sh).x;
    dvec r2(nblk), t2(nblk);
    for (int j = 0; j < nblk; ++j) {
      const double rr = sh[j] - (A(j, 0) * coef[0] + A(j, 1) * coef[1]);
      r2[j] = rr * rr;
      t2[j] = tcc[j] * tcc[j];
    }
    const double sig_slope = std::sqrt(np::sum(r2) / std::max(nblk - 2, 1) / np::sum(t2));
    const double slope = coef[0] * (nblk - 1) / nblk;
    sig = rate * sig_slope * (nblk - 1) / nblk;
    rate = rate / (1 + slope);
  }
  return {rate, sig};
}

SyncTip sync_tip_freq(const SpectralChannelizer& ch, double fs, double f_guess, double rate, double t0, double half) {
  const Band yb = slice_band(ch, f_guess, half, std::min(fs, 2.5 * half));
  const cdvec& y = yb.y;
  const double ofs = yb.fs;
  const double T = 1 / rate;
  const size_t nseg = static_cast<size_t>(3.3e-6 * ofs);
  const int64_t nl = static_cast<int64_t>(static_cast<double>(y.size()) / ofs / T) - 1;
  std::vector<int64_t> starts;
  for (int64_t i = 0; i < nl; ++i) {
    const int64_t s = static_cast<int64_t>((t0 + 0.7e-6 + static_cast<double>(i) * T) * ofs);
    if (s + static_cast<int64_t>(nseg) < static_cast<int64_t>(y.size())) starts.push_back(s);
  }
  const size_t nfft = 512;
  const dvec win = np::hanning(nseg);
  dvec acc(nfft, 0.0);
  cdvec seg(nfft);
  for (int64_t s : starts) {
    std::fill(seg.begin(), seg.end(), cd(0, 0));
    for (size_t j = 0; j < nseg; ++j) seg[j] = cd(y[s + j].real() * win[j], y[s + j].imag() * win[j]);
    const cdvec X = fft::fft(seg);
    for (size_t j = 0; j < nfft; ++j) acc[(j + nfft / 2) % nfft] += np::cabs2(X[j]);
  }
  dvec p(nfft);
  for (size_t j = 0; j < nfft; ++j) p[j] = acc[j] / static_cast<double>(starts.size());
  const size_t k = np::argmax(p);
  double d = 0.0;
  if (0 < k && k < nfft - 1)
    d = parabolic(std::log(p[k - 1] + 1e-30), std::log(p[k] + 1e-30), std::log(p[k + 1] + 1e-30));
  const double frk = static_cast<double>(static_cast<int64_t>(k) - static_cast<int64_t>(nfft / 2)) / nfft * ofs;
  SyncTip r;
  r.prom = p[k] / np::median(p);
  r.f = f_guess + frk + d * ofs / nfft;
  return r;
}

Timing refine_timing(const SpectralChannelizer& ch, double fs, double f_sync, double rate, double score,
                     double eclip) {
  auto env = [&](double fq, double hb, double& efs) {
    const Band yb = slice_band(ch, fq, hb, std::min(fs, 8e6));
    efs = yb.fs;
    dvec e(yb.y.size());
    for (size_t i = 0; i < e.size(); ++i) e[i] = np::cabs2(yb.y[i]);
    const double m = np::median(e);
    for (auto& v : e) v = std::min(v / m, eclip);
    return e;
  };
  Timing r;
  double efs = 0;
  dvec e = env(f_sync, 1.5e6, efs);
  rate = grid_line_rate(e, efs, rate);
  double t0 = sync_start(fold_prof(phase_bins(time_axis(e.size(), efs), rate, 1024), e.data(), e.size(), 1024),
                         rate, 1024);
  const SyncTip st = sync_tip_freq(ch, fs, f_sync, rate, t0);
  if (st.prom > 4 && std::fabs(st.f - f_sync) < 1.5e6) {
    f_sync = st.f;
    r.f_sync_np = true;
    e = env(f_sync, 1.0e6, efs);
    rate = grid_line_rate(e, efs, rate, 8.0);
  }
  const int nblk = score > 100 ? 8 : 4;
  auto [r2, sig] = drift_line_rate(e, efs, rate, nblk, 1024, 2, 16);
  const double sig_ppm = sig / rate * 1e6;
  if (std::isfinite(sig) && sig_ppm < 3.0 && std::fabs(r2 / rate - 1) * 1e6 < 8.0) rate = r2;
  else sig = rate * (score > 30 ? 10e-6 : 40e-6);
  const dvec prof = fold_prof(phase_bins(time_axis(e.size(), efs), rate, 1024), e.data(), e.size(), 1024);
  t0 = sync_start(prof, rate, 1024);
  r.rate = rate;
  r.sig = sig;
  r.t0 = t0;
  r.f_sync = f_sync;
  r.sync_prom = st.prom;
  r.e = std::move(e);
  r.efs = efs;
  return r;
}

std::vector<Vbi> field_timing(const dvec& e, double efs, double rate, double t0, const std::string& std_name) {
  const double T = 1 / rate;
  const int64_t nh = static_cast<int64_t>((static_cast<double>(e.size()) / efs - t0) / (T / 2));
  if (nh < 20) return {};
  const int64_t on0 = static_cast<int64_t>(6e-6 * efs), on1 = static_cast<int64_t>(26e-6 * efs);
  const int64_t g0 = static_cast<int64_t>(28.2e-6 * efs), g1 = static_cast<int64_t>(31.3e-6 * efs);
  dvec cs(e.size() + 1, 0.0);
  for (size_t i = 0; i < e.size(); ++i) cs[i + 1] = cs[i] + e[i];
  dvec hp;
  for (int64_t i = 0; i < nh; ++i) {
    const double hs = t0 + static_cast<double>(i) * T / 2;
    const int64_t a = static_cast<int64_t>(hs * efs);
    if (!(a + g1 < static_cast<int64_t>(e.size()))) continue;
    hp.push_back((cs[a + on1] - cs[a + on0]) / static_cast<double>(on1 - on0) -
                 (cs[a + g1] - cs[a + g0]) / static_cast<double>(g1 - g0));
  }
  const int run = std_name == "pal" ? 5 : 6;
  const dvec sm = np::convolve(hp, dvec(run, 1.0 / run), np::ConvMode::Same);
  const double med = np::median(sm);
  dvec ad(sm.size());
  for (size_t i = 0; i < sm.size(); ++i) ad[i] = std::fabs(sm[i] - med);
  const double mad = 1.4826 * np::median(ad) + 1e-30;
  dvec zz(sm.size());
  for (size_t i = 0; i < sm.size(); ++i) zz[i] = (sm[i] - med) / mad;
  const std::vector<int64_t> order = np::argsort_quick(zz);
  std::vector<int64_t> pk;
  for (auto it = order.rbegin(); it != order.rend(); ++it) {
    const int64_t k = *it;
    if (zz[k] < 8) break;
    bool far = true;
    for (int64_t j : pk) far &= std::llabs(k - j) > 100;
    if (far) pk.push_back(k);
  }
  if (pk.empty()) return {};
  const double lpf2 = nominal(std_name).lines;
  const int64_t k0 = pk[0];
  std::vector<Vbi> out;
  for (int64_t k : pk) {
    const double r = static_cast<double>(k - k0) / lpf2;
    if (std::fabs(r - std::nearbyint(r)) * lpf2 <= 4) out.push_back({k, zz[k], t0 + static_cast<double>(k) * T / 2});
  }
  std::sort(out.begin(), out.end(), [](const Vbi& a, const Vbi& b) { return a.k < b.k; });
  return out;
}

namespace {
// shared tail of spectrum_features / line_locked_spectrum: centre, edges, rms of seg over fb
template <class R>
void band_moments(const std::vector<R>& seg, const dvec& fb, double df, Spec& sp, R& tot_out) {
  const R tot = np::sum(seg) + static_cast<R>(1e-30);
  tot_out = tot;
  dvec sf(seg.size());
  for (size_t i = 0; i < seg.size(); ++i) sf[i] = static_cast<double>(seg[i]) * fb[i];
  sp.centre = np::sum(sf) / static_cast<double>(tot);
  std::vector<R> cdf = cumsum_seq(seg);
  for (auto& v : cdf) v = v / tot;
  const size_t n = fb.size();
  sp.e_lo = fb[std::min(searchsorted_d(cdf, 0.005), n - 1)] - df / 2;
  sp.e_hi = fb[std::min(searchsorted_d(cdf, 0.995), n - 1)] + df / 2;
  for (size_t i = 0; i < seg.size(); ++i) {
    const double d = fb[i] - sp.centre;
    sf[i] = static_cast<double>(seg[i]) * (d * d);
  }
  sp.rms = std::sqrt(std::max(np::sum(sf) / static_cast<double>(tot) - df * df / 12, 0.0));
}

// the band-growing loop shared by both spectrum estimators (on: bins above threshold)
void grow(const std::vector<uint8_t>& on, int k_sync, int a0, int b0, int& lo, int& hi) {
  auto any = [&](int a, int b) {
    for (int i = a; i < b; ++i)
      if (on[i]) return true;
    return false;
  };
  lo = hi = k_sync;
  while (lo > a0 && any(std::max(lo - 4, a0), lo)) --lo;
  while (hi < b0 && any(hi + 1, std::min(hi + 5, b0 + 1))) ++hi;
  while (lo < k_sync && !on[lo]) ++lo;
  while (hi > k_sync && !on[hi]) --hi;
}
}  // namespace

Spec spectrum_features(const fmat& P, const dvec& f, int k_sync, double max_half_hz) {
  const int nb = static_cast<int>(f.size());
  const double df = f[1] - f[0];
  const int W = static_cast<int>(std::ceil(max_half_hz / df));
  const int a0 = std::max(k_sync - W, 0), b0 = std::min(k_sync + W, nb - 1);
  fvec bp(P.rows);
  for (size_t i = 0; i < P.rows; ++i) bp[i] = np::sum(P.row(i) + a0, static_cast<size_t>(b0 - a0 + 1));
  const float med = np::median(bp);
  fvec ad(bp.size());
  for (size_t i = 0; i < bp.size(); ++i) ad[i] = std::fabs(bp[i] - med);
  const float mad = 1.4826f * np::median(ad);
  std::vector<uint8_t> keep(P.rows);
  size_t nkeep = 0;
  for (size_t i = 0; i < P.rows; ++i) nkeep += keep[i] = bp[i] <= med + 6.0f * mad;
  const fvec m = mean_rows(P, &keep);
  fmat Pk(nkeep, P.cols);
  for (size_t i = 0, j = 0; i < P.rows; ++i)
    if (keep[i]) std::copy(P.row(i), P.row(i) + P.cols, Pk.row(j++));
  fvec flo = np::percentile_axis0(Pk, 25);
  for (auto& v : flo) v = v / 0.2877f;
  const float p20 = np::percentile(flo, 20);
  bool n0_double = false;
  fvec n0f(nb, p20);
  dvec n0d;
  const double sig = 1 / std::sqrt(static_cast<double>(nkeep) / 2);
  int lo = k_sync, hi = k_sync;
  for (int it = 0; it < 2; ++it) {
    std::vector<uint8_t> on(nb);
    double thr;
    if (!n0_double) {
      fvec ex(nb);
      for (int i = 0; i < nb; ++i) ex[i] = m[i] / n0f[i] - 1.0f;
      const fvec exs = ndi::uniform_filter1d(ex, 3);
      const float pc = np::percentile(exs.data() + a0, b0 - a0 + 1, 90);
      thr = std::max(4 * sig / std::sqrt(3.0), static_cast<double>(0.03f * pc));
      for (int i = 0; i < nb; ++i) on[i] = static_cast<double>(exs[i]) > thr;
    } else {
      dvec ex(nb);
      for (int i = 0; i < nb; ++i) ex[i] = static_cast<double>(m[i]) / n0d[i] - 1;
      const dvec exs = ndi::uniform_filter1d(ex, 3);
      const double pc = np::percentile(exs.data() + a0, b0 - a0 + 1, 90);
      thr = std::max(4 * sig / std::sqrt(3.0), 0.03 * pc);
      for (int i = 0; i < nb; ++i) on[i] = exs[i] > thr;
    }
    grow(on, k_sync, a0, b0, lo, hi);
    std::vector<int> side;
    for (int k = lo - 4; k < lo - 1; ++k)
      if (k >= 0 && k < nb) side.push_back(k);
    for (int k = hi + 2; k < hi + 5; ++k)
      if (k >= 0 && k < nb) side.push_back(k);
    if (side.size() >= 2) {
      bool below = false, above = false, at = false;
      for (int k : side) {
        below |= k < k_sync;
        above |= k > k_sync;
        at |= k == k_sync;
      }
      const int deg = (static_cast<int>(below) + static_cast<int>(above) + static_cast<int>(at)) > 1 ? 1 : 0;
      dvec xs, ys;
      for (int k : side) {
        xs.push_back(f[k]);
        ys.push_back(static_cast<double>(flo[k]));
      }
      const dvec c = np::polyfit(xs, ys, deg);
      const double fl = static_cast<double>(0.5f * p20);
      n0d.assign(nb, 0.0);
      for (int i = 0; i < nb; ++i) n0d[i] = std::max(np::polyval(c, f[i]), fl);
      n0_double = true;
    }
  }
  Spec sp;
  const dvec fb(f.begin() + lo, f.begin() + hi + 1);
  if (!n0_double) {
    fvec exa(hi - lo + 1);
    for (int i = lo; i <= hi; ++i) exa[i - lo] = std::max(m[i] - n0f[i], 0.0f);
    float tot;
    band_moments(exa, fb, df, sp, tot);
    const float r = tot / np::sum(n0f.data() + lo, static_cast<size_t>(hi - lo + 1));
    sp.snr_db = r > 1e-6f ? static_cast<double>(10.0f * std::log10(r)) : 10 * std::log10(1e-6);
  } else {
    dvec exa(hi - lo + 1);
    for (int i = lo; i <= hi; ++i) exa[i - lo] = std::max(static_cast<double>(m[i]) - n0d[i], 0.0);
    double tot;
    band_moments(exa, fb, df, sp, tot);
    sp.snr_db = 10 * std::log10(std::max(tot / np::sum(n0d.data() + lo, static_cast<size_t>(hi - lo + 1)), 1e-6));
  }
  sp.lo = lo;
  sp.hi = hi;
  sp.edge_cut = lo <= 1 || hi >= nb - 2;
  return sp;
}

Spec line_locked_spectrum(const fmat& P, const dvec& t, const dvec& f, double rate, double t0, int k_sync,
                          double max_half_hz, double clip) {
  const double T = 1 / rate;
  const int nb = static_cast<int>(f.size());
  const double df = f[1] - f[0];
  const size_t nfr = P.rows;
  std::vector<uint8_t> core(nfr), act(nfr), mid(nfr);
  size_t ncore = 0;
  for (size_t i = 0; i < nfr; ++i) {
    const double ph = np::mod(t[i] - t0, T);
    core[i] = ph > 1.2e-6 && ph < 3.5e-6;
    act[i] = ph > 12.5e-6 && ph < T - 2.5e-6;
    mid[i] = ph > 20e-6 && ph < T - 10e-6;
    ncore += core[i];
  }
  fvec med = np::median_axis0(P);
  for (auto& v : med) v = v + 1e-30f;
  const float clipf = static_cast<float>(clip);
  fmat Pn(nfr, P.cols);
  for (size_t i = 0; i < nfr; ++i)
    for (size_t k = 0; k < P.cols; ++k) Pn(i, k) = std::min(P(i, k) / med[k], clipf);
  const fvec mc = mean_rows(Pn, &core), mm = mean_rows(Pn, &mid), ma = mean_rows(Pn, nullptr),
             mact = mean_rows(Pn, &act);
  fvec S(nb), Sv(nb), rel(nb);
  for (int k = 0; k < nb; ++k) {
    const float b = std::min(mc[k], mm[k]);   // np.minimum (NaN propagates)
    const float bb = (std::isnan(mc[k]) || std::isnan(mm[k])) ? std::numeric_limits<float>::quiet_NaN() : b;
    S[k] = (ma[k] - bb) * med[k];
    Sv[k] = (mact[k] - bb) * med[k];
    rel[k] = S[k] / med[k];
  }
  const int W = static_cast<int>(std::ceil(max_half_hz / df));
  const int a0 = std::max(k_sync - W, 0), b0 = std::min(k_sync + W, nb - 1);
  const fvec vc = var_rows(Pn, &core), va = var_rows(Pn, nullptr);
  const double half_all = static_cast<double>(nfr) / 2;
  dvec sig(nb);
  const bool strong = ncore >= 2;   // nc = core.sum() (numpy int64: strong scalar) -> float64 division
  for (int k = 0; k < nb; ++k) {
    const float tb = static_cast<float>(static_cast<double>(va[k]) / half_all);
    if (strong) sig[k] = std::sqrt(static_cast<double>(vc[k]) / (static_cast<double>(ncore) / 2) + static_cast<double>(tb)) + 1e-9;
    else sig[k] = static_cast<double>(std::sqrt(vc[k] / 1.0f + tb) + 1e-9f);
  }
  const fvec rs = ndi::uniform_filter1d(rel, 3);
  const float pc = np::percentile(rs.data() + a0, static_cast<size_t>(b0 - a0 + 1), 90);
  const double pcs = static_cast<double>(0.03f * pc);
  std::vector<uint8_t> on(nb);
  const double s3 = std::sqrt(3.0);
  for (int k = 0; k < nb; ++k) {
    const double s4 = strong ? 4 * sig[k] : static_cast<double>(4.0f * static_cast<float>(sig[k]));
    const double thr = std::max(s4 / s3, pcs);
    on[k] = static_cast<double>(rs[k]) > thr;
  }
  int lo, hi;
  grow(on, k_sync, a0, b0, lo, hi);
  if (hi - lo < 3) {
    dvec sf(b0 - a0 + 1);
    for (int i = a0; i <= b0; ++i) sf[i - a0] = static_cast<double>(S[i]) * (f[i] - f[k_sync]);
    const double sm = np::sum(sf);
    const double side = sm > 0 ? 1 : (sm < 0 ? -1 : 1);
    const int span = static_cast<int>(std::ceil(2.6 * 10e6 / df));
    if (side > 0) {
      lo = std::max(k_sync - 2, 0);
      hi = std::min(k_sync + span, nb - 1);
    } else {
      lo = std::max(k_sync - span, 0);
      hi = std::min(k_sync + 2, nb - 1);
    }
    fvec seg(S.begin() + lo, S.begin() + hi + 1);
    seg = ndi::uniform_filter1d(seg, 3);
    for (auto& v : seg) v = std::max(v, 0.0f);
    fvec cs = cumsum_seq(seg);
    const float last = cs.back() + 1e-30f;
    for (auto& v : cs) v = v / last;
    const int a = static_cast<int>(searchsorted_d(cs, 0.01)), z = static_cast<int>(searchsorted_d(cs, 0.99));
    const int lo0 = lo;
    lo = lo0 + a;
    hi = lo0 + std::min(z, hi - lo0);
  }
  Spec sp;
  fvec seg(hi - lo + 1);
  for (int i = lo; i <= hi; ++i) seg[i - lo] = std::max(S[i], 0.0f);
  dvec fb(f.begin() + lo, f.begin() + hi + 1);
  float tot;
  band_moments(seg, fb, df, sp, tot);
  fvec v(hi - lo + 1);
  for (int i = lo; i <= hi; ++i) v[i - lo] = std::max(Sv[i], 0.0f);
  dvec fv = fb;
  if (f[k_sync] > sp.centre) {
    std::reverse(v.begin(), v.end());
    std::reverse(fv.begin(), fv.end());
  }
  fvec cv = cumsum_seq(v);
  const float vs = np::sum(v) + 1e-30f;
  for (auto& x : cv) x = x / vs;
  sp.f_white = fv[std::min(searchsorted_d(cv, 0.995), fv.size() - 1)];
  fvec flo = np::percentile_axis0(P, 25);
  for (auto& x : flo) x = x / 0.2877f;
  const float n0 = np::percentile(flo, 20);
  const float r = tot / (n0 * static_cast<float>(seg.size()));
  sp.snr_db = r > 1e-6f ? static_cast<double>(10.0f * std::log10(r)) : 10 * std::log10(1e-6);
  sp.lo = lo;
  sp.hi = hi;
  sp.edge_cut = lo <= 1 || hi >= nb - 2;
  return sp;
}

Disc fm_discriminate(const SpectralChannelizer& ch, double f_off, double half_bw, double fs) {
  const Band yb = slice_band(ch, f_off, half_bw, std::min(fs, std::max(2.4 * half_bw, 4e6)));
  Disc r;
  r.fs = yb.fs;
  const size_t n = yb.y.size();
  r.d.resize(n > 0 ? n - 1 : 0);
  // numpy temporary elision: in `y[1:] * np.conj(y[:-1])` the conj() temporary, when it holds
  // >= 256 KiB, receives the product in place, i.e. numpy computes conj(prev) * cur; the
  // fma-based complex multiply is not symmetric in its last bit, so the order is mirrored.
  const bool elide = (n - 1) * sizeof(cd) >= 256 * 1024;
  for (size_t i = 1; i < n; ++i) {
    const cd z = elide ? np::cmul(std::conj(yb.y[i - 1]), yb.y[i]) : np::cmul(yb.y[i], std::conj(yb.y[i - 1]));
    r.d[i - 1] = std::atan2(z.imag(), z.real()) * yb.fs / (2 * kPi);
  }
  return r;
}

// ============================================================== SECAM (python secam_test; research, ITU-R BT.470)
// No burst: the back porch carries the FM subcarrier at its rest frequency, Db 4.250 MHz / Dr 4.40625 MHz on
// alternate lines -> per-line r = (e_Db - e_Dr) / (e_Db + e_Dr) alternates; rest frequencies on a 5 kHz grid.
Json secam_test(const dvec& d, double ofs, double rate, double t0) {
  constexpr double W0 = 5.8e-6, W1 = 10.0e-6, kAlt = 0.25;
  const double T = 1.0 / rate;
  Json out = Json::object();
  const int64_t k0 = static_cast<int64_t>(std::ceil((-t0 - W0) / T));
  const int64_t k1 = static_cast<int64_t>(std::floor((static_cast<double>(d.size()) / ofs - t0 - W1) / T));
  const int64_t m = static_cast<int64_t>((W1 - W0) * ofs);
  if (k1 - k0 + 1 < 32 || m < 8) return out;
  std::vector<int64_t> i0, lines;
  for (int64_t k = k0; k <= k1; ++k) {
    const int64_t v = static_cast<int64_t>(std::ceil((t0 + static_cast<double>(k) * T + W0) * ofs));
    if (v >= 0 && v + m <= static_cast<int64_t>(d.size())) i0.push_back(v), lines.push_back(k);
  }
  if (i0.size() < 32) return out;
  const size_t nl = i0.size(), M = static_cast<size_t>(m);
  std::vector<dvec> seg(nl, dvec(M));
  for (size_t r = 0; r < nl; ++r) {
    for (size_t j = 0; j < M; ++j) seg[r][j] = d[static_cast<size_t>(i0[r]) + j];
    const double mu = np::mean(seg[r]);
    for (auto& v : seg[r]) v -= mu;
  }
  // energies at frs: |sum_n seg[n] e^{-j 2 pi f n / ofs}|^2 per line (the line's start phase drops out)
  auto energies = [&](const dvec& frs) {
    std::vector<dvec> E(nl, dvec(frs.size()));
    for (size_t q = 0; q < frs.size(); ++q) {
      dvec c(M), s(M);
      for (size_t j = 0; j < M; ++j) {
        const double ph = 2 * kPi * (static_cast<double>(j) / ofs) * frs[q];
        c[j] = std::cos(ph), s[j] = std::sin(ph);
      }
      for (size_t r = 0; r < nl; ++r) {
        double a = 0, b = 0;
        for (size_t j = 0; j < M; ++j) a += seg[r][j] * c[j], b += seg[r][j] * s[j];
        E[r][q] = a * a + b * b;
      }
    }
    return E;
  };
  const auto E0 = energies({4.25e6, 4.40625e6, 2.6e6, 3.1e6, 5.0e6, 5.4e6});
  dvec a(nl), alt(nl);
  for (size_t r = 0; r < nl; ++r) {
    const double eb = E0[r][0], er = E0[r][1];
    alt[r] = pymod(lines[r], int64_t(2)) == 0 ? 1.0 : -1.0;
    a[r] = (eb - er) / (eb + er + 1e-30) * alt[r];
  }
  const double mu = np::mean(a);
  double var = 0;
  for (double v : a) var += (v - mu) * (v - mu);
  const double sd = std::sqrt(var / static_cast<double>(nl));
  const double tstat = std::fabs(mu) / (sd / std::sqrt(static_cast<double>(nl)) + 1e-12);
  dvec rm(4);
  for (size_t q = 0; q < 4; ++q) {
    dvec col(nl);
    for (size_t r = 0; r < nl; ++r) col[r] = E0[r][q + 2];
    rm[q] = np::mean(col);
  }
  const double ref = np::median(rm) + 1e-30;
  dvec ebr(nl);
  for (size_t r = 0; r < nl; ++r) ebr[r] = E0[r][0] + E0[r][1];
  const double prom = 10 * std::log10(np::mean(ebr) / 2 / ref);
  out["alternation"] = pyr(std::fabs(mu), 3);
  out["t_stat"] = pyr(tstat, 1);
  out["prominence_db"] = pyr(prom, 1);
  if (std::fabs(mu) >= kAlt) {
    dvec grid;
    for (int k = 0; 4.10e6 + k * 5e3 < 4.56e6 + 1; ++k) grid.push_back(4.10e6 + k * 5e3);
    const auto E = energies(grid);
    double bb = -1, br = -1;
    size_t jb = 0, jr = 0;
    for (size_t q = 0; q < grid.size(); ++q) {
      double sb = 0, sr = 0;
      size_t nb = 0, nr = 0;
      for (size_t r = 0; r < nl; ++r) {
        const bool pb = mu > 0 ? alt[r] > 0 : alt[r] < 0;
        if (pb) sb += E[r][q], ++nb; else sr += E[r][q], ++nr;
      }
      sb /= std::max<size_t>(nb, 1), sr /= std::max<size_t>(nr, 1);
      if (sb > bb) bb = sb, jb = q;
      if (sr > br) br = sr, jr = q;
    }
    out["db_hz"] = grid[jb];
    out["dr_hz"] = grid[jr];
  }
  return out;
}

// ============================================================== colour standard (python chroma_standard; research:
// research/drone_transmitters/01_analog_video.md). Line-locked burst phasors at each standard's fsc of the line family.
Json chroma_standard(const dvec& d, double ofs, double rate, double t0, const std::string& std_name) {
  struct Cand { const char* name; double k, off; };
  static const Cand kChroma[] = {{"ntsc", 227.5, 0.0}, {"pal_m", 227.25, 0.0}, {"pal", 283.75, 25.0}, {"pal_n", 229.25, 25.0}};
  constexpr double W0 = 5.7e-6, W1 = 7.5e-6, kSearch = 0.06, kCoh = 0.35, kProm = 6.0;
  const double T = 1.0 / rate;
  const int64_t k0 = static_cast<int64_t>(std::ceil((-t0 - W0) / T));
  const int64_t k1 = static_cast<int64_t>(std::floor((static_cast<double>(d.size()) / ofs - t0 - W1) / T));
  const int64_t m = static_cast<int64_t>((W1 - W0) * ofs);
  Json out = Json::object();
  if (k1 - k0 + 1 < 32 || m < 4) return out;
  std::vector<int64_t> i0;
  for (int64_t k = k0; k <= k1; ++k) {
    const int64_t v = static_cast<int64_t>(std::ceil((t0 + static_cast<double>(k) * T + W0) * ofs));
    if (v >= 0 && v + m <= static_cast<int64_t>(d.size())) i0.push_back(v);
  }
  if (i0.size() < 32) return out;
  const size_t nl = i0.size();
  const size_t nfft = size_t(1) << static_cast<int>(std::ceil(std::log2(static_cast<double>(nl * 8))));
  const size_t kmax = static_cast<size_t>(std::ceil(kSearch * static_cast<double>(nfft)));
  std::vector<size_t> idx;
  for (size_t i = 0; i <= kmax; ++i) idx.push_back(i);
  for (size_t i = nfft - kmax; i < nfft; ++i) idx.push_back(i);
  // per line, remove a straight line (porch level + sync-edge tail; python: seg - mean - outer(slope, u))
  const size_t M = static_cast<size_t>(m);
  dvec u(M), uu(M);
  for (size_t j = 0; j < M; ++j) u[j] = static_cast<double>(j) - static_cast<double>(m - 1) / 2, uu[j] = u[j] * u[j];
  const double suu = np::sum(uu);
  std::vector<dvec> seg(nl, dvec(M));
  {
    dvec row(M), su(M);
    for (size_t r = 0; r < nl; ++r) {
      for (size_t j = 0; j < M; ++j) row[j] = d[static_cast<size_t>(i0[r]) + j], su[j] = row[j] * u[j];
      const double mu = np::mean(row), sl = np::sum(su) / suu;
      for (size_t j = 0; j < M; ++j) seg[r][j] = (row[j] - mu) - sl * u[j];
    }
  }
  auto phasors = [&](double fr) {
    const double w = (-2.0 * kPi) * fr;
    cdvec zl(nl), row(M);
    for (size_t r = 0; r < nl; ++r) {
      for (size_t j = 0; j < M; ++j) {
        const int64_t g = i0[r] + static_cast<int64_t>(j);
        const double th = w * (static_cast<double>(g) / ofs);
        row[j] = cd(seg[r][j] * std::cos(th), seg[r][j] * std::sin(th));
      }
      zl[r] = np::sum(row);
    }
    return zl;
  };
  auto mean_pow = [](const cdvec& z) {
    dvec p(z.size());
    for (size_t i = 0; i < z.size(); ++i) p[i] = np::cabs2(z[i]);
    return np::mean(p);
  };
  // references below 0.45 ofs and >= 1 MHz from every candidate subcarrier (python CHROMA_REF_STEP / _GAP)
  dvec refs;
  {
    dvec fscs;
    for (const auto& c : kChroma) fscs.push_back(c.k * rate + c.off * rate / 15625.0);
    const double stop = 0.45 * ofs;
    const int64_t nstep = static_cast<int64_t>(std::ceil((stop - 1.5e6) / 0.25e6));   // np.arange(1.5e6, stop, 0.25e6)
    for (int64_t i = 0; i < nstep; ++i) {
      const double fr = 1.5e6 + static_cast<double>(i) * 0.25e6;
      double dmin = 1e300;
      for (double f0 : fscs) dmin = std::min(dmin, std::fabs(fr - f0));
      if (dmin >= 1.0e6) refs.push_back(mean_pow(phasors(fr)));
    }
    if (refs.empty()) refs.push_back(mean_pow(phasors(2.0e6)));
  }
  const double ref = np::median(refs) + 1e-30;
  const bool fam525 = std_name == "ntsc";
  Json cands = Json::object();
  for (const auto& c : kChroma) {
    const std::string name = c.name;
    const bool in525 = name == "ntsc" || name == "pal_m";
    if (in525 != fam525) continue;
    const double fsc = c.k * rate + c.off * rate / 15625.0;
    const cdvec zl = phasors(fsc);
    dvec pw(nl);
    for (size_t i = 0; i < nl; ++i) pw[i] = np::cabs2(zl[i]);
    const double e = np::sum(pw) + 1e-30;
    cdvec za(nl);
    for (size_t i = 0; i < nl; ++i) za[i] = (i % 2 ? -1.0 : 1.0) * zl[i];
    const cdvec Z = fft::fft(zl, nfft), Zh = fft::fft(za, nfft);
    size_t jb = 0;
    double cb = -1;
    for (size_t q = 0; q < idx.size(); ++q) {
      const double v = (np::cabs2(Z[idx[q]]) + np::cabs2(Zh[idx[q]])) / (static_cast<double>(nl) * e);
      if (v > cb) cb = v, jb = q;
    }
    const size_t ij = idx[jb];
    const double nu = ij <= nfft / 2 ? static_cast<double>(ij) / nfft : static_cast<double>(ij) / nfft - 1;
    const double swing = std::atan2(std::abs(Zh[ij]), std::abs(Z[ij])) * 180.0 / kPi;
    cands[name] = {{"coherence", pyr(cb, 3)},
                   {"prominence_db", pyr(10 * std::log10(e / static_cast<double>(nl) / ref), 1)},
                   {"swing_deg", pyr(swing, 1)},
                   {"fsc_hz", npr(fsc + nu * rate, 1)}};
  }
  out["chroma_candidates"] = cands;
  std::string best;
  double bp = 0, bc = 0;
  for (auto it = cands.begin(); it != cands.end(); ++it) {
    const Json& v = it.value();
    const double coh = v["coherence"].get<double>(), prom = v["prominence_db"].get<double>(), sw = v["swing_deg"].get<double>();
    const bool swing_ok = it.key() != "ntsc" ? sw > 25.0 : sw < 20.0;
    if (coh >= kCoh && prom >= kProm && swing_ok && (best.empty() || prom > bp || (prom == bp && coh > bc)))
      best = it.key(), bp = prom, bc = coh;
  }
  if (best.empty() && std_name != "ntsc") {
    const Json sec = secam_test(d, ofs, rate, t0);
    if (!sec.empty()) {
      out["secam_test"] = sec;
      if (sec["alternation"].get<double>() >= 0.25 && sec["t_stat"].get<double>() >= 6.0 &&
          sec["prominence_db"].get<double>() >= kProm) {
        out["color"] = true;
        out["color_standard"] = "secam";
        out["secam_db_hz"] = sec["db_hz"];
        out["secam_dr_hz"] = sec["dr_hz"];
        return out;
      }
    }
  }
  if (!best.empty()) {
    const Json& v = cands[best];
    out["color"] = true;
    out["color_standard"] = best;
    out["chroma_fsc_hz"] = v["fsc_hz"];
    out["burst_swing_deg"] = v["swing_deg"];
    out["chroma_coherence"] = v["coherence"];
  } else {
    out["color"] = false;
  }
  return out;
}

// python audio_features(d, ofs, fa): the audio FM subcarrier sliced +-100 kHz and FM-discriminated (research, A5)
Json audio_features(const dvec& d, double ofs, double fa) {
  constexpr double half = 100e3;
  const size_t n = d.size();
  const double mu = np::mean(d);
  cdvec xin(n);
  for (size_t i = 0; i < n; ++i) xin[i] = cd(d[i] - mu, 0.0);
  const cdvec X = fft::fft(xin);
  const int64_t k0 = py_round_int(fa / ofs * static_cast<double>(n));
  const int64_t kk = static_cast<int64_t>(half / ofs * static_cast<double>(n));
  const int64_t M = int64_t(1) << static_cast<int>(std::ceil(std::log2(static_cast<double>(2 * kk + 1))));
  cdvec Y(static_cast<size_t>(M), cd(0, 0));
  for (int64_t m = -kk; m <= kk; ++m) Y[static_cast<size_t>(pymod(m, M))] = X[static_cast<size_t>(pymod(k0 + m, static_cast<int64_t>(n)))];
  const double fs2 = static_cast<double>(M) * ofs / static_cast<double>(n);
  const cdvec z = fft::ifft(Y);
  dvec a(z.size() - 1);
  for (size_t i = 0; i + 1 < z.size(); ++i) a[i] = std::arg(z[i + 1] * std::conj(z[i])) * fs2 / (2 * kPi);
  const double off = np::median(a);
  const double am = np::mean(a);
  for (auto& v : a) v -= am;
  const size_t la = a.size();
  const size_t nps = std::min(la, size_t(1) << static_cast<int>(std::floor(std::log2(static_cast<double>(la)))));
  const auto W = sig::welch(a, fs2, nps);
  dvec pb;
  for (size_t i = 0; i < W.f.size(); ++i)
    if (W.f[i] >= 50 && W.f[i] <= 8000) pb.push_back(W.Pxx[i]);
  Json out = Json::object();
  if (pb.empty()) return out;
  dvec lg(pb.size());
  for (size_t i = 0; i < pb.size(); ++i) lg[i] = std::log(pb[i] + 1e-30);
  const double flat = std::exp(np::mean(lg)) / (np::mean(pb) + 1e-30);
  const size_t nf = 8 * (size_t(1) << static_cast<int>(std::ceil(std::log2(static_cast<double>(la)))));
  const dvec hw = np::hanning(la);
  dvec aw(la);
  for (size_t i = 0; i < la; ++i) aw[i] = a[i] * hw[i];
  const cdvec F = fft::rfft(aw, nf);
  size_t j = 0;
  double best = -2;
  for (size_t i = 0; i < F.size(); ++i) {
    const double f2 = static_cast<double>(i) * fs2 / static_cast<double>(nf);
    const double v = (f2 >= 50 && f2 <= 8000) ? np::cabs2(F[i]) : -1.0;
    if (v > best) best = v, j = i;
  }
  const double l0 = std::log(np::cabs2(F[j - 1]) + 1e-30), l1 = std::log(np::cabs2(F[j]) + 1e-30),
               l2 = std::log(np::cabs2(F[j + 1]) + 1e-30);
  const double dj = (l0 - 2 * l1 + l2) < 0 ? 0.5 * (l0 - l2) / (l0 - 2 * l1 + l2) : 0.0;
  const double fpk = (static_cast<double>(j) + dj) * fs2 / static_cast<double>(nf);
  size_t jp = 0;
  for (size_t i = 1; i < W.f.size(); ++i)
    if (std::fabs(W.f[i] - fpk) < std::fabs(W.f[jp] - fpk)) jp = i;
  const dvec lp = sig::sosfiltfilt(sig::butter_sos(4, 15e3, sig::BType::Lowpass, fs2), a);
  const double lmu = np::mean(lp);
  double var = 0;
  for (double v : lp) var += (v - lmu) * (v - lmu);
  out["audio_subcarrier_hz"] = npr(fa + off + (static_cast<double>(k0) / static_cast<double>(n) * ofs - fa), 0);
  out["audio_dev_rms_khz"] = pyr(std::sqrt(var / static_cast<double>(lp.size())) / 1e3, 2);
  out["audio_peak_hz"] = pyr(fpk, 1);
  out["audio_peak_db"] = pyr(10 * std::log10(W.Pxx[jp] / np::median(pb)), 1);
  out["audio_flatness"] = pyr(flat, 3);
  return out;
}

Json demod_features(const dvec& d, double ofs, double rate, double t0, const std::string& std_name, int polarity) {
  Json out = Json::object();
  const double T = 1 / rate;
  const size_t n = d.size();
  const int nbf = 256;
  dvec tt(n), ph(n);
  std::vector<int64_t> b(n);
  for (size_t i = 0; i < n; ++i) {
    tt[i] = static_cast<double>(i) / ofs;
    ph[i] = np::mod(tt[i] - t0, T);
    b[i] = static_cast<int64_t>(ph[i] / T * nbf);
  }
  const dvec line = fold_prof(b, d.data(), n, nbf);
  dvec s1, s2;
  const double a_us = std_name == "pal" ? 8.5 : 8.0, b_us = std_name == "pal" ? 10.0 : 9.2;
  for (int j = 0; j < nbf; ++j) {
    const double us = static_cast<double>(j) / nbf * T * 1e6;
    if (us > 1.0 && us < 3.7) s1.push_back(line[j]);
    if (us > a_us && us < b_us) s2.push_back(line[j]);
  }
  const double sync = np::median(s1), blank = np::median(s2);
  out["sync_level_mhz"] = sync / 1e6;
  out["blank_level_mhz"] = blank / 1e6;
  dvec sm = ndi::uniform_filter1d(d, std::max<size_t>(static_cast<size_t>(ofs * 0.3e-6), 1));
  for (auto& v : sm) v = v * polarity;
  int64_t lmin = std::numeric_limits<int64_t>::max(), lmax = std::numeric_limits<int64_t>::min();
  std::vector<int64_t> ln(n);
  for (size_t i = 0; i < n; ++i) {
    ln[i] = static_cast<int64_t>(std::floor((tt[i] - t0) / T));
    if (ph[i] > 12e-6 && ph[i] < T - 3e-6) {
      lmin = std::min(lmin, ln[i]);
      lmax = std::max(lmax, ln[i]);
    }
  }
  dvec mx(static_cast<size_t>(lmax - lmin + 1), -kInf);
  for (size_t i = 0; i < n; ++i)
    if (ph[i] > 12e-6 && ph[i] < T - 3e-6) {
      double& m = mx[ln[i] - lmin];
      if (sm[i] > m || std::isnan(sm[i])) m = sm[i];
    }
  dvec fin;
  for (double v : mx)
    if (std::isfinite(v)) fin.push_back(v);
  const double white = polarity * np::percentile(fin, 90);
  out["white_level_mhz"] = white / 1e6;
  out["deviation_demod_mhz"] = std::fabs(white - sync) / 2e6;
  // colour burst: gated back porch
  dvec seg, tb;
  for (size_t i = 0; i < n; ++i)
    if (ph[i] > 5.6e-6 && ph[i] < 7.6e-6) {
      seg.push_back(d[i]);
      tb.push_back(tt[i]);
    }
  auto tone = [&](double fr) {
    const double w = (-2.0 * kPi) * fr;
    std::vector<int64_t> l(seg.size());
    int64_t l0 = std::numeric_limits<int64_t>::max();
    for (size_t i = 0; i < seg.size(); ++i) {
      l[i] = static_cast<int64_t>(std::floor((tb[i] - t0) / T));
      l0 = std::min(l0, l[i]);
    }
    int64_t l1 = 0;
    for (auto& v : l) l1 = std::max(l1, v -= l0);
    dvec re(l1 + 1, 0.0), im(l1 + 1, 0.0);
    for (size_t i = 0; i < seg.size(); ++i) {
      const double th = w * tb[i];
      re[l[i]] += seg[i] * std::cos(th);
      im[l[i]] += seg[i] * std::sin(th);
    }
    dvec p(re.size());
    for (size_t i = 0; i < re.size(); ++i) p[i] = np::cabs2(cd(re[i], im[i]));
    return np::mean(p);
  };
  if (ofs > 9.5e6 && !seg.empty()) {
    const double e_pal = tone(4.43361875e6), e_ntsc = tone(3.579545e6);
    dvec refs;
    for (double fr : {2.6e6, 3.1e6, 4.0e6, 5.0e6, 5.4e6}) refs.push_back(tone(fr));
    const double ref = np::median(refs);
    const double r_pal = 10 * std::log10(e_pal / ref), r_ntsc = 10 * std::log10(e_ntsc / ref);
    out["burst_db"] = {{"f4_43", npr(r_pal, 1)}, {"f3_58", npr(r_ntsc, 1)}};
    if (std::max(r_pal, r_ntsc) > 6) out["color_burst_mhz"] = r_pal > r_ntsc ? 4.434 : 3.580;
  }
  if (ofs > 9.5e6) {
    // both H-sync edges (the sync time can lock on the trailing edge at low SNR); keep the reading in which a
    // colour standard passes, else the stronger coherence (python demod_features)
    std::vector<Json> rs;
    for (double sh : {0.0, -kSyncS}) {
      Json r = chroma_standard(d, ofs, rate, t0 + sh, std_name);
      if (!r.empty()) rs.push_back(std::move(r));
    }
    if (!rs.empty()) {
      auto key = [](const Json& r) {
        return std::make_pair(r.value("color", false) ? 1 : 0, r.contains("chroma_coherence") ? r["chroma_coherence"].get<double>() : 0.0);
      };
      size_t bi = 0;
      for (size_t i = 1; i < rs.size(); ++i)
        if (key(rs[i]) > key(rs[bi])) bi = i;
      Json r = rs[bi];
      if (bi == rs.size() - 1 && rs.size() == 2 && r.value("color", false))
        r["chroma_gate_note"] = "burst found 4.7 us before the measured sync time (trailing-edge lock)";
      for (auto it = r.begin(); it != r.end(); ++it) out[it.key()] = it.value();
    }
  }
  if (ofs > 15e6) {
    const auto W = sig::welch(d, ofs, 8192);
    bool have = false;
    double prom = 0, fa_best = 0;
    for (double fa : {5.5e6, 6.0e6, 6.5e6, 7.0e6}) {
      if (fa > ofs / 2 * 0.9) continue;
      double pmax = -kInf;
      dvec ring;
      for (size_t i = 0; i < W.f.size(); ++i) {
        const double a = std::fabs(W.f[i] - fa);
        if (a < 100e3) pmax = std::max(pmax, W.Pxx[i]);
        if (a > 300e3 && a < 700e3) ring.push_back(W.Pxx[i]);
      }
      const double v = 10 * std::log10(pmax / np::median(ring));
      if (!have || v > prom || (v == prom && fa > fa_best)) {
        prom = v;
        fa_best = fa;
        have = true;
      }
    }
    if (have) {
      out["audio_prominence_db"] = npr(prom, 1);
      if (prom > 10) {
        out["audio_subcarrier_mhz"] = fa_best / 1e6;
        if (prom >= 15.0) {
          const Json af = audio_features(d, ofs, fa_best);
          for (auto it = af.begin(); it != af.end(); ++it) out[it.key()] = it.value();
        }
      }
    }
  }
  return out;
}

// ============================================================== main entry
namespace {

// python _capture_locked: LO leakage modulated by a periodic process of the capture chain -- line rate within
// 2 ppm of fs / N, "sync tip" within 100 kHz of the LO, no measured field sync
bool capture_locked(const Emission& em, double fs, double fc) {
  const Json& F = em.features;
  const double r = F["line_rate_hz"].get<double>();
  const double lock_ppm = (r * static_cast<double>(py_round_int(fs / r)) / fs - 1) * 1e6;
  const bool vbi = F.contains("field_rate_source") && F["field_rate_source"] == "measured_vbi";
  return std::fabs(lock_ppm) < 2.0 && std::fabs(F["sync_tip_mhz"].get<double>() * 1e6 - fc) < 100e3 && !vbi;
}

bool is_part_of(const ScanBest& best, const dvec& f, double fc, const std::vector<Emission>& found) {
  const double T = 1 / best.rate;
  const double ph_us = best.b0 / static_cast<double>(kNB) * T * 1e6;
  const double fk = (fc + f[best.k0]) / 1e6;
  for (const auto& o : found) {
    const Json& F = o.features;
    if (F["standard"].get<std::string>() != best.std_name ||
        std::fabs(F["line_rate_ppm"].get<double>() - best.ppm) > 150)
      continue;
    const double d = np::mod(ph_us - F["sync_phase_us"].get<double>(), T * 1e6);
    const double fsm = (f[1] - f[0]) * static_cast<double>(f.size()) / 1e6;
    const double dsync = std::fabs(np::mod((fk - F["sync_tip_mhz"].get<double>()) + fsm / 2, fsm) - fsm / 2);
    const bool near = (o.f_lo_mhz - 5 <= fk && fk <= o.f_hi_mhz + 5) || dsync < 25;
    if ((d < 14 || d > T * 1e6 - 4) && near) return true;
  }
  return false;
}

Emission characterise(double fs, double fc, const SpectralChannelizer& ch, const TfPower& tf, const ScanBest& best,
                      std::pair<int, int> cols, double dur_us) {
  const fmat& P = tf.P;
  const dvec& t = tf.t;
  const dvec& f = tf.f;
  const std::string std_name = best.std_name;
  double rate = best.rate;
  const double df = f[1] - f[0];
  const int nbins = static_cast<int>(f.size());
  fvec pm = np::median_axis0(P);
  for (auto& v : pm) v = v + 1e-30f;
  fmat Pn(P.rows, P.cols);
  for (size_t i = 0; i < P.rows; ++i)
    for (size_t k = 0; k < P.cols; ++k) Pn(i, k) = std::min(P(i, k) / pm[k], 50.0f);
  const Fold Fr = fold(Pn, t, rate);
  const int nb_box = static_cast<int>(py_round_int(kSyncS * rate * kNB));
  const dmat Fb = ndi::convolve1d(Fr.F, dvec(nb_box, 1.0 / nb_box), 0, ndi::Mode::Wrap);
  const dvec mFr = np::median_axis0(Fr.F), mFb = np::median_axis0(Fb);
  dmat exc(kNB, nbins), ad(kNB, nbins);
  for (int b = 0; b < kNB; ++b)
    for (int k = 0; k < nbins; ++k) {
      exc(b, k) = Fb(b, k) - mFr[k];
      ad(b, k) = std::fabs(Fb(b, k) - mFb[k]);
    }
  const dvec mad = np::median_axis0(ad);
  const int k_lo = std::max(cols.first - 3, 0), k_hi = std::min(cols.second + 3, nbins - 1);
  dvec zmax(k_hi - k_lo + 1, -kInf);
  for (int k = k_lo; k <= k_hi; ++k)
    for (int b = 0; b < kNB; ++b) zmax[k - k_lo] = std::max(zmax[k - k_lo], exc(b, k) / (1.4826 * mad[k] + 1e-12));
  const double zm = *std::max_element(zmax.begin(), zmax.end());
  std::vector<int> cand;
  for (int k = k_lo; k <= k_hi; ++k)
    if (zmax[k - k_lo] >= 0.4 * zm) cand.push_back(k);
  std::vector<int64_t> bph(t.size());
  for (size_t i = 0; i < t.size(); ++i) bph[i] = static_cast<int64_t>(mod1(t[i] * rate) * kNB);
  const int half = std::max(nb_box / 2 - 1, 1);
  auto circ = [](int64_t d) { return imod(d + kNB / 2, kNB) - kNB / 2; };
  std::vector<std::tuple<double, int, int>> cinfo;
  for (int k : cand) {
    size_t bk = 0;
    for (int b = 1; b < kNB; ++b)
      if (exc(b, k) > exc(bk, k)) bk = b;
    std::vector<uint8_t> ins(t.size());
    for (size_t i = 0; i < t.size(); ++i) ins[i] = std::llabs(circ(bph[i] - static_cast<int64_t>(bk))) <= half;
    cinfo.emplace_back(static_cast<double>(median_col(P, k, &ins)), k, static_cast<int>(bk));
  }
  const auto [lmax, kref, bref] = *std::max_element(cinfo.begin(), cinfo.end());
  std::vector<std::tuple<int64_t, int, int>> loud;
  const int64_t lim = static_cast<int64_t>(8e-6 * rate * kNB);
  for (const auto& [lv, k, bk] : cinfo) {
    const int64_t c = circ(bk - bref);
    if (lv >= 0.5 * lmax && std::llabs(c) <= lim) loud.emplace_back(c, k, bk);
  }
  const auto [c0_, k0, b0] = *std::min_element(loud.begin(), loud.end());
  double wsum = 0, wf = 0;
  dvec wv, wfv;
  for (int k = std::max(k0 - 1, 0); k < std::min(k0 + 2, nbins); ++k) {
    const double w = std::max(exc(b0, k), 0.0);
    wv.push_back(w);
    wfv.push_back(w * f[k]);
  }
  wsum = np::sum(wv);
  wf = np::sum(wfv);
  double f_sync = wsum > 0 ? wf / (wsum + 1e-12) : f[k0];
  std::vector<uint8_t> insync(t.size());
  for (size_t i = 0; i < t.size(); ++i) insync[i] = std::llabs(circ(bph[i] - b0)) <= half;
  const float lvl = median_col(P, k0, &insync) / (median_col(P, k0, nullptr) + 1e-30f);
  const double clip_v = std::max(6.0, 2.0 * static_cast<double>(lvl));

  const Timing tm = refine_timing(ch, fs, f_sync, rate, best.score, clip_v);
  rate = tm.rate;
  const double sig_rate = tm.sig, t0 = tm.t0;
  f_sync = tm.f_sync;
  const bool fsnp = tm.f_sync_np;
  const Nominal& nom = nominal(std_name);
  const double ppm = (rate / nom.fl - 1) * 1e6;
  int k_s = 0;
  for (int k = 1; k < nbins; ++k)
    if (std::fabs(f[k] - f_sync) < std::fabs(f[k_s] - f_sync)) k_s = k;
  Spec sp = line_locked_spectrum(P, t, f, rate, t0, k_s, 25e6, clip_v);
  sp.source = "line_locked";
  if (sp.hi - sp.lo < 3) {
    Spec sp2 = spectrum_features(P, f, k_s);
    sp2.f_white = f_sync < sp2.centre ? sp2.e_hi : sp2.e_lo;
    sp2.source = "psd_excess";
    sp = sp2;
  }
  const bool spec_ok = sp.hi - sp.lo >= 3;
  const double centre = sp.centre, e_lo = sp.e_lo, e_hi = sp.e_hi, rms_dev = sp.rms, snr_db = sp.snr_db;
  const int polarity = f_sync < centre ? 1 : -1;
  const double dev = std::max(std::fabs(sp.f_white - f_sync) - kDevCorrHz, 0.0) / 2;
  auto fr = [&](double x, int n) { return fsnp ? npr(x, n) : pyr(x, n); };   // f_sync-derived values
  Json F = Json::object();
  F["standard"] = std_name;
  F["line_rate_hz"] = npr(rate, 4);
  F["line_rate_ppm"] = npr(ppm, 2);
  F["line_rate_ppm_sigma"] = std::isfinite(sig_rate) ? Json(npr(sig_rate / nom.fl * 1e6, 3)) : Json(nullptr);
  F["line_period_us"] = npr(1e6 / rate, 6);
  F["sync_phase_us"] = npr(t0 * 1e6, 3);
  F["sync_tip_mhz"] = fr((fc + f_sync) / 1e6, 3);
  F["sync_offset_mhz"] = fr((f_sync - centre) / 1e6, 3);
  F["polarity"] = polarity > 0 ? "sync_low" : "sync_high";
  F["rms_dev_mhz"] = pyr(rms_dev / 1e6, 3);
  F["deviation_mhz"] = fr(dev / 1e6, 2);
  F["bw_mhz"] = pyr((e_hi - e_lo) / 1e6, 2);
  F["snr_db"] = pyr(snr_db, 1);
  F["score"] = pyr(best.score, 1);
  F["band_edge_truncated"] = sp.edge_cut;
  F["scan_bin_mhz"] = npr(df / 1e6, 3);
  F["ppm_search_best"] = pyr(best.ppm, 1);
  const std::vector<Vbi> vbi = field_timing(tm.e, tm.efs, rate, t0, std_name);
  Json vu = Json::array();
  for (const auto& v : vbi) vu.push_back(npr(v.t * 1e6, 1));
  F["vbi_us"] = vu;
  const double lpf_nom = nom.lines / 2.0;
  double lpf;
  if (vbi.size() >= 2) {
    const int64_t dk = vbi[1].k - vbi[0].k;
    const int64_t nf = std::max<int64_t>(py_round_int(static_cast<double>(dk) / nom.lines), 1);
    lpf = static_cast<double>(dk) / 2 / static_cast<double>(nf);
    F["lines_per_field"] = lpf;
    F["scan"] = std::fabs(lpf - std::nearbyint(lpf)) > 0.25 ? "interlaced" : "progressive";
    F["field_rate_source"] = "measured_vbi";
  } else {
    lpf = lpf_nom;
    F["field_rate_source"] = "standard";
  }
  F["field_rate_hz"] = npr(rate / lpf, 5);
  if (snr_db > 4) {
    double hb = std::max(std::fabs(e_hi - centre), std::fabs(centre - e_lo)) + 1.0e6;
    hb = std::min(hb, fs / 2 - std::fabs(centre));
    const Disc dd = fm_discriminate(ch, centre, hb, fs);
    const Json dm = demod_features(dd.d, dd.fs, rate, t0, std_name, polarity);
    for (auto it = dm.begin(); it != dm.end(); ++it) {
      const std::string& k = it.key();
      if (it.value().is_number_float())
        F[k] = (k == "color_burst_mhz" || k == "audio_subcarrier_mhz") ? pyr(it.value().get<double>(), 3)
                                                                          : npr(it.value().get<double>(), 3);
      else F[k] = it.value();
    }
    if (dm.contains("color") && dm["color"].is_boolean() && !dm["color"].get<bool>() && snr_db < 10) {
      // no burst found, but below ~10 dB a burst is lost in FM-threshold noise: not measurable, not B/W
      F["color"] = nullptr;
      F["color_note"] = fmt("colour burst not measurable at %.1f dB SNR", snr_db);
    }
    if (snr_db > 6 && dm.contains("deviation_demod_mhz") && dm["deviation_demod_mhz"].get<double>() != 0) {
      F["deviation_mhz"] = npr(dm["deviation_demod_mhz"].get<double>(), 2);
      F["deviation_source"] = "discriminator";
    }
  }
  if (!F.contains("deviation_source")) F["deviation_source"] = "spectrum";
  if (!F.contains("audio_subcarrier_mhz")) F["audio_subcarrier_mhz"] = nullptr;
  std::string label = std_name;
  if (F.contains("lines_per_field") && std::fabs(F["lines_per_field"].get<double>() - lpf_nom) > 1.0)
    label = "analog_fm_video";
  if (std::fabs(ppm) > 0.9 * kPpmMax) label = "analog_fm_video";
  // colour standard (research): PAL-M rides 525-line timing, PAL-N 625-line
  if (F.contains("color_standard") && F["color_standard"].is_string()) {
    const std::string cs = F["color_standard"].get<std::string>();
    const std::string fam = (cs == "ntsc" || cs == "pal_m") ? "ntsc" : "pal";
    if (fam != std_name) {
      F["color_note"] = "colour burst of " + cs + " on " + std_name + " line timing: colour standard not used";
      F.erase("color_standard");
    } else if (label == std_name && (cs == "pal_m" || cs == "pal_n" || cs == "secam")) {
      label = cs;
    }
  }
  const double cen_mhz = (fc + centre) / 1e6;
  const auto cands = spec_ok ? match_channel(cen_mhz, 3.0) : std::vector<std::pair<std::string, int>>{};
  F["spectrum_source"] = spec_ok ? sp.source : "unmeasurable";
  if (!spec_ok) {
    F["deviation_mhz"] = nullptr;
    F["bw_mhz"] = nullptr;
    F["rms_dev_mhz"] = nullptr;
  }
  if (!cands.empty()) {
    const int ch_mhz = cands[0].second;
    F["channel"] = cands[0].first;
    Json cc = Json::array();
    for (const auto& c : cands) cc.push_back(c.first);
    F["channel_candidates"] = cc;
    for (const auto& c : cands)
      if (c.first.rfind("LB", 0) == 0) {
        F["channel_note"] = "LB = Lowband 40 MHz-step table (vendor tables, research)";
        break;
      }
    F["nominal_mhz"] = ch_mhz;
    const double cfo = pyr((cen_mhz - ch_mhz) * 1e6, -3);
    F["cfo_hz"] = cfo;
    F["cfo_ppm"] = pyr(cfo / (ch_mhz * 1e6) * 1e6, 1);
  } else {
    F["channel"] = nullptr;
  }
  F["center_mhz"] = pyr(cen_mhz, 3);
  const double conf = std::min(std::max(1 - kThreshold / (2 * best.score), 0.3), 0.99);
  Emission em;
  em.family = "analog_vtx";
  em.label = label;
  em.start_us = 0.0;
  em.end_us = pyr(dur_us, 1);
  em.f_lo_mhz = pyr((fc + e_lo) / 1e6, 3);
  em.f_hi_mhz = pyr((fc + e_hi) / 1e6, 3);
  em.center_mhz = pyr(cen_mhz, 3);
  em.snr_db = pyr(snr_db, 1);
  em.confidence = conf;
  em.decoded = nullptr;
  em.features = std::move(F);
  return em;
}

}  // namespace

std::vector<Emission> detect(const cfvec& xin, double fs, double fc) {
  const size_t n = std::min(xin.size(), kMaxSamples);
  if (n < 1024) return {};
  const cfvec x(xin.begin(), xin.begin() + static_cast<std::ptrdiff_t>(n));
  const double dur_us = static_cast<double>(n) / fs * 1e6;
  const TfPower tf = tf_power(x.data(), n, fs);
  std::vector<std::pair<double, double>> notch;
  for (const Spur& sp : find_spurs(x, fs, fc)) notch.emplace_back(sp.offset_hz, 50e3);
  cdvec xd(n);
  for (size_t i = 0; i < n; ++i) xd[i] = cd(x[i].real(), x[i].imag());
  const SpectralChannelizer ch(xd, fs, notch);
  xd = cdvec();
  const Scanner scanner(tf);
  const dvec& f = tf.f;
  std::vector<uint8_t> mask(f.size(), 0);
  std::vector<Emission> found;
  for (int it = 0; it < kMaxEmitters; ++it) {
    const bool any = std::any_of(mask.begin(), mask.end(), [](uint8_t v) { return v != 0; });
    const ScanBest best = scanner.best(any ? &mask : nullptr);
    if (best.score < kThreshold) break;
    const int nbw = static_cast<int>(py_round_int(10e-6 * best.rate * kNB));
    const auto cols = emitter_columns(best.zb, best.b0, best.k0, nbw);
    for (int k = std::max(cols.first - 1, 0); k < std::min<int>(cols.second + 2, static_cast<int>(mask.size())); ++k)
      mask[k] = 1;
    if (is_part_of(best, f, fc, found)) continue;
    Emission em = characterise(fs, fc, ch, tf, best, cols, dur_us);
    if (capture_locked(em, fs, fc)) continue;
    const int64_t lo_k = static_cast<int64_t>(searchsorted_d(f, (em.f_lo_mhz * 1e6 - fc) - 1e6));
    const int64_t hi_k = static_cast<int64_t>(searchsorted_d(f, (em.f_hi_mhz * 1e6 - fc) + 1e6));
    for (int64_t k = std::max<int64_t>(lo_k, 0); k < std::min<int64_t>(hi_k + 1, static_cast<int64_t>(mask.size())); ++k)
      mask[k] = 1;
    bool dup = false;
    for (const auto& o : found) {
      const Json& a = o.features;
      const Json& b = em.features;
      if (a["standard"] == b["standard"] &&
          std::fabs(a["line_rate_ppm"].get<double>() - b["line_rate_ppm"].get<double>()) < 5 &&
          std::fabs(np::mod((a["sync_phase_us"].get<double>() - b["sync_phase_us"].get<double>()) + 32,
                            1e6 / a["line_rate_hz"].get<double>()) - 32) < 3)
        dup = true;
    }
    if (dup) continue;
    found.push_back(std::move(em));
  }
  return found;
}

std::shared_ptr<Image> video_frame(const cfvec& x, double fs, double fc, const Emission& em, int ncols, int max_lines,
                                   const SpectralChannelizer* chp) {
  const Json& F = em.features;
  auto num = [&](const char* k) -> std::optional<double> {
    if (!F.contains(k) || !F[k].is_number()) return std::nullopt;
    return F[k].get<double>();
  };
  const auto lr = num("line_rate_hz");
  const auto sp = num("sync_phase_us");
  if (!lr || *lr == 0 || !sp) return nullptr;
  std::unique_ptr<SpectralChannelizer> own;
  if (!chp) {
    cdvec xd(x.size());
    for (size_t i = 0; i < x.size(); ++i) xd[i] = cd(x[i].real(), x[i].imag());
    own = std::make_unique<SpectralChannelizer>(xd, fs);
    chp = own.get();
  }
  double hb = std::max(em.f_hi_mhz - em.center_mhz, em.center_mhz - em.f_lo_mhz) * 1e6 + 1e6;
  hb = std::min(hb, fs / 2 - std::fabs(em.center_mhz * 1e6 - fc));
  const Disc dd = fm_discriminate(*chp, em.center_mhz * 1e6 - fc, hb, fs);
  own.reset();
  const double T = 1 / *lr;
  const double t0 = *sp * 1e-6;
  const bool pal = F.value("standard", std::string()) == "pal";
  const double act0 = pal ? 10.5e-6 : 9.4e-6, act1 = pal ? T - 1.65e-6 : T - 1.5e-6;
  const int64_t nl = std::min<int64_t>(
      static_cast<int64_t>((static_cast<double>(dd.d.size()) / dd.fs - t0) / T) - 1, max_lines);
  if (nl < 2) return nullptr;
  const dvec lin = np::linspace<double>(act0, act1, static_cast<size_t>(ncols));
  const size_t N = dd.d.size();
  const bool flip = F.value("polarity", std::string()) != "sync_low";
  dvec img(static_cast<size_t>(nl) * ncols);
  // np.interp(ts * ofs, arange(N), d) on the unit grid: d[j] + (x - j) * (d[j+1] - d[j])
  parallel_for(static_cast<size_t>(nl), [&](size_t r) {
    const double base = t0 + static_cast<double>(r) * T;
    for (int c = 0; c < ncols; ++c) {
      const double xs = (base + lin[c]) * dd.fs;
      double v;
      if (!(xs > 0)) v = dd.d[0];
      else if (xs >= static_cast<double>(N - 1)) v = dd.d[N - 1];
      else {
        const size_t j = static_cast<size_t>(xs);
        const double slope = (dd.d[j + 1] - dd.d[j]) / 1.0;
        v = slope * (xs - static_cast<double>(j)) + dd.d[j];
        if (xs == static_cast<double>(j)) v = dd.d[j];
      }
      img[r * ncols + c] = flip ? -v : v;
    }
  });
  const double lo = np::percentile(img, 1), hi = np::percentile(img, 99);
  const double den = std::max(hi - lo, 1e-12);
  auto out = std::make_shared<Image>();
  out->rows = static_cast<size_t>(nl);
  out->cols = static_cast<size_t>(ncols);
  out->pixels.resize(img.size());
  for (size_t i = 0; i < img.size(); ++i) {
    const double v = std::min(std::max((img[i] - lo) / den, 0.0), 1.0) * 255;
    out->pixels[i] = static_cast<float>(static_cast<uint8_t>(v));
  }
  return out;
}

}  // namespace si::rx::analog_vtx
