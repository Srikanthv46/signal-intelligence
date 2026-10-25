// python: receivers/analog_vtx/video.py
#include "receivers/analog_vtx/video.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "core/ndimage.hpp"
#include "core/npcompat.hpp"
#include "core/parallel.hpp"
#include "core/signal.hpp"

namespace si::rx::analog_vtx {

namespace {

// python CHROMA / LINE_GRID / BURST_MIN / COLOR_LINES
struct Chroma {
  double k, burst_us0, burst_us1, burst_per_sync, white_per_sync, black_per_sync;
};
const Chroma& chroma_standard(const std::string& s) {
  static const Chroma ntsc{227.5, 5.6, 7.6, 0.5, 2.5, 0.1875};
  static const Chroma pal{283.7516, 5.9, 8.0, 0.5, 7.0 / 3, 0.0};
  return s == "pal" ? pal : ntsc;
}
constexpr int kLineGrid = 1024;
constexpr double kBurstMin = 4.0;
constexpr int kBurstSmooth = 3;
constexpr double kBurstShiftUs[] = {0.0, 0.2, 0.4};   // burst gate offsets tried (python BURST_SHIFT_US)
constexpr double kPriorPpm = 2.0;        // the detector's line rate as a prior of this width
constexpr double kSyncCompressed = 1.3;  // picture white / the sync depth's white: compressed sync
constexpr double kNoisy = 0.09;          // sync-tip noise / (white - black): 3x3 median above
constexpr double kNoisyColor = 0.25;     // ... colour killed above

// python lowpass(x, fs, cutoff, ntaps): firwin + 'same' convolution (direct, float64)
fvec lowpass_f(const fvec& d, double fs, double cutoff, size_t ntaps) {
  const dvec h = sig::firwin(ntaps, std::min(cutoff / (fs / 2), 0.99));
  const int64_t N = static_cast<int64_t>(d.size()), H = static_cast<int64_t>(h.size()), off = (H - 1) / 2;
  fvec out(d.size());
  parallel_for((d.size() + 65535) / 65536, [&](size_t blk) {
    const int64_t a = static_cast<int64_t>(blk) * 65536, b = std::min<int64_t>(N, a + 65536);
    for (int64_t i = a; i < b; ++i) {
      double acc = 0;
      for (int64_t j = 0; j < H; ++j) {
        const int64_t k = i + off - j;
        if (k >= 0 && k < N) acc += static_cast<double>(d[k]) * h[j];
      }
      out[i] = static_cast<float>(acc);
    }
  });
  return out;
}

// np.interp(x, arange(len(d)), d) at one point
inline double interp_index(const fvec& d, double xs) {
  const size_t N = d.size();
  if (!(xs > 0)) return d[0];
  if (xs >= static_cast<double>(N - 1)) return d[N - 1];
  const size_t j = static_cast<size_t>(xs);
  return xs == static_cast<double>(j) ? static_cast<double>(d[j])
                                      : (static_cast<double>(d[j + 1]) - d[j]) * (xs - static_cast<double>(j)) + d[j];
}

// python FieldDemodulator._fold(d, fs, rate, nbf)
dvec fold_f(const fvec& d, double fs, double rate, int nbf) {
  dvec S(nbf, 0.0);
  std::vector<int64_t> N(nbf, 0);
  for (size_t i = 0; i < d.size(); ++i) {
    const double v = static_cast<double>(i) / fs * rate;
    const double ph = v - std::floor(v);   // np.mod(v, 1.0), bit-identical for finite v
    const int64_t b = static_cast<int64_t>(ph * nbf);
    S[b] += static_cast<double>(d[i]);
    ++N[b];
  }
  for (int i = 0; i < nbf; ++i) S[i] = S[i] / static_cast<double>(std::max<int64_t>(N[i], 1));
  return S;
}

// python _field_rows(vbi, n_lines)
std::vector<int64_t> field_rows(const std::vector<uint8_t>& vbi, int64_t n_lines, int need, int skip) {
  std::vector<std::pair<int64_t, int64_t>> runs;
  int64_t i = 0;
  while (i < n_lines) {
    if (vbi[i]) {
      int64_t j = i;
      while (j + 1 < n_lines && vbi[j + 1]) ++j;
      if (j - i >= 1) runs.emplace_back(i, j);
      i = j + 1;
    } else {
      ++i;
    }
  }
  auto range = [](int64_t a, int64_t b) {
    std::vector<int64_t> r;
    for (int64_t k = a; k < b; ++k) r.push_back(k);
    return r;
  };
  for (const auto& [s, e] : runs)
    if (e + skip + need <= n_lines) return range(e + skip, e + skip + need);
  for (const auto& [s, e] : runs)
    if (s - need - 2 >= 0) return range(s - need - 2, s - 2);
  return range(0, std::min<int64_t>(need, n_lines));
}

}  // namespace

const VideoStandard& video_standard(const std::string& s) {
  static const VideoStandard ntsc{10.9, 52.6, 8.0, 9.2, 12, 240, 15734.264};
  static const VideoStandard pal{12.0, 52.0, 8.5, 10.0, 17, 288, 15625.0};
  return s == "pal" ? pal : ntsc;
}

std::shared_ptr<Image> VideoFrame::to_image() const {
  auto im = std::make_shared<Image>();
  im->rows = rows;
  im->cols = cols;
  im->pixels.assign(image.begin(), image.end());
  return im;
}

FieldDemodulator::FieldDemodulator(double center_hz, double line_rate_hz, const std::string& standard,
                                   std::optional<double> sync_hz, std::optional<double> bw_hz,
                                   const std::string& polarity, bool color) {
  std::string s = standard;
  std::transform(s.begin(), s.end(), s.begin(), ::tolower);
  std_ = s == "pal" ? "pal" : "ntsc";
  p_ = video_standard(std_);
  center_hz_ = center_hz;
  sync_hz_ = sync_hz ? *sync_hz : center_hz;
  const double bw = (bw_hz && *bw_hz != 0) ? *bw_hz : 12e6;
  half_bw_ = std::max(bw / 2, std::fabs(sync_hz_ - center_hz_)) + 1.5e6;
  rate0_ = rate_ = line_rate_hz;
  polarity_ = polarity == "sync_high" ? -1.0 : 1.0;
  color_ = color;
}

FieldDemodulator FieldDemodulator::from_fields(const Json& f) {
  auto num = [&](const char* k) -> std::optional<double> {
    if (!f.contains(k) || !f[k].is_number()) return std::nullopt;
    const double v = f[k].get<double>();
    return v != 0 ? std::optional<double>(v) : std::nullopt;   // python: `if f.get(k)` (truthy)
  };
  const auto bw = num("bw_mhz");
  const auto st = num("sync_tip_mhz");
  return FieldDemodulator(f.at("center_mhz").get<double>() * 1e6, f.at("line_rate_hz").get<double>(),
                          f.contains("standard") && f["standard"].is_string() ? f["standard"].get<std::string>() : "ntsc",
                          st ? std::optional<double>(*st * 1e6) : std::nullopt,
                          bw ? std::optional<double>(*bw * 1e6) : std::nullopt,
                          f.contains("polarity") && f["polarity"].is_string() ? f["polarity"].get<std::string>()
                                                                               : "sync_low");
}

FieldDemodulator FieldDemodulator::from_emission(const Emission& e) {
  Json f = e.features;
  if (!f.contains("center_mhz") || f["center_mhz"].is_null()) f["center_mhz"] = e.center_mhz;
  return from_fields(f);
}

bool FieldDemodulator::in_band(double fs, double fc) const { return std::fabs(center_hz_ - fc) + 2e6 < fs / 2; }

Channelizer::Slice FieldDemodulator::slice(const cfvec& iq, double fs, double fc, const Channelizer* chp) const {
  const double hb = std::min(half_bw_, fs / 2 - std::fabs(center_hz_ - fc));
  std::unique_ptr<Channelizer> own;
  if (!chp) {
    own = std::make_unique<Channelizer>(iq, fs, fc);
    chp = own.get();
  }
  return chp->slice(center_hz_, hb, std::min(fs, std::max(2.4 * hb, 8e6)));
}

namespace {

// numpy temporary elision: in `y[1:] * np.conj(y[:-1])` the conj() temporary, when it holds
// >= 256 KiB, receives the product in place, i.e. numpy computes conj(prev) * cur; the
// fma-based complex multiply is not symmetric in its last bit, so the order is mirrored.
cfvec lag1_product(const cfvec& y) {
  const size_t n = y.size();
  cfvec z(n > 0 ? n - 1 : 0);
  const bool elide = (n - 1) * sizeof(cf) >= 256 * 1024;
  for (size_t i = 1; i < n; ++i)
    z[i - 1] = elide ? np::cmul(std::conj(y[i - 1]), y[i]) : np::cmul(y[i], std::conj(y[i - 1]));
  return z;
}

// python fm_discriminate (plain) * polarity
fvec discriminate_plain(const cfvec& y, double ofs, float pol) {
  const cfvec z = lag1_product(y);
  const float sc = static_cast<float>(ofs / (2 * kPi));
  fvec d(z.size());
  for (size_t i = 0; i < z.size(); ++i) d[i] = std::atan2(z[i].imag(), z[i].real()) * sc * pol;
  return d;
}

// python fm_discriminate_weighted(y, ofs) * polarity: the lag-1 product low-passed by a 7-tap
// 6 MHz FIR (symmetric, complex64) before the angle, so each phase step counts by its power
fvec discriminate_weighted(const cfvec& y, double ofs, float pol, double cutoff = 6e6, size_t ntaps = 7) {
  const cfvec z = lag1_product(y);
  const dvec hd = sig::firwin(ntaps, std::min(cutoff / (ofs / 2), 0.99));
  std::vector<float> h(hd.begin(), hd.end());
  const int64_t N = static_cast<int64_t>(z.size()), c = static_cast<int64_t>(ntaps / 2);
  const float sc = static_cast<float>(ofs / (2 * kPi));
  fvec d(z.size());
  parallel_for((z.size() + 65535) / 65536, [&](size_t blk) {
    const int64_t a = static_cast<int64_t>(blk) * 65536, b = std::min<int64_t>(N, a + 65536);
    auto at = [&](int64_t q) { return q >= 0 && q < N ? z[q] : cf(0, 0); };
    for (int64_t n = a; n < b; ++n) {
      cf acc = h[c] * z[n];
      for (int64_t i = 0; i < c; ++i) acc += h[i] * (at(n + i - c) + at(n + c - i));
      d[n] = std::atan2(acc.imag(), acc.real()) * sc * pol;
    }
  });
  return d;
}

}  // namespace

fvec FieldDemodulator::discriminate(const cfvec& iq, double fs, double fc, double& ofs, const Channelizer* chp) const {
  const Channelizer::Slice sl = slice(iq, fs, fc, chp);
  ofs = sl.fs;
  return lowpass_f(discriminate_weighted(sl.y, ofs, static_cast<float>(polarity_)), ofs, 4.8e6, 31);
}

fvec FieldDemodulator::discriminate_raw(const cfvec& iq, double fs, double fc, double& ofs,
                                        const Channelizer* chp) const {
  const Channelizer::Slice sl = slice(iq, fs, fc, chp);
  ofs = sl.fs;
  return discriminate_plain(sl.y, ofs, static_cast<float>(polarity_));
}

std::optional<VideoFrame> FieldDemodulator::process(const cfvec& iq, double fs, double fc, const Channelizer* ch) {
  if (!in_band(fs, fc)) return std::nullopt;
  const Channelizer::Slice sl = slice(iq, fs, fc, ch);
  const double ofs = sl.fs;
  const float pol = static_cast<float>(polarity_);
  auto fr = raster(lowpass_f(discriminate_weighted(sl.y, ofs, pol), ofs, 4.8e6, 31), ofs);
  if (!fr) return std::nullopt;
  fr->meta["standard"] = std_;
  fr->meta["center_mhz"] = py_round(center_hz_ / 1e6, 3);
  // colour from the plain discriminator: the weighted one's short filter attenuates the subcarrier
  if (color_) chroma(lowpass_f(discriminate_plain(sl.y, ofs, pol), ofs, 6.0e6, 63), ofs, *fr);
  return fr;
}

std::optional<VideoFrame> FieldDemodulator::raster(const fvec& d, double ofs) {
  const double T_nom = 1 / rate_;
  if (static_cast<double>(d.size()) / ofs < 4 * T_nom) return std::nullopt;
  // coarse copy (~4 MSps) for the sync search
  const size_t k = std::max<size_t>(static_cast<size_t>(std::floor(ofs / 4e6)), 1);
  fvec dc(d.size() / k);
  for (size_t i = 0; i < dc.size(); ++i) dc[i] = np::mean(d.data() + i * k, k);
  const double cfs = ofs / static_cast<double>(k);
  const int nbf = 512;
  // sync phase from the fold at the detector's rate, then the line clock from every line's edge
  double T = 1 / rate0_;
  dvec prof = fold_f(dc, cfs, rate0_, nbf);
  const int w = std::max<int>(static_cast<int>(py_round_int(4.7e-6 / T * nbf)), 1);
  dvec ext(prof);
  ext.insert(ext.end(), prof.begin(), prof.begin() + w);
  const dvec box = np::convolve(ext, dvec(w, 1.0 / w), np::ConvMode::Valid);
  double t0 = static_cast<double>(np::argmin(box.data(), static_cast<size_t>(nbf))) / nbf * T;
  std::optional<EdgeFit> tm;
  for (int it = 0; it < 2; ++it) {
    tm = edge_fit(d, ofs, t0, T);
    if (!tm) break;
    t0 = tm->t0;
    T = tm->T;
  }
  if (tm) prior_blend(t0, T, tm->se, static_cast<double>(d.size()) / ofs);
  rate_ = 1 / T;
  prof = fold_f(dc, cfs, rate_, nbf);
  dvec tip, porch;
  for (int i = 0; i < nbf; ++i) {
    const double ph_us =
        np::mod((static_cast<double>(i) + 0.5) / nbf * T - (t0 - static_cast<double>(k - 1) / 2 / ofs), T) * 1e6;
    if (ph_us > 0.8 && ph_us < 4.0) tip.push_back(prof[i]);
    if (ph_us > p_.porch_us0 && ph_us < p_.porch_us1) porch.push_back(prof[i]);
  }
  const double sync = np::median(tip);
  const double black = np::median(porch);
  const int64_t n_lines = static_cast<int64_t>((static_cast<double>(d.size()) / ofs - t0) / T) - 1;
  if (n_lines < 20) return std::nullopt;
  dvec u(kVideoWidth);
  for (int j = 0; j < kVideoWidth; ++j) u[j] = (p_.active_us0 + j * (p_.active_len_us / kVideoWidth)) * 1e-6;
  dvec ras(static_cast<size_t>(n_lines) * kVideoWidth);
  parallel_for(static_cast<size_t>(n_lines), [&](size_t r) {
    const double base = t0 + static_cast<double>(r) * T;
    for (int c = 0; c < kVideoWidth; ++c) ras[r * kVideoWidth + c] = interp_index(d, (base + u[c]) * ofs);
  });
  // broad-pulse (VBI) lines: most of the line at the sync level
  const double mid = 0.5 * (sync + black);
  std::vector<uint8_t> vbi(n_lines);
  bool any_vbi = false;
  for (int64_t r = 0; r < n_lines; ++r) {
    int cnt = 0;
    for (int c = 0; c < kVideoWidth; ++c) cnt += ras[r * kVideoWidth + c] < mid;
    vbi[r] = static_cast<double>(cnt) / kVideoWidth > 0.3;
    any_vbi |= vbi[r];
  }
  const std::vector<int64_t> rows = field_rows(vbi, n_lines, p_.lines, p_.skip);
  const double noise = sync_noise(d, ofs, t0, T, n_lines);
  Mat<double> field(rows.size(), kVideoWidth);
  for (size_t i = 0; i < rows.size(); ++i)
    std::copy(ras.begin() + rows[i] * kVideoWidth, ras.begin() + (rows[i] + 1) * kVideoWidth, field.row(i));
  bool any_act = false;
  for (int64_t r : rows) any_act |= !vbi[r];
  auto white_of = [&](const Mat<double>& f) {
    dvec pool;
    for (size_t i = 0; i < rows.size(); ++i)
      if (!any_act || !vbi[rows[i]]) pool.insert(pool.end(), f.row(i), f.row(i) + kVideoWidth);
    return np::percentile(pool, 99.5);
  };
  double white = white_of(field);
  if (white - black < 1e3) return std::nullopt;
  double rel = noise / (white - black);
  const bool denoise = rel > kNoisy;
  if (denoise) {                                   // FM clicks: 3x3 median; they do not set white
    field = ndi::median_filter(field, 3, 3, ndi::Mode::Nearest);
    white = white_of(field);
    if (white - black < 1e3) return std::nullopt;
    rel = noise / (white - black);
  }
  timing_ = Timing{t0, T, sync, black, white, rows, rel};
  VideoFrame fr;
  fr.rows = rows.size();
  fr.cols = kVideoWidth;
  fr.image.resize(fr.rows * fr.cols);
  for (size_t i = 0; i < field.data.size(); ++i) {
    const double v = std::min(std::max((field.data[i] - black) / (white - black), 0.0), 1.0);
    fr.image[i] = static_cast<uint8_t>(v * 255);
  }
  const double ppm = (rate_ / p_.fl - 1) * 1e6;
  fr.meta["line_rate_hz"] = np::round(rate_, 3);
  fr.meta["line_rate_ppm"] = np::round(ppm, 1);
  fr.meta["vsync"] = any_vbi;
  fr.meta["lines"] = static_cast<int64_t>(rows.size());
  fr.meta["noise"] = py_round(rel, 3);
  fr.meta["denoise"] = rel > kNoisy;
  if (tm) fr.meta["sync_jitter_us"] = py_round(tm->jitter_us, 3);
  return fr;
}

// python _edge_fit: every line's sync leading edge (box-placed, step-timed) -> robust line fit
std::optional<FieldDemodulator::EdgeFit> FieldDemodulator::edge_fit(const fvec& d, double ofs, double t0,
                                                                    double T) const {
  const int64_t N = static_cast<int64_t>(d.size());
  dvec c(d.size() + 1, 0.0);
  for (size_t i = 0; i < d.size(); ++i) c[i + 1] = c[i] + static_cast<double>(d[i]);
  const int64_t sw = py_round_int(4.7e-6 * ofs), a = py_round_int(0.5e-6 * ofs), b = py_round_int(2.0e-6 * ofs);
  const int64_t w = py_round_int(3.0e-6 * ofs), r = py_round_int(0.6e-6 * ofs);
  const int64_t n_lines = static_cast<int64_t>((static_cast<double>(N) / ofs - t0) / T) - 1;
  dvec kv, te;
  for (int64_t k = 0; k < n_lines; ++k) {
    const int64_t cen = static_cast<int64_t>(np::round((t0 + static_cast<double>(k) * T) * ofs, 0));
    if (!(cen - w - r - a >= 0 && cen + w + r + sw < N)) continue;
    int64_t j1 = cen - w;
    double bmin = std::numeric_limits<double>::infinity();
    for (int64_t i = cen - w; i <= cen + w; ++i) {
      const double v = c[i + sw] - c[i];
      if (v < bmin) { bmin = v; j1 = i; }
    }
    const int64_t M = 2 * r + 1;
    dvec m(M);
    for (int64_t q = 0; q < M; ++q) {
      const int64_t i = j1 - r + q;
      m[q] = (c[i] - c[i - a]) / static_cast<double>(a) - (c[i + b] - c[i]) / static_cast<double>(b);
    }
    const int64_t j = std::min<int64_t>(std::max<int64_t>(static_cast<int64_t>(np::argmax(m)), 1), M - 2);
    const double y0 = m[j - 1], y1 = m[j], y2 = m[j + 1], den = y0 - 2 * y1 + y2;
    const double frac = den < 0 ? 0.5 * (y0 - y2) / den : 0.0;
    kv.push_back(static_cast<double>(k));
    te.push_back((static_cast<double>(j1 - r + j) + frac) / ofs);
  }
  if (kv.size() < 20) return std::nullopt;
  std::vector<uint8_t> keep(kv.size(), 1);
  dvec p, res(kv.size());
  double s = 0;
  for (int it = 0; it < 4; ++it) {   // VBI half-line pulses, click-hit lines: outliers
    dvec xk, yk;
    for (size_t i = 0; i < kv.size(); ++i)
      if (keep[i]) { xk.push_back(kv[i]); yk.push_back(te[i]); }
    p = np::polyfit(xk, yk, 1);
    dvec ar;
    for (size_t i = 0; i < kv.size(); ++i) {
      res[i] = te[i] - np::polyval(p, kv[i]);
      if (keep[i]) ar.push_back(std::fabs(res[i]));
    }
    s = 1.4826 * np::median(ar);
    size_t nk = 0;
    for (size_t i = 0; i < kv.size(); ++i) nk += keep[i] = std::fabs(res[i]) < std::max(3 * s, 0.05e-6);
    if (nk < 20) return std::nullopt;
  }
  double km = 0, rm = 0;
  size_t nk = 0;
  for (size_t i = 0; i < kv.size(); ++i)
    if (keep[i]) { km += kv[i]; rm += res[i]; ++nk; }
  km /= static_cast<double>(nk);
  rm /= static_cast<double>(nk);
  double sk = 0, sr = 0;
  for (size_t i = 0; i < kv.size(); ++i)
    if (keep[i]) { sk += (kv[i] - km) * (kv[i] - km); sr += (res[i] - rm) * (res[i] - rm); }
  const double se = std::sqrt(sr / static_cast<double>(nk)) / std::sqrt(std::max(sk, 1.0));
  return EdgeFit{p[1], p[0], se, s * 1e6};
}

// python _prior_blend: the block's line period weighed against the detector's (+-kPriorPpm),
// unless they disagree by more than 4 combined sigmas; sync phase held at mid-block
void FieldDemodulator::prior_blend(double& t0, double& T, double se, double dur) const {
  const double T0 = 1 / rate0_, sp = kPriorPpm * 1e-6 * T0;
  if (std::fabs(T - T0) > 4 * std::hypot(se, sp)) return;
  const double Tn = se > 0 ? (T / (se * se) + T0 / (sp * sp)) / (1 / (se * se) + 1 / (sp * sp)) : T;
  const double km = (dur / 2 - t0) / T;
  t0 = t0 + km * (T - Tn);
  T = Tn;
}

// python _sync_noise: robust std (Hz) of the sync tip, each line's own level removed
double FieldDemodulator::sync_noise(const fvec& d, double ofs, double t0, double T, int64_t n_lines) const {
  constexpr int K = 24;
  dvec u = np::linspace(1.0e-6, 3.8e-6, K), dev;
  dev.reserve(static_cast<size_t>(n_lines) * K);
  dvec row(K);
  for (int64_t l = 0; l < n_lines; ++l) {
    for (int q = 0; q < K; ++q) row[q] = interp_index(d, (t0 + static_cast<double>(l) * T + u[q]) * ofs);
    const double med = np::median(row);
    for (int q = 0; q < K; ++q) dev.push_back(std::fabs(row[q] - med));
  }
  return 1.4826 * np::median(dev);
}

void FieldDemodulator::chroma(const fvec& d, double ofs, VideoFrame& fr) const {
  if (!timing_ || timing_->rows.empty()) return;
  const Timing& tm = *timing_;
  const Chroma& C = chroma_standard(std_);
  const double T = tm.T;
  double depth = tm.black - tm.sync;                           // luma "black" = blanking (porch)
  // a picture whiter than kSyncCompressed x the sync depth's white: compressed sync (weak VTX near
  // the FM threshold): the depth from the picture's white instead
  const double wspan = C.white_per_sync - C.black_per_sync;
  const bool compressed = depth > 0 && tm.white - tm.black > kSyncCompressed * wspan * depth;
  if (compressed) depth = (tm.white - tm.black) / wspan;
  const double black = tm.black + C.black_per_sync * depth;
  const double span = (C.white_per_sync - C.black_per_sync) * depth;
  std::vector<int64_t> lines;                                  // + the line before (PAL)
  lines.push_back(tm.rows[0] > 0 ? tm.rows[0] - 1 : tm.rows[0]);
  lines.insert(lines.end(), tm.rows.begin(), tm.rows.end());
  const size_t L = lines.size(), G = kLineGrid;
  dvec tau(G);
  for (size_t j = 0; j < G; ++j) tau[j] = static_cast<double>(j) * (T / kLineGrid);
  auto tl = [&](size_t r, size_t j) { return (tm.t0 + static_cast<double>(lines[r]) * T) + tau[j]; };
  if (!(tl(L - 1, G - 1) * ofs < static_cast<double>(d.size()) - 1) || depth <= 0) {
    fr.meta["color"] = false;
    return;
  }
  if (tm.noise > kNoisyColor) {                                // chroma is noise: grey beats colour
    fr.meta["color"] = false;
    fr.meta["color_killer"] = "noise";
    return;
  }
  const double w = (-2 * kPi) * (C.k / T), w_off = (-2 * kPi) * (C.k / T - 1.0e6);
  dvec R(L * G);
  std::vector<std::complex<double>> lo(L * G), Cx(L * G);
  const dvec h = sig::firwin(95, 1.0e6, static_cast<double>(kLineGrid) / T);
  const int64_t H = static_cast<int64_t>(h.size()), hoff = (H - 1) / 2, Gi = static_cast<int64_t>(G);
  dvec tau_us(G);
  for (size_t j = 0; j < G; ++j) tau_us[j] = tau[j] * 1e6;
  constexpr size_t NS = sizeof(kBurstShiftUs) / sizeof(kBurstShiftUs[0]);
  std::vector<size_t> bgs[NS];
  for (size_t q = 0; q < NS; ++q)
    for (size_t j = 0; j < G; ++j)
      if (tau_us[j] >= C.burst_us0 + kBurstShiftUs[q] && tau_us[j] < C.burst_us1 + kBurstShiftUs[q]) bgs[q].push_back(j);
  std::vector<std::complex<double>> zbs[NS], zofs[NS];
  for (size_t q = 0; q < NS; ++q) {
    zbs[q].resize(L);
    zofs[q].resize(L);
  }
  parallel_for(L, [&](size_t r) {
    double* Rr = &R[r * G];
    std::complex<double>* lr = &lo[r * G];
    std::vector<std::complex<double>> m(G);
    for (size_t j = 0; j < G; ++j) {
      const double t = tl(r, j);
      Rr[j] = interp_index(d, t * ofs);
      const double ph = w * t;
      lr[j] = {std::cos(ph), std::sin(ph)};
      m[j] = 2 * Rr[j] * lr[j];
    }
    for (int64_t i = 0; i < Gi; ++i) {                          // 'same' convolution along the line
      std::complex<double> acc = 0;
      for (int64_t k = 0; k < H; ++k) {
        const int64_t q = i + hoff - k;
        if (q >= 0 && q < Gi) acc += m[q] * h[k];
      }
      Cx[r * G + i] = acc;
    }
    // burst phasor from the gate alone; reference: the same gate 1 MHz below the subcarrier
    for (size_t q = 0; q < NS; ++q) {
      const std::vector<size_t>& bg = bgs[q];
      double mean = 0;
      for (size_t j : bg) mean += Rr[j];
      mean /= static_cast<double>(bg.size());
      std::complex<double> z = 0, zo = 0;
      for (size_t j : bg) {
        const double g = Rr[j] - mean;
        const double po = w_off * tl(r, j);
        z += 2 * g * lr[j];
        zo += 2 * g * std::complex<double>(std::cos(po), std::sin(po));
      }
      zbs[q][r] = z / static_cast<double>(bg.size());
      zofs[q][r] = zo / static_cast<double>(bg.size());
    }
  });
  const size_t nr = L - 1;
  // field burst: |mean zb[l] conj zb[l-2]| is the burst power with the noise averaged out; the
  // reference floored at its noise expectation; the gate offset with the strongest burst kept
  const double nb = static_cast<double>(L >= 2 ? L - 2 : 0);
  double field = -1;
  size_t qb = 0;
  for (size_t q = 0; q < NS; ++q) {
    std::complex<double> pb = 0, po = 0;
    double p2 = 0;
    for (size_t r = 2; r < L; ++r) {
      pb += zbs[q][r] * std::conj(zbs[q][r - 2]);
      po += zofs[q][r] * std::conj(zofs[q][r - 2]);
    }
    for (size_t r = 0; r < L; ++r) p2 += std::norm(zofs[q][r]);
    const double pof = std::max(std::abs(po) / nb, p2 / static_cast<double>(L) / std::sqrt(std::max(nb, 1.0)));
    const double f = std::sqrt((std::abs(pb) / nb) / std::max(pof, 1e-18));
    if (f > field) {
      field = f;
      qb = q;
    }
  }
  const std::vector<std::complex<double>>& zb = zbs[qb];
  dvec ratio(L);
  for (size_t r = 0; r < L; ++r) ratio[r] = std::abs(zb[r]) / std::max(std::abs(zofs[qb][r]), 1e-9);
  fr.meta["color"] = false;
  fr.meta["burst_ratio"] = py_round(np::median(dvec(ratio.begin() + 1, ratio.end())), 2);
  fr.meta["burst_field"] = py_round(field, 2);
  fr.meta["burst_gate_us"] = kBurstShiftUs[qb];
  fr.meta["sync_compressed"] = compressed;
  if (!(field >= kBurstMin)) return;
  // phase reference: the burst averaged over same-parity neighbours (+-kBurstSmooth)
  std::vector<std::complex<double>> zs(L);
  for (size_t par = 0; par < 2; ++par) {
    std::vector<std::complex<double>> c(1, 0);
    for (size_t r = par; r < L; r += 2) c.push_back(c.back() + zb[r]);
    const int64_t n = static_cast<int64_t>(c.size()) - 1;
    for (int64_t i = 0; i < n; ++i) {
      const int64_t a = std::max<int64_t>(i - kBurstSmooth, 0), b = std::min<int64_t>(i + kBurstSmooth + 1, n);
      zs[par + 2 * static_cast<size_t>(i)] = (c[b] - c[a]) / static_cast<double>(b - a);
    }
  }
  // burst sin(ph + theta) mixes to e^{j(theta - pi/2)}: NTSC theta = pi, PAL pi +- pi/4
  const bool pal = std_ == "pal";
  std::vector<std::complex<double>> zn(L);
  dvec zs_abs;
  for (size_t r = 0; r < L; ++r) {
    zn[r] = zs[r] / std::max(std::abs(zs[r]), 1e-9);
    if (r > 0) zs_abs.push_back(std::abs(zs[r]));
  }
  std::vector<std::complex<double>> rot(nr);
  dvec sw(nr, 1.0);
  for (size_t i = 0; i < nr; ++i) {
    const std::complex<double> ref = pal ? zn[i + 1] + zn[i] : zn[i + 1];
    if (pal) sw[i] = std::arg(zn[i + 1] * std::conj(ref)) >= 0 ? 1.0 : -1.0;
    rot[i] = std::exp(std::complex<double>(0, -1) * (std::arg(ref) - kPi / 2));
  }
  // ACC: the burst is burst_per_sync x the sync depth; chroma to luma units (white - black = 1)
  const double acc = C.burst_per_sync * depth / np::median(zs_abs) / span;
  const double a_us = p_.active_us0, len_us = p_.active_len_us;
  dvec ju(kVideoWidth), jidx(G);
  for (size_t j = 0; j < G; ++j) jidx[j] = static_cast<double>(j);
  for (int c = 0; c < kVideoWidth; ++c) ju[c] = np::interp(a_us + c * (len_us / kVideoWidth), tau_us, jidx);
  Mat<double> Yg(nr, kVideoWidth), Ug(nr, kVideoWidth), Vg(nr, kVideoWidth);
  parallel_for(nr, [&](size_t i) {
    const size_t r = i + 1;
    const std::complex<double> s = rot[i] * acc;
    const std::complex<double> sp = (i == 0 ? rot[0] : rot[i - 1]) * acc;
    const double swp = i == 0 ? -sw[0] : sw[i - 1];
    dvec Y(G), U(G), V(G);
    for (size_t j = 0; j < G; ++j) {
      const std::complex<double> cn = Cx[r * G + j] * s;
      double u = -cn.imag(), v = cn.real() * sw[i];
      if (pal) {                                                // delay line
        const std::complex<double> cp = Cx[(r - 1) * G + j] * sp;
        u = 0.5 * (u - cp.imag());
        v = 0.5 * (v + cp.real() * swp);
      }
      U[j] = u;
      V[j] = v;
      Y[j] = (R[r * G + j] - (Cx[r * G + j] * std::conj(lo[r * G + j])).real() - black) / span;
    }
    for (int c = 0; c < kVideoWidth; ++c) {
      const size_t j0 = std::min<size_t>(static_cast<size_t>(ju[c]), G - 2);
      const double f = ju[c] - static_cast<double>(j0);
      auto at = [&](const dvec& A) { return A[j0] * (1 - f) + A[j0 + 1] * f; };
      Yg(i, c) = at(Y);
      Ug(i, c) = at(U);
      Vg(i, c) = at(V);
    }
  });
  if (tm.noise > kNoisy) {   // below the FM threshold: clicks in Y, noise in the chroma band
    Yg = ndi::median_filter(Yg, 3, 3, ndi::Mode::Nearest);
    Ug = ndi::median_filter(Ug, 3, 5, ndi::Mode::Nearest);
    Vg = ndi::median_filter(Vg, 3, 5, ndi::Mode::Nearest);
  }
  fr.rgb.assign(nr * kVideoWidth * 3, 0);
  for (size_t i = 0; i < Yg.data.size(); ++i) {
    const double y = Yg.data[i], u = Ug.data[i], v = Vg.data[i];
    const double rgb[3] = {y + 1.140 * v, y - 0.395 * u - 0.581 * v, y + 2.032 * u};
    for (int q = 0; q < 3; ++q) fr.rgb[i * 3 + q] = static_cast<uint8_t>(std::min(std::max(rgb[q], 0.0), 1.0) * 255);
  }
  fr.meta["color"] = true;
}

}  // namespace si::rx::analog_vtx
