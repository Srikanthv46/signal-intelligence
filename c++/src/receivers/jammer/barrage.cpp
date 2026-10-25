// Barrage-noise and noise-FM jammers (research, signal gaps S12; python: receivers/jammer/barrage.py): plateau,
// Gaussianity, no cyclic-prefix / cyclic structure, stationarity.
#include "receivers/jammer/barrage.hpp"

#include <algorithm>
#include <cmath>

#include "core/bands_india.hpp"
#include "core/fft.hpp"
#include "core/npcompat.hpp"

namespace si::rx::jammer::barrage {

namespace {

constexpr double kMinBw = 1e6;
constexpr size_t kNSub = 64;
constexpr double kFlatDb = 3.0, kFlatFrac = 0.9, kEdgeDb = 6.0, kLineDb = 6.0;
constexpr size_t kMaxN = size_t(1) << 17;
constexpr double kCpMin = 0.025, kCpLagBw = 8.0;
constexpr double kCycMargin = 9.0;
constexpr size_t kCycBlock = 512;
constexpr double kStatDb = 1.0, kMinJamBw = 10e6, kMatchHz = 2e6;
constexpr double kPi = 3.14159265358979323846;

double db(double x) { return 10 * std::log10(x + 1e-30); }

// python Channelizer._bins(f_off, half_bw)
std::vector<size_t> bins(const Channelizer& ch, double f_off, double half) {
  const size_t L = ch.L();
  const double df = ch.fs() / static_cast<double>(L);
  int64_t lo = static_cast<int64_t>(std::floor((f_off - half) / df)) - 2;
  int64_t hi = static_cast<int64_t>(std::ceil((f_off + half) / df)) + 2;
  lo = std::max<int64_t>(lo, -static_cast<int64_t>(L / 2) - 1);
  hi = std::min<int64_t>(hi, static_cast<int64_t>((L - 1) / 2) + 1);
  std::vector<size_t> idx;
  for (int64_t k = lo; k <= hi; ++k) idx.push_back(static_cast<size_t>(((k % static_cast<int64_t>(L)) + L) % L));
  std::sort(idx.begin(), idx.end());
  idx.erase(std::unique(idx.begin(), idx.end()), idx.end());
  std::vector<size_t> out;
  for (size_t i : idx)
    if (std::fabs(ch.f()[i] - f_off) <= half) out.push_back(i);
  return out;
}

dvec power_at(const Channelizer& ch, const std::vector<size_t>& k) {
  dvec p(k.size());
  for (size_t i = 0; i < k.size(); ++i) {
    const double a = static_cast<double>(np::cabs(ch.X()[k[i]]));
    p[i] = a * a;
  }
  return p;
}

// python plateau(ch, f0, f1): adds flatness, edge_db, line_db, plateau_ok to out
void plateau(const Channelizer& ch, double f0, double f1, Json& out) {
  const double bw = f1 - f0;
  const dvec P = power_at(ch, bins(ch, (f0 + f1) / 2 - ch.fc(), 0.45 * bw));
  dvec sub;   // np.array_split: the first len % 64 pieces one longer
  const size_t q = P.size() / kNSub, r = P.size() % kNSub;
  size_t pos = 0;
  for (size_t i = 0; i < kNSub; ++i) {
    const size_t len = q + (i < r ? 1 : 0);
    sub.push_back(np::mean(P.data() + pos, len));
    pos += len;
  }
  const double med = np::median(sub);
  size_t nflat = 0;
  for (double v : sub) nflat += std::fabs(10 * std::log10(v / med + 1e-30)) <= kFlatDb;
  const double flat = static_cast<double>(nflat) / static_cast<double>(kNSub);
  Json edges = Json::array();
  double emin = 1e300;
  for (const auto& [a, b] : {std::make_pair(f0 - 0.15 * bw, f0 - 0.05 * bw), std::make_pair(f1 + 0.05 * bw, f1 + 0.15 * bw)}) {
    if (a < ch.fc() - 0.45 * ch.fs() || b > ch.fc() + 0.45 * ch.fs()) continue;
    const double e = py_round(db(med / np::mean(power_at(ch, bins(ch, (a + b) / 2 - ch.fc(), (b - a) / 2)))), 1);
    edges.push_back(e);
    emin = std::min(emin, e);
  }
  const size_t m = P.size() / 16;
  dvec sm;
  if (m) {
    for (size_t i = 0; i < m; ++i) sm.push_back(np::mean(P.data() + 16 * i, 16));
  } else {
    sm = P;
  }
  const double line = db(*std::max_element(sm.begin(), sm.end()) / np::median(sm));
  out["flatness"] = py_round(flat, 3);
  out["edge_db"] = edges;
  out["line_db"] = py_round(line, 1);
  out["plateau_ok"] = flat >= kFlatFrac && !edges.empty() && emin >= kEdgeDb && line < kLineDb;
}

// python _cyclic_db(v, fs, f_lo, f_hi) on a complex input (a real one has zero imaginary parts)
double cyclic_db(const cdvec& v, double fs, double f_lo, double f_hi) {
  const size_t n = v.size();
  const cd mu = np::mean(v);
  cdvec z(n);
  for (size_t i = 0; i < n; ++i) z[i] = v[i] - mu;
  const cdvec F = fft::fft(z);
  const dvec fr = fft::fftfreq(n, 1 / fs);
  const size_t h = n / 2;   // fftshift: index i of the shifted array is (i + n - n/2) % n ... numpy: roll by n // 2
  dvec s;
  for (size_t i = 0; i < n; ++i) {
    const size_t j = (i + n - h) % n;
    const double f = std::fabs(fr[j]);
    if (f >= f_lo && f <= f_hi) {
      const double a = np::cabs(F[j]);
      s.push_back(a * a);
    }
  }
  const size_t m = s.size() / kCycBlock;
  if (m < 1) return -99.0;
  double best = 0.0;
  dvec blk(kCycBlock);
  for (size_t b = 0; b < m; ++b) {
    std::copy(s.begin() + static_cast<ptrdiff_t>(b * kCycBlock), s.begin() + static_cast<ptrdiff_t>((b + 1) * kCycBlock),
              blk.begin());
    const double mx = *std::max_element(blk.begin(), blk.end());
    const double ratio = mx / (np::median(blk) / std::log(2.0));
    best = b == 0 ? ratio : std::max(best, ratio);
  }
  return db(best / (std::log(static_cast<double>(m * kCycBlock)) + kCycMargin));
}

cdvec as_complex(const dvec& x) {
  cdvec z(x.size());
  for (size_t i = 0; i < x.size(); ++i) z[i] = cd(x[i], 0.0);
  return z;
}

}  // namespace

std::vector<std::string> covers_band(double f0, double f1) {
  std::vector<std::string> out;
  for (const auto& b : bands::bands())
    if (std::fabs(f0 - b.lo) <= kMatchHz && std::fabs(f1 - b.hi) <= kMatchHz) out.push_back(b.key);
  return out;
}

std::optional<Json> analyse(const Channelizer& ch, double f0, double f1) {
  const double bw = f1 - f0;
  if (bw < kMinBw) return std::nullopt;
  Json out = Json::object();
  out["plateau_bw_mhz"] = py_round(bw / 1e6, 2);
  plateau(ch, f0, f1, out);
  const auto sl = ch.slice((f0 + f1) / 2, 0.45 * bw, std::min(1.25 * bw, ch.fs()));
  const double fs = sl.fs;
  if (sl.y.size() < 4096) return std::nullopt;
  dvec pall(sl.y.size());
  for (size_t i = 0; i < sl.y.size(); ++i) {
    const double a = static_cast<double>(np::cabs(sl.y[i]));
    pall[i] = a * a;
  }
  const size_t n = std::min(sl.y.size(), kMaxN);   // the FFT tests on the first 2^17 samples (latency)
  cdvec y(n);
  for (size_t i = 0; i < n; ++i) y[i] = cd(sl.y[i].real(), sl.y[i].imag());
  dvec p(n), p2(n), r2(n), r4(n), i2(n), i4(n);
  for (size_t i = 0; i < n; ++i) {
    const double a = np::cabs(y[i]);
    p[i] = a * a;
    p2[i] = p[i] * p[i];
    r2[i] = y[i].real() * y[i].real();
    r4[i] = std::pow(y[i].real(), 4.0);
    i2[i] = y[i].imag() * y[i].imag();
    i4[i] = std::pow(y[i].imag(), 4.0);
  }
  const double mp = np::mean(p);
  const double k = np::mean(p2) / (mp * mp);
  const double mr = np::mean(r2), mi = np::mean(i2);
  const double ki = np::mean(r4) / (mr * mr), kq = np::mean(i4) / (mi * mi);
  out["envelope_k"] = py_round(k, 3);
  out["kurtosis_i"] = py_round(ki, 3);
  out["kurtosis_q"] = py_round(kq, 3);
  // autocorrelation of the Hann-tapered spectrum (no sinc sidelobes)
  const cdvec Y = fft::fft(y, 2 * n);
  const dvec fy = fft::fftfreq(2 * n, 1 / fs);
  cdvec S(2 * n);
  for (size_t i = 0; i < 2 * n; ++i) {
    const double a = np::cabs(Y[i]);
    const double c = std::cos(kPi * fy[i] / (0.9 * bw));
    const double w = std::fabs(fy[i]) < 0.45 * bw ? c * c : 0.0;
    S[i] = cd(a * a * w, 0.0);
  }
  const cdvec R = fft::ifft(S);
  const double r0 = np::cabs(R[0]);
  const int64_t l0 = static_cast<int64_t>(std::ceil(std::max(1e-6, kCpLagBw / bw) * fs));
  const int64_t l1 = std::min<int64_t>(static_cast<int64_t>(1000e-6 * fs), static_cast<int64_t>(n / 4));
  if (l1 > l0) {
    int64_t j = l0;
    double best = -1;
    for (int64_t i = l0; i < l1; ++i) {
      const double v = np::cabs(R[static_cast<size_t>(i)]) / r0;
      if (v > best) {
        best = v;
        j = i;
      }
    }
    out["cp_corr"] = py_round(best, 4);
    out["cp_lag_us"] = py_round(static_cast<double>(j) / fs * 1e6, 2);
  } else {
    out["cp_corr"] = 0.0;
    out["cp_lag_us"] = nullptr;
  }
  const double cp_thr = std::max(kCpMin, 5 / std::sqrt(static_cast<double>(n)));
  out["cyclic_env_db"] = py_round(cyclic_db(as_complex(p), fs, 1e3, bw), 1);
  cdvec ysq(n);
  for (size_t i = 0; i < n; ++i) {
    const double a = y[i].real(), b = y[i].imag();
    ysq[i] = cd(a * a - b * b, a * b + b * a);
  }
  out["cyclic_conj_db"] = py_round(cyclic_db(ysq, fs, 0.0, fs / 2), 1);
  const size_t ms = static_cast<size_t>(1e-3 * fs);
  const size_t nb = ms ? pall.size() / ms : 0;
  if (nb >= 4) {   // stationarity over the whole block
    dvec pw(nb);
    for (size_t b = 0; b < nb; ++b) pw[b] = np::mean(pall.data() + b * ms, ms);
    const double m = np::mean(pw);
    double dev = 0;
    for (double v : pw) dev = std::max(dev, std::fabs(10 * std::log10(v / m)));
    out["stationarity_db"] = py_round(dev, 2);
  } else {
    out["stationarity_db"] = 0.0;
  }
  const double cpv = out["cp_corr"].get<double>(), env = out["cyclic_env_db"].get<double>(),
               conj = out["cyclic_conj_db"].get<double>(), stat = out["stationarity_db"].get<double>();
  const bool structure_free = cpv < cp_thr && env < 0 && conj < 0 && stat <= kStatDb;
  const auto cb = covers_band(f0, f1);
  out["covers_band"] = cb;
  Json cls = nullptr;
  if (1.85 <= k && k <= 2.15 && 2.7 <= ki && ki <= 3.3 && 2.7 <= kq && kq <= 3.3 && structure_free &&
      out["plateau_ok"].get<bool>() && (bw >= kMinJamBw || !cb.empty())) {
    cls = "barrage";
  } else if (k <= 1.2) {
    // numpy temporary elision (see voice_fm): a conj() temporary >= 256 KiB receives the product in place
    const bool elide = (n - 1) * sizeof(cd) >= 256 * 1024;
    dvec d(n - 1);
    const double sc = fs / (2 * kPi);
    for (size_t i = 1; i < n; ++i) {
      const cd z = elide ? np::cmul(std::conj(y[i - 1]), y[i]) : np::cmul(y[i], std::conj(y[i - 1]));
      d[i - 1] = std::atan2(z.imag(), z.real()) * sc;
    }
    const double sd = np::stddev(d);
    const double dm = np::mean(d);
    dvec d4(d.size());
    for (size_t i = 0; i < d.size(); ++i) d4[i] = std::pow(d[i] - dm, 4.0);
    const double var = np::var(d);
    const double kd = np::mean(d4) / (var * var);
    out["if_std_mhz"] = py_round(sd / 1e6, 3);
    out["if_kurtosis"] = py_round(kd, 3);
    out["cyclic_if_db"] = py_round(cyclic_db(as_complex(d), fs, 1e3, fs / 2), 1);
    if (sd >= 0.15 * bw && 2.2 <= kd && kd <= 4.5 && out["cyclic_if_db"].get<double>() < 0 && cpv < cp_thr &&
        stat <= kStatDb)
      cls = "noise_fm";
  }
  out["class"] = cls;
  return out;
}

}  // namespace si::rx::jammer::barrage
