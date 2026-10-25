// python receivers/gnss/codes.py + acquire.py (see there for the method and the sources).
#include "receivers/gnss/acquire.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

#include "core/fft.hpp"
#include "core/npcompat.hpp"
#include "core/parallel.hpp"

namespace si::rx::gnss {

namespace {
constexpr double kPpdMinFs = 30e6, kPpdWindowS = 400e-6, kPpdRidgeDb = 6.0, kPpdCoherence = 0.6;
constexpr int kPpdNfft = 32, kPpdHop = 16;
}  // namespace

namespace {

constexpr int kChips = 1023;
constexpr double kRate = 1.023e6;
constexpr double kFsDec = 2.048e6;
constexpr double kHalfBw = 1.024e6;
constexpr int kTCohMs = 4;
constexpr int kMGate = 4;
constexpr double kBinHz = 125.0;
constexpr double kPfaPrn = 1e-3;
constexpr double kSecondDb = 15.0;
constexpr double kCn0High = 45.0, kCn0Implausible = 55.0;
constexpr double kLImplDb = 1.8;
constexpr double kXcorrDb = 18.0, kXcorrNearDb = 15.0;

const int kTaps[33][2] = {{0, 0}, {2, 6}, {3, 7}, {4, 8}, {5, 9}, {1, 9}, {2, 10}, {1, 8}, {2, 9}, {3, 10},
                          {2, 3}, {3, 4}, {5, 6}, {6, 7}, {7, 8}, {8, 9}, {9, 10}, {1, 4}, {2, 5}, {3, 6},
                          {4, 7}, {5, 8}, {6, 9}, {1, 3}, {4, 6}, {5, 7}, {6, 8}, {7, 9}, {8, 10}, {1, 6},
                          {2, 7}, {3, 8}, {4, 9}};

cdvec replica_fft(int prn, size_t n, double fs_eff) {
  const auto code = ca_code(prn);
  cdvec r(n);
  for (size_t i = 0; i < n; ++i) {
    const int64_t k = static_cast<int64_t>(std::floor(static_cast<double>(i) * (kRate / fs_eff))) % kChips;
    r[i] = cd(1.0 - 2.0 * static_cast<double>(code[k]), 0.0);
  }
  cdvec F = fft::fft(r);
  for (cd& v : F) v = std::conj(v);
  return F;
}

// python _search: S[prn][doppler][code sample]
std::vector<dmat> search(const cfvec& y, size_t step, double fs, double fs_eff, const std::vector<int>& prns,
                         const dvec& dop, int m) {
  const size_t n = static_cast<size_t>(py_round_int(fs * 1e-3));
  const size_t k = kTCohMs;
  const size_t len = static_cast<size_t>(m) * k * n;
  cdvec x(len);
  for (size_t i = 0; i < len; ++i) x[i] = cd(y[i * step].real(), y[i * step].imag());
  std::vector<cdvec> R;
  for (int p : prns) R.push_back(replica_fft(p, n, fs_eff));
  std::vector<dmat> out(prns.size(), dmat(dop.size(), n));
  parallel_for(dop.size(), [&](size_t j) {
    const double fd = dop[j];
    const double w = -2 * kPi * fd;
    std::vector<cdvec> F(m);
    for (int s = 0; s < m; ++s) {
      cdvec seg(n, cd(0, 0));
      for (size_t q = 0; q < k; ++q)
        for (size_t i = 0; i < n; ++i) {
          const size_t a = (static_cast<size_t>(s) * k + q) * n + i;
          const double th = w * (static_cast<double>(a) / fs);
          seg[i] += np::cmul(x[a], cd(std::cos(th), std::sin(th)));
        }
      F[s] = fft::fft(seg);
    }
    cdvec prod(n);
    dvec acc(n);
    for (size_t p = 0; p < prns.size(); ++p) {
      std::fill(acc.begin(), acc.end(), 0.0);
      for (int s = 0; s < m; ++s) {
        for (size_t i = 0; i < n; ++i) prod[i] = np::cmul(F[s][i], R[p][i]);
        const cdvec c = fft::ifft(prod);
        for (size_t i = 0; i < n; ++i) {
          const double a = np::cabs(c[i]);
          acc[i] += a * a;
        }
      }
      for (size_t i = 0; i < n; ++i) out[p](j, i) = acc[i] / static_cast<double>(m);
    }
  });
  return out;
}

struct Peak {
  size_t j, i;
  double pk, noise, second;
};

Peak peak(const dmat& g, double chip_samples, int64_t bin_excl) {
  size_t bj = 0, bi = 0;
  for (size_t j = 0; j < g.rows; ++j)
    for (size_t i = 0; i < g.cols; ++i)
      if (g(j, i) > g(bj, bi)) bj = j, bi = i;
  const int64_t n = static_cast<int64_t>(g.cols);
  dvec rest;
  double second = -std::numeric_limits<double>::infinity();
  for (size_t j = 0; j < g.rows; ++j)
    for (size_t i = 0; i < g.cols; ++i) {
      const int64_t d = ((static_cast<int64_t>(i) - static_cast<int64_t>(bi) + n / 2) % n + n) % n - n / 2;
      const bool region = std::llabs(static_cast<long long>(static_cast<int64_t>(j) - static_cast<int64_t>(bj))) <= bin_excl &&
                          static_cast<double>(std::llabs(d)) <= 1.5 * chip_samples;
      if (region) continue;
      rest.push_back(g(j, i));
      second = std::max(second, g(j, i));
    }
  return {bj, bi, g(bj, bi), np::mean(rest), second};
}


double band_power(const Channelizer& ch, double f0, double half, bool& ok) {
  dvec v;
  for (size_t i = 0; i < ch.L(); ++i)
    if (std::fabs(ch.f()[i] + ch.fc() - f0) <= half && std::fabs(ch.f()[i]) < 0.45 * ch.fs()) {
      const double a = static_cast<double>(np::cabs(ch.X()[i]));
      v.push_back(a * a);
    }
  ok = v.size() >= 16;
  return ok ? np::mean(v) : 0.0;
}

Json interference(const Channelizer& ch, double f_l1, const dvec& exclude_hz) {
  bool have_in;
  const double inband = band_power(ch, f_l1, 1.023e6, have_in);
  dvec use;
  for (size_t i = 0; i < ch.L(); ++i) {
    const double f = ch.f()[i] + ch.fc();
    if (std::fabs(ch.f()[i]) < 0.45 * ch.fs() && (f < 1559e6 || f > 1610e6)) {
      const double a = static_cast<double>(np::cabs(ch.X()[i]));
      use.push_back(a * a);
    }
  }
  const bool have_ref = use.size() >= 4096;
  const double ref = have_ref ? np::median(use) / std::log(2.0) : 0.0;
  Json out = Json::object();
  out["jn_db"] = (have_in && inband != 0 && have_ref && ref != 0) ? Json(py_round(10 * std::log10(inband / ref), 2))
                                                                   : Json(nullptr);
  std::vector<size_t> sel;
  for (size_t i = 0; i < ch.L(); ++i)
    if (std::fabs(ch.f()[i] + ch.fc() - f_l1) <= 1.023e6 && std::fabs(ch.f()[i]) < 0.45 * ch.fs()) sel.push_back(i);
  if (!sel.empty()) {
    dvec p(sel.size()), neg(sel.size());
    for (size_t q = 0; q < sel.size(); ++q) {
      const double a = static_cast<double>(np::cabs(ch.X()[sel[q]]));
      p[q] = a * a;
      neg[q] = -p[q];
    }
    const double med = np::median(p) + 1e-30;
    std::vector<std::pair<double, double>> lines;
    const auto order = np::argsort_stable(neg);
    for (size_t r = 0; r < order.size() && r < 8; ++r) {
      const size_t k = sel[order[r]];
      const double a = static_cast<double>(np::cabs(ch.X()[k]));
      const double pk = a * a;
      if (10 * std::log10(pk / med) < 25) break;
      const double fk = ch.f()[k] + ch.fc();
      bool skip = false;
      for (double e : exclude_hz) skip = skip || std::fabs(fk - e) < 2e3;
      for (const auto& [g, _] : lines) skip = skip || std::fabs(fk - g) < 5e3;
      if (skip) continue;
      lines.emplace_back(fk, py_round(10 * std::log10(pk / med), 1));
    }
    Json L = Json::array();
    for (const auto& [a, b] : lines) L.push_back({{"freq_mhz", py_round(a / 1e6, 6)}, {"db_over_band", b}});
    out["lines"] = L;
  }
  if (auto sw = swept_ppd(ch, f_l1)) out["swept"] = *sw;
  return out;
}

Json indicators(const std::vector<Json>& sats, const Json& itf) {
  dvec strong;
  double mx = -std::numeric_limits<double>::infinity();
  for (const Json& s : sats) {
    const double c = s["cn0_dbhz"].get<double>();
    mx = std::max(mx, c);
    if (c > kCn0High) strong.push_back(c);
  }
  const bool has_jn = !itf["jn_db"].is_null();
  const double jn = has_jn ? itf["jn_db"].get<double>() : 0.0;
  const size_t nlines = itf.contains("lines") ? itf["lines"].size() : 0;
  int64_t multi = 0;
  for (const Json& s : sats) multi += s["multipeak"].get<bool>();
  Json sp = Json::object();
  sp["s1_cn0_implausible"] = !sats.empty() && mx > kCn0Implausible;
  sp["s2_uniform_cn0"] = strong.size() >= 4 && np::stddev(strong) < 1.5;
  sp["s3_power_rise_with_tracking"] = has_jn && jn >= 1.5 && strong.size() >= 4;
  sp["s4_multipeak"] = multi >= 2;
  int64_t n = 0, others = 0;
  for (auto it = sp.begin(); it != sp.end(); ++it) {
    n += it.value().get<bool>();
    if (it.key() != "s3_power_rise_with_tracking") others += it.value().get<bool>();
  }
  const bool spoofing = sp["s3_power_rise_with_tracking"].get<bool>() || (sats.size() >= 4 && others >= 2);
  const bool jam_floor = has_jn && jn >= 6.0;
  const bool swept = itf.contains("swept");
  const bool jamming = swept || ((jam_floor || nlines > 0) && !spoofing && strong.size() < 4);
  sp["n"] = n;
  sp["verdict"] = spoofing;
  Json jm = Json::object();
  jm["floor_rise"] = jam_floor;
  jm["narrowband_lines"] = static_cast<int64_t>(nlines);
  jm["swept"] = swept;
  jm["verdict"] = jamming;
  return Json{{"spoofing", sp}, {"jamming", jm}};
}

}  // namespace

// python swept_ppd(ch, f_l1): a fast swept carrier (GNSS personal privacy device) around L1
std::optional<Json> swept_ppd(const Channelizer& ch, double f_l1) {
  if (ch.fs() < kPpdMinFs) return std::nullopt;
  const double half = std::min(25e6, 0.45 * ch.fs() - std::fabs(f_l1 - ch.fc()));
  if (half < 2e6) return std::nullopt;
  const auto sl = ch.slice(f_l1, half, ch.fs());
  const double fs = sl.fs;
  const size_t n = std::min(sl.y.size(), static_cast<size_t>(kPpdWindowS * fs));
  if (n < kPpdNfft) return std::nullopt;
  const size_t nf = (n - kPpdNfft) / kPpdHop + 1;
  if (nf < 64) return std::nullopt;
  const dvec win = np::hanning(kPpdNfft);
  dvec ridge(nf), ratio_db(nf);
  for (size_t r = 0; r < nf; ++r) {
    cdvec fr(kPpdNfft);
    for (int i = 0; i < kPpdNfft; ++i) {
      const cf v = sl.y[r * kPpdHop + static_cast<size_t>(i)];
      fr[i] = cd(static_cast<double>(v.real()) * win[i], static_cast<double>(v.imag()) * win[i]);
    }
    const cdvec F = fft::fft(fr);
    dvec P(kPpdNfft);
    for (int i = 0; i < kPpdNfft; ++i) {
      const double a = np::cabs(F[i]);
      P[i] = a * a;
    }
    int k = 0;
    for (int i = 1; i < kPpdNfft; ++i)
      if (P[i] > P[k]) k = i;
    const double pk = P[k];
    ratio_db[r] = 10 * std::log10(pk / (np::median(P) + 1e-30) + 1e-30);
    const double a = std::log(P[(k - 1 + kPpdNfft) % kPpdNfft] + 1e-30), b = std::log(pk + 1e-30),
                 c = std::log(P[(k + 1) % kPpdNfft] + 1e-30);
    const double den = a - 2 * b + c;
    const double frac = den < 0 ? 0.5 * (a - c) / den : 0.0;
    const double kk = k < kPpdNfft / 2 ? k : k - kPpdNfft;
    ridge[r] = (kk + frac) * (fs / kPpdNfft);
  }
  if (np::median(ratio_db) < kPpdRidgeDb) return std::nullopt;
  const double frate = fs / kPpdHop;
  const double mu = np::mean(ridge);
  dvec x(nf), x2(nf);
  for (size_t i = 0; i < nf; ++i) {
    x[i] = ridge[i] - mu;
    x2[i] = x[i] * x[i];
  }
  const double var = np::mean(x2);
  const size_t m = nf;
  const int64_t lo = std::max<int64_t>(2, py_round_int(0.5e-6 * frate));
  const int64_t hi = std::min<int64_t>(static_cast<int64_t>(50e-6 * frate), static_cast<int64_t>(m / 4));
  if (var <= 0 || hi <= lo + 2) return std::nullopt;
  const cdvec X = fft::rfft(x, 2 * m);
  cdvec X2(X.size());
  for (size_t i = 0; i < X.size(); ++i) {
    const double a = np::cabs(X[i]);
    X2[i] = cd(a * a, 0.0);
  }
  const dvec rr = fft::irfft(X2.data(), X2.size(), 2 * m);
  double top = -std::numeric_limits<double>::infinity();
  for (int64_t q = lo + 1; q < hi - 1; ++q) top = std::max(top, rr[q]);
  int64_t j = -1;
  for (int64_t q = lo + 1; q < hi - 1; ++q)
    if (rr[q] > rr[q - 1] && rr[q] >= rr[q + 1] && rr[q] >= 0.85 * top) {
      j = q;
      break;
    }
  if (j < 0 || rr[j] / rr[0] < 0.5) return std::nullopt;
  double T = static_cast<double>(j);
  {
    const double ra = rr[j - 1], rb = rr[j], rc = rr[j + 1];
    const double d2 = ra - 2 * rb + rc;
    if (d2 < 0) T += 0.5 * (ra - rc) / d2;
  }
  const int64_t nb = std::max<int64_t>(4, py_round_int(T));
  std::vector<int64_t> ph(m);
  dvec cnt(nb, 0.0), sum(nb, 0.0);
  for (size_t i = 0; i < m; ++i) {
    const double v = std::fmod(static_cast<double>(i), T);
    ph[i] = static_cast<int64_t>(std::floor(v / T * static_cast<double>(nb))) % nb;
    cnt[ph[i]] += 1.0;
    sum[ph[i]] += x[i];
  }
  dvec fold(nb);
  for (int64_t b = 0; b < nb; ++b) fold[b] = sum[b] / std::max(cnt[b], 1.0);
  dvec f2(m);
  for (size_t i = 0; i < m; ++i) f2[i] = fold[ph[i]] * fold[ph[i]];
  const double coh = np::mean(f2) / var;
  const double fmax = *std::max_element(fold.begin(), fold.end()), fmin = *std::min_element(fold.begin(), fold.end());
  const double span = fmax - fmin;
  if (coh < kPpdCoherence || span < 1e6) return std::nullopt;
  const double centre = mu + 0.5 * (fmax + fmin);
  const int64_t lag = std::max<int64_t>(1, nb / 10);
  int64_t ups = 0;
  for (int64_t b = 0; b < nb; ++b) ups += (fold[(b + lag) % nb] - fold[b]) > 0;
  const double up = static_cast<double>(ups) / static_cast<double>(nb);
  Json o = Json::object();
  o["period_us"] = py_round(T / frate * 1e6, 3);
  o["span_mhz"] = py_round(span / 1e6, 2);
  o["centre_offset_mhz"] = py_round(centre / 1e6, 2);
  o["covers_l1"] = centre - span / 2 <= 0 && 0 <= centre + span / 2;
  o["shape"] = up >= 0.85 ? "sawtooth up" : up <= 0.15 ? "sawtooth down" : "triangle";
  o["coherence"] = py_round(coh, 3);
  return o;
}

std::vector<uint8_t> ca_code(int prn) {
  const int s1 = kTaps[prn][0], s2 = kTaps[prn][1];
  int g1[10], g2[10];
  for (int i = 0; i < 10; ++i) g1[i] = g2[i] = 1;
  std::vector<uint8_t> out(kChips);
  for (int n = 0; n < kChips; ++n) {
    out[n] = static_cast<uint8_t>(g1[9] ^ g2[s1 - 1] ^ g2[s2 - 1]);
    const int f1 = g1[2] ^ g1[9];
    const int f2 = g2[1] ^ g2[2] ^ g2[5] ^ g2[7] ^ g2[8] ^ g2[9];
    for (int i = 9; i > 0; --i) g1[i] = g1[i - 1], g2[i] = g2[i - 1];
    g1[0] = f1;
    g2[0] = f2;
  }
  return out;
}

double chi2_sf_even(double x, int dof) {
  const double h = x / 2;
  double term = 1.0, s = 1.0;
  for (int k = 1; k < dof / 2; ++k) {
    term *= h / k;
    s += term;
  }
  return std::exp(-h) * s;
}

double threshold(int m, double p_cell) {
  double lo = 0.0, hi = 1000.0;
  for (int it = 0; it < 100; ++it) {
    const double mid = 0.5 * (lo + hi);
    if (chi2_sf_even(2 * m * mid, 2 * m) > p_cell) lo = mid;
    else hi = mid;
  }
  return hi;
}

std::optional<Json> analyse(const Channelizer& ch, std::optional<double> ppm, const dvec& exclude_hz) {
  const double eps = (ppm ? *ppm : 0.0) * 1e-6;
  const double f_meas = kFL1 / (1 + eps);
  const double usable = 0.45 * ch.fs();
  if (ch.fs() < 2.1e6 || std::fabs(f_meas - ch.fc()) + 1.1e6 > usable ||
      static_cast<double>(ch.L()) / ch.fs() < 0.016)
    return std::nullopt;
  const auto sl = ch.slice(f_meas, kHalfBw, kFsDec);
  const double fs = sl.fs, fs_eff = fs * (1 + eps);
  std::vector<int> prns;
  for (int p = 1; p <= 32; ++p) prns.push_back(p);
  const double span = 5e3 + (ppm ? 3.2e3 : 20e3);
  dvec dop;
  for (int64_t i = 0;; ++i) {   // np.arange(-span, span + BIN/2, BIN)
    const double v = -span + static_cast<double>(i) * kBinHz;
    if (!(v < span + kBinHz / 2)) break;
    dop.push_back(v);
  }
  // gate at 1 sample / chip (every second sample)
  const double fg = fs / 2;
  if ((sl.y.size() + 1) / 2 < static_cast<size_t>(kMGate * kTCohMs * py_round_int(fg * 1e-3))) return std::nullopt;
  const auto G = search(sl.y, 2, fg, fs_eff / 2, prns, dop, kMGate);
  const size_t n_g = G[0].cols;
  const double gamma = threshold(kMGate, kPfaPrn / static_cast<double>(n_g * dop.size()));
  struct Hit {
    int prn;
    double fd;
  };
  std::vector<Hit> hits;
  for (size_t q = 0; q < prns.size(); ++q) {
    const Peak pk = peak(G[q], fg / kRate, 2);
    if (pk.pk / pk.noise >= gamma) hits.push_back({prns[q], dop[pk.j]});
  }
  const int m_full = std::max<int>(1, static_cast<int>(sl.y.size() / (kTCohMs * static_cast<size_t>(py_round_int(fs * 1e-3)))));
  std::vector<Json> sats;
  for (const Hit& h : hits) {
    dvec fine(7);
    for (int i = 0; i < 7; ++i) fine[i] = h.fd + static_cast<double>(i - 3) * (kBinHz / 2);
    const dmat g = search(sl.y, 1, fs, fs_eff, {h.prn}, fine, m_full)[0];
    const Peak pk = peak(g, fs / kRate, 4);
    const double gam = threshold(m_full, kPfaPrn / static_cast<double>(g.cols * fine.size()));
    const double s = pk.pk / pk.noise;
    if (s < gam) continue;
    const size_t n = g.cols;
    const double a0 = std::sqrt(std::max(g(pk.j, (pk.i + n - 1) % n), 0.0)), a1 = std::sqrt(pk.pk),
                 a2 = std::sqrt(std::max(g(pk.j, (pk.i + 1) % n), 0.0));
    const double den = a0 - 2 * a1 + a2;
    const double dlt = den < 0 ? 0.5 * (a0 - a2) / den : 0.0;
    const double a_pk = a1 - 0.25 * (a0 - a2) * dlt;
    const double s_corr = (a_pk * a_pk) / pk.noise;
    const double cn0 = 10 * std::log10(std::max(s_corr - 1, 1e-9) / (kTCohMs * 1e-3)) + kLImplDb;
    double cp = std::fmod(-(static_cast<double>(pk.i) + dlt) * kRate / fs_eff, static_cast<double>(kChips));
    if (cp < 0) cp += kChips;
    const bool multipeak = pk.second >= gam * pk.noise && 10 * std::log10(pk.second / pk.pk) >= -kSecondDb;
    Json o = Json::object();
    o["prn"] = h.prn;
    o["cn0_dbhz"] = py_round(cn0, 1);
    o["doppler_hz"] = py_round_int(fine[pk.j] * (1 + eps));
    o["code_phase_chips"] = py_round(cp, 2);
    o["peak_over_noise"] = py_round(s, 1);
    o["threshold"] = py_round(gam, 2);
    o["second_peak_db"] = py_round(10 * std::log10(std::max(pk.second, 1e-30) / pk.pk), 1);
    o["multipeak"] = multipeak;
    sats.push_back(o);
  }
  std::stable_sort(sats.begin(), sats.end(), [](const Json& a, const Json& b) {
    return a["cn0_dbhz"].get<double>() > b["cn0_dbhz"].get<double>();
  });
  std::vector<Json> kept;   // python _xcorr_guard
  for (const Json& s : sats) {
    bool ghost = false;
    for (const Json& b : kept) {
      const double d = s["doppler_hz"].get<double>() - b["doppler_hz"].get<double>() + 500;
      const double m = d - 1000 * std::floor(d / 1000);   // python % (non-negative)
      const double gap = b["cn0_dbhz"].get<double>() - s["cn0_dbhz"].get<double>();
      ghost = ghost || (gap >= kXcorrDb && std::fabs(m - 500) <= 100) ||
              (gap >= kXcorrNearDb && s["peak_over_noise"].get<double>() < 1.5 * s["threshold"].get<double>());
    }
    if (!ghost) kept.push_back(s);
  }
  Json offset = nullptr;
  if (!ppm && !kept.empty()) {
    dvec d;
    for (const Json& s : kept) d.push_back(s["doppler_hz"].get<double>());
    offset = static_cast<int64_t>(np::median(d));
  }
  const Json itf = interference(ch, f_meas, exclude_hz);
  Json r = Json::object();
  r["system"] = "GPS L1 C/A";
  r["satellites"] = Json(kept);
  r["n_satellites"] = static_cast<int64_t>(kept.size());
  r["receiver_offset_hz_estimate"] = offset;
  Json sr = Json::object();
  sr["prns"] = static_cast<int64_t>(prns.size());
  sr["doppler_span_hz"] = span;
  sr["doppler_bin_hz"] = kBinHz;
  sr["coherent_ms"] = kTCohMs;
  sr["gate_segments"] = kMGate;
  sr["refine_segments"] = m_full;
  sr["receiver_ppm_applied"] = ppm ? Json(*ppm) : Json(nullptr);
  r["search"] = sr;
  r["interference"] = itf;
  r["indicators"] = indicators(kept, itf);
  return r;
}

}  // namespace si::rx::gnss
