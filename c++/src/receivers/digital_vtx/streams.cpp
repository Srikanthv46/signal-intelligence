// python: receivers/digital_vtx/streams.py (same steps, same order; numpy rounding = half to even).
#include "receivers/digital_vtx/streams.hpp"

#include <algorithm>
#include <cmath>

#include "core/fft.hpp"
#include "core/npcompat.hpp"
#include "core/pyfmt.hpp"
#include "receivers/digital_vtx/common.hpp"

namespace si::rx::digital_vtx::streams {

namespace {

constexpr int MIN_SEG = 20;
constexpr int MIN_SC = 200;
constexpr double OCC_MIN = 0.75, CP_MIN = 0.5, BAND_FRAC = 0.9, K1_MAX = 1.43, PAIR_MIN = 0.10, PAIR_NULL = 0.07;

using Mat = std::vector<dvec>;  // rows = symbols
using detail::trunc_i;

int64_t np_round(double v) { return static_cast<int64_t>(std::nearbyint(v)); }  // FE_TONEAREST: half to even

// p / mean(p, axis=1) per row
Mat row_norm(const Mat& P) {
  Mat p = P;
  for (auto& r : p) {
    const double m = np::mean(r);
    for (double& v : r) v /= m;
  }
  return p;
}

dvec col_mean(const Mat& p) {
  dvec m(p[0].size(), 0.0);
  for (const auto& r : p)
    for (size_t k = 0; k < r.size(); ++k) m[k] += r[k];
  for (double& v : m) v /= static_cast<double>(p.size());
  return m;
}

double kurtosis(const Mat& P) {
  const Mat p = row_norm(P);
  const size_t n = p.size(), K = p[0].size();
  dvec s1(K, 0.0), s2(K, 0.0);
  for (const auto& r : p)
    for (size_t k = 0; k < K; ++k) s1[k] += r[k], s2[k] += r[k] * r[k];
  dvec kk(K);
  for (size_t k = 0; k < K; ++k) {
    const double m1 = s1[k] / n, m2 = s2[k] / n;
    kk[k] = m2 / (m1 * m1);
  }
  return np::median(kk);
}

// median over the other axis of the correlation of vectors a_j, b_j (python _rho)
double rho_rows(const Mat& q, size_t r0, size_t r1, size_t step, size_t count) {  // rows (r0 + i*step, r1 + i*step)
  dvec out;
  const size_t K = q[0].size();
  for (size_t i = 0; i < count; ++i) {
    const dvec& a = q[r0 + i * step];
    const dvec& b = q[r1 + i * step];
    const double ma = np::mean(a), mb = np::mean(b);
    dvec ab(K), aa(K), bb(K);
    for (size_t k = 0; k < K; ++k) {
      const double x = a[k] - ma, y = b[k] - mb;
      ab[k] = x * y, aa[k] = x * x, bb[k] = y * y;
    }
    out.push_back(np::sum(ab) / (std::sqrt(np::sum(aa) * np::sum(bb)) + 1e-30));
  }
  return np::median(out);
}

double rho_cols(const Mat& p, size_t c0, size_t c1, size_t count) {  // columns (c0 + 2j, c1 + 2j)
  dvec out;
  const size_t n = p.size();
  for (size_t j = 0; j < count; ++j) {
    const size_t ka = c0 + 2 * j, kb = c1 + 2 * j;
    double ma = 0, mb = 0;
    for (size_t i = 0; i < n; ++i) ma += p[i][ka], mb += p[i][kb];
    ma /= n, mb /= n;
    double ab = 0, aa = 0, bb = 0;
    for (size_t i = 0; i < n; ++i) {
      const double x = p[i][ka] - ma, y = p[i][kb] - mb;
      ab += x * y, aa += x * x, bb += y * y;
    }
    out.push_back(ab / (std::sqrt(aa * bb) + 1e-30));
  }
  return np::median(out);
}

std::pair<double, double> pairs(const Mat& P) {
  const Mat p = row_norm(P);
  const dvec mc = col_mean(p);
  Mat q = p;
  for (auto& r : q)
    for (size_t k = 0; k < r.size(); ++k) r[k] /= mc[k];
  const size_t m = p.size() / 2 * 2;
  const double dt = rho_rows(q, 0, 1, 2, m / 2) - rho_rows(q, 1, 2, 2, (m - 2) / 2);
  // blocks of live subcarriers (a null DC restarts the pairs)
  const double med = np::median(mc);
  std::vector<int> live(mc.size());
  for (size_t k = 0; k < mc.size(); ++k) live[k] = mc[k] > 0.1 * med;
  dvec df;
  size_t k = 0;
  while (k < live.size()) {
    if (!live[k]) {
      ++k;
      continue;
    }
    size_t b = k;
    while (b < live.size() && live[b]) ++b;
    if (b - k >= 64) {
      const size_t n = (b - k) / 2 * 2;
      const double d = rho_cols(p, k, k + 1, n / 2) - rho_cols(p, k + 1, k + 2, (n - 2) / 2);
      df.push_back(std::fabs(d));
    }
    k = b;
  }
  return {dt, df.empty() ? 0.0 : -*std::min_element(df.begin(), df.end())};
}

dvec movsum_r(const dvec& v, size_t k) {
  dvec c(v.size() + 1, 0.0);
  for (size_t i = 0; i < v.size(); ++i) c[i + 1] = c[i] + v[i];
  dvec o(c.size() - k);
  for (size_t i = 0; i < o.size(); ++i) o[i] = c[i + k] - c[i];
  return o;
}

cdvec movsum_c(const cdvec& v, size_t k) {
  cdvec c(v.size() + 1, cd(0, 0));
  for (size_t i = 0; i < v.size(); ++i) c[i + 1] = c[i] + v[i];
  cdvec o(c.size() - k);
  for (size_t i = 0; i < o.size(); ++i) o[i] = c[i + k] - c[i];
  return o;
}

struct Seg {
  std::vector<cdvec> X;  // per symbol, natural FFT order
  std::vector<uint8_t> ok;
};

// python segments(): symbol FFTs of y[a:b] on the folded CP grid; empty when too short / too few CP symbols
bool segments(const cdvec& y, int64_t N, int64_t ncp, double Lf, int64_t a, int64_t b, Seg& out) {
  a = std::max<int64_t>(0, a);
  b = std::min<int64_t>(static_cast<int64_t>(y.size()), b);
  if (b <= a) return false;
  const cd* seg = y.data() + a;
  const int64_t L = b - a;
  if (static_cast<double>(L) < (MIN_SEG + 1) * Lf + N) return false;
  const int64_t nc = L - N;
  cdvec c(static_cast<size_t>(nc));
  dvec e(static_cast<size_t>(nc));
  for (int64_t i = 0; i < nc; ++i) {
    c[i] = seg[i] * std::conj(seg[i + N]);
    e[i] = (np::cabs2(seg[i]) + np::cabs2(seg[i + N])) / 2;
  }
  const cdvec cs = movsum_c(c, static_cast<size_t>(ncp));
  const dvec es = movsum_r(e, static_cast<size_t>(ncp));
  const int64_t lm = static_cast<int64_t>(cs.size());
  dvec m(cs.size());
  for (size_t i = 0; i < m.size(); ++i) m[i] = std::abs(cs[i]) / (es[i] + 1e-30);
  const int64_t n_sym = static_cast<int64_t>(std::floor(static_cast<double>(lm - 1) / Lf));
  const int64_t nph = static_cast<int64_t>(std::ceil(Lf));
  auto grid = [&](int64_t ph, int64_t i) { return std::min(np_round(static_cast<double>(ph) + static_cast<double>(i) * Lf), lm - 1); };
  double best = -1;
  int64_t phase = 0;
  dvec row(static_cast<size_t>(n_sym));
  for (int64_t ph = 0; ph < nph; ++ph) {  // symbol grid: fold the CP metric at the symbol period
    for (int64_t i = 0; i < n_sym; ++i) row[i] = m[grid(ph, i)];
    const double v = np::sum(row);
    if (v > best) best = v, phase = ph;
  }
  const int64_t off = ncp / 2;  // mid-CP window: only |X|^2 is used
  std::vector<int64_t> st;
  for (int64_t i = 0; i < n_sym; ++i) {
    const int64_t s = grid(phase, i);
    if (s + off + N <= L) st.push_back(s);
  }
  std::vector<uint8_t> ok(st.size());
  int nok = 0;
  cdvec acc;
  for (size_t j = 0; j < st.size(); ++j) {
    ok[j] = m[st[j]] >= CP_MIN;
    if (ok[j]) ++nok, acc.push_back(cs[st[j]]);
  }
  if (nok < MIN_SEG) return false;
  const double eps = std::arg(np::sum(acc)) / (2 * M_PI * static_cast<double>(N));
  out.X.clear();
  cdvec w(static_cast<size_t>(N));
  for (int64_t s : st) {
    for (int64_t n = 0; n < N; ++n) {
      const int64_t t = s + off + n;
      w[n] = seg[t] * std::exp(cd(0, 2 * M_PI * eps * static_cast<double>(t)));
    }
    out.X.push_back(fft::fft(w));
  }
  out.ok = ok;
  return true;
}

}  // namespace

Json tx_streams(const cdvec& y, double ofs, double tu, double ts, double bw, const std::vector<blind::Run>& runs,
                double snr_db, double dur_s) {
  Json out = Json::object();
  out["verdict"] = "unknown";
  out["tx_antennas"] = nullptr;
  if (snr_db < SNR_MIN) {
    out["note"] = fmt("link SNR %.1f dB < %.0f dB: noise masks the constellation statistics", snr_db, SNR_MIN);
    return out;
  }
  const int64_t N = np_round(tu * ofs);
  const int64_t ncp = np_round((ts - tu) * ofs);
  const double Lf = ts * ofs;
  if (ncp < 4) {
    out["note"] = "cyclic prefix too short at this rate";
    return out;
  }
  std::vector<size_t> band;  // natural FFT bins of the band, in fftshifted (ascending frequency) order
  for (int64_t k = -(N / 2); k < N - N / 2; ++k)
    if (std::fabs(static_cast<double>(k)) <= BAND_FRAC * bw / 2 * tu) band.push_back(static_cast<size_t>((k + N) % N));
  if (static_cast<int>(band.size()) < MIN_SC) {
    out["note"] = fmt("%d subcarriers < %d: too few for the statistics", static_cast<int>(band.size()), MIN_SC);
    return out;
  }
  struct Res {
    size_t n;
    double kappa, dt, df;
    bool edge;
  };
  std::vector<Res> res;
  for (const auto& r : runs) {
    const int64_t a = trunc_i(r.t0 * ofs) - trunc_i(Lf), b = trunc_i(r.t1 * ofs) + trunc_i(Lf);
    const bool edge = r.t0 < 2 * ts || r.t1 > dur_s - 2 * ts;  // cut off by the capture: parity unknown
    Seg sg;
    if (!segments(y, N, ncp, Lf, a, b, sg)) continue;
    const size_t ns = sg.X.size(), K = band.size();
    Mat P(ns, dvec(K));
    for (size_t i = 0; i < ns; ++i)
      for (size_t k = 0; k < K; ++k) P[i][k] = np::cabs2(sg.X[i][band[k]]);
    dvec ref(K, 0.0);
    size_t nok = 0;
    for (size_t i = 0; i < ns; ++i)
      if (sg.ok[i]) {
        ++nok;
        for (size_t k = 0; k < K; ++k) ref[k] += P[i][k];
      }
    for (double& v : ref) v = v / static_cast<double>(nok) + 1e-30;
    std::vector<size_t> idx;
    for (size_t i = 0; i < ns; ++i) {
      size_t cnt = 0;
      for (size_t k = 0; k < K; ++k) cnt += P[i][k] > 0.1 * ref[k];
      if (sg.ok[i] && static_cast<double>(cnt) / static_cast<double>(K) >= OCC_MIN) idx.push_back(i);
    }
    size_t g0 = 0;
    while (g0 < idx.size()) {  // contiguous runs of kept symbols
      size_t g1 = g0 + 1;
      while (g1 < idx.size() && idx[g1] == idx[g1 - 1] + 1) ++g1;
      if (g1 - g0 >= static_cast<size_t>(MIN_SEG)) {
        Mat Pg(P.begin() + static_cast<ptrdiff_t>(idx[g0]), P.begin() + static_cast<ptrdiff_t>(idx[g1 - 1] + 1));
        const auto [dt, df] = pairs(Pg);
        res.push_back({g1 - g0, kurtosis(Pg), dt, df, edge || (idx[g0] == 0 && r.t0 < 2 * ts)});
      }
      g0 = g1;
    }
  }
  if (res.empty()) {
    out["note"] = fmt("no run of >= %d contiguous full-band data symbols", MIN_SEG);
    return out;
  }
  dvec ks;
  for (const auto& r : res) ks.push_back(r.kappa);
  const double kap = np::median(ks);
  dvec dts, dfs;
  size_t n_inner = 0, nsym = 0;
  for (const auto& r : res) nsym += r.n;
  for (const auto& r : res)
    if (!r.edge) dts.push_back(r.dt), dfs.push_back(r.df), ++n_inner;
  if (!n_inner) {
    for (const auto& r : res) dts.push_back(r.dt), dfs.push_back(r.df);
    n_inner = res.size();
  }
  const double dt = np::median(dts), df = np::median(dfs);
  out["kurtosis"] = py_round(kap, 3);
  out["time_pair_delta"] = py_round(dt, 3);
  out["freq_pair_delta"] = py_round(df, 3);
  out["segments"] = res.size();
  out["symbols"] = nsym;
  out["subcarriers"] = band.size();
  bool all_neg = true, all_pos = true, f_all = true;
  for (double d : dts) all_neg = all_neg && d <= -PAIR_MIN / 2, all_pos = all_pos && d >= PAIR_MIN / 2;
  for (double d : dfs) f_all = f_all && d <= -PAIR_MIN;
  const bool paired_t = std::fabs(dt) >= PAIR_MIN && (all_neg || all_pos);
  if (paired_t && std::fabs(dt) >= std::fabs(df)) {
    out["verdict"] = "2tx_stbc_like";
    out["tx_antennas"] = 2;
    out["note"] = fmt("Alamouti over OFDM symbol pairs (%s pairs from the run start): pair power correlation delta "
                      "%+.2f on %d run(s), kurtosis %.2f",
                      dt < 0 ? "even" : "odd", dt, static_cast<int>(n_inner), kap);
  } else if (f_all && std::fabs(df) >= PAIR_MIN) {
    out["verdict"] = "2tx_sfbc_like";
    out["tx_antennas"] = 2;
    out["note"] = fmt("Alamouti over subcarrier pairs: pair power correlation delta %+.2f, kurtosis %.2f", df, kap);
  } else if (kap <= K1_MAX && std::fabs(dt) < PAIR_NULL && std::fabs(df) < PAIR_NULL) {
    out["verdict"] = "1_stream_like";
    out["note"] = fmt("kurtosis %.2f <= %g: one stream per subcarrier (one transmit antenna, or several with cyclic "
                      "delay / switching), no Alamouti pairing",
                      kap, K1_MAX);
  } else {
    out["note"] = fmt("kurtosis %.2f and pair deltas %+.2f / %+.2f match no verdict: above a single stream without "
                      "the Alamouti pairing (spatial multiplexing, a varying allocation or non-QAM data)",
                      kap, dt, df);
  }
  return out;
}

}  // namespace si::rx::digital_vtx::streams
