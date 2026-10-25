// python receivers/rc_link/gfsk_coherent.py (see the header there for the method).
#include "receivers/rc_link/gfsk_coherent.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

#include "core/fft.hpp"
#include "core/npcompat.hpp"
#include "core/signal.hpp"

namespace si::rx::rc_link::gfsk {
namespace {

constexpr double kPi = 3.14159265358979323846;
const double kC0[32] = {
    1.5493647186e-08, 1.2515359415e-05, 1.1692026245e-04, 7.2607693471e-04, 3.3765011577e-03, 1.2201251491e-02,
    3.5143701476e-02, 8.2708419697e-02, 1.6331527480e-01, 2.7794917074e-01, 4.1790894564e-01, 5.6710940193e-01,
    7.0748264018e-01, 8.2411817347e-01, 9.0807272557e-01, 9.5695397637e-01, 9.7267765168e-01, 9.5676223260e-01,
    9.0767404231e-01, 8.2351009211e-01, 7.0669872824e-01, 5.6622299350e-01, 4.1702312740e-01, 2.7717237872e-01,
    1.6272616787e-01, 8.2330007259e-02, 3.4943224682e-02, 1.2116025896e-02, 3.3481098154e-03, 7.1881877861e-04,
    1.1552896953e-04, 1.2323470268e-05};
constexpr int kMfDelay = 9;
const double kTab[2][2] = {{0.6064282430011425, 1.3034506592017188}, {0.8962552450020563, 1.5725184171543038}};
const double kHGrid[9] = {0.42, 0.44, 0.46, 0.48, 0.5, 0.52, 0.54, 0.56, 0.58};
constexpr int64_t kW = 16;
const cd kAxis[4] = {cd(0, 1), cd(-1, 0), cd(0, -1), cd(1, 0)};
const double kSign[2] = {-1.0, 1.0};

double p1() {
  // np.correlate(C0, C0, "full") at lag 0 and lag SPS
  double a0 = 0, a1 = 0;
  for (int i = 0; i < 32; ++i) a0 += kC0[i] * kC0[i];
  for (int i = 0; i + kSps < 32; ++i) a1 += kC0[i] * kC0[i + kSps];
  return a1 / a0;
}
const double kP1 = p1();

const dvec& lpf() {
  static const dvec h = sig::firwin(33, 0.6e6, 8e6);
  return h;
}

// python _km(v): 2-means on float64 values
std::pair<double, double> km(const dvec& v, int it = 20) {
  double a = np::percentile(v, 25), b = np::percentile(v, 75);
  dvec x, y;
  for (int k = 0; k < it; ++k) {
    x.clear(); y.clear();
    for (double t : v) (std::fabs(t - a) < std::fabs(t - b) ? x : y).push_back(t);
    if (y.empty() || x.empty()) break;
    a = np::pairwise_sum(x.data(), x.size()) / static_cast<double>(x.size());
    b = np::pairwise_sum(y.data(), y.size()) / static_cast<double>(y.size());
  }
  return {std::min(a, b), std::max(a, b)};
}

inline cd cexpj(double x) { return std::exp(cd(0.0, x)); }

// python a // b for ints
inline int64_t fdiv(int64_t a, int64_t b) {
  int64_t q = a / b;
  if ((a % b != 0) && ((a < 0) != (b < 0))) --q;
  return q;
}

// np.linspace(0, n, m + 1).astype(int)
std::vector<int64_t> edges(int64_t n, int m) {
  std::vector<int64_t> e(static_cast<size_t>(m) + 1);
  const double step = static_cast<double>(n) / m;
  for (int i = 0; i < m; ++i) e[i] = static_cast<int64_t>(static_cast<double>(i) * step);
  e[m] = n;
  return e;
}

std::pair<std::vector<int>, double> first_pass(const cdvec& p, int64_t s0, int64_t n, const std::vector<int>* a_in) {
  cdvec q0(n);
  dvec d(n);
  for (int64_t k = 0; k < n; ++k) {
    q0[k] = p[s0 + k * kSps];
    d[k] = std::arg(q0[k]);
  }
  std::vector<int> a(n);
  if (a_in) {
    a = *a_in;
  } else {
    const auto [lo, hi] = km(d);
    const double thr = (lo + hi) / 2;
    for (int64_t k = 0; k < n; ++k) a[k] = d[k] > thr ? 1 : -1;
  }
  dvec t(n);
  for (int64_t k = 0; k < n; ++k) {
    const int prv = k ? a[k - 1] : a[0], nxt = k + 1 < n ? a[k + 1] : a[n - 1];
    t[k] = kTab[prv == a[k]][nxt == a[k]] * a[k];
  }
  const auto e = edges(n, 8);
  double sc[9];
  for (int j = 0; j < 9; ++j) {
    const double hh = kHGrid[j] / 0.5;
    double tot = 0;
    for (int i = 0; i < 8; ++i) {
      cd s(0, 0);
      for (int64_t k = e[i]; k < e[i + 1]; ++k) s += q0[k] * cexpj(-hh * t[k]);
      tot += std::abs(s);
    }
    sc[j] = tot;
  }
  int j = static_cast<int>(std::max_element(sc, sc + 9) - sc);
  double h = kHGrid[j];
  if (j > 0 && j < 8) {
    const double den = sc[j - 1] - 2 * sc[j] + sc[j + 1];
    if (den < 0) h += std::clamp(0.5 * (sc[j - 1] - sc[j + 1]) / den, -0.5, 0.5) * (kHGrid[1] - kHGrid[0]);
  }
  return {a, h};
}

// python _tone_freq(s, pad=8) -> (rad / sample, peak magnitude)
std::pair<double, double> tone_freq(const cd* s, size_t len) {
  const size_t nf = size_t(1) << static_cast<int>(std::ceil(std::log2(static_cast<double>(std::max<size_t>(len, 16) * 8))));
  const cdvec F = fft::fft(s, len, nf);
  dvec m(nf);
  for (size_t i = 0; i < nf; ++i) m[i] = std::abs(F[i]);
  const size_t j = static_cast<size_t>(std::max_element(m.begin(), m.end()) - m.begin());
  const double a = m[(j + nf - 1) % nf], b = m[j], c = m[(j + 1) % nf];
  const double dj = 0.5 * (a - c) / (a - 2 * b + c + 1e-30);
  double x = (static_cast<double>(j) + dj) / static_cast<double>(nf) + 0.5;
  x = x - std::floor(x);
  return {2 * kPi * (x - 0.5), b};
}

// np.polyfit(x, y, deg, w) for deg <= 2 (scaled least squares, as numpy), evaluated at 0..n-1
dvec polyfit_eval(const dvec& x, const dvec& y, const dvec& w, int deg, int64_t n) {
  const size_t m = x.size(), c = static_cast<size_t>(deg) + 1;
  std::vector<dvec> A(m, dvec(c));
  dvec rhs(m);
  for (size_t i = 0; i < m; ++i) {
    for (size_t k = 0; k < c; ++k) A[i][k] = std::pow(x[i], static_cast<double>(deg - static_cast<int>(k))) * w[i];
    rhs[i] = y[i] * w[i];
  }
  dvec scale(c, 0.0);
  for (size_t k = 0; k < c; ++k) {
    for (size_t i = 0; i < m; ++i) scale[k] += A[i][k] * A[i][k];
    scale[k] = std::sqrt(scale[k]);
    for (size_t i = 0; i < m; ++i) A[i][k] /= scale[k];
  }
  // normal equations on the scaled system (c <= 3)
  double N[3][3] = {}, r[3] = {};
  for (size_t i = 0; i < m; ++i)
    for (size_t k = 0; k < c; ++k) {
      r[k] += A[i][k] * rhs[i];
      for (size_t l = 0; l < c; ++l) N[k][l] += A[i][k] * A[i][l];
    }
  double sol[3] = {};
  // Gaussian elimination with partial pivoting
  int idx[3] = {0, 1, 2};
  for (size_t k = 0; k < c; ++k) {
    size_t piv = k;
    for (size_t l = k + 1; l < c; ++l)
      if (std::fabs(N[l][k]) > std::fabs(N[piv][k])) piv = l;
    if (piv != k) {
      for (size_t l = 0; l < c; ++l) std::swap(N[k][l], N[piv][l]);
      std::swap(r[k], r[piv]);
      std::swap(idx[k], idx[piv]);
    }
    for (size_t l = k + 1; l < c; ++l) {
      const double f = N[l][k] / N[k][k];
      for (size_t q = k; q < c; ++q) N[l][q] -= f * N[k][q];
      r[l] -= f * r[k];
    }
  }
  for (int k = static_cast<int>(c) - 1; k >= 0; --k) {
    double s = r[k];
    for (size_t l = static_cast<size_t>(k) + 1; l < c; ++l) s -= N[k][l] * sol[l];
    sol[k] = s / N[k][k];
  }
  dvec coef(c);
  for (size_t k = 0; k < c; ++k) coef[k] = sol[k] / scale[k];
  dvec out(static_cast<size_t>(n));
  for (int64_t i = 0; i < n; ++i) {
    double v = 0;
    for (size_t k = 0; k < c; ++k) v = v * static_cast<double>(i) + coef[k];
    out[i] = v;
  }
  return out;
}

dvec track(const cdvec& s, std::initializer_list<int> stages) {
  const int64_t n = static_cast<int64_t>(s.size());
  dvec ph(n, 0.0);
  cdvec res(n);
  for (int ns : stages) {
    for (int64_t k = 0; k < n; ++k) res[k] = s[k] * cexpj(-ph[k]);
    const auto e = edges(n, ns);
    dvec fr, wt, kc;
    for (int i = 0; i < ns; ++i) {
      if (e[i + 1] - e[i] < 8) continue;
      const auto [w, m] = tone_freq(res.data() + e[i], static_cast<size_t>(e[i + 1] - e[i]));
      fr.push_back(w); wt.push_back(std::sqrt(m)); kc.push_back(static_cast<double>(e[i] + e[i + 1] - 1) / 2);
    }
    if (fr.empty()) break;
    const int deg = std::min<int>(2, static_cast<int>(fr.size()) - 1);
    const dvec fk = deg > 0 ? polyfit_eval(kc, fr, wt, deg, n) : dvec(static_cast<size_t>(n), fr[0]);
    double cs = 0;
    for (int64_t k = 0; k < n; ++k) {
      cs += fk[k];
      ph[k] = ph[k] + cs - fk[0];
    }
  }
  return ph;
}

// np.unwrap
dvec unwrap(const dvec& p) {
  dvec out(p.size());
  if (p.empty()) return out;
  out[0] = p[0];
  double corr = 0;
  for (size_t i = 1; i < p.size(); ++i) {
    const double dd = p[i] - p[i - 1];
    double m = np::mod(dd + kPi, 2 * kPi) - kPi;
    if (m == -kPi && dd > 0) m = kPi;
    double pc = m - dd;
    if (std::fabs(dd) < kPi) pc = 0;
    corr += pc;
    out[i] = p[i] + corr;
  }
  return out;
}

struct Coh {
  std::vector<int> alpha;
  dvec rel;
  double h;
  int dt;
  dvec ph;
};

Coh coherent(const cdvec& r, const cdvec& p, int64_t s0, int64_t n, const std::vector<int>* a_in) {
  const auto [a, h] = first_pass(p, s0, n, a_in);
  dvec walk(n);
  double cw = 0;
  for (int64_t k = 0; k < n; ++k) {
    cw += a[k];
    walk[k] = kPi * (h - 0.5) * cw;
  }
  cdvec isi(n);
  for (int64_t k = 0; k < n; ++k) {
    const int nx = k + 1 < n ? a[k + 1] : a[n - 1];
    isi[k] = cd(1.0, kP1 * static_cast<double>(nx - a[k]));
  }
  const int64_t nr = static_cast<int64_t>(r.size());
  double best_m = -1;
  int best_dt = 0;
  cdvec best_rk, best_s;
  dvec best_ph;
  for (int dt = -3; dt <= 3; ++dt) {
    cdvec rk(n), s(n);
    for (int64_t k = 0; k < n; ++k) {
      const int64_t i = std::clamp<int64_t>(s0 - kMfDelay + dt + k * kSps, 0, nr - 1);
      rk[k] = r[i] * cexpj(-walk[k]);
      const cd u = rk[k] / isi[k];
      s[k] = u * u * kSign[k % 2];
    }
    const dvec ph = track(s, {8, 4});
    cdvec res(n);
    for (int64_t k = 0; k < n; ++k) res[k] = s[k] * cexpj(-ph[k]);
    const double m = tone_freq(res.data(), res.size()).second;
    if (best_m < 0 || m > best_m) {
      best_m = m; best_dt = dt; best_rk = std::move(rk); best_s = std::move(s); best_ph = ph;
    }
  }
  dvec ph = best_ph;
  {
    cdvec res(n);
    for (int64_t k = 0; k < n; ++k) res[k] = best_s[k] * cexpj(-ph[k]);
    const dvec p2 = track(res, {1});
    for (int64_t k = 0; k < n; ++k) ph[k] += p2[k];
  }
  cdvec cs(static_cast<size_t>(n) + 1);
  cs[0] = 0;
  for (int64_t k = 0; k < n; ++k) cs[k + 1] = cs[k] + best_s[k] * cexpj(-ph[k]);
  dvec ang(n);
  for (int64_t k = 0; k < n; ++k) {
    const int64_t hi = std::clamp<int64_t>(k + kW + 1, 0, n), lo = std::clamp<int64_t>(k - kW, 0, n);
    ang[k] = std::arg(cs[hi] - cs[lo]);
  }
  dvec ph2 = unwrap(ang);
  for (int64_t k = 0; k < n; ++k) ph2[k] += ph[k];
  Coh out;
  out.rel.resize(n);
  cdvec bh(n);
  for (int64_t k = 0; k < n; ++k) {
    const cd v = best_rk[k] * cexpj(-0.5 * ph2[k]);
    const cd ax = kAxis[k % 4];
    const double c = std::real(v * std::conj(ax));
    bh[k] = (c >= 0 ? 1.0 : -1.0) * ax;
    out.rel[k] = std::fabs(c);
  }
  out.alpha.resize(n);
  out.alpha[0] = a[0];
  for (int64_t k = 1; k < n; ++k) out.alpha[k] = std::imag(bh[k] * std::conj(bh[k - 1])) > 0 ? 1 : -1;
  out.h = h;
  out.dt = best_dt;
  out.ph = std::move(ph);
  return out;
}

// ---- CRC-24
constexpr int64_t kCrcMask = (1 << 1) | (1 << 3) | (1 << 4) | (1 << 6) | (1 << 9) | (1 << 10);
constexpr size_t kNmax = 16 + 8 * 255 + 24;

inline int64_t clock0(int64_t S) {
  const int64_t fb = (S >> 23) & 1;
  S = ((S << 1) & 0xFFFFFF) | fb;
  return fb ? S ^ kCrcMask : S;
}

struct CrcTables {
  std::vector<int64_t> M, I;
  CrcTables() : M(kNmax), I(kNmax + 1) {
    int64_t S = 1 ^ kCrcMask;
    for (size_t j = 0; j < kNmax; ++j) { M[j] = S; S = clock0(S); }
    S = 0;
    for (int p = 0; p < 24; ++p)
      if ((0x555555 >> p) & 1) S |= int64_t(1) << p;
    for (size_t N = 0; N <= kNmax; ++N) { I[N] = S; S = clock0(S); }
  }
};
const CrcTables& crc_tables() {
  static const CrcTables t;
  return t;
}

}  // namespace

cdvec diff_products(const cfvec& y) {
  const dvec& h = lpf();
  const int64_t L = static_cast<int64_t>(y.size()) + kSps, H = static_cast<int64_t>(h.size()), half = (H - 1) / 2;
  cdvec z(L);
  for (int64_t i = 0; i < L; ++i) {
    cd acc(0, 0);
    for (int64_t j = 0; j < H; ++j) {
      const int64_t m = i + half - j;
      if (m >= 0 && m < static_cast<int64_t>(y.size())) acc += cd(y[m].real(), y[m].imag()) * h[j];
    }
    z[i] = acc;
  }
  cdvec p(static_cast<size_t>(L - kSps));
  for (int64_t i = 0; i + kSps < L; ++i) p[i] = z[i + kSps] * std::conj(z[i]);
  return p;
}

std::vector<std::pair<int64_t, double>> find_access_address(const cdvec& p, const std::vector<uint8_t>& aa_bits,
                                                            int top, double min_ncorr) {
  const int64_t nt = static_cast<int64_t>(aa_bits.size());
  const int64_t L = kSps * (nt - 1) + 1, np_ = static_cast<int64_t>(p.size());
  std::vector<std::pair<int64_t, double>> out;
  if (np_ < L) return out;
  dvec d(np_);
  for (int64_t i = 0; i < np_; ++i) d[i] = std::arg(p[i]);
  const int64_t m = np_ - L + 1;
  dvec nc(m);
  for (int64_t i = 0; i < m; ++i) {
    double num = 0, den = 0;
    for (int64_t k = 0; k < nt; ++k) {
      const double v = d[i + k * kSps];
      num += (aa_bits[k] ? 1.0 : -1.0) * v;
      den += std::fabs(v);
    }
    nc[i] = num / (den + 1e-12);
  }
  std::vector<int64_t> order(m);
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](int64_t x, int64_t y) { return -nc[x] < -nc[y]; });
  for (int64_t i : order) {
    if (nc[i] < min_ncorr || static_cast<int>(out.size()) >= top) break;
    bool far = true;
    for (const auto& o : out)
      if (std::llabs(i - o.first) < 32 * kSps) far = false;
    if (far) out.emplace_back(i, nc[i]);
  }
  return out;
}

std::optional<Demod> demodulate(const cfvec& y, int64_t s0, int64_t n, int iters) {
  const int64_t ny = static_cast<int64_t>(y.size());
  n = std::min<int64_t>(n, fdiv(ny - s0, kSps) - 1);
  if (n < 48 || s0 < 0) return std::nullopt;
  // r[m] = sum_j y[m + j] C0[j], m in [0, len(y) + 2 SPS + 32)
  const int64_t nr = ny + 2 * kSps + 32;
  cdvec r(nr);
  for (int64_t m = 0; m < nr; ++m) {
    cd acc(0, 0);
    for (int j = 0; j < 32; ++j) {
      const int64_t i = m + j;
      if (i < ny) acc += cd(y[i].real(), y[i].imag()) * kC0[j];
    }
    r[m] = acc;
  }
  const cdvec p = diff_products(y);
  const int64_t np_ = static_cast<int64_t>(p.size());
  {
    int64_t best_o = 0;
    double best = -1;
    bool any = false;
    for (int o = -kSps / 2; o <= kSps / 2; ++o) {
      if (s0 + o < 0 || s0 + o + (n - 1) * kSps >= np_) continue;
      dvec v(n);
      for (int64_t k = 0; k < n; ++k) v[k] = std::fabs(std::imag(p[s0 + o + k * kSps]));
      const double mval = np::pairwise_sum(v.data(), v.size()) / static_cast<double>(n);
      if (!any || mval > best) { best = mval; best_o = o; any = true; }
    }
    s0 += best_o;
  }
  Coh c;
  std::vector<int> a_in;
  bool have = false;
  int dt = 0;
  for (int it = 0; it < iters; ++it) {
    c = coherent(r, p, s0, n, have ? &a_in : nullptr);
    a_in = c.alpha;
    have = true;
    dt = c.dt;
    if (dt && s0 + dt >= 0 && fdiv(ny - s0 - dt, kSps) - 1 >= n) {
      s0 += dt;
      dt = 0;
    }
  }
  Demod out;
  out.bits.resize(n);
  for (int64_t k = 0; k < n; ++k) out.bits[k] = c.alpha[k] > 0;
  out.rel = c.rel;
  out.h = c.h;
  const int64_t mid = n / 2;
  const double w2 = c.ph[std::min(mid + 1, n - 1)] - c.ph[mid];
  out.cfo_hz = w2 / 2 / (2 * kPi) * 1e6;
  out.dt = dt;
  return out;
}

std::pair<int64_t, std::vector<int64_t>> crc_syndrome_parts(const std::vector<uint8_t>& dw, size_t n_pdu) {
  const auto& T = crc_tables();
  std::vector<int64_t> contrib(n_pdu + 24);
  for (size_t i = 0; i < n_pdu; ++i) contrib[i] = T.M[n_pdu - 1 - i];
  for (size_t t = 0; t < 24; ++t) contrib[n_pdu + t] = int64_t(1) << (23 - t);
  int64_t syn = 0;
  for (size_t i = 0; i < n_pdu + 24 && i < dw.size(); ++i)
    if (dw[i]) syn ^= contrib[i];
  return {syn ^ T.I[n_pdu], contrib};
}

std::optional<std::pair<std::vector<uint8_t>, int>> chase(const std::vector<uint8_t>& dw, size_t n_pdu,
                                                          const dvec& rel, int64_t off, size_t plo, size_t phi,
                                                          int k_least, int max_units) {
  const auto [syn, contrib] = crc_syndrome_parts(dw, n_pdu);
  if (syn == 0) return std::make_pair(dw, 0);
  const int64_t npos = static_cast<int64_t>(n_pdu) + 24;
  std::vector<int64_t> ks(npos);
  std::iota(ks.begin(), ks.end(), off);
  // the last b has no successor (ISI model, TX ramp-down): always tried
  auto rk = [&](int64_t x) { return x == off + npos - 1 ? -1.0 : rel[x]; };
  std::stable_sort(ks.begin(), ks.end(), [&](int64_t x, int64_t y) { return rk(x) < rk(y); });
  struct Unit { std::vector<int64_t> pos; int64_t v; };
  std::vector<Unit> units;
  for (int64_t k : ks) {
    Unit u{{}, 0};
    for (int64_t q : {k - off, k - off + 1})
      if (q >= 0 && q < npos) u.pos.push_back(q);
    bool prot = false;
    for (int64_t q : u.pos)
      if (q >= static_cast<int64_t>(plo) && q < static_cast<int64_t>(phi)) prot = true;
    if (prot) continue;
    for (int64_t q : u.pos) u.v ^= contrib[q];
    units.push_back(std::move(u));
    if (static_cast<int>(units.size()) >= k_least) break;
  }
  const int nu = static_cast<int>(units.size());
  // itertools.combinations order, sizes 1..max_units
  for (int m = 1; m <= max_units && m <= nu; ++m) {
    std::vector<int> idx(m);
    std::iota(idx.begin(), idx.end(), 0);
    while (true) {
      int64_t v = 0;
      for (int u : idx) v ^= units[u].v;
      if (v == syn) {
        std::vector<uint8_t> out = dw;
        for (int u : idx)
          for (int64_t q : units[u].pos) out[q] ^= 1;
        return std::make_pair(out, m);
      }
      int i = m - 1;
      while (i >= 0 && idx[i] == nu - m + i) --i;
      if (i < 0) break;
      ++idx[i];
      for (int j = i + 1; j < m; ++j) idx[j] = idx[j - 1] + 1;
    }
  }
  return std::nullopt;
}

}  // namespace si::rx::rc_link::gfsk
