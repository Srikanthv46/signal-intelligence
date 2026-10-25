// Emitter association (python: receivers/emitters/association.py). Same method, same
// operation order as numpy (see association.hpp and the Python module docstring).
#include "receivers/emitters/association.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <set>
#include <tuple>

#include "core/fft.hpp"
#include "core/npcompat.hpp"
#include "core/pyfmt.hpp"
#include "core/signal.hpp"

namespace si::rx::emitters {
namespace {

// 15 kHz-SCS OFDM at a 15.36 MSps working rate
constexpr double R15 = 15.36e6;
constexpr size_t N15 = 1024;
constexpr size_t CP15 = 72;
constexpr size_t SYM15 = N15 + CP15;
constexpr double PI2 = 2 * M_PI;   // python 2 * np.pi
constexpr double INF = std::numeric_limits<double>::infinity();
const double NaN = std::numeric_limits<double>::quiet_NaN();

constexpr double EPS_OUTLIER = 0.03;
constexpr double EPS_TIMING_FINE = 0.001;
constexpr double EPS_BROAD = 0.003;
constexpr double EPS_BROAD_FINE = 3e-4;
constexpr double LLR_VETO = -12.0;

// Python max(a, b) / min(a, b): the first argument unless the second compares greater / less
inline double pymax(double a, double b) { return b > a ? b : a; }
inline double pymin(double a, double b) { return b < a ? b : a; }

// np.exp(1j * im) (libm cexp; on Apple that is __sincos)
inline cd expj(double im) {
#ifdef __APPLE__
  double s, c;
  __sincos(im, &s, &c);
  return {c, s};
#else
  return {std::cos(im), std::sin(im)};
#endif
}

inline int64_t pymod(int64_t a, int64_t b) {
  const int64_t r = a % b;
  return r < 0 ? r + b : r;
}

double mean_or_nan(const dvec& v) { return v.empty() ? NaN : np::mean(v); }
double median_or_nan(const dvec& v) { return v.empty() ? NaN : np::median(v); }

// ------------------------------------------------------------------ feature access
const Json& feats(const Emission& e) {
  static const Json empty = Json::object();
  return e.features.is_object() ? e.features : empty;
}
const Json* jget(const Json& o, const char* k) {
  if (!o.is_object()) return nullptr;
  const auto it = o.find(k);
  return it == o.end() ? nullptr : &*it;
}
bool truthy(const Json* v) {
  if (!v || v->is_null()) return false;
  if (v->is_boolean()) return v->get<bool>();
  if (v->is_number()) return v->get<double>() != 0;
  if (v->is_string()) return !v->get<std::string>().empty();
  return !v->empty();
}
bool not_none(const Json* v) { return v && !v->is_null(); }
double num(const Json& f, const char* k, double dflt) {   // float(f.get(k, dflt))
  const Json* v = jget(f, k);
  return v ? v->get<double>() : dflt;
}
std::string tsrc(const Emission& e) {
  const Json* v = jget(feats(e), "timing_source");
  return v && v->is_string() ? v->get<std::string>() : std::string();
}
bool both_fine(const Emission& a, const Emission& b) {
  return tsrc(a) == "reference_symbols" && tsrc(b) == "reference_symbols";
}

bool is_ofdm15(const Emission& e) {
  const Json* v = jget(feats(e), "scs_khz");
  const double scs = not_none(v) ? v->get<double>() : 0.0;
  return std::fabs(scs - 15.0) < 0.5 || e.family == "dji_ocusync" || e.family == "dji_droneid";
}

// =================================================================== measurement
struct CpItem {
  double f_slice = 0, t_cp = 0, sig_t = 0, carrier = 0, sig_nu = 0, snr_db = 0, edge_db = 0;
  int64_t nsym = 0, n_occ = 0, cp_peaks = 0;
  bool dc_null = false;
  // template stage
  std::optional<double> ref_corr;
  bool ref_ok = false;
  double t_fine = 0, sig_t_fine = 0, carrier_fine = 0, sig_nu_fine = 0;
};

double parabolic(const dvec& v, int64_t i) {
  if (i <= 0 || i >= static_cast<int64_t>(v.size()) - 1) return 0.0;
  const double a = v[i - 1], b = v[i], c = v[i + 1];
  const double d = a - 2 * b + c;
  return d == 0 ? 0.0 : 0.5 * (a - c) / d;
}

struct Centre {
  double k = 0;
  int64_t W = 0;
  bool dc = false;
};

std::optional<Centre> integer_centre(const dvec& P, const dvec& kk, double bw_bins, double lev_in, double lev_out,
                                     double k_live) {
  const auto o = np::argsort_stable(kk);
  dvec k, lp;
  for (const int64_t i : o)
    if (std::fabs(kk[i]) <= k_live) {
      k.push_back(kk[i]);
      lp.push_back(10 * std::log10(P[i] + 1e-30));
    }
  const int64_t W = py_round_int(bw_bins);
  const double hw = static_cast<double>(W) / 2;
  double best = -INF;
  std::optional<double> c_best;
  const auto lim = static_cast<int64_t>(0.35 * bw_bins);
  for (const double c : np::arange<double>(static_cast<double>(-lim), static_cast<double>(lim) + 0.01, 0.5)) {
    dvec in, sh;
    for (size_t i = 0; i < k.size(); ++i) {
      const double r = std::fabs(k[i] - c);
      if (r <= hw) in.push_back(lp[i]);
      if (r > hw && r <= hw + 6) sh.push_back(lp[i]);
    }
    if (sh.size() < 4) continue;
    const double sc = mean_or_nan(in) - np::mean(sh);
    if (sc > best) {
      best = sc;
      c_best = c;
    }
  }
  if (!c_best) return std::nullopt;
  dvec depth(k.size(), -INF);
  for (size_t i = 4; i + 4 < k.size(); ++i)
    if (std::fabs(k[i]) <= 0.3 * bw_bins && lp[i] <= lp[i - 1] && lp[i] <= lp[i + 1])
      depth[i] = np::median(dvec{lp[i - 4], lp[i - 3], lp[i - 2], lp[i + 2], lp[i + 3], lp[i + 4]}) - lp[i];
  const auto o2 = np::argsort_quick(depth);
  const int64_t i0 = o2.back();
  const double second = o2.size() >= 2 ? depth[o2[o2.size() - 2]] : -INF;
  const double edge_db = 10 * std::log10(lev_in / pymax(lev_out, 1e-30));
  if (depth[i0] >= pymin(7.0, 0.5 * edge_db) && depth[i0] - pymax(second, 0) >= 2.5) return Centre{k[i0], W, true};
  return Centre{*c_best, W, false};
}

std::optional<CpItem> cp_stage(const SpectralChannelizer& ch, double fc, const Emission& e, double half_bw) {
  const double R = R15;
  const size_t N = N15, cp = CP15, sym = SYM15;
  const double f0 = e.center_mhz * 1e6 - fc;
  const cdvec y = ch.slice(f0, half_bw, R);
  const int64_t ny = static_cast<int64_t>(y.size());
  const int64_t a = std::max<int64_t>(0, static_cast<int64_t>((e.start_us - 150) * 1e-6 * R));
  const int64_t b = std::min<int64_t>(ny, static_cast<int64_t>((e.end_us + 150) * 1e-6 * R));
  if (b - a < static_cast<int64_t>(3 * sym)) return std::nullopt;
  const cd* seg = y.data() + a;
  const size_t L = static_cast<size_t>(b - a) - N;
  // _cp_metric; y[:-N] * conj(y[N:]): numpy elides the conj temporary on >= 256 KiB (operands swap)
  const bool elided = L * sizeof(cd) >= 256 * 1024;
  cdvec c(L);
  dvec en(L);
  for (size_t i = 0; i < L; ++i) {
    const cd v = std::conj(seg[i + N]);
    c[i] = elided ? np::cmul(v, seg[i]) : np::cmul(seg[i], v);
    en[i] = (np::cabs2(seg[i]) + np::cabs2(seg[i + N])) / 2;
  }
  const cdvec C = np::convolve(c, cdvec(cp, cd(1, 0)), np::ConvMode::Valid);
  const dvec den = np::convolve(en, dvec(cp, 1.0), np::ConvMode::Valid);
  dvec m(C.size());
  for (size_t i = 0; i < m.size(); ++i) m[i] = np::cabs(C[i]) / (den[i] + 1e-30);
  const double hi = pymax(0.25, 0.5 * m[np::argmax(m)]);
  sig::PeakOpts opt;
  opt.height_min = hi;
  opt.distance = static_cast<double>(static_cast<int64_t>(0.8 * sym));
  const double s_rel = (e.start_us * 1e-6 * R) - a;
  const double e_rel = (e.end_us * 1e-6 * R) - a;
  std::vector<int64_t> pk;
  for (const int64_t p : sig::find_peaks(m, opt).peaks) {
    const double pd = static_cast<double>(p);
    if (pd >= s_rel - 0.6 * sym && pd <= e_rel - 0.5 * sym) pk.push_back(p);
  }
  if (pk.size() < 2) return std::nullopt;
  const double S = static_cast<double>(sym);
  dvec fr(pk.size()), ph(pk.size());
  for (size_t i = 0; i < pk.size(); ++i) fr[i] = static_cast<double>(pk[i]) + parabolic(m, pk[i]);
  for (size_t i = 0; i < pk.size(); ++i) {
    const double v = np::mod(fr[i] - fr[0], S);
    ph[i] = v > S / 2 ? v - S : v;
  }
  double base = fr[0] + np::median(ph);
  dvec kq(pk.size()), res(pk.size());
  std::vector<size_t> good;
  for (size_t i = 0; i < pk.size(); ++i) {
    kq[i] = std::nearbyint((fr[i] - base) / S);
    res[i] = fr[i] - (base + kq[i] * S);
    if (std::fabs(res[i]) < 0.08 * S) good.push_back(i);
  }
  if (good.size() < 2) return std::nullopt;
  dvec rg;
  for (const size_t i : good) rg.push_back(res[i]);
  base = base + np::mean(rg);
  for (size_t j = 0; j < good.size(); ++j) rg[j] = fr[good[j]] - (base + kq[good[j]] * S);
  const auto kmin = static_cast<int64_t>(std::nearbyint((s_rel - base) / S));
  CpItem it;
  it.t_cp = (static_cast<double>(a) + base + static_cast<double>(kmin * static_cast<int64_t>(sym))) / R;
  int64_t nsym = static_cast<int64_t>(
      std::nearbyint((e_rel - (base + static_cast<double>(kmin * static_cast<int64_t>(sym)))) / S));
  const double ng = static_cast<double>(good.size());
  it.sig_t = pymax(good.size() > 2 ? np::stddev(rg) : 0.25 * cp, 0.05 * cp) / std::sqrt(ng) / R;
  // fine CFO from the CP correlation phase (mod 15 kHz)
  cdvec cc, cn;
  for (const size_t i : good) cc.push_back(C[pk[i]]);
  for (const cd& v : cc) cn.push_back(np::cdiv(v, cd(np::cabs(v) + 1e-30, 0)));
  const double phase = std::arg(np::sum(cc));
  const double nu = -phase / (PI2 * N / R);
  const double circ = np::cabs(np::sum(cn)) / static_cast<double>(cc.size());
  const double sig_phi = std::sqrt(pymax(-2 * std::log(pymax(circ, 1e-6)), 1e-4)) / std::sqrt(static_cast<double>(cc.size()));
  it.sig_nu = sig_phi / (PI2 * N / R);
  // integer subcarrier: occupied bins after removing the fine CFO
  const cdvec y2 = ch.slice(f0 + nu, half_bw * 1.6, R);
  nsym = std::max<int64_t>(nsym, 1);
  dvec P(N, 0.0);
  for (int64_t j = 0; j < nsym; ++j) {
    const int64_t s = py_round_int(it.t_cp * R) + j * static_cast<int64_t>(sym) + static_cast<int64_t>(cp / 2);
    if (s >= 0 && s + static_cast<int64_t>(N) <= static_cast<int64_t>(y2.size())) {
      const cdvec X = fft::fft(y2.data() + s, N);
      for (size_t i = 0; i < N; ++i) {
        const double v = np::cabs(X[i]);
        P[i] += v * v;
      }
    }
  }
  dvec kk = fft::fftfreq(N);
  for (double& v : kk) v *= static_cast<double>(N);
  const double bw_bins = (e.f_hi_mhz - e.f_lo_mhz) * 1e6 / (R / N);
  dvec pin, pout;
  for (size_t i = 0; i < N; ++i) {
    const double ak = std::fabs(kk[i]);
    if (ak < 0.3 * bw_bins) pin.push_back(P[i]);
    if (ak > 0.75 * bw_bins && ak < 1.5 * half_bw / (R / N)) pout.push_back(P[i]);
  }
  const double lev_in = median_or_nan(pin);
  const double lev_out = pout.empty() ? lev_in * 1e-3 : np::median(pout);
  it.snr_db = 10 * std::log10(pymax(lev_in / pymax(lev_out, 1e-30) - 1, 1e-3));
  const auto kc = integer_centre(P, kk, bw_bins, lev_in, lev_out, 1.6 * half_bw / (R / N) - 3);
  if (!kc) return std::nullopt;
  it.carrier = fc + f0 + nu + kc->k * R / N;
  it.f_slice = f0 + nu + kc->k * R / N;
  it.nsym = nsym;
  it.n_occ = kc->W;
  it.dc_null = kc->dc;
  it.cp_peaks = static_cast<int64_t>(good.size());
  it.edge_db = 10 * std::log10(lev_in / pymax(lev_out, 1e-30));
  return it;
}

// FFT of each symbol (window starting mid-CP), occupied bins |k| <= kmax only
std::vector<cdvec> symbols(const SpectralChannelizer& ch, double f_slice, double t_cp, int64_t nsym, double half_bw,
                           int64_t kmax, dvec* k_out) {
  const double R = R15;
  const size_t N = N15;
  const cdvec y = ch.slice(f_slice, half_bw, R);
  const double s0 = t_cp * R;
  dvec kk = fft::fftfreq(N);
  for (double& v : kk) v *= static_cast<double>(N);
  std::vector<size_t> sel;
  for (size_t i = 0; i < N; ++i)
    if (std::fabs(kk[i]) <= static_cast<double>(kmax)) sel.push_back(i);
  if (k_out) {
    k_out->clear();
    for (const size_t i : sel) k_out->push_back(kk[i]);
  }
  std::vector<cdvec> S;
  for (int64_t j = 0; j < nsym; ++j) {
    const double s = s0 + static_cast<double>(j * static_cast<int64_t>(SYM15)) + CP15 / 2.0;
    const auto a = static_cast<int64_t>(std::floor(s));
    const double frac = s - std::floor(s);
    if (a + static_cast<int64_t>(N) + 1 > static_cast<int64_t>(y.size()) || a < 0) break;
    const cdvec X = fft::fft(y.data() + a, N);
    cdvec row;
    row.reserve(sel.size());
    for (const size_t i : sel) row.push_back(np::cmul(X[i], expj(PI2 * kk[i] * frac / N)));
    S.push_back(std::move(row));
  }
  return S;
}

// argmax of |ifft(Z)| over +-6 us: (delay in samples @R, peak index)
std::pair<double, size_t> coarse_delay(const cdvec& r, size_t up) {
  const size_t L = N15 * up;
  const auto lim = static_cast<size_t>(6e-6 * R15 * static_cast<double>(up));
  size_t best = 0;
  double bv = -1;
  auto visit = [&](size_t i) {
    const double v = np::cabs(r[i]);
    if (v > bv) {
      bv = v;
      best = i;
    }
  };
  for (size_t i = 0; i < lim; ++i) visit(i);
  for (size_t i = L - lim; i < L; ++i) visit(i);
  const double d = (best < L / 2 ? static_cast<double>(best) : static_cast<double>(best) - static_cast<double>(L)) /
                   static_cast<double>(up);
  return {d, best};
}

cdvec spread_bins(const cdvec& v, const dvec& k, size_t up) {
  const auto L = static_cast<int64_t>(N15 * up);
  cdvec Z(static_cast<size_t>(L), cd(0, 0));
  for (size_t i = 0; i < k.size(); ++i) Z[pymod(static_cast<int64_t>(k[i]), L)] += v[i];
  return Z;
}

double delay_xcorr_rho(const cdvec& A, const cdvec& B, const dvec& k) {
  const size_t up = 16;
  cdvec prod(A.size());
  for (size_t i = 0; i < A.size(); ++i) prod[i] = np::cmul(A[i], std::conj(B[i]));
  const cdvec r = fft::ifft(spread_bins(prod, k, up));
  const size_t i = coarse_delay(r, up).second;
  return np::cabs(r[i]) * static_cast<double>(N15) * static_cast<double>(up) / (np::norm(A) * np::norm(B) + 1e-30);
}

struct Fit {
  double d = 0, dnu = 0, sig_d = 0, sig_dnu = 0, rho = 0;
};

cdvec ramp_of(const dvec& k, double d) {
  cdvec r(k.size());
  for (size_t i = 0; i < k.size(); ++i) r[i] = expj(PI2 * k[i] * d / N15);
  return r;
}

double pw_sq_dev(const dvec& v) {
  const double mu = np::mean(v);
  dvec t(v.size());
  for (size_t i = 0; i < v.size(); ++i) t[i] = (v[i] - mu) * (v[i] - mu);
  return np::sum(t);
}

Fit fit_delay_phase(const std::vector<cdvec>& prods, const dvec& k, const dvec& tj) {
  const double N = static_cast<double>(N15);
  const size_t nk = k.size();
  cdvec P = prods[0];
  for (size_t r = 1; r < prods.size(); ++r)
    for (size_t i = 0; i < nk; ++i) P[i] += prods[r][i];
  double d = coarse_delay(fft::ifft(spread_bins(P, k, 16)), 16).first;
  cdvec Q(nk), q(nk);
  auto derot = [&]() {
    const cdvec ramp = ramp_of(k, d);
    for (size_t i = 0; i < nk; ++i) Q[i] = np::cmul(P[i], ramp[i]);
    const cd e = expj(-std::arg(np::sum(Q)));
    for (size_t i = 0; i < nk; ++i) q[i] = np::cmul(Q[i], e);
  };
  for (int it = 0; it < 2; ++it) {
    derot();
    dvec w(nk), ang(nk);
    for (size_t i = 0; i < nk; ++i) {
      w[i] = np::cabs(q[i]) + 1e-30;
      ang[i] = std::arg(q[i]);
    }
    const double ws = np::sum(w);
    dmat A(nk, 2);
    dvec bb(nk);
    for (size_t i = 0; i < nk; ++i) {
      const double sw = std::sqrt(w[i] / ws);
      A(i, 0) = k[i] * sw;
      A(i, 1) = 1.0 * sw;
      bb[i] = ang[i] * sw;
    }
    const double sl = np::lstsq(A, bb).x[0];
    d -= sl * N / PI2;
  }
  derot();
  dvec aq(nk), wr(nk);
  for (size_t i = 0; i < nk; ++i) {
    const double resid = std::arg(q[i]);
    aq[i] = np::cabs(q[i]);
    wr[i] = aq[i] * (resid * resid);
  }
  Fit f;
  f.d = d;
  const double sig_sc = std::sqrt(np::sum(wr) / np::sum(aq));
  f.sig_d = sig_sc / (PI2 / N * std::sqrt(pw_sq_dev(k)));
  const cdvec ramp = ramp_of(k, d);
  dvec ph;
  dvec absall;
  for (const cdvec& p : prods) {
    cdvec pr(nk);
    for (size_t i = 0; i < nk; ++i) {
      pr[i] = np::cmul(p[i], ramp[i]);
      absall.push_back(np::cabs(p[i]));
    }
    ph.push_back(std::arg(np::sum(pr)));
  }
  f.rho = np::cabs(np::sum(Q)) / (np::sum(absall) + 1e-30);
  if (prods.size() >= 2) {
    const dvec phu = np::unwrap(ph);
    f.dnu = np::polyfit(tj, phu, 1)[0] / PI2;
    const double sig_ph = sig_sc / std::sqrt(static_cast<double>(nk)) * std::sqrt(static_cast<double>(prods.size()));
    f.sig_dnu = sig_ph / (PI2 * std::sqrt(pw_sq_dev(tj)));
  } else {
    f.dnu = 0.0;
    f.sig_dnu = INF;
  }
  return f;
}

Json round_list(const dvec& v, int nd) {
  Json out = Json::array();
  for (const double x : v) out.push_back(np::round(x, nd));
  return out;
}

// Learn the symbols that repeat burst-to-burst and re-time every burst against them.
Json template_stage(const SpectralChannelizer& ch, std::vector<CpItem>& items, double half_bw,
                    double min_ref_corr = 0.5, double min_burst_corr = 0.35) {
  const double R = R15, N = static_cast<double>(N15);
  if (items.size() < 2) return Json{{"ref_symbols", Json::array()}, {"note", "single burst: no learned reference"}};
  dvec nocc, ns;
  for (const auto& it : items) {
    nocc.push_back(static_cast<double>(it.n_occ));
    ns.push_back(static_cast<double>(it.nsym));
  }
  const auto kmax = static_cast<int64_t>(std::floor(np::median(nocc) / 2));
  const auto nsym = static_cast<int64_t>(pymin(np::median(ns), *std::min_element(ns.begin(), ns.end())));
  if (nsym < 2) return Json{{"ref_symbols", Json::array()}, {"note", "bursts too short"}};
  const size_t nb = items.size();
  std::vector<std::vector<cdvec>> S;
  dvec k;
  for (const auto& it : items) S.push_back(symbols(ch, it.f_slice, it.t_cp, nsym, half_bw, kmax, &k));
  for (const auto& s : S)
    if (static_cast<int64_t>(s.size()) < nsym)
      return Json{{"ref_symbols", Json::array()}, {"note", "burst runs off the capture"}};
  dvec neg;
  for (const auto& it : items) neg.push_back(-it.snr_db);
  const auto order = np::argsort_quick(neg);
  std::vector<std::pair<int64_t, int64_t>> pairs;
  for (size_t j = 1; j < order.size(); ++j) pairs.emplace_back(order[0], order[j]);
  for (size_t j = 2; j < std::min<size_t>(order.size(), 6); ++j) pairs.emplace_back(order[1], order[j]);
  dvec rho(static_cast<size_t>(nsym));
  for (int64_t jj = 0; jj < nsym; ++jj) {
    dvec v;
    for (const auto& [a, b] : pairs) v.push_back(delay_xcorr_rho(S[a][jj], S[b][jj], k));
    rho[jj] = np::median(v);
  }
  std::vector<int64_t> J;
  for (int64_t j = 0; j < nsym; ++j)
    if (rho[j] > min_ref_corr) J.push_back(j);
  if (J.empty())
    return Json{{"ref_symbols", Json::array()}, {"symbol_repeat_corr", round_list(rho, 2)}, {"note", "no repeating symbol"}};
  const size_t nJ = J.size(), nk = k.size();
  dvec tj(nJ);
  for (size_t i = 0; i < nJ; ++i) tj[i] = static_cast<double>(J[i] * static_cast<int64_t>(SYM15)) / R;
  std::vector<cdvec> T;
  for (const int64_t j : J) T.push_back(S[order[0]][j]);
  auto products = [&](const std::vector<cdvec>& Sb) {
    std::vector<cdvec> prods(nJ, cdvec(nk));
    for (size_t i = 0; i < nJ; ++i)
      for (size_t m = 0; m < nk; ++m) prods[i][m] = np::cmul(Sb[J[i]][m], std::conj(T[i][m]));
    return prods;
  };
  // integer-subcarrier repair: re-cut a non-correlating burst +-1..3 subcarriers off
  int64_t n_fixed = 0;
  for (size_t b = 0; b < nb; ++b) {
    const double r0 = fit_delay_phase(products(S[b]), k, tj).rho;
    if (r0 >= 0.5) continue;
    double br = r0;
    int bm = 0;
    std::vector<cdvec> bs;
    for (const int m : {-1, 1, -2, 2, -3, 3}) {
      auto sm = symbols(ch, items[b].f_slice + m * R / N, items[b].t_cp, nsym, half_bw, kmax, nullptr);
      if (static_cast<int64_t>(sm.size()) < nsym) continue;
      const double r = fit_delay_phase(products(sm), k, tj).rho;
      if (r > br) {
        br = r;
        bm = m;
        bs = std::move(sm);
      }
    }
    if (bm && br > 0.5 && br > 1.5 * r0) {
      items[b].f_slice += bm * R / N;
      items[b].carrier += bm * R / N;
      S[b] = std::move(bs);
      ++n_fixed;
    }
  }
  std::vector<Fit> fits;
  const double tjm = np::mean(tj);
  for (int it = 0; it < 3; ++it) {
    fits.clear();
    std::vector<cdvec> acc(nJ, cdvec(nk, cd(0, 0)));
    for (size_t b = 0; b < nb; ++b) {
      const auto prods = products(S[b]);
      const Fit f = fit_delay_phase(prods, k, tj);
      fits.push_back(f);
      const cdvec ramp = ramp_of(k, f.d);
      cdvec flat;
      flat.reserve(nJ * nk);
      for (const auto& p : prods)
        for (size_t m = 0; m < nk; ++m) flat.push_back(np::cmul(p[m], ramp[m]));
      const double ph0 = std::arg(np::sum(flat));
      for (size_t i = 0; i < nJ; ++i) {
        const cd e = expj(-(ph0 + PI2 * f.dnu * (tj[i] - tjm)));
        for (size_t m = 0; m < nk; ++m) acc[i][m] += np::cmul(np::cmul(S[b][J[i]][m], ramp[m]), e);
      }
    }
    for (size_t i = 0; i < nJ; ++i)
      for (size_t m = 0; m < nk; ++m) T[i][m] = np::cdiv(acc[i][m], cd(static_cast<double>(nb), 0));
  }
  // reference the template timing / CFO to the group mean of the CP estimates
  dvec t_new(nb), dt_cp(nb), w(nb), car_new(nb), wf(nb), dcar(nb);
  for (size_t b = 0; b < nb; ++b) {
    t_new[b] = items[b].t_cp + fits[b].d / R;
    w[b] = 1 / (items[b].sig_t * items[b].sig_t);
    dt_cp[b] = w[b] * (t_new[b] - items[b].t_cp);
    car_new[b] = items[b].carrier + fits[b].dnu;
    wf[b] = 1 / (items[b].sig_nu * items[b].sig_nu);
    dcar[b] = wf[b] * (car_new[b] - items[b].carrier);
  }
  const double bias_t = np::sum(dt_cp) / np::sum(w);
  const double bias_f = np::sum(dcar) / np::sum(wf);
  dvec tdiff(nb);
  for (size_t b = 0; b < nb; ++b) {
    CpItem& it = items[b];
    tdiff[b] = t_new[b] - it.t_cp;
    const double dt = t_new[b] - bias_t - it.t_cp;
    const double df = car_new[b] - bias_f - it.carrier;
    const bool ok = fits[b].rho > min_burst_corr && std::fabs(dt) < pymax(5 * it.sig_t, 2e-6) &&
                    std::fabs(df) < pymax(5 * it.sig_nu, 1200.0);
    it.ref_corr = py_round(fits[b].rho, 3);
    if (ok) {
      it.t_fine = t_new[b] - bias_t;
      it.sig_t_fine = std::hypot(fits[b].sig_d / R, 3e-9);
      if (fits[b].rho >= 0.6) {
        it.carrier_fine = car_new[b] - bias_f;
        it.sig_nu_fine = std::hypot(fits[b].sig_dnu, 2.0);
      } else {
        it.carrier_fine = it.carrier;
        it.sig_nu_fine = it.sig_nu;
      }
    }
    it.ref_ok = ok;
  }
  return Json{{"ref_symbols", J},
              {"symbol_repeat_corr", round_list(rho, 2)},
              {"integer_fixes", n_fixed},
              {"cp_vs_ref_timing_rms_us", py_round(np::stddev(tdiff) * 1e6, 3)}};
}

// =================================================================== physics models
Model model_for(const Emission& e, const Models* models) {
  const ModelKey key{e.family, e.label};
  Model m;
  const Models& src = models ? *models : MODELS();
  if (auto it = src.find(key); it != src.end()) m = it->second;
  else if (auto it2 = MODELS().find(key); it2 != MODELS().end()) m = it2->second;
  const Json* fp = jget(feats(e), "frame_period_ms");
  if (truthy(fp) && m.lattices.empty()) m.lattices = {Lattice{fp->get<double>() * 1e3, {0.0}, 1.0}};
  return m;
}

std::pair<double, double> t_of(const Emission& e) {
  const Json& f = feats(e);
  if (jget(f, "t_start_us")) return {f["t_start_us"].get<double>(), num(f, "t_sigma_us", 0.5)};
  if (jget(f, "frame_phase_us") && jget(f, "t_sigma_us")) return {e.start_us, f["t_sigma_us"].get<double>()};
  return {e.start_us, 5.0};
}

std::optional<std::pair<double, double>> carrier_of(const Emission& e) {
  const Json& f = feats(e);
  if (jget(f, "carrier_hz")) return std::make_pair(f["carrier_hz"].get<double>(), num(f, "carrier_sigma_hz", 200.0));
  return std::nullopt;
}

struct Ppm {
  std::optional<double> p, s;
  Json src = nullptr;
};

Ppm ppm_of(const Emission& e, const Model& m) {
  const Json& f = feats(e);
  const auto c = carrier_of(e);
  const bool refined = c && m.raster && tsrc(e) == "reference_symbols";
  if (not_none(jget(f, "cfo_ppm")) && !refined)
    return {f["cfo_ppm"].get<double>(), num(f, "cfo_ppm_sigma", 0.05), "detector"};
  if (c && m.raster) {
    const auto [off, step] = *m.raster;
    const double nom = off + static_cast<double>(py_round_int((c->first - off) / step)) * step;
    return {(c->first - nom) / nom * 1e6, c->second / nom * 1e6, fmt("raster %g MHz", step / 1e6)};
  }
  const Json* ref = jget(f, "cfo_ref");
  if (not_none(jget(f, "cfo_hz")) && ref && ref->is_string() && ref->get<std::string>() == "nominal") {
    const double cen = e.center_mhz * 1e6;
    return {f["cfo_hz"].get<double>() / cen * 1e6, num(f, "cfo_sigma_hz", 200.0) / cen * 1e6, "detector cfo_hz"};
  }
  return {};
}

Json ids_of(const Emission& e) {
  Json out = Json::object();
  for (const Json* src : {&e.decoded, &feats(e)}) {
    if (!src->is_object()) continue;
    for (const char* k : {"serial", "serial_number", "drone_serial", "mac", "uuid"})
      if (truthy(jget(*src, k))) {
        std::string s = py_str((*src)[k]);
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char ch) { return std::tolower(ch); });
        out["id"] = s;
      }
    for (const char* k : {"zc_root", "pilot_id", "sync_word", "bind_id"})
      if (not_none(jget(*src, k))) out[k] = (*src)[k];
  }
  return out;
}

double bw_of(const Emission& e) { return num(feats(e), "bw_mhz", e.f_hi_mhz - e.f_lo_mhz); }

// ------------------------------------------------------------------ cue likelihoods
const double SQ2PI = std::sqrt(2 * M_PI);

double npdf(double r, double s) {
  const double z = r / s;
  return std::exp(-0.5 * (z * z)) / (s * SQ2PI);
}
double lnpdf(double r, double s) {
  const double z = r / s;
  return -0.5 * (z * z) - std::log(s * SQ2PI);
}
struct Term {
  double w, mu, sg;
};
double lmix(double d, const std::vector<Term>& terms) {
  dvec v;
  for (const auto& t : terms)
    if (t.w > 0) v.push_back(std::log(t.w) + lnpdf(d - t.mu, t.sg));
  const double m = *std::max_element(v.begin(), v.end());
  double s = 0;
  for (const double x : v) s += std::exp(x - m);
  return m + std::log(s);
}

struct LatticeFit {
  double ratio, res, d;
};
LatticeFit lattice_ratio(double dt, double sig, double frame, const dvec& offsets) {
  std::set<double> Dset;
  for (const double a : offsets)
    for (const double b : offsets) Dset.insert(np::mod(b - a, frame));
  const double ph = np::mod(dt, frame);
  double tot = 0.0, br = INF, bd = 0.0;
  for (const double d : Dset)
    for (int kk = -3; kk <= 3; ++kk) {
      const double r = ph - d + kk * frame;
      tot += npdf(r, sig);
      if (std::fabs(r) < std::fabs(br)) {
        br = r;
        bd = d;
      }
    }
  return {tot * frame / static_cast<double>(Dset.size()), br, bd};
}

std::optional<Cue> cue_timing(const Emission& a, const Emission& b, const Model& ma, std::optional<double> ppm_pair) {
  if (ma.lattices.empty()) return std::nullopt;
  const auto [ta, sa] = t_of(a);
  const auto [tb, sb] = t_of(b);
  const double dt = tb - ta;
  const double e = ppm_pair ? *ppm_pair : 0.0, se = ppm_pair ? 1.0 : 25.0;
  const double dfm = std::fabs(a.center_mhz - b.center_mhz);
  const bool same_ch = dfm < 0.5 * pymax(bw_of(a), 0.2);
  const double cdl = ma.chan_delay_us_per_mhz * dfm, dr = std::fabs(dt) * se * 1e-6;
  const double sig = std::sqrt(sa * sa + sb * sb + ma.jitter_us * ma.jitter_us + (same_ch ? 0 : cdl * cdl) + dr * dr);
  const double eps = both_fine(a, b) ? EPS_TIMING_FINE : EPS_OUTLIER;
  double ratio = 0.0;
  std::optional<LatticeFit> best;
  double best_frame = 0;
  for (const auto& L : ma.lattices) {
    const double fr = L.frame_us * (1 - e * 1e-6);
    dvec offs;
    for (const double o : L.offsets_us) offs.push_back(o * (1 - e * 1e-6));
    const LatticeFit lf = lattice_ratio(dt, sig, fr, offs);
    ratio += L.weight * ((1 - eps) * lf.ratio + eps);
    if (!best || L.frame_us < best_frame) {
      best = lf;
      best_frame = L.frame_us;
    }
  }
  const double llr = std::log(pymax(ratio, 1e-300));
  const double frame = best_frame, res = best->res, d = best->d;
  const int64_t n = frame != 0 ? py_round_int((dt - d) / frame) : 0;
  double W = INF;
  for (const auto& L : ma.lattices) W = pymin(W, L.frame_us / static_cast<double>(L.offsets_us.size()));
  Cue c{"timing_lattice", llr,
        fmt("dt %+.4f ms = %lld x %g ms + %g ms %+.0f ns (sigma %.0f ns%s); an independent radio's phase is uniform "
            "over ~%g ms",
            dt / 1e3, static_cast<long long>(n), frame / 1e3, d / 1e3, res * 1e3, sig * 1e3,
            same_ch ? "" : ", cross-channel", W / 1e3)};
  c.extra = Json{{"residual_us", res}, {"sigma_us", sig}};
  return c;
}

std::optional<Cue> cue_oscillator(const Emission& a, const Emission& b, const Model& ma) {
  const Ppm pa = ppm_of(a, ma), pb = ppm_of(b, ma);
  const auto ca = carrier_of(a), cb = carrier_of(b);
  const double dtau = std::fabs(t_of(b).first - t_of(a).first) * 1e-6;
  double d, s;
  std::string how;
  if (!pa.p && ca && cb && std::fabs(ca->first - cb->first) < 100e3) {
    d = (cb->first - ca->first) / ca->first * 1e6;
    s = std::hypot(ca->second, cb->second) / ca->first * 1e6;
    how = "same-channel carrier difference";
  } else if (pa.p && pb.p) {
    d = *pb.p - *pa.p;
    s = std::hypot(*pa.s, *pb.s);
    how = "carrier ppm vs " + py_str(pa.src);
  } else {
    return std::nullopt;
  }
  const double dr = ma.ppm_drift_per_s * dtau;
  const double s1 = std::sqrt(s * s + ma.ppm_h1 * ma.ppm_h1 + dr * dr);
  const double s0 = std::sqrt(s * s + 2 * (ma.ppm_pop * ma.ppm_pop));
  const Json* sv = jget(feats(a), "scs_khz");
  const double scs = (truthy(sv) ? sv->get<double>() : 0.0) * 1e3;
  std::vector<Term> out{{both_fine(a, b) ? EPS_BROAD_FINE : EPS_BROAD, 0.0, 10.0 * ma.ppm_pop}};
  if (scs != 0) {
    const double step = scs / (a.center_mhz * 1e6) * 1e6;
    for (const int k : {-2, -1, 1, 2}) out.push_back({EPS_OUTLIER / 4, k * step, s1});
  }
  std::vector<Term> h1{{1 - EPS_OUTLIER, 0.0, s1}}, h0{{1 - EPS_OUTLIER, 0.0, s0}};
  h1.insert(h1.end(), out.begin(), out.end());
  h0.insert(h0.end(), out.begin(), out.end());
  const double llr = lmix(d, h1) - lmix(d, h0);
  std::string det = pa.p ? fmt("%s: %+.3f vs %+.3f ppm, ", how.c_str(), *pa.p, *pb.p) : how + ": ";
  det += fmt("diff %+.3f ppm (%+.0f Hz); one crystal => 0 +- %.3f, independent crystals => 0 +- %.2f ppm", d,
             d * a.center_mhz, s1, s0);
  Cue c{"oscillator_ppm", llr, det};
  c.extra = Json{{"dppm", d}};
  return c;
}

std::optional<Cue> cue_clock(const Emission& a, const Emission& b, const Model& ma) {
  const Json &fa = feats(a), &fb = feats(b);
  if (!not_none(jget(fa, "clock_ppm")) || !not_none(jget(fb, "clock_ppm"))) return std::nullopt;
  if (a.family.rfind("rc_", 0) == 0 || b.family.rfind("rc_", 0) == 0) return std::nullopt;
  const double d = fb["clock_ppm"].get<double>() - fa["clock_ppm"].get<double>();
  auto sig = [](const Json& f) {   // float(f.get("clock_ppm_sigma", f.get("clock_ppm_err", 2.0)) or 2.0)
    const Json* v = jget(f, "clock_ppm_sigma");
    if (!v) v = jget(f, "clock_ppm_err");
    return truthy(v) ? v->get<double>() : 2.0;
  };
  const double s = std::hypot(sig(fa), sig(fb));
  const double s1 = std::hypot(s, 0.2), s0 = std::sqrt(s * s + 2 * (ma.ppm_pop * ma.ppm_pop));
  const Term out{EPS_OUTLIER, 0.0, 10 * ma.ppm_pop};
  return Cue{"sample_clock_ppm",
             lmix(d, {{1 - EPS_OUTLIER, 0.0, s1}, out}) - lmix(d, {{1 - EPS_OUTLIER, 0.0, s0}, out}),
             fmt("symbol-clock ppm diff %+.2f (sigma %.2f)", d, s)};
}

std::optional<Cue> cue_overlap(const Emission& a, const Emission& b, const Model& ma) {
  const double ov = std::min(a.end_us, b.end_us) - std::max(a.start_us, b.start_us);
  const double shrt = pymax(pymin(a.end_us - a.start_us, b.end_us - b.start_us), 1.0);
  const bool sep = std::fabs(a.center_mhz - b.center_mhz) > 0.5 * (bw_of(a) + bw_of(b)) / 2;
  if (ma.single_tx && ov > 0.2 * shrt && sep)
    return Cue{"simultaneous", -20.0,
               fmt("on air at the same time (%.0f us overlap) on channels %.2f MHz apart: one transmit chain cannot",
                   ov, std::fabs(a.center_mhz - b.center_mhz))};
  return std::nullopt;
}

std::optional<Cue> cue_power(const Emission& a, const Emission& b, const Model& ma) {
  if (!a.snr_db || !b.snr_db) return std::nullopt;
  const double d = *b.snr_db - *a.snr_db;
  const double s1 = ma.snr_h1_db, s0 = 10.0;
  const Term out{EPS_OUTLIER, 0.0, 25.0};
  return Cue{"power", lmix(d, {{1 - EPS_OUTLIER, 0.0, s1}, out}) - lmix(d, {{1 - EPS_OUTLIER, 0.0, s0}, out}),
             fmt("SNR %.1f vs %.1f dB (hop-to-hop fading ~%g dB; different ranges ~%g dB)", *a.snr_db, *b.snr_db, s1,
                 s0)};
}

std::optional<Cue> cue_shape(const Emission& a, const Emission& b, const Model& ma) {
  const double ba = bw_of(a), bb = bw_of(b);
  double z = std::fabs(ba - bb) / (0.08 * pymax(ba, bb) + 0.05);
  const double da = a.end_us - a.start_us, db = b.end_us - b.start_us;
  if (!ma.variable_length) z = pymax(z, std::fabs(da - db) / (0.1 * pymax(da, db) + 20.0));
  if (z <= 2) return std::nullopt;
  return Cue{"shape", -pymin(0.5 * (z * z - 4), 8.0),
             fmt("bandwidth %.2f vs %.2f MHz, duration %.0f vs %.0f us", ba, bb, da, db)};
}

std::vector<Cue> cue_identity(const Emission& a, const Emission& b) {
  const Json ia = ids_of(a), ib = ids_of(b);
  std::vector<Cue> out;
  if (ia.contains("id") && ib.contains("id")) {
    const bool same = ia["id"] == ib["id"];
    out.push_back({"decoded_id", same ? 15.0 : -15.0,
                   "decoded identity " + py_str(ia["id"]) + " vs " + py_str(ib["id"])});
  }
  for (const char* k : {"zc_root", "pilot_id", "sync_word", "bind_id"})
    if (ia.contains(k) && ib.contains(k) && ia[k] != ib[k])
      out.push_back({k, -8.0,
                     std::string(k) + " " + py_str(ia[k]) + " vs " + py_str(ib[k]) + " (radio configuration differs)"});
  return out;
}

double sigmoid(double v) { return 1 / (1 + std::exp(-pymax(pymin(v, 700), -700))); }

// ------------------------------------------------------------------ grouping
double mean_block(const dmat& M, const std::vector<int64_t>& r, const std::vector<int64_t>& c, double* mn) {
  dvec v;
  v.reserve(r.size() * c.size());
  for (const int64_t i : r)
    for (const int64_t j : c) v.push_back(M(i, j));
  if (mn) *mn = *std::min_element(v.begin(), v.end());
  return np::mean(v);
}

std::vector<std::vector<int64_t>> cluster(const std::vector<int64_t>& idx, const dmat& M, double prior_logit) {
  std::vector<std::vector<int64_t>> cl;
  for (const int64_t i : idx) cl.push_back({i});
  while (cl.size() > 1) {
    double best = 0.0;
    std::optional<std::pair<size_t, size_t>> bij;
    for (size_t i = 0; i < cl.size(); ++i)
      for (size_t j = i + 1; j < cl.size(); ++j) {
        double mn;
        const double mu = mean_block(M, cl[i], cl[j], &mn);
        if (mn <= LLR_VETO) continue;
        const double s = mu + prior_logit;
        if (s > best) {
          best = s;
          bij = {i, j};
        }
      }
    if (!bij) break;
    const auto [i, j] = *bij;
    cl[i].insert(cl[i].end(), cl[j].begin(), cl[j].end());
    std::sort(cl[i].begin(), cl[i].end());
    cl.erase(cl.begin() + static_cast<ptrdiff_t>(j));
  }
  return cl;
}

struct FrameFit {
  double t0_us = 0, frame_phase_us = 0, frame_period_ms = 0, residual_rms_ns = 0, lattice_grid_ms = 0, frame_ms = 0;
  std::optional<double> clock_ppm, clock_ppm_sigma;
  dvec residual_ns, slot_offsets_in_frame_ms;
};

std::optional<FrameFit> frame_fit(const std::vector<const Emission*>& es, const Model& ma, std::optional<double> ppm) {
  if (ma.lattices.empty() || es.size() < 2) return std::nullopt;
  const double frame = ma.lattices[0].frame_us;
  double grid = INF;
  bool any1 = false;
  for (const auto& L : ma.lattices)
    if (L.offsets_us.size() == 1) {
      grid = any1 ? std::min(grid, L.frame_us) : L.frame_us;
      any1 = true;
    }
  if (!any1) grid = frame;
  const size_t n = es.size();
  dvec t(n), s(n), f(n);
  for (size_t i = 0; i < n; ++i) {
    std::tie(t[i], s[i]) = t_of(*es[i]);
    f[i] = es[i]->center_mhz;
  }
  const double e0 = (ppm && *ppm != 0 ? *ppm : 0.0) * 1e-6;
  const double fmed = np::median(f);
  dvec u(n), w(n);
  for (size_t i = 0; i < n; ++i) {
    u[i] = std::nearbyint((t[i] - t[0]) / (grid * (1 - e0))) * grid;
    const double cdl = ma.chan_delay_us_per_mhz * (f[i] - fmed);
    w[i] = 1 / (s[i] * s[i] + ma.jitter_us * ma.jitter_us + cdl * cdl);
  }
  const double ptp = *std::max_element(u.begin(), u.end()) - *std::min_element(u.begin(), u.end());
  double p0, p1;
  FrameFit ff;
  std::optional<double> eps;
  if (n >= 3 && ptp > 0) {
    dmat AtWA(2, 2, 0.0);
    double b0 = 0, b1 = 0;
    for (size_t i = 0; i < n; ++i) {
      AtWA(0, 0) += w[i];
      AtWA(0, 1) += u[i] * w[i];
      AtWA(1, 1) += u[i] * u[i] * w[i];
      b0 += w[i] * t[i];
      b1 += u[i] * (w[i] * t[i]);
    }
    AtWA(1, 0) = AtWA(0, 1);
    const dmat cov = np::inv(AtWA);
    p0 = cov(0, 0) * b0 + cov(0, 1) * b1;
    p1 = cov(1, 0) * b0 + cov(1, 1) * b1;
    eps = -(p1 - 1) * 1e6;
    ff.clock_ppm = np::round(*eps, 3);
    ff.clock_ppm_sigma = py_round(std::sqrt(cov(1, 1)) * 1e6, 3);
  } else {
    p0 = t[0] - u[0] * (1 - e0);
    p1 = 1 - e0;
  }
  dvec res(n), r2(n);
  for (size_t i = 0; i < n; ++i) {
    res[i] = t[i] - (p0 + u[i] * p1);
    r2[i] = res[i] * res[i];
    ff.residual_ns.push_back(py_round(res[i] * 1e3, 1));
  }
  std::set<double> sub;
  for (const double v : u) sub.insert(py_round(np::mod(v, frame) / 1e3, 3));
  ff.slot_offsets_in_frame_ms.assign(sub.begin(), sub.end());
  ff.t0_us = p0;
  ff.frame_phase_us = py_round(np::mod(p0, frame), 3);
  ff.frame_period_ms = frame * (1 - (eps ? *eps : (ppm ? *ppm : 0.0)) * 1e-6) / 1e3;
  ff.residual_rms_ns = py_round(std::sqrt(np::mean(r2)) * 1e3, 1);
  ff.lattice_grid_ms = grid / 1e3;
  ff.frame_ms = frame / 1e3;
  return ff;
}

Json dvec_json(const dvec& v) {
  Json a = Json::array();
  for (const double x : v) a.push_back(x);
  return a;
}

Emitter summarise(const std::string& eid, const std::vector<int64_t>& cl, const std::vector<Emission>& all,
                  const dmat& M, const Model& ma, double prior_logit) {
  std::vector<int64_t> idx = cl;
  std::stable_sort(idx.begin(), idx.end(),
                   [&](int64_t x, int64_t y) { return t_of(all[x]).first < t_of(all[y]).first; });
  std::vector<const Emission*> es;
  for (const int64_t i : idx) es.push_back(&all[i]);
  const size_t n = es.size();
  dvec t(n);
  for (size_t i = 0; i < n; ++i) t[i] = t_of(*es[i]).first;
  std::set<double> fset;
  for (const auto* e : es) {
    const auto c = carrier_of(*e);
    const double hz = (c && c->first != 0) ? c->first : e->center_mhz * 1e6;
    fset.insert(py_round(hz / 1e6, 3));
  }
  dvec chans;
  for (const double f : fset)
    if (chans.empty() || f - chans.back() > 0.2) chans.push_back(f);
  std::vector<Ppm> pp;
  dvec pv, ps;
  for (const auto* e : es) {
    pp.push_back(ppm_of(*e, ma));
    if (pp.back().p) {
      pv.push_back(*pp.back().p);
      ps.push_back(pymax(*pp.back().s, 1e-4));
    }
  }
  Json fp = Json::object();
  std::optional<double> ppm;
  double spread = 0;
  if (!pv.empty()) {
    dvec w(pv.size()), wp(pv.size()), wd(pv.size());
    for (size_t i = 0; i < pv.size(); ++i) {
      w[i] = 1 / (ps[i] * ps[i] + ma.ppm_h1 * ma.ppm_h1);
      wp[i] = w[i] * pv[i];
    }
    const double sw = np::sum(w);
    ppm = np::sum(wp) / sw;
    for (size_t i = 0; i < pv.size(); ++i) wd[i] = w[i] * ((pv[i] - *ppm) * (pv[i] - *ppm));
    spread = py_round(std::sqrt(np::sum(wd) / sw), 4);
    Json per = Json::array();
    for (size_t i = 0; i < n; ++i)
      per.push_back(pp[i].p ? Json(py_round(*pp[i].p * es[i]->center_mhz, 0)) : Json(nullptr));
    fp["carrier_ppm"] = py_round(*ppm, 4);
    fp["carrier_ppm_sigma"] = py_round(1 / std::sqrt(sw), 4);
    fp["carrier_ppm_spread"] = spread;
    fp["carrier_ppm_source"] = pp[0].src;
    fp["carrier_offset_hz_per_hop"] = per;
  }
  const auto ff = frame_fit(es, ma, ppm);
  bool consistent = false;
  if (ff) {
    fp["timing_residual_rms_ns"] = ff->residual_rms_ns;
    fp["timing_residual_ns"] = dvec_json(ff->residual_ns);
    if (ff->clock_ppm) {
      fp["clock_ppm"] = *ff->clock_ppm;
      fp["clock_ppm_sigma"] = *ff->clock_ppm_sigma;
      if (ppm) {
        const double z = (*ff->clock_ppm - *ppm) / std::hypot(*ff->clock_ppm_sigma, 1.0);
        consistent = std::fabs(z) < 3;
        fp["clock_minus_carrier_ppm"] = np::round(*ff->clock_ppm - *ppm, 3);
        fp["clock_carrier_z"] = np::round(z, 2);
        fp["clock_carrier_consistent"] = consistent;
      }
    }
  }
  Json ts = Json::object();
  for (const auto* e : es) {
    const Json* v = jget(feats(*e), "timing_source");
    const std::string k = v ? py_str(*v) : "energy_edge";
    ts[k] = ts.value(k, 0) + 1;
  }
  Json snr = Json::array();
  dvec bws, durs;
  for (const auto* e : es) {
    snr.push_back(e->snr_db ? Json(py_round(*e->snr_db, 1)) : Json(nullptr));
    bws.push_back(bw_of(*e));
    durs.push_back(e->end_us - e->start_us);
  }
  fp["timing_source"] = ts;
  fp["snr_db"] = snr;
  fp["bw_mhz"] = py_round(np::median(bws), 3);
  fp["burst_us"] = py_round(np::median(durs), 1);
  Json seq = Json::array(), iv = Json::array();
  for (size_t i = 0; i < n; ++i) seq.push_back(Json::array({py_round(t[i] / 1e3, 4), py_round(es[i]->center_mhz, 3)}));
  for (size_t i = 1; i < n; ++i) iv.push_back(py_round((t[i] - t[i - 1]) / 1e3, 4));
  Json hop{{"n_bursts", n}, {"n_channels", chans.size()}, {"sequence", seq}, {"intervals_ms", iv}};
  if (ff) {
    hop["frame_ms"] = ff->frame_ms;
    hop["slot_offsets_in_frame_ms"] = dvec_json(ff->slot_offsets_in_frame_ms);
    hop["frame_phase_us"] = ff->frame_phase_us;
  }
  Emitter em;
  em.emitter_id = eid;
  em.family = es[0]->family;
  em.label = es[0]->label;
  em.members = idx;
  em.frequencies = chans;
  em.hop_pattern = hop;
  if (ff) em.frame_period_ms = np::round(ff->frame_period_ms, 7);
  if (n >= 2) {
    // M[np.ix_(cl, cl)] rows in cluster (index) order
    dvec coh;
    dvec upper;
    for (size_t r = 0; r < cl.size(); ++r) {
      dvec row;
      for (size_t c = 0; c < cl.size(); ++c) {
        if (c != r) row.push_back(M(cl[r], cl[c]));
        if (c > r) upper.push_back(M(cl[r], cl[c]));
      }
      coh.push_back(np::mean(row));
    }
    const double cmin = *std::min_element(coh.begin(), coh.end());
    em.confidence = py_round(sigmoid(cmin + prior_logit), 4);
    em.evidence.push_back(fmt("cohesion: weakest member's mean LLR to the rest %+.1f nats (mean pairwise %+.1f)", cmin,
                              np::mean(upper)));
    if (ff) {
      double W = INF;
      for (const auto& L : ma.lattices) W = pymin(W, L.frame_us / static_cast<double>(L.offsets_us.size()));
      const double tptp = *std::max_element(t.begin(), t.end()) - *std::min_element(t.begin(), t.end());
      em.evidence.push_back(fmt("timing: %zu burst starts on ONE frame clock (%g ms frame, slots at %s ms, phase %.4f "
                                "ms): residual rms %.0f ns over %.1f ms; an independent radio's start is uniform over "
                                "~%g ms",
                                n, ff->frame_ms, py_repr(dvec_json(ff->slot_offsets_in_frame_ms)).c_str(),
                                ff->frame_phase_us / 1e3, ff->residual_rms_ns, tptp / 1e3, W / 1e3));
    }
    if (pv.size() >= 2)
      em.evidence.push_back(fmt("oscillator: carrier offset %+.3f ppm on all %zu hops (weighted spread %.3f ppm) "
                                "across %.1f-%.1f MHz: one crystal; two radios differ by ~%.1f ppm rms",
                                *ppm, pv.size(), spread, chans.front(), chans.back(), std::sqrt(2.0) * ma.ppm_pop));
    if (fp.contains("clock_ppm") && ppm)
      em.evidence.push_back(fmt("clock/carrier: frame clock %+.2f +- %.2f ppm vs carrier %+.3f ppm (%s with one "
                                "reference crystal)",
                                *ff->clock_ppm, *ff->clock_ppm_sigma, *ppm, consistent ? "consistent" : "INCONSISTENT"));
  } else {
    em.evidence.push_back("single emission: identity rests on the separation evidence");
  }
  em.fingerprint = fp;
  return em;
}

}  // namespace

// =================================================================== public API
const Models& MODELS() {
  static const Models m = [] {
    Model up;
    up.lattices = {Lattice{10000.0, {0.0, 4000.0}, 0.6}, Lattice{500.0, {0.0}, 0.4}};
    up.jitter_us = 0.03;
    up.chan_delay_us_per_mhz = 0.02;
    up.raster = std::make_pair(0.0, 0.5e6);
    up.ppm_h1 = 0.06;
    up.ppm_drift_per_s = 0.05;
    up.ppm_pop = 2.0;
    up.snr_h1_db = 4.0;
    up.single_tx = true;
    Model did = up;
    did.lattices = {};
    did.jitter_us = 1.0;
    Model video = did;
    video.variable_length = true;
    Models out;
    out[{"dji_ocusync", "ocusync_uplink"}] = up;
    out[{"dji_droneid", "droneid"}] = did;
    out[{"dji_ocusync", "ocusync_video"}] = video;
    for (const char* lab : {"video_10MHz", "video_20MHz", "video_o4_10MHz", "video_o4_20MHz"})
      out[{"dji_ocusync", lab}] = video;
    out[{"dji_ocusync", "uplink"}] = up;
    return out;
  }();
  return m;
}

Json Cue::to_json() const {
  Json j{{"cue", cue}, {"llr", llr}, {"detail", detail}};
  for (const auto& [k, v] : extra.items()) j[k] = v;
  return j;
}

Json Emitter::to_dict() const {
  Json sep = separation, lk = linked;
  return Json{{"emitter_id", emitter_id},
              {"family", family},
              {"label", label},
              {"members", members},
              {"frequencies", dvec_json(frequencies)},
              {"hop_pattern", hop_pattern},
              {"frame_period_ms", frame_period_ms ? Json(*frame_period_ms) : Json(nullptr)},
              {"fingerprint", fingerprint},
              {"confidence", confidence ? Json(*confidence) : Json(nullptr)},
              {"evidence", evidence},
              {"separation", sep},
              {"linked", lk}};
}

Json measure_fingerprints(const cfvec& x, double fs, double fc, std::vector<Emission>& emissions,
                          const SpectralChannelizer* ch, const std::vector<std::pair<double, double>>& notch_hz) {
  std::unique_ptr<SpectralChannelizer> own;
  if (!ch) {   // notch CW spurs / LO leakage: a tone correlates at every lag
    std::vector<std::pair<double, double>> notch = notch_hz;
    if (notch.empty())
      for (const auto& sp : find_spurs(x, fs, fc)) notch.emplace_back(sp.offset_hz, 30e3);
    cdvec xd(x.size());
    for (size_t i = 0; i < x.size(); ++i) xd[i] = cd(x[i].real(), x[i].imag());
    own = std::make_unique<SpectralChannelizer>(xd, fs, notch);
    ch = own.get();
  }
  std::vector<std::pair<std::tuple<std::string, std::string, int64_t>, std::vector<size_t>>> groups;
  for (size_t i = 0; i < emissions.size(); ++i) {
    const Emission& e = emissions[i];
    if (!is_ofdm15(e)) continue;
    const auto key = std::make_tuple(e.family, e.label, py_round_int((e.f_hi_mhz - e.f_lo_mhz) / 2));
    auto it = std::find_if(groups.begin(), groups.end(), [&](const auto& g) { return g.first == key; });
    if (it == groups.end()) groups.push_back({key, {i}});
    else it->second.push_back(i);
  }
  Json info = Json::object();
  for (const auto& [key, idx] : groups) {
    dvec bws;
    for (const size_t i : idx) bws.push_back(emissions[i].f_hi_mhz - emissions[i].f_lo_mhz);
    const double bw = np::median(bws) * 1e6;
    const double half_bw = pymin(bw / 2 * 1.08, R15 / 2 * 0.95);
    std::vector<CpItem> items;
    std::vector<size_t> ids;
    for (const size_t i : idx)
      if (auto it = cp_stage(*ch, fc, emissions[i], half_bw)) {
        items.push_back(*it);
        ids.push_back(i);
      }
    const Json tinfo = items.size() >= 2 ? template_stage(*ch, items, half_bw) : Json::object();
    info[fmt("('%s', '%s', %lld)", std::get<0>(key).c_str(), std::get<1>(key).c_str(),
             static_cast<long long>(std::get<2>(key)))] = tinfo;
    for (size_t n = 0; n < ids.size(); ++n) {
      Emission& e = emissions[ids[n]];
      const CpItem& it = items[n];
      if (!e.features.is_object()) e.features = Json::object();
      Json& f = e.features;
      const bool fine = it.ref_ok;
      f["t_start_us"] = (fine ? it.t_fine : it.t_cp) * 1e6;
      f["t_sigma_us"] = (fine ? it.sig_t_fine : it.sig_t) * 1e6;
      f["carrier_hz"] = fine ? it.carrier_fine : it.carrier;
      f["carrier_sigma_hz"] = fine ? it.sig_nu_fine : it.sig_nu;
      f["t_cp_us"] = it.t_cp * 1e6;
      f["carrier_cp_hz"] = it.carrier;
      f["timing_source"] = fine ? "reference_symbols" : "cp_train";
      f["ref_corr"] = it.ref_corr ? Json(*it.ref_corr) : Json(nullptr);
      f["n_occ_sc"] = it.n_occ;
      f["cp_peaks"] = it.cp_peaks;
      if (!e.snr_db) e.snr_db = py_round(it.snr_db, 1);
      f["snr_meas_db"] = py_round(it.snr_db, 1);
    }
  }
  return info;
}

std::pair<double, std::vector<Cue>> same_emitter_score(const Emission& a, const Emission& b, const Models* models,
                                                       double prior_same) {
  std::vector<Cue> ev;
  double llr;
  if (a.family != b.family || a.label != b.label) {
    ev.push_back({"type", -30.0,
                  a.family + "/" + a.label + " vs " + b.family + "/" + b.label + ": different transmitter types"});
    llr = -30.0;
  } else {
    const Model ma = model_for(a, models);
    for (auto& c : cue_identity(a, b)) ev.push_back(std::move(c));
    if (auto c = cue_overlap(a, b, ma)) ev.push_back(*c);
    if (auto c = cue_oscillator(a, b, ma)) ev.push_back(*c);
    const auto pa = ppm_of(a, ma).p, pb = ppm_of(b, ma).p;
    std::optional<double> ppm_pair;
    if (pa && pb) ppm_pair = 0.5 * (*pa + *pb);
    for (auto c : {cue_timing(a, b, ma, ppm_pair), cue_clock(a, b, ma), cue_power(a, b, ma), cue_shape(a, b, ma)})
      if (c) ev.push_back(*c);
    llr = 0;
    for (const auto& c : ev) llr += c.llr;
  }
  const double lp = llr + std::log(prior_same / (1 - prior_same));
  const double post = 1 / (1 + std::exp(-pymax(pymin(lp, 700), -700)));
  Cue tot{"total", llr, fmt("sum LLR %+.1f nats -> P(same | prior %g) = %.4g", llr, prior_same, post)};
  tot.extra = Json{{"p_same", post}};
  ev.push_back(tot);
  for (auto& c : ev) c.llr = py_round(c.llr, 2);
  return {llr, ev};
}

Association associate_full(const std::vector<Emission>& emissions, const Models* models, double prior_same) {
  const size_t n = emissions.size();
  const double prior_logit = std::log(prior_same / (1 - prior_same));
  Association A;
  A.M = dmat(n, n, -30.0);
  for (size_t i = 0; i < n; ++i) {
    A.M(i, i) = 0.0;
    for (size_t j = i + 1; j < n; ++j) {
      const Emission &a = emissions[i], &b = emissions[j];
      if (a.family != b.family || a.label != b.label) continue;
      auto [llr, ev] = same_emitter_score(a, b, models, prior_same);
      A.M(i, j) = A.M(j, i) = llr;
      A.EV[{static_cast<int64_t>(i), static_cast<int64_t>(j)}] = std::move(ev);
    }
  }
  std::vector<std::pair<ModelKey, std::vector<int64_t>>> classes;
  for (size_t i = 0; i < n; ++i) {
    const ModelKey key{emissions[i].family, emissions[i].label};
    auto it = std::find_if(classes.begin(), classes.end(), [&](const auto& c) { return c.first == key; });
    if (it == classes.end()) classes.push_back({key, {static_cast<int64_t>(i)}});
    else it->second.push_back(static_cast<int64_t>(i));
  }
  auto& out = A.emitters;
  for (const auto& [key, idx] : classes) {
    const Model ma = model_for(emissions[idx[0]], models);
    auto cls = cluster(idx, A.M, prior_logit);
    auto tmin = [&](const std::vector<int64_t>& c) {
      double m = INF;
      for (const int64_t i : c) m = pymin(m, t_of(emissions[i]).first);
      return m;
    };
    std::stable_sort(cls.begin(), cls.end(), [&](const auto& x, const auto& y) { return tmin(x) < tmin(y); });
    for (const auto& cl : cls)
      out.push_back(summarise(fmt("E%zu", out.size() + 1), cl, emissions, A.M, ma, prior_logit));
  }
  // separation between emitters of the same type
  for (auto& a : out) {
    for (const auto& b : out) {
      if (&a == &b || a.family != b.family || a.label != b.label) continue;
      const double mean = mean_block(A.M, a.members, b.members, nullptr);
      std::vector<std::pair<std::string, dvec>> cues;
      for (const int64_t i : a.members)
        for (const int64_t j : b.members) {
          const auto it = A.EV.find({std::min(i, j), std::max(i, j)});
          if (it == A.EV.end()) continue;
          for (const auto& c : it->second) {
            if (c.cue == "total") continue;
            auto ci = std::find_if(cues.begin(), cues.end(), [&](const auto& p) { return p.first == c.cue; });
            if (ci == cues.end()) cues.push_back({c.cue, {c.llr}});
            else ci->second.push_back(c.llr);
          }
        }
      Json cm = Json::object();
      for (const auto& [k, v] : cues) cm[k] = py_round(np::mean(v), 2);
      a.separation.push_back(Json{{"other", b.emitter_id},
                                  {"mean_llr", py_round(mean, 2)},
                                  {"p_same", py_round(sigmoid(mean + prior_logit), 6)},
                                  {"verdict", mean + prior_logit < 0 ? "DIFFERENT radio" : "ambiguous"},
                                  {"cue_mean_llr", cm}});
    }
    if (!a.confidence && !a.separation.empty()) {   // singleton: confidence it is distinct
      double mx = -INF;
      for (const auto& s : a.separation) mx = pymax(mx, s["p_same"].get<double>());
      a.confidence = py_round(1 - mx, 4);
    }
  }
  // same drone system: different-type emitters sharing the oscillator
  auto prefix = [](const std::string& s) { return s.substr(0, s.find('_')); };
  for (auto& a : out)
    for (const auto& b : out) {
      if (&a == &b || prefix(a.family) != prefix(b.family) || a.label == b.label) continue;
      if (!a.fingerprint.contains("carrier_ppm") || !b.fingerprint.contains("carrier_ppm")) continue;
      const double pa = a.fingerprint["carrier_ppm"].get<double>(), pb = b.fingerprint["carrier_ppm"].get<double>();
      const double sa = std::hypot(a.fingerprint["carrier_ppm_sigma"].get<double>(), 0.1);
      const double sb = std::hypot(b.fingerprint["carrier_ppm_sigma"].get<double>(), 0.1);
      const double z = (pa - pb) / std::hypot(sa, sb);
      if (std::fabs(z) < 3)
        a.linked.push_back(Json{{"other", b.emitter_id},
                                {"label", b.label},
                                {"dppm", py_round(pa - pb, 4)},
                                {"note", "same reference oscillator within error: likely the same drone system (RC "
                                         "uplink locks to the drone)"}});
    }
  return A;
}

}  // namespace si::rx::emitters
