#include "core/signal.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

#include "core/fft.hpp"

namespace si::sig {

using np::cabs;
using np::cdiv;
using np::cmul;
using np::is_cplx;
using np::real_t;

// ------------------------------------------------------------------ windows
double i0(double x) {
  // Cephes i0 (scipy.special.i0 / numpy's np.i0 coefficient tables)
  static const double A[] = {
      -4.41534164647933937950E-18, 3.33079451882223809783E-17, -2.43127984654795469359E-16,
      1.71539128555513303061E-15,  -1.16853328779934516808E-14, 7.67618549860493561688E-14,
      -4.85644678311192946090E-13, 2.95505266312963983461E-12, -1.72682629144155570723E-11,
      9.67580903537323691224E-11,  -5.18979560163526290666E-10, 2.65982372468238665035E-9,
      -1.30002500998624804212E-8,  6.04699502254191894932E-8,  -2.67079385394061173391E-7,
      1.11738753912010371815E-6,   -4.41673835845875056359E-6, 1.64484480707288970893E-5,
      -5.75419501008210370398E-5,  1.88502885095841655729E-4,  -5.76375574538582365885E-4,
      1.63947561694133579842E-3,   -4.32430999505057594430E-3, 1.05464603945949983183E-2,
      -2.37374148058994688156E-2,  4.93052842396707084878E-2,  -9.49010970480476444210E-2,
      1.71620901522208775349E-1,   -3.04682672343198398683E-1, 6.76795274409476084995E-1};
  static const double B[] = {
      -7.23318048787475395456E-18, -4.83050448594418207126E-18, 4.46562142029675999901E-17,
      3.46122286769746109310E-17,  -2.82762398051658348494E-16, -3.42548561967721913462E-16,
      1.77256013305652638360E-15,  3.81168066935262242075E-15,  -9.55484669882830764870E-15,
      -4.15056934728722208663E-14, 1.54008621752140982691E-14,  3.85277838274214270114E-13,
      7.18012445138366623367E-13,  -1.79417853150680611778E-12, -1.32158118404477131188E-11,
      -3.14991652796324136454E-11, 1.18891471078464383424E-11,  4.94060238822496958910E-10,
      3.39623202570838634515E-9,   2.26666899049817806459E-8,   2.04891858946906374183E-7,
      2.89137052083475648297E-6,   6.88975834691682398426E-5,   3.36911647825569408990E-3,
      8.04490411014108831608E-1};
  auto chbevl = [](double v, const double* p, int n) {
    double b0 = *p++, b1 = 0.0, b2 = 0.0;
    int i = n - 1;
    do {
      b2 = b1;
      b1 = b0;
      b0 = v * b1 - b2 + *p++;
    } while (--i);
    return 0.5 * (b0 - b2);
  };
  if (x < 0) x = -x;
  if (x <= 8.0) return std::exp(x) * chbevl(x / 2.0 - 2.0, A, 30);
  return std::exp(x) * chbevl(32.0 / x - 2.0, B, 25) / std::sqrt(x);
}

namespace {
// scipy.signal.windows.general_cosine(M, a, sym)
dvec general_cosine(size_t M, const dvec& a, bool sym) {
  if (M <= 1) return dvec(M, 1.0);
  const size_t Mx = sym ? M : M + 1;
  const dvec fac = np::linspace<double>(-kPi, kPi, Mx, true);
  dvec w(Mx, 0.0);
  for (size_t k = 0; k < a.size(); ++k)
    for (size_t i = 0; i < Mx; ++i) w[i] += a[k] * std::cos(static_cast<double>(k) * fac[i]);
  w.resize(M);
  return w;
}
}  // namespace

dvec get_window(const Window& win, size_t N, bool fftbins) {
  const bool sym = !fftbins;
  switch (win.kind) {
    case Window::Hann: return general_cosine(N, {0.5, 1.0 - 0.5}, sym);
    case Window::Hamming: return general_cosine(N, {0.54, 1.0 - 0.54}, sym);
    case Window::Boxcar: return dvec(N, 1.0);
    case Window::Kaiser: {
      if (N <= 1) return dvec(N, 1.0);
      const size_t M = sym ? N : N + 1;
      const double alpha = static_cast<double>(M - 1) / 2.0;
      const double den = i0(win.beta);
      dvec w(M);
      for (size_t n = 0; n < M; ++n) {
        const double r = (static_cast<double>(n) - alpha) / alpha;
        w[n] = i0(win.beta * std::sqrt(1.0 - r * r)) / den;
      }
      w.resize(N);
      return w;
    }
  }
  return {};
}

// ------------------------------------------------------------------ find_peaks
namespace {
struct LocalMax {
  std::vector<int64_t> mid, left, right;
};
LocalMax local_maxima_1d(const dvec& x) {
  LocalMax r;
  const int64_t n = static_cast<int64_t>(x.size());
  int64_t i = 1, i_max = n - 1;
  while (i < i_max) {
    if (x[i - 1] < x[i]) {
      int64_t ahead = i + 1;
      while (ahead < i_max && x[ahead] == x[i]) ++ahead;
      if (x[ahead] < x[i]) {
        r.left.push_back(i);
        r.right.push_back(ahead - 1);
        r.mid.push_back((i + ahead - 1) / 2);
        i = ahead;
      }
    }
    ++i;
  }
  return r;
}

template <class V>
void keep_by(std::vector<V>& v, const std::vector<uint8_t>& keep) {
  if (v.empty()) return;
  size_t k = 0;
  for (size_t i = 0; i < v.size(); ++i)
    if (keep[i]) v[k++] = v[i];
  v.resize(k);
}

void apply_keep(Peaks& P, const std::vector<uint8_t>& keep) {
  keep_by(P.peaks, keep);
  keep_by(P.peak_heights, keep);
  keep_by(P.left_thresholds, keep);
  keep_by(P.right_thresholds, keep);
  keep_by(P.prominences, keep);
  keep_by(P.left_bases, keep);
  keep_by(P.right_bases, keep);
  keep_by(P.widths, keep);
  keep_by(P.width_heights, keep);
  keep_by(P.left_ips, keep);
  keep_by(P.right_ips, keep);
}

std::vector<uint8_t> select_by_property(const dvec& v, const std::optional<double>& lo, const std::optional<double>& hi) {
  std::vector<uint8_t> keep(v.size(), 1);
  for (size_t i = 0; i < v.size(); ++i) {
    if (lo && !(*lo <= v[i])) keep[i] = 0;
    if (hi && !(v[i] <= *hi)) keep[i] = 0;
  }
  return keep;
}
}  // namespace

Peaks find_peaks(const dvec& x, const PeakOpts& opt) {
  Peaks P;
  if (opt.distance && *opt.distance < 1) throw std::invalid_argument("`distance` must be greater or equal to 1");
  P.peaks = local_maxima_1d(x).mid;
  if (opt.height_min || opt.height_max) {
    P.peak_heights.resize(P.peaks.size());
    for (size_t i = 0; i < P.peaks.size(); ++i) P.peak_heights[i] = x[P.peaks[i]];
    apply_keep(P, select_by_property(P.peak_heights, opt.height_min, opt.height_max));
  }
  if (opt.threshold_min || opt.threshold_max) {
    const size_t n = P.peaks.size();
    P.left_thresholds.resize(n);
    P.right_thresholds.resize(n);
    std::vector<uint8_t> keep(n, 1);
    for (size_t i = 0; i < n; ++i) {
      const int64_t p = P.peaks[i];
      const double l = x[p] - x[p - 1], r = x[p] - x[p + 1];
      P.left_thresholds[i] = l;
      P.right_thresholds[i] = r;
      if (opt.threshold_min && !(*opt.threshold_min <= std::min(l, r))) keep[i] = 0;
      if (opt.threshold_max && !(std::max(l, r) <= *opt.threshold_max)) keep[i] = 0;
    }
    apply_keep(P, keep);
  }
  if (opt.distance) {
    const size_t n = P.peaks.size();
    const double dist = std::ceil(*opt.distance);
    std::vector<uint8_t> keep(n, 1);
    dvec prio(n);
    for (size_t i = 0; i < n; ++i) prio[i] = x[P.peaks[i]];
    const std::vector<int64_t> order = np::argsort_quick(prio);
    for (int64_t i = static_cast<int64_t>(n) - 1; i >= 0; --i) {
      const int64_t j = order[i];
      if (!keep[j]) continue;
      int64_t k = j - 1;
      while (k >= 0 && static_cast<double>(P.peaks[j] - P.peaks[k]) < dist) { keep[k] = 0; --k; }
      k = j + 1;
      while (k < static_cast<int64_t>(n) && static_cast<double>(P.peaks[k] - P.peaks[j]) < dist) { keep[k] = 0; ++k; }
    }
    apply_keep(P, keep);
  }
  const bool want_prom = opt.prominence || opt.prominence_min || opt.prominence_max;
  const bool want_width = opt.width || opt.width_min || opt.width_max;
  if (want_prom || want_width) {
    int64_t wlen = -1;
    if (opt.wlen > 1) wlen = static_cast<int64_t>(std::ceil(opt.wlen));
    else if (opt.wlen >= 0) throw std::invalid_argument("`wlen` must be larger than 1");
    const size_t n = P.peaks.size();
    P.prominences.resize(n);
    P.left_bases.resize(n);
    P.right_bases.resize(n);
    const int64_t nx = static_cast<int64_t>(x.size());
    for (size_t pn = 0; pn < n; ++pn) {
      const int64_t peak = P.peaks[pn];
      int64_t i_min = 0, i_max = nx - 1;
      if (2 <= wlen) {
        i_min = std::max(peak - wlen / 2, i_min);
        i_max = std::min(peak + wlen / 2, i_max);
      }
      int64_t i = P.left_bases[pn] = peak;
      double left_min = x[peak];
      while (i_min <= i && x[i] <= x[peak]) {
        if (x[i] < left_min) { left_min = x[i]; P.left_bases[pn] = i; }
        --i;
      }
      i = P.right_bases[pn] = peak;
      double right_min = x[peak];
      while (i <= i_max && x[i] <= x[peak]) {
        if (x[i] < right_min) { right_min = x[i]; P.right_bases[pn] = i; }
        ++i;
      }
      P.prominences[pn] = x[peak] - std::max(left_min, right_min);
    }
  }
  if (want_prom) apply_keep(P, select_by_property(P.prominences, opt.prominence_min, opt.prominence_max));
  if (want_width) {
    const size_t n = P.peaks.size();
    P.widths.resize(n);
    P.width_heights.resize(n);
    P.left_ips.resize(n);
    P.right_ips.resize(n);
    for (size_t p = 0; p < n; ++p) {
      const int64_t i_min = P.left_bases[p], i_max = P.right_bases[p], peak = P.peaks[p];
      const double height = P.width_heights[p] = x[peak] - P.prominences[p] * opt.rel_height;
      int64_t i = peak;
      while (i_min < i && height < x[i]) --i;
      double left_ip = static_cast<double>(i);
      if (x[i] < height) left_ip += (height - x[i]) / (x[i + 1] - x[i]);
      i = peak;
      while (i < i_max && height < x[i]) ++i;
      double right_ip = static_cast<double>(i);
      if (x[i] < height) right_ip -= (height - x[i]) / (x[i - 1] - x[i]);
      P.widths[p] = right_ip - left_ip;
      P.left_ips[p] = left_ip;
      P.right_ips[p] = right_ip;
    }
    apply_keep(P, select_by_property(P.widths, opt.width_min, opt.width_max));
  }
  return P;
}

Peaks find_peaks(const fvec& x, const PeakOpts& opt) { return find_peaks(dvec(x.begin(), x.end()), opt); }

// ------------------------------------------------------------------ firwin
dvec firwin(size_t numtaps, const dvec& cutoff_in, std::optional<double> fs, const Window& window, bool pass_zero,
            bool scale) {
  const double nyq = 0.5 * (fs ? *fs : 2.0);
  if (cutoff_in.empty()) throw std::invalid_argument("At least one cutoff frequency must be given.");
  dvec cutoff(cutoff_in.size());
  for (size_t i = 0; i < cutoff.size(); ++i) cutoff[i] = cutoff_in[i] / nyq;
  for (size_t i = 0; i < cutoff.size(); ++i) {
    if (cutoff[i] <= 0 || cutoff[i] >= 1)
      throw std::invalid_argument("Invalid cutoff frequency: frequencies must be greater than 0 and less than fs/2.");
    if (i && !(cutoff[i] > cutoff[i - 1])) throw std::invalid_argument("Invalid cutoff frequencies: the frequencies must be strictly increasing.");
  }
  const bool pass_nyquist = static_cast<bool>(cutoff.size() & 1) ^ pass_zero;
  if (pass_nyquist && numtaps % 2 == 0)
    throw std::invalid_argument("A filter with an even number of coefficients must have zero response at the Nyquist frequency.");
  dvec edges;
  if (pass_zero) edges.push_back(0.0);
  edges.insert(edges.end(), cutoff.begin(), cutoff.end());
  if (pass_nyquist) edges.push_back(1.0);
  const double alpha = 0.5 * static_cast<double>(numtaps - 1);
  dvec m(numtaps), h(numtaps, 0.0);
  for (size_t i = 0; i < numtaps; ++i) m[i] = static_cast<double>(i) - alpha;
  for (size_t b = 0; b + 1 < edges.size(); b += 2) {
    const double left = edges[b], right = edges[b + 1];
    for (size_t i = 0; i < numtaps; ++i) h[i] += right * np::sinc(right * m[i]);
    for (size_t i = 0; i < numtaps; ++i) h[i] -= left * np::sinc(left * m[i]);
  }
  const dvec win = get_window(window, numtaps, false);
  for (size_t i = 0; i < numtaps; ++i) h[i] *= win[i];
  if (scale) {
    const double left = edges[0], right = edges[1];
    double sf;
    if (left == 0) sf = 0.0;
    else if (right == 1) sf = 1.0;
    else sf = 0.5 * (left + right);
    dvec hc(numtaps);
    for (size_t i = 0; i < numtaps; ++i) hc[i] = h[i] * std::cos(kPi * m[i] * sf);
    const double s = np::pairwise_sum(hc.data(), numtaps);
    for (auto& v : h) v /= s;
  }
  return h;
}

// ------------------------------------------------------------------ butter (zpk -> sos)
namespace {
struct ZPK {
  cdvec z, p;
  double k = 1.0;
};

cd cprod(const cdvec& v) {
  cd acc(1.0, 0.0);
  for (const cd& x : v) acc = cmul(acc, x);
  return acc;
}

ZPK buttap(int N) {
  ZPK s;
  for (int m = -N + 1; m < N; m += 2) {
    // -exp(1j*pi*m/(2*N)): (1j*pi)*m -> (0, pi*m) ; / (2N) -> Smith division
    const cd num = cmul(cd(0.0, kPi), cd(static_cast<double>(m), 0.0));
    const cd arg = cdiv(num, cd(2.0 * N, 0.0));
    const double e = std::exp(arg.real());
    s.p.push_back(-cd(e * std::cos(arg.imag()), e * std::sin(arg.imag())));
  }
  s.k = 1.0;
  return s;
}

ZPK lp2lp(const ZPK& in, double wo) {
  ZPK o;
  const int degree = static_cast<int>(in.p.size()) - static_cast<int>(in.z.size());
  for (const cd& z : in.z) o.z.push_back(cmul(cd(wo, 0), z));
  for (const cd& p : in.p) o.p.push_back(cmul(cd(wo, 0), p));
  o.k = in.k * std::pow(wo, degree);
  return o;
}

ZPK lp2hp(const ZPK& in, double wo) {
  ZPK o;
  const int degree = static_cast<int>(in.p.size()) - static_cast<int>(in.z.size());
  for (const cd& z : in.z) o.z.push_back(cdiv(cd(wo, 0), z));
  for (const cd& p : in.p) o.p.push_back(cdiv(cd(wo, 0), p));
  for (int i = 0; i < degree; ++i) o.z.emplace_back(0.0, 0.0);
  cdvec nz, np_;
  for (const cd& z : in.z) nz.push_back(-z);
  for (const cd& p : in.p) np_.push_back(-p);
  o.k = in.k * cdiv(cprod(nz), cprod(np_)).real();
  return o;
}

ZPK lp2bp(const ZPK& in, double wo, double bw) {
  ZPK o;
  const int degree = static_cast<int>(in.p.size()) - static_cast<int>(in.z.size());
  auto shift = [&](const cdvec& v, cdvec& out) {
    cdvec lp, s;
    for (const cd& x : v) lp.push_back(cdiv(cmul(x, cd(bw, 0)), cd(2.0, 0)));
    for (const cd& x : lp) {
      const cd sq = cmul(x, x);
      s.push_back(np::csqrt(cd(sq.real() - wo * wo, sq.imag())));
    }
    for (size_t i = 0; i < lp.size(); ++i) out.push_back(lp[i] + s[i]);
    for (size_t i = 0; i < lp.size(); ++i) out.push_back(lp[i] - s[i]);
  };
  shift(in.z, o.z);
  shift(in.p, o.p);
  for (int i = 0; i < degree; ++i) o.z.emplace_back(0.0, 0.0);
  o.k = in.k * std::pow(bw, degree);
  return o;
}

ZPK bilinear(const ZPK& in, double fs) {
  ZPK o;
  const int degree = static_cast<int>(in.p.size()) - static_cast<int>(in.z.size());
  const double fs2 = 2.0 * fs;
  cdvec nz, np_;
  const cd F(fs2, 0.0);
  for (const cd& z : in.z) {
    o.z.push_back(cdiv(F + z, F - z));
    nz.push_back(F - z);
  }
  for (const cd& p : in.p) {
    o.p.push_back(cdiv(F + p, F - p));
    np_.push_back(F - p);
  }
  for (int i = 0; i < degree; ++i) o.z.emplace_back(-1.0, 0.0);
  const cd num = in.z.empty() ? cd(1.0, 0.0) : cprod(nz);
  o.k = in.k * cdiv(num, cprod(np_)).real();
  return o;
}

inline bool isreal(const cd& v) { return v.imag() == 0.0; }

// scipy.signal._filter_design._cplxreal: (complex representatives with imag > 0, reals)
void cplxreal(cdvec z, cdvec& zc, cdvec& zr) {
  zc.clear();
  zr.clear();
  if (z.empty()) return;
  const double tol = 100 * std::numeric_limits<double>::epsilon();
  // lexsort((abs(z.imag), z.real)): by real, then |imag| (stable)
  std::vector<size_t> ord(z.size());
  std::iota(ord.begin(), ord.end(), 0);
  std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) {
    if (z[a].real() != z[b].real()) return z[a].real() < z[b].real();
    return std::fabs(z[a].imag()) < std::fabs(z[b].imag());
  });
  cdvec zs;
  for (size_t i : ord) zs.push_back(z[i]);
  cdvec rest;
  for (const cd& v : zs) {
    if (std::fabs(v.imag()) <= tol * cabs(v)) zr.push_back(cd(v.real(), 0.0));
    else rest.push_back(v);
  }
  if (rest.empty()) return;
  cdvec zp, zn;
  for (const cd& v : rest) {
    if (v.imag() > 0) zp.push_back(v);
    else if (v.imag() < 0) zn.push_back(v);
  }
  if (zp.size() != zn.size()) throw std::runtime_error("Array contains complex value with no matching conjugate.");
  // runs of (approximately) equal real part, each sorted by |imag|
  std::vector<int> same(zp.size() > 0 ? zp.size() - 1 : 0);
  for (size_t i = 0; i + 1 < zp.size(); ++i) same[i] = (zp[i + 1].real() - zp[i].real()) <= tol * cabs(zp[i]);
  std::vector<int> padded;
  padded.push_back(0);
  padded.insert(padded.end(), same.begin(), same.end());
  padded.push_back(0);
  std::vector<size_t> starts, stops;
  for (size_t i = 0; i + 1 < padded.size(); ++i) {
    const int d = padded[i + 1] - padded[i];
    if (d > 0) starts.push_back(i);
    if (d < 0) stops.push_back(i);
  }
  for (size_t r = 0; r < starts.size(); ++r) {
    const size_t a = starts[r], b = stops[r] + 1;
    for (cdvec* chunk : {&zp, &zn}) {
      std::stable_sort(chunk->begin() + a, chunk->begin() + b,
                       [](const cd& u, const cd& v) { return std::fabs(u.imag()) < std::fabs(v.imag()); });
    }
  }
  for (size_t i = 0; i < zp.size(); ++i) {
    const cd c = std::conj(zn[i]);
    if (cabs(cd(zp[i].real() - c.real(), zp[i].imag() - c.imag())) > tol * cabs(zn[i]))
      throw std::runtime_error("Array contains complex value with no matching conjugate.");
    const cd s(zp[i].real() + c.real(), zp[i].imag() + c.imag());
    zc.push_back(cdiv(s, cd(2.0, 0.0)));
  }
}

size_t idx_worst(const cdvec& p) {
  dvec d(p.size());
  for (size_t i = 0; i < p.size(); ++i) d[i] = std::fabs(1.0 - cabs(p[i]));
  return np::argmin(d);
}

size_t nearest_idx(const cdvec& fro, const cd& to, int which) {  // 0 any, 1 real, 2 complex
  dvec d(fro.size());
  for (size_t i = 0; i < fro.size(); ++i) d[i] = cabs(cd(fro[i].real() - to.real(), fro[i].imag() - to.imag()));
  const std::vector<int64_t> order = np::argsort_quick(d);
  if (which == 0) return static_cast<size_t>(order[0]);
  for (int64_t o : order) {
    const bool r = isreal(fro[o]);
    if ((which == 1 && r) || (which == 2 && !r)) return static_cast<size_t>(o);
  }
  throw std::runtime_error("zpk2sos: no matching zero");
}

// np.poly(roots) (real result: the receivers' roots come in conjugate pairs)
dvec poly_real(const cdvec& roots) {
  cdvec a{cd(1.0, 0.0)};
  for (const cd& z : roots) {
    cdvec nxt(a.size() + 1, cd(0.0, 0.0));
    const cd mz = -z;
    for (size_t i = 0; i < nxt.size(); ++i) {
      cd s(0.0, 0.0);
      if (i < a.size()) s += a[i];                  // a[i] * 1
      if (i >= 1) s += cmul(a[i - 1], mz);          // a[i-1] * (-z)
      nxt[i] = s;
    }
    a.swap(nxt);
  }
  dvec r(a.size());
  for (size_t i = 0; i < a.size(); ++i) r[i] = a[i].real();
  return r;
}

void single_zpksos(const cdvec& z, const cdvec& p, double* row) {
  for (int i = 0; i < 6; ++i) row[i] = 0.0;
  const dvec b = poly_real(z), a = poly_real(p);
  for (size_t i = 0; i < b.size(); ++i) row[3 - b.size() + i] = b[i];
  for (size_t i = 0; i < a.size(); ++i) row[6 - a.size() + i] = a[i];
}

dmat zpk2sos(ZPK s) {
  cdvec z = s.z, p = s.p;
  if (z.empty() && p.empty()) {
    dmat out(1, 6, 0.0);
    out(0, 0) = s.k; out(0, 3) = 1.0;
    return out;
  }
  while (p.size() < z.size()) p.emplace_back(0.0, 0.0);
  while (z.size() < p.size()) z.emplace_back(0.0, 0.0);
  const size_t n_sections = (std::max(p.size(), z.size()) + 1) / 2;
  if (p.size() % 2 == 1) {
    p.emplace_back(0.0, 0.0);
    z.emplace_back(0.0, 0.0);
  }
  cdvec zc, zr;
  cplxreal(z, zc, zr);
  z = zc; z.insert(z.end(), zr.begin(), zr.end());
  cplxreal(p, zc, zr);
  p = zc; p.insert(p.end(), zr.begin(), zr.end());
  dmat sos(n_sections, 6, 0.0);
  auto count_real = [](const cdvec& v) { size_t c = 0; for (const cd& x : v) c += isreal(x); return c; };
  for (size_t si_ = n_sections; si_-- > 0;) {
    const size_t p1_idx = idx_worst(p);
    const cd p1 = p[p1_idx];
    p.erase(p.begin() + p1_idx);
    double* row = sos.row(si_);
    if (isreal(p1) && count_real(p) == 0) {
      const size_t z1_idx = nearest_idx(z, p1, 1);
      const cd z1 = z[z1_idx];
      z.erase(z.begin() + z1_idx);
      single_zpksos({z1, cd(0, 0)}, {p1, cd(0, 0)}, row);
    } else if (p.size() + 1 == z.size() && !isreal(p1) && count_real(p) == 1 && count_real(z) == 1) {
      const size_t z1_idx = nearest_idx(z, p1, 2);
      const cd z1 = z[z1_idx];
      z.erase(z.begin() + z1_idx);
      single_zpksos({z1, std::conj(z1)}, {p1, std::conj(p1)}, row);
    } else {
      cd p2;
      if (isreal(p1)) {
        std::vector<size_t> ridx;
        cdvec preal;
        for (size_t i = 0; i < p.size(); ++i)
          if (isreal(p[i])) { ridx.push_back(i); preal.push_back(p[i]); }
        const size_t p2_idx = ridx[idx_worst(preal)];
        p2 = p[p2_idx];
        p.erase(p.begin() + p2_idx);
      } else {
        p2 = std::conj(p1);
      }
      if (!z.empty()) {
        const size_t z1_idx = nearest_idx(z, p1, 0);
        const cd z1 = z[z1_idx];
        z.erase(z.begin() + z1_idx);
        if (!isreal(z1)) {
          single_zpksos({z1, std::conj(z1)}, {p1, p2}, row);
        } else if (!z.empty()) {
          const size_t z2_idx = nearest_idx(z, p1, 1);
          const cd z2 = z[z2_idx];
          z.erase(z.begin() + z2_idx);
          single_zpksos({z1, z2}, {p1, p2}, row);
        } else {
          single_zpksos({z1}, {p1, p2}, row);
        }
      } else {
        single_zpksos({}, {p1, p2}, row);
      }
    }
  }
  for (int i = 0; i < 3; ++i) sos(0, i) *= s.k;
  return sos;
}
}  // namespace

dmat butter_sos(int N, const dvec& Wn_in, BType btype, std::optional<double> fs) {
  if (N < 1) throw std::invalid_argument("butter: order must be >= 1");
  dvec Wn = Wn_in;
  if (fs)
    for (auto& w : Wn) w = w / (*fs / 2);
  for (double w : Wn)
    if (!(w > 0) || !(w < 1)) throw std::invalid_argument("Digital filter critical frequencies must be 0 < Wn < 1");
  ZPK s = buttap(N);
  const double fsd = 2.0;
  dvec warped(Wn.size());
  for (size_t i = 0; i < Wn.size(); ++i) warped[i] = 2 * fsd * std::tan(kPi * Wn[i] / fsd);
  if (btype == BType::Lowpass || btype == BType::Highpass) {
    if (Wn.size() != 1) throw std::invalid_argument("Must specify a single critical frequency Wn for lowpass or highpass filter");
    s = btype == BType::Lowpass ? lp2lp(s, warped[0]) : lp2hp(s, warped[0]);
  } else {
    if (Wn.size() != 2) throw std::invalid_argument("Wn must specify start and stop frequencies for bandpass or bandstop filter");
    const double bw = warped[1] - warped[0];
    const double wo = std::sqrt(warped[0] * warped[1]);
    s = lp2bp(s, wo, bw);
  }
  s = bilinear(s, fsd);
  return zpk2sos(s);
}

// ------------------------------------------------------------------ IIR filtering
dvec sosfilt(const dmat& sos, const dvec& x, dmat* zi) {
  if (sos.cols != 6) throw std::invalid_argument("sos array must be shape (n_sections, 6)");
  const size_t ns = sos.rows;
  dmat z(ns, 2, 0.0);
  if (zi) {
    if (zi->rows != ns || zi->cols != 2) throw std::invalid_argument("Invalid zi shape");
    z = *zi;
  }
  dvec y(x.size());
  for (size_t n = 0; n < x.size(); ++n) {
    double x_cur = x[n];
    for (size_t s = 0; s < ns; ++s) {
      const double* c = sos.row(s);
      const double x_new = c[0] * x_cur + z(s, 0);
      z(s, 0) = c[1] * x_cur - c[4] * x_new + z(s, 1);
      z(s, 1) = c[2] * x_cur - c[5] * x_new;
      x_cur = x_new;
    }
    y[n] = x_cur;
  }
  if (zi) *zi = z;
  return y;
}

dvec lfilter_zi(const dvec& b_in, const dvec& a_in) {
  dvec b = b_in, a = a_in;
  while (a.size() > 1 && a[0] == 0.0) a.erase(a.begin());
  if (a.empty()) throw std::invalid_argument("There must be at least one nonzero `a` coefficient.");
  if (a[0] != 1.0) {
    const double a0 = a[0];
    for (auto& v : b) v = v / a0;
    for (auto& v : a) v = v / a0;
  }
  const size_t n = std::max(a.size(), b.size());
  a.resize(n, 0.0);
  b.resize(n, 0.0);
  if (n < 2) return {};
  // IminusA = eye(n-1) - companion(a).T ; companion first row -a[1:]/a[0]
  dmat M(n - 1, n - 1, 0.0);
  for (size_t i = 0; i < n - 1; ++i) M(i, i) = 1.0;
  const double a0 = 1.0 * a[0];
  dmat C(n - 1, n - 1, 0.0);
  for (size_t j = 0; j < n - 1; ++j) C(0, j) = -a[j + 1] / a0;
  for (size_t i = 1; i < n - 1; ++i) C(i, i - 1) = 1.0;
  for (size_t i = 0; i < n - 1; ++i)
    for (size_t j = 0; j < n - 1; ++j) M(i, j) = M(i, j) - C(j, i);
  dvec B(n - 1);
  for (size_t i = 0; i < n - 1; ++i) B[i] = b[i + 1] - a[i + 1] * b[0];
  return np::solve(M, B);
}

dmat sosfilt_zi(const dmat& sos) {
  dmat zi(sos.rows, 2, 0.0);
  double scale = 1.0;
  for (size_t s = 0; s < sos.rows; ++s) {
    const dvec b{sos(s, 0), sos(s, 1), sos(s, 2)}, a{sos(s, 3), sos(s, 4), sos(s, 5)};
    const dvec z = lfilter_zi(b, a);
    zi(s, 0) = scale * z[0];
    zi(s, 1) = scale * z[1];
    scale *= np::pairwise_sum(b.data(), 3) / np::pairwise_sum(a.data(), 3);
  }
  return zi;
}

template <class T>
dvec sosfiltfilt(const dmat& sos, const std::vector<T>& x, PadType padtype, long padlen) {
  const size_t ns = sos.rows;
  size_t nz2 = 0, nz5 = 0;
  for (size_t s = 0; s < ns; ++s) { nz2 += sos(s, 2) == 0; nz5 += sos(s, 5) == 0; }
  const size_t ntaps = 2 * ns + 1 - std::min(nz2, nz5);
  size_t edge = padtype == PadType::None ? 0 : (padlen < 0 ? ntaps * 3 : static_cast<size_t>(padlen));
  if (x.size() <= edge)
    throw std::invalid_argument("The length of the input vector x must be greater than padlen, which is " + std::to_string(edge) + ".");
  std::vector<T> ext;
  if (padtype != PadType::None && edge > 0) {
    const size_t n = x.size();
    ext.reserve(n + 2 * edge);
    for (size_t i = edge; i >= 1; --i) {  // left: x[edge:0:-1]
      if (padtype == PadType::Odd) ext.push_back(T(2) * x[0] - x[i]);
      else if (padtype == PadType::Even) ext.push_back(x[i]);
      else ext.push_back(x[0]);
    }
    ext.insert(ext.end(), x.begin(), x.end());
    for (size_t i = 0; i < edge; ++i) {  // right: x[-2:-(edge+2):-1]
      const T v = x[n - 2 - i];
      if (padtype == PadType::Odd) ext.push_back(T(2) * x[n - 1] - v);
      else if (padtype == PadType::Even) ext.push_back(v);
      else ext.push_back(x[n - 1]);
    }
  } else {
    ext = x;
  }
  const dmat zi = sosfilt_zi(sos);
  dvec e(ext.begin(), ext.end());
  dmat z0 = zi;
  for (auto& v : z0.data) v = v * e[0];
  dvec y = sosfilt(sos, e, &z0);
  const double y0 = y.back();
  std::reverse(y.begin(), y.end());
  dmat z1 = zi;
  for (auto& v : z1.data) v = v * y0;
  y = sosfilt(sos, y, &z1);
  std::reverse(y.begin(), y.end());
  if (edge > 0) y = dvec(y.begin() + edge, y.end() - edge);
  return y;
}
template dvec sosfiltfilt(const dmat&, const std::vector<float>&, PadType, long);
template dvec sosfiltfilt(const dmat&, const std::vector<double>&, PadType, long);

dvec lfilter(const dvec& b_in, const dvec& a_in, const dvec& x, dvec* zi) {
  if (a_in.empty() || b_in.empty()) throw std::invalid_argument("lfilter: empty coefficients");
  if (a_in.size() == 1) {  // FIR path: np.convolve(b / a[0], x)
    dvec b = b_in;
    for (auto& v : b) v /= a_in[0];
    if (x.empty()) {
      if (zi) *zi = dvec(*zi);
      return {};
    }
    dvec full = np::convolve(b, x, np::ConvMode::Full);
    if (zi) {
      if (zi->size() != b.size() - 1) throw std::invalid_argument("Unexpected shape for zi");
      for (size_t i = 0; i < zi->size(); ++i) full[i] += (*zi)[i];
    }
    dvec y(full.begin(), full.begin() + x.size());
    if (zi) *zi = dvec(full.begin() + x.size(), full.end());
    return y;
  }
  const size_t nal = std::max(a_in.size(), b_in.size());
  dvec a(nal, 0.0), b(nal, 0.0);
  std::copy(a_in.begin(), a_in.end(), a.begin());
  std::copy(b_in.begin(), b_in.end(), b.begin());
  const double a0 = a[0];
  if (a0 == 0.0) throw std::invalid_argument("BUG: filter coefficient a[0] == 0 not supported yet");
  for (size_t n = 0; n < nal; ++n) { b[n] /= a0; a[n] /= a0; }
  dvec Z(nal - 1, 0.0);
  if (zi) {
    if (zi->size() != nal - 1) throw std::invalid_argument("The number of initial conditions must be max([len(a),len(b)]) - 1");
    Z = *zi;
  }
  dvec y(x.size());
  for (size_t k = 0; k < x.size(); ++k) {
    const double xn = x[k];
    const double yn = Z[0] + b[0] * xn;
    for (size_t n = 0; n + 2 < nal; ++n) Z[n] = Z[n + 1] + xn * b[n + 1] - yn * a[n + 1];
    Z[nal - 2] = xn * b[nal - 1] - yn * a[nal - 1];
    y[k] = yn;
  }
  if (zi) *zi = Z;
  return y;
}

// ------------------------------------------------------------------ FFT convolution
size_t next_fast_len(size_t n, bool real) { return real ? fft::good_size_real(n) : fft::good_size_cmplx(n); }

namespace {
// scipy.fft transforms (native precision) through core/fft, the only pocketfft user.
template <class R>
std::vector<std::complex<R>> pf_rfft(const R* x, size_t len, size_t n) {
  return fft::rfft(x, len, n);
}
template <class R>
std::vector<R> pf_irfft(const std::vector<std::complex<R>>& X, size_t n) {
  return fft::irfft(X.data(), X.size(), n);
}
// scipy.fft.fft(x, n) for real or complex x (real: r2c + conjugate-symmetric fill, as
// pypocketfft's c2c_sym)
template <class T>
std::vector<std::complex<real_t<T>>> pf_fft_full(const T* x, size_t len, size_t n) {
  using R = real_t<T>;
  if constexpr (is_cplx<T>::value) {
    return fft::fft(x, len, n);
  } else {
    std::vector<std::complex<R>> h = pf_rfft(x, len, n);
    std::vector<std::complex<R>> y(n);
    std::copy(h.begin(), h.end(), y.begin());
    for (size_t i = n / 2 + 1; i < n; ++i) y[i] = std::conj(y[n - i]);
    return y;
  }
}
template <class R>
std::vector<std::complex<R>> pf_ifft(std::vector<std::complex<R>> X, size_t n) {
  return fft::ifft(X.data(), std::min(X.size(), n), n);
}

template <class O, class I>
inline O promote(const I& v) {
  if constexpr (is_cplx<O>::value) {
    using R = real_t<O>;
    if constexpr (is_cplx<I>::value) return O(static_cast<R>(v.real()), static_cast<R>(v.imag()));
    else return O(static_cast<R>(v), R(0));
  } else {
    return static_cast<O>(v);
  }
}

// one linear convolution of length n1+n2-1 through FFTs of size N (scipy
// _freq_domain_conv for one block)
template <class O, class A, class B>
std::vector<O> freq_conv(const A* a, size_t na, const B* b, size_t nb, size_t N) {
  using R = real_t<O>;
  std::vector<O> out;
  if constexpr (!is_cplx<O>::value) {
    const auto Sa = pf_rfft<real_t<A>>(a, na, N);
    const auto Sb = pf_rfft<real_t<B>>(b, nb, N);
    std::vector<std::complex<R>> P(Sa.size());
    for (size_t i = 0; i < P.size(); ++i)
      P[i] = cmul(promote<std::complex<R>>(Sa[i]), promote<std::complex<R>>(Sb[i]));
    out = pf_irfft<R>(P, N);
  } else {
    const auto Sa = pf_fft_full(a, na, N);
    const auto Sb = pf_fft_full(b, nb, N);
    std::vector<std::complex<R>> P(N);
    for (size_t i = 0; i < N; ++i) P[i] = cmul(promote<std::complex<R>>(Sa[i]), promote<std::complex<R>>(Sb[i]));
    out = pf_ifft<R>(std::move(P), N);
  }
  return out;
}

template <class O>
std::vector<O> centered(const std::vector<O>& r, size_t newlen) {
  const size_t start = (r.size() - newlen) / 2;
  return std::vector<O>(r.begin() + start, r.begin() + start + newlen);
}

template <class O>
std::vector<O> apply_conv_mode(const std::vector<O>& ret, size_t s1, size_t s2, np::ConvMode mode, bool has_axis) {
  switch (mode) {
    case np::ConvMode::Full: return ret;
    case np::ConvMode::Same: return centered(ret, s1);
    case np::ConvMode::Valid: return centered(ret, has_axis ? s1 - s2 + 1 : ret.size());
  }
  return ret;
}

// in1 * in2 with broadcasting of a length-1 operand (scipy: no FFT axis)
template <class O, class A, class B>
std::vector<O> bcast_mul(const std::vector<A>& a, const std::vector<B>& b) {
  const size_t n = std::max(a.size(), b.size());
  std::vector<O> r(n);
  for (size_t i = 0; i < n; ++i) {
    const O x = promote<O>(a[a.size() == 1 ? 0 : i]), y = promote<O>(b[b.size() == 1 ? 0 : i]);
    if constexpr (is_cplx<O>::value) r[i] = cmul(x, y);
    else r[i] = x * y;
  }
  return r;
}

double lambertw_m1(double z) {  // W_{-1}(z), z in [-1/e, 0)
  double w;
  if (z < -0.25) {
    const double p = -std::sqrt(2.0 * (std::exp(1.0) * z + 1.0));
    w = -1.0 + p - p * p / 3.0 + 11.0 / 72.0 * p * p * p;
  } else {
    const double l1 = std::log(-z), l2 = std::log(-l1);
    w = l1 - l2 + l2 / l1;
  }
  for (int i = 0; i < 200; ++i) {
    const double ew = std::exp(w), f = w * ew - z;
    const double wn = w - f / (ew * (w + 1.0) - (w + 2.0) * f / (2.0 * w + 2.0));
    if (std::fabs(wn - w) <= 1e-16 * std::fabs(wn)) { w = wn; break; }
    w = wn;
  }
  return w;
}

struct OaLens {
  size_t block, overlap, step1, step2;
  bool fallback;
};
OaLens calc_oa_lens(size_t s1, size_t s2) {
  OaLens fb{s1 + s2 - 1, 0, s1, s2, true};
  if (s1 == s2 || s1 == 1 || s2 == 1) return fb;
  bool swapped = false;
  if (s2 > s1) { std::swap(s1, s2); swapped = true; }
  if (static_cast<double>(s2) >= static_cast<double>(s1) / 2) { if (swapped) std::swap(s1, s2); return OaLens{s1 + s2 - 1, 0, s1, s2, true}; }
  const size_t overlap = s2 - 1;
  const double e = 2.718281828459045;
  const double opt = -static_cast<double>(overlap) * lambertw_m1(-1.0 / (2 * e * static_cast<double>(overlap)));
  const size_t block = next_fast_len(static_cast<size_t>(std::ceil(opt)), false);
  if (block >= s1) { if (swapped) std::swap(s1, s2); return OaLens{s1 + s2 - 1, 0, s1, s2, true}; }
  if (!swapped) return {block, overlap, block - s2 + 1, s2, false};
  return {block, overlap, s2, block - s2 + 1, false};
}
}  // namespace

template <class A, class B>
std::vector<conv_out_t<A, B>> fftconvolve(const std::vector<A>& in1_, const std::vector<B>& in2_, np::ConvMode mode) {
  using O = conv_out_t<A, B>;
  if (in1_.empty() || in2_.empty()) return {};
  // 'valid' swaps so that in1 is the longer input (convolution commutes)
  if (mode == np::ConvMode::Valid && in2_.size() > in1_.size()) {
    if constexpr (std::is_same_v<A, B>) return fftconvolve(in2_, in1_, mode);
    else {
      auto r = fftconvolve(in2_, in1_, mode);
      return std::vector<O>(r.begin(), r.end());
    }
  }
  const size_t s1 = in1_.size(), s2 = in2_.size();
  if (s1 == 1 || s2 == 1) return apply_conv_mode(bcast_mul<O>(in1_, in2_), s1, s2, mode, false);
  const size_t shape = s1 + s2 - 1;
  const size_t N = next_fast_len(shape, !is_cplx<O>::value);
  std::vector<O> ret = freq_conv<O>(in1_.data(), s1, in2_.data(), s2, N);
  ret.resize(shape);
  return apply_conv_mode(ret, s1, s2, mode, true);
}

template <class A, class B>
std::vector<conv_out_t<A, B>> oaconvolve(const std::vector<A>& in1_, const std::vector<B>& in2_, np::ConvMode mode) {
  using O = conv_out_t<A, B>;
  if (in1_.empty() || in2_.empty()) return {};
  if (in1_.size() == in2_.size()) return fftconvolve(in1_, in2_, mode);
  if (mode == np::ConvMode::Valid && in2_.size() > in1_.size()) {
    auto r = oaconvolve(in2_, in1_, mode);
    return std::vector<O>(r.begin(), r.end());
  }
  const size_t s1 = in1_.size(), s2 = in2_.size();
  if (s1 == 1 || s2 == 1) return apply_conv_mode(bcast_mul<O>(in1_, in2_), s1, s2, mode, false);
  const size_t shape_final = s1 + s2 - 1;
  const OaLens L = calc_oa_lens(s1, s2);
  if (L.step1 == s1 && L.step2 == s2) return fftconvolve(in1_, in2_, mode);
  auto nsteps = [&](size_t s, size_t step, size_t& pad) {
    if (s > step) {
      size_t n = (s + 1 + step - 1) / step;  // ceil((s+1)/step)
      if ((L.block - L.overlap) * n < shape_final) n += 1;
      pad = n * step - s;
      return n;
    }
    pad = 0;
    return size_t(1);
  };
  size_t pad1, pad2;
  const size_t n1 = nsteps(s1, L.step1, pad1), n2 = nsteps(s2, L.step2, pad2);
  std::vector<A> a(in1_);
  a.resize(s1 + pad1, A(0));
  std::vector<B> b(in2_);
  b.resize(s2 + pad2, B(0));
  const size_t ns = std::max(n1, n2), blk = L.block;
  std::vector<std::vector<O>> rows(ns);
  for (size_t k = 0; k < ns; ++k) {
    const A* pa = a.data() + (n1 == 1 ? 0 : k * L.step1);
    const B* pb = b.data() + (n2 == 1 ? 0 : k * L.step2);
    rows[k] = freq_conv<O>(pa, L.step1, pb, L.step2, blk);
  }
  // overlap-add: block k's first `overlap` outputs += block k-1's last `overlap`
  const size_t ov = L.overlap, keep = blk - ov;
  for (size_t k = 1; k < ns; ++k)
    for (size_t i = 0; i < ov; ++i) rows[k][i] += rows[k - 1][keep + i];
  std::vector<O> ret;
  ret.reserve(ns * keep);
  for (size_t k = 0; k < ns; ++k) ret.insert(ret.end(), rows[k].begin(), rows[k].begin() + keep);
  ret.resize(shape_final);
  return apply_conv_mode(ret, s1, s2, mode, true);
}

#define SI_CONV(A, B)                                                                                                \
  template std::vector<conv_out_t<A, B>> fftconvolve(const std::vector<A>&, const std::vector<B>&, np::ConvMode); \
  template std::vector<conv_out_t<A, B>> oaconvolve(const std::vector<A>&, const std::vector<B>&, np::ConvMode);
SI_CONV(float, float)
SI_CONV(double, double)
SI_CONV(cf, cf)
SI_CONV(cd, cd)
SI_CONV(cd, double)
SI_CONV(double, cd)
SI_CONV(cf, double)
SI_CONV(double, cf)
SI_CONV(cf, float)
SI_CONV(float, cf)
#undef SI_CONV

// ------------------------------------------------------------------ resample_poly
namespace {
int64_t output_len(int64_t len_h, int64_t in_len, int64_t up, int64_t down) {
  const int64_t in_len_copy = in_len + (len_h + ((up - len_h % up) % up)) / up - 1;
  const int64_t nt = in_len_copy * up;
  int64_t need = nt / down;
  if (nt % down > 0) need += 1;
  return need;
}
}  // namespace

template <class T>
std::vector<T> resample_poly(const std::vector<T>& x, int up_, int down_, const Window& window) {
  using R = real_t<T>;
  if (up_ < 1 || down_ < 1) throw std::invalid_argument("up and down must be >= 1");
  const int64_t g = std::gcd(up_, down_);
  const int64_t up = up_ / g, down = down_ / g;
  if (up == 1 && down == 1) return x;
  const int64_t n_in = static_cast<int64_t>(x.size());
  int64_t n_out = n_in * up;
  n_out = n_out / down + (n_out % down ? 1 : 0);
  const int64_t max_rate = std::max(up, down);
  const double f_c = 1.0 / static_cast<double>(max_rate);
  const int64_t half_len = 10 * max_rate;
  const dvec hd = firwin(static_cast<size_t>(2 * half_len + 1), f_c, std::nullopt, window);
  std::vector<R> h0(hd.size());
  for (size_t i = 0; i < hd.size(); ++i) h0[i] = static_cast<R>(hd[i]) * static_cast<R>(up);
  const int64_t n_pre_pad = down - half_len % down;
  int64_t n_post_pad = 0;
  const int64_t n_pre_remove = (half_len + n_pre_pad) / down;
  while (output_len(static_cast<int64_t>(h0.size()) + n_pre_pad + n_post_pad, n_in, up, down) < n_out + n_pre_remove)
    ++n_post_pad;
  std::vector<R> h(static_cast<size_t>(n_pre_pad), R(0));
  h.insert(h.end(), h0.begin(), h0.end());
  h.resize(h.size() + static_cast<size_t>(n_post_pad), R(0));
  // upfirdn: _pad_h (transposed, flipped phases) then scipy's _apply_impl
  const int64_t len_h0 = static_cast<int64_t>(h.size());
  const int64_t h_padlen = len_h0 + ((up - len_h0 % up) % up);
  const int64_t hpp = h_padlen / up;  // h_per_phase
  std::vector<R> htf(static_cast<size_t>(h_padlen), R(0));
  for (int64_t ph = 0; ph < up; ++ph)
    for (int64_t j = 0; j < hpp; ++j) {
      const int64_t src = (hpp - 1 - j) * up + ph;
      htf[ph * hpp + j] = src < len_h0 ? h[src] : R(0);
    }
  const int64_t len_out = output_len(len_h0, n_in, up, down);
  std::vector<T> out(static_cast<size_t>(len_out), T(0));
  auto mac = [&](int64_t y, const T& xv, R hv) {
    if constexpr (is_cplx<T>::value) out[y] = T(out[y].real() + xv.real() * hv, out[y].imag() + xv.imag() * hv);
    else out[y] = out[y] + xv * hv;
  };
  const int64_t padded_len = n_in + hpp - 1;
  int64_t x_idx = 0, y_idx = 0, t = 0;
  bool done = len_out == 0;
  while (!done && x_idx < n_in) {
    int64_t h_idx = t * hpp;
    int64_t xc = x_idx - hpp + 1;
    if (xc < 0) { h_idx -= xc; xc = 0; }
    for (; xc < x_idx + 1; ++xc) mac(y_idx, x[xc], htf[h_idx++]);
    ++y_idx;
    if (y_idx >= len_out) { done = true; break; }
    t += down;
    x_idx += t / up;
    t = t % up;
  }
  while (!done && x_idx < padded_len) {
    int64_t h_idx = t * hpp;
    for (int64_t xc = x_idx - hpp + 1; xc < x_idx + 1; ++xc) {
      const T xv = (xc >= n_in || xc < 0) ? T(0) : x[xc];
      mac(y_idx, xv, htf[h_idx++]);
    }
    ++y_idx;
    if (y_idx >= len_out) break;
    t += down;
    x_idx += t / up;
    t = t % up;
  }
  return std::vector<T>(out.begin() + n_pre_remove, out.begin() + n_pre_remove + n_out);
}
template std::vector<float> resample_poly(const std::vector<float>&, int, int, const Window&);
template std::vector<double> resample_poly(const std::vector<double>&, int, int, const Window&);
template std::vector<cf> resample_poly(const std::vector<cf>&, int, int, const Window&);
template std::vector<cd> resample_poly(const std::vector<cd>&, int, int, const Window&);

// ------------------------------------------------------------------ welch
template <class T>
Welch<real_t<T>> welch(const std::vector<T>& x, double fs, size_t nperseg, long noverlap_, bool return_onesided) {
  using R = real_t<T>;
  using C = std::complex<R>;
  constexpr bool single = std::is_same_v<R, float>;
  constexpr bool cplx = is_cplx<T>::value;
  Welch<R> out;
  if (x.empty()) return out;
  if (x.size() < nperseg) nperseg = x.size();
  const size_t nfft = nperseg;
  const size_t noverlap = noverlap_ < 0 ? nperseg / 2 : static_cast<size_t>(noverlap_);
  if (noverlap >= nperseg) throw std::invalid_argument("noverlap must be less than nperseg.");
  const size_t nstep = nperseg - noverlap;
  const dvec wd = get_window(Window::hann(), nperseg, true);
  std::vector<R> w(nperseg);
  for (size_t i = 0; i < nperseg; ++i) w[i] = static_cast<R>(wd[i]);  // complex64 cast in single precision
  R scale;
  {
    std::vector<R> ww(nperseg);
    for (size_t i = 0; i < nperseg; ++i) ww[i] = w[i] * w[i];
    if constexpr (single) {
      std::vector<C> wc(nperseg);
      for (size_t i = 0; i < nperseg; ++i) wc[i] = C(ww[i], 0);
      const C s = np::pairwise_sum(wc.data(), nperseg);
      scale = R(1) / (static_cast<R>(fs) * s.real());
    } else {
      scale = 1.0 / (fs * np::pairwise_sum(ww.data(), nperseg));
    }
  }
  const bool onesided = return_onesided && !cplx;
  out.f = onesided ? fft::rfftfreq(nfft, 1.0 / fs) : fft::fftfreq(nfft, 1.0 / fs);
  const size_t nf = onesided ? nfft / 2 + 1 : nfft;
  const size_t nseg = (x.size() - nperseg) / nstep + 1;
  std::vector<R> acc(nf, R(0));
  std::vector<T> seg(nperseg);
  for (size_t s = 0; s < nseg; ++s) {
    const T* p = x.data() + s * nstep;
    const T m = np::mean(p, nperseg);
    for (size_t i = 0; i < nperseg; ++i) seg[i] = p[i] - m;
    std::vector<C> F;
    if (onesided) {
      std::vector<R> r(nperseg);
      if constexpr (!cplx)
        for (size_t i = 0; i < nperseg; ++i) r[i] = w[i] * seg[i];
      F = pf_rfft<R>(r.data(), nperseg, nfft);
    } else if constexpr (cplx) {
      std::vector<C> c(nperseg);
      for (size_t i = 0; i < nperseg; ++i) c[i] = C(w[i] * seg[i].real(), w[i] * seg[i].imag());
      F = pf_fft_full(c.data(), nperseg, nfft);
    } else if constexpr (single) {  // float32 two-sided: complex64 product -> c2c
      std::vector<C> c(nperseg);
      for (size_t i = 0; i < nperseg; ++i) c[i] = C(w[i] * seg[i], R(0));
      F = pf_fft_full(c.data(), nperseg, nfft);
    } else {  // float64 two-sided: real product -> scipy.fft.fft of a real array
      std::vector<R> r(nperseg);
      for (size_t i = 0; i < nperseg; ++i) r[i] = w[i] * seg[i];
      F = pf_fft_full(r.data(), nperseg, nfft);
    }
    for (size_t k = 0; k < nf; ++k) {
      R v = cmul(std::conj(F[k]), F[k]).real() * scale;
      if (onesided && k >= 1 && (nfft % 2 ? true : k < nf - 1)) v = v * R(2);
      acc[k] = (s == 0) ? v : acc[k] + v;
    }
  }
  out.Pxx.resize(nf);
  for (size_t k = 0; k < nf; ++k) out.Pxx[k] = nseg > 1 ? acc[k] / static_cast<R>(nseg) : acc[k];
  return out;
}
template Welch<float> welch(const std::vector<float>&, double, size_t, long, bool);
template Welch<double> welch(const std::vector<double>&, double, size_t, long, bool);
template Welch<float> welch(const std::vector<cf>&, double, size_t, long, bool);
template Welch<double> welch(const std::vector<cd>&, double, size_t, long, bool);

// ------------------------------------------------------------------ stft / istft
template <class C>
Stft<C> stft(const std::vector<C>& x, double fs, size_t nperseg, long noverlap_) {
  using R = real_t<C>;
  constexpr bool single = std::is_same_v<R, float>;
  Stft<C> out;
  if (x.size() < nperseg) nperseg = x.size();
  if (nperseg < 1) throw std::invalid_argument("nperseg must be a positive integer");
  const size_t noverlap = noverlap_ < 0 ? nperseg / 2 : static_cast<size_t>(noverlap_);
  if (noverlap >= nperseg) throw std::invalid_argument("noverlap must be less than nperseg.");
  const size_t nstep = nperseg - noverlap, nfft = nperseg;
  const dvec wd = get_window(Window::hann(), nperseg, true);
  dvec w(nperseg);
  double scale;
  if constexpr (single) {  // window cast to complex64; scale = sqrt(1 / sum(w)**2) in complex64
    std::vector<cf> wc(nperseg);
    for (size_t i = 0; i < nperseg; ++i) { wc[i] = cf(static_cast<float>(wd[i]), 0.f); w[i] = static_cast<float>(wd[i]); }
    const float s = np::pairwise_sum(wc.data(), nperseg).real();
    const float s2 = s * s;
    scale = static_cast<double>(std::sqrt(1.0f / s2));
  } else {
    w = wd;
    const double s = np::pairwise_sum(wd.data(), nperseg);
    scale = std::sqrt(1.0 / (s * s));
  }
  // zero_ext by nperseg//2 each side, then pad to a whole number of segments (complex128)
  const size_t half = nperseg / 2;
  const int64_t L0 = static_cast<int64_t>(x.size() + 2 * half);
  const int64_t ns = static_cast<int64_t>(nstep), np_ = static_cast<int64_t>(nperseg);
  int64_t nadd = (-(L0 - np_)) % ns;
  if (nadd < 0) nadd += ns;
  nadd = nadd % np_;
  cdvec xe(static_cast<size_t>(L0 + nadd), cd(0, 0));
  for (size_t i = 0; i < x.size(); ++i) xe[half + i] = cd(x[i].real(), x[i].imag());
  const size_t nseg = (xe.size() - nperseg) / nstep + 1;
  out.f = fft::fftfreq(nfft, 1.0 / fs);
  out.Z = Mat<C>(nfft, nseg);
  cdvec seg(nperseg);
  for (size_t s = 0; s < nseg; ++s) {
    for (size_t i = 0; i < nperseg; ++i) {
      const cd v = xe[s * nstep + i];
      seg[i] = cd(w[i] * v.real(), w[i] * v.imag());
    }
    const cdvec F = pf_fft_full(seg.data(), nperseg, nfft);
    for (size_t k = 0; k < nfft; ++k) {
      const cd v(F[k].real() * scale, F[k].imag() * scale);
      out.Z(k, s) = C(static_cast<R>(v.real()), static_cast<R>(v.imag()));
    }
  }
  const dvec tt = np::arange<double>(static_cast<double>(nperseg) / 2,
                                     static_cast<double>(xe.size()) - static_cast<double>(nperseg) / 2 + 1,
                                     static_cast<double>(nstep));
  out.t.resize(tt.size());
  for (size_t i = 0; i < tt.size(); ++i) out.t[i] = tt[i] / fs - (static_cast<double>(nperseg) / 2) / fs;
  return out;
}
template Stft<cf> stft(const std::vector<cf>&, double, size_t, long);
template Stft<cd> stft(const std::vector<cd>&, double, size_t, long);

template <class C>
Istft<C> istft(const Mat<C>& Z, double fs, size_t nperseg, long noverlap_) {
  using R = real_t<C>;
  constexpr bool single = std::is_same_v<R, float>;
  Istft<C> out;
  const size_t nfft = Z.rows, nseg = Z.cols;
  if (nperseg == 0) nperseg = nfft;
  if (nfft < nperseg) throw std::invalid_argument("nfft must be greater than or equal to nperseg.");
  const size_t noverlap = noverlap_ < 0 ? nperseg / 2 : static_cast<size_t>(noverlap_);
  if (noverlap >= nperseg) throw std::invalid_argument("noverlap must be less than nperseg.");
  const size_t nstep = nperseg - noverlap;
  const dvec wd = get_window(Window::hann(), nperseg, true);
  std::vector<R> w(nperseg);
  for (size_t i = 0; i < nperseg; ++i) w[i] = static_cast<R>(wd[i]);
  R wsum;
  if constexpr (single) {
    std::vector<C> wc(nperseg);
    for (size_t i = 0; i < nperseg; ++i) wc[i] = C(w[i], 0);
    wsum = np::pairwise_sum(wc.data(), nperseg).real();
  } else {
    wsum = np::pairwise_sum(w.data(), nperseg);
  }
  const size_t outlen = nperseg + (nseg - 1) * nstep;
  std::vector<C> xs(outlen, C(0, 0));
  std::vector<R> norm(outlen, R(0));
  std::vector<C> col(nfft);
  for (size_t s = 0; s < nseg; ++s) {
    for (size_t k = 0; k < nfft; ++k) col[k] = Z(k, s);
    const std::vector<C> sub = pf_ifft<R>(col, nfft);
    for (size_t i = 0; i < nperseg; ++i) {
      const C v(sub[i].real() * wsum, sub[i].imag() * wsum);
      const size_t j = s * nstep + i;
      xs[j] = C(xs[j].real() + v.real() * w[i], xs[j].imag() + v.imag() * w[i]);
      norm[j] = norm[j] + w[i] * w[i];
    }
  }
  const size_t h = nperseg / 2;
  const size_t lo = h, hi = outlen - h;
  out.x.resize(hi > lo ? hi - lo : 0);
  const R thr = static_cast<R>(1e-10);
  for (size_t i = lo; i < hi; ++i) {
    const R nv = norm[i] > thr ? norm[i] : R(1);
    const C d = cdiv(xs[i], C(nv, R(0)));
    out.x[i - lo] = d;
  }
  out.t.resize(out.x.size());
  for (size_t i = 0; i < out.t.size(); ++i) out.t[i] = static_cast<double>(i) / fs;
  return out;
}
template Istft<cf> istft(const Mat<cf>&, double, size_t, long);
template Istft<cd> istft(const Mat<cd>&, double, size_t, long);

}  // namespace si::sig
