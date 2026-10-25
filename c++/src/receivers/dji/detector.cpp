// DJI OcuSync video / RC uplink / DroneID (python: receivers/dji/detector.py). Line-by-line
// port; numpy exactness notes in README.md.
#include "receivers/dji/detector.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <tuple>

#include "core/bursts.hpp"
#include "core/fft.hpp"
#include "core/ndimage.hpp"
#include "core/npcompat.hpp"
#include "core/ofdm15k.hpp"
#include "core/parallel.hpp"
#include "core/signal.hpp"
#include "core/spectral.hpp"

namespace si::rx::dji {

const std::vector<VideoProfile>& video_profiles() {
  static const std::vector<VideoProfile> p = {
      {"15k_10MHz", 15, 15.36e6, 1024, 72, 9.0e6, 1.0e6, 601},
      {"15k_20MHz", 15, 30.72e6, 2048, 144, 18.0e6, 2.0e6, 1201},
      {"30k_10MHz", 30, 15.36e6, 512, 36, 9.0e6, 1.0e6, 0},
      {"30k_20MHz", 30, 30.72e6, 1024, 72, 18.0e6, 2.0e6, 0},
  };
  return p;
}

namespace {

constexpr double SLOT_S = 0.5e-3;
constexpr double PI = 3.141592653589793;
constexpr size_t UP_L = 73;
constexpr size_t DID_N = 1024;
const size_t DID_CP[9] = {80, 72, 72, 72, 72, 72, 72, 72, 80};
constexpr size_t DID_LEN = 9 * 1024 + 80 * 2 + 72 * 7;
// numpy elides the temporary of `a * np.conj(b)` (>= 256 KiB) and then computes conj(b) * a
constexpr size_t ELIDE_CPLX = 256 * 1024 / 16;

// np.exp(0 + 1j*im) (numpy -> libm cexp; on Apple that is __sincos, which can differ from
// sin() in the last bit)
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
  int64_t r = a % b;
  return r < 0 ? r + b : r;
}

// Python float // float
double py_floordiv(double vx, double wx) {
  double mod = std::fmod(vx, wx);
  double div = (vx - mod) / wx;
  if (mod != 0.0) {
    if ((wx < 0) != (mod < 0)) div -= 1.0;
  }
  double fl;
  if (div != 0.0) {
    fl = std::floor(div);
    if (div - fl > 0.5) fl += 1.0;
  } else {
    fl = std::copysign(0.0, vx / wx);
  }
  return fl;
}

inline int64_t npround_int(double x) { return static_cast<int64_t>(std::nearbyint(x)); }

// fftshift(fft(y[s:s+N])) (complex128)
cdvec sym_spectrum(const cd* y, size_t N) { return fft::fftshift(fft::fft(y, N)); }

// y * exp(-2j*pi*shift*n/ofs) as numpy computes it (elided operand order when large)
cdvec mix(const cdvec& y, double shift_hz, double ofs) {
  const double c = (-2.0 * PI) * shift_hz, inv = 1.0 / ofs;
  const bool swap = y.size() >= ELIDE_CPLX;
  cdvec out(y.size());
  for (size_t n = 0; n < y.size(); ++n) {
    const cd e = expj((c * static_cast<double>(n)) * inv);
    out[n] = swap ? np::cmul(e, y[n]) : np::cmul(y[n], e);
  }
  return out;
}

// np.convolve(v, ones(L), "valid")
dvec box_valid(const dvec& v, size_t L) { return np::convolve(v, dvec(L, 1.0), np::ConvMode::Valid); }

// sum of a selection, numpy pairwise
cd sum_sel(const cdvec& v, const std::vector<int64_t>& idx) {
  cdvec s;
  s.reserve(idx.size());
  for (int64_t i : idx) s.push_back(v[static_cast<size_t>(i)]);
  return np::sum(s);
}

// ============================================================== in-band SNR (python _SNR)
class Snr {
 public:
  Snr(const cfvec& x, double fs, double fc, size_t nfft = 1024) {
    SpectralSpectrogram sp = spectrogram(x, fs, fc, nfft);
    const dvec w = np::hanning(nfft);
    dvec w2(nfft);
    for (size_t i = 0; i < nfft; ++i) w2[i] = w[i] * w[i];
    const double sw = np::sum(w2);
    P_ = std::move(sp.P);
    for (auto& v : P_.data) v = v / sw;
    t_ = std::move(sp.t_us);
    f_ = std::move(sp.f_mhz);
    fl_ = noise_floor_psd(x, nfft).floor_power_per_sample;
  }
  std::optional<double> operator()(double t0, double t1, double f0, double f1) const {
    std::vector<size_t> ti, fi;
    for (size_t i = 0; i < t_.size(); ++i)
      if (t_[i] >= t0 && t_[i] <= t1) ti.push_back(i);
    for (size_t i = 0; i < f_.size(); ++i)
      if (f_[i] >= f0 && f_[i] <= f1) fi.push_back(i);
    if (ti.empty() || fi.empty()) return std::nullopt;
    dvec sig(fi.size()), fl(fi.size(), fl_);
    for (size_t j = 0; j < fi.size(); ++j) {
      double s = P_(ti[0], fi[j]);
      for (size_t r = 1; r < ti.size(); ++r) s += P_(ti[r], fi[j]);
      sig[j] = s / static_cast<double>(ti.size()) - fl_;
    }
    const double v = 10 * std::log10(std::max(np::mean(sig), 1e-12) / np::mean(fl));
    return py_round(v, 1);
  }

 private:
  dmat P_;
  dvec t_, f_;
  double fl_ = 0;
};

Json cfo_feats(double carrier_hz, Json f, double raster = 100e3) {
  const double nom = std::nearbyint(carrier_hz / raster) * raster;
  const double cfo = carrier_hz - nom;
  f["carrier_mhz"] = np::round(carrier_hz / 1e6, 5);
  f["cfo_hz"] = py_round(cfo, 1);
  f["cfo_ppm"] = py_round(cfo / carrier_hz * 1e6, 3);
  return f;
}

}  // namespace

// ============================================================== small DSP helpers
namespace detail {

int64_t FastSlicer::slice_len(double out_fs) const {
  return py_round_int(static_cast<double>(X_.size()) * out_fs / fs_);
}

cdvec FastSlicer::slice(double f_off, double half_bw, double out_fs) const {
  const int64_t L = static_cast<int64_t>(X_.size());
  const int64_t M = py_round_int(static_cast<double>(L) * out_fs / fs_);
  const double df = fs_ / static_cast<double>(L);
  const int64_t lo = static_cast<int64_t>(std::ceil((f_off - half_bw) / df - 1e-9));
  const int64_t hi = static_cast<int64_t>(std::floor((f_off + half_bw) / df + 1e-9));
  const int64_t k0 = py_round_int(f_off / fs_ * static_cast<double>(L));
  cdvec Y(static_cast<size_t>(M), cd(0.0, 0.0));
  for (int64_t k = lo; k <= hi; ++k) Y[static_cast<size_t>(pymod(k - k0, M))] = X_[static_cast<size_t>(pymod(k, L))];
  cdvec y = fft::ifft(Y);
  const double s = static_cast<double>(M) / static_cast<double>(L);
  for (auto& v : y) v = cd(v.real() * s, v.imag() * s);
  return y;
}

CpOut cp_metric(const cd* y, size_t len, size_t N, size_t cp, const double* a2) {
  CpOut r;
  if (len <= N) return r;
  const size_t n = len - N;
  const bool swap = n >= ELIDE_CPLX;
  cdvec c(n + 1);
  dvec ce(n + 1);
  c[0] = cd(0.0, 0.0);
  ce[0] = 0.0;
  // products and energies first (vectorised), then the running sums (sequential, as np.cumsum)
  for (size_t i = 0; i < n; ++i) {
    c[i + 1] = swap ? np::cmul(std::conj(y[i + N]), y[i]) : np::cmul(y[i], std::conj(y[i + N]));
    ce[i + 1] = a2 ? (a2[i] + a2[i + N]) / 2 : (np::cabs2(y[i]) + np::cabs2(y[i + N])) / 2;
  }
  for (size_t i = 0; i < n; ++i) {
    c[i + 1] = c[i] + c[i + 1];
    ce[i + 1] = ce[i] + ce[i + 1];
  }
  if (n + 1 <= cp) return r;
  const size_t nv = n + 1 - cp;
  r.m.resize(nv);
  r.cs.resize(nv);
  const double* pc = reinterpret_cast<const double*>(c.data());
  double* ps = reinterpret_cast<double*>(r.cs.data());
  const double* pe = ce.data();
  double* pm = r.m.data();
  for (size_t j = 0; j < nv; ++j) {
    const double re = pc[2 * (j + cp)] - pc[2 * j], im = pc[2 * (j + cp) + 1] - pc[2 * j + 1];
    ps[2 * j] = re;
    ps[2 * j + 1] = im;
    pm[j] = np::cabs(cd(re, im)) / ((pe[j + cp] - pe[j]) + 1e-20);
  }
  return r;
}

FoldIdx fold_index(size_t n, double P) {
  FoldIdx f;
  f.idx.resize(n);
  for (size_t i = 0; i < n; ++i) f.idx[i] = static_cast<int64_t>(std::floor(np::mod(static_cast<double>(i), P)));
  f.cnt = np::bincount(f.idx, static_cast<size_t>(std::ceil(P)));
  return f;
}

dvec fold(const dvec& v, const FoldIdx& f) {
  dvec s(f.cnt.size(), 0.0);
  for (size_t i = 0; i < v.size(); ++i) s[static_cast<size_t>(f.idx[i])] += v[i];
  for (size_t i = 0; i < s.size(); ++i) s[i] = s[i] / static_cast<double>(std::max<int64_t>(f.cnt[i], 1));
  return s;
}

dvec fold(const dvec& v, double P) { return fold(v, fold_index(v.size(), P)); }

dvec robust_z(const dvec& prof) {
  const double med = np::median(prof);
  dvec d(prof.size());
  for (size_t i = 0; i < prof.size(); ++i) d[i] = std::fabs(prof[i] - med);
  const double sig = 1.4826 * np::median(d) + 1e-12;
  dvec z(prof.size());
  for (size_t i = 0; i < prof.size(); ++i) z[i] = (prof[i] - med) / sig;
  return z;
}

Segments segment(const dvec& m, int64_t phase, double P, size_t n_avg, double k_sigma) {
  Segments S;
  const double Kf = py_floordiv(static_cast<double>(m.size()) - static_cast<double>(phase) - P, P);
  const int64_t K = Kf > 0 ? static_cast<int64_t>(Kf) : 0;
  S.t.resize(static_cast<size_t>(K));
  for (int64_t k = 0; k < K; ++k) S.t[static_cast<size_t>(k)] = static_cast<double>(phase) + P * static_cast<double>(k);
  auto samp = [&](double dt) {
    dvec o(S.t.size());
    const int64_t n = static_cast<int64_t>(m.size());
    for (size_t i = 0; i < S.t.size(); ++i) {
      const int64_t r = npround_int(S.t[i] + dt);
      const int64_t a = std::max<int64_t>(r - 2, 0), b = std::min<int64_t>(r + 3, n);
      double mx = m[static_cast<size_t>(a)];
      for (int64_t j = a + 1; j < b; ++j) mx = std::max(mx, m[static_cast<size_t>(j)]);
      o[i] = mx;
    }
    return o;
  };
  S.on = samp(0);
  S.off = samp(P / 2);
  const dvec q = samp(P / 4);
  const size_t n = S.t.size();
  if (n == 0) return S;
  dvec nd(n), dd(n);
  for (size_t i = 0; i < n; ++i) {
    nd[i] = S.off[i] - q[i];
    dd[i] = S.on[i] - S.off[i];
  }
  const dvec nref = ndi::uniform_filter1d(nd, n_avg);
  const double mn = np::median(nref);
  dvec ad(n);
  for (size_t i = 0; i < n; ++i) ad[i] = std::fabs(nref[i] - mn);
  const double sig = 1.4826 * np::median(ad) + 1e-9;
  const dvec d = ndi::uniform_filter1d(dd, n_avg);
  std::vector<uint8_t> mask(n), dil(n), clo(n);
  for (size_t i = 0; i < n; ++i) mask[i] = d[i] > k_sigma * sig;
  // binary_closing(mask, ones(3)): dilation then erosion, border 0 (edges erode)
  for (size_t i = 0; i < n; ++i) dil[i] = mask[i] || (i > 0 && mask[i - 1]) || (i + 1 < n && mask[i + 1]);
  for (size_t i = 0; i < n; ++i) clo[i] = dil[i] && i > 0 && dil[i - 1] && i + 1 < n && dil[i + 1];
  for (size_t i = 0; i < n;) {
    if (!clo[i]) {
      ++i;
      continue;
    }
    size_t a = i, b = i;
    while (b < n && clo[b]) ++b;
    i = b;
    while (a < b && S.on[a] - S.off[a] < sig) ++a;
    while (b > a && S.on[b - 1] - S.off[b - 1] < sig) --b;
    if (b - a >= 3) S.runs.emplace_back(a, b);
  }
  return S;
}

std::pair<cd, double> zc_dd(const cdvec& X, size_t L) {
  const size_t N = X.size(), h = (L - 1) / 2;
  cd acc(0.0, 0.0);
  double mag = 0.0;
  for (int s = 0; s < 2; ++s) {
    const cd* seg = s == 0 ? X.data() + (N / 2 - h) : X.data() + (N / 2 + 1);
    cdvec D(h - 1), DD(h - 2);
    for (size_t i = 0; i + 1 < h; ++i) D[i] = np::cmul(seg[i + 1], std::conj(seg[i]));
    dvec a(h - 2);
    for (size_t i = 0; i + 2 < h; ++i) {
      DD[i] = np::cmul(D[i + 1], std::conj(D[i]));
      a[i] = np::cabs(DD[i]);
    }
    acc += np::sum(DD);
    mag += np::sum(a);
  }
  return {acc, mag};
}

int64_t root_from_dd(cd acc, size_t L) {
  const double v = -std::atan2(acc.imag(), acc.real()) * static_cast<double>(L) / (2 * PI);
  return pymod(npround_int(v), static_cast<int64_t>(L));
}

dvec droneid_cp_score(const cdvec& y) {
  const size_t N = DID_N;
  if (y.size() <= DID_LEN) return {};
  const size_t n = y.size() - N, nmax = y.size() - DID_LEN;
  const bool swap = n >= ELIDE_CPLX;
  cdvec cp(n + 1);
  dvec ca(n + 1), cb(n + 1);
  cp[0] = cd(0.0, 0.0);
  ca[0] = cb[0] = 0.0;
  for (size_t i = 0; i < n; ++i) {
    const cd p = swap ? np::cmul(std::conj(y[i]), y[i + N]) : np::cmul(y[i + N], std::conj(y[i]));
    cp[i + 1] = cp[i] + p;
    ca[i + 1] = ca[i] + np::cabs2(y[i]);
    cb[i + 1] = cb[i] + np::cabs2(y[i + N]);
  }
  dvec score(nmax + 1, 0.0);
  size_t off = 0;
  for (int i = 0; i < 9; ++i) {
    const size_t c = DID_CP[i];
    if (i) {
      for (size_t j = 0; j <= nmax; ++j) {
        const size_t k = off + j;
        const double v = np::cabs(cd(cp[k + c] - cp[k])) /
                         std::sqrt(std::max((ca[k + c] - ca[k]) * (cb[k + c] - cb[k]), 1e-30));
        score[j] += v;
      }
    }
    off += N + c;
  }
  for (auto& v : score) v = v / 8;
  return score;
}

std::optional<FramePeriod> frame_period(const dvec& starts_ms) {
  dvec s = starts_ms;
  std::sort(s.begin(), s.end());
  if (s.size() < 3) return std::nullopt;
  if (auto T = ofdm15k::cycle(s)) return FramePeriod{*T, true};
  const dvec iv = np::diff(s);
  const double miv = np::median(iv);
  std::optional<FramePeriod> best;
  for (double Tc : np::arange(1.0, 10.01, 0.5)) {
    size_t hit = 0;
    for (double v : s) {
      const double r = np::mod(v - s[0] + Tc / 2, Tc) - Tc / 2;
      hit += std::fabs(r) < 0.1;
    }
    if (static_cast<double>(hit) / static_cast<double>(s.size()) >= 0.8 && Tc <= miv + 0.1) best = FramePeriod{Tc, false};
  }
  return best;
}

}  // namespace detail

using namespace detail;

cdvec zc_freq(int64_t root, size_t L) {
  const double c = -PI * static_cast<double>(root), inv = 1.0 / static_cast<double>(L);
  const size_t h = L / 2;
  cdvec out(L, cd(0.0, 0.0));
  for (size_t n = 0; n < L; ++n) {
    if (n == h) continue;
    const double nn = static_cast<double>(n);
    out[n] = expj(((c * nn) * (nn + 1.0)) * inv);
  }
  return out;
}

std::pair<double, int64_t> zc_root_search(const cdvec& X, size_t L, const std::vector<int64_t>& roots_in) {
  std::vector<int64_t> roots = roots_in;
  if (roots.empty())
    for (size_t r = 1; r < L; ++r) roots.push_back(static_cast<int64_t>(r));
  const double nx = np::norm(X) + 1e-20;
  cdvec Xn(X.size());
  for (size_t i = 0; i < X.size(); ++i) Xn[i] = np::cdiv(X[i], cd(nx, 0.0));
  const cdvec F = fft::fft(Xn.data(), Xn.size(), 2 * L);
  std::pair<double, int64_t> best{0.0, -1};
  cdvec G(2 * L);
  for (int64_t r : roots) {
    const cdvec z = zc_freq(r, L);
    const cdvec Z = fft::fft(z.data(), z.size(), 2 * L);
    for (size_t i = 0; i < 2 * L; ++i) G[i] = np::cmul(F[i], std::conj(Z[i]));
    const cdvec g = fft::ifft(G);
    double mx = np::cabs(g[0]);
    for (size_t i = 1; i < g.size(); ++i) mx = std::max(mx, np::cabs(g[i]));
    const double cc = mx / np::norm(z);
    if (cc > best.first) best = {cc, r};
  }
  return best;
}

namespace {

// ============================================================== video
struct Cand {
  const VideoProfile* prof;
  double f_off, z;  // z: np.float64
  int64_t phase;
  double P;
};

struct Evidence {
  std::string kind;
  double z = 0;
  bool has_ok = false, ok = true;
  int64_t zc_root = -1;  // -1 = None
  int bin_shift = 0;
  double zc_corr = 0, zc_ref = 0;  // rounded to 3 decimals, as python's evidence dict
  bool okay() const { return !has_ok || ok; }
};

struct Link {
  Cand c;
  double f_off = 0, carrier = 0;
  cdvec y;
  dvec m;
  std::vector<std::pair<size_t, size_t>> runs;
  dvec t;
  std::vector<size_t> syms;
  std::vector<size_t> ev_syms;  // symbols the full-band test finds on its own (evidence, clock fit)
  dvec S;
  double S_off = 0;       // slice centre of S
  bool retimed = false;   // one symbol grid per burst (no capture-level grid)
  Evidence ev;
  bool image = false;
  int64_t parent = -1;    // a further symbol grid of links[parent] (regrid)
  Evidence grid_ev;       // that grid's own DJI evidence (informative)
};

std::vector<Cand> scan_video(const FastSlicer& ch, double usable_half, double z_thr) {
  struct Group {
    double ofs, occ, step;
    std::vector<const VideoProfile*> profs;
  };
  std::vector<Group> groups;
  for (const auto& p : video_profiles()) {
    auto it = std::find_if(groups.begin(), groups.end(),
                           [&](const Group& g) { return g.ofs == p.ofs && g.occ == p.occ && g.step == p.step; });
    if (it == groups.end()) groups.push_back({p.ofs, p.occ, p.step, {&p}});
    else it->profs.push_back(&p);
  }
  std::vector<Cand> cands;
  for (const auto& g : groups) {
    const dvec offs = np::arange(-usable_half + g.occ / 2, usable_half - g.occ / 2 + 1, g.step);
    using Z = std::vector<std::pair<double, int64_t>>;
    // every slice of the group has the same length: one fold index per profile (python _FOLD_IDX)
    const size_t M = static_cast<size_t>(ch.slice_len(g.ofs));
    std::vector<FoldIdx> fidx;
    for (const VideoProfile* p : g.profs)
      fidx.push_back(fold_index(M > p->N + p->cp ? M - p->N - p->cp + 1 : 0, SLOT_S * g.ofs / (p->scs == 15 ? 7 : 14)));
    const std::vector<Z> res = pmap(offs, [&](double off) {
      const cdvec y = ch.slice(off, g.occ / 2 * 1.05, g.ofs);
      dvec a2(y.size());
      for (size_t i = 0; i < y.size(); ++i) a2[i] = np::cabs2(y[i]);
      Z out;
      for (size_t pi = 0; pi < g.profs.size(); ++pi) {
        const VideoProfile* p = g.profs[pi];
        const CpOut co = cp_metric(y.data(), y.size(), p->N, p->cp, a2.data());
        const dvec z = robust_z(co.m.size() == fidx[pi].idx.size() ? fold(co.m, fidx[pi])
                                                                   : fold(co.m, SLOT_S * g.ofs / (p->scs == 15 ? 7 : 14)));
        const size_t k = np::argmax(z);
        out.emplace_back(z[k], static_cast<int64_t>(k));
      }
      return out;
    });
    for (size_t pi = 0; pi < g.profs.size(); ++pi) {
      const VideoProfile* p = g.profs[pi];
      const double P = SLOT_S * g.ofs / (p->scs == 15 ? 7 : 14);
      const size_t n = offs.size();
      for (size_t i = 0; i < n; ++i) {
        const double zi = res[i][pi].first;
        if (zi < z_thr || zi < res[i ? i - 1 : 0][pi].first || zi < res[std::min(i + 1, n - 1)][pi].first) continue;
        cands.push_back({p, offs[i], zi, res[i][pi].second, P});
      }
    }
  }
  return cands;
}

std::vector<size_t> run_syms(const std::vector<std::pair<size_t, size_t>>& runs) {
  std::vector<size_t> s;
  for (const auto& [a, b] : runs)
    for (size_t k = a; k < b; ++k) s.push_back(k);
  return s;
}

struct CarrierOut {
  double car;
  dvec S;
};

CarrierOut find_carrier(const FastSlicer& ch, double f_off, double ofs, size_t N, size_t cp, const dvec& t_syms,
                        size_t L_used) {
  const double scs = ofs / static_cast<double>(N);
  const cdvec y = ch.slice(f_off, ofs / 2 * 0.98, ofs);
  const CpOut co = cp_metric(y, N, cp);
  std::vector<int64_t> idx;
  for (double t : t_syms) {
    const int64_t r = npround_int(t);
    if (r < static_cast<int64_t>(co.cs.size())) idx.push_back(r);
  }
  const cd s = sum_sel(co.cs, idx);
  const double frac = -std::atan2(s.imag(), s.real()) / (2 * PI) * scs;
  // python _med_spectrum: per-bin MEDIAN over the symbols (an interferer in fewer than half of
  // them, Wi-Fi 14 dB above the video on D_20260515-192835, must not move the band position)
  const cdvec yc = mix(y, frac, ofs);
  std::vector<dvec> rows;
  for (int64_t t : idx) {
    const int64_t st = t + static_cast<int64_t>(cp);
    if (st >= 0 && st + static_cast<int64_t>(N) <= static_cast<int64_t>(yc.size())) {
      const cdvec X = sym_spectrum(yc.data() + st, N);
      dvec r(N);
      for (size_t k = 0; k < N; ++k) r[k] = np::cabs2(X[k]);
      rows.push_back(std::move(r));
    }
  }
  dvec S(N, 0.0);
  if (!rows.empty()) {
    dvec col(rows.size());
    for (size_t k = 0; k < N; ++k) {
      for (size_t i = 0; i < rows.size(); ++i) col[i] = rows[i][k];
      S[k] = np::median(col);
    }
  }
  const dvec w = box_valid(S, L_used);
  const int64_t c0 = static_cast<int64_t>(np::argmax(w)) + static_cast<int64_t>(L_used / 2);
  const dvec loc = ndi::uniform_filter1d(S, 7);
  int64_t dip = -1;
  double best = 0;
  for (int64_t k = std::max<int64_t>(c0 - 3, 1); k < std::min<int64_t>(c0 + 4, static_cast<int64_t>(N) - 1); ++k) {
    const double v = S[static_cast<size_t>(k)] / (loc[static_cast<size_t>(k)] + 1e-20);
    if (dip < 0 || v < best) dip = k, best = v;
  }
  const int64_t k_int = dip - static_cast<int64_t>(N / 2);
  return {f_off + static_cast<double>(k_int) * scs + frac, S};
}

// python _segment with a list of metrics (full band + sub-bands, same timing): a symbol is on
// when ANY of them shows the CP structure. Returns the runs, the grid and the smoothed z of the
// first (full-band) metric. exclude: sample intervals of another known emission (in-band hops).
struct SegMulti {
  std::vector<std::pair<size_t, size_t>> runs;
  dvec t, z_full, z1;  // z1: per-symbol z of the full-band metric, -inf where excluded
};

SegMulti segment_multi(const std::vector<dvec>& ms, int64_t phase, double P, size_t min_len,
                       const std::vector<std::pair<double, double>>& exclude, size_t n_avg = 4, double k_sigma = 4.0) {
  SegMulti S;
  const double Kf = py_floordiv(static_cast<double>(ms[0].size()) - static_cast<double>(phase) - P, P);
  const int64_t K = Kf > 0 ? static_cast<int64_t>(Kf) : 0;
  S.t.resize(static_cast<size_t>(K));
  for (int64_t k = 0; k < K; ++k) S.t[static_cast<size_t>(k)] = static_cast<double>(phase) + P * static_cast<double>(k);
  const size_t n = S.t.size();
  if (n == 0) return S;
  const bool multi = ms.size() > 1;
  dvec z(n, -std::numeric_limits<double>::infinity()), zraw(n, -std::numeric_limits<double>::infinity());
  for (size_t mi = 0; mi < ms.size(); ++mi) {
    const dvec& m = ms[mi];
    auto samp = [&](double dt) {
      dvec o(n);
      const int64_t nm = static_cast<int64_t>(m.size());
      for (size_t i = 0; i < n; ++i) {
        const int64_t r = npround_int(S.t[i] + dt);
        const int64_t a = std::max<int64_t>(r - 2, 0), b = std::min<int64_t>(r + 3, nm);
        double mx = m[static_cast<size_t>(a)];
        for (int64_t j = a + 1; j < b; ++j) mx = std::max(mx, m[static_cast<size_t>(j)]);
        o[i] = mx;
      }
      return o;
    };
    const dvec on = samp(0), off = samp(P / 2), q4 = samp(P / 4);
    dvec q(n), dd(n);
    for (size_t i = 0; i < n; ++i) {
      q[i] = off[i] - q4[i];
      dd[i] = on[i] - off[i];
    }
    const dvec nref = ndi::uniform_filter1d(q, n_avg);
    const double mn = np::median(nref);
    dvec ad(n);
    for (size_t i = 0; i < n; ++i) ad[i] = std::fabs(nref[i] - mn);
    const double sig = 1.4826 * np::median(ad) + 1e-9;
    double sr = sig;
    if (multi) {
      const double mq = np::median(q);
      for (size_t i = 0; i < n; ++i) ad[i] = std::fabs(q[i] - mq);
      sr = 1.4826 * np::median(ad) + 1e-9;
    }
    const dvec d = ndi::uniform_filter1d(dd, n_avg);
    if (mi == 0) {
      S.z_full.resize(n);
      S.z1.resize(n);
      for (size_t i = 0; i < n; ++i) {
        S.z_full[i] = d[i] / sig;
        S.z1[i] = dd[i] / sr;
      }
    }
    for (size_t i = 0; i < n; ++i) {
      z[i] = std::max(z[i], d[i] / sig);
      zraw[i] = std::max(zraw[i], dd[i] / sr);
    }
  }
  for (size_t i = 0; i < n; ++i)
    for (const auto& [s_, e_] : exclude)
      if (S.t[i] < e_ && S.t[i] + P > s_) {
        z[i] = -std::numeric_limits<double>::infinity();
        zraw[i] = -std::numeric_limits<double>::infinity();
        S.z1[i] = -std::numeric_limits<double>::infinity();
      }
  const double trim = multi ? 3.0 : 1.0;
  std::vector<uint8_t> mask(n), dil(n), clo(n);
  for (size_t i = 0; i < n; ++i) mask[i] = z[i] > k_sigma;
  for (size_t i = 0; i < n; ++i) dil[i] = mask[i] || (i > 0 && mask[i - 1]) || (i + 1 < n && mask[i + 1]);
  for (size_t i = 0; i < n; ++i) clo[i] = dil[i] && i > 0 && dil[i - 1] && i + 1 < n && dil[i + 1];
  for (size_t i = 0; i < n;) {
    if (!clo[i]) {
      ++i;
      continue;
    }
    size_t a = i, b = i;
    while (b < n && clo[b]) ++b;
    i = b;
    while (a < b && zraw[a] < trim) ++a;
    while (b > a && zraw[b - 1] < trim) --b;
    if (multi && b > a) {
      while (a > 0 && zraw[a - 1] >= trim) --a;
      while (b < n && zraw[b] >= trim) ++b;
    }
    if (b - a >= min_len) {
      if (!S.runs.empty() && a <= S.runs.back().second) S.runs.back().second = std::max(b, S.runs.back().second);
      else S.runs.emplace_back(a, b);
    }
  }
  return S;
}

constexpr size_t N_SUB = 4;
constexpr double MAX_BRIDGE_S = 0.6e-3;

// python _sub_slices: N_SUB equal sub-bands of the occupied band, same rate and timing
std::vector<cdvec> sub_slices(const FastSlicer& ch, double car, double ofs, double occ) {
  const double w = occ / static_cast<double>(N_SUB);
  std::vector<cdvec> out;
  for (size_t i = 0; i < N_SUB; ++i)
    out.push_back(ch.slice(car + (static_cast<double>(i) - static_cast<double>(N_SUB - 1) / 2) * w, w / 2, ofs));
  return out;
}

// python _bridge: join runs across a short gap in which every sub-band carries at least the
// power of the neighbouring link symbols (a stronger signal covers the whole band)
std::vector<std::pair<size_t, size_t>> bridge(const std::vector<std::pair<size_t, size_t>>& runs, const dvec& t,
                                              double P, const std::vector<cdvec>& subs, size_t max_gap) {
  if (runs.size() < 2) return runs;
  std::vector<dvec> cums;
  for (const auto& y : subs) {
    dvec c(y.size() + 1, 0.0);
    for (size_t i = 0; i < y.size(); ++i) c[i + 1] = c[i] + np::cabs2(y[i]);
    cums.push_back(std::move(c));
  }
  const int64_t Pi = py_round_int(P);
  auto pw = [&](const dvec& c, const std::vector<size_t>& ks) {
    dvec o;
    const int64_t last = static_cast<int64_t>(c.size()) - 1;
    for (size_t k : ks) {
      const int64_t r = npround_int(t[k]);
      o.push_back((c[static_cast<size_t>(std::min(r + Pi, last))] - c[static_cast<size_t>(std::min(r, last))]) /
                  static_cast<double>(Pi));
    }
    return o;
  };
  std::vector<std::pair<size_t, size_t>> out{runs[0]};
  for (size_t ri = 1; ri < runs.size(); ++ri) {
    const auto [a, b] = runs[ri];
    const auto [a0, b0] = out.back();
    const size_t glen = a > b0 ? a - b0 : 0;
    if (glen > 0 && glen <= max_gap) {
      std::vector<size_t> gap, near;
      for (size_t k = b0; k < a; ++k) gap.push_back(k);
      for (size_t k = std::max<size_t>(a0, b0 >= 8 ? b0 - 8 : 0); k < b0; ++k) near.push_back(k);
      for (size_t k = a; k < std::min(b, a + 8); ++k) near.push_back(k);
      bool masked = true;
      for (const auto& c : cums) masked = masked && np::median(pw(c, gap)) >= np::median(pw(c, near));
      if (masked) {
        out.back().second = b;
        continue;
      }
    }
    out.emplace_back(a, b);
  }
  return out;
}

// python _frame_complete: low-SNR burst completion from the TDD frame structure (per-symbol
// full-band evidence folded over 2.5 / 5 ms within the link's active span)
std::vector<std::pair<size_t, size_t>> frame_complete(const std::vector<std::pair<size_t, size_t>>& runs,
                                                      const dvec& z1, double P, double ofs) {
  using RV = std::vector<std::pair<size_t, size_t>>;
  if (runs.size() < 2) return runs;
  const size_t K = z1.size();
  dvec e(K);
  for (size_t i = 0; i < K; ++i) e[i] = std::isfinite(z1[i]) ? std::clamp(z1[i], -3.0, 6.0) : 0.0;
  struct Best {
    double con;
    size_t F;
    dvec zt;
    size_t a0, a1;
  };
  std::optional<Best> best;
  for (double T : {2.5e-3, 5e-3}) {
    const int64_t Fi = py_round_int(T / (P / ofs));
    if (Fi <= 0) continue;
    const size_t F = static_cast<size_t>(Fi);
    if (std::fabs(static_cast<double>(F) * P / ofs - T) > 0.02 * T || 3 * F > K) continue;
    const size_t a0 = runs.front().first / F;
    const int64_t a1i = std::min<int64_t>(static_cast<int64_t>(runs.back().second / F), static_cast<int64_t>(K / F) - 1);
    if (a1i - static_cast<int64_t>(a0) + 1 < 3) continue;
    const size_t a1 = static_cast<size_t>(a1i), nf = a1 - a0 + 1;
    dvec zt(F, 0.0);
    for (size_t r = 0; r < nf; ++r)
      for (size_t j = 0; j < F; ++j) zt[j] += e[(a0 + r) * F + j];
    for (size_t j = 0; j < F; ++j) zt[j] = zt[j] / static_cast<double>(nf) * std::sqrt(static_cast<double>(nf));
    dvec srt = zt;
    std::sort(srt.begin(), srt.end());
    const size_t q = std::max<size_t>(F / 4, 1);
    const double con = np::mean(dvec(srt.end() - static_cast<int64_t>(q), srt.end()));
    if (!best || con > best->con) best = Best{con, F, zt, a0, a1};
  }
  if (!best || best->con < 6.0) return runs;
  const size_t F = best->F;
  const dvec sm = ndi::uniform_filter1d(best->zt, 5, ndi::Mode::Wrap);
  const size_t off = np::argmin(sm);
  std::vector<std::pair<size_t, size_t>> regions;
  {
    size_t i = 0;
    while (i < F) {
      if (!(best->zt[(i + off) % F] > 3.0)) {
        ++i;
        continue;
      }
      size_t j = i;
      while (j < F && best->zt[(j + off) % F] > 3.0) ++j;
      if (j - i >= 3) regions.emplace_back(i, j);
      i = j;
    }
  }
  RV add, present;
  for (size_t f = best->a0; f <= best->a1; ++f) {
    for (const auto& [ra, rb] : regions) {
      const size_t lo = f * F + off + ra, hi = std::min(f * F + off + rb, K);
      if (hi <= lo || hi - lo < 3) continue;
      const dvec seg(e.begin() + static_cast<int64_t>(lo), e.begin() + static_cast<int64_t>(hi));
      if (np::mean(seg) * std::sqrt(static_cast<double>(seg.size())) < 4.0) continue;
      present.emplace_back(lo, hi);
      const dvec smo = ndi::uniform_filter1d(seg, 4);
      std::vector<uint8_t> on(seg.size());
      std::vector<size_t> idx;
      for (size_t i = 0; i < seg.size(); ++i) {
        on[i] = smo[i] >= 1.0;
        if (on[i]) idx.push_back(i);
      }
      for (size_t i = 0; i + 1 < idx.size(); ++i)
        if (idx[i + 1] - idx[i] > 1 && idx[i + 1] - idx[i] <= 4)
          for (size_t k = idx[i]; k < idx[i + 1]; ++k) on[k] = 1;
      // the contiguous stretch with the most evidence (max of (sum, start, stop), as python)
      bool have = false;
      double bs = 0;
      size_t ba = 0, bb = 0;
      for (size_t i = 0; i < on.size();) {
        if (!on[i]) {
          ++i;
          continue;
        }
        size_t j = i;
        while (j < on.size() && on[j]) ++j;
        if (j - i >= 3) {
          const double sum = np::sum(dvec(seg.begin() + static_cast<int64_t>(i), seg.begin() + static_cast<int64_t>(j)));
          if (!have || std::make_tuple(sum, i, j) > std::make_tuple(bs, ba, bb)) bs = sum, ba = i, bb = j, have = true;
        }
        i = j;
      }
      if (!have) continue;
      while (ba < bb && seg[ba] < 1.0) ++ba;
      while (bb > ba && seg[bb - 1] < 1.0) --bb;
      if (bb - ba >= 3) add.emplace_back(lo + ba, lo + bb);
    }
  }
  RV all = runs;
  all.insert(all.end(), add.begin(), add.end());
  std::sort(all.begin(), all.end());
  RV out;
  for (const auto& [a, b] : all) {
    bool same = false;
    if (!out.empty() && a <= out.back().second + 6)
      for (const auto& [lo, hi] : present) same = same || (lo <= out.back().second && a < hi);
    if (!out.empty() && (a <= out.back().second + 1 || same)) out.back().second = std::max(b, out.back().second);
    else out.emplace_back(a, b);
  }
  return out;
}

std::optional<Link> video_link(const FastSlicer& ch, const Cand& c, double fc, const std::vector<Emission>& hops) {
  const VideoProfile& p = *c.prof;
  cdvec y = ch.slice(c.f_off, p.occ / 2 * 1.05, p.ofs);
  dvec m = cp_metric(y, p.N, p.cp).m;
  Segments sg = segment(m, c.phase, c.P);
  if (sg.runs.empty()) return std::nullopt;
  std::vector<size_t> syms = run_syms(sg.runs);
  const size_t L_used = static_cast<size_t>(py_round_int(p.occ / (p.ofs / static_cast<double>(p.N)))) + 1;
  dvec ts;
  for (size_t k : syms) ts.push_back(sg.t[k]);
  CarrierOut co = find_carrier(ch, c.f_off, p.ofs, p.N, p.cp, ts, L_used);
  double S_off = c.f_off;
  // the strongest scan centre can be MHz off the band centre and its slice then cuts the band
  // edge (D_20260515-192148): re-centre on the estimate until it stops moving
  for (int it = 0; it < 3; ++it) {
    CarrierOut c2 = find_carrier(ch, co.car, p.ofs, p.N, p.cp, ts, L_used);
    const bool moved = std::fabs(c2.car - co.car) > 1.5 * p.ofs / static_cast<double>(p.N);
    S_off = co.car;
    co = std::move(c2);
    if (!moved) break;
  }
  Link lk;
  lk.c = c;
  lk.f_off = lk.carrier = co.car;
  lk.S = std::move(co.S);
  lk.S_off = S_off;
  lk.y = ch.slice(co.car, p.occ / 2 * 1.05, p.ofs);
  lk.m = cp_metric(lk.y, p.N, p.cp).m;
  // full band + sub-bands (Wi-Fi over 3/4 of the band must not cut a burst), in-band uplink hops
  // excluded, masked gaps bridged; a run needs full-band evidence somewhere (D_20260515-192835)
  const std::vector<cdvec> subs = sub_slices(ch, co.car, p.ofs, p.occ);
  std::vector<dvec> ms{lk.m};
  for (const auto& u : subs) ms.push_back(cp_metric(u, p.N, p.cp).m);
  std::vector<std::pair<double, double>> ex;
  for (const auto& u : hops)
    if (std::fabs((u.center_mhz * 1e6 - fc) - co.car) <= p.occ / 2) ex.emplace_back(u.start_us * 1e-6 * p.ofs, u.end_us * 1e-6 * p.ofs);
  SegMulti sm = segment_multi(ms, c.phase, c.P, 1, ex);
  const size_t max_gap = static_cast<size_t>(MAX_BRIDGE_S / (c.P / p.ofs));
  std::vector<std::pair<size_t, size_t>> runs;
  for (const auto& [a, b] : bridge(sm.runs, sm.t, c.P, subs, max_gap)) {
    bool anchored = false;
    for (size_t k = a; k < b; ++k) anchored = anchored || sm.z_full[k] > 4.0;
    if (b - a >= 3 && anchored) runs.emplace_back(a, b);
  }
  runs = frame_complete(runs, sm.z1, c.P, p.ofs);
  lk.runs = bridge(runs, sm.t, c.P, subs, max_gap);
  if (lk.runs.empty()) return std::nullopt;
  lk.t = std::move(sm.t);
  lk.syms = run_syms(lk.runs);
  // DJI evidence and the clock fit use only the full-band symbols (bridged / sub-band-only ones
  // are masked or weak and dilute the slot statistics)
  lk.ev_syms = run_syms(segment(lk.m, c.phase, c.P).runs);
  if (lk.ev_syms.empty()) lk.ev_syms = lk.syms;
  return lk;
}

// One symbol for the exact DJI-ZC test: DD coherence, fftshifted spectrum, DD sum.
struct ZcSym {
  double p;
  const cdvec* X;
  cd dd;
};

struct ZcExact {
  int64_t root;
  double corr, ref, top_coh;
  int shift;
  bool ok;
};

// python _zc_exact: exact DJI-ZC correlation on the strongest symbols (best: strongest first).
ZcExact zc_exact(const std::vector<ZcSym>& best, cd acc, size_t L, size_t N) {
  const size_t h = (L - 1) / 2;
  const int64_t Li = static_cast<int64_t>(L);
  const int64_t r0 = root_from_dd(acc, L);
  dvec hist(L, 0.0);
  for (const auto& b : best) hist[static_cast<size_t>(root_from_dd(b.dd, L))] += b.p;
  dvec ext;
  ext.insert(ext.end(), hist.end() - 8, hist.end());
  ext.insert(ext.end(), hist.begin(), hist.end());
  ext.insert(ext.end(), hist.begin(), hist.begin() + 8);
  const dvec hs = box_valid(ext, 17);
  std::vector<int64_t> modes;
  const std::vector<int64_t> order = np::argsort_quick(hs);
  for (auto it = order.rbegin(); it != order.rend(); ++it) {
    const int64_t i = *it;
    bool far = true;
    for (int64_t m_ : modes)
      if (std::min(pymod(i - m_, Li), pymod(m_ - i, Li)) <= 12) far = false;
    if (far) modes.push_back(i);
    if (modes.size() == 3) break;
  }
  std::set<int64_t> cand;
  modes.push_back(r0);
  for (int64_t m_ : modes)
    for (int64_t d = -10; d <= 10; ++d) cand.insert(pymod(m_ + d, Li));
  cand.erase(0);
  // lags r_[0:25, N-24:N]
  std::vector<size_t> lags;
  for (size_t l = 0; l < 25; ++l) lags.push_back(l);
  for (size_t l = N - 24; l < N; ++l) lags.push_back(l);
  // Xfull rows normalised by the norm of the central 2h+1 bins
  const size_t R = best.size();
  std::vector<cdvec> Xf(R);
  for (size_t r = 0; r < R; ++r) {
    Xf[r] = *best[r].X;
    dvec s2(2 * h + 1);
    for (size_t k = 0; k < 2 * h + 1; ++k) {
      const cd v = Xf[r][N / 2 - h + k];
      s2[k] = np::cmul(std::conj(v), v).real();
    }
    const double nrm = std::sqrt(np::sum(s2)) + 1e-20;
    for (auto& v : Xf[r]) v = np::cdiv(v, cd(nrm, 0.0));
  }
  const std::vector<int64_t> roots(cand.begin(), cand.end());
  const double sq = std::sqrt(static_cast<double>(L - 1));
  struct Sc {
    double v = 0;
    int sh = 0;
    bool set = false;
  };
  const std::vector<Sc> sc = pmap(roots, [&](int64_t r) {
    const cdvec zf = zc_freq(r, L);
    cdvec Zc(N, cd(0.0, 0.0));
    for (size_t k = 0; k < L; ++k) Zc[N / 2 - h + k] = std::conj(np::cdiv(zf[k], cd(sq, 0.0)));
    Sc out;
    cdvec prod(N);
    dvec qmax(R);
    for (int dsh = -2; dsh <= 2; ++dsh) {
      for (size_t r_ = 0; r_ < R; ++r_) {
        const cdvec& X = Xf[r_];
        // ifftshift(roll(X, -dsh) * conj(Z)): element j <- product at (j + N/2) % N
        for (size_t j = 0; j < N; ++j) {
          const size_t a = (j + N / 2) % N;
          prod[j] = np::cmul(X[static_cast<size_t>(pymod(static_cast<int64_t>(a) + dsh, static_cast<int64_t>(N)))], Zc[a]);
        }
        const cdvec g = fft::ifft(prod);
        double mx = -1;
        for (size_t l : lags) mx = std::max(mx, np::cabs(g[l]) * static_cast<double>(N));
        qmax[r_] = mx;
      }
      const double v = np::mean(qmax) * (dsh == 0 ? 1.0 : 0.92);
      if (v > (out.set ? out.v : 0.0)) out.v = v, out.sh = dsh, out.set = true;
    }
    return out;
  });
  int64_t zc_root = -1;
  double zc_corr = 0;
  for (size_t i = 0; i < roots.size(); ++i)
    if (sc[i].set && (zc_root < 0 || sc[i].v > zc_corr)) zc_root = roots[i], zc_corr = sc[i].v;
  int shift = 0;
  for (size_t i = 0; i < roots.size(); ++i)
    if (roots[i] == zc_root) shift = sc[i].sh;
  dvec oth;
  for (size_t i = 0; i < roots.size(); ++i)
    if (sc[i].set && std::min(pymod(roots[i] - zc_root, Li), pymod(zc_root - roots[i], Li)) > 3) oth.push_back(sc[i].v);
  const double ref = oth.empty() ? 0.0 : np::median(oth);
  dvec tops;
  for (const auto& b : best) tops.push_back(b.p);
  const double top_coh = np::mean(tops);
  const bool ok = zc_corr - ref > 0.08 && zc_corr > 0.5 * std::sqrt(std::max(top_coh, 1e-6));
  return {zc_root, zc_corr, ref, top_coh, shift, ok};
}

Evidence dji_evidence(const Link& lk) {
  const VideoProfile& p = *lk.c.prof;
  const size_t N = p.N, cp = p.cp, L = p.L;
  const cdvec& y = lk.y;
  const dvec& t = lk.t;
  const auto& syms = lk.ev_syms;
  const size_t ncls = p.scs == 15 ? 7 : 14;
  Evidence ev;
  if (L) {
    ev.kind = "zc_slot_symbol";
    ev.has_ok = true;
    // _symbol_spectra(.., shift 0.0): the mix by exp(-0j) is the identity
    std::vector<cdvec> spec(syms.size());
    std::vector<uint8_t> have(syms.size(), 0);
    parallel_for(syms.size(), [&](size_t i) {
      const int64_t s = npround_int(t[syms[i]]) + static_cast<int64_t>(cp);
      if (s + static_cast<int64_t>(N) <= static_cast<int64_t>(y.size())) {
        spec[i] = sym_spectrum(y.data() + s, N);
        have[i] = 1;
      }
    });
    std::vector<cd> acc(ncls, cd(0.0, 0.0));
    dvec per(syms.size());
    std::vector<std::pair<cd, double>> dds(syms.size());
    parallel_for(syms.size(), [&](size_t i) {
      if (have[i]) dds[i] = zc_dd(spec[i], L);
    });
    for (size_t i = 0; i < syms.size(); ++i) {
      if (!have[i]) {
        per[i] = 0.0;
        continue;
      }
      const auto& [a, mg] = dds[i];
      acc[syms[i] % ncls] += a;
      per[i] = std::hypot(a.real(), a.imag()) / (mg + 1e-20);
    }
    dvec means(ncls, 0.0);
    for (size_t c_ = 0; c_ < ncls; ++c_) {
      dvec v;
      for (size_t i = 0; i < syms.size(); ++i)
        if (syms[i] % ncls == c_) v.push_back(per[i]);
      if (!v.empty()) means[c_] = np::mean(v);
    }
    const size_t cls = np::argmax(means);
    dvec others;
    size_t n_cls = 0;
    for (size_t i = 0; i < syms.size(); ++i) {
      if (syms[i] % ncls == cls) ++n_cls;
      else others.push_back(per[i]);
    }
    const double sd = np::stddev(others) + 1e-9;
    ev.z = (means[cls] - np::median(others)) / (sd / std::sqrt(static_cast<double>(std::max<size_t>(n_cls, 1))));
    // best: (p_, i) of the slot position, stable sort by -p_, first 48
    std::vector<std::pair<double, size_t>> best;
    for (size_t i = 0; i < syms.size(); ++i)
      if (have[i] && syms[i] % ncls == cls) best.emplace_back(per[i], i);
    std::stable_sort(best.begin(), best.end(), [](const auto& a, const auto& b) { return -a.first < -b.first; });
    if (best.size() > 48) best.resize(48);
    std::vector<ZcSym> bs;
    for (const auto& [p_, i] : best) bs.push_back({p_, &spec[i], dds[i].first});
    const ZcExact zx = zc_exact(bs, acc[cls], L, N);
    ev.ok = zx.ok;
    ev.zc_root = zx.ok ? zx.root : -1;
    ev.bin_shift = zx.ok ? zx.shift : 0;
    ev.zc_corr = py_round(zx.corr, 3);
    ev.zc_ref = py_round(zx.ref, 3);
    return ev;
  }
  // O4: repeated reference symbol 4 symbols later
  ev.kind = "repeated_ref_symbol_4";
  const std::set<size_t> on(syms.begin(), syms.end());
  auto vdot = [](const cd* a, const cd* b, size_t n) {
    cd s(0.0, 0.0);
    for (size_t i = 0; i < n; ++i) s += np::cmul(std::conj(a[i]), b[i]);
    return s;
  };
  auto corr = [&](size_t k, size_t lag, double& out) {
    const int64_t s = npround_int(t[k]) + static_cast<int64_t>(cp);
    const int64_t s2 = npround_int(t[k + lag]) + static_cast<int64_t>(cp);
    if (s2 + static_cast<int64_t>(N) + 3 > static_cast<int64_t>(y.size())) return false;
    const cd* a = y.data() + s;
    const double ea = vdot(a, a, N).real();
    double mx = 0;
    for (int d = -3; d <= 3; ++d) {
      const cd* b = y.data() + s2 + d;
      const cd v = vdot(a, b, N);
      const double r = std::hypot(v.real(), v.imag()) / std::sqrt(ea * vdot(b, b, N).real() + 1e-20);
      mx = d == -3 ? r : std::max(mx, r);
    }
    out = mx;
    return true;
  };
  std::vector<size_t> ks;
  for (size_t k : syms)
    if (on.count(k + 3) && on.count(k + 4) && on.count(k + 5)) ks.push_back(k);
  struct V {
    bool ok;
    double v3, v4, v5;
  };
  const std::vector<V> vs = pmap(ks, [&](size_t k) {
    V v{true, 0, 0, 0};
    v.ok = corr(k, 3, v.v3) && corr(k, 4, v.v4) && corr(k, 5, v.v5);
    return v;
  });
  dvec d4, d35;
  std::vector<size_t> kk;
  for (size_t i = 0; i < vs.size(); ++i) {
    const auto& v = vs[i];
    if (!v.ok) continue;
    kk.push_back(ks[i]);
    d4.push_back(v.v4 - (v.v3 + v.v5) / 2);
    d35.push_back(v.v3 - v.v5);
  }
  if (d4.size() < 30) return ev;
  const double sd = std::sqrt(1.5 * np::var(d35) / 2 / static_cast<double>(d4.size())) + 1e-9;
  const double z_all = np::mean(d4) / sd;
  // python: the repeat sits at two slot positions 8 apart (a comb-4 pilot pair); best of the 14
  // rotations of {r, r+8}
  const double sd1 = std::sqrt(1.5 * np::var(d35) / 2) + 1e-9;
  double zpos = 0.0;
  for (size_t r = 0; r < ncls; ++r) {
    dvec sel;
    for (size_t i = 0; i < kk.size(); ++i)
      if (kk[i] % ncls == r || kk[i] % ncls == (r + 8) % ncls) sel.push_back(d4[i]);
    if (sel.size() >= 4) {
      const double zr = np::mean(sel) / (sd1 / std::sqrt(static_cast<double>(sel.size())));
      if (zr > zpos) zpos = zr;
    }
  }
  ev.z = std::max(z_all, zpos);
  return ev;
}

struct ClockPpm {
  double ppm, sd;
};

std::optional<ClockPpm> clock_ppm(const dvec& m, const std::vector<std::pair<size_t, double>>& t_syms, double P,
                                  size_t ncls) {
  std::vector<int64_t> ks;
  dvec ts;
  const double med = np::median(m);
  for (const auto& [k, t] : t_syms) {
    const int64_t a = npround_int(t);
    const int64_t lo = std::max<int64_t>(a - 6, 0), hi = std::min<int64_t>(a + 7, static_cast<int64_t>(m.size()));
    if (hi - lo < 13) continue;
    dvec w(13), wa(13);
    for (size_t i = 0; i < 13; ++i) {
      w[i] = std::max(m[static_cast<size_t>(lo) + i] - med, 0.0);
      wa[i] = w[i] * static_cast<double>(i);
    }
    const double ws = np::sum(w);
    if (ws <= 0) continue;
    ks.push_back(static_cast<int64_t>(k));
    ts.push_back(static_cast<double>(a - 6) + np::sum(wa) / ws);
  }
  if (ks.size() < 20 || *std::max_element(ks.begin(), ks.end()) - *std::min_element(ks.begin(), ks.end()) < 50)
    return std::nullopt;
  const size_t n = ks.size();
  dmat A(n, ncls + 1, 0.0);
  for (size_t i = 0; i < n; ++i) {
    A(i, 0) = static_cast<double>(ks[i]);
    A(i, 1 + static_cast<size_t>(ks[i] % static_cast<int64_t>(ncls))) = 1;
  }
  const dvec sol = np::lstsq(A, ts).x;
  dvec r2(n), kd(n);
  dvec kf(ks.begin(), ks.end());
  const double km = np::mean(kf);
  for (size_t i = 0; i < n; ++i) {
    double s = 0;
    for (size_t j = 0; j <= ncls; ++j) s += A(i, j) * sol[j];
    const double r = ts[i] - s;
    r2[i] = r * r;
    kd[i] = (kf[i] - km) * (kf[i] - km);
  }
  const double sk = std::sqrt(np::sum(r2) / static_cast<double>(std::max<int64_t>(static_cast<int64_t>(n) - static_cast<int64_t>(ncls) - 1, 1)) /
                              np::sum(kd));
  return ClockPpm{(sol[0] / P - 1) * 1e6, sk / P * 1e6};
}

// python _local_bursts: symbol starts found one by one (CP-metric peaks ~P apart), grouped into
// bursts of >= 3 symbols, for a link whose symbol grid moves between bursts.
std::vector<std::vector<int64_t>> local_bursts(const dvec& m, double P) {
  const double thr = 4 * np::median(m);
  const int64_t Pi = py_round_int(P), n = static_cast<int64_t>(m.size());
  std::vector<int64_t> starts;
  int64_t i = 0;
  while (i + Pi < n) {
    int64_t k = 0;
    for (int64_t j = 1; j < Pi; ++j)
      if (m[static_cast<size_t>(i + j)] > m[static_cast<size_t>(i + k)]) k = j;
    if (m[static_cast<size_t>(i + k)] > thr) {
      starts.push_back(i + k);
      i += k + static_cast<int64_t>(0.9 * P);
    } else {
      i += Pi;
    }
  }
  std::vector<std::vector<int64_t>> groups;
  if (starts.empty()) return groups;
  std::vector<int64_t> cur{starts[0]};
  for (size_t j = 1; j < starts.size(); ++j) {
    const double d = static_cast<double>(starts[j] - cur.back());
    const int64_t nn = py_round_int(d / P);
    if (nn >= 1 && nn <= 4 && std::fabs(d - static_cast<double>(nn) * P) <= 24) {
      cur.push_back(starts[j]);
    } else {
      groups.push_back(cur);
      cur = {starts[j]};
    }
  }
  groups.push_back(cur);
  // a burst's first symbol (the ZC symbol of the 2484 MHz link) can stand up to P + 256 off
  std::vector<std::vector<int64_t>> merged;
  for (auto& g : groups) {
    if (!merged.empty() && merged.back().size() == 1) {
      const double d = static_cast<double>(g[0] - merged.back()[0]);
      if (P - 12 <= d && d <= P + 256) {
        merged.back().insert(merged.back().end(), g.begin(), g.end());
        continue;
      }
    }
    merged.push_back(std::move(g));
  }
  std::vector<std::vector<int64_t>> out;
  for (auto& g : merged)
    if (g.size() >= 3) out.push_back(std::move(g));
  return out;
}

// python _retimed_link: the link re-timed burst by burst, carrier re-estimated on those symbols.
std::optional<Link> retimed_link(const FastSlicer& ch, const Link& lk) {
  const VideoProfile& p = *lk.c.prof;
  const auto groups = local_bursts(lk.m, lk.c.P);
  if (groups.empty()) return std::nullopt;
  Link rt;
  rt.c = lk.c;
  rt.ev = lk.ev;
  size_t a = 0;
  for (const auto& g : groups) {
    for (int64_t s : g) rt.t.push_back(static_cast<double>(s));
    rt.runs.emplace_back(a, a + g.size());
    a += g.size();
  }
  for (size_t k = 0; k < rt.t.size(); ++k) rt.syms.push_back(k);
  rt.ev_syms = rt.syms;
  const size_t L_used = static_cast<size_t>(py_round_int(p.occ / (p.ofs / static_cast<double>(p.N)))) + 1;
  CarrierOut co = find_carrier(ch, lk.carrier, p.ofs, p.N, p.cp, rt.t, L_used);
  rt.f_off = rt.carrier = co.car;
  rt.S = std::move(co.S);
  rt.S_off = lk.carrier;
  rt.y = ch.slice(co.car, p.occ / 2 * 1.05, p.ofs);
  rt.m = cp_metric(rt.y, p.N, p.cp).m;
  rt.retimed = true;
  return rt;
}

// python _sparse_zc_evidence: symbols whose DD coherence stands > 8 sigma above the link's own,
// then the exact-sequence test on them (no slot position across re-timed bursts).
Evidence sparse_zc_evidence(const Link& lk) {
  const VideoProfile& p = *lk.c.prof;
  const size_t N = p.N, cp = p.cp, L = p.L;
  Evidence ev;
  ev.kind = "zc_any_symbol";
  ev.has_ok = true;
  ev.ok = false;
  std::vector<cdvec> spec;
  for (double tt : lk.t) {
    const int64_t s = npround_int(tt) + static_cast<int64_t>(cp);
    if (s + static_cast<int64_t>(N) <= static_cast<int64_t>(lk.y.size())) spec.push_back(sym_spectrum(lk.y.data() + s, N));
  }
  if (spec.size() < 20) return ev;
  std::vector<std::pair<cd, double>> dd(spec.size());
  parallel_for(spec.size(), [&](size_t i) { dd[i] = zc_dd(spec[i], L); });
  dvec per(spec.size());
  for (size_t i = 0; i < spec.size(); ++i) per[i] = std::hypot(dd[i].first.real(), dd[i].first.imag()) / (dd[i].second + 1e-20);
  const double med = np::median(per);
  dvec dev(per.size());
  for (size_t i = 0; i < per.size(); ++i) dev[i] = std::fabs(per[i] - med);
  const double sd = 1.4826 * np::median(dev) + 1e-9;
  const std::vector<int64_t> order = np::argsort_quick(per);
  std::vector<size_t> top;
  for (size_t j = 0; j < order.size() && j < 48; ++j) {  // argsort(per)[::-1][:48]
    const size_t i = static_cast<size_t>(order[order.size() - 1 - j]);
    if (per[i] > med + 8 * sd) top.push_back(i);
  }
  if (top.size() < 2) return ev;
  dvec pt;
  for (size_t i : top) pt.push_back(per[i]);
  const double mt = np::mean(pt);
  ev.z = (mt - med) / (sd / std::sqrt(static_cast<double>(top.size())));
  std::vector<ZcSym> bs;
  cd acc(0.0, 0.0);
  for (size_t i : top) {
    bs.push_back({per[i], &spec[i], dd[i].first});
    acc += dd[i].first;
  }
  const ZcExact zx = zc_exact(bs, acc, L, N);
  ev.ok = zx.ok;
  ev.zc_root = zx.ok ? zx.root : -1;
  ev.bin_shift = zx.ok ? zx.shift : 0;
  ev.zc_corr = py_round(zx.corr, 3);
  ev.zc_ref = py_round(zx.ref, 3);
  return ev;
}

// python _regrid: bursts of a confirmed link on another symbol grid. A link can re-time mid-capture
// (D_20260609-151042: the O4 grid jumps by 316 samples, 10.3 us, at ~16.8 ms); the capture-level fold
// keeps one grid. Further peaks of the link's own CP fold (same carrier and profile, more than a CP
// from the grids already used) are segmented on their own phase; runs that do not overlap the link
// in time are kept, as segments sharing the parent's identity.
constexpr size_t MAX_GRIDS = 2;
constexpr size_t MIN_GRID_SYMS = 20;

std::vector<Link> regrid(const FastSlicer& ch, const Link& lk, double fc, const std::vector<Emission>& hops,
                         double z_thr) {
  const double P = lk.c.P;
  const int64_t cp = static_cast<int64_t>(lk.c.prof->cp);
  const dvec f = robust_z(fold(lk.m, P));
  const int64_t nb = static_cast<int64_t>(f.size());
  std::vector<int64_t> taken{lk.c.phase};
  std::vector<std::pair<double, double>> busy;
  for (const auto& [a, b] : lk.runs) busy.emplace_back(lk.t[a], lk.t[b - 1] + P);
  std::vector<Link> out;
  for (size_t it = 0; it < MAX_GRIDS; ++it) {
    int64_t k = -1;
    for (int64_t i = 0; i < nb; ++i) {
      bool far = true;
      for (int64_t p : taken) far = far && std::min(std::llabs(i - p), nb - std::llabs(i - p)) > cp + 16;
      if (far && (k < 0 || f[static_cast<size_t>(i)] > f[static_cast<size_t>(k)])) k = i;
    }
    if (k < 0 || f[static_cast<size_t>(k)] < z_thr) break;
    taken.push_back(k);
    Cand c2 = lk.c;
    c2.phase = k;
    c2.z = f[static_cast<size_t>(k)];
    c2.f_off = lk.carrier;
    std::optional<Link> sg = video_link(ch, c2, fc, hops);
    if (!sg) continue;
    const dvec& t = sg->t;
    std::vector<std::pair<size_t, size_t>> runs;
    size_t n = 0;
    for (const auto& [a, b] : sg->runs) {
      bool hit = false;
      for (const auto& [s_, e_] : busy) hit = hit || (t[a] < e_ && t[b - 1] + P > s_);
      if (!hit) runs.emplace_back(a, b), n += b - a;
    }
    if (n < MIN_GRID_SYMS) continue;
    sg->grid_ev = dji_evidence(*sg);
    sg->runs = std::move(runs);
    sg->ev = lk.ev;
    sg->carrier = lk.carrier;
    sg->S = lk.S;
    sg->S_off = lk.S_off;
    sg->image = false;
    for (const auto& [a, b] : sg->runs) busy.emplace_back(t[a], t[b - 1] + P);
    out.push_back(std::move(*sg));
  }
  return out;
}

std::vector<Emission> detect_video(const FastSlicer& ch, double fc, const Snr& snr, double usable_half,
                                   const std::vector<Emission>& hops, const bursts::TFPower* tf = nullptr,
                                   double z_thr = 8.0, double ev_thr = 8.0) {
  std::vector<Cand> cands = scan_video(ch, usable_half, z_thr);
  std::stable_sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) { return -a.z < -b.z; });
  std::vector<Link> links;
  auto strong = [](const Link& l) { return l.ev.okay() && l.ev.z > 8.0; };
  // a link that failed the DJI test hides only another hypothesis of itself (same profile, same
  // refined carrier): hiding all same-profile candidates within half the band lost D_20260515-192148
  auto blocks = [&](const Link& l, const VideoProfile& cp, double c_f_off, double lf) {
    if (strong(l))
      return l.c.prof->scs == cp.scs && std::fabs(lf - c_f_off) < std::max(l.c.prof->occ, cp.occ) / 2;
    return l.c.prof->name == cp.name && std::fabs(l.carrier - c_f_off) < 0.5e6;
  };
  for (const Cand& c : cands) {
    if (links.size() >= 8) break;
    bool blk = false;
    for (const auto& l : links) blk = blk || blocks(l, *c.prof, c.f_off, l.f_off);
    if (blk) continue;
    std::optional<Link> lk = video_link(ch, c, fc, hops);
    if (!lk) continue;
    for (const auto& l : links) blk = blk || blocks(l, *c.prof, lk->carrier, l.carrier);
    if (blk) continue;
    lk->ev = dji_evidence(*lk);
    if (c.prof->L && !lk->ev.okay()) {
      // no ZC on the capture-level grid: the grid may move between bursts
      std::optional<Link> rt = retimed_link(ch, *lk);
      if (rt) {
        Evidence ev = sparse_zc_evidence(*rt);
        if (ev.ok && ev.z > ev_thr) {
          lk = std::move(rt);
          lk->ev = ev;
        }
      }
    }
    links.push_back(std::move(*lk));
  }
  std::vector<uint8_t> drop(links.size(), 0);
  for (size_t i = 0; i < links.size(); ++i)
    for (size_t j = 0; j < links.size(); ++j)
      if (j != i && strong(links[j]) && !strong(links[i]) && links[j].c.prof->scs == links[i].c.prof->scs &&
          std::fabs(links[j].carrier - links[i].carrier) < 9e6)
        drop[i] = 1;
  std::vector<Link> kept;
  for (size_t i = 0; i < links.size(); ++i)
    if (!drop[i]) kept.push_back(std::move(links[i]));
  links = std::move(kept);
  for (auto& lk : links) {
    lk.image = false;
    for (const auto& o : links) {
      if (&o == &lk || o.c.prof->scs != lk.c.prof->scs) continue;
      if (std::fabs(o.carrier + lk.carrier) < 0.3e6 && o.c.z > lk.c.z * 1.5) lk.image = true;
    }
  }
  // tiers (python): 1 = DJI (z > ev_thr, or an overwhelming exact-sequence match), 2 = "-like"
  // (z in [EV_LIKE, ev_thr), structure present), one per band and never next to another tier
  constexpr double EV_LIKE = 5.0;
  std::vector<int> tier(links.size(), 0);
  for (size_t i = 0; i < links.size(); ++i) {
    const Evidence& ev = links[i].ev;
    if (ev.okay() && !links[i].image && (ev.z > ev_thr || ev.zc_corr - ev.zc_ref > 0.3)) tier[i] = 1;
  }
  std::vector<size_t> order(links.size());
  for (size_t i = 0; i < order.size(); ++i) order[i] = i;
  std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return -links[a].ev.z < -links[b].ev.z; });
  for (size_t i : order) {
    const Evidence& ev = links[i].ev;
    if (tier[i] || !ev.okay() || links[i].image || !(ev.z >= EV_LIKE)) continue;
    bool near = false;
    for (size_t j = 0; j < links.size(); ++j)
      near = near || (tier[j] && links[j].c.prof->scs == links[i].c.prof->scs &&
                      std::fabs(links[j].carrier - links[i].carrier) < 9e6);
    if (!near) tier[i] = 2;
  }
  // further symbol grids of every reported link (a link re-timed mid-capture)
  const size_t n_links = links.size();
  std::vector<size_t> lid(n_links);
  for (size_t i = 0; i < n_links; ++i) lid[i] = i;
  for (size_t i = 0; i < n_links; ++i) {
    if (!tier[i] || links[i].retimed) continue;
    std::vector<Link> segs = regrid(ch, links[i], fc, hops, z_thr);
    for (auto& sg : segs) {
      sg.parent = static_cast<int64_t>(i);
      tier.push_back(tier[i]);
      lid.push_back(lid[i]);
      links.push_back(std::move(sg));
    }
  }
  // per reported link (its further grids follow it): runs after the in-band hop exclusion
  std::vector<size_t> gorder;                       // link ids in order of first appearance
  std::map<size_t, std::vector<size_t>> groups;     // link id -> indices into links
  std::vector<std::vector<std::pair<size_t, size_t>>> runs_out(links.size());
  for (size_t li = 0; li < links.size(); ++li) {
    Link& lk = links[li];
    const VideoProfile& p = *lk.c.prof;
    const double ofs = p.ofs, occ = p.occ;
    const size_t N = p.N;
    const Evidence& ev = lk.ev;
    if (!tier[li]) continue;
    if (lk.parent >= 0) lk.carrier = links[static_cast<size_t>(lk.parent)].carrier;   // its link's corrected carrier
    else if (ev.bin_shift) lk.carrier += static_cast<double>(ev.bin_shift) * ofs / static_cast<double>(N);
    const dvec& t = lk.t;
    std::vector<std::pair<size_t, size_t>> runs = lk.runs;
    if (!hops.empty()) {
      // the RC uplink is TDD with the video and has the same numerology: a video symbol
      // overlapping an in-band hop is the hop
      const double lo_ = (fc + lk.carrier - occ / 2) / 1e6, hi_ = (fc + lk.carrier + occ / 2) / 1e6;
      std::vector<std::pair<double, double>> hs;
      for (const auto& u : hops)
        if (lo_ <= u.center_mhz && u.center_mhz <= hi_) hs.emplace_back(u.start_us, u.end_us);
      const double Pus = lk.c.P / ofs * 1e6;
      std::vector<std::pair<size_t, size_t>> nr;
      for (const auto& [a, b] : runs) {
        bool have = false;
        size_t c0 = 0, c1 = 0;
        for (size_t k = a; k < b; ++k) {
          const double tk = t[k] / ofs * 1e6;
          bool hit = false;
          for (const auto& [s_, e_] : hs) hit = hit || (tk < e_ && tk + Pus > s_);
          if (!hit) {
            if (have) c1 = k + 1;
            else c0 = k, c1 = k + 1, have = true;
            continue;
          }
          if (have && c1 - c0 >= 3) nr.emplace_back(c0, c1);
          have = false;
        }
        if (have && c1 - c0 >= 3) nr.emplace_back(c0, c1);
      }
      runs = std::move(nr);
    }
    runs_out[li] = std::move(runs);
    if (!groups.count(lid[li])) gorder.push_back(lid[li]);
    groups[lid[li]].push_back(li);
  }
  std::vector<Emission> out;
  for (size_t gid : gorder) {
    const std::vector<size_t>& grp = groups[gid];
    const size_t li = grp[0];
    Link& lk = links[li];
    const VideoProfile& p = *lk.c.prof;
    const double ofs = p.ofs, occ = p.occ, P = lk.c.P;
    const size_t N = p.N, cp = p.cp, L = p.L;
    const Evidence& ev = lk.ev;
    const size_t ncls = p.scs == 15 ? 7 : 14;
    std::vector<std::pair<size_t, double>> tk;
    for (size_t k : lk.ev_syms) tk.emplace_back(k, lk.t[k]);
    // the symbol-clock fit needs one continuous grid; a re-timed link has one per burst
    const auto cppm = lk.retimed ? std::nullopt : clock_ppm(lk.m, tk, P, ncls);
    const dvec ps = ndi::uniform_filter1d(lk.S, 9);
    const double scs_hz = ofs / static_cast<double>(N);
    int64_t k0 = npround_int(static_cast<double>(N / 2) + (lk.carrier - lk.S_off) / scs_hz);
    k0 = std::clamp<int64_t>(k0, 0, static_cast<int64_t>(N) - 1);
    dvec inb;
    for (size_t k = 0; k < N; ++k)
      if (std::fabs(static_cast<double>(static_cast<int64_t>(k) - k0)) < 0.25 * occ / scs_hz) inb.push_back(ps[k]);
    const double thr = (np::median(inb) + np::percentile(ps, 5.0)) / 2;
    // dips inside the band (O4's pilot/null pattern, a faded stretch) must not end the walk:
    // D_20260609-151042, an 18 MHz O4 link, measured 7-11 MHz and was labelled 10 MHz
    std::vector<uint8_t> occm(N);
    for (size_t k = 0; k < N; ++k) occm[k] = ps[k] > thr;
    occm = ndi::binary_closing(occm, ndi::ones1(std::max<size_t>(3, static_cast<size_t>(0.6e6 / scs_hz))));
    int64_t lo = k0, hi = k0;
    while (lo > 0 && occm[static_cast<size_t>(lo - 1)]) --lo;
    while (hi < static_cast<int64_t>(N) - 1 && occm[static_cast<size_t>(hi + 1)]) ++hi;
    const int64_t h2 = static_cast<int64_t>(N / 2);
    double f_lo = lk.S_off + static_cast<double>(lo - h2) * ofs / static_cast<double>(N);
    double f_hi = lk.S_off + static_cast<double>(hi + 1 - h2) * ofs / static_cast<double>(N);
    // bursts on every grid of the link, in seconds
    bursts::Runs seeds;
    for (size_t gi : grp)
      for (const auto& [a, b] : runs_out[gi])
        seeds.emplace_back(links[gi].t[a] / ofs, (links[gi].t[b - 1] + P) / ofs);
    std::sort(seeds.begin(), seeds.end());
    const bool zc_bw = L && ev.zc_root >= 0;
    std::string band_note = "symbol-synchronous spectrum edges";
    if (tf && !seeds.empty()) {
      // 2026-10-13, D_20260609-151042 ch0: band from in-band power on vs off, bursts from in-band power
      // seeded by the CP evidence (python detect_video)
      const double span = f_hi > f_lo ? f_hi - f_lo : occ;
      if (!zc_bw) {
        if (const auto fb = bursts::fit_band(*tf, seeds, f_lo, f_hi, std::max(f_lo - 0.5 * span, -usable_half),
                                             std::min(f_hi + 0.5 * span, usable_half))) {
          f_lo = fb->lo, f_hi = fb->hi;
          band_note = "in-band power on vs off (rectangle fit)";
        }
      }
      const double hb = static_cast<double>(L ? L - 1 : 0) * ofs / static_cast<double>(N) / 2;
      const double elo = zc_bw ? lk.carrier - hb : f_lo, ehi = zc_bw ? lk.carrier + hb : f_hi;
      auto [er, einfo] = bursts::energy_bursts(*tf, seeds, elo, ehi, usable_half, P / ofs);
      seeds = std::move(er);
      if (einfo.used && !zc_bw) {
        if (const auto fb = bursts::fit_band(*tf, seeds, f_lo, f_hi, std::max(f_lo - 0.5 * (f_hi - f_lo), -usable_half),
                                             std::min(f_hi + 0.5 * (f_hi - f_lo), usable_half))) {
          f_lo = fb->lo, f_hi = fb->hi;
          band_note = "in-band power on vs off (rectangle fit)";
        }
      }
    }
    const double bw_meas = f_hi - f_lo;
    double bw = bw_meas;
    const double cen = (fc + lk.carrier) / 1e6;
    double b_lo, b_hi;
    if (zc_bw) {
      bw = static_cast<double>(L - 1) * ofs / static_cast<double>(N);
      b_lo = cen - bw / 2e6, b_hi = cen + bw / 2e6;
    } else {
      b_lo = (fc + f_lo) / 1e6, b_hi = (fc + f_hi) / 1e6;
    }
    const std::string bw_label = bw < 13e6 ? "10MHz" : (bw < 26e6 ? "20MHz" : "40MHz");
    dvec starts_ms;
    for (const auto& [a, b] : seeds) starts_ms.push_back(a * 1e3);
    const auto T = frame_period(starts_ms);
    std::vector<CpOut> csums;
    for (size_t gi : grp) csums.push_back(cp_metric(links[gi].y, N, cp));
    for (const auto& [t0s, t1s] : seeds) {
      // per-burst residual carrier from the CP phase, on the grid of the link that covers the burst
      size_t best = 0;
      double bov = -1;
      for (size_t j = 0; j < grp.size(); ++j) {
        const Link& g = links[grp[j]];
        double ov = 0;
        for (const auto& [a, b] : runs_out[grp[j]])
          ov += std::max(0.0, std::min(t1s, (g.t[b - 1] + P) / ofs) - std::max(t0s, g.t[a] / ofs));
        if (ov > bov) bov = ov, best = j;
      }
      const Link& g = links[grp[best]];
      const CpOut& cs = csums[best];
      std::vector<int64_t> idx;
      for (double v : g.t) {
        const int64_t r = npround_int(v);
        if (t0s <= v / ofs && v / ofs < t1s && r < static_cast<int64_t>(cs.cs.size())) idx.push_back(r);
      }
      const double t0 = t0s * 1e6, t1 = t1s * 1e6;
      Json f = Json::object();
      f["scs_khz"] = p.scs;
      f["symbol_us"] = py_round(P / ofs * 1e6, 3);
      f["cp_us"] = py_round(static_cast<double>(cp) / ofs * 1e6, 2);
      f["tu_us"] = py_round(static_cast<double>(N) / ofs * 1e6, 2);
      f["bw_mhz"] = py_round(bw / 1e6, 2);
      f["bw_measured_mhz"] = py_round(bw_meas / 1e6, 2);
      f["burst_us"] = py_round(t1 - t0, 1);
      f["n_symbols"] = py_round_int((t1s - t0s) / (P / ofs));
      f["link_id"] = "video" + std::to_string(gid);
      f["cp_fold_z"] = np::round(lk.c.z, 1);
      f["evidence"] = ev.kind;
      f["evidence_z"] = py_round(ev.z, 1);
      f["evidence_tier"] = tier[li] == 1 ? "dji" : "like";
      f["band_from"] = band_note;
      f = cfo_feats(fc + lk.carrier, std::move(f));
      if (!idx.empty()) {
        const cd s = sum_sel(cs.cs, idx);
        f["cfo_burst_hz"] = py_round(f["cfo_hz"].get<double>() - std::atan2(s.imag(), s.real()) / (2 * PI) * scs_hz, 1);
      }
      if (T) {
        f["frame_period_ms"] = T->from_cycle ? py_round(T->T, 3) : np::round(T->T, 3);
        f["frame_phase_us"] = py_round(np::mod(t0, T->T * 1e3), 1);
      }
      if (cppm) {
        f["clock_ppm"] = py_round(cppm->ppm, 2);
        f["clock_ppm_sd"] = py_round(cppm->sd, 2);
      }
      if (ev.zc_root >= 0) f["zc_root"] = ev.zc_root;
      Emission e;
      e.family = "dji_ocusync";
      e.label = (p.scs == 15 ? "video_" : "video_o4_") + bw_label;
      e.start_us = py_round(t0, 1);
      e.end_us = py_round(t1, 1);
      e.f_lo_mhz = zc_bw ? np::round(b_lo, 3) : py_round(b_lo, 3);
      e.f_hi_mhz = zc_bw ? np::round(b_hi, 3) : py_round(b_hi, 3);
      e.center_mhz = zc_bw ? np::round(cen, 4) : py_round((b_lo + b_hi) / 2, 4);
      const double bc = (b_lo + b_hi) / 2, bh = (b_hi - b_lo) / 2.4;
      e.snr_db = snr(t0, t1, bc - bh, bc + bh);
      e.confidence = std::min(tier[li] == 1 ? 0.99 : 0.6, 0.5 + ev.z / 60);
      e.features = std::move(f);
      out.push_back(std::move(e));
    }
  }
  return out;
}

// ============================================================== RC uplink
std::vector<Emission> detect_uplink(const FastSlicer& fch, size_t nx, double fs, double fc, const Snr& snr,
                                    double usable_half, double k_sigma = 5.0) {
  const ofdm15k::Slicer ch = [&](double off, double half, double ofs) { return fch.slice(off, half, ofs); };
  const std::vector<ofdm15k::Profile> prof = {ofdm15k::profiles()[0]};
  const ofdm15k::Scan scan = ofdm15k::ofdm15k_scan(ch, prof, 0.5e-3, k_sigma, usable_half);
  const auto refined = pmap(scan.hits, [&](const ofdm15k::Hit& h) { return ofdm15k::refine_ofdm15k(ch, h); });
  std::vector<ofdm15k::Burst> rb;
  for (const auto& b : refined)
    if (b) rb.push_back(*b);
  std::vector<ofdm15k::Burst> bs = ofdm15k::merge_ofdm15k(rb);
  std::vector<ofdm15k::Burst> tr =
      ofdm15k::track_hops(ch, bs, scan.stats, static_cast<double>(nx) / fs, 4.0, usable_half);
  bs.insert(bs.end(), tr.begin(), tr.end());
  bs = ofdm15k::merge_ofdm15k(bs);
  const auto& P0 = ofdm15k::profiles()[0];
  const double ofs = P0.ofs;
  const size_t N = P0.N, cp = P0.cp;
  const double fcm = fc / 1e6;
  std::vector<std::optional<Emission>> res = pmap(bs, [&](const ofdm15k::Burst& b) -> std::optional<Emission> {
    const double dur = (b.end_s - b.start_s) * 1e6;
    if (b.bw_hz > 2.5e6 || !(250 <= dur && dur <= 1200)) return std::nullopt;
    const cdvec y = fch.slice(b.f_off, 0.9e6, ofs);
    const int64_t s0 = static_cast<int64_t>(b.start_s * ofs), e0 = static_cast<int64_t>(b.end_s * ofs);
    const int64_t a = std::max<int64_t>(0, s0 - static_cast<int64_t>(N));
    const int64_t e = std::min<int64_t>(static_cast<int64_t>(y.size()), e0 + static_cast<int64_t>(N));
    const CpOut co = cp_metric(y.data() + a, static_cast<size_t>(std::max<int64_t>(e - a, 0)), N, cp);
    sig::PeakOpts po;
    po.height_min = 0.3;
    po.distance = static_cast<double>(static_cast<int64_t>(0.8 * static_cast<double>(N + cp)));
    const sig::Peaks pk = sig::find_peaks(co.m, po);
    struct Z {
      double cc;
      int64_t r;
      double t_us;
      cdvec X;
    };
    std::vector<Z> zcs;
    for (int64_t p : pk.peaks) {
      const int64_t s = a + p + static_cast<int64_t>(cp);
      if (s + static_cast<int64_t>(N) > static_cast<int64_t>(y.size())) continue;
      cdvec X = sym_spectrum(y.data() + s, N);
      const cdvec sub(X.begin() + (N / 2 - 36), X.begin() + (N / 2 + 37));
      const auto [cc, r] = zc_root_search(sub, UP_L);
      zcs.push_back({cc, r, static_cast<double>(a + p) / ofs * 1e6, std::move(X)});
    }
    std::stable_sort(zcs.begin(), zcs.end(), [](const Z& u, const Z& v) { return -u.cc < -v.cc; });
    const bool ok = zcs.size() >= 2 && zcs[0].cc > 0.45 && zcs[1].cc > 0.33 && zcs[0].r == zcs[1].r &&
                    std::fabs(std::fabs(zcs[0].t_us - zcs[1].t_us) - 428.6) < 30;
    if (!ok) return std::nullopt;
    int64_t kd = -1;
    double best = 0;
    for (size_t k = N / 2 - 4; k < N / 2 + 5; ++k) {
      const double v = np::cabs2(zcs[0].X[k]) + np::cabs2(zcs[1].X[k]);
      if (kd < 0 || v < best) kd = static_cast<int64_t>(k), best = v;
    }
    const int64_t k_dc = kd - static_cast<int64_t>(N / 2);
    const double thr = np::percentile(co.m, 85.0);
    cdvec sel;
    for (size_t i = 0; i < co.m.size(); ++i)
      if (co.m[i] > thr) sel.push_back(co.cs[i]);
    const cd ss = np::sum(sel);
    const double frac = -std::atan2(ss.imag(), ss.real()) / (2 * PI) * (ofs / static_cast<double>(N));
    const double car = fc + b.f_off + static_cast<double>(k_dc) * ofs / static_cast<double>(N) + frac;
    const double cen = car / 1e6;
    const double t0 = b.start_s * 1e6, t1 = b.end_s * 1e6;
    Json f = Json::object();
    f["scs_khz"] = 15;
    f["symbol_us"] = b.period_us;
    f["cp_us"] = py_round(static_cast<double>(cp) / ofs * 1e6, 2);
    f["tu_us"] = 66.67;
    f["bw_mhz"] = py_round(b.bw_hz / 1e6, 2);
    f["burst_us"] = py_round(dur, 1);
    f["n_symbols"] = b.nsym;
    f["zc_root"] = zcs[0].r;
    f["zc_corr"] = py_round(zcs[0].cc, 2);
    f["tracked"] = b.tracked;
    f = cfo_feats(car, std::move(f));
    Emission em;
    em.family = "dji_ocusync";
    em.label = "uplink";
    em.start_us = py_round(t0, 1);
    em.end_us = py_round(t1, 1);
    em.f_lo_mhz = np::round(cen - b.bw_hz / 2e6, 3);
    em.f_hi_mhz = np::round(cen + b.bw_hz / 2e6, 3);
    em.center_mhz = np::round(cen, 4);
    em.snr_db = snr(t0, t1, cen - 0.5, cen + 0.5);
    em.confidence = 0.8;
    em.features = std::move(f);
    return em;
  });
  std::vector<Emission> out;
  for (auto& e : res)
    if (e) out.push_back(std::move(*e));
  auto sv = [](const Emission& e) { return e.snr_db && *e.snr_db != 0 ? *e.snr_db : -99.0; };
  std::vector<Emission> kept;
  for (size_t i = 0; i < out.size(); ++i) {
    bool img = false;
    for (size_t j = 0; j < out.size(); ++j)
      if (j != i && std::fabs(out[j].center_mhz + out[i].center_mhz - 2 * fcm) < 0.15 &&
          std::fabs(out[j].start_us - out[i].start_us) < 50 && sv(out[j]) > sv(out[i]) + 10)
        img = true;
    if (!img) kept.push_back(out[i]);
  }
  if (kept.size() >= 3) {
    dvec st;
    for (const auto& e : kept) st.push_back(e.start_us / 1e3);
    if (const auto T = frame_period(st)) {
      for (auto& e : kept) {
        e.features["frame_period_ms"] = T->from_cycle ? py_round(T->T, 3) : np::round(T->T, 3);
        e.features["frame_phase_us"] = py_round(np::mod(e.start_us, T->T * 1e3), 1);
      }
    }
  }
  return kept;
}

// ============================================================== DroneID
struct DidZc {
  double z600 = 0, z147 = 0;
  int64_t kb = 0;
};

DidZc droneid_zc(const cdvec& y, int64_t start) {
  size_t offs[9];
  offs[0] = 0;
  for (int i = 1; i < 9; ++i) offs[i] = offs[i - 1] + DID_N + DID_CP[i - 1];
  std::vector<cdvec> spectra;
  for (int i : {3, 5}) {
    const int64_t s = start + static_cast<int64_t>(offs[i]) + static_cast<int64_t>(DID_CP[i]);
    if (s + static_cast<int64_t>(DID_N) > static_cast<int64_t>(y.size()) || s < 0) return {};
    spectra.push_back(sym_spectrum(y.data() + s, DID_N));
  }
  const size_t L = 601;
  dvec Pz(DID_N);
  for (size_t k = 0; k < DID_N; ++k) Pz[k] = np::cabs2(spectra[0][k]) + np::cabs2(spectra[1][k]);
  const dvec w = box_valid(Pz, L);
  // DC bin from the power summed over all the frame's symbols (every symbol leaves DC empty; PLAN D4)
  dvec P2(DID_N, 0.0);
  for (int i = 0; i < 9; ++i) {
    const int64_t s = start + static_cast<int64_t>(offs[i]) + static_cast<int64_t>(DID_CP[i]);
    if (s < 0 || s + static_cast<int64_t>(DID_N) > static_cast<int64_t>(y.size())) continue;
    const cdvec X = sym_spectrum(y.data() + s, DID_N);
    for (size_t k = 0; k < DID_N; ++k) P2[k] += np::cabs2(X[k]);
  }
  const size_t b0 = np::argmax(w);
  DidZc r;
  const int64_t rts[2] = {600, 147};
  for (int i = 0; i < 2; ++i) {
    const cdvec sub(spectra[i].begin() + b0, spectra[i].begin() + b0 + L);
    const double cc = zc_root_search(sub, L, {rts[i]}).first;
    (i == 0 ? r.z600 : r.z147) = cc;
  }
  const int64_t c0 = static_cast<int64_t>(b0 + L / 2);
  int64_t dc = -1;
  double best = 0;
  for (int64_t k = std::max<int64_t>(c0 - 4, 0); k < std::min<int64_t>(c0 + 5, DID_N); ++k)
    if (dc < 0 || P2[static_cast<size_t>(k)] < best) dc = k, best = P2[static_cast<size_t>(k)];
  r.kb = dc - static_cast<int64_t>(DID_N / 2);
  return r;
}

std::vector<Emission> detect_droneid(const FastSlicer& ch, double fc, const Snr& snr, double usable_half,
                                     bool decode, double cp_thr = 0.2, double step = 3e6) {
  const double ofs = 15.36e6;
  const dvec offs = np::arange(-usable_half + 4.5e6, usable_half - 4.5e6 + 1, step);
  using C = std::tuple<double, double, double, int64_t, double, double, int64_t>;
  const std::vector<std::vector<C>> per = pmap(offs, [&](double off) {
    std::vector<C> cs;
    const cdvec y = ch.slice(off, 5.0e6, ofs);
    const dvec sc = droneid_cp_score(y);
    if (sc.empty()) return cs;
    sig::PeakOpts po;
    po.height_min = cp_thr;
    po.distance = static_cast<double>(DID_LEN / 2);
    for (int64_t p : sig::find_peaks(sc, po).peaks) {
      // the CP score is ambiguous by one symbol: the ZC pair is tested at p and one symbol either side (PLAN D4)
      DidZc z;
      int64_t q = p;
      bool first = true;
      for (int64_t qq : {p, p + static_cast<int64_t>(DID_N + DID_CP[0]), p - static_cast<int64_t>(DID_N + DID_CP[1])}) {
        const DidZc zz = droneid_zc(y, qq);
        if (first || zz.z600 + zz.z147 > z.z600 + z.z147) z = zz, q = qq;
        first = false;
      }
      if (z.z600 > 0.25 && z.z147 > 0.2)
        cs.emplace_back(z.z600 + z.z147, sc[static_cast<size_t>(p)], off, q, z.z600, z.z147, z.kb);
    }
    return cs;
  });
  std::vector<C> cands;
  for (const auto& v : per) cands.insert(cands.end(), v.begin(), v.end());
  std::stable_sort(cands.begin(), cands.end(), [](const C& a, const C& b) { return b < a; });
  std::vector<Emission> out;
  dvec seen;
  size_t offs_s[9];
  offs_s[0] = 0;
  for (int i = 1; i < 9; ++i) offs_s[i] = offs_s[i - 1] + DID_N + DID_CP[i - 1];
  for (const auto& [sum_, score, off, p, z600, z147, kb] : cands) {
    const double t0 = static_cast<double>(p) / ofs;
    bool dup = false;
    for (double s : seen) dup = dup || std::fabs(t0 - s) < 0.5e-3;
    if (dup) continue;
    const cdvec y = ch.slice(off, 5.0e6, ofs);
    seen.push_back(t0);
    // pp = y[N:] * conj(y[:-N]) (elided: conj(y[:-N]) * y[N:])
    const bool swap = y.size() - DID_N >= ELIDE_CPLX;
    cd cs(0.0, 0.0);
    for (int i = 1; i < 9; ++i) {
      cdvec seg(DID_CP[i]);
      for (size_t j = 0; j < DID_CP[i]; ++j) {
        const size_t q = static_cast<size_t>(p) + offs_s[i] + j;
        seg[j] = swap ? np::cmul(std::conj(y[q]), y[q + DID_N]) : np::cmul(y[q + DID_N], std::conj(y[q]));
      }
      cs += np::sum(seg);
    }
    const double frac = std::atan2(cs.imag(), cs.real()) / (2 * PI) * (ofs / static_cast<double>(DID_N));
    double car = fc + off + static_cast<double>(kb) * ofs / static_cast<double>(DID_N) + frac;
    // re-slice on the coarse carrier (a 3 MHz-step slice can cut the 9 MHz band) and take the integer bin again
    for (int it = 0; it < 3; ++it) {  // until the band sits in the slice centre (|bin| <= 1)
      const double off2 = std::nearbyint((car - fc) / 1e3) * 1e3;
      const int64_t kb2 = droneid_zc(ch.slice(off2, 5.0e6, ofs), p).kb;
      car = fc + off2 + static_cast<double>(kb2) * ofs / static_cast<double>(DID_N) + frac;
      if (std::llabs(kb2) <= 1) break;
    }
    // 9 symbols (OcuSync 2/3) or 8 (legacy: the first symbol dropped): an OFDM symbol in the first slot?
    double c0 = 0.0;
    {
      const int64_t cp = static_cast<int64_t>(DID_CP[0]), N = static_cast<int64_t>(DID_N);
      if (p >= 0 && p + cp + N <= static_cast<int64_t>(y.size())) {
        cd acc(0, 0);
        double ea = 0, eb = 0;
        dvec va(static_cast<size_t>(cp)), vb(static_cast<size_t>(cp));
        cdvec pr(static_cast<size_t>(cp));
        for (int64_t j = 0; j < cp; ++j) {
          const cd a = y[static_cast<size_t>(p + j)], b = y[static_cast<size_t>(p + N + j)];
          pr[static_cast<size_t>(j)] = b * std::conj(a);
          va[static_cast<size_t>(j)] = np::cabs2(a);
          vb[static_cast<size_t>(j)] = np::cabs2(b);
        }
        acc = np::sum(pr);
        ea = np::sum(va), eb = np::sum(vb);
        c0 = std::abs(acc) / std::sqrt(ea * eb + 1e-30);
      }
    }
    const int nsym = c0 >= std::max(0.3, 0.5 * score) ? 9 : 8;
    const double t0f = nsym == 8 ? static_cast<double>(p + static_cast<int64_t>(DID_N + DID_CP[0])) / ofs : t0;
    const double t1 = t0f + static_cast<double>(nsym == 9 ? DID_LEN : DID_LEN - DID_N - DID_CP[1]) / ofs;
    const double cen = car / 1e6;
    Json f = Json::object();
    f["scs_khz"] = 15;
    f["symbol_us"] = py_round(static_cast<double>(DID_N + 72) / ofs * 1e6, 2);
    f["cp_us"] = py_round(72 / ofs * 1e6, 2);
    f["tu_us"] = 66.67;
    f["bw_mhz"] = 9.0;
    f["burst_us"] = py_round((t1 - t0f) * 1e6, 1);
    f["n_symbols"] = nsym;
    f["zc_root"] = Json::array({600, 147});
    f["zc_corr"] = Json::array({py_round(z600, 2), py_round(z147, 2)});
    f["cp_score"] = py_round(score, 2);
    f["droneid_frame"] = nsym == 9 ? "9-symbol (OcuSync 2/3 generation: Mini 2, Air 2 and later)"
                                   : "8-symbol (first symbol absent: Mavic Pro / Mavic 2 family, and measured on the DJI Mini 4 Pro (O4); payload not decoded)";
    f["symbol0_cp"] = py_round(c0, 2);
    f = cfo_feats(car, std::move(f));
    // The payload decoder (droneid-decoder-reference) is not ported (no licence, no
    // verification vector): status as Python reports it when that decoder is absent.
    if (decode && nsym == 9) f["decode_status"] = "import_failed: No module named 'droneid_api'";
    Emission e;
    e.family = "dji_droneid";
    e.label = "droneid";
    e.start_us = py_round(t0f * 1e6, 1);
    e.end_us = py_round(t1 * 1e6, 1);
    e.f_lo_mhz = np::round(cen - 4.5, 3);
    e.f_hi_mhz = np::round(cen + 4.5, 3);
    e.center_mhz = np::round(cen, 4);
    e.snr_db = snr(t0f * 1e6, t1 * 1e6, cen - 4, cen + 4);
    e.confidence = 0.85;
    e.features = std::move(f);
    out.push_back(std::move(e));
  }
  return out;
}

}  // namespace

// ============================================================== entry point
std::vector<Emission> detect(const cfvec& x, double fs, double fc, bool decode_droneid) {
  // notch CW tones / LO leakage first (a tone correlates at every lag)
  std::vector<std::pair<double, double>> notch;
  for (const auto& s : find_spurs(x, fs, fc)) notch.emplace_back(s.offset_hz, 30e3);
  cdvec xd(x.size());
  for (size_t i = 0; i < x.size(); ++i) xd[i] = cd(x[i].real(), x[i].imag());
  const SpectralChannelizer chz(xd, fs, notch);
  cdvec().swap(xd);
  const FastSlicer ch(chz.X(), fs);
  const double usable_half = std::min(fs, 56e6) / 2;
  const Snr snr(x, fs, fc);
  std::vector<Emission> up = detect_uplink(ch, x.size(), fs, fc, snr, usable_half);
  const bursts::TFPower tf(x, fs);
  std::vector<Emission> video = detect_video(ch, fc, snr, usable_half, up, &tf);
  std::vector<Emission> em = std::move(video);
  for (auto& e : up) em.push_back(std::move(e));
  for (auto& e : detect_droneid(ch, fc, snr, usable_half, decode_droneid)) em.push_back(std::move(e));
  std::stable_sort(em.begin(), em.end(), [](const Emission& a, const Emission& b) {
    return a.start_us != b.start_us ? a.start_us < b.start_us : a.center_mhz < b.center_mhz;
  });
  return em;
}

}  // namespace si::rx::dji
