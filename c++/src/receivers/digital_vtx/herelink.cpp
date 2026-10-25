// Herelink half of the digital_vtx detector (python: receivers/digital_vtx/detector.py,
// pss_fold_search ... _occupied): blind LTE PSS fold over the band, carrier/root search,
// 5 ms grid, SSS -> PCI, PSS purity, subframe activity -> control / video / uplink bursts.
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <tuple>

#include "core/fft.hpp"
#include "core/ndimage.hpp"
#include "core/npcompat.hpp"
#include "core/parallel.hpp"
#include "core/pyfmt.hpp"
#include "core/seq3gpp.hpp"
#include "receivers/digital_vtx/common.hpp"
#include "receivers/digital_vtx/detector.hpp"

namespace si::rx::digital_vtx {
namespace detail {
namespace {

constexpr std::array<int, 3> kRoots = {25, 29, 34};
inline int nid2_of(int u) { return u == 25 ? 0 : u == 29 ? 1 : 2; }

constexpr double kPssFoldThreshold = 7.0;
constexpr double kPssPurityMin = 0.5;
constexpr double kPssInstanceQ = 6.0;
const double kPssOffsetUs = (160 + 2048 + 5 * 2192 + 144) / 30.72;  // HERELINK["pss_offset_us"]
const double kHalfFrameMs = 5.0;

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// python _corr_q(seg, tpl): q (Exp(1) in noise) and r at every lag of seg; e = sliding
// energy over the template (precomputed by the caller when shared)
void corr_q(const cd* seg, size_t n, const cdvec& tpl, double etpl, const double* e, dvec& q, cdvec& r) {
  const size_t nt = tpl.size(), nl = n - nt + 1;
  q.resize(nl);
  r.resize(nl);
  for (size_t k = 0; k < nl; ++k) {
    r[k] = dot_conj(seg + k, tpl.data(), nt);
    q[k] = (static_cast<double>(nt) * np::cabs2(r[k])) / (etpl * e[k] + 1e-30);
  }
}

dvec abs2(const cd* y, size_t n) {
  dvec a(n);
  for (size_t i = 0; i < n; ++i) a[i] = np::cabs2(y[i]);
  return a;
}

double sum_abs2(const cdvec& v) { return np::sum(abs2(v.data(), v.size())); }

// mean over blocks of k samples of |y|^2 (np.abs(y[a:a+nb*k])**2).reshape(-1, k).mean(1)
dvec block_power(const cdvec& y, size_t a, size_t nb, size_t k) {
  dvec out(nb);
  dvec tmp(k);
  for (size_t b = 0; b < nb; ++b) {
    for (size_t i = 0; i < k; ++i) tmp[i] = np::cabs2(y[a + b * k + i]);
    out[b] = np::mean(tmp);
  }
  return out;
}

// np.abs(z)**2 for complex64 (np::cabs squared), branch-free for finite z so the loop vectorises
inline float cabs2_finite(const cf& z) {
  const float re = z.real(), im = z.imag();
  if (!std::isfinite(re) || !std::isfinite(im)) return np::cabs2(z);
  const float ar = std::fabs(re), ai = std::fabs(im);
  const float l = std::max(ar, ai), s = std::min(ar, ai);
  const float r = s / (l == 0.0f ? 1.0f : l);
  const float a = std::sqrt(std::fma(r, r, 1.0f)) * l;
  return a * a;
}

struct Cand {
  double z;
  int root;
  double f_off;
  int64_t tau;
};

// python _pss_candidates(sr, thr)
std::vector<Cand> pss_candidates(const FoldSearch& sr, double thr, size_t max_cands = 4, double ridge = 400e3) {
  std::vector<std::tuple<double, int, size_t>> rows;
  for (int u : kRoots) {
    const dvec& z = sr.Z.at(u);
    for (size_t i = 0; i < z.size(); ++i)
      if (z[i] > thr) rows.emplace_back(z[i], u, i);
  }
  std::sort(rows.begin(), rows.end(), [](const auto& a, const auto& b) { return a > b; });
  std::vector<Cand> out;
  for (const auto& [z, u, i] : rows) {
    const double f = sr.offs[i];
    bool near = false;
    for (const auto& o : out)
      if (std::fabs(f - o.f_off) <= ridge) near = true;
    if (near) continue;
    out.push_back({z, u, f, sr.T.at(u)[i]});
    if (out.size() >= max_cands) break;
  }
  return out;
}

struct Inst {
  double t = 0, q = 0;
  cd r;
  bool t_np = true;
};

// python _instances(y, ofs, tpl, t_guess, win)
std::vector<std::optional<Inst>> instances(const cdvec& y, const cdvec& tpl, const dvec& guess, int64_t win) {
  const int64_t nt = static_cast<int64_t>(tpl.size());
  const double etpl = sum_abs2(tpl);
  std::vector<std::optional<Inst>> out;
  dvec q;
  cdvec r;
  for (double g : guess) {
    const int64_t a = py_round_int(g) - win;
    if (a < 0 || a + 2 * win + nt > static_cast<int64_t>(y.size())) {
      out.push_back(std::nullopt);
      continue;
    }
    const size_t n = static_cast<size_t>(2 * win + nt);
    const dvec p2 = abs2(y.data() + a, n);
    const dvec e = box_valid(p2.data(), n, static_cast<size_t>(nt));
    corr_q(y.data() + a, n, tpl, etpl, e.data(), q, r);
    const size_t i = np::argmax(q);
    double di = 0.0;
    if (i > 0 && i < q.size() - 1) {
      const double l = np::cabs(r[i - 1]), c = np::cabs(r[i]), rr = np::cabs(r[i + 1]);
      const double den = (l - 2 * c) + rr;
      di = den != 0 ? (0.5 * (l - rr)) / den : 0.0;
    }
    out.push_back(Inst{static_cast<double>(a + static_cast<int64_t>(i)) + di, q[i], r[i], true});
  }
  return out;
}

struct Sss {
  int nid1 = 0, pci = 0, parity = 0;
  double z = 0, z_second = 0;
  int n = 0;
};

// python _sss_decode(y, ofs, nid2, pss_useful_idx)
std::optional<Sss> sss_decode(const cdvec& y, double ofs, int nid2, const dvec& pss_idx) {
  const int64_t N = py_round_int(ofs / 15e3);
  const int64_t sym = py_round_int(2192.0 / 2048 * static_cast<double>(N));
  const cdvec d = seq3gpp::lte_pss_freq(kRoots[nid2]);
  const auto& bank = seq3gpp::sss_bank()[nid2];
  const auto& sc = seq3gpp::sc62();
  std::vector<dvec> acc(2, dvec(168, 0.0));
  int used = 0;
  for (size_t k = 0; k < pss_idx.size(); ++k) {
    const int64_t p = py_round_int(pss_idx[k]);
    const int64_t s = p - sym;
    if (s < 0 || p + N > static_cast<int64_t>(y.size())) continue;
    const cdvec Fp = fft::fft(y.data() + p, static_cast<size_t>(N));
    const cdvec Fs = fft::fft(y.data() + s, static_cast<size_t>(N));
    cdvec H(62), E(62);
    dvec h2(62);
    for (size_t i = 0; i < 62; ++i) {
      const size_t b = static_cast<size_t>(pymod(sc[i], N));
      H[i] = np::cdiv(Fp[b], d[i]);
      h2[i] = np::cabs2(H[i]);
    }
    const double mh = np::mean(h2);
    for (size_t i = 0; i < 62; ++i) {
      const size_t b = static_cast<size_t>(pymod(sc[i], N));
      E[i] = np::cdiv(np::cmul(Fs[b], std::conj(H[i])), cd(h2[i] + mh * 0.05, 0.0));
    }
    double scv[2][168];
    for (int sb = 0; sb < 2; ++sb)
      for (int n1 = 0; n1 < 168; ++n1) {
        const dvec& bv = bank[sb][n1];
        double a = 0.0;
        for (size_t i = 0; i < 62; ++i) a += bv[i] * E[i].real() - 0.0 * E[i].imag();
        scv[sb][n1] = a;
      }
    for (int par = 0; par < 2; ++par) {
      const int sb = (static_cast<int64_t>(k) + par) % 2 == 0 ? 0 : 1;
      for (int n1 = 0; n1 < 168; ++n1) acc[par][n1] += scv[sb][n1];
    }
    ++used;
  }
  if (used == 0) return std::nullopt;
  dvec flat;
  flat.reserve(336);
  for (int par = 0; par < 2; ++par) flat.insert(flat.end(), acc[par].begin(), acc[par].end());
  const size_t am = np::argmax(flat);
  std::sort(flat.begin(), flat.end());
  const double med = np::median(flat);
  dvec dev(flat.size());
  for (size_t i = 0; i < flat.size(); ++i) dev[i] = std::fabs(flat[i] - med);
  const double noise = 1.4826 * np::median(dev) + 1e-12;
  Sss r;
  r.parity = static_cast<int>(am / 168);
  r.nid1 = static_cast<int>(am % 168);
  r.pci = 3 * r.nid1 + nid2;
  r.z = (flat[flat.size() - 1] - med) / noise;
  r.z_second = (flat[flat.size() - 2] - med) / noise;
  r.n = used;
  return r;
}

// python _cfo_cp(y, ofs, sym_starts) -> fractional CFO (None when no symbol fits)
std::optional<double> cfo_cp(const cdvec& y, double ofs, const dvec& starts) {
  const int64_t N = py_round_int(ofs / 15e3);
  const int64_t cp = py_round_int(144.0 / 2048 * static_cast<double>(N));
  cd acc(0.0, 0.0);
  bool any = false;
  for (double sv : starts) {
    const int64_t s = py_round_int(sv);
    if (s < 0 || s + N + cp > static_cast<int64_t>(y.size())) continue;
    cd v(0.0, 0.0);
    for (int64_t i = 0; i < cp; ++i) v += np::cmul(std::conj(y[s + N + i]), y[s + i]);
    acc += v;
    any = true;
  }
  if (!any) return std::nullopt;
  return -std::atan2(acc.imag(), acc.real()) / (2 * kPi * static_cast<double>(N) / ofs);
}

// python _sym_starts(t_sf, ofs)
void sym_starts(double t_sf, double ofs, dvec& out) {
  const double u = ofs / 30.72e6;
  const double a = t_sf * ofs;
  for (int sl = 0; sl < 2; ++sl)
    for (int s = 0; s < 7; ++s)
      out.push_back(a + u * static_cast<double>(sl * 15360 + (s == 0 ? 0 : 2208 + (s - 1) * 2192)));
}

// python _pss_purity(X, fs, f_c, u, t_s) -> (purity, window SNR)
std::pair<double, double> pss_purity(const cdvec& X, double fs, double f_c, int u, const dvec& t_s) {
  double ofs = 0;
  const cdvec y = slice(X, fs, f_c, 0.96e6, 1.92e6, &ofs);
  const cdvec tpl = pss_time(u, ofs);
  const size_t nt = tpl.size();
  const dvec p2 = abs2(y.data(), y.size());
  const dvec e = box_valid(p2.data(), p2.size(), nt);
  const double n0 = np::percentile(e, 20);
  const double etpl = sum_abs2(tpl);
  dvec qs, ss, q;
  cdvec r;
  for (double t : t_s) {
    const int64_t a = py_round_int(t * ofs) - 3;
    if (a < 0 || a + static_cast<int64_t>(nt) + 6 > static_cast<int64_t>(y.size())) continue;
    corr_q(y.data() + a, nt + 6, tpl, etpl, e.data() + a, q, r);
    const size_t i = np::argmax(q);
    qs.push_back(q[i]);
    ss.push_back(std::max(e[static_cast<size_t>(a) + i] / n0 - 1, 1e-3));
  }
  if (qs.empty()) return {1.0, 0.0};
  const double s_med = np::median(ss);
  return {np::median(qs) / (static_cast<double>(nt) * s_med / (1 + s_med)), s_med};
}

// python _herelink_bw(y, ofs, t_list)
double herelink_bw(const cdvec& y, double ofs, const dvec& t_list) {
  const size_t nf = 2048;
  std::vector<dvec> per;
  for (size_t ti = 0; ti < t_list.size() && ti < 6; ++ti) {
    const int64_t a0 = trunc_i(t_list[ti] * ofs);
    const int64_t b0 = a0 + trunc_i(0.9e-3 * ofs);
    const int64_t n = static_cast<int64_t>(y.size());
    const int64_t a = std::min(std::max<int64_t>(a0, 0), n), b = std::min(std::max<int64_t>(b0, 0), n);
    const int64_t len = std::max<int64_t>(b - a, 0);
    if (len < 4096) continue;
    const size_t rows = static_cast<size_t>(len) / nf;
    dvec m(nf, 0.0);
    for (size_t r = 0; r < rows; ++r) {
      const cdvec F = fft::fft(y.data() + a + r * nf, nf);
      for (size_t i = 0; i < nf; ++i) m[i] = r == 0 ? np::cabs2(F[i]) : m[i] + np::cabs2(F[i]);
    }
    for (auto& v : m) v = v / static_cast<double>(rows);
    per.push_back(std::move(m));
  }
  if (per.empty()) return 18.0;
  dvec S(nf);
  for (size_t i = 0; i < nf; ++i) {
    double s = per[0][i];
    for (size_t k = 1; k < per.size(); ++k) s += per[k][i];
    S[i] = s / static_cast<double>(per.size());
  }
  const dvec f = fft::fftfreq(nf, 1 / ofs);
  dvec in, ring, edge;
  for (size_t i = 0; i < nf; ++i) {
    const double af = std::fabs(f[i]);
    if (af < 4e6) in.push_back(S[i]);
    if (af > 5e6 && af < 8.8e6) ring.push_back(S[i]);
    if (af > 9.6e6 && af < 12e6) edge.push_back(S[i]);
  }
  const double inner = np::mean(in), rg = np::mean(ring), ed = np::mean(edge);
  return (rg - ed) > 0.3 * (inner - ed) ? 18.0 : 9.0;
}

struct Edge {
  double t;
  bool np_;
};

// python _edges(yb, ofs, t_a, t_b, n0, dur_s, pad)
std::pair<Edge, Edge> edges(const cdvec& yb, double ofs, Edge t_a, Edge t_b, double n0, double dur_s, double pad) {
  const int64_t a = trunc_i(std::max(t_a.t - pad, 0.0) * ofs);
  const int64_t b = trunc_i(std::min(t_b.t + pad, dur_s) * ofs);
  const int64_t k = trunc_i(10e-6 * ofs);
  auto clamp = [&]() {
    Edge ea = t_a.t < 0 ? Edge{0.0, false} : t_a;
    Edge eb = dur_s < t_b.t ? Edge{dur_s, false} : t_b;
    return std::make_pair(ea, eb);
  };
  if (b - a < 4 * k) return clamp();
  const size_t nb = static_cast<size_t>((b - a) / k);
  const dvec e = block_power(yb, static_cast<size_t>(a), nb, static_cast<size_t>(k));
  std::vector<size_t> w;
  for (size_t i = 0; i < nb; ++i)
    if (e[i] > 4 * n0) w.push_back(i);
  if (w.size() < 3) return clamp();
  const double ta = static_cast<double>(a + static_cast<int64_t>(w.front()) * k) / ofs;
  const double tb = static_cast<double>(a + static_cast<int64_t>(w.back() + 1) * k) / ofs;
  return {Edge{ta, true}, dur_s < tb ? Edge{dur_s, false} : Edge{tb, true}};
}

// python _occupied(yb, ofs, t_a, t_b, f_c, fc, n0, lo, hi)
std::pair<double, double> occupied(const cdvec& yb, double ofs, double t_a, double t_b, double f_c, double fc,
                                   double n0, double lo, double hi) {
  const int64_t n = static_cast<int64_t>(yb.size());
  int64_t a = trunc_i(t_a * ofs), b = trunc_i(t_b * ofs);
  a = std::min(std::max<int64_t>(a, 0), n);
  b = std::min(std::max<int64_t>(b, 0), n);
  const size_t nf = 1024;
  const int64_t len = std::max<int64_t>(b - a, 0);
  if (len < static_cast<int64_t>(4 * nf)) return {lo, hi};
  const size_t rows = static_cast<size_t>(len) / nf;
  const dvec w = np::hanning(nf);
  dvec w2(nf);
  for (size_t i = 0; i < nf; ++i) w2[i] = w[i] * w[i];
  const double sw = np::sum(w2);
  dvec S(nf, 0.0);
  cdvec buf(nf);
  for (size_t r = 0; r < rows; ++r) {
    for (size_t i = 0; i < nf; ++i) buf[i] = np::cmul(yb[a + r * nf + i], cd(w[i], 0.0));
    const cdvec F = fft::fftshift(fft::fft(buf));
    for (size_t i = 0; i < nf; ++i) S[i] = r == 0 ? np::cabs2(F[i]) : S[i] + np::cabs2(F[i]);
  }
  for (auto& v : S) v = v / static_cast<double>(rows) / sw;
  S = ndi::uniform_filter1d(S, 9);
  std::vector<size_t> m;
  for (size_t i = 0; i < nf; ++i)
    if (S[i] > 2 * n0) m.push_back(i);
  if (m.size() < 5) return {lo, hi};
  auto f = [&](size_t i) { return (static_cast<double>(static_cast<int64_t>(i) - 512) / 1024.0) * ofs; };
  return {(fc + f_c + f(m.front())) / 1e6, (fc + f_c + f(m.back())) / 1e6};
}

// python _herelink_refine(x, fs, fc, X, cand, sr) -> (emissions, accepted f_off)
std::pair<std::vector<Emission>, std::optional<double>> herelink_refine(size_t n_samples, double fs, double fc,
                                                                        const cdvec& X, const Cand& cand,
                                                                        const FoldSearch& sr) {
  std::vector<Emission> out;
  const double ofs1 = sr.ofs, P = sr.P;
  const double dur_s = static_cast<double>(n_samples) / fs;
  const int64_t nhf = trunc_i(dur_s / 5e-3) + 2;
  dvec t_inst;
  for (int64_t k = 0; k < nhf; ++k) {
    const double t = (static_cast<double>(cand.tau) + static_cast<double>(k) * P) / ofs1;
    if (t + 70e-6 < dur_s) t_inst.push_back(t);
  }
  // ---- exact carrier and root, then the 30.72 MSps working slice
  const CarrierHit ch = carrier_search(X, fs, cand.f_off, t_inst);
  double f_c = ch.f_c;
  const int u = ch.root;
  double ofs = 0;
  const cdvec y = slice(X, fs, f_c, 9.3e6, 30.72e6, &ofs);
  const cdvec tpl = pss_time(u, ofs);
  dvec guess;
  for (double t : t_inst) guess.push_back(t * ofs);
  const auto inst = instances(y, tpl, guess, trunc_i(45e-6 * ofs));
  std::vector<std::pair<int64_t, Inst>> ok, good;
  for (size_t k = 0; k < inst.size(); ++k)
    if (inst[k]) ok.emplace_back(static_cast<int64_t>(k), *inst[k]);
  if (ok.size() < 2) return {out, std::nullopt};
  const double Pn = 5e-3 * ofs;
  dvec ph;
  for (const auto& [k, d] : ok)
    if (d.q > kPssInstanceQ) ph.push_back(d.t - static_cast<double>(k) * Pn);
  if (ph.empty())
    for (const auto& [k, d] : ok) ph.push_back(d.t - static_cast<double>(k) * Pn);
  const double phase = np::median(ph);
  for (const auto& [k, d] : ok)
    if (std::fabs(d.t - (phase + static_cast<double>(k) * Pn)) < 0.6e-6 * ofs && d.q > kPssInstanceQ)
      good.emplace_back(k, d);
  Json feats = Json::object();
  double Tp, t0;
  const bool t_np = good.size() >= 3;  // polyfit -> numpy floats (round() semantics downstream)
  if (good.size() >= 3) {
    dvec kk, tt;
    for (const auto& [k, d] : good) {
      kk.push_back(static_cast<double>(k));
      tt.push_back(d.t / ofs);
    }
    const dvec c = np::polyfit(kk, tt, 1);
    Tp = c[0];
    t0 = c[1];
    dvec r2(kk.size());
    for (size_t i = 0; i < kk.size(); ++i) {
      const double r = tt[i] - (t0 + Tp * kk[i]);
      r2[i] = r * r;
    }
    feats["timing_rms_us"] = py_round(std::sqrt(np::mean(r2)) * 1e6, 3);
    feats["frame_period_ms"] = np::round(Tp * 1e3, 6);
    feats["clock_ppm"] = np::round((kHalfFrameMs * 1e-3 / Tp - 1) * 1e6, 2);
  } else {
    Tp = 5e-3;
    t0 = phase / ofs;
  }
  dvec pss_idx;
  for (int64_t kk = 0; kk < nhf + 1; ++kk) {
    const double p = (t0 + static_cast<double>(kk) * Tp) * ofs;
    if (2192 < p && p + 2048 < static_cast<double>(y.size())) pss_idx.push_back(p);
  }
  // ---- SSS -> cell id and 10 ms frame phase
  const auto sss = sss_decode(y, ofs, nid2_of(u), pss_idx);
  const bool sss_ok = sss && sss->z > 8 && sss->z - sss->z_second > 3;
  // ---- purity: a real PSS owns its 1.08 MHz, another 15 kHz system's ZC sync does not
  dvec ts;
  for (const auto& [k, d] : (good.empty() ? ok : good)) ts.push_back(d.t / ofs);
  const auto [pur, s_med] = pss_purity(X, fs, f_c, u, ts);
  if (s_med > 0.15 && pur < kPssPurityMin) return {out, std::nullopt};
  if (!(sss_ok || good.size() >= 3)) return {out, std::nullopt};
  const double f_acc = f_c;
  const double t_sfp = t0 - kPssOffsetUs * 1e-6;
  double t_sf0;
  if (sss_ok) {
    t_sf0 = sss->parity == 0 ? t_sfp : t_sfp - 5e-3;
    feats["frame_phase_us"] = py_round(np::mod(t_sf0, 10e-3) * 1e6, 2);
    feats["pci"] = sss->pci;
    feats["nid1"] = sss->nid1;
  } else {
    t_sf0 = t_sfp;
    feats["frame_phase_us"] = py_round(np::mod(t_sf0, 5e-3) * 1e6, 2);
    feats["frame_phase_ref"] = "5 ms (SSS not decoded: subframe 0/5 ambiguous)";
  }
  feats["nid2"] = nid2_of(u);
  feats["pss_root"] = u;
  // ---- subframe grid (index j mod 5 relative to a PSS-bearing subframe)
  const double sf_len = 1e-3 * Tp / 5e-3;
  const int64_t j0 = trunc_i(std::floor(-t_sf0 / sf_len)) - 1;
  struct Sf {
    int64_t j;
    double t;
  };
  std::vector<Sf> sfs;
  const int64_t nsf = trunc_i(dur_s / sf_len) + 4;
  for (int64_t j = j0; j < j0 + nsf; ++j) {
    const double t = t_sf0 + static_cast<double>(j) * sf_len;
    if (t + sf_len > 0 && t < dur_s) sfs.push_back({pymod(j, 5), t});
  }
  dvec tbw;
  for (const auto& s : sfs)
    if ((s.j == 3 || s.j == 4) && 0 <= s.t && s.t < dur_s - 1e-3) tbw.push_back(s.t);
  const double bw = herelink_bw(y, ofs, tbw);
  const cdvec yb = slice(X, fs, f_c, bw / 2 * 1e6, 30.72e6);
  // per-subframe in-band power (trimmed mean over 10 us blocks)
  const int64_t kb10 = trunc_i(10e-6 * ofs);
  dvec pw;
  for (const auto& s : sfs) {
    const int64_t a = trunc_i(std::max(s.t, 0.0) * ofs), b = trunc_i(std::min(s.t + sf_len, dur_s) * ofs);
    if (static_cast<double>(b - a) < 0.2e-3 * ofs) {
      pw.push_back(kNaN);
      continue;
    }
    dvec e = block_power(yb, static_cast<size_t>(a), static_cast<size_t>((b - a) / kb10), static_cast<size_t>(kb10));
    std::sort(e.begin(), e.end());
    const size_t m = std::max<size_t>(1, static_cast<size_t>(trunc_i(static_cast<double>(e.size()) * 0.95)));
    pw.push_back(np::mean(e.data(), m));
  }
  std::vector<bool> full;
  for (const auto& s : sfs) full.push_back(s.t >= 0 && s.t + sf_len <= dur_s);
  dvec prof(5, kNaN);
  for (int64_t j = 0; j < 5; ++j) {
    dvec v;
    for (size_t i = 0; i < sfs.size(); ++i)
      if (sfs[i].j == j && full[i]) v.push_back(pw[i]);
    if (v.empty()) continue;
    dvec nn;
    for (double x : v)
      if (!std::isnan(x)) nn.push_back(x);
    prof[j] = nn.empty() ? kNaN : np::median(nn);
  }
  double pmin = kNaN;
  for (double p : prof)
    if (!std::isnan(p) && (std::isnan(pmin) || p < pmin)) pmin = p;
  const double dof = sf_len * bw * 1e6;
  dvec ref;
  for (size_t i = 0; i < sfs.size(); ++i)
    if ((sfs[i].j == 1 || sfs[i].j == 2) && full[i]) ref.push_back(pw[i]);
  const double n0 = ref.empty() ? pmin : np::median(ref);
  double sig = n0 / std::sqrt(dof);
  if (ref.size() >= 4) {
    const double mr = np::median(ref);
    dvec dv(ref.size());
    for (size_t i = 0; i < ref.size(); ++i) dv[i] = std::fabs(ref[i] - mr);
    sig = std::max(sig, 1.4826 * np::median(dv));
  }
  dvec zsf(pw.size());
  for (size_t i = 0; i < pw.size(); ++i) zsf[i] = (pw[i] - n0) / sig;
  // ---- fine CFO: CP phase over DL symbols
  dvec starts;
  for (const auto& s : sfs)
    if ((s.j == 0 || s.j == 3 || s.j == 4) && 0 <= s.t && s.t + sf_len < dur_s) sym_starts(s.t, ofs, starts);
  const auto cfo_frac = cfo_cp(y, ofs, starts);
  if (cfo_frac && std::fabs(*cfo_frac) < 3e3) f_c = f_c + *cfo_frac;
  const double f_abs = fc + f_c;
  const double nominal = std::nearbyint(f_abs / 100e3) * 100e3;
  feats["cfo_hz"] = np::round(f_abs - nominal, 1);
  feats["cfo_ppm"] = np::round((f_abs - nominal) / nominal * 1e6, 2);
  feats["nominal_carrier_mhz"] = nominal / 1e6;
  feats["carrier_mhz"] = np::round(f_abs / 1e6, 5);
  feats["cfo_ref"] = "nearest 100 kHz raster (Herelink channel plan not public)";
  feats["symbol_us"] = py_round(2192 / 30.72, 3);
  feats["tu_us"] = py_round(2048 / 30.72, 3);
  feats["cp_us"] = py_round(144 / 30.72, 3);
  feats["scs_khz"] = 15.0;
  feats["bw_mhz"] = bw;
  feats["numerology"] = "lte_like";
  feats["modulation"] = "OFDM (LTE Rel-8 downlink PHY)";
  feats["pss_instances"] = good.size();
  if (!feats.contains("frame_period_ms")) feats["frame_period_ms"] = kHalfFrameMs;
  // ---- emissions
  const double cen = f_abs / 1e6;
  const double lo = cen - bw / 2, hi = cen + bw / 2;
  auto snr_of = [&](double p) { return 10 * std::log10(std::max(p - n0, n0 * 1e-3) / n0); };
  const double conf_link = std::min(std::max(0.55 + 0.03 * cand.z + (sss_ok ? 0.2 : 0.0), 0.0), 0.99);
  std::map<int64_t, Inst> grid_q;
  for (const auto& [k, d] : good) grid_q[py_round_int((d.t / ofs - t0) / Tp)] = d;
  auto mk = [&](const std::string& label, double s_us, double e_us, double flo, double fhi, double cm,
                std::optional<double> snr, double conf, Json f) {
    Emission e;
    e.family = "herelink";
    e.label = label;
    e.start_us = s_us;
    e.end_us = e_us;
    e.f_lo_mhz = flo;
    e.f_hi_mhz = fhi;
    e.center_mhz = cm;
    e.snr_db = snr;
    e.confidence = conf;
    e.features = std::move(f);
    return e;
  };
  // start / end of a subframe-derived time: max(t, 0) -> int 0, min(.., dur_s) -> Python float
  auto r_start = [&](double t) { return t < 0 ? 0.0 : rnd(t * 1e6, 1, t_np); };
  auto r_end = [&](double t) { return dur_s < t ? py_round(dur_s * 1e6, 1) : rnd(t * 1e6, 1, t_np); };
  for (size_t i = 0; i < sfs.size(); ++i) {  // control: PSS/SSS-bearing subframes
    if (sfs[i].j != 0) continue;
    const double t = sfs[i].t;
    const int64_t k = py_round_int((t + kPssOffsetUs * 1e-6 - t0) / Tp);
    const auto it = grid_q.find(k);
    const Inst* d = it == grid_q.end() ? nullptr : &it->second;
    if (!d && !(std::isfinite(zsf[i]) && zsf[i] > 6)) continue;
    Json f = feats;
    f["burst_us"] = rnd(sf_len * 1e6, 1, t_np);
    f["pss_q"] = d ? Json(py_round(d->q, 1)) : Json(nullptr);
    if (sss_ok)
      f["subframe"] = (k + sss->parity) % 2 == 0 ? 0 : 5;
    else
      f["subframe"] = "0|5";
    f["content"] = "PSS+SSS+CRS (+PBCH in subframe 0), no PDSCH observed";
    if (d) f["pss_t_us"] = np::round(d->t / ofs * 1e6, 3);
    out.push_back(mk("control", r_start(t), r_end(t + sf_len), np::round(lo, 3), np::round(hi, 3),
                     np::round(cen, 4),
                     std::isfinite(pw[i]) ? std::optional<double>(py_round(snr_of(pw[i]), 1)) : std::nullopt,
                     d ? conf_link : 0.5, f));
  }
  // video = runs of active non-sync DL subframes; uplink = activity in idle positions
  struct Act {
    size_t i;
    int64_t j;
    double t;
  };
  std::vector<std::vector<Act>> runs;
  std::vector<Act> cur;
  auto idle = [](int64_t j) { return j == 1 || j == 2; };
  for (size_t i = 0; i < sfs.size(); ++i) {
    if (sfs[i].j == 0 || !(std::isfinite(zsf[i]) && zsf[i] > 6)) continue;
    const Act a{i, sfs[i].j, sfs[i].t};
    if (!cur.empty() && (a.i != cur.back().i + 1 || idle(a.j) != idle(cur.back().j))) {
      runs.push_back(cur);
      cur.clear();
    }
    cur.push_back(a);
  }
  if (!cur.empty()) runs.push_back(cur);
  for (const auto& r : runs) {
    const bool is_ul = idle(r[0].j);
    auto [ea, eb] = edges(yb, ofs, Edge{r[0].t, t_np}, Edge{r.back().t + sf_len, t_np}, n0, dur_s,
                          is_ul ? 0.15e-3 : 0.03e-3);
    const double lim_a = r[0].t - (is_ul ? 150e-6 : 12e-6);
    if (lim_a > ea.t) ea = Edge{lim_a, t_np};
    const double lim_b = r.back().t + sf_len + 3e-6;
    if (lim_b < eb.t) eb = Edge{lim_b, t_np};
    const double t_a = ea.t, t_b = eb.t;
    dvec pr;
    for (const auto& a : r) pr.push_back(pw[a.i]);
    const double p = np::mean(pr);
    double flo, fhi;
    Json f = Json::object();
    if (is_ul) {
      std::tie(flo, fhi) = occupied(yb, ofs, t_a, t_b, f_c, fc, n0, lo, hi);
      const int64_t a_ = trunc_i(t_a * ofs), b_ = trunc_i(t_b * ofs);
      const int64_t nb = b_ - a_ >= 0 ? (b_ - a_) / kb10 : -((a_ - b_ + kb10 - 1) / kb10);
      double duty = 0.0;
      if (nb > 0) {
        const dvec blk = block_power(yb, static_cast<size_t>(a_), static_cast<size_t>(nb), static_cast<size_t>(kb10));
        size_t c = 0;
        for (double v : blk) c += v > 1.5 * n0;
        duty = static_cast<double>(c) / static_cast<double>(blk.size());
      }
      if (t_b - t_a < 300e-6 || fhi - flo < 3.0 || duty < 0.6) continue;
      f["frame_period_ms"] = feats["frame_period_ms"];
      f["frame_phase_us"] = feats["frame_phase_us"];
      f["pci_of_downlink"] = feats.contains("pci") ? feats["pci"] : Json(nullptr);
      f["uplink_unverified"] = true;
      f["duty"] = py_round(duty, 2);
      f["note"] =
          "TDD-aligned emission in the air unit's silent subframes (UpPTS+U of an LTE TDD config-2 frame); "
          "presumed ground-unit uplink, PHY not characterised";
    } else {
      f = feats;
      dvec st;
      for (const auto& a : r) sym_starts(a.t, ofs, st);
      const auto cb = cfo_cp(y, ofs, st);
      if (cb && std::fabs(*cb) < 3e3)
        f["cfo_burst_hz"] = np::round(f_abs - (cfo_frac ? *cfo_frac : 0.0) + *cb - nominal, 1);
      flo = lo;
      fhi = hi;
    }
    const bool both_np = ea.np_ && eb.np_;
    f["burst_us"] = rnd((t_b - t_a) * 1e6, 1, both_np || ea.np_ || eb.np_);
    Json sub = Json::array();
    for (const auto& a : r) sub.push_back(a.j);
    f["subframes"] = sub;
    f["truncated"] = t_a <= 1e-6 || t_b >= dur_s - 1e-6;
    out.push_back(mk(is_ul ? "uplink" : "video", rnd(t_a * 1e6, 1, ea.np_), rnd(t_b * 1e6, 1, eb.np_),
                     np::round(flo, 3), np::round(fhi, 3), np::round((flo + fhi) / 2, 4), py_round(snr_of(p), 1),
                     is_ul ? 0.4 : conf_link, f));
  }
  // link solid but bursts below per-burst sensitivity: the expected video subframes
  bool any_video = false;
  for (const auto& e : out) any_video |= e.label == "video";
  if (!any_video) {
    for (const auto& s : sfs) {
      if (!(s.j == 3 && s.t >= 0 && s.t + 2 * sf_len <= dur_s)) continue;
      Json f = feats;
      f["burst_us"] = rnd(2 * sf_len * 1e6, 1, t_np);
      f["predicted"] = true;
      f["note"] = "below per-burst sensitivity; position from PSS timing + measured frame";
      out.push_back(mk("video", rnd(s.t * 1e6, 1, t_np), rnd((s.t + 2 * sf_len) * 1e6, 1, t_np), np::round(lo, 3),
                       np::round(hi, 3), np::round(cen, 4), std::nullopt, 0.3, f));
    }
  }
  return {out, f_acc};
}

}  // namespace

// ============================================================== PSS fold search
FoldSearch pss_fold_search(const cdvec& X, double fs, double usable_half, double step_hz, double period_s,
                           double clip, std::optional<std::pair<double, double>> f_range) {
  const size_t batch = 16;  // rows per transform batch (python: 48; per-row results are the same)
  const int64_t L = static_cast<int64_t>(X.size());
  const int64_t M = py_round_int(static_cast<double>(L) * 1.92e6 / fs);
  const double ofs = static_cast<double>(M) * fs / static_cast<double>(L);
  const double P = period_s * ofs;
  const int64_t Pint = static_cast<int64_t>(std::ceil(P));
  const bool integral = std::fabs(P - std::nearbyint(P)) < 1e-9;
  std::vector<int64_t> ph(static_cast<size_t>(M));
  for (int64_t n = 0; n < M; ++n) ph[n] = static_cast<int64_t>(std::floor(np::mod(static_cast<double>(n), P)));
  dvec cnt(static_cast<size_t>(Pint), 0.0);
  {
    const auto c = np::bincount(ph, static_cast<size_t>(Pint));
    cnt.assign(c.size(), 0.0);
    for (size_t i = 0; i < c.size(); ++i) cnt[i] = c[i] == 0 ? 1.0 : static_cast<double>(c[i]);
  }
  std::map<int, cfvec> Pc;
  cdvec tpl0;
  for (int u : kRoots) {
    const cdvec tpl = pss_time(u, ofs);
    if (u == kRoots[0]) tpl0 = tpl;
    const cdvec F = fft::fftshift(fft::fft(tpl, static_cast<size_t>(M)));
    cfvec c(F.size());
    for (size_t i = 0; i < F.size(); ++i)
      c[i] = cf(static_cast<float>(F[i].real()), static_cast<float>(-F[i].imag()));
    Pc[u] = std::move(c);
  }
  const size_t nt = tpl0.size();
  const double Ep = sum_abs2(tpl0);
  const int64_t sb = std::max<int64_t>(1, py_round_int(step_hz / fs * static_cast<double>(L)));
  const int64_t kb_a = trunc_i(-usable_half / fs * static_cast<double>(L));
  const int64_t kb_b = trunc_i(usable_half / fs * static_cast<double>(L)) + 1;
  std::vector<int64_t> kb;
  for (int64_t k = kb_a; k < kb_b; k += sb) {
    const double f = static_cast<double>(k) * fs / static_cast<double>(L);
    if (!f_range || (f >= f_range->first && f <= f_range->second)) kb.push_back(k);   // python f_range filter
  }
  FoldSearch sr;
  sr.ofs = ofs;
  sr.M = M;
  sr.P = P;
  sr.nt = nt;
  for (int64_t k : kb) sr.offs.push_back(static_cast<double>(k) * fs / static_cast<double>(L));
  if (kb.empty()) {
    for (int u : kRoots) sr.Z[u] = dvec(), sr.T[u] = std::vector<int64_t>();
    return sr;
  }
  const cdvec Xsh = fft::fftshift(X);
  cfvec Xs(static_cast<size_t>(L + M));
  for (int64_t i = 0; i < L + M; ++i) {
    const cd v = Xsh[static_cast<size_t>(i % L)];
    Xs[i] = cf(static_cast<float>(v.real()), static_cast<float>(v.imag()));
  }
  const int64_t lo = L / 2 - M / 2;
  const int64_t K = static_cast<int64_t>(py_floordiv(static_cast<double>(M), P));
  const size_t H = kb.size();
  std::map<int, dvec> pk, md;
  for (int u : kRoots) {
    pk[u] = dvec(H);
    md[u] = dvec(H);
    sr.T[u] = std::vector<int64_t>(H);
  }
  const float clipf = static_cast<float>(clip);
  const float invc = static_cast<float>(static_cast<double>(nt) / Ep);
  const float eps = static_cast<float>(1e-3);
  const size_t nb = (H + batch - 1) / batch;
  // batches of rows through pocketfft's multi-row (SIMD lane) transform, as scipy's
  // ifft(axis=1) on the python batch: per-row results are identical to a 1-D ifft
  parallel_for(nb, [&](size_t bi) {
    const size_t Ms = static_cast<size_t>(M);
    const size_t h0 = bi * batch, B = std::min(H, (bi + 1) * batch) - h0;
    std::vector<double> cs(Ms + nt);
    fvec Ey(Ms), q(Ms), a2(Ms), Qf;
    fmat INV(B, Ms);
    cfmat E(B, Ms), R(B, Ms);
    dvec Q, Q2;
    std::vector<const cf*> Yb(B);
    for (size_t b = 0; b < B; ++b) {
      Yb[b] = Xs.data() + pymod(lo + kb[h0 + b], L);
      std::copy(Yb[b], Yb[b] + Ms, E.row(b));
    }
    fft::fft_rows(E, false, 1);
    for (size_t b = 0; b < B; ++b) {
      const cf* r0 = E.row(b);
      for (size_t j = 0; j < Ms; ++j) a2[j] = cabs2_finite(r0[j]);
      double run = 0.0;
      for (size_t j = 0; j < Ms + nt; ++j) {
        run += static_cast<double>(a2[j < Ms ? j : j - Ms]);
        cs[j] = run;
      }
      for (size_t j = 0; j < Ms; ++j) Ey[j] = static_cast<float>(cs[nt - 1 + j] - (j ? cs[j - 1] : 0.0));
      const float add = eps * np::mean(Ey);
      float* inv = INV.row(b);
      for (size_t j = 0; j < Ms; ++j) inv[j] = invc / (Ey[j] + add);
    }
    for (int u : kRoots) {
      const cfvec& pc = Pc[u];
      for (size_t b = 0; b < B; ++b) {
        np::cmul_n(Yb[b], pc.data(), R.row(b), Ms);
      }
      fft::fft_rows(R, false, 1);
      for (size_t b = 0; b < B; ++b) {
        const size_t h = h0 + b;
        const cf* r = R.row(b);
        const float* inv = INV.row(b);
        for (size_t j = 0; j < Ms; ++j) {
          const float re = r[j].real(), im = r[j].imag();
          q[j] = std::min((re * re + im * im) * inv[j], clipf);
        }
        size_t nq;
        if (integral) {
          const int64_t Pi = static_cast<int64_t>(std::nearbyint(P));
          nq = static_cast<size_t>(Pi);
          Qf.assign(q.begin(), q.begin() + Pi);
          for (int64_t k = 1; k < K; ++k)
            for (int64_t c = 0; c < Pi; ++c) Qf[c] = Qf[c] + q[k * Pi + c];
          for (int64_t c = 0; c < M - K * Pi; ++c) Qf[c] = Qf[c] + q[K * Pi + c];
          Q.resize(nq);
          for (size_t c = 0; c < nq; ++c) Q[c] = static_cast<double>(Qf[c]) / cnt[c];
        } else {
          nq = static_cast<size_t>(Pint);
          Q.assign(nq, 0.0);
          for (size_t j = 0; j < Ms; ++j) Q[ph[j]] += static_cast<double>(q[j]);
          for (size_t c = 0; c < nq; ++c) Q[c] = Q[c] / cnt[c];
        }
        Q2.resize(nq);
        // max(Q, roll(Q, 1), roll(Q, -1)): the wrap-around ends apart so the interior vectorises
        auto max3 = [&](size_t c, size_t l, size_t r) { Q2[c] = std::max(Q[c], std::max(Q[l], Q[r])); };
        max3(0, nq - 1, nq > 1 ? 1 : 0);
        for (size_t c = 1; c + 1 < nq; ++c) Q2[c] = std::max(Q[c], std::max(Q[c - 1], Q[c + 1]));
        if (nq > 1) max3(nq - 1, nq - 2, 0);
        const size_t i = np::argmax(Q2);
        pk[u][h] = Q2[i];
        sr.T[u][h] = static_cast<int64_t>(i);
        md[u][h] = np::median_inplace(Q2.data(), nq);  // Q2 is rebuilt for the next row
      }
    }
  });
  dvec all;
  for (int u : kRoots) all.insert(all.end(), md[u].begin(), md[u].end());
  sr.global_median = np::median(all);
  for (int u : kRoots) {
    dvec z(H);
    for (size_t h = 0; h < H; ++h) z[h] = pk[u][h] / std::max(md[u][h], sr.global_median);
    sr.Z[u] = std::move(z);
  }
  return sr;
}

// ============================================================== carrier / root search
CarrierHit carrier_search(const cdvec& X, double fs, double f0, const dvec& t_inst_s) {
  const int n_sc = 28;
  double ofs = 0;
  const cdvec y = slice(X, fs, f0, 0.96e6 + n_sc * 15e3, 3.84e6, &ofs);
  std::map<int, cdvec> tpls;
  for (int u : kRoots) tpls[u] = pss_time(u, ofs);
  const size_t nt = tpls[25].size();
  const int64_t w = trunc_i(45e-6 * ofs);
  struct Win {
    size_t a, n;
    dvec e;
  };
  std::vector<Win> wins;
  for (double t : t_inst_s) {
    const int64_t c = py_round_int(t * ofs);
    const int64_t a = c - w, b = c + w + static_cast<int64_t>(nt);
    if (a < 0 || b > static_cast<int64_t>(y.size())) continue;
    Win wi{static_cast<size_t>(a), static_cast<size_t>(b - a), {}};
    const dvec p2 = abs2(y.data() + a, wi.n);
    wi.e = box_valid(p2.data(), wi.n, nt);
    wins.push_back(std::move(wi));
  }
  if (wins.empty()) return {f0, 25, 0.0};
  // python: tpls[u] * np.exp(2j*np.pi*df*np.arange(nt)/ofs) -> theta = ((2pi*df)*n)*(1/ofs)
  auto score = [&](int u, double df) {
    const cdvec& tp0 = tpls[u];
    const double c = 2 * kPi * df, scl = 1.0 / ofs;
    cdvec tp(nt);
    for (size_t n = 0; n < nt; ++n) tp[n] = np::cmul(tp0[n], expj((c * static_cast<double>(n)) * scl));
    const double et = sum_abs2(tp);
    double s = 0.0;
    dvec q;
    cdvec r;
    for (const auto& wi : wins) {
      corr_q(y.data() + wi.a, wi.n, tp, et, wi.e.data(), q, r);
      double m = q[0];
      for (double v : q) m = std::max(m, v);
      s += m;
    }
    return s;
  };
  struct Hyp {
    double df;
    int u;
  };
  std::vector<Hyp> hyps;
  for (int di = 0; di < 6; ++di) {
    const double d = -7.5e3 + di * 2.5e3;
    for (int k = -n_sc; k <= n_sc; ++k)
      for (int u : kRoots) hyps.push_back({k * 15e3 + d, u});
  }
  const dvec s = pmap(hyps, [&](const Hyp& h) { return score(h.u, h.df); });
  double best_s = -1.0, best_f = f0;
  int best_u = 25;
  for (size_t i = 0; i < hyps.size(); ++i)
    if (s[i] > best_s) {
      best_s = s[i];
      best_f = f0 + hyps[i].df;
      best_u = hyps[i].u;
    }
  // fine: +-1.25 kHz in 250 Hz steps at the winner
  const double f1 = best_f;
  const int u = best_u;
  std::vector<double> fine;
  for (int i = 0; i < 11; ++i) fine.push_back(-1.25e3 + i * 250.0);
  const dvec sf = pmap(fine, [&](double d) { return score(u, f1 + d - f0); });
  for (size_t i = 0; i < fine.size(); ++i)
    if (sf[i] > best_s) {
      best_s = sf[i];
      best_f = f1 + fine[i];
      best_u = u;
    }
  return {best_f, best_u, best_s};
}

std::vector<Emission> herelink_detect(size_t n_samples, double fs, double fc, const cdvec& X) {
  const double usable = std::min(fs, 56e6) / 2 - 0.6e6;
  // latency (PLAN P1; python HERELINK_RF): Herelink is 2.4 GHz only, so only carrier hypotheses inside
  // 2370-2520 MHz are searched (none at all elsewhere)
  const double lo = std::max(-usable, 2370e6 - fc), hi = std::min(usable, 2520e6 - fc);
  if (lo > hi) return {};
  const bool full = lo <= -usable && hi >= usable;
  const FoldSearch sr = pss_fold_search(X, fs, usable, 10e3, 5e-3, 24.0,
                                        full ? std::nullopt : std::optional<std::pair<double, double>>({lo, hi}));
  if (sr.offs.empty()) return {};
  std::vector<Emission> out;
  std::vector<std::pair<double, double>> accepted;
  for (const auto& c : pss_candidates(sr, kPssFoldThreshold)) {
    bool skip = false;
    for (const auto& [a, z] : accepted)
      if (std::fabs(c.f_off - a) < 9.5e6) skip = true;
    if (skip) continue;
    // receiver IQ imbalance mirrors a strong carrier to -f (conjugated: root 25 <-> 34)
    for (const auto& [a, z] : accepted)
      if (std::fabs(c.f_off + a) < 0.3e6 && c.z < 0.8 * z) skip = true;
    if (skip) continue;
    auto [em, f_off] = herelink_refine(n_samples, fs, fc, X, c, sr);
    if (f_off) {
      accepted.emplace_back(*f_off, c.z);
      out.insert(out.end(), em.begin(), em.end());
    }
  }
  return out;
}

}  // namespace detail
}  // namespace si::rx::digital_vtx
