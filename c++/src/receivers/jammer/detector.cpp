// Swept-carrier jammer detector (python: receivers/jammer/detector.py). dtypes as numpy: the occupancy, its
// log and cyclic spectrum are float32 (scipy single-precision FFTs); folds accumulate float32 in frame order (scipy
// sparse product) and continue in float64; the scalars are Python floats.
#include "receivers/jammer/detector.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>

#include "core/fft.hpp"
#include "core/npcompat.hpp"
#include "core/pyfmt.hpp"

namespace si::rx::jammer {

namespace {

constexpr double kTvLines[] = {15625.0, 15734.264};

// python _interp: parabolic interpolation on log y (float64) of a peak at k -> fractional index
double interp_peak(const double* y, size_t n, size_t k) {
  if (k > 0 && k + 1 < n) {
    const double a = std::log(y[k - 1] + 1e-30), b = std::log(y[k] + 1e-30), c = std::log(y[k + 1] + 1e-30);
    const double d = a - 2 * b + c;
    if (d < 0) return static_cast<double>(k) + 0.5 * (a - c) / d;
  }
  return static_cast<double>(k);
}

// python _local_median(S, w=64)
dvec local_median(const dvec& S, size_t w = 64) {
  dvec xp, fp, buf;
  for (size_t i = 0; i < S.size(); i += w) {
    const size_t a = i >= w ? i - w : 0, b = std::min(S.size(), i + w + 1);
    xp.push_back(static_cast<double>(i));
    fp.push_back(np::median_buf(S.data() + a, b - a, buf));
  }
  dvec x(S.size());
  for (size_t i = 0; i < x.size(); ++i) x[i] = static_cast<double>(i);
  dvec out = np::interp(x, xp, fp);
  for (auto& v : out) v += 1e-30;
  return out;
}

// python round(x) for a Python float (ties to even)
double pyround0(double x) { return std::nearbyint(x); }

// python candidates(): refined line rates and the strongest line contrast (dB)
std::pair<dvec, double> candidates(const fmat& P, double fr, double T) {
  const size_t m = P.rows, nb = P.cols;
  // med = np.median(P[::8], axis=0) + 1e-20 (float32)
  std::vector<float> med(nb), col;
  for (size_t b = 0; b < nb; ++b) {
    col.clear();
    for (size_t r = 0; r < m; r += 8) col.push_back(P(r, b));
    med[b] = np::median_inplace(col.data(), col.size()) + 1e-20f;
  }
  fmat Z(m, nb);
  for (size_t r = 0; r < m; ++r)
    for (size_t b = 0; b < nb; ++b) Z(r, b) = std::log(std::max(P(r, b) / med[b], 1e-3f));
  const std::vector<float> zm = np::mean_axis0(Z);
  for (size_t r = 0; r < m; ++r)
    for (size_t b = 0; b < nb; ++b) Z(r, b) -= zm[b];
  const size_t n = fft::good_size_real(m);
  const size_t nr = n / 2 + 1;
  const double amin = std::max(kMinRateHz, 8.0 / T), amax = fr / 8;
  const auto lo = static_cast<int64_t>(std::ceil(amin * n / fr));
  const int64_t hi = std::min(static_cast<int64_t>(nr) - 3, static_cast<int64_t>(amax * n / fr));
  if (hi - lo < 16) return {{}, 0.0};
  // A = R.real**2 + R.imag**2 per column (float32), whitened by its mean over [lo, hi), summed over bins
  fmat A(nr, nb);
  std::vector<float> c(m);
  for (size_t b = 0; b < nb; ++b) {
    for (size_t r = 0; r < m; ++r) c[r] = Z(r, b);
    const cfvec R = fft::rfft(c.data(), m, n);
    for (size_t r = 0; r < nr; ++r) {
      const float re2 = R[r].real() * R[r].real();
      const float im2 = R[r].imag() * R[r].imag();
      A(r, b) = re2 + im2;
    }
  }
  std::vector<float> am(nb, 0.0f);
  for (int64_t r = lo; r < hi; ++r)
    for (size_t b = 0; b < nb; ++b) am[b] += A(r, b);
  for (auto& v : am) v = v / static_cast<float>(hi - lo) + 1e-30f;
  dvec S(nr);
  std::vector<float> row(nb);
  for (size_t r = 0; r < nr; ++r) {
    for (size_t b = 0; b < nb; ++b) row[b] = A(r, b) / am[b];
    S[r] = static_cast<double>(np::pairwise_sum(row.data(), nb));
  }
  const dvec lm = local_median(S);
  dvec C(nr);
  for (size_t r = 0; r < nr; ++r) C[r] = S[r] / lm[r];
  const dvec Cw(C.begin() + lo, C.begin() + hi);
  const auto order = np::argsort_quick(Cw);
  std::vector<int64_t> peaks;
  for (auto it = order.rbegin(); it != order.rend(); ++it) {
    const int64_t k = lo + *it;
    if (C[k] < std::pow(10.0, kTrigDb / 10) || peaks.size() >= kNCand) break;
    bool ok = true;
    for (int64_t j : peaks) {
      const double kj = static_cast<double>(k) / j, jk = static_cast<double>(j) / k;
      if (!(std::llabs(k - j) > 3 && std::fabs(kj - pyround0(kj)) * j > 3 && std::fabs(jk - pyround0(jk)) * k > 3)) {
        ok = false;
        break;
      }
    }
    if (ok) peaks.push_back(k);
  }
  dvec rates;
  for (int64_t k : peaks) {
    double best = interp_peak(S.data(), nr, static_cast<size_t>(k)) * fr / n;
    for (int h = 2; h < 40; ++h) {
      if (h * best > amax) break;
      const auto j = static_cast<int64_t>(pyround0(h * best * n / fr));
      if (j + 3 > static_cast<int64_t>(nr)) break;
      const int64_t jj = j - 2 + static_cast<int64_t>(np::argmax(S.data() + (j - 2), 5));
      if (C[jj] >= std::pow(10.0, kHarmDb / 10)) best = interp_peak(S.data(), nr, static_cast<size_t>(jj)) * fr / n / h;
    }
    rates.push_back(best);
  }
  return {rates, 10 * std::log10(*std::max_element(Cw.begin(), Cw.end()))};
}

// python fold(): mean occupancy per (phase, bin) over all periods; rows [r0, r1) of P, the first at time t0
dmat fold(const fmat& P, size_t r0, size_t r1, double fr, double period, size_t nph, double t0 = 0.0) {
  const size_t nb = P.cols;
  fmat acc(nph, nb, 0.0f);
  std::vector<int64_t> cnt(nph, 0);
  for (size_t i = r0; i < r1; ++i) {
    const double t = t0 + (static_cast<double>(i - r0) + 0.5) / fr;
    auto ph = static_cast<int64_t>(static_cast<int32_t>(std::fmod(t / period, 1.0) * static_cast<double>(nph)));
    ph = std::min<int64_t>(ph, static_cast<int64_t>(nph) - 1);
    ++cnt[ph];
    float* a = acc.row(ph);
    const float* p = P.row(i);
    for (size_t b = 0; b < nb; ++b) a[b] += p[b];
  }
  dmat out(nph, nb);
  for (size_t r = 0; r < nph; ++r)
    for (size_t b = 0; b < nb; ++b)
      out(r, b) = static_cast<double>(acc(r, b)) / static_cast<double>(std::max<int64_t>(cnt[r], 1));
  return out;
}

// python _norm: F / (F.mean(axis=0) + 1e-20) - 1.0
void norm(dmat& F) {
  const dvec mu = np::mean_axis0(F);
  for (size_t r = 0; r < F.rows; ++r)
    for (size_t b = 0; b < F.cols; ++b) F(r, b) = F(r, b) / (mu[b] + 1e-20) - 1.0;
}

size_t phases(double p, double fr) {
  return static_cast<size_t>(std::min(kMaxPhases, std::max(8.0, p * fr)));
}

// python _period(): the true period from a cyclic line's period, rate-corrected; false when refused
bool true_period(const fmat& P, double fr, double T, double& period, size_t& nph, Json& info) {
  const double amin = std::max(kMinRateHz, 8.0 / T);
  std::map<int64_t, std::pair<double, size_t>> cache;
  auto contrast = [&](double p) {
    const auto key = static_cast<int64_t>(pyround0(p * 1e9));
    auto it = cache.find(key);
    if (it == cache.end()) {
      const size_t n = phases(p, fr);
      dmat F = fold(P, 0, P.rows, fr, p, n);
      norm(F);
      dvec sq(F.size());
      for (size_t i = 0; i < sq.size(); ++i) sq[i] = F.data[i] * F.data[i];
      it = cache.emplace(key, std::make_pair(np::mean(sq), n)).first;
    }
    return it->second.first;
  };
  std::vector<int> ds;
  for (int d = 1; d <= kMaxDiv; ++d)
    if (1.0 / (d * period) >= amin) ds.push_back(d);
  if (ds.empty()) {
    info["reject"] = "sweep period too long for the block";
    return false;
  }
  std::map<int, double> v;
  double vmax = -1;
  for (int d : ds) vmax = std::max(vmax, v[d] = contrast(d * period));
  int d = ds.front();
  for (int dd : ds)
    if (v[dd] >= kContrastFrac * vmax) {
      d = dd;
      break;
    }
  period *= d;
  nph = cache[static_cast<int64_t>(pyround0(period * 1e9))].second;
  double longer = -1;
  bool any_longer = false;
  for (int k = 2; k < 6; ++k)
    if (1.0 / (k * period) >= amin) {
      longer = std::max(longer, contrast(k * period));
      any_longer = true;
    }
  info["contrast"] = py_round(v[d], 4);
  if (any_longer && longer * kContrastFrac > v[d]) {
    info["reject"] = "the pattern repeats only on a longer period (pulsed, not a continuous sweep)";
    return false;
  }
  // rate correction: lag of each segment's fold against the whole-block fold, fitted linearly over time
  dmat ref = fold(P, 0, P.rows, fr, period, nph);
  norm(ref);
  const auto nseg = static_cast<size_t>(std::min(kNSeg, std::max(2.0, T / period / kSegPeriods)));
  const size_t q = P.rows / nseg;
  dvec lags, mids, prod(ref.size()), cc(nph);
  for (size_t i = 0; i < nseg; ++i) {
    dmat B = fold(P, i * q, (i + 1) * q, fr, period, nph, static_cast<double>(i * q) / fr);
    norm(B);
    for (size_t s = 0; s < nph; ++s) {   // np.sum(ref * np.roll(B, s, axis=0))
      for (size_t r = 0; r < nph; ++r) {
        const double* br = B.row((r + nph - s) % nph);
        const double* rr = ref.row(r);
        for (size_t b = 0; b < ref.cols; ++b) prod[r * ref.cols + b] = rr[b] * br[b];
      }
      cc[s] = np::sum(prod);
    }
    const size_t s0 = np::argmax(cc);
    const double cmin = *std::min_element(cc.begin(), cc.end());
    const double c3[3] = {cc[(s0 + nph - 1) % nph] - cmin + 1e-9, cc[s0] - cmin + 1e-9, cc[(s0 + 1) % nph] - cmin + 1e-9};
    const double sf = static_cast<double>(s0) + interp_peak(c3, 3, 1) - 1;
    lags.push_back(sf > nph / 2.0 ? sf - static_cast<double>(nph) : sf);
    mids.push_back((i + 0.5) * static_cast<double>(q) / fr);
  }
  const double slope = np::polyfit(mids, lags, 1)[0];   // phases per second the pattern runs early
  if (std::fabs(slope) * static_cast<double>(P.rows) / fr <= nph / 4.0) period *= 1 - slope * period / nph;
  return true;
}

// python _confirm(): the sweep on the fold, or found = false
Sweep confirm(const fmat& P, double fr, double T, double line_hz, double df, Json& info) {
  Sweep sw;
  for (double tv : kTvLines)
    for (int h = 1; h < 17; ++h)
      if (std::fabs(line_hz / (h * tv) - 1) < 2e-3) {
        info["reject"] = "analog video line rate";
        return sw;
      }
  double period = 1.0 / line_hz;
  size_t nph = 0;
  if (!true_period(P, fr, T, period, nph, info)) return sw;
  info["rate_hz"] = 1.0 / period;
  const size_t nb = P.cols;
  const dmat Pf = fold(P, 0, P.rows, fr, period, nph);
  const dvec floor = np::percentile_axis0(Pf, 10.0);
  dmat E(nph, nb);
  for (size_t r = 0; r < nph; ++r)
    for (size_t b = 0; b < nb; ++b) E(r, b) = Pf(r, b) / (floor[b] + 1e-20);
  std::vector<size_t> k(nph);
  dvec pk_db(nph);
  for (size_t r = 0; r < nph; ++r) {
    k[r] = np::argmax(E.row(r), nb);
    pk_db[r] = 10 * std::log10(E(r, k[r]));
  }
  dvec sorted = pk_db;
  std::sort(sorted.begin(), sorted.end());
  const double thr_db = std::max(kActiveDb, np::median(dvec(sorted.begin() + nph / 2, sorted.end())) - kActiveRelDb);
  std::vector<size_t> act;
  std::vector<uint8_t> active(nph);
  for (size_t r = 0; r < nph; ++r)
    if ((active[r] = pk_db[r] >= thr_db)) act.push_back(r);
  const double afrac = static_cast<double>(act.size()) / static_cast<double>(nph);
  info["active_frac"] = py_round(afrac, 3);
  if (afrac < kMinActive) {
    info["reject"] = fmt("carrier in band for %.0f%% of the period (< 30%%)", afrac * 100);
    return sw;
  }
  // the trajectory, smoothed by a circular median over kSmooth phases (none below 40 phases per period): within a
  // fast transit (smeared by the sweep-to-sweep jitter) a single phase's peak can land on a clipping product
  const size_t h = std::min<size_t>(kSmooth / 2, nph / 40);
  dvec fk0(nph), fk(nph), ks(nph);
  for (size_t r = 0; r < nph; ++r) fk0[r] = interp_peak(E.row(r), nb, k[r]);
  for (size_t r = 0; r < nph; ++r) {
    dvec wf;
    std::vector<int64_t> wk;
    for (size_t o = 0; o <= 2 * h; ++o) {
      const size_t i = (r + nph + o - h) % nph;
      wf.push_back(fk0[i]);
      wk.push_back(static_cast<int64_t>(k[i]));
    }
    std::sort(wf.begin(), wf.end());
    std::sort(wk.begin(), wk.end());
    fk[r] = wf[h];                                   // odd count: the middle value
    ks[r] = static_cast<double>(wk[h]);
  }
  double fmin = fk[act[0]], fmax = fk[act[0]];
  auto lo_b = static_cast<size_t>(ks[act[0]]), hi_b = lo_b;
  for (size_t r : act) {
    fmin = std::min(fmin, fk[r]);
    fmax = std::max(fmax, fk[r]);
    lo_b = std::min(lo_b, static_cast<size_t>(ks[r]));
    hi_b = std::max(hi_b, static_cast<size_t>(ks[r]));
  }
  ++hi_b;
  const double span = (fmax - fmin) * df;
  dvec dom(nph), buf;
  for (size_t r = 0; r < nph; ++r)
    dom[r] = pk_db[r] - 10 * std::log10(np::median_buf(E.row(r) + lo_b, hi_b - lo_b, buf) + 1e-20);
  // a jump between consecutive phases is a flyback (or a hop) only when some bin between the two peaks is dark in
  // both phases (below the weaker peak - kBridgeDb); a jitter-smeared fast transit keeps them lit
  const size_t L = act.size();
  std::vector<uint8_t> gap(L);
  int flyback = 0, exits = 0;
  const double jump = std::max(4.0, 0.3 * span / df);
  for (size_t j = 0; j < L; ++j) {
    const size_t a = act[j], b = act[(j + 1) % L];
    gap[j] = ((b + nph - a) % nph) > 1;
    if (gap[j]) {
      ++exits;
      continue;
    }
    if (std::fabs(ks[b] - ks[a]) <= jump) continue;
    const auto lo_k = static_cast<size_t>(std::min(ks[a], ks[b])), hi_k = static_cast<size_t>(std::max(ks[a], ks[b]));
    const double floor_lin = std::pow(10.0, std::max(kActiveDb, std::min(pk_db[a], pk_db[b]) - kBridgeDb) / 10);
    for (size_t c = lo_k; c <= hi_k; ++c)
      if (std::max(E(a, c), E(b, c)) < floor_lin) {
        ++flyback;
        break;
      }
  }
  dvec dact;
  for (size_t r : act) dact.push_back(dom[r]);
  const double dom75 = np::percentile(dact, 75.0);
  info["span_hz"] = py_round_int(span);
  info["flyback"] = flyback;
  info["exits"] = exits;
  info["dom_db"] = py_round(dom75, 1);
  if (span < kMinSpanHz) {
    info["reject"] = fmt("sweep span %.1f MHz (< 10)", span / 1e6);
    return sw;
  }
  if (flyback > kMaxFlyback || exits > kMaxExits) {
    info["reject"] = fmt("%d frequency jumps, %d band exits per period (hopping, not a sweep)", flyback, exits);
    return sw;
  }
  if (dom75 < kMinDomDb) {
    info["reject"] = fmt("no single dominant carrier per phase (%.1f dB < 4.0)", dom75);
    return sw;
  }
  // occupied extent: from the trajectory's range outward over bins lit at some phase, at most kMaxOvershoot bins
  std::vector<uint8_t> lit(nb, 0);
  const double thr_lin = std::pow(10.0, thr_db / 10);
  for (size_t b = 0; b < nb; ++b) {
    double mx = E(0, b);
    for (size_t r = 1; r < nph; ++r) mx = std::max(mx, E(r, b));
    lit[b] = mx >= thr_lin;
  }
  const size_t lo_lim = lo_b > kMaxOvershoot ? lo_b - kMaxOvershoot : 0;
  const size_t hi_lim = std::min(nb, hi_b + kMaxOvershoot);
  while (lo_b > lo_lim && lit[lo_b - 1]) --lo_b;
  while (hi_b < hi_lim && lit[hi_b]) ++hi_b;
  // where the carrier leaves the band: the trajectory on both sides of each gap is near one edge
  std::set<std::string> beyond;
  for (size_t j = 0; j < L; ++j)
    if (gap[j]) beyond.insert(ks[act[j]] + ks[act[(j + 1) % L]] < nb - 1 ? "below" : "above");
  sw.found = true;
  sw.period = period;
  sw.nph = nph;
  sw.fk = std::move(fk);
  sw.active = std::move(active);
  sw.pk_db = std::move(pk_db);
  sw.span = span;
  sw.dom75 = dom75;
  sw.flyback = flyback;
  sw.exits = exits;
  sw.beyond.assign(beyond.begin(), beyond.end());
  sw.occupied_hz = {static_cast<double>(lo_b), static_cast<double>(hi_b)};   // bins; made absolute in measure()
  return sw;
}

}  // namespace

fmat occupancy(const cfvec& x) {
  const size_t h = kNfft / 2;
  const size_t m = x.size() / kNfft - 1;
  const dvec hd = np::hanning(kNfft);
  std::vector<float> w(kNfft);
  for (size_t i = 0; i < kNfft; ++i) w[i] = static_cast<float>(hd[i]);
  cfmat F(2 * m, kNfft);
  for (size_t r = 0; r < m; ++r)
    for (size_t i = 0; i < kNfft; ++i) {
      F(2 * r, i) = cf(x[r * kNfft + i].real() * w[i], x[r * kNfft + i].imag() * w[i]);
      F(2 * r + 1, i) = cf(x[h + r * kNfft + i].real() * w[i], x[h + r * kNfft + i].imag() * w[i]);
    }
  fft::fft_rows(F, true, 1);
  fmat P(m, kNfft);
  for (size_t r = 0; r < m; ++r)
    for (size_t i = 0; i < kNfft; ++i) {
      const cf a = F(2 * r, i), b = F(2 * r + 1, i);
      const float ar = a.real() * a.real();
      const float ai = a.imag() * a.imag();
      const float br = b.real() * b.real();
      const float bi = b.imag() * b.imag();
      const float pa = ar + ai;
      const float pb = br + bi;
      P(r, i) = pa + pb;
    }
  return P;
}

Sweep measure(const cfvec& x, double fs, double fc) {
  Sweep none;
  if (x.size() < 4 * kNfft) return none;
  const fmat Q = occupancy(x);
  const double fr = fs / kNfft, T = static_cast<double>(x.size()) / fs, df = fs / kNfft;
  const size_t edge = std::max<size_t>(1, kNfft / 16), nb = kNfft - 2 * edge;
  fmat P(Q.rows, nb);   // fftshift, outer 1/16 each side dropped: bins low -> high frequency
  for (size_t r = 0; r < Q.rows; ++r)
    for (size_t j = 0; j < nb; ++j) P(r, j) = Q(r, (j + edge + kNfft / 2) % kNfft);
  const auto [rates, top_db] = candidates(P, fr, T);
  none.top_line_db = py_round(top_db, 1);
  for (double r : rates) {
    Json info = {{"line_hz", py_round(r, 1)}};
    Sweep sw = confirm(P, fr, T, r, df, info);
    none.tried.push_back(info);
    if (sw.found) {
      sw.freqs.resize(sw.nph);
      const double off = static_cast<double>(edge) - static_cast<double>(kNfft / 2);
      for (size_t i = 0; i < sw.nph; ++i) sw.freqs[i] = fc + (sw.fk[i] + edge - static_cast<double>(kNfft / 2)) * df;
      sw.occupied_hz = {fc + (sw.occupied_hz.first + off - 0.5) * df, fc + (sw.occupied_hz.second + off - 0.5) * df};
      sw.top_line_db = none.top_line_db;
      sw.tried = none.tried;
      return sw;
    }
  }
  return none;
}

std::vector<Emission> detect(const cfvec& x, double fs, double fc) {
  const Sweep s = measure(x, fs, fc);
  if (!s.found) return {};
  const double T = static_cast<double>(x.size()) / fs;
  std::vector<size_t> act;
  for (size_t i = 0; i < s.nph; ++i)
    if (s.active[i]) act.push_back(i);
  double f_min = s.freqs[act[0]], f_max = f_min;
  dvec pk;
  for (size_t i : act) {
    f_min = std::min(f_min, s.freqs[i]);
    f_max = std::max(f_max, s.freqs[i]);
    pk.push_back(s.pk_db[i]);
  }
  const bool below = std::find(s.beyond.begin(), s.beyond.end(), "below") != s.beyond.end();
  const bool above = std::find(s.beyond.begin(), s.beyond.end(), "above") != s.beyond.end();
  const double lo = below ? fc - fs / 2 : s.occupied_hz.first;
  const double hi = above ? fc + fs / 2 : s.occupied_hz.second;
  const size_t step = std::max<size_t>(1, act.size() / kTrajPoints);
  Json traj = Json::array();
  for (size_t j = 0; j < act.size(); j += step) {
    const size_t i = act[j];
    traj.push_back({py_round(static_cast<double>(i) / static_cast<double>(s.nph) * s.period * 1e6, 2),
                    py_round(s.freqs[i] / 1e6, 3)});
  }
  const double fold_peak = py_round(np::median(pk), 1);
  Json F = Json::object();
  F["sweep_rate_hz"] = py_round(1.0 / s.period, 2);
  F["sweep_period_us"] = py_round(s.period * 1e6, 3);
  F["sweep_low_mhz"] = py_round(f_min / 1e6, 3);
  F["sweep_high_mhz"] = py_round(f_max / 1e6, 3);
  F["sweep_span_mhz"] = py_round((f_max - f_min) / 1e6, 3);
  F["beyond_capture"] = s.beyond;
  F["shape"] = !s.flyback ? "up-down" : !s.exits ? "sawtooth" : "undetermined";
  F["in_band_fraction"] = py_round(static_cast<double>(act.size()) / static_cast<double>(s.nph), 3);
  F["flybacks"] = s.flyback;
  F["band_exits"] = s.exits;
  F["dominance_db"] = py_round(s.dom75, 1);
  F["fold_peak_db"] = fold_peak;
  F["trajectory_us_mhz"] = traj;
  F["modulation"] = "swept CW (FM by a periodic tuning voltage)";
  Emission e;
  e.family = "jammer";
  e.label = "swept_carrier";
  e.start_us = 0.0;
  e.end_us = T * 1e6;
  e.f_lo_mhz = lo / 1e6;
  e.f_hi_mhz = hi / 1e6;
  e.center_mhz = (lo + hi) / 2e6;
  e.snr_db = fold_peak;
  e.confidence = s.dom75 >= 6.0 && s.span >= 20e6 ? 0.9 : 0.75;
  e.features = F;
  return {e};
}

}  // namespace si::rx::jammer
