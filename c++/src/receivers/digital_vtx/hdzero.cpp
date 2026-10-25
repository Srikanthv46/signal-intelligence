// python receivers/digital_vtx/hdzero.py (constants, slice_plan, _fold_phase, _coherence, _preambles, hdzero_links).
#include "receivers/digital_vtx/hdzero.hpp"

#include <algorithm>
#include <cmath>

#include "core/fft.hpp"
#include "core/npcompat.hpp"
#include "receivers/digital_vtx/common.hpp"

namespace si::rx::digital_vtx::hdzero {

const int CONTINUAL[N_CONTINUAL] = {-444, -400, -360, -312, -280, -240, -224, -211, -200, -168, -136, -96, -56,
                                    0,    56,   104,  128,  160,  184,  216,  264,  280,  312,  364,  404,  444};

namespace {

using detail::kPi;
using detail::expj;

constexpr double MAX_PPM = 25.0;
constexpr double PILOT_MIN = 0.2, PILOT_CONTRAST = 4.0;
constexpr int64_t MIN_SYMBOLS = 24;
constexpr double PREAMBLE_RHO = 0.3;
constexpr int64_t FRAME_LO = 300, FRAME_HI = 800;
constexpr int64_t NW_LO = 1000, NW_HI = 1048;
constexpr double NW_TOL = 0.02;

struct Channel { const char* name; double mhz; };
const std::vector<Channel>& channels() {
  static const std::vector<Channel> c = [] {
    std::vector<Channel> v = {{"R1", 5658}, {"R2", 5695}, {"R3", 5732}, {"R4", 5769}, {"R5", 5806}, {"R6", 5843},
                              {"R7", 5880}, {"R8", 5917}, {"E1", 5705}, {"F1", 5740}, {"F2", 5760}, {"F4", 5800}};
    static const char* L[] = {"L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8"};
    for (int i = 0; i < 8; ++i) v.push_back({L[i], 5362.0 + 37 * i});
    return v;
  }();
  return c;
}

bool is_continual(int64_t k) {
  for (int c : CONTINUAL)
    if (c == k) return true;
  return false;
}

// python _REF: data carriers +-2 from the pilots, not pilots, |k| <= HALF, not on a scattered-pilot position
const std::vector<int64_t>& ref_carriers() {
  static const std::vector<int64_t> r = [] {
    std::vector<int64_t> v;
    for (int c : CONTINUAL)
      for (int d : {-2, 2}) v.push_back(c + d);
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    std::vector<int64_t> o;
    for (int64_t k : v)
      if (!is_continual(k) && std::llabs(k) <= HALF && detail::pymod(k - 4, 8) != 0) o.push_back(k);
    return o;
  }();
  return r;
}

double np_rint(double v) { return std::nearbyint(v); }

// mean over cols of |mean_n q| / (mean_n |q| + 1e-30), q = S[n+1, c] conj(S[n, c])
double coherence(const std::vector<cdvec>& S, const std::vector<int64_t>& cols) {
  const size_t M = S.size();
  dvec coh(cols.size());
  cdvec q(M - 1);
  dvec a(M - 1);
  for (size_t j = 0; j < cols.size(); ++j) {
    const size_t c = static_cast<size_t>(cols[j]);
    for (size_t n = 0; n + 1 < M; ++n) {
      q[n] = np::cmul(S[n + 1][c], std::conj(S[n][c]));
      a[n] = np::cabs(q[n]);
    }
    coh[j] = np::cabs(np::mean(q)) / (np::mean(a) + 1e-30);
  }
  return np::mean(coh);
}

}  // namespace

bool smooth(int64_t n) {
  for (int64_t p : {2, 3, 5, 7, 11, 13})
    while (n % p == 0) n /= p;
  return n == 1;
}

std::pair<int64_t, int64_t> slice_plan(int64_t n_samples, double fs) {
  const double m0 = static_cast<double>(n_samples) * FS / fs;
  std::vector<int64_t> nws;
  for (int64_t v = NW_LO; v <= NW_HI; ++v) nws.push_back(v);
  std::stable_sort(nws.begin(), nws.end(), [](int64_t a, int64_t b) { return std::llabs(a - N) < std::llabs(b - N); });
  std::optional<std::pair<int64_t, int64_t>> best;
  for (int64_t nw : nws) {
    if (best && std::llabs(nw - N) > std::llabs(best->second - N)) break;
    const double m = m0 * static_cast<double>(nw) / static_cast<double>(N);
    const int64_t lo = static_cast<int64_t>(std::ceil(m * (1 - NW_TOL / N)));
    const int64_t hi = static_cast<int64_t>(std::floor(m * (1 + NW_TOL / N)));
    for (int64_t c = lo; c <= hi; ++c)
      if (smooth(c)) {
        best = std::make_pair(c, nw);
        break;
      }
  }
  return best ? *best : std::make_pair(py_round_int(m0), N);
}

std::vector<Link> hdzero_links(const Slicer& sl, double fs, double fc, double usable, int64_t n_samples,
                               std::vector<LogEntry>* log) {
  std::vector<Link> out;
  const double span = (HALF + 0.5) * SCS;
  const auto [m1, nw] = slice_plan(n_samples, fs);
  const double out_fs = static_cast<double>(m1) * fs / static_cast<double>(n_samples);
  for (const auto& ch : channels()) {
    const double off = ch.mhz * 1e6 - fc;
    if (std::fabs(off) + span > usable) continue;
    double rate = 0;
    cdvec y = sl(off, span + 0.15e6, out_fs, &rate);
    const double tsf = static_cast<double>(TS) * rate / FS;
    const int64_t lag = nw;
    const double rho = static_cast<double>(nw) / (static_cast<double>(N) * rate / FS);
    auto col = [&](int64_t k, int64_t s) {
      return detail::pymod(static_cast<int64_t>(np_rint(static_cast<double>(k) * rho)) + s, nw);
    };
    const int64_t L = static_cast<int64_t>(y.size());
    const int64_t n_sym = static_cast<int64_t>(detail::py_floordiv(static_cast<double>(L - lag - 2 * TS), tsf));
    if (n_sym < MIN_SYMBOLS) continue;
    // ---- _fold_phase
    const int64_t nc = L - lag;
    cdvec cs(static_cast<size_t>(nc + 1), cd(0, 0));
    for (int64_t i = 0; i < nc; ++i)
      cs[static_cast<size_t>(i + 1)] = cs[static_cast<size_t>(i)] +
                                       np::cmul(y[static_cast<size_t>(i + lag)], std::conj(y[static_cast<size_t>(i)]));
    const int64_t W = static_cast<int64_t>(tsf);
    cdvec f(static_cast<size_t>(W), cd(0, 0));
    for (int64_t k = 0; k < n_sym; ++k) {
      const int64_t b = static_cast<int64_t>(np_rint(static_cast<double>(k) * tsf));
      for (int64_t j = 0; j < W; ++j) f[static_cast<size_t>(j)] += cs[static_cast<size_t>(b + j + CP)] - cs[static_cast<size_t>(b + j)];
    }
    dvec af(f.size());
    for (size_t j = 0; j < f.size(); ++j) af[j] = np::cabs(f[j]);
    const int64_t ph = static_cast<int64_t>(std::max_element(af.begin(), af.end()) - af.begin());
    dvec tmp = af;
    const double contrast = af[static_cast<size_t>(ph)] / (np::median_inplace(tmp.data(), tmp.size()) + 1e-30);
    const cd acc = f[static_cast<size_t>(ph)];
    const double cfo_frac = std::atan2(acc.imag(), acc.real()) / (2 * kPi * static_cast<double>(lag) / rate);
    for (int64_t i = 0; i < L; ++i)
      y[static_cast<size_t>(i)] = np::cmul(y[static_cast<size_t>(i)],
                                           expj(-2 * kPi * cfo_frac * static_cast<double>(i) / rate));
    // ---- per-symbol spectra on the fractional grid, window centred in the CP
    const int64_t nk = static_cast<int64_t>(detail::py_floordiv(static_cast<double>(L - ph - CP - nw), tsf));
    std::vector<cdvec> S;
    S.reserve(static_cast<size_t>(std::max<int64_t>(nk, 0)));
    for (int64_t k = 0; k < nk; ++k) {
      const double exact = static_cast<double>(ph + CP / 2) + static_cast<double>(k) * tsf;
      const int64_t st = static_cast<int64_t>(np_rint(exact));
      cdvec row = fft::fft(y.data() + st, static_cast<size_t>(nw));
      const double d = exact - static_cast<double>(st);
      for (int64_t b = 0; b < nw; ++b) {
        const int64_t fb = b < (nw + 1) / 2 ? b : b - nw;    // np.fft.fftfreq(nw) * nw
        row[static_cast<size_t>(b)] = np::cmul(row[static_cast<size_t>(b)],
                                               expj(2 * kPi * d * static_cast<double>(fb) / static_cast<double>(nw)));
      }
      S.push_back(std::move(row));
    }
    if (S.size() < 2) continue;
    const int64_t smax = static_cast<int64_t>(std::ceil(MAX_PPM * 1e-6 * ch.mhz * 1e6 / SCS)) + 1;
    int64_t shift = 0;
    double coh = -1, ref = 0;
    std::vector<int64_t> cols;
    for (int64_t s = -smax; s <= smax; ++s) {
      cols.clear();
      for (int c : CONTINUAL) cols.push_back(col(c, s));
      const double cp_ = coherence(S, cols);
      cols.clear();
      for (int64_t c : ref_carriers()) cols.push_back(col(c, s));
      const double rf = coherence(S, cols);
      if (coh < 0 || cp_ > coh) shift = s, coh = cp_, ref = rf;
    }
    const bool ok = coh >= PILOT_MIN && coh >= PILOT_CONTRAST * ref;
    if (log) log->push_back({ch.name, contrast, coh, ref, shift, ok});
    if (!ok) continue;
    const double cfo = static_cast<double>(shift) * rate / static_cast<double>(nw) + cfo_frac;
    dvec p(static_cast<size_t>(L));
    for (int64_t i = 0; i < L; ++i) p[static_cast<size_t>(i)] = np::cabs2(y[static_cast<size_t>(i)]);
    const double p_band = np::mean(p);
    // ---- _preambles
    std::vector<int64_t> keepc;
    for (int64_t k = -HALF; k <= HALF; ++k)
      if (!is_continual(k) && detail::pymod(k - 4, 8) != 0) keepc.push_back(col(k, shift));
    const size_t M = S.size();
    std::vector<cdvec> X(M, cdvec(keepc.size()));
    for (size_t n = 0; n < M; ++n) {
      double nn = 0;
      for (size_t j = 0; j < keepc.size(); ++j) {
        X[n][j] = S[n][static_cast<size_t>(keepc[j])];
        nn += np::cabs2(X[n][j]);
      }
      const double nr = std::sqrt(nn) + 1e-30;
      for (auto& v : X[n]) v = cd(v.real() / nr, v.imag() / nr);
    }
    std::vector<uint8_t> rep(M, 0);
    for (size_t d = 1; d <= 4; ++d)
      for (size_t n = 0; n + d < M; ++n) {
        cd s(0, 0);
        for (size_t j = 0; j < keepc.size(); ++j) s += np::cmul(X[n + d][j], std::conj(X[n][j]));
        if (np::cabs(s) > PREAMBLE_RHO) rep[n] = rep[n + d] = 1;
      }
    std::vector<std::pair<int64_t, int64_t>> blocks;
    int64_t bs = -1, bp = -1;
    for (size_t i = 0; i < M; ++i) {
      if (!rep[i]) continue;
      const int64_t ii = static_cast<int64_t>(i);
      if (bs < 0) {
        bs = bp = ii;
      } else {
        if (ii - bp > 6) {
          blocks.emplace_back(bs, bp);
          bs = ii;
        }
        bp = ii;
      }
    }
    if (bs >= 0) blocks.emplace_back(bs, bp);
    std::vector<std::pair<int64_t, int64_t>> bl;
    for (const auto& b : blocks)
      if (b.second - b.first >= 3) bl.push_back(b);
    Link lk;
    lk.channel = ch.name;
    lk.channel_mhz = ch.mhz;
    lk.off = off;
    lk.cfo_hz = cfo;
    lk.ppm = cfo / (ch.mhz * 1e6) * 1e6;
    lk.lo = off + cfo - span;
    lk.hi = off + cfo + span;
    lk.pilot_coherence = coh;
    lk.data_coherence = ref;
    lk.fold_contrast = contrast;
    lk.p_band = p_band;
    lk.half_bw = span;
    lk.n_symbols = static_cast<int64_t>(M);
    const double t_sym = tsf / rate, t0 = static_cast<double>(ph + CP) / rate;
    for (const auto& b : bl) {
      lk.preamble_ms.push_back(py_round((t0 + static_cast<double>(b.first) * t_sym) * 1e3, 3));
      lk.preamble_symbols.push_back(b.second - b.first + 1);
    }
    if (bl.size() >= 2) {
      dvec d;
      for (size_t i = 1; i < bl.size(); ++i) {
        const int64_t v = bl[i].first - bl[i - 1].first;
        if (v >= FRAME_LO && v <= FRAME_HI) d.push_back(static_cast<double>(v));
      }
      if (!d.empty()) lk.frame_symbols = static_cast<int64_t>(np::median_inplace(d.data(), d.size()));
    }
    if (lk.frame_symbols) lk.frame_ms = static_cast<double>(*lk.frame_symbols) * t_sym * 1e3;
    lk.t0 = 0.0;
    lk.t1 = static_cast<double>(n_samples) / fs;
    out.push_back(std::move(lk));
  }
  // one link per band (python: sort by -pilot_coherence, keep when > span from every kept link)
  std::stable_sort(out.begin(), out.end(), [](const Link& a, const Link& b) { return a.pilot_coherence > b.pilot_coherence; });
  std::vector<Link> keep;
  for (auto& d : out) {
    bool ok = true;
    for (const auto& k : keep) ok = ok && std::fabs(d.off + d.cfo_hz - k.off - k.cfo_hz) > 2 * span * 0.5;
    if (ok) keep.push_back(std::move(d));
  }
  return keep;
}

}  // namespace si::rx::digital_vtx::hdzero
