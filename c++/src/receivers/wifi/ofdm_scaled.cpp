// Port of python receivers/wifi/ofdm_scaled.py (see ofdm_scaled.hpp).
#include "receivers/wifi/ofdm_scaled.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>

#include "core/fft.hpp"
#include "core/ndimage.hpp"
#include "core/npcompat.hpp"
#include "core/parallel.hpp"
#include "core/pyfmt.hpp"
#include "core/signal.hpp"

namespace si::rx::wifi::scaled {

namespace {

constexpr double kPi = 3.141592653589793;
constexpr double kTuLo = 2.4e-6, kTuHi = 26e-6;
constexpr double kMinWidth = 1.2e6, kMaxWidth = 45e6;
using np::cabs;
using np::cabs2;
using np::cmul;

inline cd expj(double x) { return {std::cos(x), std::sin(x)}; }

// Python float floor division a // b (b > 0)
double py_floordiv(double a, double b) {
  const double mod = std::fmod(a, b);
  double div = (a - mod) / b;
  if (mod != 0 && ((b < 0) != (mod < 0))) div -= 1.0;
  if (div != 0) {
    double fl = std::floor(div);
    if (div - fl > 0.5) fl += 1.0;
    return fl;
  }
  return std::copysign(0.0, a / b);
}

// python _snap(clock) for a numpy float64 clock
double snap(double clock) {
  const double g = static_cast<double>(static_cast<int64_t>(std::nearbyint(clock / 0.5e6))) * 0.5e6;
  return (g != 0 && std::fabs(g - clock) / clock < 5e-3) ? g : np::round(clock, -3);
}

}  // namespace

Json Numerology::json() const {
  Json j = Json::object();
  j["center_hz"] = center_hz;
  j["occupied_hz"] = occupied_hz;
  j["tu_s"] = tu_s;
  j["ts_s"] = ts_s;
  j["guard"] = guard;
  j["cp_corr"] = cp_corr;
  j["cp_corr_over_median"] = cp_corr_over_median;
  j["symbol_line_over_median"] = symbol_line_over_median;
  if (clock_hz) j["clock_hz"] = *clock_hz;
  if (clock_measured_hz) j["clock_measured_hz"] = *clock_measured_hz;
  if (occupied_over_clock) j["occupied_over_clock"] = *occupied_over_clock;
  return j;
}

// ------------------------------------------------------------------ bands
std::pair<std::vector<Band>, std::vector<Band>> find_bands(const cfvec& x, double fs, double fc, double thresh_db) {
  size_t n = size_t(1) << static_cast<int>(std::nearbyint(std::log2(fs / 60e3)));
  n = std::max<size_t>(64, std::min<size_t>(n, 4096));
  const size_t m = x.size() / n;
  const dvec win = np::hanning(n);
  Mat<cd> X(m, n);
  for (size_t r = 0; r < m; ++r)
    for (size_t c = 0; c < n; ++c) {
      const cf v = x[r * n + c];
      X(r, c) = cmul(cd(v.real(), v.imag()), cd(win[c], 0));
    }
  fft::fft_rows(X, true);
  dmat P(m, n);
  const size_t h = n / 2;
  for (size_t r = 0; r < m; ++r)
    for (size_t c = 0; c < n; ++c) P(r, (c + h) % n) = cabs2(X(r, c));
  X = Mat<cd>();
  const double bin_hz = fs / static_cast<double>(n);
  dvec floor = np::percentile_axis0(P, 20.0);
  const size_t k = std::max<size_t>(3, static_cast<size_t>(3e6 / bin_hz));
  floor = ndi::maximum_filter1d(ndi::minimum_filter1d(floor, k, ndi::Mode::Nearest), k, ndi::Mode::Nearest);
  const double fac = std::pow(10.0, thresh_db / 10);
  std::vector<uint8_t> occ0(m * n), occ(m * n);
  const size_t edge = std::max<size_t>(1, n / 50);
  for (size_t r = 0; r < m; ++r)
    for (size_t c = 0; c < n; ++c) occ0[r * n + c] = (c >= edge && c < n - edge) && P(r, c) > floor[c] * fac;
  for (size_t r = 0; r < m; ++r) {
    const uint8_t* o = &occ0[r * n];
    for (size_t c = 0; c < n; ++c) {
      auto at = [&](int64_t d) { return o[((static_cast<int64_t>(c) + d) % static_cast<int64_t>(n) + n) % n]; };
      occ[r * n + c] = o[c] | (at(-1) & at(1)) | (at(-2) & at(2));
    }
  }
  const double lo_w = kMinWidth / bin_hz, hi_w = kMaxWidth / bin_hz;
  std::map<std::pair<int64_t, int64_t>, size_t> key_ix;
  std::vector<std::pair<std::pair<int64_t, int64_t>, std::vector<std::array<int64_t, 3>>>> edges;
  for (size_t t = 0; t < m; ++t) {
    const uint8_t* o = &occ[t * n];
    int64_t s = -1;
    for (size_t c = 0; c <= n; ++c) {
      const bool on = c < n && o[c];
      if (on && s < 0) s = static_cast<int64_t>(c);
      if (!on && s >= 0) {
        const int64_t e = static_cast<int64_t>(c);
        const double w = static_cast<double>(e - s);
        if (lo_w <= w && w <= hi_w) {
          const auto key = std::make_pair(s / 2, e / 2);
          auto it = key_ix.find(key);
          if (it == key_ix.end()) {
            key_ix[key] = edges.size();
            edges.push_back({key, {}});
            it = key_ix.find(key);
          }
          edges[it->second].second.push_back({static_cast<int64_t>(t), s, e});
        }
        s = -1;
      }
    }
  }
  const int64_t min_frames = std::max<int64_t>(4, static_cast<int64_t>(150e-6 * fs / static_cast<double>(n)));
  std::stable_sort(edges.begin(), edges.end(),
                   [](const auto& a, const auto& b) { return a.second.size() > b.second.size(); });
  std::vector<uint8_t> used(edges.size(), 0);
  std::vector<Band> bands;
  for (size_t ci = 0; ci < edges.size(); ++ci) {
    if (used[ci]) continue;
    const auto key = edges[ci].first;
    std::vector<std::array<int64_t, 3>> group;
    for (size_t cj = 0; cj < edges.size(); ++cj) {
      const auto k2 = edges[cj].first;
      if (!used[cj] && std::llabs(k2.first - key.first) <= 2 && std::llabs(k2.second - key.second) <= 2) {
        group.insert(group.end(), edges[cj].second.begin(), edges[cj].second.end());
        used[cj] = 1;
      }
    }
    if (static_cast<int64_t>(group.size()) < min_frames) continue;
    dvec ss, ee;
    for (const auto& g : group) {
      ss.push_back(static_cast<double>(g[1]));
      ee.push_back(static_cast<double>(g[2]));
    }
    const int64_t s = static_cast<int64_t>(np::median(ss)), e = static_cast<int64_t>(np::median(ee));
    Band b;
    b.frames.assign(m, 0);
    for (const auto& g : group) b.frames[g[0]] = 1;
    b.f_lo = fc - fs / 2 + static_cast<double>(s) * bin_hz;
    b.f_hi = fc - fs / 2 + static_cast<double>(e) * bin_hz;
    b.n_frames = static_cast<int64_t>(group.size());
    b.frame_len = n;
    bands.push_back(std::move(b));
  }
  std::stable_sort(bands.begin(), bands.end(), [](const Band& a, const Band& b) { return a.n_frames > b.n_frames; });
  std::vector<Band> unions;
  for (size_t i = 0; i < bands.size(); ++i)
    for (size_t j = i + 1; j < bands.size(); ++j) {
      const Band& a = bands[i];
      const Band& b = bands[j];
      const Band& lo = b.f_lo < a.f_lo ? b : a;
      const Band& hi = b.f_lo < a.f_lo ? a : b;
      int64_t both = 0;
      for (size_t t = 0; t < m; ++t) both += a.frames[t] & b.frames[t];
      if (-0.3e6 < hi.f_lo - lo.f_hi && hi.f_lo - lo.f_hi < 0.6e6 && hi.f_hi - lo.f_lo <= kMaxWidth &&
          static_cast<double>(both) >= 0.5 * static_cast<double>(std::min(a.n_frames, b.n_frames))) {
        Band u;
        u.frames.resize(m);
        int64_t cnt = 0;
        for (size_t t = 0; t < m; ++t) cnt += (u.frames[t] = a.frames[t] | b.frames[t]);
        u.f_lo = lo.f_lo;
        u.f_hi = hi.f_hi;
        u.n_frames = cnt;
        u.frame_len = n;
        u.merged = true;
        unions.push_back(std::move(u));
      }
    }
  std::stable_sort(unions.begin(), unions.end(), [](const Band& a, const Band& b) { return a.n_frames > b.n_frames; });
  return {bands, unions};
}

// ------------------------------------------------------------------ numerology
std::optional<Numerology> numerology(const cfvec& x, double fs, double fc, const Band& band) {
  double cen = (band.f_lo + band.f_hi) / 2;
  double bw = band.f_hi - band.f_lo;
  const int64_t dec = std::max<int64_t>(1, static_cast<int64_t>(py_floordiv(fs, 1.6 * bw)));
  const double fsd = fs / static_cast<double>(dec);
  const dvec h = sig::firwin(129, 0.55 * bw, fs);
  const double b = (-2.0 * kPi * (cen - fc)) / fs;
  cdvec mixed(x.size());
  for (size_t i = 0; i < x.size(); ++i)
    mixed[i] = cmul(cd(x[i].real(), x[i].imag()), expj(b * static_cast<double>(i)));
  cdvec yf = sig::oaconvolve(mixed, h, np::ConvMode::Same);
  mixed = cdvec();
  // gate to the active frames, closed over gaps of up to 3 frames
  const size_t nfr = band.frames.size();
  std::vector<uint8_t> fr_on(nfr, 0);
  for (size_t i = 0; i < nfr; ++i)
    for (int64_t d = -3; d <= 3; ++d) {
      const int64_t j = static_cast<int64_t>(i) + d;
      if (j >= 0 && j < static_cast<int64_t>(nfr) && band.frames[j]) fr_on[i] = 1;
    }
  const size_t ngate_full = std::min(nfr * band.frame_len, x.size());
  const size_t ngate = (ngate_full + dec - 1) / dec;
  const size_t ny_full = (yf.size() + dec - 1) / dec;
  const size_t ny = std::min(ngate, ny_full);
  cdvec y(ny);
  for (size_t i = 0; i < ny; ++i) {
    const size_t src = i * dec;
    y[i] = fr_on[src / band.frame_len] ? yf[src] : cd(0, 0);
  }
  yf = cdvec();
  dvec a2(ny);
  for (size_t i = 0; i < ny; ++i) a2[i] = cabs2(y[i]);
  const double e = np::sum(a2);
  if (!(e > 0)) return std::nullopt;
  // occupied band at -10 dB from the in-band level
  const size_t nseg = 256;
  cdvec on;
  for (size_t i = 0; i < ny; ++i)
    if (cabs(y[i]) > 0) on.push_back(y[i]);
  if (on.size() >= 4 * nseg) {
    const size_t rows = on.size() / nseg;
    const dvec hw = np::hanning(nseg);
    Mat<cd> S(rows, nseg);
    for (size_t r = 0; r < rows; ++r)
      for (size_t c = 0; c < nseg; ++c) S(r, c) = cmul(on[r * nseg + c], cd(hw[c], 0));
    fft::fft_rows(S, true);
    dmat A(rows, nseg);
    for (size_t i = 0; i < S.size(); ++i) A.data[i] = cabs2(S.data[i]);
    const dvec psd = fft::fftshift(np::mean_axis0(A));
    const dvec f = fft::fftshift(fft::fftfreq(nseg, 1 / fsd));
    dvec inner;
    for (size_t i = 0; i < nseg; ++i)
      if (std::fabs(f[i]) < 0.3 * bw) inner.push_back(psd[i]);
    const double lvl = np::median(inner);
    std::vector<uint8_t> hot(nseg);
    for (size_t i = 0; i < nseg; ++i) hot[i] = psd[i] > lvl / 10;
    dvec af(nseg);
    for (size_t i = 0; i < nseg; ++i) af[i] = std::fabs(f[i]);
    const size_t c0 = np::argmin(af);
    size_t aa = c0, bb = c0;
    while (aa > 0 && (hot[aa - 1] || (aa > 1 && hot[aa - 2]))) --aa;
    while (bb < nseg - 1 && (hot[bb + 1] || (bb < nseg - 2 && hot[bb + 2]))) ++bb;
    const double lo = f[aa] - 0.5 * fsd / nseg, hi = f[bb] + 0.5 * fsd / nseg;
    if (hi - lo > 0.3 * bw) {
      cen = cen + (lo + hi) / 2;
      bw = hi - lo;
    }
  }
  const int64_t lmin = static_cast<int64_t>(kTuLo * fsd);
  int64_t lmax = static_cast<int64_t>(std::ceil(kTuHi * fsd));
  lmax = std::min<int64_t>(lmax, static_cast<int64_t>(ny / 4));
  if (lmax <= lmin + 4) return std::nullopt;
  size_t nf = size_t(1) << static_cast<int>(std::ceil(std::log2(2.0 * static_cast<double>(ny))));
  dvec r(lmax + 2);
  {
    cdvec F = fft::fft(y, nf);
    for (auto& z : F) z = cd(cabs2(z), 0);
    const cdvec ac = fft::ifft(F);
    for (int64_t i = 0; i < lmax + 2; ++i) r[i] = cabs(ac[i]) / e;
  }
  const dvec seg(r.begin() + lmin, r.begin() + lmax + 1);
  const int64_t i = static_cast<int64_t>(np::argmax(seg)) + lmin;
  const double peak = r[i], base = np::median(seg);
  if (peak < 0.06 || peak < 5 * base) return std::nullopt;
  double lag;
  {
    const double a = r[i - 1], bq = r[i], c = r[i + 1];
    const double den = a - 2 * bq + c;
    lag = den != 0 ? static_cast<double>(i) + 0.5 * (a - c) / den : static_cast<double>(i);
  }
  const double tu = lag / fsd;
  const int64_t L = static_cast<int64_t>(std::nearbyint(lag));
  const int64_t w = std::max<int64_t>(2, L / 4);
  const int64_t ncc = static_cast<int64_t>(ny) - L;
  const int64_t nm = ncc - w + 1;
  if (nm < 2) return std::nullopt;
  dvec mm(nm);
  {
    cdvec cc(ncc);
    for (int64_t k = 0; k < ncc; ++k) cc[k] = cmul(y[k], std::conj(y[k + L]));
    for (int64_t k = 0; k < nm; ++k) {
      cd s = cc[k];
      for (int64_t j = 1; j < w; ++j) s += cc[k + j];
      mm[k] = cabs(s);
    }
    const double mean = np::mean(mm);
    for (auto& v : mm) v = v - mean;
  }
  nf = size_t(1) << static_cast<int>(std::ceil(std::log2(static_cast<double>(nm) * 4)));
  const cdvec Sc = fft::rfft(mm, nf);
  dvec S(Sc.size());
  for (size_t q = 0; q < Sc.size(); ++q) S[q] = cabs(Sc[q]);
  const dvec fr = fft::rfftfreq(nf, 1 / fsd);
  const double frame_rate = fs / static_cast<double>(band.frame_len);
  std::vector<size_t> selix;
  for (size_t q = 0; q < fr.size(); ++q) {
    bool sel = fr[q] > 0.6 / tu && fr[q] < 0.97 / tu;
    const double harm = std::nearbyint(fr[q] / frame_rate);
    if (harm >= 1 && std::fabs(fr[q] - harm * frame_rate) < 0.01 * frame_rate) sel = false;
    if (sel) selix.push_back(q);
  }
  if (selix.empty()) return std::nullopt;
  dvec ssel;
  for (size_t q : selix) ssel.push_back(S[q]);
  const size_t j = selix[np::argmax(ssel)];
  const double line = S[j] / (np::median(ssel) + 1e-20);
  double fj;
  {
    const double a = S[j - 1], bq = S[j], c = S[j + 1];
    const double den = a - 2 * bq + c;
    fj = fr[j] + (den != 0 ? 0.5 * (a - c) / den : 0.0) * (fr[1] - fr[0]);
  }
  const double ts = 1 / fj;
  const double gi = ts / tu - 1;
  Numerology out;
  out.center_hz = cen;
  out.occupied_hz = bw;
  out.tu_s = tu;
  out.ts_s = ts;
  out.guard = py_round(gi, 3);
  out.cp_corr = py_round(peak, 3);
  out.cp_corr_over_median = py_round(peak / std::max(base, 1e-9), 1);
  out.symbol_line_over_median = py_round(line, 1);
  if (!(0.2 <= gi && gi <= 0.3 && line > 3)) {
    out.reject = fmt("guard %.3f of the symbol (802.11 long guard is 0.25)", gi);
    return out;
  }
  const double clock = std::fabs(gi - 0.25) < 0.02 ? 80 / ts : 64 / tu;
  out.clock_hz = snap(clock);
  out.clock_measured_hz = clock;
  const double ratio = bw / *out.clock_hz;
  out.occupied_over_clock = py_round(ratio, 3);
  if (!(0.7 <= ratio && ratio <= 1.0))
    out.reject = fmt("occupied width %.2f MHz does not fit 52-56 of 64 subcarriers", bw / 1e6);
  return out;
}

// ------------------------------------------------------------------ decode
namespace {

struct Decoded {
  std::vector<ofdm::Packet> pk;
  size_t n_rejects = 0;
};

Decoded decode(const cfvec& x, const cfvec& xc, double fs, double fc, double center_hz, double clock_hz, bool inverted) {
  const double k = 20e6 / clock_hz;
  const double fsv = fs * k;
  const double off = (center_hz - fc) * k * (inverted ? -1 : 1);
  std::vector<ofdm::Channel> chans;
  for (double d : {-1.0e6, 0.0, 1.0e6}) chans.push_back({"scaled", std::nullopt, off + d, off + d, fsv / 2, false});
  ofdm::DetectorConfig cfg;
  cfg.usable_bw = fsv;
  cfg.stf_thresh = 0.6;
  cfg.power_margin_db = 2.5;
  cfg.min_plateau = 40;
  ofdm::WifiDetector det(fsv, 0.0, cfg, chans);
  Decoded out;
  out.pk = det.run(inverted ? xc : x, 1);
  out.n_rejects = det.n_rejects();
  return out;
}

int evidence(const std::vector<ofdm::Packet>& pk) {
  int s = 0;
  for (const auto& p : pk) s += (p.sig_crc_ok && *p.sig_crc_ok) ? 2 : 1;
  return s;
}

void to_real(ofdm::Packet& p, double fs, double fc, double k, bool inverted, double clock) {
  const double sgn = inverted ? -1 : 1;
  const double cen = fc + sgn * p.center_mhz * 1e6 / k;  // numpy float64
  const double half = (p.fmt == "HT" ? 8.75e6 : 8.3e6) / k;
  p.center_mhz = np::round(cen / 1e6, 4);
  p.cfo_khz = py_round(sgn * p.cfo_khz / k, 3);
  p.start_us = py_round(static_cast<double>(p.start_sample) / fs * 1e6, 3);
  p.duration_us = py_round(p.duration_us * k, 2);
  p.end_us = py_round(p.start_us + p.duration_us, 3);
  if (p.energy_end_us) p.energy_end_us = np::round(*p.energy_end_us * k, 2);
  const double lo_e = fc - fs / 2, hi_e = fc + fs / 2;
  p.f_lo_mhz = lo_e > cen - half ? py_round(lo_e / 1e6, 4) : np::round((cen - half) / 1e6, 4);
  p.f_hi_mhz = hi_e < cen + half ? py_round(hi_e / 1e6, 4) : np::round((cen + half) / 1e6, 4);
  p.bandwidth_mhz = clock / 1e6;
  p.channel = fmt("%.2fMHz@%sMHz", cen / 1e6, fmt("%g", clock / 1e6).c_str());
}

bool standard(double clock_hz, double center_hz) {
  if (std::fabs(clock_hz - 20e6) > 0.06e6) return false;
  for (const auto& c : ofdm::all_channels())
    if (std::fabs(c.freq_hz - center_hz) < 1e6) return true;
  return false;
}

}  // namespace

std::vector<Link> scan(const cfvec& x, double fs, double fc, size_t max_bands, size_t max_unions) {
  std::vector<Link> links;
  if (fs < 4e6) return links;
  auto [bands, unions] = find_bands(x, fs, fc);
  std::vector<const Band*> cand;
  for (size_t i = 0; i < std::min(max_bands, bands.size()); ++i) cand.push_back(&bands[i]);
  for (size_t i = 0; i < std::min(max_unions, unions.size()); ++i) cand.push_back(&unions[i]);
  // numerology of every candidate (independent; two at a time: each holds a few large FFTs)
  const auto nms = pmap(cand, [&](const Band* b) { return numerology(x, fs, fc, *b); }, std::min(2, threads()));
  cfvec xc;  // the mirrored capture, made once when needed
  std::vector<std::pair<double, double>> seen;
  for (size_t bi = 0; bi < cand.size(); ++bi) {
    if (!nms[bi] || nms[bi]->reject || !nms[bi]->clock_hz) continue;
    Numerology nm = *nms[bi];
    if (standard(*nm.clock_hz, nm.center_hz) ||
        (std::fabs(*nm.clock_measured_hz - 20e6) < 0.3e6 && standard(20e6, nm.center_hz)))
      continue;
    bool dup_seen = false;
    for (const auto& [c, k] : seen)
      if (std::fabs(nm.center_hz - c) < 0.25 * *nm.clock_hz && std::fabs(*nm.clock_hz - k) < 0.01 * k) dup_seen = true;
    if (dup_seen) continue;
    seen.push_back({nm.center_hz, *nm.clock_hz});
    const double meas = *nm.clock_measured_hz;
    std::vector<double> grid;
    for (int64_t g = static_cast<int64_t>(meas * 0.985 / 0.5e6); g < static_cast<int64_t>(meas * 1.015 / 0.5e6) + 2; ++g)
      if (std::fabs(static_cast<double>(g) * 0.5e6 - meas) <= 0.015 * meas) grid.push_back(static_cast<double>(g) * 0.5e6);
    std::stable_sort(grid.begin(), grid.end(),
                     [&](double a, double b) { return std::fabs(a - meas) < std::fabs(b - meas); });
    grid.push_back(py_round(meas, -3));
    std::vector<double> clocks;
    for (double c : grid)
      if (std::find(clocks.begin(), clocks.end(), c) == clocks.end()) clocks.push_back(c);
    if (xc.empty()) {
      xc.resize(x.size());
      for (size_t i = 0; i < x.size(); ++i) xc[i] = std::conj(x[i]);
    }
    std::vector<std::pair<double, bool>> trials;
    for (double c : clocks)
      for (bool inv : {false, true}) trials.push_back({c, inv});
    const auto res = pmap(trials, [&](const std::pair<double, bool>& t) {
      return decode(x, xc, fs, fc, nm.center_hz, t.first, t.second);
    });
    int best = -1;
    for (size_t t = 0; t < trials.size(); ++t)
      if (!res[t].pk.empty() && (best < 0 || evidence(res[t].pk) > evidence(res[best].pk))) best = static_cast<int>(t);
    if (best < 0) continue;
    Link ln;
    ln.clock_hz = trials[best].first;
    ln.inverted = trials[best].second;
    nm.clock_hz = ln.clock_hz;
    ln.numerology = nm;
    ln.packets = res[best].pk;
    ln.n_rejects = res[best].n_rejects;
    const double k = 20e6 / ln.clock_hz;
    dvec cens;
    for (auto& p : ln.packets) {
      to_real(p, fs, fc, k, ln.inverted, ln.clock_hz);
      cens.push_back(p.center_mhz);
    }
    ln.center_hz = np::median(cens) * 1e6;
    Link* dup = nullptr;
    for (auto& o : links)
      if (std::fabs(o.center_hz - ln.center_hz) < 0.25 * ln.clock_hz && o.clock_hz == ln.clock_hz) {
        dup = &o;
        break;
      }
    if (!dup) links.push_back(std::move(ln));
    else if (evidence(ln.packets) > evidence(dup->packets)) *dup = std::move(ln);
  }
  return links;
}

}  // namespace si::rx::wifi::scaled
