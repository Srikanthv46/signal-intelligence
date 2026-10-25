// python receivers/navaids/detector.py (see there for the method and the sources).
#include "receivers/navaids/detector.hpp"

#include <algorithm>
#include <cmath>

#include "core/fft.hpp"
#include "core/npcompat.hpp"

namespace si::rx::navaids {

namespace {

const std::pair<double, double> kIlsPairs[40] = {
    {108.10, 334.70}, {108.15, 334.55}, {108.30, 334.10}, {108.35, 333.95}, {108.50, 329.90}, {108.55, 329.75},
    {108.70, 330.50}, {108.75, 330.35}, {108.90, 329.30}, {108.95, 329.15}, {109.10, 331.40}, {109.15, 331.25},
    {109.30, 332.00}, {109.35, 331.85}, {109.50, 332.60}, {109.55, 332.45}, {109.70, 333.20}, {109.75, 333.05},
    {109.90, 333.80}, {109.95, 333.65}, {110.10, 334.40}, {110.15, 334.25}, {110.30, 335.00}, {110.35, 334.85},
    {110.50, 329.60}, {110.55, 329.45}, {110.70, 330.20}, {110.75, 330.05}, {110.90, 330.80}, {110.95, 330.65},
    {111.10, 331.70}, {111.15, 331.55}, {111.30, 332.30}, {111.35, 332.15}, {111.50, 332.90}, {111.55, 332.75},
    {111.70, 333.50}, {111.75, 333.35}, {111.90, 331.10}, {111.95, 330.95}};
constexpr double kLineDb = 15.0;
constexpr double kSearchHz = 35e3;
constexpr double kIlsFs = 6e3, kVorFs = 32e3;
constexpr double kSubHz = 9960.0, kSubBw = 700.0;

dvec scales() {   // np.linspace(0.97, 1.03, 61)
  dvec s(61);
  const double step = (1.03 - 0.97) / 60;
  for (int i = 0; i < 61; ++i) s[i] = 0.97 + i * step;
  s[60] = 1.03;
  return s;
}

// python _lines(ch, f_meas): carrier lines (Hz, dB over the window median), strongest first
std::optional<std::vector<std::pair<double, double>>> lines(const Channelizer& ch, double f_meas) {
  const double f_off = f_meas - ch.fc();
  const size_t L = ch.L();
  const double df = ch.fs() / static_cast<double>(L);
  int64_t lo = static_cast<int64_t>(std::floor((f_off - kSearchHz) / df)) - 2;
  int64_t hi = static_cast<int64_t>(std::ceil((f_off + kSearchHz) / df)) + 2;
  lo = std::max<int64_t>(lo, -static_cast<int64_t>(L / 2) - 1);
  hi = std::min<int64_t>(hi, static_cast<int64_t>((L - 1) / 2) + 1);
  std::vector<size_t> idx;
  for (int64_t k = lo; k <= hi; ++k) idx.push_back(static_cast<size_t>(((k % static_cast<int64_t>(L)) + L) % L));
  std::sort(idx.begin(), idx.end());
  idx.erase(std::unique(idx.begin(), idx.end()), idx.end());
  std::vector<size_t> k;
  for (size_t i : idx)
    if (std::fabs(ch.f()[i] - f_off) <= kSearchHz) k.push_back(i);
  if (k.size() < 64) return std::nullopt;
  std::vector<std::pair<double, double>> fp;   // (f, p) sorted by f
  for (size_t i : k) {
    const double a = static_cast<double>(np::cabs(ch.X()[i]));
    fp.emplace_back(ch.f()[i] + ch.fc(), a * a);
  }
  std::sort(fp.begin(), fp.end());
  dvec p(fp.size()), f(fp.size());
  for (size_t i = 0; i < fp.size(); ++i) f[i] = fp[i].first, p[i] = fp[i].second;
  const double floor = np::median(p) + 1e-30;
  dvec db(p.size()), neg(p.size());
  for (size_t i = 0; i < p.size(); ++i) {
    db[i] = 10 * std::log10(p[i] / floor + 1e-30);
    neg[i] = -db[i];
  }
  std::vector<std::pair<double, double>> out;
  for (int64_t i : np::argsort_stable(neg)) {
    if (db[i] < kLineDb || out.size() >= 2) break;
    bool far = true;
    for (const auto& [g, _] : out) far = far && std::fabs(f[i] - g) >= 2e3;
    if (far) out.emplace_back(f[i], db[i]);
  }
  return out;
}

// python _fit(a, fs, tones): least squares on [1, cos, sin per tone] -> coefficients, residual power
std::pair<dvec, double> fit(const dvec& a, double fs, const dvec& tones) {
  const size_t n = a.size(), m = 1 + 2 * tones.size();
  dmat H(n, m);
  for (size_t i = 0; i < n; ++i) {
    const double t = static_cast<double>(i) / fs;
    H(i, 0) = 1.0;
    for (size_t j = 0; j < tones.size(); ++j) {
      H(i, 1 + 2 * j) = std::cos(2 * kPi * tones[j] * t);
      H(i, 2 + 2 * j) = std::sin(2 * kPi * tones[j] * t);
    }
  }
  dmat A(m, m);
  dvec b(m, 0.0);
  for (size_t r = 0; r < m; ++r) {
    for (size_t c = 0; c < m; ++c) {
      double s = 0;
      for (size_t i = 0; i < n; ++i) s += H(i, r) * H(i, c);
      A(r, c) = s;
    }
    double s = 0;
    for (size_t i = 0; i < n; ++i) s += H(i, r) * a[i];
    b[r] = s;
  }
  const dvec c = np::solve(A, b);
  dvec r2(n);
  for (size_t i = 0; i < n; ++i) {
    double y = 0;
    for (size_t j = 0; j < m; ++j) y += H(i, j) * c[j];
    r2[i] = (a[i] - y) * (a[i] - y);
  }
  return {c, np::mean(r2)};
}

// python _whole(n, fs, rate)
size_t whole(size_t n, double fs, double rate = 30.0) {
  const double per = fs / rate;
  const int64_t k = static_cast<int64_t>(static_cast<double>(n) / per + 1e-6);
  return k >= 1 ? static_cast<size_t>(py_round_int(static_cast<double>(k) * per)) : 0;
}

struct Ils {
  double m90, m150, ddm, sdm, scale, sigma, window_s;
  bool edge;
};

std::optional<Ils> measure_ils(const Channelizer& ch, double f_line) {
  const auto sl = ch.slice(f_line, 1.5e3, kIlsFs);
  const size_t n = whole(sl.y.size(), sl.fs);
  if (n < 64) return std::nullopt;
  dvec a(n);
  for (size_t i = 0; i < n; ++i) a[i] = static_cast<double>(np::cabs(sl.y[i]));
  const dvec S = scales();
  dvec res;
  for (double s : S) res.push_back(fit(a, sl.fs, {90.0 * s, 150.0 * s}).second);
  const size_t i = static_cast<size_t>(std::min_element(res.begin(), res.end()) - res.begin());
  double s = S[i];
  const double step = S[1] - S[0];
  if (0 < i && i < S.size() - 1) {
    const double r0 = res[i - 1], r1 = res[i], r2 = res[i + 1];
    const double den = r0 - 2 * r1 + r2;
    if (den > 0) s += 0.5 * (r0 - r2) / den * step;
  }
  const auto [c, rp] = fit(a, sl.fs, {90.0 * s, 150.0 * s});
  const double a0 = c[0];
  Ils o;
  o.m90 = std::hypot(c[1], c[2]) / a0;
  o.m150 = std::hypot(c[3], c[4]) / a0;
  o.ddm = o.m90 - o.m150;
  o.sdm = o.m90 + o.m150;
  o.scale = s;
  o.sigma = std::sqrt(2.0 / (static_cast<double>(n) * (a0 * a0 / (rp + 1e-30))));
  o.edge = i == 0 || i == S.size() - 1;
  o.window_s = static_cast<double>(n) / sl.fs;
  return o;
}

struct Vor {
  double am_depth, deviation_hz, subcarrier_db, bearing_deg, window_s;
};

std::optional<Vor> measure_vor(const Channelizer& ch, double f_line) {
  const auto sl = ch.slice(f_line, 12.5e3, kVorFs);
  const double fs = sl.fs;
  const size_t n = whole(sl.y.size(), fs);
  if (n < 256) return std::nullopt;
  dvec a(n);
  for (size_t i = 0; i < n; ++i) a[i] = static_cast<double>(np::cabs(sl.y[i]));
  const double mu = np::mean(a);
  for (double& v : a) v = v / mu;
  const auto [c, _r] = fit(a, fs, {30.0});
  const double am_depth = std::hypot(c[1], c[2]), am_phase = std::atan2(c[2], c[1]);
  cdvec ac(n);
  for (size_t i = 0; i < n; ++i) ac[i] = cd(a[i], 0.0);
  const cdvec A = fft::fft(ac);
  const double val = 1.0 / (static_cast<double>(n) * (1 / fs));
  std::vector<size_t> sub;
  double ps = 0, pn = 0;
  size_t nside = 0;
  for (size_t k = 0; k < n; ++k) {
    const int64_t kk = k < (n + 1) / 2 ? static_cast<int64_t>(k) : static_cast<int64_t>(k) - static_cast<int64_t>(n);
    const double fr = static_cast<double>(kk) * val;
    const double d = std::fabs(fr - kSubHz);
    if (d <= kSubBw) sub.push_back(k);
    if (d > 1.5 * kSubBw && d <= 3 * kSubBw) ++nside;
  }
  if (sub.size() < 8 || nside < 8) return std::nullopt;
  {
    dvec vs, vn;
    for (size_t k = 0; k < n; ++k) {
      const int64_t kk = k < (n + 1) / 2 ? static_cast<int64_t>(k) : static_cast<int64_t>(k) - static_cast<int64_t>(n);
      const double d = std::fabs(static_cast<double>(kk) * val - kSubHz);
      const double m = np::cabs(A[k]);
      if (d <= kSubBw) vs.push_back(m * m);
      if (d > 1.5 * kSubBw && d <= 3 * kSubBw) vn.push_back(m * m);
    }
    ps = np::mean(vs);
    pn = np::mean(vn);
  }
  const double sub_db = 10 * std::log10(ps / (pn + 1e-30) + 1e-30);
  cdvec Z(n, cd(0, 0));
  const int64_t k0 = py_round_int(kSubHz * static_cast<double>(n) / fs);
  for (size_t k : sub) {
    const int64_t j = ((static_cast<int64_t>(k) - k0) % static_cast<int64_t>(n) + static_cast<int64_t>(n)) % static_cast<int64_t>(n);
    Z[j] = A[k];
  }
  const cdvec z = fft::ifft(Z);
  dvec fi(n - 1);
  // numpy temporary elision (see voice_fm): a conj() temporary >= 256 KiB receives the product in place
  const bool elide = (n - 1) * sizeof(cd) >= 256 * 1024;
  for (size_t i = 1; i < n; ++i) {
    const cd q = elide ? np::cmul(std::conj(z[i - 1]), z[i]) : np::cmul(z[i], std::conj(z[i - 1]));
    fi[i - 1] = std::atan2(q.imag(), q.real()) * fs / (2 * kPi);
  }
  // fit at the midpoints t = (i + 0.5) / fs
  const size_t m = fi.size();
  dmat A3(3, 3);
  dvec b3(3, 0.0);
  dmat H(m, 3);
  for (size_t i = 0; i < m; ++i) {
    const double t = (static_cast<double>(i) + 0.5) / fs;
    H(i, 0) = 1.0;
    H(i, 1) = std::cos(2 * kPi * 30 * t);
    H(i, 2) = std::sin(2 * kPi * 30 * t);
  }
  for (size_t r = 0; r < 3; ++r) {
    for (size_t cc = 0; cc < 3; ++cc) {
      double s = 0;
      for (size_t i = 0; i < m; ++i) s += H(i, r) * H(i, cc);
      A3(r, cc) = s;
    }
    double s = 0;
    for (size_t i = 0; i < m; ++i) s += H(i, r) * fi[i];
    b3[r] = s;
  }
  const dvec cf3 = np::solve(A3, b3);
  Vor o;
  o.am_depth = am_depth;
  o.deviation_hz = std::hypot(cf3[1], cf3[2]);
  const double fm_phase = std::atan2(cf3[2], cf3[1]);
  double deg = (am_phase - fm_phase) * (180.0 / kPi);   // np.degrees(x) % 360.0
  deg = std::fmod(deg, 360.0);
  if (deg < 0) deg += 360.0;
  o.bearing_deg = deg;
  o.subcarrier_db = sub_db;
  o.window_s = static_cast<double>(n) / fs;
  return o;
}

}  // namespace

std::vector<Channel> channels() {
  std::vector<Channel> out;
  for (const auto& [l, g] : kIlsPairs) out.push_back({"ils_gs", g, l});
  for (const auto& [l, g] : kIlsPairs) out.push_back({"ils_loc", l, g});
  for (int k = 0; k < 200; ++k) {
    const double f = py_round(108.0 + 0.05 * k, 2);
    if (f < 112.0 && (py_round_int(f * 100) / 10) % 2) continue;
    out.push_back({"vor", f, std::nullopt});
  }
  return out;
}

std::vector<Emission> detect(const Channelizer& ch, double ppm) {
  const double lo = ch.fc() - 0.45 * ch.fs(), hi = ch.fc() + 0.45 * ch.fs();
  const double dur_us = static_cast<double>(ch.L()) / ch.fs() * 1e6;
  std::vector<Emission> out;
  for (const Channel& c : channels()) {
    const double f_nom = c.mhz * 1e6;
    const double f_meas = f_nom / (1 + ppm * 1e-6);
    if (!(lo + kSearchHz < f_meas && f_meas < hi - kSearchHz)) continue;
    const auto ln = lines(ch, f_meas);
    if (!ln || ln->empty()) continue;
    const auto [f0, db0] = (*ln)[0];
    std::optional<int64_t> two;
    if (ln->size() > 1 && 4e3 <= std::fabs((*ln)[1].first - f0) && std::fabs((*ln)[1].first - f0) <= 32e3 &&
        (*ln)[1].second >= db0 - 15)
      two = py_round_int(std::fabs((*ln)[1].first - f0));
    if (two && *two == 0) two.reset();   // python: `if two:`
    Json F = Json::object();
    F["kind"] = c.kind;
    F["channel_mhz"] = c.mhz;
    F["carrier_mhz"] = py_round(f0 / 1e6, 6);
    F["offset_hz"] = py_round_int(f0 * (1 + ppm * 1e-6) - f_nom);
    F["line_db"] = py_round(db0, 1);
    F["receiver_ppm_applied"] = ppm;
    if (two) F["two_frequency_hz"] = *two;
    bool ok;
    double half;
    if (c.kind == "vor") {
      const auto m = measure_vor(ch, f0);
      if (!m) continue;
      ok = m->subcarrier_db >= 10 && 0.15 <= m->am_depth && m->am_depth <= 0.45 && 300 <= m->deviation_hz &&
           m->deviation_hz <= 700;
      if (m->subcarrier_db < 6) continue;
      F["radial_deg"] = py_round(m->bearing_deg, 1);
      F["am_depth"] = py_round(m->am_depth, 3);
      F["deviation_hz"] = py_round_int(m->deviation_hz);
      F["subcarrier_db"] = py_round(m->subcarrier_db, 1);
      F["window_s"] = py_round(m->window_s, 4);
      half = 12.5e3;
    } else {
      const auto m = measure_ils(ch, f0);
      if (!m) continue;
      const bool tones = std::min(m->m90, m->m150) >= 5 * m->sigma;
      if (!tones && std::max(m->m90, m->m150) < 5 * m->sigma) continue;
      const bool sdm_ok = c.kind == "ils_gs" ? (0.6 <= m->sdm && m->sdm <= 1.0) : (0.25 <= m->sdm && m->sdm <= 0.65);
      ok = tones && sdm_ok && !m->edge;
      const std::string side = c.kind == "ils_gs" ? (m->ddm > 0 ? "above path" : "below path")
                                                  : (m->ddm > 0 ? "left of course" : "right of course");
      F["paired_mhz"] = *c.paired;
      F["m90"] = py_round(m->m90, 4);
      F["m150"] = py_round(m->m150, 4);
      F["ddm"] = py_round(m->ddm, 4);
      F["sdm"] = py_round(m->sdm, 4);
      F["side"] = side;
      F["tone_scale"] = py_round(m->scale, 4);
      F["depth_sigma"] = py_round(m->sigma, 5);
      F["window_s"] = py_round(m->window_s, 4);
      half = static_cast<double>(two ? *two : 0) / 2 + 2e3;
    }
    F["verified_structure"] = ok;
    Emission e;
    e.family = "navaids";
    e.label = c.kind;
    e.start_us = 0.0;
    e.end_us = dur_us;
    if (two) {
      const double f1 = (*ln)[1].first;
      e.f_lo_mhz = (std::min(f0, f1) - 2e3) / 1e6;
      e.f_hi_mhz = (std::max(f0, f1) + 2e3) / 1e6;
    } else {
      e.f_lo_mhz = (f0 - half) / 1e6;
      e.f_hi_mhz = (f0 + half) / 1e6;
    }
    e.center_mhz = f0 / 1e6;
    e.snr_db = py_round(db0, 1);
    e.confidence = ok ? 0.8 : 0.55;
    e.features = std::move(F);
    out.push_back(std::move(e));
  }
  return out;
}

}  // namespace si::rx::navaids
