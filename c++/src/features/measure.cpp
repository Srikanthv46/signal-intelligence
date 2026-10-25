#include "features/measure.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

#include "core/deinterleave.hpp"
#include "core/bands_india.hpp"
#include "core/fft.hpp"
#include "core/npcompat.hpp"
#include "core/pyfmt.hpp"

namespace si::features {

const size_t kMaxSamples = size_t(1) << 19;

Json timing(const Signal& sig, double capture_s) {
  std::vector<Box> b = sig.bursts;
  std::stable_sort(b.begin(), b.end(), [](const Box& x, const Box& y) { return x.t0 < y.t0; });
  double on = 0.0;
  dvec durs;
  bool cont = false;
  for (const auto& x : b) {
    on += x.t1 - x.t0;  // python sum(): sequential
    durs.push_back(x.t1 - x.t0);
    cont = cont || x.kind == "continuous";
  }
  Json out = {{"n_bursts", b.size()},
              {"duty_cycle", py_round(std::min(1.0, on / capture_s), 3)},
              {"burst_dur_s_median", np::median(durs.data(), durs.size())},
              {"burst_period_s", nullptr},
              {"continuous", cont}};
  if (b.size() >= 3) {
    dvec d(b.size() - 1);
    for (size_t i = 0; i + 1 < b.size(); ++i) d[i] = b[i + 1].t0 - b[i].t0;
    const double med = np::median(d.data(), d.size());
    if (med > 0 && np::stddev(d.data(), d.size()) / med < 0.15) out["burst_period_s"] = py_round(med, 6);
  }
  return out;
}

constexpr double kPulseThrDb = 10.0;   // 1 us envelope above the band's median level
constexpr double kPulseGapS = 20e-6;   // gaps shorter than this join one pulse group (DME pair, radar short + long)
constexpr double kPriMinS = 100e-6;    // group-to-group intervals below this are fragments of one group
constexpr double kPriTol = 0.15;       // intervals within +-15 % of the shortest cluster are one PRI (stagger)
constexpr double kPriMinFrac = 0.5;    // ... and that cluster holds at least half of the intervals
constexpr double kLfmMinS = 10e-6;     // chirp test on pulses at least this long
constexpr double kOfdmMinLags = 4;     // a CP-OFDM autocorrelation lag fits >= 4 times into its bursts
constexpr double kLfmR2 = 0.8;         // linear fit of the instantaneous frequency (1 us blocks)
constexpr size_t kDeintMin = 6;      // pulses before the deinterleaver runs (core/deinterleave)
constexpr size_t kDeintMax = 2000;   // ... and the most it takes (the first ones)
constexpr double kPulseMinBw = 100e3;  // a pulse of width tau occupies >~ 1/tau

// python pulse_train(): pulse groups over the whole capture on the signal's band, 1 us envelope; PRI; LFM fit
Json pulse_train(const Signal& sig, const Channelizer& ch) {
  Json out = Json::object();
  if (sig.bw() < kPulseMinBw) return out;
  const Channelizer::Slice sl = ch.slice(sig.fc(), 0.5 * sig.bw(), std::min(ch.fs(), std::max(1.25 * sig.bw(), 1e6)));
  const cfvec& y = sl.y;
  const double ofs = sl.fs;
  const int64_t m = std::max<int64_t>(1, py_round_int(ofs * 1e-6));
  const int64_t k = static_cast<int64_t>(y.size()) / m;
  if (k < 64) return out;
  fvec e(static_cast<size_t>(k)), row(static_cast<size_t>(m));
  for (int64_t r = 0; r < k; ++r) {
    for (int64_t c = 0; c < m; ++c) {
      const float a = np::cabs(y[static_cast<size_t>(r * m + c)]);
      row[static_cast<size_t>(c)] = a * a;
    }
    e[static_cast<size_t>(r)] = np::mean(row.data(), row.size());
  }
  const float thr = np::median(e) * static_cast<float>(std::pow(10.0, kPulseThrDb / 10));
  std::vector<int64_t> idx;
  for (int64_t i = 0; i < k; ++i)
    if (e[static_cast<size_t>(i)] > thr) idx.push_back(i);
  if (idx.empty()) return out;
  const double dt = static_cast<double>(m) / ofs;
  std::vector<int64_t> st{idx[0]}, en;
  for (size_t i = 0; i + 1 < idx.size(); ++i)
    if (static_cast<double>(idx[i + 1] - idx[i]) * dt > kPulseGapS) {
      en.push_back(idx[i]);
      st.push_back(idx[i + 1]);
    }
  en.push_back(idx.back());
  dvec w(st.size());
  for (size_t i = 0; i < st.size(); ++i) w[i] = static_cast<double>(en[i] - st[i] + 1) * dt;
  out["pulse_n"] = st.size();
  out["pulse_width_s_median"] = py_round(np::median(w), 7);
  out["pulse_duty"] = py_round(np::sum(w) / (static_cast<double>(k) * dt), 4);
  Json pt = Json::array();  // pulse starts (detect::merge_pulse_trains)
  for (int64_t v : st) pt.push_back(static_cast<double>(v) * dt);
  out["_pulse_t_s"] = pt;
  if (st.size() >= kDeintMin) {   // emitters by PRI (research, signal gaps S10)
    const size_t nd = std::min(st.size(), kDeintMax);
    std::vector<double> toa(nd), pw(w.begin(), w.begin() + static_cast<ptrdiff_t>(nd));
    for (size_t i = 0; i < nd; ++i) toa[i] = static_cast<double>(st[i]) * dt;
    const Json de = deint::deinterleave(toa, pw);
    out["pulse_emitters"] = de["emitters"];
    out["pulse_unassigned"] = de["unassigned"];
  }
  dvec d;
  for (size_t i = 0; i + 1 < st.size(); ++i) {
    const double v = static_cast<double>(st[i + 1] - st[i]) * dt;
    if (v >= kPriMinS) d.push_back(v);
  }
  if (d.size() >= 2) {
    // weak pulses below the threshold leave 2 x PRI gaps: the PRI is the shortest interval cluster
    const double p0 = np::percentile(d, 25.0);
    dvec c;
    for (double v : d)
      if (std::fabs(v - p0) <= kPriTol * p0) c.push_back(v);
    // random intervals also leave a few in the window: the cluster holds at least half of them
    if (static_cast<double>(c.size()) >= std::max(2.0, kPriMinFrac * static_cast<double>(d.size()))) {
      const double pri = np::median(c);
      out["pri_s"] = py_round(pri, 7);
      out["pri_spread"] = py_round((*std::max_element(c.begin(), c.end()) - *std::min_element(c.begin(), c.end())) / pri, 3);
    }
  }
  const size_t j = static_cast<size_t>(std::max_element(w.begin(), w.end()) - w.begin());
  if (w[j] >= kLfmMinS) {
    int64_t a = st[j] * m, b = (en[j] + 1) * m;
    const int64_t a2 = a + (b - a) / 10, b2 = b - (b - a) / 10;  // central 80 %: edges ride on the rise/fall
    a = a2;
    b = b2;
    // instantaneous frequency per 1 us block (phase step summed over the block)
    const int64_t nq = b - a - 1, nb = nq / m;
    if (nb >= 8) {
      cfvec q(static_cast<size_t>(nb * m));
      for (int64_t i = 0; i < nb * m; ++i)
        q[static_cast<size_t>(i)] = np::cmul(y[static_cast<size_t>(a + 1 + i)], std::conj(y[static_cast<size_t>(a + i)]));
      const float kofs = static_cast<float>(ofs), two_pi = static_cast<float>(2 * kPi);
      fvec fi(static_cast<size_t>(nb));
      dvec t(static_cast<size_t>(nb));
      for (int64_t r = 0; r < nb; ++r) {
        const cf s = np::sum(q.data() + r * m, static_cast<size_t>(m));
        fi[static_cast<size_t>(r)] = std::atan2(s.imag(), s.real()) * kofs / two_pi;
        t[static_cast<size_t>(r)] = (static_cast<double>(r) + 0.5) * static_cast<double>(m) / ofs;
      }
      const double tm = np::mean(t);
      const float fm = np::mean(fi);
      dvec tt(t.size()), tf(t.size());
      fvec ff(fi.size());
      for (size_t i = 0; i < t.size(); ++i) {
        const double tc = t[i] - tm;
        const float fc = fi[i] - fm;
        tt[i] = tc * tc;
        ff[i] = fc * fc;
        tf[i] = tc * static_cast<double>(fc);
      }
      const double stt = np::sum(tt), sff = static_cast<double>(np::sum(ff)) + 1e-30, stf = np::sum(tf);
      const double r2 = stf * stf / (stt * sff);
      const double sweep = stf / stt * static_cast<double>(b - a) / ofs;
      if (r2 >= kLfmR2 && std::fabs(sweep) > 2 * ofs / static_cast<double>(b - a)) {
        out["lfm_sweep_hz"] = py_round_int(sweep);
        out["lfm_pulse_s"] = py_round(w[j], 7);
        out["lfm_r2"] = py_round(r2, 3);
      }
    }
  }
  return out;
}

namespace {

// top-k (frequency, peak/median) of spec inside [fmin, fmax] (python _peaks)
std::vector<std::pair<double, double>> peaks(const fvec& spec, const dvec& freqs, double fmin, double fmax, size_t k = 3) {
  fvec v;
  dvec f;
  for (size_t i = 0; i < spec.size(); ++i)
    if (freqs[i] >= fmin && freqs[i] <= fmax) {
      v.push_back(spec[i]);
      f.push_back(freqs[i]);
    }
  if (v.size() < 8) return {};
  fvec tmp = v;
  const float med = np::median(tmp.data(), tmp.size()) + 1e-30f;
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  // np.argsort(v)[::-1][:k]: the k largest, ties: higher index first
  std::stable_sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return v[a] < v[b]; });
  std::vector<std::pair<double, double>> out;
  for (size_t j = 0; j < k && j < idx.size(); ++j) {
    const size_t i = idx[idx.size() - 1 - j];
    out.push_back({f[i], static_cast<double>(v[i] / med)});
  }
  return out;
}

size_t argmax_f(const fvec& v) { return static_cast<size_t>(std::max_element(v.begin(), v.end()) - v.begin()); }

}  // namespace

std::pair<Json, std::vector<std::string>> iq_features(const Signal& sig, const Channelizer& ch, size_t n) {
  std::vector<std::string> notes;
  const double bw = std::max(sig.bw(), 5e3);
  const double ofs_req = std::min(ch.fs(), std::max(1.25 * bw, 20e3));
  Channelizer::Slice sl = ch.slice(sig.fc(), 0.55 * bw, ofs_req);
  const double ofs = sl.fs;
  const std::vector<uint8_t> m = burst_mask(n, ch.fs(), sig.bursts);
  cfvec y;
  y.reserve(sl.y.size());
  for (size_t i = 0; i < sl.y.size(); ++i) {
    int64_t j = static_cast<int64_t>(static_cast<double>(i) * ch.fs() / ofs);
    j = std::clamp<int64_t>(j, 0, static_cast<int64_t>(n) - 1);
    if (m[static_cast<size_t>(j)]) y.push_back(sl.y[i]);
  }
  if (y.size() > kMaxSamples) y.resize(kMaxSamples);
  Json out = {{"analysis_rate_hz", py_round_int(ofs)}, {"analysis_samples", y.size()}};
  if (y.size() < 64) return {out, {"too few on-burst samples for IQ features"}};
  const size_t N = y.size();

  fvec a(N), a2(N), a4(N);
  for (size_t i = 0; i < N; ++i) {
    a[i] = np::cabs(y[i]);
    a2[i] = a[i] * a[i];
    a4[i] = std::pow(a[i], 4.0f);  // np.power(float32, 4) -> powf (only ** 2 has a fast path, np.square)
  }
  const float amean = np::mean(a.data(), N);
  const float astd = np::stddev(a.data(), N);
  const float cv = astd / (amean + 1e-12f);
  out["envelope_cv"] = py_round(static_cast<double>(cv), 3);
  const float m2 = np::mean(a2.data(), N), m4 = np::mean(a4.data(), N);
  out["envelope_kurtosis"] = py_round(static_cast<double>(m4 / (m2 * m2 + 1e-30f)), 2);
  notes.push_back(
      "envelope CV/kurtosis on on-burst samples (Gaussian noise/OFDM: CV~0.52, kurt~2; constant envelope: CV->0, "
      "kurt->1)");

  // cyclic lines of |y|^2 -> symbol / chip / frame clocks of linear modulations
  fvec e(N);
  for (size_t i = 0; i < N; ++i) e[i] = a2[i] - m2;
  const cfvec Ec = fft::np_rfft(e.data(), N);
  fvec E(Ec.size());
  for (size_t i = 0; i < Ec.size(); ++i) E[i] = np::cabs(Ec[i]);
  const dvec fr = fft::rfftfreq(N, 1.0 / ofs);
  auto pk = peaks(E, fr, std::max(2e3, 4 * ofs / static_cast<double>(N)), 0.98 * ofs / 2);
  if (!pk.empty() && pk[0].second > 25) {
    out["cyclic_line_hz"] = py_round(pk[0].first, 1);
    out["cyclic_line_ratio"] = py_round(pk[0].second, 1);
  }
  // carrier lines of y^2 / y^4 (BPSK / QPSK signature)
  const dvec freqs = fft::fftfreq(N, 1.0 / ofs);
  cfvec yp(N);
  for (int p : {2, 4}) {
    for (size_t i = 0; i < N; ++i) {
      const cf s = np::cmul(y[i], y[i]);  // numpy y**2 / y**4 (fma-based complex multiply)
      yp[i] = p == 2 ? s : np::cmul(s, s);
    }
    const cfvec Z = fft::np_fft(yp);
    fvec Za(N);
    for (size_t i = 0; i < N; ++i) Za[i] = np::cabs(Z[i]);
    const size_t k = argmax_f(Za);
    fvec tmp = Za;
    const float r = Za[k] / (np::median(tmp.data(), N) + 1e-30f);
    if (r > 60) {
      out["x" + std::to_string(p) + "_line_ratio"] = py_round(static_cast<double>(r), 1);
      out["x" + std::to_string(p) + "_line_hz"] = py_round(freqs[k], 1);
    }
  }
  notes.push_back(
      "cyclic lines: FFT of |x|^2 (symbol clock), x^2 and x^4 (BPSK/QPSK carrier lines); peak/median > 25 (|x|^2) or "
      "> 60 (x^p) reported");

  // OFDM: cyclic-prefix autocorrelation at lag Tu, searched beyond 3/BW
  cfvec Y = fft::np_fft(y, 2 * N);
  for (auto& v : Y) {
    const float m = np::cabs(v);
    v = cf(m * m, 0.0f);
  }
  const cfvec Rc = fft::np_ifft(Y);
  fvec R(N);
  const float r0 = np::cabs(Rc[0]) + 1e-30f;
  for (size_t i = 0; i < N; ++i) R[i] = np::cabs(Rc[i]) / r0;
  const double lag_lo = std::max(3 / bw, 2e-6), lag_hi = std::min(4e-3, static_cast<double>(N) / ofs / 4);
  fvec Rs;
  dvec ls;
  for (size_t i = 0; i < N; ++i) {
    const double lag = static_cast<double>(i) / ofs;
    if (lag > lag_lo && lag < lag_hi) {
      Rs.push_back(R[i]);
      ls.push_back(lag);
    }
  }
  if (Rs.size() > 16) {
    const size_t k = argmax_f(Rs);
    const double r = Rs[k], lag = ls[k];
    fvec tmp = Rs;
    const double noise = static_cast<double>(np::median(tmp.data(), tmp.size()) + 1e-12f);
    out["acf_peak_lag_s"] = py_round(lag, 9);
    out["acf_peak"] = py_round(r, 4);
    out["acf_peak_over_median"] = py_round(r / noise, 1);
    // a real CP-OFDM signal: strong, isolated peak and >= 12 subcarriers inside the band
    if (r / noise > 15 && r > 0.06 && bw * lag >= 12 && bw >= 150e3) out["ofdm_scs_hz"] = py_round(1 / lag, 1);
  }
  notes.push_back(
      "OFDM test: normalised autocorrelation peak at lag > 3/BW; a CP-OFDM signal correlates at its useful symbol "
      "time Tu (SCS = 1/Tu)");

  // FSK: instantaneous-frequency clusters (only meaningful for constant envelope)
  if (static_cast<double>(py_round(static_cast<double>(cv), 3)) < 0.35) {
    const float k_ofs = static_cast<float>(ofs), two_pi = static_cast<float>(2 * kPi);
    dvec ifr(N - 1);
    for (size_t i = 0; i + 1 < N; ++i) {
      const cf d = y[i + 1] * std::conj(y[i]);
      ifr[i] = static_cast<double>(std::atan2(d.imag(), d.real()) * k_ofs / two_pi);
    }
    dvec tmp = ifr;
    const double lo = np::percentile(tmp.data(), tmp.size(), 10.0);
    tmp = ifr;
    const double hi = np::percentile(tmp.data(), tmp.size(), 90.0);
    out["inst_freq_p10_p90_hz"] = Json::array({py_round_int(lo), py_round_int(hi)});
    const auto h = np::histogram<double>(ifr, 64, lo - 0.2 * (hi - lo), hi + 0.2 * (hi - lo)).counts;
    const int64_t mid_min = *std::min_element(h.begin() + 28, h.begin() + 36);
    const int64_t lmax = *std::max_element(h.begin(), h.begin() + 32);
    const int64_t rmax = *std::max_element(h.begin() + 32, h.end());
    out["inst_freq_bimodal"] = lmax > 1.5 * mid_min && rmax > 1.5 * mid_min;
    notes.push_back("instantaneous-frequency 10/90th percentiles and bimodality (2-FSK test)");
  }
  return {out, notes};
}

std::pair<std::string, std::vector<std::string>> infer_modulation(const Signal& sig, const Json& f, double df_fine) {
  auto num = [&](const char* k) { return f.contains(k) && f[k].is_number() ? f[k].get<double>() : 0.0; };
  auto has = [&](const char* k) {  // python truthiness of f.get(k)
    if (!f.contains(k)) return false;
    const Json& v = f[k];
    if (v.is_null()) return false;
    if (v.is_boolean()) return v.get<bool>();
    if (v.is_number()) return v.get<double>() != 0;
    return true;
  };
  std::vector<std::string> why;
  if (sig.bw() <= 3 * df_fine + 1)
    return {"CW (unmodulated carrier)", {"occupied bandwidth " + fmt("%.0f", sig.bw()) + " Hz <= 3 fine bins"}};
  if (sig.snr_db < 3)
    return {"unknown (SNR too low)",
            {"in-band SNR " + py_float_repr(sig.snr_db) + " dB: envelope/cyclic/CP statistics not reliable"}};
  if (has("lfm_sweep_hz")) {
    why.push_back("pulses with a linear frequency sweep (" + fmt("%+.0f", num("lfm_sweep_hz") / 1e3) + " kHz over " +
                  fmt("%.0f", num("lfm_pulse_s") * 1e6) + " us, fit R^2 " + py_str(f["lfm_r2"]) +
                  ") -> LFM chirp (pulse compression)");
    return {"LFM chirp pulses", why};
  }
  const double pw = num("pulse_n") >= 2 ? num("pulse_width_s_median") : 0.0;
  if (has("ofdm_scs_hz") && pw != 0 && pw < kOfdmMinLags * num("acf_peak_lag_s")) {
    why.push_back("autocorrelation peak at " + fmt("%.2f", num("acf_peak_lag_s") * 1e6) + " us but pulses only " +
                  fmt("%.1f", pw * 1e6) + " us long -> not OFDM (pulse train)");
  } else if (has("ofdm_scs_hz")) {
    why.push_back("autocorrelation peak " + py_str(f["acf_peak"]) + " (" + py_str(f["acf_peak_over_median"]) +
                  "x median) at " + fmt("%.2f", num("acf_peak_lag_s") * 1e6) + " us -> OFDM, SCS " +
                  fmt("%.2f", num("ofdm_scs_hz") / 1e3) + " kHz");
    return {"OFDM", why};
  }
  const bool cl = has("cyclic_line_hz");
  const double clv = num("cyclic_line_hz");
  if (cl && std::fabs(clv - 1e6) < 2e3 && sig.bw() > 8e6) {
    why.push_back("|x|^2 cyclic line at " + fmt("%.4f", clv / 1e6) + " MHz (x" + py_str(f["cyclic_line_ratio"]) +
                  ") -> 1 Msym/s Barker-spread DSSS clock");
    return {"DSSS (DBPSK/DQPSK, Barker-11)", why};
  }
  if (f.contains("envelope_cv") && f["envelope_cv"].is_number() && f["envelope_cv"].get<double>() < 0.25 &&
      sig.snr_db > 10) {
    const std::string cv = py_str(f["envelope_cv"]);
    if (has("inst_freq_bimodal")) {
      why.push_back("constant envelope (CV " + cv + ") and bimodal instantaneous frequency -> 2-FSK/GFSK");
      return {"GFSK/2-FSK", why};
    }
    why.push_back("constant envelope (CV " + cv + ") -> FM / FSK / GMSK family");
    return {"constant-envelope (FM/FSK/GMSK)", why};
  }
  if (has("x2_line_ratio") && sig.snr_db > 10) {
    why.push_back("x^2 carrier line (x" + py_str(f["x2_line_ratio"]) + ") -> BPSK-like");
    return {"BPSK-like", why};
  }
  if (has("x4_line_ratio") && sig.snr_db > 10) {
    why.push_back("x^4 carrier line (x" + py_str(f["x4_line_ratio"]) + ") -> QPSK-like");
    return {"QPSK-like", why};
  }
  if (cl) {
    why.push_back("|x|^2 cyclic line at " + fmt("%.3f", clv / 1e3) +
                  " kHz -> single-carrier linear modulation, symbol rate ~" + fmt("%.1f", clv / 1e3) +
                  " ksym/s (or a multiple)");
    return {"single-carrier linear (PSK/QAM)", why};
  }
  return {"unknown", {"no decisive envelope, cyclostationary or OFDM signature"}};
}

void extract(Signal& sig, const Channelizer& ch, size_t n, double capture_s, double df_fine) {
  Json f = timing(sig, capture_s);
  auto [iq, notes] = iq_features(sig, ch, n);
  for (auto it = iq.begin(); it != iq.end(); ++it) f[it.key()] = it.value();
  if (!f["continuous"].get<bool>() && bands::pulsed_band(sig.f0, sig.f1, sig.fc())) {
    const Json pt = pulse_train(sig, ch);
    for (auto it = pt.begin(); it != pt.end(); ++it) f[it.key()] = it.value();
    notes.push_back("pulse train: 1 us envelope of the signal's band over the whole capture, pulses > median + " +
                    fmt("%g", kPulseThrDb) + " dB, gaps < " + fmt("%g", kPulseGapS * 1e6) +
                    " us joined; PRI = shortest interval cluster holding >= half of the intervals; LFM = linear "
                    "instantaneous-frequency fit on the widest pulse; pulse_emitters: core.deinterleave (PW groups, "
                    "CDIF + sequence search, stagger, jitter mode; research)");
  }
  if (f.contains("cyclic_line_hz") && f["cyclic_line_hz"].is_number() && f["cyclic_line_hz"].get<double>() != 0)
    f["symbol_rate"] = f["cyclic_line_hz"];
  auto [mod, why] = infer_modulation(sig, f, df_fine);
  for (auto it = f.begin(); it != f.end(); ++it) sig.features[it.key()] = it.value();
  sig.modulation = mod;
  if (!sig.features.contains("symbol_rate")) sig.features["symbol_rate"] = nullptr;
  sig.detection["feature_methods"] = notes;
  sig.classify_trace.insert(sig.classify_trace.end(), why.begin(), why.end());
}

}  // namespace si::features
