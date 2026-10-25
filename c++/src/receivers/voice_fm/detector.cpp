// python: receivers/voice_fm/detector.py (float32 where numpy keeps the complex64 channel's
// dtype; Python-float thresholds compared in float32 as NEP 50 does).
#include "receivers/voice_fm/detector.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <numeric>
#include <set>

#include "core/fft.hpp"
#include "core/npcompat.hpp"
#include "core/pyfmt.hpp"
#include "core/signal.hpp"
#include "receivers/voice_fm/afsk.hpp"
#include "receivers/voice_fm/channels.hpp"
#include "receivers/voice_fm/subaudio.hpp"

namespace si::rx::voice_fm {

namespace {

constexpr double kRes = 1000.0;
constexpr double kMid = 48000.0;
constexpr double kFrameS = 0.005;
constexpr double kMaxDev = 6000.0;
constexpr double kKeyedDb = 10.0;
constexpr double kFskRates[2] = {2400.0, 4800.0};  // 4FSK digital voice/data symbol rates
constexpr double kFskLine = 6.0;                     // symbol-rate line over its neighbourhood: data, not speech
constexpr double kEnvCvFloor = 0.12;                 // FM envelope cv allowed whatever the SNR
constexpr double kEnvCvExcess = 1.6;                 // ... and up to this multiple of what the noise alone gives
constexpr double kAmFrame = 0.05;                    // AM: a 20 ms frame with speech (envelope rms / carrier)
constexpr double kAmMinRms = 0.04;
constexpr double kAmSpeechFrac = 0.31;               // share of the slice's envelope noise inside 300-3400 Hz

// python _symbol_line(dl, fs)["ratio"]: strongest FSK symbol-rate line in the squared derivative of the
// low-passed discriminator output over the median of its +-40..400 Hz neighbourhood
double symbol_line(const double* dl, size_t n, double fs) {
  if (n < static_cast<size_t>(0.02 * fs) || n < 2) return 0.0;
  dvec e(n - 1);
  for (size_t i = 0; i + 1 < n; ++i) {
    const double d = dl[i + 1] - dl[i];
    e[i] = d * d;
  }
  const double mu = np::mean(e);
  for (double& v : e) v = v - mu;
  const cdvec E = fft::rfft(e);
  const double val = 1.0 / (static_cast<double>(e.size()) * (1.0 / fs));
  double best = 0.0;
  for (double r : kFskRates) {
    double on = -1.0;
    dvec bg;
    for (size_t k = 0; k < E.size(); ++k) {
      const double f = static_cast<double>(k) * val;
      const double a = np::cabs(E[k]);
      if (std::fabs(f - r) <= 15) on = std::max(on, a);
      if (std::fabs(f - r) > 40 && std::fabs(f - r) < 400) bg.push_back(a);
    }
    if (on >= 0 && !bg.empty()) {
      const double q = on / (np::median(bg) + 1e-30);
      if (q > best) best = q;
    }
  }
  return best;
}

// python _noise_cv(snr_db): envelope cv of a constant-envelope carrier from noise alone, ~1/sqrt(2 SNR)
double noise_cv(double snr_db) { return 1.0 / std::sqrt(2 * std::pow(10.0, snr_db / 10)); }

// python _runs(m): (start, stop) of the True runs
std::vector<std::pair<int64_t, int64_t>> runs(const std::vector<uint8_t>& m) {
  std::vector<std::pair<int64_t, int64_t>> out;
  int64_t s = -1;
  for (size_t i = 0; i <= m.size(); ++i) {
    const bool on = i < m.size() && m[i];
    if (on && s < 0) s = static_cast<int64_t>(i);
    if (!on && s >= 0) {
      out.emplace_back(s, static_cast<int64_t>(i));
      s = -1;
    }
  }
  return out;
}

double kurtosis(dvec v) {
  const double m = np::mean(v);
  dvec v2(v.size()), v4(v.size());
  for (size_t i = 0; i < v.size(); ++i) {
    v[i] = v[i] - m;
    v2[i] = v[i] * v[i];
    v4[i] = std::pow(v[i], 4.0);
  }
  const double s = np::mean(v2);
  return np::mean(v4) / (s * s + 1e-30);
}

std::optional<Json> artifact_of(const Measurement& m, const std::vector<Measurement>& meas, size_t n_stronger,
                                double fc) {
  const double fo = m.f_hz - fc;
  for (size_t i = 0; i < n_stronger; ++i) {
    const Measurement& s = meas[i];
    const double so = s.f_hz - fc;
    if (std::fabs(so) < 1e3) continue;
    for (int k : {-1, 3, -3, 5, -5, 7, -7, 9, -9, 11}) {
      if (std::fabs(fo - k * so) >= 300 + 60 * std::abs(k)) continue;
      bool ok;
      if (k == -1) {
        ok = s.peak_db - m.peak_db >= 15;
      } else {
        const double r = m.deviation_hz / std::max(s.deviation_hz, 1.0);
        ok = 0.6 * std::abs(k) <= r && r <= 1.4 * std::abs(k);
      }
      if (ok) return Json{{"of_mhz", py_round(s.f_hz / 1e6, 5)}, {"k", k}};
    }
  }
  return std::nullopt;
}

bool truthy(const Json& j, const char* k) {
  if (!j.contains(k)) return false;
  const Json& v = j[k];
  if (v.is_null()) return false;
  if (v.is_boolean()) return v.get<bool>();
  if (v.is_number()) return v.get<double>() != 0;
  if (v.is_string()) return !v.get<std::string>().empty();
  return !v.empty();
}

// python _noise_fm_rms(snr_db): speech-band discriminator rms of an unmodulated carrier from noise alone
double noise_fm_rms(double snr_db) {
  const double nd = 1.0 / (2 * std::pow(10.0, snr_db / 10)) / 20e3;
  return std::sqrt(2 * nd * (std::pow(3400.0, 3) - std::pow(300.0, 3)) / 3);
}

// python _am_emission(m, ppm): full-carrier AM voice in the aeronautical bands (research, S13)
std::optional<Emission> am_emission(const Measurement& m, double ppm) {
  const Json air = air_channel(m.f_hz, ppm);
  if (air.empty() || !m.snr_db || *m.snr_db < 10 || m.fsk_line_ratio >= kFskLine) return std::nullopt;
  const double snr = *m.snr_db;
  const double noise_am = noise_cv(snr) * std::sqrt(kAmSpeechFrac);
  if (m.am_voice_rms < std::max(kAmMinRms, 3 * noise_am) || m.am_activity < 0.1) return std::nullopt;
  if (m.voice_rms_hz > 3 * noise_fm_rms(snr) + 150) return std::nullopt;   // the phase moves: not AM
  double conf = 0.7;
  Json why = Json::array({"carrier with speech-band envelope modulation",
                          fmt("envelope speech rms %.0f %% of the carrier (noise alone %.1f %%)", m.am_voice_rms * 100,
                              noise_am * 100),
                          fmt("phase quiet: discriminator speech-band rms %.0f Hz", m.voice_rms_hz),
                          fmt("speech in %.0f %% of 20 ms frames", m.am_activity * 100), fmt("keyed %.2f s", m.keyed_s)});
  if (air.contains("emergency"))
    why.push_back(air["emergency"]);
  else
    why.push_back(air["service"].get<std::string>() + ", channel " + py_str(air["channel_mhz"]) + " MHz (" +
                  fmt("%+lld", static_cast<long long>(air["offset_hz"].get<int64_t>())) + " Hz)");
  if (m.am_activity >= 0.2) conf += 0.1;
  Json F = Json::object();
  F["center_mhz"] = py_round(m.f_hz / 1e6, 6);
  F["modulation_rms"] = py_round(m.am_voice_rms, 3);
  F["channel_bw_khz"] = 25.0;
  F["keyed_s"] = np::round(m.keyed_s, 3);
  F["keyed_fraction"] = np::round(m.keyed_s / m.duration_s, 3);
  F["envelope_cv"] = py_round(m.envelope_cv, 3);
  F["voice_activity"] = py_round(m.am_activity, 3);
  F["phase_rms_hz"] = py_round_int(m.voice_rms_hz);
  F["snr_db"] = py_round(snr, 1);
  F["evidence"] = why;
  for (auto it = air.begin(); it != air.end(); ++it) F[it.key()] = it.value();
  if (!m.artifacts.empty()) F["suppressed_artifacts"] = m.artifacts;
  Json segs = Json::array();
  for (const auto& [a, b] : m.segs)
    segs.push_back(Json::array({np::round(static_cast<double>(a) * kFrameS, 3), np::round(static_cast<double>(b) * kFrameS, 3)}));
  F["segments_s"] = segs;
  Emission e;
  e.family = "voice_fm";
  e.label = "am_voice";
  e.start_us = static_cast<double>(m.segs.front().first) * kFrameS * 1e6;
  e.end_us = static_cast<double>(m.segs.back().second) * kFrameS * 1e6;
  e.f_lo_mhz = (m.f_hz - 4e3) / 1e6;
  e.f_hi_mhz = (m.f_hz + 4e3) / 1e6;
  e.center_mhz = m.f_hz / 1e6;
  e.snr_db = F["snr_db"].get<double>();
  e.confidence = py_round(std::min(conf, 0.9), 2);
  e.decoded = nullptr;
  e.features = std::move(F);
  return e;
}

// python _packet_emission(m): AX.25 / APRS packets (CRC-valid) carried by the FM channel (research, S19)
std::optional<Emission> packet_emission(const Measurement& m) {
  if (m.packets.empty()) return std::nullopt;
  bool aprs = false;
  std::set<std::string> srcs;
  for (const Json& p : m.packets) {
    aprs = aprs || p.contains("aprs");
    srcs.insert(p["src"].get<std::string>());
  }
  Json st = Json::array();
  for (const auto& x : srcs) st.push_back(x);
  Json F = Json::object();
  F["center_mhz"] = py_round(m.f_hz / 1e6, 6);
  F["deviation_hz"] = py_round_int(m.deviation_hz);
  F["channel_bw_khz"] = m.deviation_hz <= 3200 ? 12.5 : 25.0;
  F["keyed_s"] = np::round(m.keyed_s, 3);
  F["snr_db"] = m.snr_db ? Json(py_round(*m.snr_db, 1)) : Json(nullptr);
  F["stations"] = st;
  F["packets"] = Json(m.packets);
  F["evidence"] = Json::array({fmt("%zu AX.25 UI frame(s) with a valid CRC-16/X.25 (Bell 202 AFSK 1200 bit/s)",
                                   m.packets.size()),
                               "stations " + py_repr(st)});
  Json segs = Json::array();
  for (const auto& [a, b] : m.segs)
    segs.push_back(Json::array({np::round(static_cast<double>(a) * kFrameS, 3), np::round(static_cast<double>(b) * kFrameS, 3)}));
  F["segments_s"] = segs;
  const double half = m.deviation_hz <= 3200 ? 6.25e3 : 12.5e3;
  Emission e;
  e.family = "voice_fm";
  e.label = aprs ? "aprs" : "ax25_packet";
  e.start_us = static_cast<double>(m.segs.front().first) * kFrameS * 1e6;
  e.end_us = static_cast<double>(m.segs.back().second) * kFrameS * 1e6;
  e.f_lo_mhz = (m.f_hz - half) / 1e6;
  e.f_hi_mhz = (m.f_hz + half) / 1e6;
  e.center_mhz = m.f_hz / 1e6;
  if (!F["snr_db"].is_null()) e.snr_db = F["snr_db"].get<double>();
  e.confidence = 0.97;
  e.decoded = Json{{"crc_ok", true}, {"packets", static_cast<int64_t>(m.packets.size())}};
  e.features = std::move(F);
  return e;
}

std::optional<Emission> emission(const Measurement& m, double ppm) {
  if (auto pk = packet_emission(m)) return pk;
  if (auto am = am_emission(m, ppm)) return am;
  const Json sub = m.sub ? *m.sub : Json::object();
  if (m.envelope_cv > 0.35 || m.deviation_hz > kMaxDev) return std::nullopt;
  const bool has_sub = !sub.empty();
  if (!has_sub && m.snr_db && m.envelope_cv > std::max(kEnvCvFloor, kEnvCvExcess * noise_cv(*m.snr_db)))
    return std::nullopt;                                    // envelope varies beyond noise: two tones / AM
  const bool modulated = has_sub || m.voice_activity >= 0.05 || m.deviation_hz >= 300;
  if (!modulated) return std::nullopt;
  if (!has_sub && m.voice_kurtosis < 2.2 && m.voice_rms_hz > 500) return std::nullopt;
  if (m.fsk_line_ratio >= kFskLine) return std::nullopt;   // FSK symbol-rate line: digital, not analog
  const std::optional<Json> chan = match_channel(m.f_hz, ppm);
  double conf = 0.5;
  Json why = Json::array({"constant-envelope narrowband FM carrier", fmt("peak deviation %.0f Hz", m.deviation_hz),
                          fmt("keyed %.2f s", m.keyed_s)});
  const std::string type = sub.contains("type") ? sub["type"].get<std::string>() : "";
  if (type == "DCS") {
    conf += 0.3;
    why.push_back("DCS " + sub["dcs_code"].get<std::string>() + " (" + py_str(sub["dcs_words_verified"]) +
                  " Golay-verified words)");
  } else if (truthy(sub, "ctcss_hz")) {
    conf += 0.25;
    why.push_back("CTCSS " + py_str(sub["ctcss_hz"]) + " Hz (+" + py_str(sub["ctcss_snr_db"]) + " dB)");
  } else if (has_sub) {
    conf += 0.1;
    why.push_back("sub-audible tone " + py_str(sub["ctcss_measured_hz"]) + " Hz (non-standard)");
  }
  if (m.voice_activity >= 0.1) {
    conf += 0.1;
    why.push_back(fmt("speech-band activity %.0f%% of 20 ms frames", m.voice_activity * 100));
  }
  if (chan) {
    conf += 0.1;
    why.push_back((*chan)["service"].get<std::string>() + " channel " + py_str((*chan)["channel"]) + " (" +
                  fmt("%+lld", static_cast<long long>((*chan)["offset_hz"].get<int64_t>())) + " Hz)");
  }
  if (m.keyed_s >= 0.3) conf += 0.05;
  conf = std::min(conf, 0.97);
  double chan_bw = chan && truthy(*chan, "channel_bw_khz") ? (*chan)["channel_bw_khz"].get<double>()
                                                            : (m.deviation_hz <= 3200 ? 12.5 : 25.0);
  const double half = std::min(chan_bw * 1e3 / 2, m.deviation_hz + 3400);
  Json F = Json::object();
  F["center_mhz"] = py_round(m.f_hz / 1e6, 6);
  F["deviation_hz"] = py_round_int(m.deviation_hz);
  F["channel_bw_khz"] = chan_bw;
  F["keyed_s"] = np::round(m.keyed_s, 3);
  F["keyed_fraction"] = np::round(m.keyed_s / m.duration_s, 3);
  F["envelope_cv"] = py_round(m.envelope_cv, 3);
  F["voice_activity"] = py_round(m.voice_activity, 3);
  F["voice_rms_hz"] = py_round_int(m.voice_rms_hz);
  F["voice_kurtosis"] = py_round(m.voice_kurtosis, 2);
  F["snr_db"] = m.snr_db ? Json(py_round(*m.snr_db, 1)) : Json(nullptr);
  F["evidence"] = why;
  if (chan)
    for (auto it = chan->begin(); it != chan->end(); ++it) F[it.key()] = it.value();
  for (auto it = sub.begin(); it != sub.end(); ++it) F[it.key()] = it.value();
  if (!m.artifacts.empty()) F["suppressed_artifacts"] = m.artifacts;
  Json segs = Json::array();
  for (const auto& [a, b] : m.segs)
    segs.push_back(Json::array({np::round(static_cast<double>(a) * kFrameS, 3), np::round(static_cast<double>(b) * kFrameS, 3)}));
  F["segments_s"] = segs;
  std::string label = "nbfm_voice";
  if (chan) {
    std::string svc = (*chan)["service"].get<std::string>();
    for (auto& c : svc) c = c == '/' ? '_' : static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    label += "_" + svc + "_ch" + py_str((*chan)["channel"]);
  }
  Emission e;
  e.family = "voice_fm";
  e.label = label;
  e.start_us = static_cast<double>(m.segs.front().first) * kFrameS * 1e6;
  e.end_us = static_cast<double>(m.segs.back().second) * kFrameS * 1e6;
  e.f_lo_mhz = (m.f_hz - half) / 1e6;
  e.f_hi_mhz = (m.f_hz + half) / 1e6;
  e.center_mhz = m.f_hz / 1e6;
  if (!F["snr_db"].is_null()) e.snr_db = F["snr_db"].get<double>();
  e.confidence = py_round(conf, 2);
  e.decoded = has_sub ? sub : Json(nullptr);
  e.features = std::move(F);
  return e;
}

}  // namespace

const std::vector<std::pair<double, double>>& voice_ranges() {
  static const std::vector<std::pair<double, double>> r = {{118e6, 174e6}, {225e6, 520e6}, {806e6, 870e6}};
  return r;
}

Json air_channel(double f_hz, double ppm) {
  struct Band {
    double a, b;
    const char* name;
  };
  static const Band bands[] = {{118e6, 137e6, "VHF air band"}, {225e6, 400e6, "UHF military air band"}};
  const double f = f_hz * (1 + ppm * 1e-6);
  Json out = Json::object();
  for (const Band& bd : bands)
    if (bd.a <= f && f < bd.b) {
      out["service"] = bd.name;
      const double nom = static_cast<double>(py_round_int(f / 25e3)) * 25e3;
      out["channel_mhz"] = py_round(nom / 1e6, 4);
      out["offset_hz"] = py_round_int(f - nom);
      break;
    }
  if (std::fabs(f - 121.5e6) <= 5e3) out["emergency"] = "aeronautical emergency 121.5 MHz";
  if (std::fabs(f - 243.0e6) <= 5e3) out["emergency"] = "military emergency 243.0 MHz";
  return out;
}

std::vector<Candidate> candidates(const Channelizer& ch, size_t max_n) {
  const double fs = ch.fs(), fc = ch.fc();
  const size_t L = ch.L();
  const size_t g = static_cast<size_t>(std::max<int64_t>(1, py_round_int(kRes * static_cast<double>(L) / fs)));
  const size_t start = (L + 1) / 2;   // np.argsort(fftfreq): negative frequencies first
  fvec P(L);
  dvec f(L);
  for (size_t i = 0; i < L; ++i) {
    const size_t j = (start + i) % L;
    P[i] = np::cabs2(ch.X()[j]);
    f[i] = ch.f()[j];
  }
  const size_t n = L / g;
  fvec Pg(n);
  dvec fg(n);
  std::vector<uint8_t> inr(n);
  size_t ninr = 0;
  fvec inner;
  for (size_t j = 0; j < n; ++j) {
    Pg[j] = np::mean(P.data() + j * g, g);
    fg[j] = np::mean(f.data() + j * g, g);
    const double fa = fc + fg[j];
    bool in = false;
    for (const auto& [a, b] : voice_ranges()) in |= fa >= a && fa <= b;
    const bool centre = std::fabs(fg[j]) < 0.45 * fs;
    inr[j] = in && centre;
    ninr += inr[j];
    if (centre) inner.push_back(Pg[j]);
  }
  if (ninr < 20) return {};
  const double floor = static_cast<double>(np::median(inner));
  const float floorf = static_cast<float>(floor);
  fvec Ldb(n);
  for (size_t j = 0; j < n; ++j) Ldb[j] = 10.0f * std::log10(Pg[j] / floorf + 1e-30f);
  sig::PeakOpts o;
  o.height_min = kKeyedDb;
  o.prominence = true;
  o.prominence_min = kKeyedDb;
  o.width = true;
  o.width_max = 40e3 / kRes;
  o.rel_height = 0.5;
  o.distance = std::max(1.0, std::floor(3e3 / kRes));
  const sig::Peaks pk = sig::find_peaks(Ldb, o);
  std::vector<Candidate> out;
  for (size_t q = 0; q < pk.peaks.size(); ++q) {
    const int64_t i = pk.peaks[q];
    if (!inr[i]) continue;
    const int64_t a = std::max<int64_t>(0, i - 4), b = std::min<int64_t>(static_cast<int64_t>(n), i + 5);
    dvec fp;
    fvec p;
    for (int64_t k = a; k < b; ++k) {
      fp.push_back(fg[k] * static_cast<double>(Pg[k]));
      p.push_back(Pg[k]);
    }
    out.push_back({fc + np::sum(fp) / static_cast<double>(np::sum(p)), static_cast<double>(Ldb[i]),
                   pk.widths[q] * kRes, floor / (fs * static_cast<double>(L))});
  }
  std::stable_sort(out.begin(), out.end(), [](const Candidate& x, const Candidate& y) { return -x.peak_db < -y.peak_db; });
  if (out.size() > max_n) out.resize(max_n);
  return out;
}

std::optional<Measurement> measure(const Channelizer& ch, const Candidate& c) {
  const Channelizer::Slice sl = ch.slice(c.f_hz, 10e3, kMid);
  const cfvec& y = sl.y;
  const double fs_y = sl.fs;
  const size_t n_fr = std::max<size_t>(1, static_cast<size_t>(kFrameS * fs_y));
  const size_t nf = y.size() / n_fr;
  if (nf < 4) return std::nullopt;
  fvec p(nf);
  {
    fvec a(n_fr);
    for (size_t i = 0; i < nf; ++i) {
      for (size_t j = 0; j < n_fr; ++j) a[j] = np::cabs2(y[i * n_fr + j]);
      p[i] = np::mean(a);
    }
  }
  const double noise = c.noise_density * 20e3;
  const float thr = static_cast<float>(noise * std::pow(10.0, kKeyedDb / 10));
  std::vector<uint8_t> keyed(nf), notk(nf);
  for (size_t i = 0; i < nf; ++i) {
    keyed[i] = p[i] > thr;
    notk[i] = !keyed[i];
  }
  std::vector<uint8_t> k = keyed;
  for (const auto& [a, b] : runs(notk))
    if (0 < a && b < static_cast<int64_t>(nf) && static_cast<double>(b - a) * kFrameS < 0.03)
      for (int64_t i = a; i < b; ++i) k[i] = 1;
  std::vector<std::pair<int64_t, int64_t>> segs;
  for (const auto& r : runs(k))
    if (static_cast<double>(r.second - r.first) * kFrameS >= 0.01) segs.push_back(r);
  const double dur = static_cast<double>(y.size()) / fs_y;
  int64_t nk = 0;
  for (const auto& [a, b] : segs) nk += b - a;
  const double keyed_s = static_cast<double>(nk) * kFrameS;
  if (segs.empty() || keyed_s < std::min(0.04, 0.5 * dur)) return std::nullopt;
  std::vector<size_t> idx;
  for (const auto& [a, b] : segs)
    for (int64_t i = a * static_cast<int64_t>(n_fr); i < b * static_cast<int64_t>(n_fr); ++i) idx.push_back(i);
  fvec env(idx.size());
  for (size_t i = 0; i < idx.size(); ++i) env[i] = np::cabs(y[idx[i]]);
  const float cv = np::stddev(env) / (np::mean(env) + 1e-30f);
  const float sc = static_cast<float>(fs_y / (2 * kPi));
  fvec d(y.size());
  // numpy temporary elision: in `y[1:] * np.conj(y[:-1])` the conj() temporary, when it holds
  // >= 256 KiB, receives the product in place, i.e. numpy computes conj(prev) * cur; the
  // fma-based complex multiply is not symmetric in its last bit, so the order is mirrored.
  const bool elide = (y.size() - 1) * sizeof(cf) >= 256 * 1024;
  for (size_t i = 1; i < y.size(); ++i) {
    const cf z = elide ? np::cmul(std::conj(y[i - 1]), y[i]) : np::cmul(y[i], std::conj(y[i - 1]));
    d[i] = std::atan2(z.imag(), z.real()) * sc;
  }
  if (y.size() > 1) d[0] = d[1];
  const dvec dl = sig::sosfiltfilt(sig::butter_sos(4, 4000.0, sig::BType::Lowpass, fs_y), d);
  dvec dk(idx.size());
  for (size_t i = 0; i < idx.size(); ++i) dk[i] = dl[idx[i]];
  const double mean_off = np::mean(dk);
  dvec ad(dk.size());
  for (size_t i = 0; i < dk.size(); ++i) ad[i] = std::fabs(dk[i] - mean_off);
  const double dev = np::percentile(ad, 99.5);
  const dvec vb = sig::sosfiltfilt(sig::butter_sos(4, dvec{300, 3400}, sig::BType::Bandpass, fs_y), d);
  dvec v(idx.size());
  for (size_t i = 0; i < idx.size(); ++i) v[i] = vb[idx[i]];
  const size_t nfrm = static_cast<size_t>(0.02 * fs_y);
  const size_t nv = nfrm ? v.size() / nfrm : 0;
  double activity = 0.0;
  if (nv) {
    size_t loud = 0;
    dvec sq(nfrm);
    for (size_t i = 0; i < nv; ++i) {
      for (size_t j = 0; j < nfrm; ++j) sq[j] = v[i * nfrm + j] * v[i * nfrm + j];
      loud += std::sqrt(np::mean(sq)) > 250;
    }
    activity = static_cast<double>(loud) / static_cast<double>(nv);
  }
  // AM: the envelope's speech band relative to the carrier (python ea / ae / arms)
  const double env_mean = static_cast<double>(np::mean(env));
  double am_rms = 0.0, am_act = 0.0;
  {
    const float em = static_cast<float>(env_mean);
    dvec e0(y.size());
    for (size_t i = 0; i < y.size(); ++i) e0[i] = static_cast<double>(np::cabs(y[i]) - em);
    const dvec eb = sig::sosfiltfilt(sig::butter_sos(4, dvec{300, 3400}, sig::BType::Bandpass, fs_y), e0);
    dvec ae(idx.size()), a2(idx.size());
    for (size_t i = 0; i < idx.size(); ++i) {
      ae[i] = eb[idx[i]] / (env_mean + 1e-30);
      a2[i] = ae[i] * ae[i];
    }
    am_rms = std::sqrt(np::mean(a2));
    if (nv) {
      size_t loud = 0;
      dvec sq(nfrm);
      for (size_t i = 0; i < nv; ++i) {
        for (size_t j = 0; j < nfrm; ++j) sq[j] = ae[i * nfrm + j] * ae[i * nfrm + j];
        loud += std::sqrt(np::mean(sq)) > kAmFrame;
      }
      am_act = static_cast<double>(loud) / static_cast<double>(nv);
    }
  }
  size_t lg = 0;
  for (size_t i = 1; i < segs.size(); ++i)
    if (segs[i].second - segs[i].first > segs[lg].second - segs[lg].first) lg = i;
  const fvec dseg(d.begin() + segs[lg].first * static_cast<int64_t>(n_fr), d.begin() + segs[lg].second * static_cast<int64_t>(n_fr));
  Measurement m;
  m.sub = subaudible(dseg, fs_y);
  {
    const size_t a0 = static_cast<size_t>(segs[lg].first) * n_fr, a1 = static_cast<size_t>(segs[lg].second) * n_fr;
    m.fsk_line_ratio = symbol_line(dl.data() + a0, std::min(a1, dl.size()) - std::min(a0, dl.size()), fs_y);
  }
  fvec pk;
  for (size_t i = 0; i < nf; ++i)
    if (keyed[i]) pk.push_back(p[i]);
  if (!pk.empty()) m.snr_db = static_cast<double>(10.0f * std::log10(np::median(pk) / static_cast<float>(noise)));
  m.f_hz = c.f_hz + mean_off;
  m.segs = segs;
  m.keyed_s = keyed_s;
  m.duration_s = dur;
  m.envelope_cv = static_cast<double>(cv);
  m.deviation_hz = dev;
  dvec v2(v.size());
  for (size_t i = 0; i < v.size(); ++i) v2[i] = v[i] * v[i];
  m.voice_rms_hz = std::sqrt(np::mean(v2));
  m.voice_kurtosis = kurtosis(v);
  m.voice_activity = activity;
  m.am_voice_rms = am_rms;
  m.am_activity = am_act;
  m.packets = afsk::decode(d, fs_y);
  m.peak_db = c.peak_db;
  return m;
}

std::vector<Emission> detect(const Channelizer& ch, double ppm, const dvec& exclude_hz) {
  const double fs = ch.fs(), fc = ch.fc();
  bool any = false;
  for (const auto& [a, b] : voice_ranges()) any |= fc - fs / 2 < b && fc + fs / 2 > a;
  if (!any) return {};
  std::vector<Measurement> measured;
  for (const Candidate& c : candidates(ch)) {
    bool skip = false;                                       // the receiver's own clock lines
    for (double f : exclude_hz) skip = skip || std::fabs(c.f_hz - f) < 2e3;
    if (skip) continue;
    if (auto m = measure(ch, c)) measured.push_back(std::move(*m));
  }
  std::vector<size_t> kept;
  for (size_t i = 0; i < measured.size(); ++i) {
    const auto art = artifact_of(measured[i], measured, i, fc);
    if (!art) {
      kept.push_back(i);
      continue;
    }
    for (size_t j = 0; j < i; ++j)
      if (py_round(measured[j].f_hz / 1e6, 5) == (*art)["of_mhz"].get<double>()) {
        measured[j].artifacts.push_back(
            {{"f_mhz", py_round(measured[i].f_hz / 1e6, 5)}, {"of_mhz", (*art)["of_mhz"]}, {"k", (*art)["k"]}});
        break;
      }
  }
  std::vector<Emission> out;
  for (size_t i : kept)
    if (auto e = emission(measured[i], ppm)) out.push_back(std::move(*e));
  return out;
}

std::vector<Emission> detect(const cfvec& x, double fs, double fc, double ppm) {
  bool any = false;
  for (const auto& [a, b] : voice_ranges()) any |= fc - fs / 2 < b && fc + fs / 2 > a;
  if (!any) return {};
  return detect(Channelizer(x, fs, fc), ppm);
}

}  // namespace si::rx::voice_fm
