#include "classify/rules.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "core/bands_india.hpp"
#include "core/pyfmt.hpp"

namespace si::classify {

using bands::band_keys;
using bands::ble_adv;
using bands::describe;
using bands::pulsed_bands;
using bands::gsm_arfcn;
using bands::lookup;

const double kMinConf = 0.5;

namespace {

constexpr double kHz = 1e3;
using bands::MHz;

struct Cand {
  std::string label;
  double score;
  std::vector<std::string> why;
};
using Rule = std::function<std::optional<Cand>(const Signal&, const Json&)>;

std::optional<double> opt(const Json& f, const char* k) {
  if (!f.contains(k) || !f[k].is_number()) return std::nullopt;
  return f[k].get<double>();
}
bool truthy(const Json& f, const char* k) {
  if (!f.contains(k)) return false;
  const Json& v = f[k];
  if (v.is_null()) return false;
  if (v.is_boolean()) return v.get<bool>();
  if (v.is_number()) return v.get<double>() != 0;
  return true;
}
double num(const Json& f, const char* k, double dflt) {
  auto v = opt(f, k);
  return v ? *v : dflt;
}

bool near(std::optional<double> v, std::initializer_list<double> targets, double tol) {
  if (!v) return false;
  for (double t : targets)
    if (std::fabs(*v - t) <= tol * t) return true;
  return false;
}

std::set<std::string> keys(const Signal& s) {
  std::set<std::string> k;
  for (const auto& x : band_keys(s.f0, s.f1)) k.insert(x);
  for (const auto* b : lookup(s.fc())) k.insert(b->key);
  return k;
}

bool any_of(const std::set<std::string>& k, std::initializer_list<const char*> names) {
  for (const char* n : names)
    if (k.count(n)) return true;
  return false;
}

std::optional<Cand> r_cw(const Signal& s, const Json&) {
  if (s.modulation.rfind("CW", 0) == 0)
    return Cand{"CW carrier (unmodulated tone)-like", 0.55,
                {"bandwidth " + fmt("%.0f", s.bw()) + " Hz: unmodulated carrier or receiver spur (sweep check decides)"}};
  return std::nullopt;
}

std::optional<Cand> r_wifi_ofdm(const Signal& s, const Json& f) {
  const auto k = keys(s);
  if (!any_of(k, {"ism_24", "unii_12", "unii_2e", "ism_58", "its_59", "wifi6e"})) return std::nullopt;
  double sc = 0.25;
  std::vector<std::string> why{"inside a delicensed Wi-Fi band"};
  if (near(opt(f, "ofdm_scs_hz"), {312.5e3}, 0.05)) {
    sc += 0.35;
    why.push_back("OFDM subcarrier spacing 312.5 kHz (802.11a/g/n/ac)");
  } else if (truthy(f, "ofdm_scs_hz")) {
    return std::nullopt;
  }
  if (near(s.bw(), {16.6e6, 17.5e6, 36e6, 38e6, 76e6, 78e6}, 0.12)) {
    sc += 0.15;
    why.push_back("occupied bandwidth " + fmt("%.1f", s.bw() / MHz) + " MHz matches 20/40/80 MHz channel");
  }
  if (!truthy(f, "continuous") && num(f, "n_bursts", 0) >= 1) {
    sc += 0.05;
    why.push_back("bursty (CSMA)");
  }
  if (sc >= 0.5) return Cand{"Wi-Fi 802.11 OFDM-like", sc, why};
  return std::nullopt;
}

std::optional<Cand> r_wifi_dsss(const Signal& s, const Json&) {
  if (keys(s).count("ism_24") && s.modulation.rfind("DSSS", 0) == 0)
    return Cand{"Wi-Fi 802.11b DSSS-like", 0.75, {"2.4 GHz ISM", "1 Msym/s Barker DSSS clock"}};
  return std::nullopt;
}

std::pair<std::string, std::string> cell_band(const std::set<std::string>& k) {
  const std::vector<std::pair<std::string, std::vector<std::pair<const char*, const char*>>>> tables = {
      {"dl", {{"b28_dl", "B28/n28"}, {"b5_dl", "B5/n5"}, {"b8_dl", "B8/n8"}, {"b3_dl", "B3/n3"}, {"b1_dl", "B1/n1"}}},
      {"tdd", {{"b40", "B40/n40"}, {"b41", "B41/n41"}, {"n78", "n78"}}},
      {"ul", {{"b28_ul", "B28"}, {"b5_ul", "B5"}, {"b8_ul", "B8"}, {"b3_ul", "B3"}, {"b1_ul", "B1"}}}};
  for (const auto& [kind, tab] : tables)
    for (const auto& [key, name] : tab)
      if (k.count(key)) return {kind, name};
  return {"", ""};
}

std::optional<Cand> r_cellular_ofdm(const Signal& s, const Json& f) {
  const auto [kind, band] = cell_band(keys(s));
  if (kind.empty()) return std::nullopt;
  const auto scs = opt(f, "ofdm_scs_hz");
  std::string K = kind;
  std::transform(K.begin(), K.end(), K.begin(), ::toupper);
  double sc = 0.25;
  std::vector<std::string> why{"inside India IMT band " + band + " (" + K + ")"};
  std::string tech;
  if (near(scs, {15e3}, 0.04)) {
    sc += 0.3;
    why.push_back("OFDM SCS 15 kHz (LTE / NR-15)");
    tech = "LTE";
  } else if (near(scs, {30e3}, 0.04)) {
    sc += 0.3;
    why.push_back("OFDM SCS 30 kHz (NR)");
    tech = "5G NR";
  } else if (truthy(f, "ofdm_scs_hz")) {
    return std::nullopt;
  }
  if (near(s.bw(), {1.4e6, 3e6, 5e6, 10e6, 15e6, 20e6, 40e6, 50e6, 60e6, 80e6, 100e6}, 0.12)) {
    sc += 0.1;
    why.push_back("bandwidth " + fmt("%.1f", s.bw() / MHz) + " MHz on a 3GPP channel size");
  }
  const auto per = opt(f, "burst_period_s");
  if (kind == "dl" && truthy(f, "continuous")) {
    sc += 0.1;
    why.push_back("always on (downlink reference signals)");
  }
  if (kind == "tdd" && per && *per != 0) {
    for (double p : {2.5e-3, 5e-3, 10e-3})
      if (std::fabs(*per - p) < 0.1e-3) {
        sc += 0.15;
        why.push_back("TDD frame periodicity " + fmt("%.2f", *per * 1e3) + " ms");
        break;
      }
  }
  if (kind == "ul" && !truthy(f, "continuous")) {
    sc += 0.05;
    why.push_back("bursty (UE uplink scheduling)");
  }
  if (tech.empty()) {
    tech = "LTE/NR";
    if (sc < 0.45) return std::nullopt;
  }
  const std::string lab = kind == "dl" ? tech + "-like downlink " + band
                          : kind == "tdd" ? tech + "-like TDD " + band
                                          : tech + "-like uplink " + band;
  return Cand{lab, sc, why};
}

std::optional<Cand> r_gsm(const Signal& s, const Json& f) {
  const auto k = keys(s);
  if (!any_of(k, {"b8_dl", "b3_dl", "b8_ul", "b3_ul"})) return std::nullopt;
  if (!(120e3 <= s.bw() && s.bw() <= 320e3)) return std::nullopt;
  double sc = 0.35;
  std::vector<std::string> why{"GSM band", "bandwidth " + fmt("%.0f", s.bw() / kHz) + " kHz (GSM 200 kHz carrier)"};
  if (s.modulation.find("constant") != std::string::npos || s.modulation.find("GMSK") != std::string::npos) {
    sc += 0.2;
    why.push_back("constant envelope (GMSK)");
  }
  const bool dl = any_of(k, {"b8_dl", "b3_dl"});
  if (truthy(f, "continuous") && dl) {
    sc += 0.1;
    why.push_back("continuous (BCCH carrier)");
  }
  const auto arfcn = gsm_arfcn(s.fc());
  return Cand{std::string("GSM-like ") + (dl ? "downlink" : "uplink") + (arfcn ? " ARFCN " + std::to_string(*arfcn) : ""),
              sc, why};
}

std::optional<Cand> r_umts(const Signal& s, const Json& f) {
  if (keys(s).count("b1_dl") && 3.5e6 <= s.bw() && s.bw() <= 5.2e6 && !truthy(f, "ofdm_scs_hz") &&
      truthy(f, "continuous")) {
    double sc = 0.5;
    std::vector<std::string> why{"B1 downlink", "bandwidth " + fmt("%.2f", s.bw() / MHz) + " MHz (3.84 Mcps WCDMA)",
                                 "continuous, no OFDM CP"};
    if (near(opt(f, "cyclic_line_hz"), {3.84e6}, 0.01)) {
      sc += 0.2;
      why.push_back("3.84 MHz chip-rate cyclic line");
    }
    return Cand{"UMTS/WCDMA-like downlink B1", sc, why};
  }
  return std::nullopt;
}

std::optional<Cand> r_ble(const Signal& s, const Json& f) {
  if (!keys(s).count("ism_24") || !(0.6e6 <= s.bw() && s.bw() <= 2.6e6)) return std::nullopt;
  if (truthy(f, "continuous") || num(f, "burst_dur_s_median", 1) > 2.2e-3) return std::nullopt;
  double sc = 0.35;
  std::vector<std::string> why{"2.4 GHz ISM", "bandwidth " + fmt("%.2f", s.bw() / MHz) + " MHz, bursts " +
                                                  fmt("%.0f", num(f, "burst_dur_s_median", 0) * 1e6) + " us"};
  if (s.modulation.find("FSK") != std::string::npos || s.modulation.find("constant") != std::string::npos) {
    sc += 0.2;
    why.push_back("GFSK / constant envelope");
  }
  for (const auto& [ch, fc] : ble_adv())
    if (std::fabs(s.fc() - fc) < 0.3e6) {
      sc += 0.1;
      why.push_back("centred on BLE advertising channel " + std::to_string(ch));
    }
  return Cand{"Bluetooth (BLE/BR)-like", sc, why};
}

std::optional<Cand> r_adsb(const Signal& s, const Json& f) {
  const double dur = num(f, "burst_dur_s_median", 1);
  if (std::fabs(s.fc() - 1090e6) < 1.5e6 && dur < 150e-6)
    return Cand{"Mode S / ADS-B (1090 MHz)-like", 0.7,
                {"1090 MHz SSR reply frequency", "short bursts (" + fmt("%.0f", dur * 1e6) + " us)"}};
  if (std::fabs(s.fc() - 1030e6) < 1.5e6 && dur < 60e-6)
    return Cand{"SSR interrogation (1030 MHz)-like", 0.65, {"1030 MHz SSR interrogation", "short pulses"}};
  return std::nullopt;
}

constexpr double kPulseTbMin = 1.0;                       // bandwidth x pulse width of a real pulse
constexpr double kDmePairS[2] = {10e-6, 45e-6};          // pulse-pair group width (X 12 us / Y 30 us mode)

// python _radar_class: radar class from the band plan, pulse width, PRI and LFM (research, signal gaps S10)
struct RadarClass {
  std::string label;
  bool cls = false;
  std::vector<std::string> extra;
};

template <class K>
RadarClass radar_class(const Signal& s, const Json& f, const K& k, bool dme, bool pulses, double dur) {
  if (dme) return {"DME/TACAN pulse pairs-like", false, {}};
  const bool has_pri = truthy(f, "pri_s");
  const double pri = num(f, "pri_s", 0);
  const double fc = s.fc();
  const bool lfm = truthy(f, "lfm_sweep_hz");
  const bool long_lfm = lfm && num(f, "lfm_pulse_s", 0) >= 10e-6;
  if (k.count("radar_s_atc")) {
    if (long_lfm)
      return {"Pulsed radar (S-band ATC, solid-state)-like", true,
              {"S-band ATC band (2.7-2.9 GHz), long LFM pulse (pulse compression)"}};
    if (pulses && dur > 1.5e-6 && has_pri && 0.77e-3 <= pri && pri <= 4e-3)
      return {"Pulsed radar (weather, S-band)-like", true, {"S-band weather-radar band (2.7-2.9 GHz)"}};
    return {"Pulsed radar (S-band ATC)-like", true, {"S-band ATC band (2.7-2.9 GHz)"}};
  }
  if (k.count("radar_s_marine")) {
    if (dur <= 1.2e-6) return {"Pulsed radar (marine S-band)-like", true, {"maritime navigation-radar band (2.9-3.1 GHz)"}};
    return {"Pulsed radar (S-band)-like", true, {"S-band radar (2.9-3.1 GHz)"}};
  }
  if (5570e6 <= fc && fc <= 5660e6)
    return {"Pulsed radar (weather, C-band)-like", true, {"C-band weather-radar band (5.6 GHz)"}};
  if (k.count("radar_5400") || k.count("unii_2e"))
    return {"Pulsed radar (C-band)-like", true, {"C-band radiolocation (5.25-5.85 GHz)"}};
  if (k.count("rnss_l2_radar") || k.count("arns_radar") || k.count("fixed_radar_1350")) {
    if (1270e6 <= fc && fc <= 1295e6 && has_pri && pri <= 250e-6)
      return {"Wind profiler radar-like", true, {"L-band wind-profiler range (~1280 MHz)"}};
    if (pulses && dur >= 1e-6 && has_pri && 1e-3 <= pri && pri <= 4e-3)
      return {"Pulsed radar (L-band surveillance)-like", true, {"L-band surveillance-radar band (1.21-1.4 GHz)"}};
    return {"Pulsed radar (L-band)-like", true, {"L-band radar (1.21-1.4 GHz)"}};
  }
  if (440e6 <= fc && fc <= 450e6 && has_pri && pri <= 250e-6)
    return {"Wind profiler radar-like", true, {"UHF wind-profiler range (440-450 MHz)"}};
  return {"Pulsed radar-like", false, {}};
}

std::optional<Cand> r_pulsed(const Signal& s, const Json& f) {
  const auto k = keys(s);
  bool in_band = false;
  for (const auto& x : k) in_band = in_band || pulsed_bands().count(x);
  if (!in_band) return std::nullopt;
  for (double f0 : {1030e6, 1090e6})
    if (std::fabs(s.fc() - f0) < 1.5e6) return std::nullopt;  // SSR interrogations / Mode S replies (r_adsb)
  // pulse widths from the 1 us envelope when measured (spectrogram burst boxes are far longer than the pulses)
  const bool pulses = num(f, "pulse_n", 0) >= 2;
  const double dur = pulses ? num(f, "pulse_width_s_median", 0) : num(f, "burst_dur_s_median", 1);
  if (truthy(f, "continuous") || dur > 200e-6) return std::nullopt;
  // DME replies come in pulse pairs, radars repeat on their PRI: one blip is not
  if ((pulses ? num(f, "pulse_n", 0) : num(f, "n_bursts", 0)) < 2) return std::nullopt;
  // a pulse of width tau occupies >~ 1/tau: shorter "pulses" are noise peaks
  if (pulses && s.bw() * dur < kPulseTbMin) return std::nullopt;
  // shared with Wi-Fi: bursts there are packets unless they repeat on a PRI (DFS radar)
  if (k.count("unii_2e") && !(truthy(f, "pri_s") || truthy(f, "burst_period_s"))) return std::nullopt;
  const bool dme = k.count("arns_l") || k.count("rnss_l5");
  // a DME reply is a pulse pair, 12 (X) or 30 (Y) us apart: one group 15-35 us
  if (dme && pulses && !(kDmePairS[0] <= dur && dur <= kDmePairS[1])) return std::nullopt;
  const RadarClass rc = radar_class(s, f, k, dme, pulses, dur);
  const std::string name = rc.label;
  double sc = 0.5;
  std::vector<std::string> why{describe(s.fc()), "short pulses (" + fmt("%.1f", dur * 1e6) + " us" +
                                                     (pulses ? ", " + fmt("%d", static_cast<int>(num(f, "pulse_n", 0))) +
                                                                   " pulse groups, 1 us envelope)"
                                                             : std::string(")"))};
  for (const auto& w : rc.extra) why.push_back(w);
  if (truthy(f, "pri_s")) {  // intervals within +-15 % of the shortest cluster
    sc += 0.15;
    const double spread = num(f, "pri_spread", 0);
    why.push_back("PRI " + fmt("%.0f", num(f, "pri_s", 0) * 1e6) + " us" +
                  (spread > 0.02 ? " (staggered, spread " + fmt("%.0f", spread * 100) + " %)" : std::string()));
  } else if (truthy(f, "burst_period_s")) {
    sc += 0.15;
    why.push_back("regular PRI " + fmt("%.0f", num(f, "burst_period_s", 0) * 1e6) + " us");
  }
  if (truthy(f, "lfm_sweep_hz")) {
    sc += 0.1;
    why.push_back("LFM chirp " + fmt("%+.0f", num(f, "lfm_sweep_hz", 0) / 1e3) + " kHz over " +
                  fmt("%.0f", num(f, "lfm_pulse_s", 0) * 1e6) + " us (pulse compression)");
  }
  if (f.contains("pulse_emitters") && f["pulse_emitters"].is_array() && !f["pulse_emitters"].empty()) {
    // core/deinterleave (research, signal gaps S10): evidence only
    std::string t;
    for (const Json& e : f["pulse_emitters"]) {
      std::string one = "PRI " + fmt("%.1f", e["pri_us"].get<double>()) + " us";
      if (e["stagger_us"].is_array()) {
        std::string st;
        for (const Json& v : e["stagger_us"]) st += (st.empty() ? "" : "/") + fmt("%.0f", v.get<double>());
        one += " staggered " + st + " us";
      }
      if (e["jitter"].get<double>() >= 0.01) one += ", jitter " + fmt("%.1f", e["jitter"].get<double>() * 100) + " %";
      one += ", " + std::to_string(e["n_pulses"].get<int64_t>()) + " pulses, PW " + fmt("%.1f", e["pw_us"].get<double>()) + " us";
      t += (t.empty() ? "" : "; ") + one;
    }
    why.push_back("deinterleaved into " + std::to_string(f["pulse_emitters"].size()) + " emitter(s): " + t);
  }
  if (rc.cls) sc = std::min(0.8, sc + 0.1);
  return Cand{name, sc, why};
}

std::optional<Cand> r_dvbt2(const Signal& s, const Json& f) {
  const auto k = keys(s);
  if (!k.count("uhf_tv") && !k.count("uhf_bcast_imt")) return std::nullopt;
  if (6.5e6 <= s.bw() && s.bw() <= 8e6 && truthy(f, "continuous")) {
    double sc = 0.45;
    std::vector<std::string> why{"UHF TV band", "bandwidth " + fmt("%.2f", s.bw() / MHz) + " MHz (8 MHz raster)",
                                 "continuous"};
    const auto scs = opt(f, "ofdm_scs_hz");
    if (near(scs, {279e3 / 64 * 1, 1116, 2232, 4464, 558, 279}, 0.05)) {
      sc += 0.35;
      why.push_back("OFDM SCS " + fmt("%.0f", *scs) + " Hz (DVB-T2 FFT mode)");
    } else if (truthy(f, "ofdm_scs_hz")) {
      sc += 0.15;
      why.push_back("OFDM SCS " + fmt("%.0f", *scs) + " Hz");
    }
    return Cand{"DVB-T2 TV (Doordarshan)-like", sc, why};
  }
  return std::nullopt;
}

std::optional<Cand> r_narrowband(const Signal& s, const Json&) {
  const auto k = keys(s);
  if (!any_of(k, {"govt_mobile_300", "govt_mobile_335", "tetra_ppdr", "govt_mobile_406", "pmr_450", "govt_mobile_435",
                  "trunk_800", "trunk_850", "metaids"}))
    return std::nullopt;
  if (!(3e3 <= s.bw() && s.bw() <= 30e3)) return std::nullopt;
  if (k.count("tetra_ppdr") && 18e3 <= s.bw() && s.bw() <= 27e3)
    return Cand{"TETRA-like (25 kHz pi/4-DQPSK)", 0.55, {"380-400 MHz public-safety band", "~25 kHz"}};
  if (k.count("metaids"))
    return Cand{"Radiosonde-like telemetry", 0.5, {"400-406 MHz MetAids", fmt("%.1f", s.bw() / kHz) + " kHz"}};
  const std::string m = (s.modulation.find("constant") != std::string::npos || s.modulation.find("FSK") != std::string::npos)
                            ? "analog FM / FSK voice-data"
                            : "narrowband";
  return Cand{"PMR land-mobile (" + m + ")-like", 0.55, {describe(s.fc()), fmt("%.1f", s.bw() / kHz) + " kHz channel"}};
}

std::optional<Cand> r_srd(const Signal& s, const Json& f) {
  const auto k = keys(s);
  if ((k.count("ism_433") || k.count("ism_865")) && s.bw() < 600e3 && !truthy(f, "continuous"))
    return Cand{"ISM short-range device (433/865 MHz)-like", 0.55,
                {describe(s.fc()), "bursty, " + fmt("%.0f", s.bw() / kHz) + " kHz"}};
  return std::nullopt;
}

std::string join(const std::vector<std::string>& v, const char* sep) {
  std::string out;
  for (size_t i = 0; i < v.size(); ++i) out += (i ? sep : "") + v[i];
  return out;
}

std::optional<Cand> r_iridium(const Signal& s, const Json& f) {
  if (!(keys(s).count("mss_iridium") && 25e3 <= s.bw() && s.bw() <= 60e3 && !truthy(f, "continuous")))
    return std::nullopt;
  // python: Iridium needs a PSK signature (x^2 / x^4 carrier line or the 25 kHz symbol line)
  std::vector<std::string> psk;
  if (truthy(f, "x2_line_ratio")) psk.push_back("x^2 carrier line");
  if (truthy(f, "x4_line_ratio")) psk.push_back("x^4 carrier line");
  if (f.contains("cyclic_line_hz") && f["cyclic_line_hz"].is_number() &&
      std::fabs(f["cyclic_line_hz"].get<double>() - 25e3) <= 0.05 * 25e3)
    psk.push_back("25 kHz symbol line");
  if (psk.empty()) return std::nullopt;
  return Cand{"Iridium-like TDMA burst", 0.55,
              {"1616-1626.5 MHz Iridium band", fmt("%.0f", s.bw() / kHz) + " kHz burst", join(psk, " + ") + " (PSK)"}};
}

const std::vector<Rule>& rules() {
  static const std::vector<Rule> r = {r_cw,  r_wifi_ofdm, r_wifi_dsss, r_cellular_ofdm, r_gsm,  r_umts,   r_ble,
                                      r_adsb, r_pulsed,   r_dvbt2,     r_narrowband,    r_srd, r_iridium};
  return r;
}

}  // namespace

void classify(Signal& s) {
  const Json& f = s.features;
  if (s.has_flag("lo_leakage")) {
    s.label = "LO leakage (receiver artifact)";
    s.confidence = 0.95;
    s.modulation = "CW";
    s.evidence = "narrow tone exactly at the tuned centre frequency (DC bin): AD936x LO feed-through";
    return;
  }
  if (s.has_flag("receiver_spur")) {
    const double fz = num(f, "spur_freq_hz", s.fc());
    if (f.contains("spur_kind") && f["spur_kind"] == "sample_clock_harmonic") {
      const int64_t k = f["spur_harmonic"].get<int64_t>();
      const double rate = num(f, "spur_sample_rate_hz", 0);
      s.label = fmt("Receiver spur (sample-clock harmonic %lld x %g MHz)", static_cast<long long>(k), rate / MHz);
      s.confidence = 0.9;
      s.evidence = fmt("narrow line at %.6f MHz = %lld x the sample rate within %.0f Hz: generated by the "
                       "receiver's own clock",
                       fz / MHz, static_cast<long long>(k), std::fabs(fz - static_cast<double>(k) * rate));
    } else {
      const double ppm = num(f, "spur_receiver_ppm", 0);
      s.label = "Receiver spur (locked to the receiver reference)";
      s.confidence = 0.8;
      s.evidence = fmt("narrow line at %.6f MHz, on the integer-MHz grid of the receiver's own frequency scale "
                       "although the receiver error (",
                       fz / MHz) +
                   py_float_repr(ppm) +
                   fmt(" ppm) would move an external line %.1f kHz off it", std::fabs(ppm) * fz * 1e-6 / 1e3);
    }
    s.modulation = "CW";
    return;
  }
  if (s.has_flag("iq_image")) {
    const double src = num(f, "image_of_f_center_hz", 0);
    s.label = "IQ image (receiver artifact)";
    s.confidence = 0.8;
    s.evidence = "mirror of a >=15 dB stronger signal at " + fmt("%.3f", src / MHz) +
                 " MHz about the LO, same burst timing: AD936x IQ-imbalance image";
    return;
  }
  std::vector<Cand> cands;
  for (const auto& r : rules())
    if (auto c = r(s, f)) cands.push_back(*c);
  std::stable_sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) { return -a.score < -b.score; });
  Json cj = Json::array();
  for (size_t i = 0; i < cands.size() && i < 4; ++i)
    cj.push_back({{"label", cands[i].label}, {"score", py_round(cands[i].score, 2)}});
  s.features["classifier_candidates"] = cj;
  const std::string band = describe(s.fc());
  std::vector<std::string> trace2(s.classify_trace.begin(),
                                  s.classify_trace.begin() + std::min<size_t>(2, s.classify_trace.size()));
  if (!cands.empty() && cands[0].score >= kMinConf) {
    s.label = cands[0].label;
    s.confidence = py_round(std::min(cands[0].score, 0.9), 2);
    std::vector<std::string> ev = cands[0].why;
    ev.insert(ev.end(), trace2.begin(), trace2.end());
    s.evidence = join(ev, "; ");
  } else {
    s.label = "unknown";
    s.confidence = cands.empty() ? 0.2 : py_round(cands[0].score, 2);
    s.evidence = "band prior: " + band + ". Features: " + s.modulation + ", bw " + fmt("%.0f", s.bw() / kHz) +
                 " kHz, duty " + py_str(f.value("duty_cycle", Json(nullptr))) +
                 (cands.empty() ? std::string("; no rule matched")
                                : "; best candidate '" + cands[0].label + "' only " + fmt("%.2f", cands[0].score)) +
                 (s.classify_trace.empty() ? std::string() : "; " + join(trace2, "; "));
  }
}

}  // namespace si::classify
