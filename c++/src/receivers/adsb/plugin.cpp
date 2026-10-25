// Pipeline plugin: ADS-B / Mode S replies as capture-level signals (python: receivers/adsb/plugin.py).
#include <memory>
#include <set>

#include "core/npcompat.hpp"
#include "core/pyfmt.hpp"
#include "decode/plugin.hpp"
#include "receivers/adsb/messages.hpp"
#include "receivers/adsb/modeac.hpp"
#include "receivers/adsb/modes.hpp"
#include "receivers/adsb/tracks.hpp"
#include "receivers/adsb/uat.hpp"

namespace si::rx::adsb {
namespace {

// (lat, lon[, alt_m]) from SigMF core:geolocation (GeoJSON point), if the capture has one
Json receiver_position(const Capture& cap) {
  if (!cap.meta.is_object() || !cap.meta.contains("global") || !cap.meta["global"].is_object()) return nullptr;
  const Json& g0 = cap.meta["global"];
  if (!g0.contains("core:geolocation")) return nullptr;
  const Json& g = g0["core:geolocation"];
  if (!g.is_object() || !g.contains("type") || g["type"] != "Point" || !g.contains("coordinates")) return nullptr;
  const Json& c = g["coordinates"];
  if (!c.is_array() || c.size() < 2) return nullptr;
  Json r = Json::array({c[1], c[0]});
  if (c.size() > 2) r.push_back(c[2]);
  return r;
}

class AdsbDecoder : public Decoder {
 public:
  std::string name() const override { return "adsb"; }
  std::string package() const override { return "adsb"; }
  std::string protocol() const override { return "Mode S / ADS-B 1090 MHz"; }
  bool capture_level() const override { return true; }

  bool can_decode(const Signal& sig, DecodeContext&) const override {
    return sig.f0 - 1e6 <= kFAdsb && kFAdsb <= sig.f1 + 1e6;
  }

  std::vector<Signal> capture_signals(DecodeContext& ctx) const override {
    const Capture& cap = ctx.cap;
    if (!(cap.fc - 0.45 * cap.fs < kFAdsb && kFAdsb < cap.fc + 0.45 * cap.fs)) return {};
    const auto sl = ctx.ch.slice(kFAdsb, 2e6, kFs);
    fvec env(sl.y.size());
    for (size_t i = 0; i < env.size(); ++i) env[i] = np::cabs(sl.y[i]);
    const float noise = np::median(env);
    const float thr = 6.0f * noise * static_cast<float>(kSpu) / 2.0f * 2.0f;   // pulses ~>6x noise amplitude
    std::vector<Json> msgs = demod(env, thr);
    fields(msgs);
    std::vector<std::pair<double, double>> busy;
    for (const Json& m : msgs) {
      const double ts = m["t_s"].get<double>();
      busy.emplace_back(ts, ts + static_cast<double>(8 + m["bits"].get<int>()) * 1e-6);
    }
    const std::vector<Json> replies = modeac::decode(env, kFs, busy, noise);
    std::vector<std::string> steps = {
        "channelise 1090 MHz +-2 MHz at 8 MSps, magnitude envelope",
        "8 us preamble correlation (pulses 0/1/3.5/4.5 us), PPM bit decisions",
        "CRC-24 (0xFFF409): DF17/18 residual 0; other DFs address/parity",
        "fields of every DF/typecode (pyModeS + DO-260B bit layouts): identity, category, altitude, "
        "velocity, intent, integrity, capabilities, Comm-B registers, ACAS and radar reply fields",
        "per-aircraft: CPR position (even/odd global, then local), BDS50/60 disambiguation, "
        "wind and air temperature, interrogation environment, plausibility flags"};
    const std::vector<std::string> steps_ac = {
        "SSR Mode A/C replies (research): 0.45 us pulses >= 12 dB over the envelope median, framing pulses "
        "F1 / F2 20.3 us apart within 3 dB, every pulse between them on the 1.45 us grid (else garbled), "
        "X empty; SPI at 24.65 us",
        "each reply read both ways: Mode A identity (octal ABCD) and Mode C Gillham altitude when the "
        "bits are a valid Gillham code; words repeated in the block are the verified ones"};
    if (msgs.empty()) {
      ctx.cache.set<std::string>("adsb_note", "no Mode S frame passed CRC in this capture");
      if (replies.empty()) return {};
      std::vector<std::string> st = {steps[0]};
      st.insert(st.end(), steps_ac.begin(), steps_ac.end());
      return {mode_ac_signal(replies, st)};
    }
    std::set<std::string> icaos;
    for (const Json& m : msgs)
      if (m.contains("icao") && m["icao"].is_string() && !m["icao"].get<std::string>().empty())
        icaos.insert(m["icao"].get<std::string>());
    const double t0 = cap.t0 ? static_cast<double>(cap.t0->epoch_s * 1000000 + cap.t0->micro) / 1e6 : 0.0;
    for (Json& m : msgs) m["t"] = t0 + m["t_s"].get<double>();
    const Json ref = receiver_position(cap);
    Json aircraft = tracks(msgs, ref);
    Json dfs = Json::object();
    for (const Json& m : msgs) {
      const std::string k = std::to_string(m["df"].get<int>());
      dfs[k] = dfs.contains(k) ? dfs[k].get<int64_t>() + 1 : 1;
    }
    Json f = Json::object();
    f["n_frames"] = msgs.size();
    f["dfs"] = dfs;
    f["icao"] = Json::array();
    for (const auto& s : icaos) f["icao"].push_back(s);
    f["aircraft"] = aircraft;
    f["receiver_position"] = ref;
    f["frames"] = Json::array();
    std::vector<std::pair<double, double>> bursts;
    for (const Json& m : msgs) {
      const double ts = m["t_s"].get<double>();
      bursts.emplace_back(ts, ts + static_cast<double>(8 + m["bits"].get<int>()) * 1e-6);
    }
    for (const Json& r : replies) {
      const double ts = r["t_s"].get<double>();
      bursts.emplace_back(ts, ts + modeac::kSpiUs * 1e-6);
    }
    if (!replies.empty()) {
      f["mode_ac"] = modeac::summarise(replies);
      steps.insert(steps.end(), steps_ac.begin(), steps_ac.end());
    }
    Json names = Json::array(), uav = Json::array();
    for (const Json& a : aircraft) {
      const bool cs = a.contains("callsign") && a["callsign"].is_string() && !a["callsign"].get<std::string>().empty();
      names.push_back(cs ? a["callsign"] : a["icao"]);
      if (a.contains("uav") && a["uav"] == true) uav.push_back(cs ? a["callsign"] : a["icao"]);
    }
    const size_t n = msgs.size();
    for (auto& m : msgs) f["frames"].push_back(std::move(m));
    DecodeResult res = result("success", "Mode S/ADS-B receiver (this package) + pyModeS", steps, f, "",
                              std::string("Mode S / ADS-B 1090 MHz"), 0.97, std::string("PPM 1 Mb/s (pulse position)"),
                              fmt("%zu Mode S frame(s) with valid parity, aircraft ", n) + py_repr(names) +
                                  (uav.empty() ? std::string() : "; UAV (emitter category B6): " + py_repr(uav)) +
                                  (replies.empty() ? std::string() : fmt("; %zu Mode A/C repl(ies)", replies.size())),
                              std::make_pair(kFAdsb - 1e6, kFAdsb + 1e6));
    return {receiver_signal(bursts, kFAdsb - 1e6, kFAdsb + 1e6, res, std::nullopt, name())};
  }

 private:
  // Mode A/C replies without any Mode S frame: their own signal (research; python _mode_ac_signal)
  Signal mode_ac_signal(const std::vector<Json>& replies, const std::vector<std::string>& steps) const {
    Json f = modeac::summarise(replies);
    std::vector<Json> rep;
    std::set<std::string> special;
    for (const Json& g : f["replies"])
      if (g["repeated"].get<bool>()) {
        rep.push_back(g);
        if (g.contains("special_code"))
          special.insert(g["code"].get<std::string>() + " " + g["special_code"].get<std::string>());
      }
    Json codes = Json::array();
    const Json src = rep.empty() ? f["replies"] : Json(rep);
    for (size_t i = 0; i < src.size() && i < 6; ++i) codes.push_back(src[i]["code"]);
    std::string label = rep.empty() ? "SSR Mode A/C replies-like" : "SSR Mode A/C replies";
    if (!special.empty()) {
      std::string j;
      for (const auto& s : special) j += (j.empty() ? "" : ", ") + s;
      label += " (squawk " + j + ")";
    }
    Json fields = Json::object();
    fields["mode_ac"] = f;
    DecodeResult res = result("partial", "Mode A/C reply decoder (this package; research)", steps, fields,
                              "Mode A/C replies carry no parity: codes are verified only by repetition", label,
                              rep.empty() ? 0.6 : 0.8, std::string("pulse code (SSR Mode A/C)"),
                              fmt("%lld reply(ies) on the 1.45 us grid between framing pulses 20.3 us apart, %lld "
                                  "distinct, %lld repeated; codes ",
                                  static_cast<long long>(f["n_replies"].get<int64_t>()),
                                  static_cast<long long>(f["n_distinct"].get<int64_t>()),
                                  static_cast<long long>(f["n_repeated"].get<int64_t>())) + py_repr(codes),
                              std::make_pair(kFAdsb - 1e6, kFAdsb + 1e6));
    std::vector<std::pair<double, double>> bursts;
    for (const Json& r : replies) {
      const double ts = r["t_s"].get<double>();
      bursts.emplace_back(ts, ts + modeac::kSpiUs * 1e-6);
    }
    return receiver_signal(bursts, kFAdsb - 1e6, kFAdsb + 1e6, res, std::nullopt, name());
  }
};

// UAT 978 MHz ADS-B downlink (research, PLAN G14): capture-level, like Mode S (python UatDecoder)
class UatDecoder : public Decoder {
 public:
  std::string name() const override { return "uat"; }
  std::string package() const override { return "adsb"; }
  std::string protocol() const override { return "UAT 978 MHz ADS-B"; }
  bool capture_level() const override { return true; }

  bool can_decode(const Signal& sig, DecodeContext&) const override {
    return sig.f0 - 1e6 <= uat::kFUat && uat::kFUat <= sig.f1 + 1e6;
  }

  std::vector<Signal> capture_signals(DecodeContext& ctx) const override {
    const Capture& cap = ctx.cap;
    if (!(cap.fc - 0.45 * cap.fs < uat::kFUat && uat::kFUat < cap.fc + 0.45 * cap.fs)) return {};
    const auto sl = ctx.ch.slice(uat::kFUat, 0.8e6, uat::kFs);
    std::vector<Json> msgs = uat::demod(sl.y);
    if (msgs.empty()) {
      ctx.cache.set<std::string>("uat_note", "no UAT downlink message passed Reed-Solomon in this capture");
      return {};
    }
    const std::vector<std::string> steps = {
        "channelise 978 MHz +-0.8 MHz at 4 samples / bit (4.1667 MSps)",
        "phase change per bit, 36-bit sync 0xEACDDA4E2 (<= 4 bit errors) at the best timing phase, slicer "
        "centred on the sync's levels (CFO)",
        "Reed-Solomon GF(256) 0x187, fcr 120: Long RS(48,34) when the payload type is not 0, else Basic "
        "RS(30,18) (as FlightAware dump978)",
        "DO-282B fields: header, state vector, Mode Status (emitter category, callsign)"};
    Json aircraft = Json::array();
    for (const Json& m : msgs) {
      const std::string addr = m["address"].get<std::string>();
      Json* a = nullptr;
      for (auto& x : aircraft)
        if (x["address"] == addr) a = &x;
      if (!a) {
        aircraft.push_back({{"address", addr}, {"qualifier", m["address_qualifier"]}, {"n", 0}});
        a = &aircraft.back();
      }
      (*a)["n"] = (*a)["n"].get<int>() + 1;
      for (const char* k : {"lat", "lon", "alt_ft", "callsign", "emitter_category", "uav"})
        if (m.contains(k) && !m[k].is_null()) (*a)[k] = m[k];
    }
    Json names = Json::array(), uav = Json::array();
    for (const Json& a : aircraft) {
      const bool cs = a.contains("callsign") && a["callsign"].is_string() && !a["callsign"].get<std::string>().empty();
      names.push_back(cs ? a["callsign"] : a["address"]);
      if (a.contains("uav") && a["uav"] == true) uav.push_back(cs ? a["callsign"] : a["address"]);
    }
    std::vector<std::pair<double, double>> bursts;
    for (const Json& m : msgs) {
      const double ts = m["t_s"].get<double>();
      bursts.emplace_back(ts, ts + (36.0 + (m["frame"] == "long" ? 8 * 48 : 8 * 30)) / uat::kRate);
    }
    Json f = Json::object();
    f["n_messages"] = msgs.size();
    f["aircraft"] = aircraft;
    f["messages"] = Json(msgs);
    DecodeResult res = result("success", "UAT 978 receiver (this package; research)", steps, f, "",
                              std::string("UAT 978 MHz ADS-B") + (uav.empty() ? "" : ", UAV"), 0.95,
                              std::string("CPFSK 1.041667 Mb/s, h 0.6"),
                              fmt("%zu RS-corrected UAT message(s), aircraft ", msgs.size()) + py_repr(names) +
                                  (uav.empty() ? std::string() : "; UAV (emitter category 14): " + py_repr(uav)),
                              std::make_pair(uat::kFUat - 0.8e6, uat::kFUat + 0.8e6));
    return {receiver_signal(bursts, uat::kFUat - 0.8e6, uat::kFUat + 0.8e6, res, std::nullopt, name())};
  }
};

}  // namespace

void register_plugins(PluginList& out) {
  out.push_back(std::make_unique<AdsbDecoder>());
  out.push_back(std::make_unique<UatDecoder>());
}

}  // namespace si::rx::adsb
