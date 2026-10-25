// Pipeline plugin: 315 / 433 MHz remotes and TPMS (research; python: receivers/ism_devices/plugin.py).
#include <algorithm>
#include <memory>
#include <optional>

#include "core/pyfmt.hpp"
#include "decode/emission_adapter.hpp"
#include "receivers/ism_devices/detector.hpp"

namespace si::rx::ism_devices {

namespace {

bool overlap(double lo, double hi) {
  for (const auto& [a, b] : kBands)
    if (lo < b && hi > a) return true;
  return false;
}

std::string decoder_of(const Signal& s) {
  return s.features.contains("decoder") && s.features["decoder"].is_string() ? s.features["decoder"].get<std::string>()
                                                                             : std::string();
}

class IsmDevicesDecoder : public EmissionReceiver {
 public:
  std::string name() const override { return "ism_devices"; }
  std::string package() const override { return "ism_devices"; }
  std::string protocol() const override { return "315 / 433 MHz remotes and TPMS"; }
  std::string key() const override { return "ism_devices"; }
  std::string method() const override {
    return "receivers.ism_devices: OOK pulse timing / FSK symbols -> EV1527, PT2262, KeeLoq, TPMS (research)";
  }
  std::vector<std::string> steps() const override {
    return {"candidates: 314-316 and 433.05-434.79 MHz at >= 2 MSps, 256-point frames; per frame the bin furthest "
            "over its median across frames (>= 12 dB); frequencies grouped at 100 kHz gaps -> up to 4 carriers",
            "per carrier +-100 kHz at 1 MSps: envelope (8-sample average) on above the geometric mean of its 10th / "
            "99.5th percentiles (>= 12 dB apart), runs < 30 us debounced -> pulses (mark, space)",
            "OOK: EV1527 / PT2262 (sync 1 + 31 Te, 24 PWM bits of 4 Te, repeats compared; PT2262 tri-state when "
            "every pair is 00 / 11 / 01); KeeLoq HCS (12 warm-up pulses, 10 Te header, 66 PWM bits of 3 Te: clear "
            "serial, buttons, battery, repeat; the encrypted hopping code is not reported)",
            "FSK: a mark >= 2 ms with a two-level discriminator (deviation >= 10 kHz) -> 52 us symbols -> TPMS "
            "Continental / VDO (Manchester, CRC-8 0x07 / 0xAA) and Toyota PMV-C210 (differential Manchester, "
            "CRC-8 0x07 / 0x80 + inverted pressure), both polarities",
            "nothing decoded: >= 16 OOK pulses of 50-2000 us or an FSK burst -> an unrecognised device ('-like')"};
  }

  std::vector<Emission> run_detector(DecodeContext& ctx) const override { return detect(ctx.ch); }
  bool in_range(double lo, double hi) const override { return overlap(lo, hi); }
  bool can_decode(const Signal& sig, DecodeContext&) const override { return overlap(sig.f0, sig.f1); }

  std::string group_key(const Emission& e) const override {
    std::string id = "None";
    if (e.decoded.is_object())
      for (const char* k : {"id", "code", "serial"})
        if (e.decoded.contains(k)) {
          id = e.decoded[k].get<std::string>();
          break;
        }
    return e.label + "|" + id + "|" + py_float_repr(py_round(e.center_mhz, 2));
  }

  // python IsmDevicesDecoder.arbitrate
  std::vector<Signal> arbitrate(std::vector<Signal> signals, DecodeContext&) const override {
    std::vector<size_t> hosts;
    for (size_t i = 0; i < signals.size(); ++i)
      if (decoder_of(signals[i]) == name() && signals[i].confidence >= 0.75) hosts.push_back(i);
    if (hosts.empty()) return signals;
    std::vector<size_t> keep;
    for (size_t i = 0; i < signals.size(); ++i) {
      const Signal& s = signals[i];
      std::optional<size_t> host;
      const std::string dec = decoder_of(s);
      if (!dec.empty() && dec != name() && s.decoded.status != "success") {
        for (size_t h : hosts) {
          const Signal& H = signals[h];
          bool t = false;
          for (const Box& a : s.bursts)
            for (const Box& b : H.bursts) t = t || std::min(a.t1, b.t1) > std::max(a.t0, b.t0);
          if (H.f0 <= s.fc() && s.fc() <= H.f1 && t) {
            host = h;
            break;
          }
        }
      }
      if (!host) {
        keep.push_back(i);
        continue;
      }
      Json& det = signals[*host].detection;
      if (!det.contains("explained_receiver_detections")) det["explained_receiver_detections"] = Json::array();
      det["explained_receiver_detections"].push_back({{"decoder", dec},
                                                      {"label", s.label},
                                                      {"f_low_hz", py_round_int(s.f0)},
                                                      {"f_high_hz", py_round_int(s.f1)},
                                                      {"note", "the same burst as this decoded device: attributed to it"}});
    }
    std::vector<Signal> out;
    for (size_t i : keep) out.push_back(std::move(signals[i]));
    return out;
  }

  DecodeResult summarise(const std::vector<Emission>& em, const Signal*) const override {
    const Emission* e = &em[0];
    for (const auto& x : em)
      if (x.confidence > e->confidence) e = &x;
    const Json d = e->decoded.is_object() ? e->decoded : Json::object();
    const Json& F = e->features;
    double lo = em[0].f_lo_mhz, hi = em[0].f_hi_mhz;
    for (const auto& x : em) {
      lo = std::min(lo, x.f_lo_mhz);
      hi = std::max(hi, x.f_hi_mhz);
    }
    const auto span = std::make_pair(lo * 1e6, hi * 1e6);
    auto with_base = [&](Json f) {
      f["carrier_mhz"] = F["carrier_mhz"];
      f["envelope_db"] = F["envelope_db"];
      return f;
    };
    if (e->label == "tpms") {
      Json fields = d;
      fields["frames"] = static_cast<int64_t>(em.size());
      fields = with_base(fields);
      const std::string proto = d["protocol"].get<std::string>();
      const std::string mod = proto == "Continental / VDO" ? "FSK, Manchester, 52 us" : "FSK, differential Manchester, 52 us";
      return result("success", method(), steps(), fields, "",
                    "TPMS (" + proto + ", ID " + d["id"].get<std::string>() + ")", e->confidence, mod,
                    std::to_string(em.size()) + " CRC-8 verified frame(s); " + py_str(d["pressure_kpa"]) + " kPa, " +
                        py_str(d["temperature_c"]) + " C",
                    span);
    }
    if (e->label == "ook_fixed_code") {
      const std::string kind =
          d.contains("tristate") ? "PT2262 tri-state " + d["tristate"].get<std::string>() : std::string("EV1527");
      return result("partial", method(), steps(), with_base(d),
                    "fixed-code remotes carry no checksum: the repeats are the evidence",
                    "OOK fixed-code remote (EV1527 / PT2262-like, ID " + d["id"].get<std::string>() + ")",
                    e->confidence, "OOK PWM",
                    py_str(d["repeats"]) + " identical 24-bit word(s) of " + py_str(d["words"]) + ", Te " +
                        py_str(d["te_us"]) + " us; fits " + kind,
                    span);
    }
    if (e->label == "keeloq") {
      return result("partial", method(), steps(), with_base(d),
                    "code-hopping remote: the encrypted 32-bit part is neither decoded nor reported",
                    "KeeLoq code-hopping remote-like (serial " + d["serial"].get<std::string>() + ")", e->confidence,
                    "OOK PWM",
                    py_str(d["packets"]) + " packet(s) of 12 warm-up pulses + 66 PWM bits, Te " + py_str(d["te_us"]) +
                        " us; " + py_str(d["matching_packets"]) + " with the same serial and button",
                    span);
    }
    const std::string kind = e->label == "ook_device" ? "OOK" : "FSK";
    const std::string ev = kind == "OOK" ? py_str(F["n_pulses"]) + " OOK pulses, marks " + py_str(F["mark_short_us"]) +
                                               "-" + py_str(F["mark_long_us"]) + " us"
                                         : py_str(F["fsk_bursts"]) + " FSK burst(s) >= 2 ms";
    return result("partial", method(), steps(), F, "no known protocol matched",
                  kind + " device-like (unrecognised protocol)", e->confidence, kind, ev, span);
  }
};

}  // namespace

void register_plugins(si::PluginList& out) { out.push_back(std::make_unique<IsmDevicesDecoder>()); }

}  // namespace si::rx::ism_devices
