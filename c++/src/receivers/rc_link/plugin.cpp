// Pipeline plugin: LoRa/ELRS, FSK RC links, BLE adverts and 2 Mb/s ESB-like links (python:
// receivers/rc_link/plugin.py).
#include <algorithm>
#include <cmath>
#include <memory>
#include <set>

#include "core/bands_india.hpp"
#include "core/drone_models.hpp"
#include "decode/emission_adapter.hpp"
#include "core/mavlink.hpp"
#include "core/npcompat.hpp"
#include "receivers/rc_link/models.hpp"
#include "receivers/rc_link/protocols.hpp"
#include "receivers/rc_link/detector.hpp"
#include "receivers/rc_link/label_check.hpp"

namespace si::rx::rc_link {
namespace {

// BLE advertising channel centre +- (receiver ppm + carrier offset), channels 2 MHz apart (python ADV_TOL_HZ)
constexpr double kAdvTolHz = 0.35e6;
// sub-GHz LoRa links outside the delicensed bands (python LORA_HZ): EU868, FCC915 drone links and off-band
// telemetry radios; a capture that touches no RC band but this range gets the LoRa scan only
constexpr double kLoraLo = 860e6, kLoraHi = 930e6;
bool lora_overlap(double lo, double hi) { return lo < kLoraHi && hi > kLoraLo; }
const char* const kEsbLabel = "2.4 GHz GFSK link (ESB-like, 2 Mb/s)";

std::string replace_all(std::string s, const std::string& a, const std::string& b) {
  for (size_t p = s.find(a); p != std::string::npos; p = s.find(a, p + b.size())) s.replace(p, a.size(), b);
  return s;
}

bool truthy(const Json& d, const char* k) {
  if (!d.is_object() || !d.contains(k)) return false;
  const Json& v = d[k];
  if (v.is_boolean()) return v.get<bool>();
  if (v.is_string()) return !v.get<std::string>().empty();
  if (v.is_number()) return v.get<double>() != 0;
  return !v.is_null() && !v.empty();
}

class RcLinkDecoder : public EmissionReceiver {
 public:
  std::string name() const override { return "rc_link"; }
  std::string package() const override { return "rc_link"; }
  std::string protocol() const override { return "LoRa/FSK/BLE/ESB"; }
  std::string key() const override { return "rc_link"; }
  std::set<std::string> bands() const override { return {"ism_24", "ism_433", "ism_865", "b8_gap"}; }
  std::string method() const override { return "receivers.rc_link (research receiver)"; }
  std::vector<std::string> steps() const override {
    return {"LoRa: blind dechirp scan per sub-band (critically sampled, K-symbol coherent FFT), 2.25-downchirp "
            "SFD verification, symbol walk -> SF/BW/preamble/payload length, ELRS channel-grid match",
            "GFSK/FSK: matched-area energy at 3 widths, narrowband flank test, constant-envelope test, "
            "inst.-frequency clusters (deviation), symbol rate from transition timing",
            "1 Mb/s h~0.5 bursts: BLE preamble + access address 0x8E89BED6, de-whitening, CRC-24",
            "2.4 GHz 0.5-3 MHz bursts: 2 Mb/s GFSK slicing, >= 16-bit 0101 preamble -> ESB-like packet; packets "
            "sharing a sync word: CRC-16/CCITT end and init found jointly (>= 3 distinct packets on one CRC)"};
  }

  bool rc_capture(const DecodeContext& ctx) const {
    return EmissionReceiver::in_range(ctx.cap.fc - ctx.cap.fs / 2, ctx.cap.fc + ctx.cap.fs / 2);
  }
  bool in_range(double lo, double hi) const override {
    return EmissionReceiver::in_range(lo, hi) || lora_overlap(lo, hi);
  }
  bool can_decode(const Signal& sig, DecodeContext& ctx) const override {
    if (rc_capture(ctx)) return EmissionReceiver::can_decode(sig, ctx);
    return lora_overlap(ctx.cap.fc - ctx.cap.fs / 2, ctx.cap.fc + ctx.cap.fs / 2) && lora_overlap(sig.f0, sig.f1);
  }

  std::vector<Emission> run_detector(DecodeContext& ctx) const override {
    const cfvec& x = ctx.cap.iq();
    const cfvec* X = ctx.ch.L() == x.size() ? &ctx.ch.X() : nullptr;
    if (rc_capture(ctx)) return detect(x, ctx.cap.fs, ctx.cap.fc, X);
    const std::pair<double, double> r{kLoraLo - ctx.cap.fc, kLoraHi - ctx.cap.fc};
    return detect(x, ctx.cap.fs, ctx.cap.fc, X, &r);
  }

  // python group_key(): ESB channels sit on the 1 MHz raster (2400 + k MHz)
  std::string group_key(const Emission& e) const override {
    if (e.family == "esb_link") return e.label + "|" + std::to_string(py_round_int(e.center_mhz));
    return EmissionReceiver::group_key(e);
  }

  // python pretty(): by family, not substring ("lora_preamble" contains "ble")
  std::string pretty(const std::string& top, const std::vector<Emission>& em) const override {
    for (const auto& e : em)
      if (e.family == "esb_link") return kEsbLabel;
    bool ble = false;
    for (const auto& e : em)
      if (e.family == "bluetooth" || truthy(e.decoded, "aa")) ble = true;
    if (ble) {
      for (const auto& e : em)
        if (truthy(e.decoded, "crc_ok")) {
          const Json& n = e.decoded.contains("pdu_name") ? e.decoded["pdu_name"] : Json(nullptr);
          const std::string pn = n.is_string() ? n.get<std::string>() : std::string("None");
          const bool lr = e.decoded.contains("phy") && e.decoded["phy"].is_string() &&
                          e.decoded["phy"].get<std::string>().rfind("LE Coded", 0) == 0;
          const std::string bl = lr ? "Bluetooth LE Long Range" : "Bluetooth LE";
          for (const auto& x : em)
            if (x.decoded.is_object() && x.decoded.contains("ad") && x.decoded["ad"].contains("remote_id"))
              return "Drone Remote ID (" + bl + ", " + pn + ")";
          for (const auto& x : em)   // item tracker (research, signal gaps S7)
            if (truthy(x.decoded, "crc_ok") && x.decoded.contains("ad") && x.decoded["ad"].contains("device") &&
                x.decoded["ad"]["device"]["tracker"].get<bool>())
              return "Bluetooth LE tracker (" + x.decoded["ad"]["device"]["label"].get<std::string>() + ")";
          return bl + " advertising (" + pn + ")";
        }
      for (const auto& e : em)   // BCH-verified BR/EDR access code (research, signal gaps S7)
        if (e.decoded.is_object() && e.decoded.contains("lap"))
          return "Bluetooth BR/EDR (LAP " + e.decoded["lap"].get<std::string>() +
                 (e.decoded["inquiry"].get<bool>() ? ", inquiry" : "") + ")";
      // advertising uses only channels 37/38/39; a 1 Mb/s GFSK hop elsewhere is a BLE data channel or a BR/EDR
      // hop (the measured modulation index does not separate h 0.5 from 0.32 reliably)
      for (const auto& e : em)
        for (const auto& [ch, f] : bands::ble_adv())
          if (std::fabs(e.center_mhz * 1e6 - f) <= kAdvTolHz) return "Bluetooth LE advertising-like";
      return "Bluetooth (BLE/BR)-like";
    }
    for (const auto& e : em)
      if (e.label == "sik_telemetry") return "SiK telemetry radio (MAVLink)";
    {
      const Json* fr = nullptr;
      for (const auto& e : em)
        if (e.label.rfind("frsky_", 0) == 0 && e.decoded.is_object() && !e.decoded.empty()) {
          if (!fr) fr = &e.decoded;
          if (e.decoded["crc_ok"].is_boolean() && e.decoded["crc_ok"].get<bool>()) { fr = &e.decoded; break; }
        }
      if (fr) return (*fr)["protocol"].get<std::string>();
    }
    {
      std::vector<const Json*> zb, ok;
      for (const auto& e : em)
        if (e.label.rfind("ieee802154", 0) == 0 && e.decoded.is_object() && !e.decoded.empty()) {
          zb.push_back(&e.decoded);
          if (e.decoded["crc_ok"].get<bool>()) ok.push_back(&e.decoded);
        }
      if (!zb.empty()) {
        if (ok.empty()) return "IEEE 802.15.4-like (O-QPSK 2.4 GHz)";
        for (const Json* d : ok)
          if (d->contains("vendor")) return "IEEE 802.15.4 frames (" + (*d)["vendor"].get<std::string>() + ")";
        return "IEEE 802.15.4 frames (XBee / Zigbee-like)";
      }
    }
    for (const auto& e : em)
      if (e.label == "toy_drone" && e.decoded.is_object() && !e.decoded.empty())
        return "Toy drone RC link: " + e.decoded["protocol"].get<std::string>();
    {
      const Json* xn = nullptr;
      for (const auto& e : em)
        if ((e.label == "bayang" || e.label == "xn297_like") && e.decoded.is_object() && !e.decoded.empty()) {
          if (!xn) xn = &e.decoded;
          if (e.decoded["checksum_ok"].is_boolean() && e.decoded["checksum_ok"].get<bool>()) { xn = &e.decoded; break; }
        }
      if (xn) return (*xn)["protocol"].get<std::string>();
    }
    {
      const Json* af = nullptr;
      for (const auto& e : em)
        if ((e.label.rfind("afhds2a", 0) == 0 || e.label.rfind("afhds1", 0) == 0) && e.decoded.is_object() &&
            !e.decoded.empty()) {
          if (!af) af = &e.decoded;
          if (e.decoded["crc_ok"].is_boolean() && e.decoded["crc_ok"].get<bool>()) { af = &e.decoded; break; }
        }
      if (af) return (*af)["protocol"].get<std::string>();
    }
    for (const auto& e : em)
      if (e.label == "elrs_gfsk") return "ExpressLRS GFSK (LR1121)";
    {
      const Json* rl = nullptr;
      for (const auto& e : em)
        if (e.label == "radiolink" && e.decoded.is_object() && !e.decoded.empty()) {
          if (!rl) rl = &e.decoded;
          if (e.decoded.value("packet_name", std::string()) == "data") { rl = &e.decoded; break; }
        }
      if (rl) return (*rl)["protocol"].get<std::string>();
    }
    for (const auto& e : em)
      if (e.label == "hubsan" && e.decoded.is_object() && !e.decoded.empty()) return "Hubsan (A7105)";
    for (const auto& e : em)   // HoTT / Hitec / S-FHSS (PLAN G12)
      if ((e.label == "hott" || e.label == "hitec" || e.label == "sfhss") && e.decoded.is_object() && !e.decoded.empty())
        return replace_all(e.decoded["protocol"].get<std::string>(), " receiver telemetry", "");
    {   // OGN tracker / ADS-L / FANET (PLAN G10)
      bool og = false, ogu = false, ad = false, adu = false, fa = false, fau = false;
      for (const auto& e : em) {
        if (!e.decoded.is_object() || e.decoded.empty()) continue;
        if (e.label == "ogn_like") og = true, ogu = ogu || e.decoded.value("uav", false);
        if (e.label == "adsl") ad = true, adu = adu || e.decoded.value("uas", false);
        if (e.decoded.value("protocol", std::string()) == "FANET" && e.decoded.value("crc_ok", false)) {
          fa = true;
          const Json& f = e.decoded.contains("fields") && e.decoded["fields"].is_object() ? e.decoded["fields"] : Json::object();
          fau = fau || f.value("uav", false);
        }
      }
      if (og) return std::string("OGN tracker-like (OGNTP)") + (ogu ? ", UAV" : "");
      if (ad) return std::string("ADS-L iConspicuity (SRD-860)") + (adu ? ", UAS" : "");
      if (fa) return std::string("FANET aircraft tracking (LoRa)") + (fau ? ", UAV" : "");
    }
    for (const auto& e : em)
      if (e.label == "flarm_like") return "FLARM-like (868 MHz collision-avoidance beacon)";
    for (const auto& e : em)
      if (e.label == "crossfire_sync") return "TBS Crossfire-like (150 Hz FSK, sync 0x3E6B3E)";
    if (top.find("mlrs") != std::string::npos)
      return replace_all(replace_all(replace_all(top, "_", " "), "mlrs", "mLRS"), "lora", "LoRa");
    if (top.find("lora") != std::string::npos || top.find("elrs") != std::string::npos)
      return replace_all(replace_all(replace_all(top, "_", " "), "elrs", "ExpressLRS"), "lora", "LoRa");
    return replace_all(top, "_", " ");
  }

  std::string modulation(const std::vector<Emission>& em) const override {
    const Json& f = em[0].features;
    const Json m = f.contains("modulation") ? f["modulation"] : Json(nullptr);
    if (m.is_string() && m.get<std::string>() == "lora") return "LoRa CSS";
    if (m.is_string() && !m.get<std::string>().empty()) return m.get<std::string>();
    return "GFSK";
  }

  // the summary keeps scalar features plus sync_syms; python also keeps elrs_channel (a
  // tuple, not a list), which the generic adapter cannot tell from a list
  DecodeResult summarise(const std::vector<Emission>& em, const Signal* sig) const override {
    DecodeResult r = EmissionReceiver::summarise(em, sig);
    {  // ExpressLRS 2.4 GHz (SX1280): IQ inverted <=> UID[5] bit 0 = 1 or binding (python summarise; research)
      int n = 0, ninv = 0;
      for (const auto& e : em) {
        const Json& f = e.features;
        if (e.label.find("elrs") == std::string::npos || e.center_mhz < 2400.0 || !f.contains("modulation") ||
            f["modulation"] != "lora" || !f.contains("iq_inverted") || !f["iq_inverted"].is_boolean())
          continue;
        ++n;
        ninv += f["iq_inverted"].get<bool>();
      }
      if (n && (ninv == 0 || ninv == n)) {
        const bool v = ninv == n;
        r.fields["elrs_identity"] = {{"uid5_bit0", v ? 1 : 0}, {"iq_inverted", v}, {"packets", n},
                                     {"note", "binding phrase UID[5] bit 0 from the LoRa IQ polarity (1 also while binding)"}};
      }
    }
    {  // ExpressLRS SYNC packets verified by their own UID bytes (python summarise; research, PLAN R2c)
      const Json* last = nullptr;
      int n = 0;
      for (const auto& e : em) {
        const Json& d = e.decoded;
        if (!d.is_object() || d.value("protocol", std::string()) != "ExpressLRS" ||
            d.value("packet_type", std::string()) != "SYNC" || !d.contains("crc_ok") || d["crc_ok"] != true)
          continue;
        last = &d;
        ++n;
      }
      if (last) {
        Json f = (*last)["sync"];
        f["packets"] = n;
        f["ota_version"] = (*last)["ota_version"];
        f["binding"] = (*last)["binding"];
        if (last->contains("hop_recovery")) f["hop_recovery"] = (*last)["hop_recovery"];   // research, PLAN R2d
        r.fields["elrs_sync"] = f;
      }
    }
    {  // OGN tracker / ADS-L / FANET (research, PLAN G10): addresses, aircraft type / category, last position
      auto pos3 = [](const Json& d, bool all3) {
        Json p = Json::object();
        for (const char* k : {"lat", "lon", "alt_m"})
          if (all3 || (d.contains(k) && !d[k].is_null())) p[k] = d.contains(k) ? d[k] : Json(nullptr);
        return p;
      };
      std::vector<const Json*> og, ad, fa;
      for (const auto& e : em) {
        if (!e.decoded.is_object() || e.decoded.empty()) continue;
        if (e.label == "ogn_like") og.push_back(&e.decoded);
        if (e.label == "adsl") ad.push_back(&e.decoded);
        if (e.decoded.value("protocol", std::string()) == "FANET" && e.decoded.value("crc_ok", false) &&
            e.decoded.contains("fields") && e.decoded["fields"].is_object() && !e.decoded["fields"].empty())
          fa.push_back(&e.decoded["fields"]);
      }
      if (!og.empty()) {
        std::set<std::string> ad_, at;
        bool uav = false, enc = false;
        const Json* last = nullptr;
        for (const Json* d : og) {
          ad_.insert((*d)["address_type"].get<std::string>() + ":" + (*d)["address"].get<std::string>());
          if (d->contains("aircraft_type")) at.insert((*d)["aircraft_type"].get<std::string>());
          uav = uav || d->value("uav", false);
          enc = enc || (*d)["encrypted"].get<bool>();
          if (d->contains("lat") && !(*d)["lat"].is_null()) last = d;
        }
        r.fields["ogn"] = {{"packets", og.size()}, {"addresses", Json(ad_)}, {"aircraft_types", Json(at)}, {"uav", uav},
                           {"encrypted", enc}, {"position", last ? pos3(*last, true) : Json(nullptr)},
                           {"fec", "LDPC not checked"}};
      }
      if (!ad.empty()) {
        std::set<std::string> addr, pt, cats;
        bool uas = false;
        const Json* last = nullptr;
        for (const Json* d : ad) {
          addr.insert((*d)["address_table"].get<std::string>() + ":" + (*d)["address"].get<std::string>());
          pt.insert((*d)["payload_type"].get<std::string>());
          if (d->contains("aircraft_category") && (*d)["aircraft_category"].is_string() &&
              !(*d)["aircraft_category"].get<std::string>().empty())
            cats.insert((*d)["aircraft_category"].get<std::string>());
          uas = uas || d->value("uas", false);
          if (d->contains("lat") && !(*d)["lat"].is_null()) last = d;
        }
        r.fields["adsl"] = {{"packets", ad.size()}, {"addresses", Json(addr)}, {"payload_types", Json(pt)},
                            {"aircraft_categories", Json(cats)}, {"uas", uas},
                            {"position", last ? pos3(*last, false) : Json(nullptr)}};
      }
      if (!fa.empty()) {
        std::set<std::string> ids, mf, mt, at, names;
        bool uav = false;
        const Json* last = nullptr;
        for (const Json* d : fa) {
          ids.insert((*d)["device_id"].get<std::string>());
          mf.insert((*d)["manufacturer"].get<std::string>());
          mt.insert((*d)["message_type"].get<std::string>());
          if (d->contains("aircraft_type") && (*d)["aircraft_type"].is_string()) at.insert((*d)["aircraft_type"].get<std::string>());
          uav = uav || d->value("uav", false);
          if (d->contains("name") && (*d)["name"].is_string() && !(*d)["name"].get<std::string>().empty())
            names.insert((*d)["name"].get<std::string>());
          if (d->contains("lat") && !(*d)["lat"].is_null()) last = d;
        }
        r.fields["fanet"] = {{"packets", fa.size()}, {"device_ids", Json(ids)}, {"manufacturers", Json(mf)},
                             {"message_types", Json(mt)}, {"aircraft_types", Json(at)}, {"uav", uav},
                             {"position", last ? pos3(*last, false) : Json(nullptr)}, {"names", Json(names)}};
      }
    }
    {  // research catalog of RC / telemetry protocols (python protocols.candidates)
      std::vector<const Emission*> fsk;
      for (const auto& e : em)
        if ((e.family == "rc_fhss" || e.family == "narrowband" || e.family == "unknown") && e.features.contains("symbol_rate") &&
            e.features["symbol_rate"].is_number() && e.features["symbol_rate"].get<double>() != 0 &&
            e.features.value("modulation", Json(nullptr)) == "gfsk")
          fsk.push_back(&e);
      if (!fsk.empty()) {
        dvec R, dv, fc, st;
        for (const auto* e : fsk) {
          R.push_back(e->features["symbol_rate"].get<double>());
          if (e->features.contains("dev_hz") && e->features["dev_hz"].is_number() && e->features["dev_hz"].get<double>() != 0)
            dv.push_back(e->features["dev_hz"].get<double>());
          fc.push_back(e->center_mhz);
          st.push_back(e->start_us);
        }
        std::sort(st.begin(), st.end());
        dvec gaps;
        for (size_t i = 1; i < st.size(); ++i)
          if (st[i] - st[i - 1] > 300.0) gaps.push_back(st[i] - st[i - 1]);
        std::optional<double> iv;
        if (st.size() >= 3 && !gaps.empty()) iv = *std::min_element(gaps.begin(), gaps.end()) / 1e3;
        const Json pc = protocols::candidates(np::median(R), dv.empty() ? std::nullopt : std::optional<double>(np::median(dv)),
                                              np::median(fc), iv);
        if (!pc.empty()) r.fields["protocol_candidates"] = pc;
      }
    }
    Json& bursts = r.fields["bursts"];
    for (size_t i = 0; i < bursts.size() && i < em.size(); ++i) {
      Json f = Json::object();
      const Json& src = em[i].features;
      for (auto it = src.begin(); it != src.end(); ++it)
        if (!(it.value().is_array() || it.value().is_object()) || it.key() == "sync_syms" ||
            it.key() == "elrs_channel")
          f[it.key()] = it.value();
      bursts[i]["features"] = f;
    }
    // Remote ID (F3411 over BLE; research capability): later messages complete / update earlier ones
    std::vector<Json> rids, ads;
    for (const auto& e : em) {
      if (!e.decoded.is_object() || !e.decoded.contains("ad")) continue;
      ads.push_back(e.decoded["ad"]);
      if (truthy(e.decoded, "crc_ok") && e.decoded["ad"].contains("remote_id")) rids.push_back(e.decoded["ad"]["remote_id"]);
    }
    if (!rids.empty()) {
      Json merged = Json::object();
      for (const auto& rr : rids)
        for (auto it = rr.begin(); it != rr.end(); ++it) {
          if (it.key() == "messages") {
            if (!merged.contains("messages")) merged["messages"] = Json::array();
            for (const auto& m : it.value())
              if (std::find(merged["messages"].begin(), merged["messages"].end(), m) == merged["messages"].end())
                merged["messages"].push_back(m);
          } else if (!it.value().is_null()) {
            merged[it.key()] = it.value();
          }
        }
      r.fields["remote_id"] = merged;
      r.fields["remote_id_adverts"] = rids.size();
      if (merged.contains("uas_id") && truthy(merged, "uas_id")) r.fields["uas_id"] = merged["uas_id"];
      r.fields["model"] = drone_models::remote_id_model(merged);   // research, PLAN M2
      std::string msgs;
      if (merged.contains("messages"))
        for (const auto& m : merged["messages"]) msgs += (msgs.empty() ? "" : ", ") + m.get<std::string>();
      r.evidence = r.evidence.value_or("") + "; Remote ID (ASTM F3411) in " + std::to_string(rids.size()) +
                   " CRC-valid advert(s): " + msgs;
    }
    {  // SiK telemetry (research): radio facts + the vehicle state from the CRC-valid MAVLink frames
      std::vector<const Json*> sk;
      for (const auto& e : em)
        if (e.label == "sik_telemetry" && truthy(e.decoded, "crc_ok")) sk.push_back(&e.decoded);
      if (!sk.empty()) {
        Json frames = Json::array();
        int nok = 0;
        for (const auto* d : sk) {
          for (const auto& f : (*d)["mavlink_frames"]) frames.push_back(f);
          nok += (*d)["mavlink_crc_ok"].get<int>();
        }
        const Json& d0 = *sk[0];
        r.fields["sik"] = {{"netid", d0["netid"]}, {"air_rate_bps", d0["air_rate_bps"]}, {"whitening", d0["whitening"]},
                           {"hw_crc_variant", d0["hw_crc_variant"]}, {"packets", sk.size()}, {"mavlink_frames_crc_ok", nok}};
        const Json mv = mavlink::summary(frames);
        if (!mv.empty()) {
          r.fields["mavlink"] = mv;
          std::string msgs;
          if (mv.contains("messages"))
            for (const auto& m : mv["messages"]) msgs += (msgs.empty() ? "" : ", ") + m.get<std::string>();
          r.evidence = r.evidence.value_or("") + "; SiK NETID " + std::to_string(d0["netid"].get<int>()) + ", MAVLink " + msgs;
        }
      }
    }
    if (!ads.empty()) {
      Json keep = Json::object();
      for (const auto& a : ads)
        for (const char* k : {"local_name", "company_id", "company", "service_uuids16", "tx_power_dbm", "appearance"})
          if (a.contains(k) && !a[k].is_null()) keep[k] = a[k];
      if (!keep.empty()) r.fields["ble_ad"] = keep;
      // trackers / notable devices (research, signal gaps S7): distinct recognitions of CRC-valid adverts
      Json devs = Json::array();
      std::set<std::pair<std::string, std::string>> seen;
      for (const auto& e : em) {
        if (!truthy(e.decoded, "crc_ok") || !e.decoded.contains("ad") || !e.decoded["ad"].contains("device")) continue;
        const Json& dv = e.decoded["ad"]["device"];
        const std::string addr = e.decoded.contains("adv_addr") && e.decoded["adv_addr"].is_string()
                                     ? e.decoded["adv_addr"].get<std::string>() : std::string("\x01none");
        if (!seen.insert({dv["label"].get<std::string>(), addr}).second) continue;
        Json d = dv;
        d["adv_addr"] = e.decoded.contains("adv_addr") ? e.decoded["adv_addr"] : Json(nullptr);
        devs.push_back(d);
      }
      if (!devs.empty()) {
        Json top = Json::array();
        for (size_t i = 0; i < devs.size() && i < 32; ++i) top.push_back(devs[i]);
        r.fields["ble_devices"] = top;
        std::set<std::string> trk;
        for (const auto& d : devs)
          if (d["tracker"].get<bool>()) trk.insert(d["label"].get<std::string>());
        if (!trk.empty()) {
          std::string j;
          for (const auto& t : trk) j += (j.empty() ? "" : ", ") + t;
          r.evidence = r.evidence.value_or("") + "; item tracker(s) in the adverts: " + j;
        }
      }
    }
    {
      // Bluetooth BR/EDR (research, signal gaps S7): LAPs of the BCH-verified access codes; UAP / NAP not on air
      std::vector<const Json*> br;
      for (const auto& e : em)
        if (e.label == "bluetooth_bredr" && e.decoded.is_object() && !e.decoded.empty()) br.push_back(&e.decoded);
      if (!br.empty()) {
        std::set<std::string> laps;
        bool inq = false;
        int64_t mx = 0;
        for (const Json* d : br) {
          laps.insert((*d)["lap"].get<std::string>());
          inq = inq || (*d)["inquiry"].get<bool>();
          mx = std::max(mx, (*d)["access_code_errors"].get<int64_t>());
        }
        Json f = Json::object();
        f["laps"] = Json(std::vector<std::string>(laps.begin(), laps.end()));
        f["packets"] = static_cast<int64_t>(br.size());
        f["inquiry"] = inq;
        f["max_access_code_errors"] = mx;
        r.fields["bredr"] = f;
        std::string j;
        for (const auto& l : laps) j += (j.empty() ? "" : ", ") + l;
        r.evidence = r.evidence.value_or("") + "; " + std::to_string(br.size()) + " BR/EDR access code(s) within " +
                     std::to_string(mx) + " bit(s) of the (64,30) BCH codeword of LAP " + j;
      }
    }
    {
      // TBS Crossfire 150 Hz (research, PLAN R9): sync-word packets, grid channels and links seen
      std::vector<const Json*> cf;
      for (const auto& e : em)
        if (e.label == "crossfire_sync" && e.decoded.is_object() && !e.decoded.empty()) cf.push_back(&e.decoded);
      if (!cf.empty()) {
        std::set<std::pair<std::string, int64_t>> chs;
        Json band = nullptr;
        for (const Json* d : cf)
          if (d->contains("grid")) {
            const Json& g = (*d)["grid"];
            if (band.is_null()) band = g["band"];
            chs.emplace(g["link"].get<std::string>(), g["channel"].get<int64_t>());
          }
        Json chj = Json::array();
        for (const auto& [l, c] : chs) chj.push_back(Json::array({l, c}));
        r.fields["crossfire"] = {{"packets", cf.size()}, {"sync", (*cf[0])["sync"]}, {"band", band},
                                 {"channels", chj}, {"note", "packets read but not CRC-verified (per-hop CRC init unknown)"}};
      }
    }
    {
      // FrSky D8 / D16 (research, PLAN R6): variant from the CRC structure + rate, TX ID from bytes 1-2
      std::vector<const Json*> fr, ok;
      for (const auto& e : em)
        if (e.label.rfind("frsky_", 0) == 0 && e.decoded.is_object() && !e.decoded.empty()) {
          fr.push_back(&e.decoded);
          if (e.decoded["crc_ok"].is_boolean() && e.decoded["crc_ok"].get<bool>()) ok.push_back(&e.decoded);
        }
      if (!fr.empty()) {
        std::set<std::string> vars, ids;
        for (const Json* d : ok.empty() ? fr : ok) vars.insert((*d)["variant"].get<std::string>());
        for (const Json* d : ok)
          if ((*d)["tx_id"].is_string()) ids.insert((*d)["tx_id"].get<std::string>());
        r.fields["frsky"] = {{"packets", fr.size()}, {"crc_ok", ok.size()}, {"variants", Json(vars)},
                             {"tx_ids", Json(ids)}, {"sync", (*fr[0])["sync"]}};
      }
    }
    {
      // FlySky AFHDS 2A (research, PLAN R7): ID-code packets; CRC-verified ones give the IDs and packet types
      std::vector<const Json*> af, ok;
      for (const auto& e : em)
        if (e.label.rfind("afhds2a", 0) == 0 && e.decoded.is_object() && !e.decoded.empty()) {
          af.push_back(&e.decoded);
          if (e.decoded["crc_ok"].is_boolean() && e.decoded["crc_ok"].get<bool>()) ok.push_back(&e.decoded);
        }
      if (!af.empty()) {
        std::set<std::string> tx, rx, types;
        bool bind = false;
        for (const Json* d : ok) {
          tx.insert((*d)["tx_id"].get<std::string>());
          rx.insert((*d)["rx_id"].get<std::string>());
          types.insert((*d)["packet_type"].get<std::string>());
          bind = bind || (*d)["binding"].get<bool>();
        }
        r.fields["afhds2a"] = {{"packets", af.size()}, {"crc_ok", ok.size()}, {"tx_ids", Json(tx)}, {"rx_ids", Json(rx)},
                               {"packet_types", Json(types)}, {"binding", bind},
                               {"convention", ok.empty() ? Json(nullptr) : (*ok[0])["convention"]}};
      }
    }
    auto collect = [&](auto pred) {
      std::vector<const Json*> v;
      for (const auto& e : em)
        if (pred(e.label) && e.decoded.is_object() && !e.decoded.empty()) v.push_back(&e.decoded);
      return v;
    };
    auto truthy = [](const Json& d, const char* k) { return d.contains(k) && d[k].is_boolean() && d[k].get<bool>(); };
    {
      // FlySky AFHDS 1 (research, PLAN G5): TX ID, hop table check, CRC convention when one passes
      const auto a1 = collect([](const std::string& l) { return l.rfind("afhds1", 0) == 0; });
      if (!a1.empty()) {
        std::vector<const Json*> ok;
        for (const Json* d : a1)
          if (truthy(*d, "crc_ok")) ok.push_back(d);
        std::set<std::string> tx, hc;
        for (const Json* d : ok.empty() ? a1 : ok) tx.insert((*d)["tx_id"].get<std::string>());
        bool bind = false;
        for (const Json* d : a1) {
          if (d->contains("hop_check")) hc.insert((*d)["hop_check"].get<std::string>());
          bind = bind || (*d)["packet_name"] == "bind";
        }
        r.fields["afhds1"] = {{"packets", a1.size()}, {"crc_ok", ok.size()}, {"tx_ids", Json(tx)},
                              {"hop_checks", Json(hc)}, {"hop_channels", (*a1[0])["hop_channels"]}, {"binding", bind},
                              {"convention", ok.empty() ? Json(nullptr) : (*ok[0])["convention"]}};
      }
    }
    {
      // ExpressLRS LR1121 GFSK (research, PLAN G3): verified SYNC packets -> UID bytes, OTA version, whitening
      const auto eg = collect([](const std::string& l) { return l == "elrs_gfsk"; });
      if (!eg.empty()) {
        std::set<std::string> wv;
        bool bind = false;
        for (const Json* d : eg) {
          wv.insert((*d)["whitening"].get<std::string>());
          bind = bind || (*d)["binding"].get<bool>();
        }
        const Json& d0 = *eg[0];
        r.fields["elrs_gfsk"] = {{"sync_packets", eg.size()}, {"uid4", d0["uid4"]}, {"uid5", d0["uid5"]},
                                 {"sync_word", d0["sync_word"]}, {"ota", d0["ota"]}, {"ota_version", d0["ota_version"]},
                                 {"whitening", Json(wv)}, {"binding", bind}};
      }
    }
    {
      // RadioLink (research, PLAN G4): sub-protocol, TX IDs, hop-table check, telemetry
      const auto rl = collect([](const std::string& l) { return l == "radiolink"; });
      if (!rl.empty()) {
        std::set<std::string> subs, ids, kinds, hc;
        Json hop = nullptr;
        bool tlm = false;
        for (const Json* d : rl) {
          subs.insert((*d)["sub_protocol"].get<std::string>());
          ids.insert((*d)["tx_id"].get<std::string>());
          kinds.insert((*d)["packet_name"].get<std::string>());
          if (d->contains("hop_check")) hc.insert((*d)["hop_check"].get<std::string>());
          if (hop.is_null() && d->contains("hop_channels")) hop = (*d)["hop_channels"];
          tlm = tlm || truthy(*d, "telemetry_request");
        }
        r.fields["radiolink"] = {{"packets", rl.size()}, {"sub_protocols", Json(subs)}, {"tx_ids", Json(ids)},
                                 {"packet_kinds", Json(kinds)}, {"hop_checks", Json(hc)}, {"hop_channels", hop},
                                 {"telemetry_requested", tlm}};
      }
    }
    {
      // HoTT / Hitec / S-FHSS (research, PLAN G12): TX IDs, packet kinds, hop table / index, model hint
      const auto cc = collect([](const std::string& l) { return l == "hott" || l == "hitec" || l == "sfhss"; });
      if (!cc.empty()) {
        std::set<std::string> pr, ids, rx, kinds, hc, mh;
        Json table = nullptr;
        for (const Json* d : cc) {
          pr.insert((*d)["protocol"].get<std::string>());
          ids.insert((*d)["tx_id"].get<std::string>());
          if (d->contains("rx_id") && (*d)["rx_id"].is_string() && !(*d)["rx_id"].get<std::string>().empty())
            rx.insert((*d)["rx_id"].get<std::string>());
          kinds.insert((*d)["packet_name"].get<std::string>());
          if (table.is_null() && d->contains("hop_table") && !(*d)["hop_table"].is_null()) table = (*d)["hop_table"];
          if (d->contains("hop_check") && (*d)["hop_check"].is_string()) hc.insert((*d)["hop_check"].get<std::string>());
          if (d->contains("model_hint") && (*d)["model_hint"].is_string()) mh.insert((*d)["model_hint"].get<std::string>());
        }
        r.fields["cc2500_rc"] = {{"packets", cc.size()}, {"protocols", Json(pr)}, {"tx_ids", Json(ids)},
                                 {"rx_ids", Json(rx)}, {"packet_kinds", Json(kinds)}, {"hop_table", table},
                                 {"hop_checks", Json(hc)}, {"model_hints", Json(mh)}};
      }
    }
    {
      // Hubsan (research, PLAN G5): session ID codes, packet kinds, channel checks
      const auto hb = collect([](const std::string& l) { return l == "hubsan"; });
      if (!hb.empty()) {
        std::set<std::string> ids, kinds, cc;
        for (const Json* d : hb) {
          ids.insert((*d)["id_code"].get<std::string>());
          kinds.insert((*d)["packet_name"].get<std::string>());
          if (d->contains("channel_check")) cc.insert((*d)["channel_check"].get<std::string>());
        }
        r.fields["hubsan"] = {{"packets", hb.size()}, {"id_codes", Json(ids)}, {"packet_kinds", Json(kinds)},
                              {"channel_checks", Json(cc)}};
      }
    }
    {
      // toy-drone protocols on nRF24 / XN297 (research, PLAN G6): protocol, addresses, TX IDs
      const auto ty = collect([](const std::string& l) { return l == "toy_drone"; });
      if (!ty.empty()) {
        std::set<std::string> pr, ad, ids;
        int64_t ok = 0;
        for (const Json* d : ty) {
          pr.insert((*d)["protocol"].get<std::string>());
          ad.insert((*d)["address"].get<std::string>());
          ids.insert((*d)["tx_id"].get<std::string>());
          ok += truthy(*d, "checksum_ok");
        }
        r.fields["toy_drone"] = {{"packets", ty.size()}, {"protocols", Json(pr)}, {"addresses", Json(ad)},
                                 {"tx_ids", Json(ids)}, {"checksum_ok", ok}};
      }
    }
    {
      // XN297 / Bayang (research, PLAN R8): addresses, Bayang TX IDs, packet kinds, CRC residuals
      std::vector<const Json*> xn, ok;
      for (const auto& e : em)
        if ((e.label == "bayang" || e.label == "xn297_like") && e.decoded.is_object() && !e.decoded.empty()) {
          xn.push_back(&e.decoded);
          if (e.decoded["checksum_ok"].is_boolean() && e.decoded["checksum_ok"].get<bool>()) ok.push_back(&e.decoded);
        }
      if (!xn.empty()) {
        std::set<std::string> addrs, ids, kinds, res;
        for (const Json* d : xn) addrs.insert((*d)["address"].get<std::string>());
        for (const Json* d : ok) {
          if (d->contains("tx_id")) ids.insert((*d)["tx_id"].get<std::string>());
          kinds.insert((*d)["packet_name"].get<std::string>());
          res.insert((*d)["crc_residual"].get<std::string>());
        }
        r.fields["xn297"] = {{"packets", xn.size()}, {"bayang_checksum_ok", ok.size()}, {"addresses", Json(addrs)},
                             {"tx_ids", Json(ids)}, {"packet_kinds", Json(kinds)}, {"crc_residuals", Json(res)}};
      }
    }
    {
      // IEEE 802.15.4 (research, PLAN T4): CRC-verified frames, channel, PANs, addresses, vendors
      std::vector<const Json*> zb, ok;
      for (const auto& e : em)
        if (e.label.rfind("ieee802154", 0) == 0 && e.decoded.is_object() && !e.decoded.empty()) {
          zb.push_back(&e.decoded);
          if (e.decoded["crc_ok"].get<bool>()) ok.push_back(&e.decoded);
        }
      if (!zb.empty()) {
        std::set<int64_t> chs;
        std::set<std::string> types, pans, addrs, vendors;
        for (const Json* d : zb) chs.insert((*d)["channel"].get<int64_t>());
        for (const Json* d : ok) {
          types.insert((*d)["frame_type"].is_string() ? (*d)["frame_type"].get<std::string>() : (*d)["frame_type"].dump());
          for (const char* k : {"dst_pan", "src_pan"})
            if (d->contains(k) && (*d)[k].is_string()) pans.insert((*d)[k].get<std::string>());
          for (const char* k : {"dst_addr", "src_addr"})
            if (d->contains(k) && (*d)[k].is_string()) addrs.insert((*d)[k].get<std::string>());
          if (d->contains("vendor")) vendors.insert((*d)["vendor"].get<std::string>());
        }
        r.fields["ieee802154"] = {{"frames", zb.size()}, {"crc_ok", ok.size()}, {"channels", Json(chs)},
                                  {"frame_types", Json(types)}, {"pan_ids", Json(pans)}, {"addresses", Json(addrs)},
                                  {"vendors", Json(vendors)}};
      }
    }
    {
      // FLARM (research, PLAN T5): clear header only; the payload is obfuscated and not decoded
      std::vector<const Json*> fl;
      for (const auto& e : em)
        if (e.label == "flarm_like" && e.decoded.is_object() && !e.decoded.empty()) fl.push_back(&e.decoded);
      if (!fl.empty()) {
        std::set<std::string> ids, types;
        for (const Json* d : fl) {
          if ((*d)["radio_id"].is_string()) ids.insert((*d)["radio_id"].get<std::string>());
          if ((*d)["address_type"].is_string()) types.insert((*d)["address_type"].get<std::string>());
        }
        r.fields["flarm"] = {{"beacons", fl.size()}, {"radio_ids", Json(ids)}, {"address_types", Json(types)}};
      }
    }
    {
      // mLRS (research, PLAN R3 / T6): air modes; CRC-verified frames give the per-binding FrameSyncWord
      int64_t n = 0, nok = 0;
      std::set<std::string> sws, modes;
      for (const auto& e : em) {
        if (e.label.rfind("mlrs", 0) != 0) continue;
        ++n;
        if (e.decoded.is_object() && e.decoded.value("crc_ok", false)) ++nok, sws.insert(e.decoded["frame_sync_word"].get<std::string>());
        if (e.features.contains("mlrs_modes"))
          for (const auto& m : e.features["mlrs_modes"]) modes.insert(m.get<std::string>());
      }
      if (n) r.fields["mlrs"] = {{"packets", n}, {"crc_ok", nok}, {"frame_sync_words", Json(sws)}, {"modes", Json(modes)}};
    }
    {
      // blind FSK frame structure (research): sync-word candidates on their first 24 bits (python Counter.most_common:
      // count descending, ties in first-seen order), preamble lengths, known sync words
      std::vector<std::pair<std::string, int64_t>> cnt;
      std::set<int64_t> pre;
      std::set<std::string> known;
      int64_t nfr = 0;
      for (const auto& e : em) {
        if (!e.features.contains("frame") || !e.features["frame"].is_object()) continue;
        const Json& f = e.features["frame"];
        ++nfr;
        const std::string k = f["sync32_hex"].get<std::string>().substr(0, 6);
        auto it = std::find_if(cnt.begin(), cnt.end(), [&](const auto& p) { return p.first == k; });
        if (it == cnt.end()) cnt.emplace_back(k, 1); else ++it->second;
        pre.insert(f["preamble_bits"].get<int64_t>());
        if (f.contains("known_sync")) known.insert(f["known_sync"]["name"].get<std::string>());
      }
      if (nfr) {
        std::stable_sort(cnt.begin(), cnt.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
        Json cj = Json::array();
        for (size_t i = 0; i < cnt.size() && i < 4; ++i) cj.push_back(Json::array({cnt[i].first, cnt[i].second}));
        r.fields["fsk_frame"] = {{"packets", nfr}, {"sync24_candidates", cj}, {"preamble_bits", Json(pre)},
                                 {"known_sync", Json(known)}};
      }
    }
    if (!r.fields.contains("model")) {   // model level (research, PLAN M5)
      Json mb = models::model(em);
      if (!mb.is_null()) r.fields["model"] = std::move(mb);
    }
    {
      // low-SNR coherent GFSK path (research, PLAN L1): adverts, CRC Chase corrections, measured modulation index
      int64_t n = 0, fixes = 0;
      dvec hs;
      for (const auto& e : em) {
        const Json& d = e.decoded;
        if (!d.is_object() || !d.contains("crc_ok") || !d["crc_ok"].is_boolean() || !d["crc_ok"].get<bool>()) continue;
        if (!d.contains("demod") || d["demod"] != "coherent") continue;
        ++n;
        fixes += d["bit_corrections"].get<int64_t>();
        hs.push_back(d["mod_index"].get<double>());
      }
      if (n) {
        Json c = Json::object();
        c["adverts"] = n;
        c["bit_corrections"] = fixes;
        c["mod_index"] = py_round(np::median(hs), 2);
        r.fields["ble_coherent"] = c;
      }
    }
    return r;
  }

  // A CRC-verified BLE advert passes only when de-whitened with its true channel index, so
  // (labelled frequency - true advertising-channel frequency) is the capture's label error.
  std::optional<Json> frequency_label_check(const std::vector<Signal>& signals) const override {
    std::vector<Emission> ems;
    for (const auto& s : signals)
      if (s.features.contains("decoder") && s.features["decoder"] == name())
        ems.insert(ems.end(), s.emissions.begin(), s.emissions.end());
    return check_frequency_label(ems);
  }
};

}  // namespace

void register_plugins(si::PluginList& out) { out.push_back(std::make_unique<RcLinkDecoder>()); }

}  // namespace si::rx::rc_link
