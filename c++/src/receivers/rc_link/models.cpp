// python receivers/rc_link/models.py (see there for the sources and the rules).
#include "receivers/rc_link/models.hpp"

#include <algorithm>
#include <set>
#include <string>

#include "core/drone_models.hpp"
#include "core/pyfmt.hpp"

namespace si::rx::rc_link::models {
namespace {
const std::set<int64_t> kSikRfd = {12, 56, 100, 125, 200, 224, 500, 750};

bool has(const Emission& e) { return e.decoded.is_object() && !e.decoded.empty(); }
bool starts(const std::string& s, const std::string& p) { return s.rfind(p, 0) == 0; }
bool truthy_bool(const Json& d, const char* k) { return d.contains(k) && d[k].is_boolean() && d[k].get<bool>(); }
}  // namespace

Json model(const std::vector<Emission>& em) {
  using drone_models::model_block;
  std::vector<const Json*> fr, af, cf, sk;
  for (const auto& e : em) {
    if (!has(e)) continue;
    if (starts(e.label, "frsky_")) fr.push_back(&e.decoded);
    if (starts(e.label, "afhds2a")) af.push_back(&e.decoded);
    if (e.label == "crossfire_sync") cf.push_back(&e.decoded);
    if (e.label == "sik_telemetry") sk.push_back(&e.decoded);
  }
  if (!fr.empty()) {
    std::set<std::string> var;
    for (const Json* d : fr)
      if (truthy_bool(*d, "crc_ok")) var.insert((*d)["variant"].get<std::string>());
    if (var.empty())
      for (const Json* d : fr) var.insert((*d)["variant"].get<std::string>());
    std::vector<std::string> ev;
    std::string fam;
    for (const auto& v : var) {
      const char* note = v == "D16 v2" ? ": ACCST v2 firmware generation" : v == "D16 LBT v1" ? ": EU-LBT region build"
                         : v == "D16 v1" ? ": ACCST v1 FCC" : "";
      ev.push_back("FrSky " + v + note);
      fam += (fam.empty() ? "" : " / ") + v;
    }
    return model_block("FrSky " + fam, {"FrSky radio / internal module", "multiprotocol module"}, ev, "decoded",
                       "no published on-air fingerprint separates FrSky radios and multiprotocol modules");
  }
  {   // ADS-L (PLAN G10)
    std::set<std::string> cats, tabs;
    bool any = false;
    for (const auto& e : em)
      if (e.label == "adsl" && has(e)) {
        any = true;
        if (e.decoded.contains("aircraft_category") && e.decoded["aircraft_category"].is_string() &&
            !e.decoded["aircraft_category"].get<std::string>().empty())
          cats.insert(e.decoded["aircraft_category"].get<std::string>());
        tabs.insert(e.decoded["address_table"].get<std::string>());
      }
    if (any) {
      auto join = [](const std::set<std::string>& v) {
        std::string o;
        for (const auto& x : v) o += (o.empty() ? "" : ", ") + x;
        return o;
      };
      const std::string c = join(cats);
      return model_block("ADS-L device", {"SkyEcho 2", "FLARM (ADS-L firmware)", "OGN tracker", "SoftRF"},
                         {"ADS-L Traffic, category " + (c.empty() ? std::string("n/a") : c) + ", address table " + join(tabs)},
                         "inferred", "ADS-L packets do not name the device");
    }
  }
  {   // FANET (PLAN G10)
    std::set<std::string> mf, ids;
    for (const auto& e : em)
      if (has(e) && e.decoded.value("protocol", std::string()) == "FANET" && truthy_bool(e.decoded, "crc_ok") &&
          e.decoded.contains("fields") && e.decoded["fields"].is_object() && !e.decoded["fields"].empty()) {
        mf.insert(e.decoded["fields"]["manufacturer"].get<std::string>());
        ids.insert(e.decoded["fields"]["manufacturer_id"].get<std::string>());
      }
    if (!mf.empty()) {
      std::string v;
      for (const auto& x : ids) v += (v.empty() ? "" : ", ") + x;
      return model_block("FANET device", std::vector<std::string>(mf.begin(), mf.end()), {"FANET manufacturer ID " + v},
                         "decoded", "the manufacturer ID names the maker, not the device model");
    }
  }
  {
    std::set<std::string> vendors;
    bool any = false;
    for (const auto& e : em)
      if (e.label.rfind("ieee802154", 0) == 0 && has(e) && truthy_bool(e.decoded, "crc_ok")) {
        any = true;
        if (e.decoded.contains("vendor")) vendors.insert(e.decoded["vendor"].get<std::string>());
      }
    if (any && !vendors.empty()) {
      std::string v;
      for (const auto& x : vendors) v += (v.empty() ? "" : ", ") + x;
      return model_block("Digi XBee (802.15.4 / Zigbee)", {"XBee 802.15.4", "XBee Zigbee 3.0", "XBee DigiMesh 2.4"},
                         {"64-bit address with the Digi OUI 00 13 A2 (" + v + ")"}, "decoded",
                         "the XBee 2.4 GHz variants share the PHY; the firmware is not visible in the MAC header");
    }
    if (any)
      return model_block("IEEE 802.15.4 device", {}, {"CRC-verified 802.15.4 frames, no vendor OUI in the addresses"},
                         "decoded", "short addresses carry no vendor");
  }
  {   // toy drones (PLAN G6)
    std::set<std::string> names, prots;
    for (const auto& e : em)
      if (e.label == "toy_drone" && has(e)) {
        const std::string p = e.decoded["protocol"].get<std::string>();
        prots.insert(p);
        names.insert(p.substr(0, p.rfind(" (")));
      }
    if (!prots.empty())
      return model_block("Toy drone (nRF24 / XN297)", std::vector<std::string>(names.begin(), names.end()),
                         std::vector<std::string>(prots.begin(), prots.end()), "decoded",
                         "several toy protocols match these packets");
  }
  {
    bool any = false, ok = false;
    for (const auto& e : em)
      if ((e.label == "bayang" || e.label == "xn297_like") && has(e)) {
        any = true;
        ok = ok || truthy_bool(e.decoded, "checksum_ok");
      }
    if (any)
      return model_block(ok ? "Bayang (XN297) toy drone" : "XN297 toy drone", {},
                         {std::string("XN297 prefix 0x710F55") + (ok ? ", Bayang checksum verified" : "")}, "decoded",
                         "many toy-drone brands share the XN297 / Bayang protocol");
  }
  {   // RadioLink (PLAN G4)
    std::set<std::string> subs;
    for (const auto& e : em)
      if (e.label == "radiolink" && has(e)) subs.insert(e.decoded["sub_protocol"].get<std::string>());
    if (!subs.empty()) {
      std::vector<std::string> models, ev;
      for (const auto& sp : subs) {
        if (sp == "air") models.insert(models.end(), {"T8FB", "T8S"});
        else if (sp == "surface") models.insert(models.end(), {"RC4GS V3", "RC6GS V3", "RC8X"});
        else if (sp == "DumboRC") models.push_back("DumboRC transmitter");
        else if (sp == "RC4G") models.push_back("RC4G");
        ev.push_back("RadioLink " + sp + " protocol (sync word / header)");
      }
      return model_block("RadioLink", models, ev, "inferred",
                         "the transmitters of one RadioLink sub-protocol share the air protocol");
    }
  }
  {   // HoTT / Hitec / S-FHSS (PLAN G12)
    const Json* first = nullptr;
    bool minima = false;
    std::set<std::string> hints, mbytes;
    for (const auto& e : em)
      if ((e.label == "hott" || e.label == "hitec" || e.label == "sfhss") && has(e)) {
        if (!first) first = &e.decoded;
        if (e.decoded.value("variant", std::string()) == "Minima") minima = true;
        if (e.decoded.contains("model_hint") && e.decoded["model_hint"].is_string())
          hints.insert(e.decoded["model_hint"].get<std::string>());
        if (e.decoded.contains("model_bytes")) mbytes.insert(e.decoded["model_bytes"].get<std::string>());
      }
    if (first) {
      const std::string p = (*first)["protocol"].get<std::string>();
      if (starts(p, "Graupner"))
        return model_block("Graupner HoTT", {"mz-12 / mz-16 / mz-24 / mz-32", "MC-series HoTT", "HoTT module"},
                           {"CC2500 MSK + FEC, sync D391, CRC-verified HoTT packet"}, "decoded",
                           "every HoTT transmitter uses the same air protocol");
      if (starts(p, "Hitec"))
        return model_block("Hitec", {minima ? "Minima receiver link" : "Aurora 9 / Flash 7 / Flash 8 (Optima)"},
                           {"CC2500 MSK, sync word = TX ID, CRC-verified"}, "decoded",
                           "Optima transmitters share the air protocol");
      std::string mb;
      for (const auto& x : mbytes) mb += (mb.empty() ? "" : ", ") + x;
      std::vector<std::string> cand(hints.begin(), hints.end());
      if (cand.empty()) cand.push_back("S-FHSS transmitter (T6J / T8J / T10J / T14SG / ...)");
      return model_block("Futaba", cand, {"S-FHSS packet, model bytes " + mb}, hints.empty() ? "inferred" : "decoded",
                         "model bytes 3-4 not in the known list");
    }
  }
  {   // Hubsan (PLAN G5)
    std::set<std::string> models;
    int64_t n = 0;
    for (const auto& e : em)
      if (e.label == "hubsan" && has(e)) {
        ++n;
        for (const auto& m : e.decoded["models"]) models.insert(m.get<std::string>());
      }
    if (n)
      return model_block("Hubsan", std::vector<std::string>(models.begin(), models.end()),
                         {"model constant bytes in " + std::to_string(n) + " packet(s)"}, "decoded",
                         "no data packet with a model's constant bytes");
  }
  for (const auto& e : em)   // FlySky AFHDS 1 (PLAN G5)
    if (starts(e.label, "afhds1") && has(e))
      return model_block("FlySky AFHDS 1", {"FS-CT6B", "FS-T6", "Turnigy 9X FlySky module", "multiprotocol module"},
                         {"A7105 ID code, 21-byte AFHDS 1 packet, hop table from the TX ID"}, "decoded",
                         "every AFHDS 1 transmitter uses the same air protocol");
  if (!af.empty()) {
    const bool ok = std::any_of(af.begin(), af.end(), [](const Json* d) { return truthy_bool(*d, "crc_ok"); });
    return model_block("FlySky AFHDS 2A", {"FS-i6", "FS-i6X", "other AFHDS 2A TX / multiprotocol module"},
                       {std::string("A7105 ID code") + (ok ? ", CRC-verified packets" : "")}, "decoded",
                       "every AFHDS 2A transmitter uses the same air protocol");
  }
  {
    bool any = false;
    const Json* ok = nullptr;
    for (const auto& e : em)
      if (e.label == "flarm_like" && has(e)) {
        any = true;
        if (!ok && e.decoded["radio_id"].is_string()) ok = &e.decoded;
      }
    if (any)
      return model_block("FLARM device", {},
                         {ok ? "FLARM radio ID " + (*ok)["radio_id"].get<std::string>() + " (" +
                                   (*ok)["address_type"].get<std::string>() + ")"
                             : std::string("FLARM sync word; header not read reliably")},
                         "decoded", "the device type is in the obfuscated payload (not decoded)");
  }
  if (!cf.empty()) {
    std::vector<std::string> ev = {"Crossfire 150 Hz sync word"};
    for (const Json* d : cf)
      if (d->contains("grid") && (*d)["grid"].is_object()) {
        ev.push_back("band plan " + (*d)["grid"]["band"].get<std::string>() + " (region build, not the model)");
        break;
      }
    return model_block("TBS Crossfire", {"Micro TX V2", "Nano TX", "full-size TX", "Tango 2 / Mambo internal"}, ev,
                       "inferred", "the same protocol on every Crossfire TX; only power differs");
  }
  if (!sk.empty()) {
    const double bps = sk[0]->contains("air_rate_bps") && (*sk[0])["air_rate_bps"].is_number()
                           ? (*sk[0])["air_rate_bps"].get<double>() : 0.0;
    const int64_t kb = py_round_int(bps / 1e3);
    if (kSikRfd.count(kb))
      return model_block("SiK telemetry (RFD 3.x)", {"RFD900x", "RFD900ux", "RFD868x"},
                         {"air rate " + std::to_string(kb) + " kb/s is only in the RFD 3.x firmware set"}, "inferred",
                         "RFD900x / ux / 868x differ only in the band edges");
    return model_block("SiK telemetry", {"100 mW SiK board (3DR / Holybro / mRo / HM-TRP)", "RFD900"},
                       {"air rate " + std::to_string(kb) + " kb/s (standard SiK set)"}, "inferred",
                       "every SiK radio with this air rate looks the same on air");
  }
  {
    std::set<std::string> modes;
    const Json* ok = nullptr;
    bool any = false;
    for (const auto& e : em)
      if (e.label.rfind("mlrs", 0) == 0) {
        any = true;
        if (e.features.contains("mlrs_modes"))
          for (const auto& m : e.features["mlrs_modes"]) modes.insert(m.get<std::string>());
        if (!ok && e.decoded.is_object() && e.decoded.value("crc_ok", false)) ok = &e.decoded;
      }
    if (any) {
      auto chip = [](const std::string& m) -> std::string {
        if (m == "mlrs_19hz_7x") return "SX1276 / LR1121 (SX127x-compatible SF6)";
        if (m == "mlrs_31hz_900" || m == "mlrs_19hz_900") return "SX126x / LLCC68 / LR1121";
        return "SX128x / LR1121";
      };
      std::vector<std::string> cands;
      std::string ml;
      for (const auto& m : modes) {
        const std::string c = chip(m) + "-based mLRS module";
        if (std::find(cands.begin(), cands.end(), c) == cands.end()) cands.push_back(c);
        ml += (ml.empty() ? "" : ", ") + m;
      }
      std::vector<std::string> ev = {"mLRS air mode " + ml};
      if (ok) ev.push_back("FrameSyncWord " + (*ok)["frame_sync_word"].get<std::string>() + " (per binding)");
      return model_block("mLRS", cands, ev, ok ? "decoded" : "inferred", "the SX126x and SX127x SF6 framings differ by 2 symbols only");
    }
  }
  for (const auto& e : em) {
    if (e.label.find("elrs") == std::string::npos || e.features.value("modulation", Json(nullptr)) != "lora") continue;
    const Json& f = e.features;
    const bool sub = e.center_mhz < 1500.0;
    const Json sf = f.contains("sf") ? f["sf"] : Json(nullptr);
    const double bw = f.contains("bw_hz") && f["bw_hz"].is_number() ? f["bw_hz"].get<double>() : 0.0;
    std::vector<std::string> ev = {"LoRa SF" + py_str(sf) + " " + fmt("%g", bw / 1e3) + " kHz " + (sub ? "sub-GHz" : "2.4 GHz")};
    if (f.contains("elrs_channel") && f["elrs_channel"].is_array() && !f["elrs_channel"].empty())
      ev.push_back("domain " + f["elrs_channel"][0].get<std::string>() + " (regulatory build)");
    if (sub && sf.is_number_integer() && sf.get<int64_t>() == 5) {
      ev.push_back("sub-GHz SF5 (200 Hz Full / 250 Hz) exists only on LR1121");
      return model_block("ExpressLRS", {"LR1121-based module"}, ev, "inferred");
    }
    return model_block("ExpressLRS", sub ? std::vector<std::string>{"SX127x-based module", "LR1121-based module"}
                                         : std::vector<std::string>{"SX1280-based module", "LR1121-based module"},
                       ev, "inferred", "this LoRa mode exists on both radio chips");
  }
  for (const auto& e : em)   // ExpressLRS LR1121 GFSK, verified SYNC (PLAN G3)
    if (e.label == "elrs_gfsk" && has(e))
      return model_block("ExpressLRS", {"LR1121-based module"},
                         {"300 kb/s GFSK SYNC packet verified (sync word " + e.decoded["sync_word"].get<std::string>() +
                              " = UID5 / UID4)",
                          "the GFSK modes exist only on LR1121 hardware"},
                         "decoded");
  bool flrc = false, gfsk = false, all24 = true;
  for (const auto& e : em) {
    flrc = flrc || e.label == "elrs_flrc_like";
    gfsk = gfsk || e.label == "elrs_gfsk_like";
    all24 = all24 && e.center_mhz >= 2400.0;
  }
  if (flrc) return model_block("ExpressLRS-like (FLRC)", {"SX1280-based module"}, {"650 kb/s FLRC-like rate"}, "inferred");
  if (gfsk && all24)
    return model_block("ExpressLRS-like (GFSK)", {"LR1121-based module"},
                       {"300 kb/s GFSK-like rate at 2.4 GHz (K1000 / DK modes: LR1121 only)"}, "inferred");
  return nullptr;
}

}  // namespace si::rx::rc_link::models
