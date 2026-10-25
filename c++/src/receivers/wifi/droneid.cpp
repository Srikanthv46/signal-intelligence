// Port of python receivers/wifi/droneid.py (see droneid.hpp).
#include "receivers/wifi/droneid.hpp"

#include <cmath>
#include <cstring>
#include <map>
#include <regex>
#include <utility>

#include "core/pyfmt.hpp"
#include "core/remote_id.hpp"
#include "receivers/wifi/pyutil.hpp"

namespace si::rx::wifi::droneid {

using detail::Bytes;
using detail::jget;
using detail::truthy;

namespace {

constexpr size_t kMsgLen = 25;
const char* kOdidOui = "fa:0b:bc";
constexpr int kOdidVendorType = 0x0D;
const double kDjiScale = 174533.0;

const char* lookup(const std::map<int, const char*>& m, int k) {
  auto it = m.find(k);
  return it == m.end() ? nullptr : it->second;
}
// TABLE.get(k, k)
Json get_or_int(const std::map<int, const char*>& m, int k) {
  const char* s = lookup(m, k);
  return s ? Json(s) : Json(k);
}

const std::map<int, const char*> kDjiModels = {
    {1, "Inspire 1"}, {2, "Phantom 3 Series"}, {3, "Phantom 3 Series"}, {4, "Phantom 3 Std"}, {5, "M100"},
    {6, "ACEONE"}, {7, "WKM"}, {8, "NAZA"}, {9, "A2"}, {10, "A3"}, {11, "Phantom 4"}, {12, "MG1"}, {14, "M600"},
    {15, "Phantom 3 4k"}, {16, "Mavic Pro"}, {17, "Inspire 2"}, {18, "Phantom 4 Pro"}, {20, "N2"},
    {21, "Spark"}, {23, "M600 Pro"}, {24, "Mavic Air"}, {25, "M200"}, {26, "Phantom 4 Series"},
    {27, "Phantom 4 Adv"}, {28, "M210"}, {30, "M210RTK"}, {31, "A3_AG"}, {32, "MG2"}, {34, "MG1A"},
    {35, "Phantom 4 RTK"}, {36, "Phantom 4 Pro V2.0"}, {38, "MG1P"}, {40, "MG1P-RTK"}, {41, "Mavic 2"},
    {44, "M200 V2 Series"}, {51, "Mavic 2 Enterprise"}, {53, "Mavic Mini"}, {58, "Mavic Air 2"},
    {59, "P4M"}, {60, "M300 RTK"}, {61, "DJI FPV"}, {63, "Mini 2"}, {64, "AGRAS T10"}, {65, "AGRAS T30"},
    {66, "Air 2S"}, {67, "M30"}, {68, "DJI Mavic 3"}, {69, "Mavic 2 Enterprise Advanced"}, {70, "Mini SE"}};

int32_t s32(const Bytes& b, size_t i) {
  uint32_t v = static_cast<uint32_t>(b[i]) | (static_cast<uint32_t>(b[i + 1]) << 8) |
               (static_cast<uint32_t>(b[i + 2]) << 16) | (static_cast<uint32_t>(b[i + 3]) << 24);
  return static_cast<int32_t>(v);
}
uint32_t u32(const Bytes& b, size_t i) { return static_cast<uint32_t>(s32(b, i)); }
int u16(const Bytes& b, size_t i) { return detail::u16le(b, i); }
int16_t s16(const Bytes& b, size_t i) { return static_cast<int16_t>(u16(b, i)); }

// b.split(b"\x00", 1)[0].decode("utf-8", "replace").strip()
std::string text(const Bytes& b) {
  size_t n = 0;
  while (n < b.size() && b.data()[n] != 0) ++n;
  return detail::py_strip(detail::decode_utf8(b.sub(0, n)));
}

// metres (0.5 m steps, -1000 m offset); 0 = not reported
}  // namespace

uint16_t crc16_kermit(const uint8_t* data, size_t n, uint16_t crc) {
  for (size_t i = 0; i < n; ++i) {
    crc ^= data[i];
    for (int k = 0; k < 8; ++k) crc = static_cast<uint16_t>((crc >> 1) ^ ((crc & 1) ? 0x8408 : 0));
  }
  return crc;
}

std::optional<Json> parse_odid(const std::vector<uint8_t>& bodyv) { return si::remote_id::parse_odid(bodyv); }
namespace {
std::optional<Json> dji_at(const Bytes& raw) {
  // struct "<BBBHH16siihhhhhhQiiiiBB20sH"
  const int version = raw[2];
  const int seq = u16(raw, 3), state = u16(raw, 5);
  const Bytes serial = raw.sub(7, 23);
  const int32_t lon_i = s32(raw, 23), lat_i = s32(raw, 27);
  const int16_t alt_i = s16(raw, 31), height_i = s16(raw, 33), vn = s16(raw, 35), ve = s16(raw, 37),
                vu = s16(raw, 39), yaw = s16(raw, 41);
  uint64_t gps = 0;
  for (int k = 7; k >= 0; --k) gps = (gps << 8) | static_cast<uint64_t>(raw[43 + k]);
  const int32_t olat_i = s32(raw, 51), olon_i = s32(raw, 55), hlon_i = s32(raw, 59), hlat_i = s32(raw, 63);
  const int model = raw[67];
  const Bytes uuid_b = raw.sub(69, 89);
  const int crc = u16(raw, 89);
  if (crc16_kermit(raw.data(), kDjiLen - 2) != crc) return std::nullopt;
  const double lat = lat_i / kDjiScale, lon = lon_i / kDjiScale;
  Json out = Json::object();
  out["standard"] = "DJI DroneID (proprietary)";
  out["transport"] = "Wi-Fi vendor element";
  out["crc_ok"] = true;
  out["version"] = version;
  out["sequence"] = seq;
  out["state_info"] = fmt("0x%04x", state);
  out["uas_id"] = text(serial);
  const char* mname = lookup(kDjiModels, model);
  out["model"] = mname ? std::string(mname) : "unknown type " + std::to_string(model);
  out["ua_type"] = "multirotor";
  if (std::fabs(lat) <= 90 && std::fabs(lon) <= 180 && (lat != 0 || lon != 0)) {
    out["lat"] = py_round(lat, 7);
    out["lon"] = py_round(lon, 7);
  }
  out["geodetic_altitude_m"] = py_round(alt_i / 3.281, 2);
  out["height_m"] = py_round(height_i / 3.281, 2);
  out["velocity_ne_up_mps"] = Json::array({vn, ve, vu});
  if (std::abs(static_cast<int>(yaw)) <= 36000) out["yaw_deg"] = py_round(yaw / 100.0, 2);
  if (gps) out["gps_time_ms"] = gps;
  const double olat = olat_i / kDjiScale, olon = olon_i / kDjiScale;
  if (std::fabs(olat) <= 90 && std::fabs(olon) <= 180 && (olat != 0 || olon != 0)) {
    out["operator_lat"] = py_round(olat, 7);
    out["operator_lon"] = py_round(olon, 7);
    out["operator_location_type"] = "live GNSS (operator moving)";
  }
  const double hlat = hlat_i / kDjiScale, hlon = hlon_i / kDjiScale;
  if (std::fabs(hlat) <= 90 && std::fabs(hlon) <= 180 && (hlat != 0 || hlon != 0)) {
    out["home_lat"] = py_round(hlat, 7);
    out["home_lon"] = py_round(hlon, 7);
  }
  std::string uuid = text(uuid_b);
  if (!uuid.empty()) out["uuid"] = uuid;
  return out;
}
}  // namespace

std::optional<Json> parse_dji(const std::vector<uint8_t>& bodyv) {
  const Bytes body(bodyv);
  const size_t last = body.size() > kDjiLen ? body.size() - kDjiLen : 0;
  for (size_t off = 0; off <= last; ++off) {
    const Bytes raw = body.sub(off, off + kDjiLen);
    if (raw.size() < kDjiLen) continue;  // struct.error
    auto got = dji_at(raw);
    if (got) {
      (*got)["payload_offset"] = off;
      return got;
    }
  }
  return std::nullopt;
}

std::optional<Json> from_nan_action(const std::vector<uint8_t>& b) {
  static const uint8_t kSid[6] = {0x88, 0x69, 0x19, 0x9D, 0x92, 0x09};   // SHA-256("org.opendroneid.remoteid")[:6]
  if (b.size() < 7 || b[0] != 4 || b[1] != 9 || b[2] != 0x50 || b[3] != 0x6F || b[4] != 0x9A || b[5] != 0x13)
    return std::nullopt;
  size_t i = 6;
  while (i + 3 <= b.size()) {
    const int aid = b[i];
    const size_t ln = static_cast<size_t>(b[i + 1] | (b[i + 2] << 8));
    const size_t v0 = i + 3, v1 = std::min(b.size(), v0 + ln);
    i = v0 + ln;
    if (aid != 0x03 || v1 - v0 < 9 || !std::equal(kSid, kSid + 6, b.begin() + static_cast<ptrdiff_t>(v0))) continue;
    const int ctrl = b[v0 + 8];
    size_t j = v0 + 9;
    if (ctrl & 0x40) j += 2;
    if ((ctrl & 0x04) && j < v1) j += 1 + b[j];
    if ((ctrl & 0x08) && j < v1) j += 1 + b[j];
    if ((ctrl & 0x10) && j < v1) {
      const size_t a = j + 1, e = std::min(v1, a + b[j]);
      if (a <= e) {
        auto rid = si::remote_id::parse_odid(std::vector<uint8_t>(b.begin() + static_cast<ptrdiff_t>(a),
                                                                  b.begin() + static_cast<ptrdiff_t>(e)),
                                             "Wi-Fi NAN service discovery");
        if (rid) return rid;
      }
    }
  }
  return std::nullopt;
}

std::optional<Json> from_vendor_ie(const std::string& oui, int vendor_type, const std::vector<uint8_t>& body) {
  const std::string o = detail::ascii_lower(oui);
  if (o == kOdidOui && vendor_type == kOdidVendorType) return parse_odid(body);
  if (o == "26:37:12") {
    std::vector<uint8_t> b;
    b.reserve(body.size() + 1);
    b.push_back(static_cast<uint8_t>(vendor_type));
    b.insert(b.end(), body.begin(), body.end());
    return parse_dji(b);
  }
  return std::nullopt;
}

std::optional<Json> scan_vendor_ies(const std::vector<uint8_t>& tlvv) {
  const Bytes tlv(tlvv);
  size_t i = 0;
  while (i + 2 <= tlv.size()) {
    const size_t ln = tlv[i + 1];
    const Bytes v = tlv.sub(i + 2, i + 2 + ln);
    if (tlv[i] == 221 && v.size() >= 4) {
      auto got = from_vendor_ie(detail::mac_str(v.sub(0, 3)), v[3], v.from(4).vec());
      if (got) return got;
    }
    i += 2 + ln;
  }
  return std::nullopt;
}

// ----------------------------------------------------------------------------- signatures
namespace {
// Python's re works on code points with Unicode \w; std::regex on bytes with ASCII \w. The
// string is lower-cased (ASCII) and every non-ASCII code point becomes one 'a' (a letter or
// digit: a word character) or ' ' (anything else), which keeps \b and ^ where Python has them.
bool cp_is_word(uint32_t c) {
  return c == 0xAA || c == 0xB5 || c == 0xBA || (c >= 0xC0 && c <= 0x24F && c != 0xD7 && c != 0xF7) ||
         (c >= 0x370 && c <= 0x3FF) || (c >= 0x400 && c <= 0x52F) || (c >= 0x5D0 && c <= 0x6FF) ||
         (c >= 0x900 && c <= 0xDFF) || (c >= 0x3040 && c <= 0x9FFF) || (c >= 0xAC00 && c <= 0xD7AF);
}

std::string regex_subject(const std::string& s) {
  std::string out;
  size_t i = 0;
  while (i < s.size()) {
    const unsigned char c = static_cast<unsigned char>(s[i]);
    if (c < 0x80) {
      out += (c >= 'A' && c <= 'Z') ? static_cast<char>(c - 'A' + 'a') : static_cast<char>(c);
      ++i;
      continue;
    }
    int n = c < 0xE0 ? 2 : c < 0xF0 ? 3 : 4;
    if (i + n > s.size()) n = 1;
    uint32_t cp = n == 2 ? (c & 0x1F) : n == 3 ? (c & 0x0F) : (c & 0x07);
    for (int k = 1; k < n; ++k) cp = (cp << 6) | (static_cast<unsigned char>(s[i + k]) & 0x3F);
    out += cp_is_word(cp) ? 'a' : ' ';
    i += n;
  }
  return out;
}

struct Pattern {
  std::regex re;
  const char* name;
};

const std::vector<Pattern>& vendor_patterns() {
  static const std::vector<Pattern> p = {
      {std::regex(R"(\bdji\b|\bsz dji\b|dji (osmo|baiwang|ronin))"), "DJI"},
      {std::regex(R"(\bparrot\b)"), "Parrot"},
      {std::regex(R"(autel (robotics|lntelligent|intelligent))"), "Autel Robotics"},
      {std::regex(R"(\bskydio\b)"), "Skydio"},
      {std::regex(R"(\byuneec\b)"), "Yuneec"},
      {std::regex(R"(\bhubsan\b)"), "Hubsan"},
      {std::regex(R"(powervision tech)"), "PowerVision"},
      {std::regex(R"(\bwalkera\b)"), "Walkera"},
      {std::regex(R"(\behang\b)"), "EHang"},
  };
  return p;
}

const std::vector<Pattern>& ssid_patterns() {
  static const std::vector<Pattern> p = {
      {std::regex(R"(^mavic[_-]?air[_-])"), "DJI Mavic Air"},
      {std::regex(R"(^mavic[_-]?mini[_-])"), "DJI Mavic Mini"},
      {std::regex(R"(^phantom3[_-])"), "DJI Phantom 3 Standard / 4K (Wi-Fi video models)"},
      {std::regex(R"(^(mavic|mavic[_-]?(air|mini|pro|2|3))[_-])"), "DJI Mavic"},
      {std::regex(R"(^phantom[_-]?\d*)"), "DJI Phantom"},
      {std::regex(R"(^spark[_-]rc)"), "DJI Spark remote controller"},
      {std::regex(R"(^spark[_-])"), "DJI Spark"},
      {std::regex(R"(^tello[_-])"), "DJI Tello"},
      {std::regex(R"(^dji[_-])"), "DJI"},
      {std::regex(R"(^(mini)[_-]?[23][_-])"), "DJI Mini"},
      {std::regex(R"(^anafi[_-])"), "Parrot Anafi"},
      {std::regex(R"(^(bebop|anafi|disco|parrot)[_-])"), "Parrot"},
      {std::regex(R"(^autel[_-])"), "Autel Robotics"},
      {std::regex(R"(^cgo3p?_)"), "Yuneec CGO3 camera video link (Typhoon Q500 4K / Typhoon H)"},
      {std::regex(R"(^(typhoon|yuneec)[_-])"), "Yuneec"},
      {std::regex(R"(^(hubsan|h501|h107))"), "Hubsan"},
      {std::regex(R"(^skydio[_-])"), "Skydio"},
      {std::regex(R"(^(fimi|xiaomi[_-]fimi)[_-])"), "FIMI"},
      {std::regex(R"(^(openhd|wfb|wifibroadcast))"), "open-source FPV link (OpenHD / wfb-ng)"},
      // PLAN W3 (research file 18): vendor-documented default hotspot names
      {std::regex(R"(^rmtt-)"), "DJI RoboMaster TT (Tello Talent)"},
      {std::regex(R"(^holystone ?fpv)"), "Holy Stone (Wi-Fi FPV: HS720 / HS100 series)"},
      {std::regex(R"(^ruko[_-]f11)"), "Ruko F11"},
      {std::regex(R"(^ruko[_-])"), "Ruko"},
      {std::regex(R"(^wifi-720p-)"), "toy Wi-Fi FPV drone (Eachine E58 class)"},
      {std::regex(R"(^expresslrs (tx|rx)$)"), "ExpressLRS radio in Wi-Fi setup mode"},
      {std::regex(R"(^(ardupilot|pixracer)$)"), "ArduPilot Wi-Fi telemetry bridge (MAVESP8266)"},
      {std::regex(R"(^cuavwlink)"), "CUAV PW-Link Wi-Fi telemetry"},
      {std::regex(R"(^dronebridge)"), "DroneBridge ESP32 telemetry"},
  };
  return p;
}

std::string first_word(const std::string& s) {
  auto w = detail::py_split_ws(s);
  return w.empty() ? std::string() : w[0];
}
}  // namespace

std::optional<Json> signature(const std::optional<std::string>& vendor_name, const std::optional<std::string>& ssid) {
  Json out = Json::object();
  if (vendor_name && !vendor_name->empty()) {
    const std::string v = regex_subject(*vendor_name);
    for (const auto& p : vendor_patterns()) {
      if (std::regex_search(v, p.re)) {
        out = Json::object();
        out["platform"] = p.name;
        out["evidence"] = "IEEE OUI registered to " + *vendor_name;
        out["confidence"] = 0.75;
        break;
      }
    }
  }
  if (ssid && !ssid->empty()) {
    const std::string s = regex_subject(*ssid);
    for (const auto& p : ssid_patterns()) {
      if (std::regex_search(s, p.re)) {
        const std::string name = p.name;
        if (!out.empty() && first_word(out["platform"].get<std::string>()) == first_word(name)) {
          const std::string cur = out["platform"].get<std::string>();
          out["platform"] = detail::len_codepoints(name) > detail::len_codepoints(cur) ? name : cur;
          out["evidence"] = out["evidence"].get<std::string>() + "; SSID '" + *ssid + "' matches " + name;
          out["confidence"] = 0.9;
        } else if (out.empty()) {
          out["platform"] = name;
          out["evidence"] = "SSID '" + *ssid + "' matches " + name;
          out["confidence"] = 0.7;
        }
        break;
      }
    }
  }
  if (out.empty()) return std::nullopt;
  return out;
}

std::string summary(const Json& rid) {
  std::vector<std::string> s = {"REMOTE ID"};
  if (truthy(jget(rid, "uas_id"))) s.push_back("UAS " + detail::pstr(rid["uas_id"]));
  if (truthy(jget(rid, "model"))) {
    s.push_back("(" + detail::pstr(rid["model"]) + ")");
  } else if (truthy(jget(rid, "ua_type")) && rid["ua_type"] != "undeclared") {
    s.push_back("(" + detail::pstr(rid["ua_type"]) + ")");
  }
  if (truthy(jget(rid, "status"))) s.push_back(detail::ascii_upper(detail::pstr(rid["status"])));
  if (!jget(rid, "lat").is_null())
    s.push_back(fmt("at %.6f,%.6f", rid["lat"].get<double>(), rid["lon"].get<double>()));
  for (const auto& [k, lab] : {std::pair<const char*, const char*>{"height_m", "AGL"}, {"geodetic_altitude_m", "alt"}}) {
    if (!jget(rid, k).is_null()) {
      s.push_back(fmt("%.0f m %s", rid[k].get<double>(), lab));
      break;
    }
  }
  if (!jget(rid, "speed_mps").is_null()) s.push_back(fmt("%.1f m/s", rid["speed_mps"].get<double>()));
  if (!jget(rid, "track_deg").is_null()) s.push_back("hdg " + detail::pstr(rid["track_deg"]));
  if (!jget(rid, "operator_lat").is_null())
    s.push_back(fmt("OPERATOR at %.6f,%.6f", rid["operator_lat"].get<double>(), rid["operator_lon"].get<double>()));
  if (truthy(jget(rid, "operator_id"))) s.push_back("operator ID " + detail::pstr(rid["operator_id"]));
  if (truthy(jget(rid, "description"))) s.push_back("'" + detail::pstr(rid["description"]) + "'");
  std::string out;
  for (size_t i = 0; i < s.size(); ++i) out += (i ? " " : "") + s[i];
  return out;
}

}  // namespace si::rx::wifi::droneid
