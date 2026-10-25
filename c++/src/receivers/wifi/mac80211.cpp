// Port of python receivers/wifi/mac80211.py (see mac80211.hpp). Python's dict semantics are
// kept: keys in insertion order (reassignment keeps the position), None values present as
// null, and IndexError raised by bytes indexing caught exactly where Python catches it.
#include "receivers/wifi/mac80211.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <optional>
#include <map>
#include <set>
#include <utility>

#include "core/mavlink.hpp"
#include "core/pyfmt.hpp"
#include "receivers/wifi/droneid.hpp"
#include "receivers/wifi/oui_table.hpp"
#include "receivers/wifi/pyutil.hpp"

namespace si::rx::wifi::mac80211 {

using detail::ascii_lower;
using detail::be16;
using detail::Bytes;
using detail::decode_ascii;
using detail::IndexError;
using detail::jget;
using detail::mac_str;
using detail::pstr;
using detail::py_strip;
using detail::truthy;
using detail::u16le;

namespace {

using IntMap = std::map<int, const char*>;

const char* lookup(const IntMap& m, int k) {
  auto it = m.find(k);
  return it == m.end() ? nullptr : it->second;
}
Json get_or_int(const IntMap& m, int k) {
  const char* s = lookup(m, k);
  return s ? Json(s) : Json(k);
}
Json get_or_null(const IntMap& m, int k) {
  const char* s = lookup(m, k);
  return s ? Json(s) : Json(nullptr);
}
std::string get_or(const IntMap& m, int k, const std::string& dflt) {
  const char* s = lookup(m, k);
  return s ? std::string(s) : dflt;
}

const IntMap kFtype = {{0, "mgmt"}, {1, "ctrl"}, {2, "data"}, {3, "ext"}};
const std::map<std::pair<int, int>, const char*> kSubtype = {
    {{0, 0}, "assoc-req"}, {{0, 1}, "assoc-resp"}, {{0, 2}, "reassoc-req"}, {{0, 3}, "reassoc-resp"},
    {{0, 4}, "probe-req"}, {{0, 5}, "probe-resp"}, {{0, 6}, "timing-adv"}, {{0, 8}, "beacon"}, {{0, 9}, "atim"},
    {{0, 10}, "disassoc"}, {{0, 11}, "auth"}, {{0, 12}, "deauth"}, {{0, 13}, "action"}, {{0, 14}, "action-no-ack"},
    {{1, 2}, "trigger"}, {{1, 4}, "bf-report-poll"}, {{1, 5}, "ndp-announcement"}, {{1, 6}, "ctrl-frame-ext"},
    {{1, 7}, "ctrl-wrapper"}, {{1, 8}, "block-ack-req"}, {{1, 9}, "block-ack"}, {{1, 10}, "ps-poll"},
    {{1, 11}, "rts"}, {{1, 12}, "cts"}, {{1, 13}, "ack"}, {{1, 14}, "cf-end"}, {{1, 15}, "cf-end+cf-ack"},
    {{2, 0}, "data"}, {{2, 1}, "data+cf-ack"}, {{2, 2}, "data+cf-poll"}, {{2, 3}, "data+cf-ack+cf-poll"},
    {{2, 4}, "null"}, {{2, 5}, "cf-ack"}, {{2, 6}, "cf-poll"}, {{2, 7}, "cf-ack+cf-poll"}, {{2, 8}, "qos-data"},
    {{2, 9}, "qos-data+cf-ack"}, {{2, 10}, "qos-data+cf-poll"}, {{2, 11}, "qos-data+cf-ack+cf-poll"},
    {{2, 12}, "qos-null"}, {{2, 14}, "qos-cf-poll"}, {{2, 15}, "qos-cf-ack+cf-poll"},
    {{3, 0}, "dmg-beacon"}, {{3, 1}, "s1g-beacon"}};
const std::vector<std::pair<int, const char*>> kCaps = {
    {0, "ESS"}, {1, "IBSS"}, {4, "privacy"}, {5, "short-preamble"}, {8, "spectrum-mgmt"}, {9, "QoS"},
    {10, "short-slot"}, {11, "APSD"}, {12, "radio-measurement"}, {13, "EPD"}, {14, "delayed-BA"},
    {15, "immediate-BA"}};
const IntMap kCipher = {{0, "use-group"}, {1, "WEP-40"}, {2, "TKIP"}, {4, "CCMP-128"}, {5, "WEP-104"},
                        {6, "BIP-CMAC-128"}, {7, "group-traffic-not-allowed"}, {8, "GCMP-128"},
                        {9, "GCMP-256"}, {10, "CCMP-256"}, {11, "BIP-GMAC-128"}, {12, "BIP-GMAC-256"},
                        {13, "BIP-CMAC-256"}};
const IntMap kAkm = {{1, "802.1X"}, {2, "PSK"}, {3, "FT-802.1X"}, {4, "FT-PSK"}, {5, "802.1X-SHA256"},
                     {6, "PSK-SHA256"}, {7, "TDLS"}, {8, "SAE"}, {9, "FT-SAE"}, {10, "AP-PeerKey"},
                     {11, "802.1X-SuiteB"}, {12, "802.1X-SuiteB-192"}, {13, "FT-802.1X-SHA384"},
                     {14, "FILS-SHA256"}, {15, "FILS-SHA384"}, {16, "FT-FILS-SHA256"}, {17, "FT-FILS-SHA384"},
                     {18, "OWE"}, {19, "FT-PSK-SHA384"}, {20, "PSK-SHA384"}, {22, "802.1X-SHA384"},
                     {24, "SAE-EXT-KEY"}, {25, "FT-SAE-EXT-KEY"}};
const IntMap kAuthAlg = {{0, "open"}, {1, "shared-key"}, {2, "FT"}, {3, "SAE"}, {4, "FILS-SK"},
                         {5, "FILS-SK-PFS"}, {6, "FILS-PK"}, {7, "PASN"}, {0xFFFF, "vendor"}};
const IntMap kReason = {
    {1, "unspecified"}, {2, "previous authentication no longer valid"}, {3, "station leaving (deauthenticated)"},
    {4, "inactivity"}, {5, "AP cannot handle all associated stations"},
    {6, "class 2 frame from non-authenticated STA"}, {7, "class 3 frame from non-associated STA"},
    {8, "station leaving BSS (disassociated)"}, {9, "not authenticated"}, {10, "power capability unacceptable"},
    {11, "supported channels unacceptable"}, {12, "BSS transition management"}, {13, "invalid element"},
    {14, "MIC failure"}, {15, "4-way handshake timeout"}, {16, "group key handshake timeout"},
    {17, "element in 4-way handshake differs"}, {18, "invalid group cipher"}, {19, "invalid pairwise cipher"},
    {20, "invalid AKMP"}, {21, "unsupported RSNE version"}, {22, "invalid RSNE capabilities"},
    {23, "802.1X authentication failed"}, {24, "cipher suite rejected"}, {25, "TDLS teardown unreachable"},
    {26, "TDLS teardown unspecified"}, {34, "excessive frame loss (poor channel)"},
    {36, "requesting STA leaving BSS"}, {39, "peer timeout"}, {45, "peer key mismatch"}, {66, "invalid PMKID"}};
const IntMap kStatus = {
    {0, "success"}, {1, "unspecified failure"}, {10, "capabilities mismatch"}, {11, "reassociation denied"},
    {12, "association denied"}, {13, "authentication algorithm not supported"}, {14, "auth sequence out of order"},
    {15, "challenge failure"}, {16, "authentication timeout"}, {17, "AP full"}, {18, "basic rates mismatch"},
    {30, "rejected temporarily (try later)"}, {31, "robust management frame policy violation"},
    {37, "request declined"}, {40, "invalid element"}, {41, "invalid group cipher"}, {42, "invalid pairwise cipher"},
    {43, "invalid AKMP"}, {53, "invalid PMKID"}, {72, "invalid SAE group"}, {76, "anti-clogging token required"},
    {77, "finite cyclic group not supported"}, {82, "rejected with suggested BSS transition"},
    {123, "SAE hash-to-element required"}, {126, "SAE public key"}};
const IntMap kActionCat = {
    {0, "spectrum-mgmt"}, {1, "QoS"}, {2, "DLS"}, {3, "block-ack"}, {4, "public"}, {5, "radio-measurement"},
    {6, "fast-BSS-transition"}, {7, "HT"}, {8, "SA-query"}, {9, "protected-dual-of-public"}, {10, "WNM"},
    {11, "unprotected-WNM"}, {12, "TDLS"}, {13, "mesh"}, {14, "multihop"}, {15, "self-protected"}, {16, "DMG"},
    {18, "fast-session-transfer"}, {19, "robust-AV-streaming"}, {20, "unprotected-DMG"}, {21, "VHT"},
    {30, "HE"}, {31, "protected-HE"}, {33, "protected-FTM"}, {36, "EHT"}, {37, "protected-EHT"},
    {126, "vendor-protected"}, {127, "vendor"}};
const std::map<std::string, const char*> kVendorIe = {
    {"00:50:f2", "Microsoft"}, {"00:10:18", "Broadcom"}, {"00:03:7f", "Qualcomm Atheros"},
    {"00:0c:43", "MediaTek (Ralink)"}, {"00:0c:e7", "MediaTek"}, {"00:e0:4c", "Realtek"}, {"00:17:f2", "Apple"},
    {"50:6f:9a", "Wi-Fi Alliance"}, {"00:90:4c", "Broadcom (Epigram)"}, {"8c:fd:f0", "Qualcomm"},
    {"00:13:92", "Ruckus"}, {"00:0b:86", "Aruba"}, {"00:40:96", "Cisco"}, {"00:a0:40", "Apple"},
    {"00:26:86", "Quantenna"}, {"f4:f5:e8", "Google"}, {"00:13:74", "Atheros"}, {"00:15:6d", "Ubiquiti"},
    {"00:25:00", "Apple"}, {"00:1d:0f", "TP-Link"}, {"ac:85:3d", "Huawei"}, {"00:e0:fc", "Huawei"},
    {"00:19:3b", "LG"}, {"00:16:32", "Samsung"}, {"00:00:f0", "Samsung"}};
const IntMap kWfaTypes = {{9, "Wi-Fi Direct (P2P)"}, {0x10, "Hotspot 2.0"}, {0x12, "OSEN"}, {0x16, "MBO/OCE"},
                          {0x1C, "OWE transition"}, {0x1A, "Wi-Fi Easy Connect (DPP)"},
                          {0x1B, "Multi-AP (EasyMesh)"}};
const IntMap kWpsAttr = {{0x1011, "device_name"}, {0x1021, "manufacturer"}, {0x1023, "model_name"},
                         {0x1024, "model_number"}, {0x1042, "serial_number"}};
const IntMap kWpsDeviceCat = {{1, "computer"}, {2, "input device"}, {3, "printer/scanner"}, {4, "camera"},
                              {5, "storage"}, {6, "network infrastructure (AP/router)"}, {7, "display"},
                              {8, "multimedia"}, {9, "gaming"}, {10, "telephone"}, {11, "audio"},
                              {12, "docking"}, {255, "other"}};
const IntMap kEthertype = {{0x0800, "IPv4"}, {0x0806, "ARP"}, {0x86DD, "IPv6"}, {0x888E, "EAPOL"},
                           {0x88C7, "802.11 pre-auth"}, {0x8100, "VLAN"}, {0x88CC, "LLDP"},
                           {0x893A, "IEEE 1905 (EasyMesh)"}, {0x890D, "802.11 TDLS"}, {0x88B4, "WAI (WAPI)"}};
const IntMap kIpProto = {{1, "ICMP"}, {2, "IGMP"}, {6, "TCP"}, {17, "UDP"}, {58, "ICMPv6"}};
const IntMap kPorts = {{53, "DNS"}, {67, "DHCP"}, {68, "DHCP"}, {123, "NTP"}, {137, "NetBIOS"}, {138, "NetBIOS"},
                       {1900, "SSDP"}, {5353, "mDNS"}, {5355, "LLMNR"}, {443, "HTTPS"}, {80, "HTTP"},
                       {546, "DHCPv6"}, {547, "DHCPv6"}, {14550, "MAVLink (GCS)"}, {14555, "MAVLink"},
                       {14540, "MAVLink (offboard)"}, {5760, "MAVLink (TCP)"}, {8889, "Tello SDK command"},
                       {8890, "Tello SDK state"}, {11111, "Tello video"}};

// wfb-ng / OpenIPC / OpenHD raw-802.11 links (python _wfb; research, PLAN G1)
constexpr int kWfbDefaultLinkId = 0x7505D6;
std::string wfb_stream(int port) {
  const int base = port & 0x7F;
  std::string s = base < 16 ? "video / data" : base < 32 ? "MAVLink" : base < 48 ? "tunnel" : base == 0x7F ? "bind" : "other";
  return s + ((port & 0x80) ? " (reverse: ground -> air)" : "");
}
Json wfb(const Bytes& a2, const Bytes& body) {
  const int link = (a2[2] << 16) | (a2[3] << 8) | a2[4], port = a2[5];
  Json out = Json::object();
  out["link_id"] = link;
  out["link_id_hex"] = fmt("0x%06X", link);
  out["default_link_domain"] = link == kWfbDefaultLinkId;
  out["radio_port"] = port;
  out["stream"] = wfb_stream(port);
  if (!body.empty()) {
    out["packet_type"] = body[0] == 1 ? Json("data") : body[0] == 2 ? Json("session") : Json(fmt("0x%02x", body[0]));
    if (body[0] == 1 && body.size() >= 9) {
      const uint64_t n = be_uint(body.sub(1, 9));
      out["block_idx"] = n >> 8;
      out["fragment_idx"] = n & 0xFF;
    }
  }
  return out;
}

// ESP-NOW + the ExpressLRS Backpack (python crc8_dvb_s2 / _msp / _espnow; research, PLAN G7)
int crc8_dvb_s2(const Bytes& d) {
  unsigned c = 0;
  for (size_t i = 0; i < d.size(); ++i) {
    c ^= static_cast<unsigned>(d[i]);
    for (int k = 0; k < 8; ++k) c = (c & 0x80) ? ((c << 1) ^ 0xD5u) & 0xFF : (c << 1) & 0xFF;
  }
  return static_cast<int>(c);
}
const IntMap kMspNames = {{0x0009, "ELRS bind"}, {0x000B, "ELRS request VTX packet"}, {0x0010, "backpack version"},
                          {0x0011, "backpack CRSF telemetry"}, {0x0058, "VTX config (get)"}, {0x0059, "VTX config (set)"},
                          {0x00B5, "ELRS set send UID"}, {0x00B6, "ELRS set OSD"}, {0x0300, "get channel index"},
                          {0x0301, "set channel index"}, {0x0302, "get frequency"}, {0x0303, "set frequency"},
                          {0x0304, "get recording state"}, {0x0305, "set recording state"}, {0x0306, "get VRX mode"},
                          {0x0307, "set VRX mode"}, {0x0308, "get RSSI"}, {0x0309, "get battery voltage"},
                          {0x030A, "get firmware"}, {0x030D, "set head tracking"}, {0x030E, "set RTC"},
                          {0x0380, "set mode"}, {0x0381, "get version"}, {0x0382, "get status"},
                          {0x0383, "set pan / tilt / roll"}};
Json msp(const Bytes& b) {
  if (b.size() < 9 || b[0] != '$' || b[1] != 'X' || (b[2] != '<' && b[2] != '>')) return nullptr;
  const int func = b[4] | (b[5] << 8), size = b[6] | (b[7] << 8);
  if (b.size() < static_cast<size_t>(9 + size)) return nullptr;
  const Bytes pay = b.sub(8, 8 + static_cast<size_t>(size));
  Json out = Json::object();
  out["direction"] = b[2] == '<' ? "command" : "response";
  out["function"] = fmt("0x%04X", func);
  const char* nm = lookup(kMspNames, func);
  out["function_name"] = nm ? nm : "?";
  out["payload_hex"] = hex(pay);
  out["crc_ok"] = crc8_dvb_s2(b.sub(3, 8 + static_cast<size_t>(size))) == b[8 + static_cast<size_t>(size)];
  int idx = -1;
  if (func == 0x0059 && pay.size() >= 2) {
    const int v = pay[0] | (pay[1] << 8);
    if (v < 48) idx = v;
    else out["vtx_frequency_mhz"] = v;
  } else if (func == 0x0301 && pay.size() >= 1) {
    idx = pay[0];
  }
  if (idx >= 0 && idx < 48) {
    out["vtx_band"] = std::string(1, "ABEFRL"[idx / 8]);
    out["vtx_channel"] = idx % 8 + 1;
  }
  return out;
}
Json espnow(const Bytes& b, Json& info) {
  Json out = Json::object();
  out["random"] = hex(b.sub(4, 8));
  size_t i = 8;
  while (i + 2 <= b.size()) {
    const int eid = b[i], ln = b[i + 1];
    const Bytes v = b.sub(i + 2, i + 2 + static_cast<size_t>(ln));
    if (eid == 221 && v.size() >= 5 && mac_str(v.sub(0, 3)) == "18:fe:34" && v[3] == 4) {
      out["version"] = v[4];
      out["body_len"] = v.size() - 5;
      Json m = msp(v.from(5));
      if (!m.is_null()) out["msp"] = m;
      break;
    }
    i += 2 + static_cast<size_t>(ln);
  }
  const Json& m = jget(out, "msp");
  const Json& a1 = jget(info, "addr1");
  const Json& a2 = jget(info, "addr2");
  if (truthy(m) && m["crc_ok"].get<bool>() && truthy(a2) &&
      (a1 == a2 || (a1 == "ff:ff:ff:ff:ff:ff" && m["function"] == "0x0009"))) {
    Json bp = Json::object();
    bp["bind_uid"] = a2;
    bp["msp"] = m["function_name"];
    for (const char* k : {"vtx_band", "vtx_channel", "vtx_frequency_mhz"})
      if (m.contains(k)) bp[k] = m[k];
    info["elrs_backpack"] = bp;
  }
  return out;
}

// python _mavlink(data): CRC-valid MAVLink frames in a datagram -> {frames, summary...} or null
Json mavlink_frames(const Bytes& data) {
  if (data.size() < 8) return nullptr;
  Json ok = Json::array();
  for (const auto& f : mavlink::scan(data.vec()))
    if (f.contains("crc_ok") && f["crc_ok"].is_boolean() && f["crc_ok"].get<bool>()) ok.push_back(f);
  if (ok.empty()) return nullptr;
  Json o = Json::object();
  o["frames"] = ok.size();
  const Json sm = mavlink::summary(ok);
  for (auto it = sm.begin(); it != sm.end(); ++it) o[it.key()] = it.value();
  return o;
}

// python _udp_app (research, PLAN W4): Tello SDK command (UDP 8889) / state (UDP 8890), and video streams (raw
// H.264 / H.265 Annex B, or RTP with a dynamic payload type carrying H.264 NAL units)
const std::set<std::string>& tello_cmds() {
  static const std::set<std::string> s = {
      "command", "takeoff", "land", "streamon", "streamoff", "emergency", "up", "down", "left", "right", "forward",
      "back", "cw", "ccw", "flip", "go", "stop", "curve", "jump", "speed", "rc", "wifi", "mon", "moff", "mdirection",
      "ap", "motoron", "motoroff", "throwfly", "reboot", "setfps", "setbitrate", "setresolution", "port", "multiwifi",
      "downvision", "battery?", "speed?", "time?", "wifi?", "sdk?", "sn?", "hardware?", "active?", "ap?", "ssid?"};
  return s;
}
const std::vector<std::pair<std::string, std::string>>& tello_state() {
  static const std::vector<std::pair<std::string, std::string>> s = {
      {"pitch", "pitch_deg"}, {"roll", "roll_deg"}, {"yaw", "yaw_deg"}, {"vgx", "vgx_dm_s"}, {"vgy", "vgy_dm_s"},
      {"vgz", "vgz_dm_s"}, {"templ", "temp_low_c"}, {"temph", "temp_high_c"}, {"tof", "tof_cm"}, {"h", "height_cm"},
      {"bat", "battery_pct"}, {"baro", "baro_m"}, {"time", "motor_time_s"}, {"agx", "agx"}, {"agy", "agy"},
      {"agz", "agz"}, {"mid", "mission_pad"}};
  return s;
}
const char* h264_nal(int t) {
  switch (t) {
    case 1: return "slice"; case 5: return "IDR slice"; case 6: return "SEI"; case 7: return "SPS";
    case 8: return "PPS"; case 9: return "AUD"; default: return nullptr;
  }
}
const char* h265_nal(int t) {
  switch (t) {
    case 1: return "slice"; case 19: case 20: return "IDR slice"; case 32: return "VPS"; case 33: return "SPS";
    case 34: return "PPS"; case 35: return "AUD"; case 39: return "SEI"; default: return nullptr;
  }
}
// python int(v) / float(v) on an ASCII value: surrounding whitespace allowed, the whole rest must parse
std::optional<Json> py_number(const std::string& v) {
  const std::string t = py_strip(v);
  if (t.empty()) return std::nullopt;
  if (v.find('.') != std::string::npos) {
    for (char c : t)
      if (!(std::isdigit(static_cast<unsigned char>(c)) || c == '.' || c == '-' || c == '+' || c == 'e' || c == 'E'))
        return std::nullopt;
    char* end = nullptr;
    const double x = std::strtod(t.c_str(), &end);
    if (end != t.c_str() + t.size()) return std::nullopt;
    return Json(x);
  }
  size_t i = (t[0] == '-' || t[0] == '+') ? 1 : 0;
  if (i == t.size()) return std::nullopt;
  for (size_t k = i; k < t.size(); ++k)
    if (!std::isdigit(static_cast<unsigned char>(t[k]))) return std::nullopt;
  return Json(std::stoll(t));
}
Json udp_app(const Bytes& d, int sp, int dp) {
  Json out = Json::object();
  if ((sp == 8889 || dp == 8889) && !d.empty() && d.size() <= 128) {
    bool ascii = true;
    for (size_t i = 0; i < d.size(); ++i) ascii = ascii && d[i] < 0x80;
    std::string t = ascii ? py_strip(std::string(d.data(), d.data() + d.size())) : std::string();
    bool printable = !t.empty();
    for (char c : t) printable = printable && c >= 0x20 && c < 0x7F;
    const std::string w = printable ? ascii_lower(t.substr(0, t.find(' '))) : std::string();
    std::string nodot = w;
    if (const size_t k = nodot.find('.'); k != std::string::npos) nodot.erase(k, 1);
    bool digits = !nodot.empty();
    for (char c : nodot) digits = digits && std::isdigit(static_cast<unsigned char>(c));
    if (tello_cmds().count(w) || (sp == 8889 && (w == "ok" || w == "error" || digits))) {
      if (w == "wifi" || w == "ap" || w == "multiwifi") t = w + " <credentials not reported>";
      out["tello"] = {{"port", "command"}, {"direction", dp == 8889 ? "to drone" : "from drone"}, {"text", t}};
    }
  }
  if ((sp == 8890 || dp == 8890) && d.size() <= 300) {
    Json kv = Json::object();
    const std::string s = py_strip(decode_ascii(d));
    size_t a = 0;
    while (a <= s.size()) {
      size_t b = s.find(';', a);
      if (b == std::string::npos) b = s.size();
      const std::string item = s.substr(a, b - a);
      const size_t c = item.find(':');
      const std::string k = item.substr(0, c), v = c == std::string::npos ? std::string() : item.substr(c + 1);
      for (const auto& [key, name] : tello_state())
        if (key == k && !v.empty()) {
          if (auto n = py_number(v)) kv[name] = *n;
        }
      a = b + 1;
    }
    if (kv.size() >= 8 && kv.contains("battery_pct") && kv.contains("height_cm") && kv.contains("pitch_deg")) {
      Json te = {{"port", "state"}, {"direction", "from drone"}};
      for (auto it = kv.begin(); it != kv.end(); ++it) te[it.key()] = it.value();
      out["tello"] = te;
    }
  }
  const size_t sc = (d.size() >= 4 && d[0] == 0 && d[1] == 0 && d[2] == 0 && d[3] == 1) ? 4
                    : (d.size() >= 3 && d[0] == 0 && d[1] == 0 && d[2] == 1)            ? 3 : 0;
  bool nal = false;
  if (sc && d.size() >= sc + 2 && !(d[sc] & 0x80)) {
    const int h = d[sc], t5 = (h >> 1) & 0x3F, t4 = h & 0x1F;
    if (!(h & 1) && h265_nal(t5) && d[sc + 1] >= 1 && d[sc + 1] <= 3) {
      out["video"] = {{"transport", "raw Annex B"}, {"codec", "H.265"}, {"nal", h265_nal(t5)}};
      nal = true;
    } else if (h264_nal(t4) && ((t4 != 7 && t4 != 8 && t4 != 5) || (h & 0x60))) {
      out["video"] = {{"transport", "raw Annex B"}, {"codec", "H.264"}, {"nal", h264_nal(t4)}};
      nal = true;
    }
  }
  if (!nal && d.size() >= 14 && (d[0] >> 6) == 2 && (d[1] & 0x7F) >= 96) {
    const size_t off = 12 + 4 * static_cast<size_t>(d[0] & 0x0F);
    if (d.size() > off + 1 && !(d[off] & 0x80)) {
      const int t = d[off] & 0x1F;
      if ((t >= 1 && t <= 23) || t == 24 || t == 28) {
        const char* n = h264_nal(t);
        out["video"] = {{"transport", "RTP"}, {"codec", "H.264-like"}, {"payload_type", d[1] & 0x7F},
                        {"nal", n ? std::string(n) : t == 24 ? "STAP-A" : t == 28 ? "FU-A" : "type " + std::to_string(t)}};
      }
    }
  }
  return out;
}

std::string uptime(uint64_t us) {
  uint64_t s = us / 1000000;
  const uint64_t d = s / 86400;
  s %= 86400;
  const uint64_t h = s / 3600;
  s %= 3600;
  const uint64_t m = s / 60;
  s %= 60;
  return fmt("%llud %02llu:%02llu:%02llu", static_cast<unsigned long long>(d), static_cast<unsigned long long>(h),
             static_cast<unsigned long long>(m), static_cast<unsigned long long>(s));
}

Json caps(int c) {
  Json out = Json::array();
  for (const auto& [b, n] : kCaps)
    if ((c >> b) & 1) out.push_back(n);
  return out;
}

// (oui string, type byte) of a 4-byte cipher / AKM suite; the type read raises IndexError
std::pair<std::string, int> suite(const Bytes& b, size_t i) {
  std::string o = mac_str(b.sub(i, i + 3));
  return {o, b[i + 3]};
}

Json rsn(const Bytes& v, const std::string& oui_expect = "00:0f:ac", bool wpa1 = false) {
  Json out = Json::object();
  out["version"] = u16le(v, 0);  // outside the try, as Python
  size_t i = 2;
  try {
    if (i + 4 <= v.size()) {
      auto [o, t] = suite(v, i);
      out["group_cipher"] = o == oui_expect ? get_or(kCipher, t, o + ":" + std::to_string(t))
                                            : o + ":" + std::to_string(t);
      i += 4;
    }
    if (i + 2 <= v.size()) {
      const size_t n = u16le(v, i);
      i += 2;
      Json pc = Json::array();
      for (size_t k = 0; k < n; ++k) pc.push_back(get_or(kCipher, suite(v, i + 4 * k).second, "?"));
      out["pairwise_ciphers"] = pc;
      i += 4 * n;
    }
    if (i + 2 <= v.size()) {
      const size_t n = u16le(v, i);
      i += 2;
      Json akm = Json::array();
      for (size_t k = 0; k < n; ++k) {
        auto [o, t] = suite(v, i + 4 * k);
        akm.push_back(o == oui_expect ? get_or(kAkm, t, std::to_string(t)) : oui_name(o) + ":" + std::to_string(t));
      }
      out["akm"] = akm;
      i += 4 * n;
    }
    if (!wpa1 && i + 2 <= v.size()) {
      const int cap = u16le(v, i);
      out["mfp_required"] = static_cast<bool>(cap & 0x40);
      out["mfp_capable"] = static_cast<bool>(cap & 0x80);
      i += 2;
      if (i + 2 <= v.size()) {
        const size_t n = u16le(v, i);
        out["n_pmkid"] = n;
        i += 2 + 16 * n;
      }
      if (i + 4 <= v.size()) out["group_mgmt_cipher"] = get_or(kCipher, v[i + 3], "?");
    }
  } catch (const IndexError&) {
    out["truncated"] = true;
  }
  return out;
}

Json wps(const Bytes& v) {
  Json out = Json::object();
  size_t i = 0;
  while (i + 4 <= v.size()) {
    const int t = (v[i] << 8) | v[i + 1];
    const size_t n = static_cast<size_t>((v[i + 2] << 8) | v[i + 3]);
    const Bytes val = v.sub(i + 4, i + 4 + n);
    if (const char* name = lookup(kWpsAttr, t)) {
      out[name] = detail::py_strip_chars(detail::decode_utf8(val), std::string("\0 ", 2));
    } else if (t == 0x1044 && n >= 1) {
      out["state"] = get_or_int({{1, "not configured"}, {2, "configured"}}, val[0]);
    } else if (t == 0x1047 && n == 16) {
      out["uuid"] = detail::hex(val);
    } else if (t == 0x1054 && n >= 8) {
      const int c = (val[0] << 8) | val[1];
      out["device_category"] = get_or_int(kWpsDeviceCat, c);
    } else if (t == 0x103B && n >= 1) {
      out["response_type"] =
          get_or_null({{0, "enrollee (open 802.1X)"}, {1, "enrollee"}, {2, "registrar"}, {3, "AP"}}, val[0]);
    } else if (t == 0x1041 && n >= 1) {
      out["selected_registrar"] = static_cast<bool>(val[0]);
    } else if (t == 0x1012 && n >= 2) {
      const int c = (val[0] << 8) | val[1];
      out["device_password_id"] = get_or_int({{0, "PIN"}, {4, "push-button"}}, c);
    }
    i += 4 + n;
  }
  return out;
}

Json p2p(const Bytes& v) {
  Json out = Json::object();
  size_t i = 0;
  while (i + 3 <= v.size()) {
    const int a = v[i];
    const size_t n = u16le(v, i + 1);
    const Bytes val = v.sub(i + 3, i + 3 + n);
    if (a == 13 && n >= 17) {
      out["device_address"] = mac_str(val.sub(0, 6));
      const size_t nsec = val[16];
      const size_t j = 17 + 8 * nsec;
      if (j + 4 <= val.size() && ((val[j] << 8) | val[j + 1]) == 0x1011) {
        const size_t ln = (val[j + 2] << 8) | val[j + 3];
        out["device_name"] = detail::decode_utf8(val.sub(j + 4, j + 4 + ln));
      }
      const int cat = (val[8] << 8) | val[9];
      out["device_category"] = get_or_int(kWpsDeviceCat, cat);
    } else if (a == 3 && n >= 6) {
      out["device_id"] = mac_str(val.sub(0, 6));
    }
    i += 3 + n;
  }
  return out;
}

// python sorted(set(list of floats))
Json sorted_unique(const std::vector<double>& v) {
  std::set<double> s(v.begin(), v.end());
  Json out = Json::array();
  for (double x : s) out.push_back(x);
  return out;
}

bool ext_bit(const Bytes& v, int b) {
  const size_t byte = static_cast<size_t>(b) / 8;
  return byte < v.size() && ((v.data()[byte] >> (b % 8)) & 1);
}

}  // namespace

// ------------------------------------------------------------------------------- registry
uint32_t crc32(const uint8_t* data, size_t n, uint32_t crc) {
  static const auto table = [] {
    std::vector<uint32_t> t(256);
    for (uint32_t i = 0; i < 256; ++i) {
      uint32_t c = i;
      for (int k = 0; k < 8; ++k) c = (c & 1) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
      t[i] = c;
    }
    return t;
  }();
  crc = ~crc;
  for (size_t i = 0; i < n; ++i) crc = table[(crc ^ data[i]) & 0xFF] ^ (crc >> 8);
  return ~crc;
}

bool fcs_ok(const std::vector<uint8_t>& m) {
  if (m.size() <= 4) return false;
  const size_t n = m.size() - 4;
  const uint32_t f = static_cast<uint32_t>(m[n]) | (static_cast<uint32_t>(m[n + 1]) << 8) |
                     (static_cast<uint32_t>(m[n + 2]) << 16) | (static_cast<uint32_t>(m[n + 3]) << 24);
  return crc32(m.data(), n) == f;
}

std::optional<std::string> vendor(const std::optional<std::string>& mac) {
  if (!mac || mac->empty()) return std::nullopt;
  const int o1 = std::stoi(mac->substr(0, 2), nullptr, 16);
  if (o1 & 0x01) return detail::ascii_lower(*mac) == "ff:ff:ff:ff:ff:ff" ? "broadcast" : "multicast";
  if (o1 & 0x02) return std::string("locally-administered (randomized/private)");
  std::string h;
  for (char c : *mac)
    if (c != ':') h += c;
  h = detail::ascii_upper(h);
  for (size_t n : {9, 7, 6}) {
    if (h.size() < n) continue;
    if (const char* v = detail::oui_registry(h.substr(0, n))) return std::string(v);
  }
  return std::string("unknown OUI");
}

std::string oui_name(const std::string& oui) {
  const std::string o = detail::ascii_lower(oui);
  auto it = kVendorIe.find(o);
  if (it != kVendorIe.end()) return it->second;
  std::string h;
  for (char c : o)
    if (c != ':') h += c;
  const char* v = detail::oui_registry(detail::ascii_upper(h));
  return v ? std::string(v) : o;
}

// ------------------------------------------------------------------------------- elements
Json parse_elements(const uint8_t* tlvp, size_t tn) {
  const Bytes tlv(tlvp, tn);
  Json out = Json::object();
  out["element_ids"] = Json::array();
  std::vector<double> rates, basic;
  Json vend = Json::array(), ext_ids = Json::array();
  size_t i = 0;
  while (i + 2 <= tlv.size()) {
    const int tag = tlv[i];
    const size_t ln = tlv[i + 1];
    const Bytes v = tlv.sub(i + 2, i + 2 + ln);
    if (v.size() < ln) {
      out["elements_truncated"] = true;
      break;
    }
    out["element_ids"].push_back(tag);
    try {
      if (tag == 0) {
        bool zero = true;
        for (size_t k = 0; k < v.size(); ++k) zero = zero && v.data()[k] == 0;
        if (ln == 0 || zero) {
          out["ssid"] = ln == 0 ? "<broadcast>" : "<hidden>";
          out["ssid_hidden"] = true;
        } else {
          out["ssid"] = detail::decode_utf8(v);
        }
      } else if (tag == 1 || tag == 50) {
        for (size_t k = 0; k < v.size(); ++k) {
          const int b = v.data()[k];
          if ((b & 0x7F) == 0x7F || (b & 0x7F) == 0x7E) continue;
          const double r = (b & 0x7F) / 2.0;
          rates.push_back(r);
          if (b & 0x80) basic.push_back(r);
        }
      } else if (tag == 3 && ln >= 1) {
        out["ds_channel"] = v[0];
      } else if (tag == 5 && ln >= 4) {
        const int bc = v[2];
        const Bytes pvb = v.from(3);
        const int offset = (bc >> 1) * 2;
        Json aids = Json::array();
        for (size_t k = 0; k < pvb.size(); ++k)
          for (int j = 0; j < 8; ++j)
            if ((pvb.data()[k] >> j) & 1) {
              const int a = 8 * (offset + static_cast<int>(k)) + j;
              if (a && aids.size() < 32) aids.push_back(a);
            }
        Json tim = Json::object();
        tim["dtim_count"] = v[0];
        tim["dtim_period"] = v[1];
        tim["multicast_buffered"] = static_cast<bool>(bc & 1);
        tim["aids_with_traffic"] = aids;
        out["tim"] = tim;
      } else if (tag == 7 && ln >= 3) {
        Json chans = Json::array();
        for (size_t k = 3; k + 2 < ln; k += 3) {
          if (v[k] >= 201) continue;
          Json c = Json::object();
          c["first"] = v[k];
          c["n"] = v[k + 1];
          c["max_dbm"] = v[k + 2] < 128 ? v[k + 2] : v[k + 2] - 256;
          chans.push_back(c);
        }
        Json c = Json::object();
        c["code"] = detail::decode_ascii(v.sub(0, 2));
        const IntMap env = {{0x20, "any"}, {0x4F, "outdoor"}, {0x49, "indoor"}, {0x58, "non-country"}};
        c["environment"] = get_or(env, v[2], fmt("0x%x", v[2]));
        c["channels"] = chans;
        out["country"] = c;
      } else if (tag == 11 && ln >= 5) {
        Json b = Json::object();
        b["stations"] = u16le(v, 0);
        b["channel_utilisation_pct"] = py_round(v[2] / 255.0 * 100, 1);
        b["available_admission_capacity"] = u16le(v, 3);
        out["bss_load"] = b;
      } else if (tag == 32 && ln >= 1) {
        out["power_constraint_db"] = v[0];
      } else if (tag == 42 && ln >= 1) {
        Json e = Json::object();
        e["non_erp_present"] = static_cast<bool>(v[0] & 1);
        e["use_protection"] = static_cast<bool>(v[0] & 2);
        e["barker_preamble_mode"] = static_cast<bool>(v[0] & 4);
        out["erp"] = e;
      } else if (tag == 45 && ln >= 26) {
        const int cap = u16le(v, 0);
        int ss = 0;
        for (size_t k = 3; k < 7; ++k) ss += v[k] != 0;
        Json h = Json::object();
        h["40mhz"] = static_cast<bool>(cap & 0x02);
        h["sgi20"] = static_cast<bool>(cap & 0x20);
        h["sgi40"] = static_cast<bool>(cap & 0x40);
        h["greenfield"] = static_cast<bool>(cap & 0x10);
        h["ldpc"] = static_cast<bool>(cap & 0x01);
        h["spatial_streams"] = ss;
        out["ht"] = h;
      } else if (tag == 61 && ln >= 2) {
        const IntMap secm = {{0, "none"}, {1, "above"}, {3, "below"}};
        const std::string sec = get_or(secm, v[1] & 3, "?");
        Json h = Json::object();
        h["primary_channel"] = v[0];
        h["secondary_channel"] = sec;
        h["width_mhz"] = ((v[1] & 4) && sec != "none") ? 40 : 20;
        out["ht_operation"] = h;
      } else if (tag == 48) {
        out["rsn"] = rsn(v);
      } else if (tag == 54 && ln >= 2) {
        out["mobility_domain"] = fmt("%04x", u16le(v, 0));
      } else if (tag == 70) {
        out["radio_measurement_802_11k"] = true;
      } else if (tag == 107) {
        out["interworking_802_11u"] = true;
      } else if (tag == 127) {
        Json ec = Json::array();
        const std::pair<int, const char*> bits[] = {
            {19, "BSS-transition (802.11v)"}, {31, "interworking"}, {62, "opmode-notification"},
            {70, "FTM-responder"}, {71, "FTM-initiator"}, {77, "TWT-requester"}, {78, "TWT-responder"}};
        for (const auto& [b, n] : bits)
          if (ext_bit(v, b)) ec.push_back(n);
        out["ext_capabilities"] = ec;
      } else if (tag == 191 && ln >= 12) {
        const uint32_t cap = static_cast<uint32_t>(detail::le_uint(v.sub(0, 4)));
        const int rx_map = u16le(v, 4);
        int ss = 0;
        for (int k = 0; k < 8; ++k) ss += ((rx_map >> (2 * k)) & 3) != 3;
        const IntMap ws = {{0, "80"}, {1, "160"}, {2, "160/80+80"}};
        Json h = Json::object();
        h["width_set"] = get_or(ws, (cap >> 2) & 3, "?");
        h["sgi80"] = static_cast<bool>(cap & 0x20);
        h["spatial_streams"] = ss;
        h["su_beamformer"] = static_cast<bool>(cap & (1u << 11));
        h["mu_beamformer"] = static_cast<bool>(cap & (1u << 19));
        out["vht"] = h;
      } else if (tag == 192 && ln >= 3) {
        const IntMap w = {{0, "20/40"}, {1, "80"}, {2, "160"}, {3, "80+80"}};
        Json h = Json::object();
        h["width"] = get_or_int(w, v[0]);
        h["center_seg0"] = v[1];
        h["center_seg1"] = v[2];
        out["vht_operation"] = h;
      } else if (tag == 221 && ln >= 4) {
        const std::string oui = mac_str(v.sub(0, 3));
        const int typ = v[3];
        std::string name = oui_name(oui);
        if (oui == "00:50:f2" && typ == 1) {
          out["wpa"] = rsn(v.from(4), "00:50:f2", true);
          name = "Microsoft WPA";
        } else if (oui == "00:50:f2" && typ == 2) {
          out["wmm"] = true;
          name = "Microsoft WMM/WME";
        } else if (oui == "00:50:f2" && typ == 4) {
          out["wps"] = wps(v.from(4));
          name = "Microsoft WPS";
        } else if (oui == "50:6f:9a") {
          const char* wt = lookup(kWfaTypes, typ);
          name = std::string("Wi-Fi Alliance ") + (wt ? std::string(wt) : std::to_string(typ));
          if (typ == 9) {
            out["p2p"] = p2p(v.from(4));
          } else if (typ == 0x1C && ln >= 11) {
            out["owe_transition_bssid"] = mac_str(v.sub(4, 10));
          } else if (typ == 0x10) {
            out["hotspot2"] = true;
          }
        } else {
          auto rid = droneid::from_vendor_ie(oui, typ, v.from(4).vec());
          if (rid) {
            name = (*rid)["standard"].get<std::string>();
            out["remote_id"] = std::move(*rid);
          }
        }
        vend.push_back(name);
      } else if (tag == 255 && ln >= 1) {
        ext_ids.push_back(v[0]);
        if (v[0] == 35) {
          out["he"] = true;
        } else if (v[0] == 36 && ln >= 7) {
          Json h = Json::object();
          h["bss_color"] = v[4] & 0x3F;
          out["he_operation"] = h;
        } else if (v[0] == 106 || v[0] == 108) {
          out["eht"] = true;
        }
      }
    } catch (const IndexError&) {
      if (!out.contains("malformed_elements")) out["malformed_elements"] = Json::array();
      out["malformed_elements"].push_back(tag);
    }
    i += 2 + ln;
  }
  if (!rates.empty()) {
    out["rates_mbps"] = sorted_unique(rates);
    out["basic_rates_mbps"] = sorted_unique(basic);
  }
  if (!vend.empty()) out["vendor_elements"] = vend;
  if (!ext_ids.empty()) out["ext_element_ids"] = ext_ids;
  return out;
}

std::optional<std::string> security(const Json& info) {
  const Json& r = jget(info, "rsn");
  const Json& w = jget(info, "wpa");
  bool priv = false;
  if (info.contains("capabilities"))
    for (const auto& c : info["capabilities"]) priv = priv || c == "privacy";
  if (!truthy(r) && !truthy(w)) {
    if (!info.contains("capabilities")) return std::nullopt;
    return std::string(priv ? "WEP" : "open");
  }
  std::vector<std::string> parts;
  if (truthy(r)) {
    std::set<std::string> akm;
    for (const auto& a : jget(r, "akm")) akm.insert(a.get<std::string>());
    std::string ci;
    for (const auto& c : jget(r, "pairwise_ciphers")) ci += (ci.empty() ? "" : "/") + c.get<std::string>();
    if (ci.empty()) ci = r.contains("group_cipher") ? r["group_cipher"].get<std::string>() : "?";
    auto any_of = [&](std::initializer_list<const char*> l) {
      for (const char* x : l)
        if (akm.count(x)) return true;
      return false;
    };
    std::string kind;
    if (akm.count("OWE")) {
      kind = "OWE (Enhanced Open)";
    } else if (any_of({"SAE", "FT-SAE", "SAE-EXT-KEY", "FT-SAE-EXT-KEY"})) {
      kind = any_of({"PSK", "FT-PSK", "PSK-SHA256"}) ? "WPA2/WPA3-Personal (transition)" : "WPA3-Personal (SAE)";
    } else if (any_of({"802.1X-SuiteB-192", "802.1X-SHA384"})) {
      kind = "WPA3-Enterprise 192-bit";
    } else if (std::any_of(akm.begin(), akm.end(), [](const std::string& a) {
                 return a.find("802.1X") != std::string::npos || a.find("FILS") != std::string::npos;
               })) {
      kind = "WPA2-Enterprise (802.1X)";
    } else if (any_of({"PSK", "FT-PSK", "PSK-SHA256"})) {
      kind = "WPA2-Personal (PSK)";
    } else {
      kind = "RSN ";
      bool first = true;
      for (const auto& a : akm) {
        kind += (first ? "" : ",") + a;
        first = false;
      }
    }
    std::string p = kind + " " + ci;
    if (std::any_of(akm.begin(), akm.end(), [](const std::string& a) { return a.rfind("FT", 0) == 0; }))
      p += " +FT";
    if (truthy(jget(r, "mfp_required"))) p += " MFP-required";
    parts.push_back(p);
  }
  if (truthy(w)) {
    std::string ci;
    bool first = true;
    for (const auto& c : jget(w, "pairwise_ciphers")) {
      ci += (first ? "" : "/") + c.get<std::string>();
      first = false;
    }
    parts.push_back("WPA1 " + ci);
  }
  std::string out;
  for (size_t i = 0; i < parts.size(); ++i) out += (i ? "; " : "") + parts[i];
  return out;
}

std::optional<std::string> phy_generation(const Json& info) {
  std::set<double> rates;
  for (const auto& r : jget(info, "rates_mbps")) rates.insert(r.get<double>());
  if (rates.empty() && !truthy(jget(info, "ht"))) return std::nullopt;
  std::vector<std::string> std_;
  for (double r : {1.0, 2.0, 5.5, 11.0})
    if (rates.count(r)) {
      std_.push_back("b");
      break;
    }
  for (double r : {6.0, 9.0, 12.0, 18.0, 24.0, 36.0, 48.0, 54.0})
    if (rates.count(r)) {
      std_.push_back("g");
      break;
    }
  for (const auto& [k, n] : {std::pair<const char*, const char*>{"ht", "n"}, {"vht", "ac"}, {"he", "ax"}, {"eht", "be"}})
    if (truthy(jget(info, k))) std_.push_back(n);
  if (std_.empty()) return std::nullopt;
  std::string out = "802.11";
  for (size_t i = 0; i < std_.size(); ++i) out += (i ? "/" : "") + std_[i];
  return out;
}

// ------------------------------------------------------------------------------- data bodies
namespace {

std::string ip4(const Bytes& b) {
  std::string s;
  for (size_t i = 0; i < b.size(); ++i) s += (i ? "." : "") + std::to_string(b.data()[i]);
  return s;
}

std::string ip6(const Bytes& b) {
  std::string s;
  for (size_t i = 0; i < 16; i += 2) s += (i ? ":" : "") + fmt("%x", be16(b, i));
  return s;
}

Json llc(const Bytes& body) {
  Json out = Json::object();
  if (body.size() < 8 || body[0] != 0xAA || body[1] != 0xAA || body[2] != 0x03) return out;
  const int et = (body[6] << 8) | body[7];
  out["ethertype"] = fmt("0x%04x", et);
  out["protocol"] = get_or(kEthertype, et, out["ethertype"].get<std::string>());
  const Bytes p = body.from(8);
  try {
    if (et == 0x888E && p.size() >= 4) {
      Json ev = Json::object();
      ev["eapol_version"] = p[0];
      ev["eapol_type"] = get_or_int({{0, "EAP"}, {1, "start"}, {2, "logoff"}, {3, "key"}}, p[1]);
      if (p[1] == 3 && p.size() >= 99) {
        const int ki = (p[5] << 8) | p[6];
        const int pair = ki & 8, inst = ki & 0x40, ack = ki & 0x80, mic = ki & 0x100, sec = ki & 0x200;
        if (pair) {
          int msg = (ack && !mic) ? 1 : (ack && mic && inst) ? 3 : (mic && sec && !ack) ? 4 : (mic && !ack) ? 2 : 0;
          ev["handshake"] = msg ? "4-way handshake message " + std::to_string(msg) : std::string("pairwise key");
        } else {
          ev["handshake"] = "group key handshake";
        }
        ev["key_descriptor_version"] = ki & 7;
        ev["replay_counter"] = detail::be_uint(p.sub(9, 17));
      } else if (p[1] == 0 && p.size() >= 9) {
        ev["eap"] = get_or_int({{1, "request"}, {2, "response"}, {3, "success"}, {4, "failure"}}, p[4]);
        if ((p[4] == 1 || p[4] == 2) && p.size() >= 9) {
          ev["eap_method"] =
              get_or_int({{1, "identity"}, {13, "EAP-TLS"}, {21, "EAP-TTLS"}, {25, "PEAP"}, {43, "EAP-FAST"},
                          {52, "EAP-PWD"}},
                         p[8]);
          if (p[8] == 1 && p[4] == 2) {
            const size_t end = 4 + static_cast<size_t>((p[6] << 8) | p[7]);
            ev["identity"] = detail::decode_utf8(p.sub(9, end));
          }
        }
      }
      out["eapol"] = ev;
    } else if (et == 0x0806 && p.size() >= 28) {
      Json a = Json::object();
      a["op"] = get_or_int({{1, "request"}, {2, "reply"}}, be16(p, 6));
      a["sender_mac"] = mac_str(p.sub(8, 14));
      a["sender_ip"] = ip4(p.sub(14, 18));
      a["target_mac"] = mac_str(p.sub(18, 24));
      a["target_ip"] = ip4(p.sub(24, 28));
      out["arp"] = a;
    } else if (et == 0x0800 && p.size() >= 20) {
      const size_t ihl = static_cast<size_t>(p[0] & 15) * 4;
      const int proto = p[9];
      Json ip = Json::object();
      ip["src"] = ip4(p.sub(12, 16));
      ip["dst"] = ip4(p.sub(16, 20));
      ip["proto"] = get_or_int(kIpProto, proto);
      ip["ttl"] = p[8];
      if ((proto == 6 || proto == 17) && p.size() >= ihl + 4) {
        const int sp = be16(p, ihl), dp = be16(p, ihl + 2);
        ip["sport"] = sp;
        ip["dport"] = dp;
        const char* svc = lookup(kPorts, dp);
        if (!svc) svc = lookup(kPorts, sp);
        if (svc) ip["service"] = svc;
        // MAVLink over an open Wi-Fi link (research, PLAN G2)
        size_t off = proto == 17 ? ihl + 8 : (p.size() > ihl + 12 ? ihl + 4 * static_cast<size_t>(p[ihl + 12] >> 4) : p.size());
        Json mv = mavlink_frames(p.from(off));
        if (!mv.is_null()) ip["mavlink"] = mv;
        if (proto == 17) {  // PLAN W4: Tello SDK, video streams (open links)
          const Json app = udp_app(p.from(off), sp, dp);
          for (auto it = app.begin(); it != app.end(); ++it) ip[it.key()] = it.value();
        }
      }
      out["ipv4"] = ip;
    } else if (et == 0x86DD && p.size() >= 40) {
      const int nh = p[6];
      Json ip = Json::object();
      ip["src"] = ip6(p.sub(8, 24));
      ip["dst"] = ip6(p.sub(24, 40));
      ip["next_header"] = get_or_int(kIpProto, nh);
      if ((nh == 6 || nh == 17) && p.size() >= 44) {
        ip["sport"] = be16(p, 40);
        ip["dport"] = be16(p, 42);
      }
      out["ipv6"] = ip;
    }
  } catch (const IndexError&) {
    out["truncated"] = true;
  }
  return out;
}

Json crypto_header(const Bytes& body) {
  Json out = Json::object();
  if (body.size() < 4) return out;
  const bool ext_iv = body[3] & 0x20;
  const int key_id = body[3] >> 6;
  if (!ext_iv) {
    out["cipher"] = "WEP";
    out["key_id"] = key_id;
    out["iv"] = detail::hex(body.sub(0, 3));
    return out;
  }
  if (body.size() < 8) {
    out["key_id"] = key_id;
    return out;
  }
  const uint64_t hi = detail::le_uint(body.sub(4, 8)) << 16;
  if (body[1] == ((body[0] | 0x20) & 0x7F)) {  // TKIP WEPSeed[1]
    out["cipher"] = "TKIP";
    out["key_id"] = key_id;
    out["pn"] = static_cast<uint64_t>(body[2]) | (static_cast<uint64_t>(body[0]) << 8) | hi;
    return out;
  }
  out["cipher"] = "CCMP/GCMP";
  out["key_id"] = key_id;
  out["pn"] = static_cast<uint64_t>(body[0]) | (static_cast<uint64_t>(body[1]) << 8) | hi;
  return out;
}

void vendors(Json& info) {
  for (const char* k : {"addr1", "addr2", "addr3", "addr4", "addr1_ra", "addr2_ta", "sa", "da", "bssid", "ra", "ta"}) {
    if (truthy(jget(info, k))) {
      auto v = vendor(info[k].get<std::string>());
      info[std::string(k) + "_vendor"] = v ? Json(*v) : Json(nullptr);
    }
  }
}

std::optional<std::string> opt_str(const Json& v) {
  if (v.is_string()) return v.get<std::string>();
  return std::nullopt;
}

void drone(Json& info) {
  const Json& rid = jget(info, "remote_id");
  if (truthy(rid)) {
    Json d = Json::object();
    if (truthy(jget(rid, "model"))) d["platform"] = rid["model"];
    else if (truthy(jget(rid, "uas_id"))) d["platform"] = rid["uas_id"];
    else d["platform"] = "UAS";
    d["source"] = "Remote ID broadcast";
    d["confidence"] = 0.99;
    d["evidence"] = droneid::summary(rid);
    info["drone"] = d;
    return;
  }
  {  // wfb-ng link (research, PLAN G1)
    const Json& w = jget(info, "wfb_ng");
    if (truthy(w)) {
      Json d = Json::object();
      d["platform"] = "wfb-ng digital FPV / telemetry link (OpenIPC / OpenHD)";
      d["source"] = "wfb-ng 802.11 injection signature";
      d["confidence"] = 0.9;
      d["evidence"] = "link_id " + w["link_id_hex"].get<std::string>() +
                      (w["default_link_domain"].get<bool>() ? " (default link domain)" : "") + ", radio port " +
                      pstr(w["radio_port"]) + " (" + w["stream"].get<std::string>() + ")";
      info["drone"] = d;
      return;
    }
  }
  {  // ExpressLRS Backpack (research, PLAN G7)
    const Json& bp = jget(info, "elrs_backpack");
    if (truthy(bp)) {
      Json d = Json::object();
      d["platform"] = "ExpressLRS Backpack (bind UID " + bp["bind_uid"].get<std::string>() + ")";
      d["source"] = "ESP-NOW MSP frame (CRC-valid)";
      d["confidence"] = 0.85;
      d["evidence"] = "MSP " + bp["msp"].get<std::string>() +
                      (bp.contains("vtx_band") ? ", VTX " + bp["vtx_band"].get<std::string>() + pstr(bp["vtx_channel"]) : "");
      info["drone"] = d;
      return;
    }
  }
  {  // a CRC-valid MAVLink vehicle over Wi-Fi (research, PLAN G2)
    const Json& mv = jget(jget(info, "ipv4"), "mavlink");
    if (truthy(mv) && mv.contains("vehicle_type") && !mv["vehicle_type"].is_null()) {
      Json d = Json::object();
      d["platform"] = "MAVLink vehicle (sysid " + pstr(jget(mv, "sysid")) + ", " + pstr(mv["vehicle_type"]) + ")";
      d["source"] = "MAVLink telemetry over Wi-Fi (CRC-valid)";
      d["confidence"] = 0.9;
      std::string msgs;
      for (const auto& m : jget(mv, "messages")) msgs += (msgs.empty() ? "" : ", ") + m.get<std::string>();
      d["evidence"] = pstr(mv["frames"]) + " MAVLink frame(s): " + msgs;
      info["drone"] = d;
      return;
    }
  }
  {  // Tello SDK over open Wi-Fi (research, PLAN W4)
    const Json& te = jget(jget(info, "ipv4"), "tello");
    if (truthy(te)) {
      const bool st = te["port"] == "state";
      Json d = Json::object();
      d["platform"] = "DJI / Ryze Tello (Tello SDK)";
      d["source"] = "Tello SDK over open Wi-Fi";
      d["confidence"] = st ? 0.9 : 0.8;
      d["evidence"] = st ? "state: battery " + pstr(jget(te, "battery_pct")) + " %, height " + pstr(jget(te, "height_cm")) + " cm"
                         : "command '" + pstr(te["text"]) + "' (" + pstr(te["direction"]) + ")";
      info["drone"] = d;
      return;
    }
  }
  const Json* tx = &jget(info, "ta");
  if (!truthy(*tx)) tx = &jget(info, "sa");
  if (!truthy(*tx)) tx = &jget(info, "bssid");
  std::optional<std::string> vend;
  if (truthy(*tx)) vend = vendor(tx->get<std::string>());
  auto sig = droneid::signature(vend, opt_str(jget(info, "ssid")));
  if (sig) {
    (*sig)["source"] = "platform signature (no Remote ID in this frame)";
    (*sig)["platform"] = (*sig)["platform"].get<std::string>() + "-like";
    info["drone"] = *sig;
  }
}

const std::set<std::string>& suffixes() {
  static const std::set<std::string> s = {
      "co", "co.", "ltd", "ltd.", "inc", "inc.", "corp", "corp.", "corporation", "pvt", "private", "limited",
      "llc", "gmbh", "company", "technology", "technologies", "systems", "communications", "electronics",
      "digital", "industrial", "international", "(india)", "india", "network", "networks", "co.,ltd",
      "co.,ltd.", "technologies,"};
  return s;
}

std::string short_name(const Json& vj) {
  if (!truthy(vj)) return "?";
  const std::string v = vj.get<std::string>();
  if (v.rfind("locally-administered", 0) == 0) return "randomized";
  auto words = detail::py_split_ws(v.substr(0, v.find(',')));
  while (words.size() > 1 && suffixes().count(detail::ascii_lower(words.back()))) words.pop_back();
  std::string out;
  for (size_t i = 0; i < words.size() && i < 3; ++i) out += (i ? " " : "") + words[i];
  return out;
}

// a or b (Python truthiness: the first truthy, else the last)
const Json& or2(const Json& a, const Json& b) { return truthy(a) ? a : b; }

}  // namespace

std::string summary(const Json& info) {
  const std::string st = info.contains("subtype") ? pstr(info["subtype"]) : "?";
  std::vector<std::string> s = {st};
  if (truthy(jget(info, "retry"))) s.push_back("(retry)");
  if (truthy(jget(info, "remote_id"))) return st + " " + droneid::summary(info["remote_id"]);
  const std::string tx = pstr(or2(jget(info, "sa"), jget(info, "ta")));
  const std::string rx = pstr(or2(jget(info, "da"), jget(info, "ra")));
  const Json& type = jget(info, "type");
  auto get_d = [&](const char* k, const char* dflt) {
    return info.contains(k) ? pstr(info[k]) : std::string(dflt);
  };
  if (type == "mgmt" && (st == "beacon" || st == "probe-resp")) {
    s.push_back("BSSID " + pstr(jget(info, "bssid")) + " [" + short_name(jget(info, "bssid_vendor")) + "]");
    s.push_back("SSID '" + get_d("ssid", "?") + "'");
    if (truthy(jget(info, "ds_channel"))) s.push_back("ch" + pstr(info["ds_channel"]));
    for (const char* k : {"security", "phy"})
      if (truthy(jget(info, k))) s.push_back(pstr(info[k]));
    if (truthy(jget(info, "ap_uptime"))) s.push_back("up " + pstr(info["ap_uptime"]));
    if (truthy(jget(info, "bss_load")))
      s.push_back(pstr(info["bss_load"]["stations"]) + " STA " + pstr(info["bss_load"]["channel_utilisation_pct"]) +
                  "% busy");
    const Json& w = jget(info, "wps");
    std::string dev;
    for (const char* k : {"manufacturer", "model_name", "device_name"})
      if (truthy(jget(w, k))) dev += (dev.empty() ? "" : " ") + pstr(w[k]);
    if (!dev.empty()) s.push_back("WPS '" + dev + "'");
  } else if (st == "probe-req") {
    const Json& v = or2(or2(jget(info, "sa_vendor"), jget(info, "ta_vendor")), jget(info, "addr2_vendor"));
    s.push_back("from " + tx + " [" + short_name(v) + "]");
    s.push_back("for '" + get_d("ssid", "?") + "'");
    if (truthy(jget(info, "phy"))) s.push_back(pstr(info["phy"]));
    const Json& name = or2(jget(jget(info, "p2p"), "device_name"), jget(jget(info, "wps"), "device_name"));
    if (truthy(name)) s.push_back("device '" + pstr(name) + "'");
  } else if (st == "deauth" || st == "disassoc") {
    s.push_back(tx + " -> " + rx + ": " + pstr(jget(info, "reason_text")));
  } else if (st == "auth") {
    s.push_back(tx + " -> " + rx + ": " + pstr(jget(info, "auth_algorithm")) + " seq " + pstr(jget(info, "auth_seq")) +
                " " + pstr(jget(info, "status_text")));
  } else if (st.rfind("assoc", 0) == 0 || st.rfind("reassoc", 0) == 0) {
    s.push_back(tx + " -> " + rx);
    if (truthy(jget(info, "ssid"))) s.push_back("SSID '" + pstr(info["ssid"]) + "'");
    if (info.contains("status_text")) s.push_back(pstr(info["status_text"]) + " AID " + pstr(jget(info, "aid")));
  } else if (type == "data") {
    s.push_back(get_d("direction", "") + " " + tx + " -> " + rx + " BSSID " + pstr(jget(info, "bssid")));
    const Json& enc = jget(info, "encryption");
    if (truthy(enc)) {
      const std::string cipher = enc.contains("cipher") ? pstr(enc["cipher"]) : "encrypted";
      if (!jget(enc, "pn").is_null()) {
        s.push_back(cipher + " PN " + pstr(enc["pn"]));
      } else {
        s.push_back(cipher);
      }
    }
    if (truthy(jget(info, "eapol"))) {
      const Json& e = info["eapol"];
      const std::string what =
          truthy(jget(e, "handshake")) ? pstr(e["handshake"]) : (e.contains("eap") ? pstr(e["eap"]) : "");
      s.push_back("EAPOL " + what);
    } else if (truthy(jget(info, "protocol"))) {
      s.push_back(pstr(info["protocol"]));
    }
    for (const char* k : {"ipv4", "ipv6"}) {
      if (truthy(jget(info, k))) {
        const Json& ip = info[k];
        std::string line = pstr(ip["src"]) + " -> " + pstr(ip["dst"]) + " " +
                           pstr(or2(jget(ip, "proto"), jget(ip, "next_header")));
        if (truthy(jget(ip, "service"))) line += " " + pstr(ip["service"]);
        if (truthy(jget(ip, "video")))
          line += " video " + pstr(ip["video"]["codec"]) + " " + pstr(ip["video"]["nal"]) + " (" +
                  pstr(ip["video"]["transport"]) + ")";
        s.push_back(line);
      }
    }
    if (truthy(jget(info, "arp"))) {
      const Json& a = info["arp"];
      s.push_back("ARP " + pstr(a["op"]) + " " + pstr(a["sender_ip"]) + " -> " + pstr(a["target_ip"]));
    }
  } else if (type == "ctrl") {
    std::string line;
    for (const char* k : {"ta", "ra"})
      if (truthy(jget(info, k))) line += (line.empty() ? "" : " -> ") + pstr(info[k]);
    s.push_back(line);
  } else if (st == "action" || st == "action-no-ack") {
    if (truthy(jget(info, "action")))
      s.push_back(tx + " -> " + rx + ": " + pstr(info["action"]["category_name"]) + " #" +
                  pstr(info["action"]["code"]));
    else
      s.push_back(tx + " -> " + rx);
  }
  if (truthy(jget(info, "drone"))) s.push_back("[" + pstr(info["drone"]["platform"]) + "]");
  std::string out;
  bool first = true;
  for (const auto& x : s) {
    if (x.empty()) continue;
    out += (first ? "" : " ") + x;
    first = false;
  }
  return out;
}

// ------------------------------------------------------------------------------- frame
Json parse(const std::vector<uint8_t>& mpdu) { return parse(mpdu.data(), mpdu.size()); }

Json parse(const uint8_t* p, size_t n) {
  const Bytes m(p, n);
  Json info = Json::object();
  if (n < 2) {
    info["parse"] = "too short";
    return info;
  }
  const int fc0 = m[0], fc1 = m[1];
  const int ftype = (fc0 >> 2) & 3, sub = (fc0 >> 4) & 0xF;
  const int to_ds = fc1 & 1, from_ds = (fc1 >> 1) & 1;
  info["type"] = get_or(kFtype, ftype, "?");
  auto sit = kSubtype.find({ftype, sub});
  info["subtype"] = sit != kSubtype.end() ? std::string(sit->second) : fmt("%d/%d", ftype, sub);
  info["protocol_version"] = fc0 & 3;
  info["to_ds"] = to_ds;
  info["from_ds"] = from_ds;
  info["more_frag"] = (fc1 >> 2) & 1;
  info["retry"] = (fc1 >> 3) & 1;
  info["pwr_mgmt"] = (fc1 >> 4) & 1;
  info["more_data"] = (fc1 >> 5) & 1;
  info["protected"] = (fc1 >> 6) & 1;
  info["order"] = (fc1 >> 7) & 1;
  info["length_bytes"] = n + 4;
  const std::string st = info["subtype"].get<std::string>();
  if (n >= 4) {
    const int d = u16le(m, 2);
    if (st == "ps-poll") info["aid"] = d & 0x3FFF;
    else info["duration_us"] = d < 0x8000 ? Json(d) : Json(nullptr);
  }
  if (ftype == 1) {  // control
    if (n >= 10) {
      info["addr1_ra"] = mac_str(m.sub(4, 10));
      info["ra"] = info["addr1_ra"];
    }
    static const std::set<std::string> with_ta = {"rts", "block-ack-req", "block-ack", "ps-poll", "bf-report-poll",
                                                  "ndp-announcement", "cf-end", "trigger"};
    if (with_ta.count(st) && n >= 16) {
      info["addr2_ta"] = mac_str(m.sub(10, 16));
      info["ta"] = info["addr2_ta"];
    }
    if (st == "ps-poll") info["bssid"] = jget(info, "ra");
    if ((st == "block-ack-req" || st == "block-ack") && n >= 20) {
      const int ctl = u16le(m, 16);
      Json ba = Json::object();
      ba["tid"] = ctl >> 12;
      ba["start_seq"] = u16le(m, 18) >> 4;
      ba["bitmap"] = st == "block-ack" ? Json(detail::hex(m.sub(20, 28))) : Json(nullptr);
      info["block_ack"] = ba;
    }
    vendors(info);
    info["summary"] = summary(info);
    return info;
  }
  if (n >= 10) info["addr1"] = mac_str(m.sub(4, 10));
  if (n >= 16) info["addr2"] = mac_str(m.sub(10, 16));
  if (n >= 22) info["addr3"] = mac_str(m.sub(16, 22));
  if (n >= 24) {
    const int sc = u16le(m, 22);
    info["seq"] = sc >> 4;
    info["frag"] = sc & 0xF;
  }
  size_t hdr = 24;
  if (ftype == 2 && to_ds && from_ds && n >= 30) {
    info["addr4"] = mac_str(m.sub(24, 30));
    hdr = 30;
  }
  const Json a1 = jget(info, "addr1"), a2 = jget(info, "addr2"), a3 = jget(info, "addr3"), a4 = jget(info, "addr4");
  info["ra"] = a1;
  info["ta"] = a2;
  if (ftype == 2) {
    if (!to_ds && !from_ds) {
      info["da"] = a1, info["sa"] = a2, info["bssid"] = a3, info["direction"] = "IBSS/direct";
    } else if (!to_ds && from_ds) {
      info["da"] = a1, info["bssid"] = a2, info["sa"] = a3, info["direction"] = "AP -> STA";
    } else if (to_ds && !from_ds) {
      info["bssid"] = a1, info["sa"] = a2, info["da"] = a3, info["direction"] = "STA -> AP";
    } else {
      info["da"] = a3, info["sa"] = a4, info["direction"] = "WDS / mesh (4-address)";
    }
    if ((sub & 8) && n >= hdr + 2) {  // QoS control
      const int q = u16le(m, hdr);
      static const char* kAck[] = {"normal", "no-ack", "no-explicit", "block-ack"};
      Json qos = Json::object();
      qos["tid"] = q & 15;
      qos["ack_policy"] = kAck[(q >> 5) & 3];
      qos["amsdu"] = static_cast<bool>(q & 0x80);
      info["qos"] = qos;
      hdr += 2;
    }
    if (truthy(info["order"]) && n >= hdr + 4) hdr += 4;  // HT control
    const Bytes body = m.from(hdr);
    info["body_bytes"] = body.size();
    {   // wfb-ng (research, PLAN G1)
      const Json& a2 = jget(info, "addr2");
      if (sub == 0 && jget(info, "addr1") == "ff:ff:ff:ff:ff:ff" && a2.is_string() &&
          a2.get<std::string>().rfind("57:42:", 0) == 0 && jget(info, "addr3") == a2)
        info["wfb_ng"] = wfb(m.sub(10, 16), body);
    }
    if (truthy(info["protected"])) {
      info["encryption"] = crypto_header(body);
    } else if (!body.empty() && !(sub & 4)) {  // not a null frame
      if (truthy(jget(jget(info, "qos"), "amsdu"))) {
        info["amsdu"] = true;
      } else {
        const Json l = llc(body);
        for (auto it = l.begin(); it != l.end(); ++it) info[it.key()] = it.value();
      }
    }
  } else if (ftype == 0) {
    info["bssid"] = jget(info, "addr3");
    if (truthy(info["protected"])) {
      info["encryption"] = crypto_header(m.from(24));
      vendors(info);
      info["summary"] = summary(info);
      return info;
    }
    const Bytes b = m.from(24);
    std::optional<Bytes> ie;
    if ((st == "beacon" || st == "probe-resp") && b.size() >= 12) {
      const uint64_t tsf = detail::le_uint(b.sub(0, 8));
      info["tsf_us"] = tsf;
      info["ap_uptime"] = uptime(tsf);
      info["beacon_interval_tu"] = u16le(b, 8);
      info["beacon_interval_ms"] = py_round(u16le(b, 8) * 1.024, 1);
      info["capabilities"] = caps(u16le(b, 10));
      ie = b.from(12);
    } else if (st == "probe-req") {
      ie = b;
    } else if ((st == "assoc-req" || st == "reassoc-req") && b.size() >= 4) {
      info["capabilities"] = caps(u16le(b, 0));
      info["listen_interval"] = u16le(b, 2);
      if (st == "reassoc-req" && b.size() >= 10) {
        info["current_ap"] = mac_str(b.sub(4, 10));
        ie = b.from(10);
      } else {
        ie = b.from(4);
      }
    } else if ((st == "assoc-resp" || st == "reassoc-resp") && b.size() >= 6) {
      info["capabilities"] = caps(u16le(b, 0));
      info["status"] = u16le(b, 2);
      info["status_text"] = get_or(kStatus, u16le(b, 2), "?");
      info["aid"] = u16le(b, 4) & 0x3FFF;
      ie = b.from(6);
    } else if (st == "auth" && b.size() >= 6) {
      const int alg = u16le(b, 0);
      info["auth_algorithm"] = get_or_int(kAuthAlg, alg);
      info["auth_seq"] = u16le(b, 2);
      info["status"] = u16le(b, 4);
      info["status_text"] = get_or(kStatus, u16le(b, 4), "?");
      if (alg == 3 && b.size() >= 8 && u16le(b, 2) == 1) info["sae_group"] = u16le(b, 6);
      if (alg != 3) ie = b.from(6);
    } else if ((st == "deauth" || st == "disassoc") && b.size() >= 2) {
      info["reason"] = u16le(b, 0);
      info["reason_text"] = get_or(kReason, u16le(b, 0), "?");
    } else if ((st == "action" || st == "action-no-ack") && b.size() >= 2) {
      Json a = Json::object();
      a["category"] = b[0];
      a["category_name"] = get_or_int(kActionCat, b[0] & 0x7F);
      if (!lookup(kActionCat, b[0] & 0x7F)) a["category_name"] = b[0];
      a["code"] = b[1];
      const bool wfa = b[0] == 4 && b[1] == 9 && b.size() >= 6 && mac_str(b.sub(2, 5)) == "50:6f:9a";
      if (wfa) a["detail"] = "Wi-Fi Alliance public action (P2P/DPP/NAN)";
      info["action"] = a;
      if (b[0] == 127 && b.size() >= 8 && mac_str(b.sub(1, 4)) == "18:fe:34") {   // ESP-NOW (research, PLAN G7)
        Json en = espnow(b, info);
        info["esp_now"] = en;
      }
      if (wfa)
        if (auto rid = droneid::from_nan_action(b.vec())) info["remote_id"] = *rid;  // Remote ID over Wi-Fi NAN (research)
    }
    if (ie) {
      const Json el = parse_elements(ie->data(), ie->size());
      for (auto it = el.begin(); it != el.end(); ++it) info[it.key()] = it.value();
      if (st != "probe-req" && jget(info, "ssid") == "<broadcast>") info["ssid"] = "<hidden>";
      auto sec = security(info);
      if (sec && !sec->empty()) info["security"] = *sec;
      auto phy = phy_generation(info);
      if (phy && !phy->empty()) info["phy"] = *phy;
    }
  }
  vendors(info);
  drone(info);
  info["summary"] = summary(info);
  return info;
}

}  // namespace si::rx::wifi::mac80211
