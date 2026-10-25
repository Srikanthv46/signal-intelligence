// python receivers/rc_link/hubsan.py (see there for the PHY, the packets and the model constants).
#include "receivers/rc_link/hubsan.hpp"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "receivers/rc_link/sik.hpp"

namespace si::rx::rc_link::hubsan {
namespace {
const char* header_name(unsigned h) {
  switch (h) {
    case 0x20: return "data";
    case 0x40: return "vtx";
    case 0xE0: case 0xE1: return "telemetry";
    case 0xE8: case 0xE9: return "keep-alive";
    case 0x01: return "bind 1";
    case 0x03: return "bind 3";
    case 0x09: return "bind 9";
    default: return nullptr;
  }
}
const char* bind_id(unsigned idc) {
  return idc == 0x55201041u ? "bind (normal)" : idc == 0xAA201041u ? "bind (plus)" : nullptr;
}
std::string hexs(const uint8_t* p, size_t n) {
  std::string s;
  char b[3];
  for (size_t i = 0; i < n; ++i) std::snprintf(b, sizeof b, "%02x", p[i]), s += b;
  return s;
}
Json model_of(const uint8_t* p) {
  Json m = Json::array();
  if (p[10] == 0x64 && p[11] == 0xDB && p[12] == 0x04 && p[13] == 0x26 && p[14] == 0x79)
    m.push_back("H107 / H107L / H107C / H107D");
  else if (p[10] == 0x18 && p[12] == 0x5C && p[14] == 0xF6)
    m.push_back("H301");
  else if (p[10] == 0x19 && p[12] == 0x5C && p[14] == 0x49)
    m.push_back("H107 Plus");
  else if ((p[11] == 0x41 || p[11] == 0x42 || p[11] == 0x43) && p[3] == 0 && p[5] == 0 && p[7] == 0)
    m.push_back("H501 / H122D / H502 class");
  return m;
}
}  // namespace

Json try_decode(const cfvec& y, double ofs, std::optional<double> f_hz, double rate) {
  for (const auto& [ph, bits] : sik::demod_bits_mid(y, ofs, rate))
    for (bool inv : {false, true}) {
      std::vector<uint8_t> b(bits.size());
      for (size_t i = 0; i < b.size(); ++i) b[i] = inv ? 1 - bits[i] : bits[i];
      // python _alternating_runs(b, 16)
      std::vector<int64_t> edges = {-1};
      for (size_t i = 0; i + 1 < b.size(); ++i)
        if (b[i] == b[i + 1]) edges.push_back(static_cast<int64_t>(i));
      edges.push_back(static_cast<int64_t>(b.size()) - 1);
      for (size_t e = 0; e + 1 < edges.size(); ++e) {
        if (edges[e + 1] - edges[e] < 16) continue;
        const int64_t r0 = edges[e] + 1, r1 = edges[e + 1] + 1;
        for (int64_t s = r0 + 8; s <= r1; ++s) {
          if (static_cast<size_t>(s) + 32 + 8 * kLen > b.size()) break;
          unsigned idc = 0;
          for (int i = 0; i < 32; ++i) idc = (idc << 1) | b[static_cast<size_t>(s + i)];
          uint8_t p[kLen];
          unsigned sum = 0;
          for (int k = 0; k < kLen; ++k) {
            unsigned v = 0;
            for (int i = 0; i < 8; ++i) v = (v << 1) | b[static_cast<size_t>(s + 32 + 8 * k + i)];
            p[k] = static_cast<uint8_t>(v);
            sum += v;
          }
          const char* kind = header_name(p[0]);
          if (!kind || (sum & 0xFF)) continue;
          const Json models = p[0] == 0x20 ? model_of(p) : Json::array();
          if (p[0] == 0x20 && models.empty()) continue;
          const std::string k = kind;
          if (k.rfind("bind", 0) == 0 && !bind_id(idc)) continue;
          char id[16], t[8];
          std::snprintf(id, sizeof id, "0x%08X", idc);
          std::snprintf(t, sizeof t, "0x%02X", p[0]);
          Json out = Json::object();
          out["protocol"] = "Hubsan (A7105)";
          out["id_code"] = id;
          out["id_kind"] = bind_id(idc) ? bind_id(idc) : "session";
          out["packet_type"] = t;
          out["packet_name"] = k;
          out["checksum_ok"] = true;
          out["polarity_inverted"] = inv;
          out["packet_hex"] = hexs(p, kLen);
          out["models"] = models;
          if (p[0] == 0x20)
            out["sticks"] = {{"throttle", p[2]}, {"rudder", 255 - p[4]}, {"elevator", 255 - p[6]}, {"aileron", p[8]}};
          if (k.rfind("bind", 0) == 0) out["tx_id"] = hexs(p + 2, 4);
          if (k == "telemetry") out["battery_v"] = py_round(p[13] / 10.0, 1);
          if (f_hz) {
            const int64_t n = static_cast<int64_t>(std::nearbyint((*f_hz - kF0) / kStep));
            out["channel"] = n;
            bool ok = false;
            for (int64_t c = 0x14; c <= 0x82; c += 10) ok = ok || n == c || n - 0x23 == c;
            out["channel_check"] = ok ? "match" : "mismatch";
          }
          return out;
        }
      }
    }
  return nullptr;
}

}  // namespace si::rx::rc_link::hubsan
