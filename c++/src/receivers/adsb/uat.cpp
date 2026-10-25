// python receivers/adsb/uat.py (see there for the sources, the PHY, the FEC and the fields).
#include "receivers/adsb/uat.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#include "core/npcompat.hpp"
#include "core/rs.hpp"

namespace si::rx::adsb::uat {
namespace {
constexpr uint64_t kSync = 0xEACDDA4E2ull;
const char* kQualifier[] = {"ADS-B ICAO", "ADS-B self-assigned", "TIS-B ICAO", "TIS-B track file", "surface vehicle",
                            "fixed beacon", "ADS-R self-assigned", "reserved"};
const char kB40[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ *??";

const rs::RS& rs_short() {
  static const rs::RS r(0x187, 120, 1, 12);
  return r;
}
const rs::RS& rs_long() {
  static const rs::RS r(0x187, 120, 1, 14);
  return r;
}

const char* category(int c) {
  switch (c) {
    case 0: return "no information"; case 1: return "light"; case 2: return "small"; case 3: return "large";
    case 4: return "high-vortex large"; case 5: return "heavy"; case 6: return "highly manoeuvrable";
    case 7: return "rotorcraft"; case 9: return "glider / sailplane"; case 10: return "lighter-than-air";
    case 11: return "parachutist / skydiver"; case 12: return "ultralight / hang-glider / paraglider";
    case 14: return "unmanned aerial vehicle"; case 15: return "space / trans-atmospheric";
    case 17: return "surface emergency vehicle"; case 18: return "surface service vehicle";
    case 19: return "point obstacle"; case 20: return "cluster obstacle"; case 21: return "line obstacle";
    default: return nullptr;
  }
}

// DO-282 bit range (byte, bit), 1-based, bit 1 = MSB, inclusive
int64_t bits(const std::vector<uint8_t>& m, int b1, int i1, int b2, int i2) {
  const int lo = (b1 - 1) * 8 + (i1 - 1), hi = (b2 - 1) * 8 + (i2 - 1);
  int64_t v = 0;
  for (int k = lo; k <= hi; ++k) v = (v << 1) | ((m[static_cast<size_t>(k / 8)] >> (7 - k % 8)) & 1);
  return v;
}
}  // namespace

Json parse(const std::vector<uint8_t>& m) {
  const int pt = m[0] >> 3, aq = m[0] & 7;
  char addr[8];
  std::snprintf(addr, sizeof addr, "%02X%02X%02X", m[1], m[2], m[3]);
  Json out = Json::object();
  out["payload_type"] = pt;
  out["address_qualifier"] = kQualifier[aq];
  out["address"] = addr;
  out["frame"] = m.size() == 34 ? "long" : "basic";
  const int64_t lat_r = bits(m, 5, 1, 7, 7), lon_r = bits(m, 7, 8, 10, 7);
  if (lat_r || lon_r) {
    double lat = static_cast<double>(lat_r) * 360 / (1 << 24), lon = static_cast<double>(lon_r) * 360 / (1 << 24);
    out["lat"] = py_round(lat > 90 ? lat - 180 : lat, 5);
    out["lon"] = py_round(lon > 180 ? lon - 360 : lon, 5);
  }
  const int64_t alt = bits(m, 11, 1, 12, 4);
  out["alt_ft"] = alt == 0 ? Json(nullptr) : Json((alt - 41) * 25);
  out["alt_type"] = bits(m, 10, 8, 10, 8) ? "geometric" : "pressure";
  out["nic"] = bits(m, 12, 5, 12, 8);
  const int64_t ag = bits(m, 13, 1, 13, 2);
  const char* ags[] = {"airborne subsonic", "airborne supersonic", "on ground", "reserved"};
  out["air_ground"] = ags[ag];
  if (ag == 0 || ag == 1) {
    const int64_t k = ag == 1 ? 4 : 1;
    const int64_t ns = bits(m, 13, 4, 13, 4), nv = bits(m, 13, 5, 14, 6);
    const int64_t es = bits(m, 14, 7, 14, 7), ev = bits(m, 14, 8, 16, 1);
    const int64_t vs = bits(m, 16, 3, 16, 3), vv = bits(m, 16, 4, 17, 4);
    if (nv) out["north_kt"] = (nv - 1) * k * (ns ? -1 : 1);
    if (ev) out["east_kt"] = (ev - 1) * k * (es ? -1 : 1);
    if (vv) out["vert_fpm"] = (vv - 1) * 64 * (vs ? -1 : 1);
  }
  if (m.size() == 34 && (pt == 1 || pt == 3)) {
    int w[3];
    for (int i = 0; i < 3; ++i) w[i] = (m[static_cast<size_t>(17 + 2 * i)] << 8) | m[static_cast<size_t>(18 + 2 * i)];
    const int cat = (w[0] / 1600) % 40;
    std::string cs;
    for (int x : w) cs += std::string(1, kB40[(x / 1600) % 40]) + kB40[(x / 40) % 40] + kB40[x % 40];
    cs = cs.substr(1);
    while (!cs.empty() && (cs.back() == ' ' || cs.back() == '*')) cs.pop_back();
    const char* cn = category(cat);
    out["emitter_category"] = cn ? std::string(cn) : "reserved " + std::to_string(cat);
    out["emitter_category_code"] = cat;
    out["uav"] = cat == 14;
    out["callsign"] = cs;
    out["emergency"] = bits(m, 24, 1, 24, 3);
    out["mops"] = bits(m, 24, 4, 24, 6);
  }
  return out;
}

std::vector<Json> demod(const cfvec& y, int max_err) {
  std::vector<double> ang(y.size());
  for (size_t i = 0; i < y.size(); ++i) ang[i] = std::atan2(static_cast<double>(y[i].imag()), static_cast<double>(y[i].real()));
  const std::vector<double> ph = np::unwrap(ang);
  int pat[36];
  for (int i = 0; i < 36; ++i) pat[i] = ((kSync >> (35 - i)) & 1) ? 1 : -1;
  std::vector<Json> out;
  std::vector<int64_t> seen;
  constexpr size_t kNLong = 8 * 48;
  for (int p = 0; p < kSpb; ++p) {
    std::vector<double> d;
    for (int64_t i = p; i < static_cast<int64_t>(ph.size()) - kSpb; i += kSpb)
      d.push_back(ph[static_cast<size_t>(i + kSpb)] - ph[static_cast<size_t>(i)]);
    if (d.size() < 36 + kNLong) continue;
    for (size_t k = 0; k + 36 <= d.size(); ++k) {
      int c = 0;
      for (int i = 0; i < 36; ++i) c += (d[k + i] > 0 ? 1 : -1) * pat[i];
      if (c < 36 - 2 * max_err) continue;
      const int64_t at = static_cast<int64_t>(k) * kSpb + p;
      bool dup = false;
      for (int64_t q : seen) dup = dup || std::llabs(at - q) < 30 * kSpb;
      if (dup) continue;
      double s1 = 0, s0 = 0;
      int n1 = 0, n0 = 0;
      for (int i = 0; i < 36; ++i)
        if (pat[i] > 0) s1 += d[k + i], ++n1;
        else s0 += d[k + i], ++n0;
      const double thr = (s1 / n1 + s0 / n0) / 2;
      int serr = 0;
      for (int i = 0; i < 36; ++i) serr += (d[k + i] > thr) != (pat[i] > 0);
      if (serr > max_err) continue;
      const size_t nb = std::min(kNLong, d.size() - (k + 36));
      if (nb < 8 * 30) continue;
      std::vector<uint8_t> raw((nb + 7) / 8, 0);
      for (size_t i = 0; i < nb; ++i)
        if (d[k + 36 + i] > thr) raw[i / 8] |= static_cast<uint8_t>(1 << (7 - i % 8));
      std::vector<uint8_t> msg;
      int nfix = -1;
      if (raw.size() >= 48) {
        auto r = rs_long().decode(std::vector<uint8_t>(raw.begin(), raw.begin() + 48));
        if (r && (r->first[0] >> 3) != 0) msg.assign(r->first.begin(), r->first.begin() + 34), nfix = r->second;
      }
      if (msg.empty()) {
        auto r = rs_short().decode(std::vector<uint8_t>(raw.begin(), raw.begin() + 30));
        if (r && (r->first[0] >> 3) == 0) msg.assign(r->first.begin(), r->first.begin() + 18), nfix = r->second;
      }
      if (msg.empty()) continue;
      seen.push_back(at);
      Json m = parse(msg);
      m["t_s"] = static_cast<double>(at) / kFs;
      m["rs_corrected_bytes"] = nfix;
      m["sync_bit_errors"] = serr;
      std::string hx;
      char b[3];
      for (uint8_t v : msg) std::snprintf(b, sizeof b, "%02x", v), hx += b;
      m["hex"] = hx;
      out.push_back(std::move(m));
    }
  }
  std::stable_sort(out.begin(), out.end(), [](const Json& a, const Json& b) { return a["t_s"].get<double>() < b["t_s"].get<double>(); });
  return out;
}

}  // namespace si::rx::adsb::uat
