// python receivers/marine/ais.py (see there for the format and the sources).
#include "receivers/marine/ais.hpp"

#include <cstdlib>
#include <string>

#include "core/pyfmt.hpp"

namespace si::rx::marine::ais {

namespace {

const int kFlag[8] = {0, 1, 1, 1, 1, 1, 1, 0};
const char* kNav[16] = {"under way using engine", "at anchor", "not under command", "restricted manoeuvrability",
                        "constrained by draught", "moored", "aground", "engaged in fishing", "under way sailing",
                        "reserved (HSC)", "reserved (WIG)", "power-driven towing astern",
                        "power-driven pushing / towing", "reserved", "AIS-SART / MOB / EPIRB active", "not defined"};

bool flag_at(const std::vector<int>& b, size_t i) {
  if (i + 8 > b.size()) return false;
  for (int k = 0; k < 8; ++k)
    if (b[i + k] != kFlag[k]) return false;
  return true;
}

using Bits = std::vector<int>;

int64_t u(const Bits& b, size_t a, size_t n) {
  int64_t v = 0;
  for (size_t i = a; i < a + n && i < b.size(); ++i) v = (v << 1) | b[i];
  return v;
}
int64_t s(const Bits& b, size_t a, size_t n) {
  const int64_t v = u(b, a, n);
  return ((v >> (n - 1)) & 1) ? v - (int64_t(1) << n) : v;
}
std::string text(const Bits& b, size_t a, size_t nchars) {
  std::string out;
  for (size_t k = 0; k < nchars; ++k) {
    if (a + 6 * k + 6 > b.size()) break;
    const int v = static_cast<int>(u(b, a + 6 * k, 6));
    out += static_cast<char>(v < 32 ? v + 64 : v);
  }
  while (!out.empty() && (out.back() == '@' || out.back() == ' ')) out.pop_back();
  std::string r;
  for (char c : out)
    if (c != '@') r += c;
  return r;
}
std::pair<Json, Json> pos(const Bits& b, size_t a_lon, size_t a_lat) {
  const double lon = static_cast<double>(s(b, a_lon, 28)) / 600000.0, lat = static_cast<double>(s(b, a_lat, 27)) / 600000.0;
  return {std::fabs(lon) > 180 ? Json(nullptr) : Json(py_round(lon, 6)),
          std::fabs(lat) > 90 ? Json(nullptr) : Json(py_round(lat, 6))};
}

}  // namespace

uint16_t crc_x25(const std::vector<uint8_t>& o) {
  uint32_t c = 0xFFFF;
  for (uint8_t b : o) {
    c ^= b;
    for (int k = 0; k < 8; ++k) c = (c & 1) ? (c >> 1) ^ 0x8408 : c >> 1;
  }
  return static_cast<uint16_t>(c ^ 0xFFFF);
}

bool residue_ok(const std::vector<uint8_t>& o) { return o.size() >= 3 && (crc_x25(o) ^ 0xFFFF) == 0xF0B8; }

std::vector<int> unnrzi(const std::vector<int>& lv) {
  std::vector<int> out;
  for (size_t i = 1; i < lv.size(); ++i) out.push_back(lv[i] == lv[i - 1] ? 1 : 0);
  return out;
}

std::vector<std::pair<int64_t, std::vector<uint8_t>>> frames(const std::vector<int>& bits) {
  std::vector<std::pair<int64_t, std::vector<uint8_t>>> out;
  const size_t n = bits.size();
  size_t i = 6;
  auto alt = [&](size_t i0) {
    const int a[6] = {0, 1, 0, 1, 0, 1};
    bool p = true, q = true;
    for (int k = 0; k < 6; ++k) {
      p = p && bits[i0 + k] == a[k];
      q = q && bits[i0 + k] == 1 - a[k];
    }
    return p || q;
  };
  while (i + 8 <= n) {
    if (!flag_at(bits, i) || !alt(i - 6)) {
      ++i;
      continue;
    }
    std::vector<int> raw;
    int run = 0;
    size_t j = i + 8;
    bool ok = false;
    while (j < n) {
      if (flag_at(bits, j)) {
        ok = true;
        break;
      }
      const int b = bits[j];
      if (run == 5) {
        if (b == 1) break;
        run = 0;
        ++j;
        continue;
      }
      raw.push_back(b);
      run = b ? run + 1 : 0;
      ++j;
      if (raw.size() > 1100) break;
    }
    if (ok && raw.size() >= 40 && raw.size() % 8 == 0) {
      std::vector<uint8_t> o;
      for (size_t k = 0; k < raw.size(); k += 8) {
        int v = 0;
        for (int q = 0; q < 8; ++q) v |= raw[k + q] << q;
        o.push_back(static_cast<uint8_t>(v));
      }
      if (residue_ok(o)) {
        o.resize(o.size() - 2);
        out.emplace_back(static_cast<int64_t>(i), o);
        i = j + 8;
        continue;
      }
    }
    ++i;
  }
  return out;
}

Json fields(const std::vector<uint8_t>& octs, std::optional<size_t> nbits) {
  Bits b;
  for (uint8_t o : octs)
    for (int i = 0; i < 8; ++i) b.push_back((o >> (7 - i)) & 1);
  if (nbits && *nbits && *nbits < b.size()) b.resize(*nbits);
  const int64_t t = u(b, 0, 6), mmsi = u(b, 8, 30);
  Json out = Json::object();
  out["msg_type"] = t;
  out["repeat"] = u(b, 6, 2);
  const std::string ms = fmt("%09lld", static_cast<long long>(mmsi));
  out["mmsi"] = ms;
  // MID: the first three digits of a ship; 00MID... base station, 99MID... aid to navigation, 111MID... SAR aircraft
  const int mid = ms.rfind("111", 0) == 0 ? std::atoi(ms.substr(3, 3).c_str())
                  : (ms.rfind("00", 0) == 0 || ms.rfind("99", 0) == 0) ? std::atoi(ms.substr(2, 3).c_str())
                                                                         : std::atoi(ms.substr(0, 3).c_str());
  if (200 <= mid && mid <= 799) {
    out["mid"] = mid;
    if (mid == 419) out["flag"] = "India";
  }
  const size_t nb = b.size();
  auto opt = [](int64_t v, int64_t na, double scale) { return v == na ? Json(nullptr) : Json(static_cast<double>(v) / scale); };
  if ((t == 1 || t == 2 || t == 3) && nb >= 168) {
    const auto [lon, lat] = pos(b, 61, 89);
    const int64_t sog = u(b, 50, 10), cog = u(b, 116, 12), hdg = u(b, 128, 9);
    out["kind"] = "class A position report";
    out["nav_status"] = kNav[u(b, 38, 4)];
    out["sog_kn"] = opt(sog, 1023, 10);
    out["accuracy"] = static_cast<bool>(b[60]);
    out["lon"] = lon;
    out["lat"] = lat;
    out["cog_deg"] = opt(cog, 3600, 10);
    out["heading_deg"] = hdg == 511 ? Json(nullptr) : Json(hdg);
    out["utc_second"] = u(b, 137, 6);
  } else if ((t == 4 || t == 11) && nb >= 168) {
    const auto [lon, lat] = pos(b, 79, 107);
    out["kind"] = t == 4 ? "base station report" : "UTC / date response";
    out["utc"] = fmt("%04lld-%02lld-%02lldT%02lld:%02lld:%02lldZ", static_cast<long long>(u(b, 38, 14)),
                     static_cast<long long>(u(b, 52, 4)), static_cast<long long>(u(b, 56, 5)),
                     static_cast<long long>(u(b, 61, 5)), static_cast<long long>(u(b, 66, 6)),
                     static_cast<long long>(u(b, 72, 6)));
    out["lon"] = lon;
    out["lat"] = lat;
  } else if (t == 5 && nb >= 420) {
    out["kind"] = "static and voyage data";
    out["imo"] = u(b, 40, 30);
    out["callsign"] = text(b, 70, 7);
    out["name"] = text(b, 112, 20);
    out["ship_type"] = u(b, 232, 8);
    out["length_m"] = u(b, 240, 9) + u(b, 249, 9);
    out["beam_m"] = u(b, 258, 6) + u(b, 264, 6);
    out["draught_m"] = static_cast<double>(u(b, 294, 8)) / 10;
    out["destination"] = text(b, 302, 20);
  } else if (t == 18 && nb >= 168) {
    const auto [lon, lat] = pos(b, 57, 85);
    const int64_t sog = u(b, 46, 10), cog = u(b, 112, 12), hdg = u(b, 124, 9);
    out["kind"] = "class B position report";
    out["sog_kn"] = opt(sog, 1023, 10);
    out["lon"] = lon;
    out["lat"] = lat;
    out["cog_deg"] = opt(cog, 3600, 10);
    out["heading_deg"] = hdg == 511 ? Json(nullptr) : Json(hdg);
  } else if (t == 21 && nb >= 272) {
    const auto [lon, lat] = pos(b, 164, 192);
    out["kind"] = "aid to navigation";
    out["aid_type"] = u(b, 38, 5);
    out["name"] = text(b, 43, 20);
    out["lon"] = lon;
    out["lat"] = lat;
    out["virtual"] = static_cast<bool>(b[269]);
  } else if (t == 24 && nb >= 160) {
    const int64_t part = u(b, 38, 2);
    out["kind"] = "class B static data";
    if (part == 0) {
      out["name"] = text(b, 40, 20);
    } else if (part == 1 && nb >= 162) {
      out["ship_type"] = u(b, 40, 8);
      out["callsign"] = text(b, 90, 7);
    }
  }
  return out;
}

}  // namespace si::rx::marine::ais
