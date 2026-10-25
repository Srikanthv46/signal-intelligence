// python core/hdlc.py (see there for the formats and the sources).
#include "core/hdlc.hpp"

#include <cmath>
#include <cstdlib>

namespace si::hdlc {

namespace {

const int kFlag[8] = {0, 1, 1, 1, 1, 1, 1, 0};

bool flag_at(const std::vector<int>& b, size_t i) {
  if (i + 8 > b.size()) return false;
  for (int k = 0; k < 8; ++k)
    if (b[i + k] != kFlag[k]) return false;
  return true;
}

bool py_space(unsigned char c) {
  return c == ' ' || (c >= 0x09 && c <= 0x0D) || (c >= 0x1C && c <= 0x1F) || c == 0x85 || c == 0xA0;
}
std::string strip(const std::string& s) {
  size_t a = 0, b = s.size();
  while (a < b && py_space(static_cast<unsigned char>(s[a]))) ++a;
  while (b > a && py_space(static_cast<unsigned char>(s[b - 1]))) --b;
  return s.substr(a, b - a);
}
std::string rstrip(const std::string& s) {
  size_t b = s.size();
  while (b > 0 && py_space(static_cast<unsigned char>(s[b - 1]))) --b;
  return s.substr(0, b);
}
bool isdigits(const std::string& s) {
  if (s.empty()) return false;
  for (unsigned char c : s)
    if (c < '0' || c > '9') return false;
  return true;
}
int ordc(const std::string& s, size_t i) { return static_cast<unsigned char>(s[i]); }

std::string call(const uint8_t* b) {
  std::string c;
  for (int i = 0; i < 6; ++i) c += static_cast<char>(b[i] >> 1);
  c = rstrip(c);
  const int ssid = (b[6] >> 1) & 15;
  return ssid ? c + "-" + std::to_string(ssid) : c;
}

// python _ddmm(s, hemi, deg_digits): nullopt on a parse error
std::optional<double> ddmm(std::string s, char hemi, size_t deg) {
  for (auto& c : s)
    if (c == ' ') c = '0';
  const std::string a = s.substr(0, deg), b = s.substr(deg);
  if (!isdigits(a)) return std::nullopt;
  // python float(): digits with an optional point
  size_t dots = 0, digs = 0;
  for (char c : b) {
    if (c == '.') ++dots;
    else if (c >= '0' && c <= '9') ++digs;
    else return std::nullopt;
  }
  if (dots > 1 || digs == 0) return std::nullopt;
  const double v = static_cast<double>(std::atoll(a.c_str())) + std::strtod(b.c_str(), nullptr) / 60;
  return py_round((hemi == 'S' || hemi == 'W') ? -v : v, 5);
}

Json uncompressed(const std::string& p) {
  Json out = Json::object();
  if (p.size() < 19 || (p[7] != 'N' && p[7] != 'S') || (p[17] != 'E' && p[17] != 'W')) return out;
  const auto lat = ddmm(p.substr(0, 7), p[7], 2), lon = ddmm(p.substr(9, 8), p[17], 3);
  if (!lat || !lon || std::fabs(*lat) > 90 || std::fabs(*lon) > 180) return out;
  out["lat"] = *lat;
  out["lon"] = *lon;
  out["symbol"] = latin1_utf8(std::string(1, p[8]) + p[18]);
  out["format"] = "uncompressed";
  std::string rest = p.substr(19);
  if (rest.size() >= 7 && rest[3] == '/' && isdigits(rest.substr(0, 3)) && isdigits(rest.substr(4, 3))) {
    out["course_deg"] = std::atoi(rest.substr(0, 3).c_str());
    out["speed_kn"] = std::atoi(rest.substr(4, 3).c_str());
    rest = rest.substr(7);
  }
  out["comment"] = latin1_utf8(rest);
  return out;
}

Json compressed(const std::string& p) {
  Json out = Json::object();
  if (p.size() < 13) return out;
  int64_t v[8];
  for (int i = 0; i < 8; ++i) {
    v[i] = ordc(p, 1 + static_cast<size_t>(i)) - 33;
    if (v[i] < 0 || v[i] > 90) return out;
  }
  const double lat = 90 - static_cast<double>(v[0] * 91 * 91 * 91 + v[1] * 91 * 91 + v[2] * 91 + v[3]) / 380926;
  const double lon = -180 + static_cast<double>(v[4] * 91 * 91 * 91 + v[5] * 91 * 91 + v[6] * 91 + v[7]) / 190463;
  out["lat"] = py_round(lat, 5);
  out["lon"] = py_round(lon, 5);
  out["symbol"] = latin1_utf8(std::string(1, p[0]) + p[9]);
  out["format"] = "compressed";
  const int c = ordc(p, 10) - 33, s = ordc(p, 11) - 33;
  if (p[10] != ' ' && 0 <= c && c <= 89) {
    out["course_deg"] = c * 4;
    out["speed_kn"] = py_round(std::pow(1.08, s) - 1, 1);
  }
  out["comment"] = latin1_utf8(p.substr(13));
  return out;
}

Json mic_e(const std::string& dest, const std::string& info) {
  Json out = Json::object();
  const std::string d = dest.substr(0, dest.find('-'));
  if (d.size() != 6 || info.size() < 9) return out;
  int dg[6];
  for (int i = 0; i < 6; ++i) {
    const char ch = d[i];
    if (ch >= '0' && ch <= '9') dg[i] = ch - 48;
    else if (ch >= 'A' && ch <= 'J') dg[i] = ch - 65;
    else if (ch >= 'P' && ch <= 'Y') dg[i] = ch - 80;
    else if (ch == 'K' || ch == 'L' || ch == 'Z') dg[i] = 0;
    else return out;
  }
  const double lat = dg[0] * 10 + dg[1] + (dg[2] * 10 + dg[3] + (dg[4] * 10 + dg[5]) / 100.0) / 60;
  const bool north = d[3] >= 'P', off100 = d[4] >= 'P', west = d[5] >= 'P';
  int deg = ordc(info, 1) - 28 + (off100 ? 100 : 0);
  if (180 <= deg && deg <= 189) deg -= 80;
  else if (190 <= deg && deg <= 199) deg -= 190;
  int mins = ordc(info, 2) - 28;
  if (mins >= 60) mins -= 60;
  const int hund = ordc(info, 3) - 28;
  const double lon = deg + (mins + hund / 100.0) / 60;
  const int sp = ordc(info, 4) - 28, dc = ordc(info, 5) - 28, se = ordc(info, 6) - 28;
  int speed = sp * 10 + (dc >= 0 ? dc / 10 : -((-dc + 9) / 10));
  if (speed >= 800) speed -= 800;
  int course = (((dc % 10) + 10) % 10) * 100 + se;
  if (course >= 400) course -= 400;
  if (!(0 <= lat && lat <= 90 && 0 <= lon && lon <= 180)) return out;
  out["lat"] = py_round(north ? lat : -lat, 5);
  out["lon"] = py_round(west ? -lon : lon, 5);
  out["symbol"] = latin1_utf8(std::string(1, info[8]) + info[7]);
  out["format"] = "Mic-E";
  out["course_deg"] = course;
  out["speed_kn"] = speed;
  out["comment"] = latin1_utf8(info.substr(9));
  return out;
}

Json position(const std::string& body) {
  if (body.empty()) return Json::object();
  return (body[0] >= '0' && body[0] <= '9') || body[0] == ' ' ? uncompressed(body) : compressed(body);
}

void update(Json& out, const Json& add) {
  for (auto it = add.begin(); it != add.end(); ++it) out[it.key()] = it.value();
}

}  // namespace

std::string latin1_utf8(const std::string& s) {
  std::string o;
  for (unsigned char c : s) {
    if (c < 0x80) {
      o += static_cast<char>(c);
    } else {
      o += static_cast<char>(0xC0 | (c >> 6));
      o += static_cast<char>(0x80 | (c & 0x3F));
    }
  }
  return o;
}

uint16_t crc_x25(const std::vector<uint8_t>& o) {
  uint32_t c = 0xFFFF;
  for (uint8_t b : o) {
    c ^= b;
    for (int k = 0; k < 8; ++k) c = (c & 1) ? (c >> 1) ^ 0x8408 : c >> 1;
  }
  return static_cast<uint16_t>(c ^ 0xFFFF);
}

std::vector<int> unnrzi(const std::vector<int>& lv) {
  std::vector<int> out;
  for (size_t i = 1; i < lv.size(); ++i) out.push_back(lv[i] == lv[i - 1] ? 1 : 0);
  return out;
}

std::vector<std::pair<int64_t, std::vector<uint8_t>>> deframe(const std::vector<int>& bits, size_t min_octets,
                                                               size_t max_bits) {
  std::vector<std::pair<int64_t, std::vector<uint8_t>>> out;
  const size_t n = bits.size();
  size_t i = 0;
  while (i + 8 <= n) {
    if (!flag_at(bits, i)) {
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
      if (raw.size() > max_bits) break;
    }
    if (ok && raw.size() >= 8 * min_octets && raw.size() % 8 == 0) {
      std::vector<uint8_t> o;
      for (size_t k = 0; k < raw.size(); k += 8) {
        int v = 0;
        for (int q = 0; q < 8; ++q) v |= raw[k + q] << q;
        o.push_back(static_cast<uint8_t>(v));
      }
      if ((crc_x25(o) ^ 0xFFFF) == 0xF0B8) {
        o.resize(o.size() - 2);
        out.emplace_back(static_cast<int64_t>(i), o);
        i = j;
        continue;
      }
    }
    i = (ok && raw.empty()) ? j : i + 1;
  }
  return out;
}

std::optional<Ax25> ax25(const std::vector<uint8_t>& o) {
  std::vector<const uint8_t*> addrs;
  size_t k = 0;
  while (k + 7 <= o.size()) {
    const uint8_t* a = o.data() + k;
    for (int i = 0; i < 6; ++i) {
      const int c = a[i] >> 1;
      if (!((0x40 <= c && c < 0x5B) || c == 0x20 || (0x30 <= c && c <= 0x39))) return std::nullopt;
    }
    addrs.push_back(a);
    k += 7;
    if (a[6] & 1) break;
  }
  if (addrs.size() < 2 || !(addrs.back()[6] & 1) || k + 2 > o.size() || o[k] != 0x03 || o[k + 1] != 0xF0)
    return std::nullopt;
  Ax25 f;
  f.dest = call(addrs[0]);
  f.src = call(addrs[1]);
  for (size_t q = 2; q < addrs.size(); ++q) f.path.push_back(call(addrs[q]) + ((addrs[q][6] & 0x80) ? "*" : ""));
  f.info.assign(o.begin() + static_cast<long>(k) + 2, o.end());
  return f;
}

Json aprs(const std::string& dest, const std::string& info) {
  Json out = Json::object();
  if (info.empty()) return out;
  const char t = info[0];
  out["data_type"] = latin1_utf8(std::string(1, t));
  if (t == '!' || t == '=') {
    update(out, position(info.substr(1)));
  } else if ((t == '/' || t == '@') && info.size() >= 8) {
    out["timestamp"] = latin1_utf8(info.substr(1, 7));
    update(out, position(info.substr(8)));
  } else if (t == '`' || t == '\'') {
    update(out, mic_e(dest, info));
  } else if (t == '>') {
    out["status"] = latin1_utf8(info.substr(1));
  } else if (t == ':' && info.size() >= 11 && info[10] == ':') {
    out["addressee"] = latin1_utf8(strip(info.substr(1, 9)));
    out["message"] = latin1_utf8(info.substr(11));
  } else if (t == ';' && info.size() >= 18) {
    out["object"] = latin1_utf8(strip(info.substr(1, 9)));
    update(out, uncompressed(info.substr(18)));
  }
  return out;
}

}  // namespace si::hdlc
