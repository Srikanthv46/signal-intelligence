#include "core/pybytes.hpp"

#include <cstdio>

#include "core/pyfmt.hpp"

namespace si::pybytes {

int u16le(const Bytes& b, size_t i) { return b[i] | (b[i + 1] << 8); }
int be16(const Bytes& b, size_t i) { return (b[i] << 8) | b[i + 1]; }

uint64_t le_uint(const Bytes& b) {
  uint64_t v = 0;
  for (size_t i = b.size(); i-- > 0;) v = (v << 8) | static_cast<uint64_t>(b[i]);
  return v;
}

uint64_t be_uint(const Bytes& b) {
  uint64_t v = 0;
  for (size_t i = 0; i < b.size(); ++i) v = (v << 8) | static_cast<uint64_t>(b[i]);
  return v;
}

std::string hex(const Bytes& b) {
  static const char* d = "0123456789abcdef";
  std::string s;
  s.reserve(2 * b.size());
  for (size_t i = 0; i < b.size(); ++i) {
    s += d[b.data()[i] >> 4];
    s += d[b.data()[i] & 15];
  }
  return s;
}

std::string mac_str(const Bytes& b) {
  static const char* d = "0123456789abcdef";
  std::string s;
  for (size_t i = 0; i < b.size(); ++i) {
    if (i) s += ':';
    s += d[b.data()[i] >> 4];
    s += d[b.data()[i] & 15];
  }
  return s;
}

namespace {
const char* kRepl = "\xEF\xBF\xBD";  // U+FFFD

// next code point of valid UTF-8 at s[i] (advances i)
uint32_t next_cp(const std::string& s, size_t& i) {
  const unsigned char c = static_cast<unsigned char>(s[i]);
  int n = c < 0x80 ? 1 : c < 0xE0 ? 2 : c < 0xF0 ? 3 : 4;
  if (i + n > s.size()) n = 1;
  uint32_t cp = n == 1 ? c : n == 2 ? (c & 0x1F) : n == 3 ? (c & 0x0F) : (c & 0x07);
  for (int k = 1; k < n; ++k) cp = (cp << 6) | (static_cast<unsigned char>(s[i + k]) & 0x3F);
  i += n;
  return cp;
}

bool py_isspace(uint32_t c) {
  return (c >= 0x09 && c <= 0x0D) || (c >= 0x1C && c <= 0x20) || c == 0x85 || c == 0xA0 || c == 0x1680 ||
         (c >= 0x2000 && c <= 0x200A) || c == 0x2028 || c == 0x2029 || c == 0x202F || c == 0x205F || c == 0x3000;
}
}  // namespace

// CPython's UTF-8 decoder with errors="replace": one U+FFFD per maximal invalid subpart.
std::string decode_utf8(const Bytes& b) {
  const uint8_t* p = b.data();
  const size_t n = b.size();
  std::string out;
  out.reserve(n);
  size_t i = 0;
  while (i < n) {
    const uint8_t c = p[i];
    if (c < 0x80) {
      out += static_cast<char>(c);
      ++i;
      continue;
    }
    int len;
    uint8_t lo = 0x80, hi = 0xBF;
    if (c >= 0xC2 && c <= 0xDF) {
      len = 2;
    } else if (c >= 0xE0 && c <= 0xEF) {
      len = 3;
      if (c == 0xE0) lo = 0xA0;
      if (c == 0xED) hi = 0x9F;
    } else if (c >= 0xF0 && c <= 0xF4) {
      len = 4;
      if (c == 0xF0) lo = 0x90;
      if (c == 0xF4) hi = 0x8F;
    } else {
      out += kRepl;
      ++i;
      continue;
    }
    size_t k = 1;
    bool bad = false;
    for (; k < static_cast<size_t>(len); ++k) {
      if (i + k >= n) {
        bad = true;
        break;
      }
      const uint8_t d = p[i + k];
      const uint8_t l = k == 1 ? lo : 0x80, h = k == 1 ? hi : 0xBF;
      if (d < l || d > h) {
        bad = true;
        break;
      }
    }
    if (bad) {
      out += kRepl;
      i += k;
      continue;
    }
    out.append(reinterpret_cast<const char*>(p + i), len);
    i += len;
  }
  return out;
}

std::string decode_ascii(const Bytes& b) {
  std::string out;
  for (size_t i = 0; i < b.size(); ++i) {
    const uint8_t c = b.data()[i];
    if (c < 0x80) out += static_cast<char>(c);
    else out += kRepl;
  }
  return out;
}

std::string py_strip(const std::string& s) {
  std::vector<std::pair<size_t, uint32_t>> cps;  // (byte offset, code point)
  size_t i = 0;
  while (i < s.size()) {
    size_t at = i;
    cps.emplace_back(at, next_cp(s, i));
  }
  size_t a = 0, e = cps.size();
  while (a < e && py_isspace(cps[a].second)) ++a;
  while (e > a && py_isspace(cps[e - 1].second)) --e;
  if (a == e) return "";
  const size_t b0 = cps[a].first, b1 = e < cps.size() ? cps[e].first : s.size();
  return s.substr(b0, b1 - b0);
}

std::string py_strip_chars(const std::string& s, const std::string& chars) {
  size_t a = 0, e = s.size();
  while (a < e && chars.find(s[a]) != std::string::npos) ++a;
  while (e > a && chars.find(s[e - 1]) != std::string::npos) --e;
  return s.substr(a, e - a);
}

std::vector<std::string> py_split_ws(const std::string& s) {
  std::vector<std::string> out;
  std::string cur;
  bool in = false;
  size_t i = 0;
  while (i < s.size()) {
    size_t at = i;
    const uint32_t cp = next_cp(s, i);
    if (py_isspace(cp)) {
      if (in) out.push_back(cur);
      cur.clear();
      in = false;
    } else {
      cur.append(s, at, i - at);
      in = true;
    }
  }
  if (in) out.push_back(cur);
  return out;
}

std::string ascii_lower(std::string s) {
  for (auto& c : s)
    if (c >= 'A' && c <= 'Z') c = static_cast<char>(c - 'A' + 'a');
  return s;
}

std::string ascii_upper(std::string s) {
  for (auto& c : s)
    if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
  return s;
}

std::string head_codepoints(const std::string& s, size_t n) {
  size_t i = 0, k = 0;
  while (i < s.size() && k < n) {
    next_cp(s, i);
    ++k;
  }
  return s.substr(0, i);
}

size_t len_codepoints(const std::string& s) {
  size_t i = 0, k = 0;
  while (i < s.size()) {
    next_cp(s, i);
    ++k;
  }
  return k;
}

bool truthy(const Json& v) {
  switch (v.type()) {
    case Json::value_t::null: return false;
    case Json::value_t::boolean: return v.get<bool>();
    case Json::value_t::number_integer: return v.get<int64_t>() != 0;
    case Json::value_t::number_unsigned: return v.get<uint64_t>() != 0;
    case Json::value_t::number_float: return v.get<double>() != 0.0;
    case Json::value_t::string: return !v.get_ref<const std::string&>().empty();
    case Json::value_t::array:
    case Json::value_t::object: return !v.empty();
    default: return true;
  }
}

const Json& jget(const Json& d, const char* k) {
  static const Json kNull;
  if (!d.is_object()) return kNull;
  auto it = d.find(k);
  return it == d.end() ? kNull : *it;
}

std::string pstr(const Json& v) { return py_str(v); }

std::string json_key(const Json& v) {
  if (v.is_string()) return v.get<std::string>();
  if (v.is_null()) return "null";
  if (v.is_boolean()) return v.get<bool>() ? "true" : "false";
  if (v.is_number_float()) return py_float_repr(v.get<double>());
  return v.dump();
}

}  // namespace si::pybytes
