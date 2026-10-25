#include "core/pyfmt.hpp"

#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace si {

std::string py_float_repr(double x) {
  if (std::isnan(x)) return "nan";
  if (std::isinf(x)) return x > 0 ? "inf" : "-inf";
  if (x == 0) return std::signbit(x) ? "-0.0" : "0.0";
  // shortest %.{p}e that round-trips -> digits and exponent
  char buf[64];
  int p = 0;
  for (p = 0; p < 17; ++p) {
    std::snprintf(buf, sizeof buf, "%.*e", p, x);
    if (std::strtod(buf, nullptr) == x) break;
  }
  // buf = [-]d.ddddde[+-]XX
  std::string s(buf);
  const bool neg = s[0] == '-';
  if (neg) s = s.substr(1);
  const size_t epos = s.find('e');
  std::string mant = s.substr(0, epos);
  const int exp10 = std::atoi(s.c_str() + epos + 1);
  std::string digits;
  for (char c : mant)
    if (c != '.') digits += c;
  while (digits.size() > 1 && digits.back() == '0') digits.pop_back();
  std::string out;
  if (exp10 < -4 || exp10 >= 16) {
    out = digits.substr(0, 1);
    if (digits.size() > 1) out += "." + digits.substr(1);
    char e[16];
    std::snprintf(e, sizeof e, "e%c%02d", exp10 < 0 ? '-' : '+', std::abs(exp10));
    out += e;
  } else if (exp10 < 0) {
    out = "0." + std::string(static_cast<size_t>(-exp10 - 1), '0') + digits;
  } else {
    const size_t ip = static_cast<size_t>(exp10) + 1;
    if (digits.size() <= ip) {
      out = digits + std::string(ip - digits.size(), '0') + ".0";
    } else {
      out = digits.substr(0, ip) + "." + digits.substr(ip);
    }
  }
  return neg ? "-" + out : out;
}

namespace {

std::string repr_str(const std::string& s) {
  const bool has_sq = s.find('\'') != std::string::npos, has_dq = s.find('"') != std::string::npos;
  const char q = (has_sq && !has_dq) ? '"' : '\'';
  std::string out(1, q);
  for (unsigned char c : s) {
    if (c == '\\') out += "\\\\";
    else if (c == static_cast<unsigned char>(q)) out += std::string("\\") + q;
    else if (c == '\n') out += "\\n";
    else if (c == '\r') out += "\\r";
    else if (c == '\t') out += "\\t";
    else if (c < 0x20 || c == 0x7f) {
      char b[8];
      std::snprintf(b, sizeof b, "\\x%02x", c);
      out += b;
    } else out += static_cast<char>(c);
  }
  out += q;
  return out;
}

}  // namespace

std::string py_repr(const Json& v) {
  switch (v.type()) {
    case Json::value_t::null: return "None";
    case Json::value_t::boolean: return v.get<bool>() ? "True" : "False";
    case Json::value_t::number_integer: return std::to_string(v.get<int64_t>());
    case Json::value_t::number_unsigned: return std::to_string(v.get<uint64_t>());
    case Json::value_t::number_float: return py_float_repr(v.get<double>());
    case Json::value_t::string: return repr_str(v.get<std::string>());
    case Json::value_t::array: {
      std::string out = "[";
      bool first = true;
      for (const auto& x : v) {
        if (!first) out += ", ";
        first = false;
        out += py_repr(x);
      }
      return out + "]";
    }
    case Json::value_t::object: {
      std::string out = "{";
      bool first = true;
      for (auto it = v.begin(); it != v.end(); ++it) {
        if (!first) out += ", ";
        first = false;
        out += repr_str(it.key()) + ": " + py_repr(it.value());
      }
      return out + "}";
    }
    default: return v.dump();
  }
}

std::string py_str(const Json& v) { return v.is_string() ? v.get<std::string>() : py_repr(v); }

std::string fmt(const char* format, ...) {
  va_list ap;
  va_start(ap, format);
  va_list ap2;
  va_copy(ap2, ap);
  const int n = std::vsnprintf(nullptr, 0, format, ap);
  va_end(ap);
  std::string out(n > 0 ? static_cast<size_t>(n) : 0, '\0');
  if (n > 0) std::vsnprintf(out.data(), static_cast<size_t>(n) + 1, format, ap2);
  va_end(ap2);
  return out;
}

}  // namespace si
