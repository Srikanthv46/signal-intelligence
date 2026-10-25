#include "core/timeutil.hpp"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "core/types.hpp"

namespace si {
namespace {

// days since 1970-01-01 of a proleptic Gregorian date (H. Hinnant's algorithm)
int64_t days_from_civil(int64_t y, unsigned m, unsigned d) {
  y -= m <= 2;
  const int64_t era = (y >= 0 ? y : y - 399) / 400;
  const unsigned yoe = static_cast<unsigned>(y - era * 400);
  const unsigned doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;
  const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
  return era * 146097 + static_cast<int64_t>(doe) - 719468;
}

void civil_from_days(int64_t z, int64_t& y, unsigned& m, unsigned& d) {
  z += 719468;
  const int64_t era = (z >= 0 ? z : z - 146096) / 146097;
  const unsigned doe = static_cast<unsigned>(z - era * 146097);
  const unsigned yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
  y = static_cast<int64_t>(yoe) + era * 400;
  const unsigned doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
  const unsigned mp = (5 * doy + 2) / 153;
  d = doy - (153 * mp + 2) / 5 + 1;
  m = mp + (mp < 10 ? 3 : -9);
  y += m <= 2;
}

}  // namespace

std::optional<UtcTime> UtcTime::parse(const std::string& s) {
  int Y, M, D, h = 0, mi = 0, sec = 0;
  if (s.size() < 10 || std::sscanf(s.c_str(), "%4d-%2d-%2d", &Y, &M, &D) != 3) return std::nullopt;
  size_t p = 10;
  int64_t micro = 0;
  int64_t offset_s = 0;
  auto two = [&](size_t at, int& v) {
    if (at + 2 > s.size() || !std::isdigit(static_cast<unsigned char>(s[at])) ||
        !std::isdigit(static_cast<unsigned char>(s[at + 1])))
      return false;
    v = (s[at] - '0') * 10 + (s[at + 1] - '0');
    return true;
  };
  if (p < s.size() && (s[p] == 'T' || s[p] == ' ')) {
    if (!two(p + 1, h) || p + 3 >= s.size() || s[p + 3] != ':' || !two(p + 4, mi)) return std::nullopt;
    p += 6;                                               // past "THH:MM"
    if (p < s.size() && s[p] == ':') {
      if (!two(p + 1, sec)) return std::nullopt;
      p += 3;
    }
    if (p < s.size() && s[p] == '.') {
      size_t q = p + 1;
      int64_t frac = 0, digits = 0;
      while (q < s.size() && std::isdigit(static_cast<unsigned char>(s[q]))) {
        if (digits < 6) { frac = frac * 10 + (s[q] - '0'); ++digits; }
        ++q;
      }
      while (digits < 6) { frac *= 10; ++digits; }
      micro = frac;
      p = q;
    }
    if (p < s.size()) {
      if (s[p] == 'Z' || s[p] == 'z') {
        offset_s = 0;
      } else if (s[p] == '+' || s[p] == '-') {
        int oh = 0, om = 0;
        if (std::sscanf(s.c_str() + p + 1, "%2d:%2d", &oh, &om) < 1) return std::nullopt;
        offset_s = (oh * 3600 + om * 60) * (s[p] == '-' ? -1 : 1);
      }
    }
  }
  UtcTime t;
  t.epoch_s = days_from_civil(Y, static_cast<unsigned>(M), static_cast<unsigned>(D)) * 86400 + h * 3600 + mi * 60 +
              sec - offset_s;
  t.micro = micro;
  return t;
}

UtcTime UtcTime::from_unix(double seconds) {
  UtcTime t;
  double fl = std::floor(seconds);
  t.epoch_s = static_cast<int64_t>(fl);
  t.micro = py_round_int((seconds - fl) * 1e6);
  if (t.micro >= 1000000) { t.epoch_s += 1; t.micro -= 1000000; }
  return t;
}

UtcTime UtcTime::plus(double s) const {
  int64_t us = py_round_int(s * 1e6);
  UtcTime t = *this;
  int64_t total = t.micro + us;
  int64_t carry = total >= 0 ? total / 1000000 : -((-total + 999999) / 1000000);
  t.epoch_s += carry;
  t.micro = total - carry * 1000000;
  return t;
}

std::string UtcTime::isoformat() const {
  int64_t days = epoch_s >= 0 ? epoch_s / 86400 : -((-epoch_s + 86399) / 86400);
  int64_t rem = epoch_s - days * 86400;
  int64_t y;
  unsigned m, d;
  civil_from_days(days, y, m, d);
  char buf[64];
  if (micro)
    std::snprintf(buf, sizeof buf, "%04lld-%02u-%02uT%02lld:%02lld:%02lld.%06lld+00:00", static_cast<long long>(y), m,
                  d, static_cast<long long>(rem / 3600), static_cast<long long>(rem % 3600 / 60),
                  static_cast<long long>(rem % 60), static_cast<long long>(micro));
  else
    std::snprintf(buf, sizeof buf, "%04lld-%02u-%02uT%02lld:%02lld:%02lld+00:00", static_cast<long long>(y), m, d,
                  static_cast<long long>(rem / 3600), static_cast<long long>(rem % 3600 / 60),
                  static_cast<long long>(rem % 60));
  return buf;
}

std::string UtcTime::iso_z() const {
  std::string s = isoformat();
  return s.substr(0, s.size() - 6) + "Z";
}

}  // namespace si
