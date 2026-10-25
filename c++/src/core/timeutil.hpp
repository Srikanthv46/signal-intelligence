// UTC timestamps as the Python reference prints them (datetime.isoformat()).
#pragma once

#include <cstdint>
#include <optional>
#include <string>

namespace si {

struct UtcTime {
  int64_t epoch_s = 0;   // seconds since 1970-01-01T00:00:00Z
  int64_t micro = 0;     // 0..999999

  // Parse ISO 8601 ("2026-10-04T08:40:21.800158+00:00", "...Z", offsets, no fraction).
  static std::optional<UtcTime> parse(const std::string& iso);
  static UtcTime from_unix(double seconds);
  // this + s seconds (rounded to the microsecond like datetime + timedelta(seconds=s))
  UtcTime plus(double s) const;
  // "YYYY-MM-DDTHH:MM:SS[.ffffff]+00:00" (python isoformat of an aware UTC datetime)
  std::string isoformat() const;
  // same with a trailing "Z" instead of "+00:00" (the response format)
  std::string iso_z() const;
};

}  // namespace si
