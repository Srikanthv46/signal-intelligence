// Report numbers and strings must read exactly as Python writes them (expected values
// printed by CPython 3.10).
#include <doctest/doctest.h>

#include "core/pyfmt.hpp"
#include "core/timeutil.hpp"
#include "core/types.hpp"

using namespace si;

TEST_CASE("py_float_repr matches repr(float)") {
  const std::pair<double, const char*> cases[] = {
      {0.1, "0.1"},
      {100.0, "100.0"},
      {1e-05, "1e-05"},
      {1e16, "1e+16"},
      {123456789.123, "123456789.123"},
      {0.0001, "0.0001"},
      {1.5e-07, "1.5e-07"},
      {1234567890123456.0, "1234567890123456.0"},
      {12345678901234567.0, "1.2345678901234568e+16"},
      {-2.5, "-2.5"},
      {5e-324, "5e-324"},
      {0.30000000000000004, "0.30000000000000004"},
      {299792458.0, "299792458.0"},
      {1.7976931348623157e308, "1.7976931348623157e+308"},
  };
  for (const auto& [v, s] : cases) CHECK(py_float_repr(v) == s);
}

TEST_CASE("py_round matches round()") {
  CHECK(py_round(2.675, 2) == 2.67);
  CHECK(py_round(0.125, 2) == 0.12);
  CHECK(py_round(1234.5678, -2) == 1200.0);
  CHECK(py_round(0.0000049, 5) == 0.0);
  CHECK(py_round_int(-1.5) == -2);
  CHECK(py_round_int(2.5) == 2);
  CHECK(py_round_int(3.5) == 4);
  CHECK(py_round_int(-2.5) == -2);
}

TEST_CASE("py_repr matches repr(dict)") {
  Json d = {{"a", 1}, {"b", Json::array({1.5, "x"})}, {"c", nullptr}, {"d", true}, {"it's", "say \"hi\""}};
  CHECK(py_repr(d) == "{'a': 1, 'b': [1.5, 'x'], 'c': None, 'd': True, \"it's\": 'say \"hi\"'}");
  CHECK(fmt("%.1f MHz", 2437.04) == "2437.0 MHz");
}

TEST_CASE("UtcTime matches datetime.isoformat") {
  auto t = UtcTime::parse("2026-10-04T08:40:21.800158+00:00");
  REQUIRE(t);
  CHECK(t->isoformat() == "2026-10-04T08:40:21.800158+00:00");
  CHECK(t->plus(0.5).isoformat() == "2026-10-04T08:40:22.300158+00:00");
  CHECK(t->plus(1e-7).isoformat() == "2026-10-04T08:40:21.800158+00:00");
  CHECK(t->plus(0.1998420).isoformat() == "2026-10-04T08:40:22+00:00");
  auto z = UtcTime::parse("2026-10-09T08:00:00Z");
  REQUIRE(z);
  CHECK(z->iso_z() == "2026-10-09T08:00:00Z");
  auto off = UtcTime::parse("2026-10-09T13:30:00+05:30");
  REQUIRE(off);
  CHECK(off->isoformat() == "2026-10-09T08:00:00+00:00");
  CHECK(UtcTime::from_unix(1760000000.25).iso_z() == "2025-10-09T08:53:20.250000Z");
  CHECK_FALSE(UtcTime::parse("yesterday"));
}
