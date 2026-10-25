#include "core/types.hpp"

#include <cfenv>
#include <cstdio>
#include <cstdlib>

namespace si {

double py_round(double x, int ndigits) {
  if (!std::isfinite(x)) return x;
  if (ndigits > 22 || ndigits < -22) return x;
  char buf[512];
  if (ndigits >= 0) {
    std::snprintf(buf, sizeof buf, "%.*f", ndigits, x);
    return std::strtod(buf, nullptr);
  }
  // negative ndigits: round to tens/hundreds..., ties to even, as CPython does
  double p = std::pow(10.0, -ndigits);
  double y = x / p;
  double r = std::nearbyint(y);  // default rounding mode: to nearest, ties to even
  return r * p;
}

int64_t py_round_int(double x) {
  const int old = std::fegetround();
  if (old != FE_TONEAREST) std::fesetround(FE_TONEAREST);
  double r = std::nearbyint(x);
  if (old != FE_TONEAREST) std::fesetround(old);
  return static_cast<int64_t>(r);
}

Json jnum(double x, int ndigits) {
  if (!std::isfinite(x)) return Json(nullptr);
  return Json(py_round(x, ndigits));
}

}  // namespace si
