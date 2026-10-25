#include "core/model.hpp"

#include <algorithm>
#include <limits>

namespace si {

double Signal::t0() const {
  double t = std::numeric_limits<double>::infinity();
  for (const auto& b : bursts) t = std::min(t, b.t0);
  return t;
}

double Signal::t1() const {
  double t = -std::numeric_limits<double>::infinity();
  for (const auto& b : bursts) t = std::max(t, b.t1);
  return t;
}

bool Signal::has_flag(const std::string& f) const {
  return std::find(flags.begin(), flags.end(), f) != flags.end();
}

}  // namespace si
