// Analysis -> response (schema/response.schema.json); python: api/response.py.
// Built from the per-block report (report/writer), so the two always agree.
#pragma once

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "analyze.hpp"
#include "core/timeutil.hpp"

namespace si::api {

struct BlockMeta {
  double center_frequency_hz = 0;
  double sample_rate_hz = 0;
  std::optional<double> bandwidth_hz;
  std::optional<double> gain_db;
  std::optional<double> receiver_ppm;
  std::optional<UtcTime> start;
  std::optional<std::string> source_id;
  std::optional<std::string> request_id;
  size_t n_samples = 0;
};

// receiver package -> its plugin names, in registry order
using Packages = std::vector<std::pair<std::string, std::vector<std::string>>>;

Json build_response(const Json& report, const Analysis& a, const BlockMeta& block, const Json& engine,
                    const Packages& packages, double total_s, bool compact);

// Same shape, empty results, the reason in errors[] (status insufficient_input / error).
Json failure_response(const BlockMeta& block, const Json& engine, const std::string& status,
                      const std::string& code, const std::string& message, double total_s);

}  // namespace si::api
