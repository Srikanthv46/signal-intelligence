// The report side of the 802.11b plugin (python plugin_dsss.py, minus the pipeline glue that
// lives in plugin.cpp): one decode result per channel with the PLCP statistics, the BSS and
// station tables and one line per frame, or a "Barker-only" partial result.
#pragma once

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "core/types.hpp"

namespace si::rx::wifi::dsss {

const std::vector<std::string>& steps();  // python plugin_dsss.STEPS

// python decode.base.result() fields, filled by the summaries below
struct Summary {
  std::string status, method;
  std::vector<std::string> steps;
  Json fields = Json::object();
  std::string reason;
  std::optional<std::string> label;
  std::optional<double> confidence;
  std::optional<std::string> modulation, evidence;
  std::optional<std::pair<double, double>> span;
};

// python WifiDsssDecoder.summarise(mine, n_all, ctx): block_s > 0 adds the field.py link profiles (PLAN W1)
Summary summarise(const std::vector<const Json*>& mine, size_t n_all, double block_s = 0);
// python WifiDsssDecoder.barker_only(regs) (regs non-empty)
Summary barker_only(const std::vector<const Json*>& regs);
// the not_possible result when a signal holds neither Barker spreading nor a PLCP header
Summary no_dsss();

Json frame_line(const Json& p);
Json bss_table(const std::vector<const Json*>& frames);
Json station_table(const std::vector<const Json*>& frames);

}  // namespace si::rx::wifi::dsss
