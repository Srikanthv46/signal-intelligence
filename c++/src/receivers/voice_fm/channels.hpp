// Walkie-talkie channel plans and the sub-audible signalling tables (CTCSS tones, DCS codes)
// (python: receivers/voice_fm/channels.py).
#pragma once

#include <optional>
#include <string>
#include <vector>

#include "core/types.hpp"

namespace si::rx::voice_fm {

struct Plan {
  std::string service;
  std::vector<std::pair<int, double>> channels;   // (channel, nominal Hz)
  double bw_hz;
  double spacing_hz;
};
const std::vector<Plan>& plans();          // PMR446, LPD433, FRS/GMRS (python dict order)
const std::vector<double>& ctcss_tones();  // EIA/TIA-603 + common extras
const std::vector<std::string>& dcs_codes();  // the 104 standard codes (octal strings)
bool is_standard_dcs(const std::string& code);

// python match_channel(f_hz, ppm): {"service", "channel", "nominal_mhz", "channel_bw_khz",
// "offset_hz", "offset_ppm", "receiver_ppm_applied"} or nullopt
std::optional<Json> match_channel(double f_hz, double ppm = 0.0);

}  // namespace si::rx::voice_fm
