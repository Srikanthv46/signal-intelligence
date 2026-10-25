// Protocol tables of the rc_link detector (python detector.py: ELRS_LORA, ELRS_DOMAINS,
// FSK_CLASSES, _nominal_channel). Insertion order as the Python dicts.
#pragma once

#include <cmath>
#include <optional>
#include <string>
#include <vector>

#include "core/types.hpp"

namespace si::rx::rc_link::detail {

// ExpressLRS air-rate modes: name, band, sf, bw, cr_den, preamble, payload bytes, hop_every, interval_us
struct ElrsMode {
  const char* name;
  const char* band;
  int sf;
  double bw;
  int cr_den, preamble, payload, hop_every, interval_us;
};
inline const std::vector<ElrsMode>& elrs_lora() {
  static const std::vector<ElrsMode> t = {
      {"2g4_500hz", "2G4", 5, 812.5e3, 6, 12, 8, 4, 2000},      {"2g4_333hz_8ch", "2G4", 5, 812.5e3, 8, 12, 13, 4, 3003},
      {"2g4_250hz", "2G4", 6, 812.5e3, 8, 14, 8, 4, 4000},      {"2g4_150hz", "2G4", 7, 812.5e3, 8, 12, 8, 4, 6666},
      {"2g4_100hz_8ch", "2G4", 7, 812.5e3, 8, 12, 13, 4, 10000}, {"2g4_50hz", "2G4", 8, 812.5e3, 8, 12, 8, 2, 20000},
      {"900_200hz", "900", 6, 500e3, 7, 8, 8, 4, 5000},         {"900_100hz_8ch", "900", 6, 500e3, 8, 8, 13, 4, 10000},
      {"900_100hz", "900", 7, 500e3, 7, 8, 8, 4, 10000},        {"900_50hz", "900", 8, 500e3, 7, 10, 8, 4, 20000},
      {"900_25hz", "900", 9, 500e3, 7, 10, 8, 2, 40000},
  };
  return t;
}

struct ElrsDomain {
  const char* name;
  double f0, f1;
  int n;
};
inline const std::vector<ElrsDomain>& elrs_domains() {
  static const std::vector<ElrsDomain> t = {
      {"ISM2G4", 2400.4e6, 2479.4e6, 80}, {"FCC915", 903.5e6, 926.9e6, 40}, {"AU915", 915.5e6, 926.9e6, 20},
      {"EU868", 863.275e6, 869.575e6, 13}, {"IN866", 865.375e6, 866.95e6, 4}, {"AU433", 433.42e6, 434.42e6, 3},
      {"EU433", 433.1e6, 434.45e6, 3},     {"US433", 433.25e6, 438.0e6, 8},   {"US433W", 423.5e6, 438.0e6, 20},
      {"TH920", 920.5e6, 924.7e6, 8},
  };
  return t;
}

struct NominalChannel {
  std::string domain;
  int64_t index;
  double f;
};
// python _nominal_channel(f_hz): nearest ELRS channel of any domain
inline std::optional<NominalChannel> nominal_channel(double f_hz) {
  std::optional<NominalChannel> best;
  for (const auto& d : elrs_domains()) {
    const double sp = (d.f1 - d.f0) / (d.n - 1);
    const int64_t i = py_round_int((f_hz - d.f0) / sp);
    if (0 <= i && i < d.n) {
      const double fn = d.f0 + static_cast<double>(i) * sp;
      if (!best || std::fabs(f_hz - fn) < std::fabs(f_hz - best->f)) best = NominalChannel{d.name, i, fn};
    }
  }
  return best;
}

// GFSK families: label, symbol-rate range, deviation range, protocol nominal rate (0 = None)
struct FskClass {
  const char* label;
  double r0, r1, d0, d1, nominal;
};
inline const std::vector<FskClass>& fsk_classes() {
  static const std::vector<FskClass> t = {
      {"elrs_flrc_like", 560e3, 740e3, 90e3, 260e3, 650e3},    {"elrs_gfsk_like", 250e3, 350e3, 60e3, 140e3, 300e3},
      {"afhds2a_like", 420e3, 580e3, 80e3, 260e3, 500e3},      {"crossfire_fsk_like", 78e3, 92e3, 30e3, 55e3, 85.1e3},
      {"frsky_d16_like", 62e3, 80e3, 40e3, 75e3, 70.0e3},      {"sik_like", 40e3, 260e3, 15e3, 160e3, 0.0},
  };
  return t;
}

}  // namespace si::rx::rc_link::detail
