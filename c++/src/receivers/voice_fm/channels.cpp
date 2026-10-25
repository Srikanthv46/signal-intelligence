#include "receivers/voice_fm/channels.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace si::rx::voice_fm {

const std::vector<Plan>& plans() {
  static const std::vector<Plan> p = [] {
    std::vector<Plan> v;
    Plan pmr{"PMR446", {}, 12.5e3, 12.5e3};
    for (int n = 1; n < 17; ++n) pmr.channels.emplace_back(n, 446.00625e6 + 12.5e3 * (n - 1));
    Plan lpd{"LPD433", {}, 25e3, 25e3};
    for (int n = 1; n < 70; ++n) lpd.channels.emplace_back(n, 433.075e6 + 25e3 * (n - 1));
    Plan frs{"FRS/GMRS", {}, 12.5e3, 12.5e3};
    const double f[] = {462.5625, 462.5875, 462.6125, 462.6375, 462.6625, 462.6875, 462.7125, 467.5625,
                        467.5875, 467.6125, 467.6375, 467.6625, 467.6875, 467.7125, 462.5500, 462.5750,
                        462.6000, 462.6250, 462.6500, 462.6750, 462.7000, 462.7250};
    for (int n = 1; n < 23; ++n) frs.channels.emplace_back(n, f[n - 1] * 1e6);
    v.push_back(pmr);
    v.push_back(lpd);
    v.push_back(frs);
    return v;
  }();
  return p;
}

const std::vector<double>& ctcss_tones() {
  static const std::vector<double> t = {
      67.0,  69.3,  71.9,  74.4,  77.0,  79.7,  82.5,  85.4,  88.5,  91.5,  94.8,  97.4,  100.0,
      103.5, 107.2, 110.9, 114.8, 118.8, 123.0, 127.3, 131.8, 136.5, 141.3, 146.2, 150.0, 151.4,
      156.7, 159.8, 162.2, 165.5, 167.9, 171.3, 173.8, 177.3, 179.9, 183.5, 186.2, 189.9, 192.8,
      196.6, 199.5, 203.5, 206.5, 210.7, 218.1, 225.7, 229.1, 233.6, 241.8, 250.3, 254.1};
  return t;
}

const std::vector<std::string>& dcs_codes() {
  static const std::vector<std::string> c = [] {
    std::istringstream s(
        "023 025 026 031 032 036 043 047 051 053 054 065 071 072 073 074 114 115 116 122 125 131 132 134 "
        "143 145 152 155 156 162 165 172 174 205 212 223 225 226 243 244 245 246 251 252 255 261 263 265 266 271 274 306 "
        "311 315 325 331 332 343 346 351 356 364 365 371 411 412 413 423 431 432 445 446 452 454 455 462 464 465 466 503 "
        "506 516 523 526 532 546 565 606 612 624 627 631 632 654 662 664 703 712 723 731 732 734 743 754");
    std::vector<std::string> v;
    std::string w;
    while (s >> w) v.push_back(w);
    return v;
  }();
  return c;
}

bool is_standard_dcs(const std::string& code) {
  const auto& c = dcs_codes();
  return std::find(c.begin(), c.end(), code) != c.end();
}

std::optional<Json> match_channel(double f_hz, double ppm) {
  const double f = f_hz * (1 + ppm * 1e-6);
  const Plan* bp = nullptr;
  int bn = 0;
  double bnom = 0;
  for (const Plan& p : plans()) {
    const double tol = std::min(0.4 * p.spacing_hz, 12e-6 * f);
    for (const auto& [n, nom] : p.channels)
      if (std::fabs(f - nom) <= tol && (!bp || std::fabs(f - nom) < std::fabs(f - bnom))) {
        bp = &p;
        bn = n;
        bnom = nom;
      }
  }
  if (!bp) return std::nullopt;
  Json j = Json::object();
  j["service"] = bp->service;
  j["channel"] = bn;
  j["nominal_mhz"] = py_round(bnom / 1e6, 5);
  j["channel_bw_khz"] = bp->bw_hz / 1e3;
  j["offset_hz"] = py_round_int(f_hz - bnom);
  j["offset_ppm"] = py_round((f_hz - bnom) / bnom * 1e6, 2);
  j["receiver_ppm_applied"] = ppm;
  return j;
}

}  // namespace si::rx::voice_fm
