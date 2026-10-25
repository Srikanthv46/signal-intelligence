// python receivers/rc_link/protocols.py (same table, order and strings).
#include "receivers/rc_link/protocols.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include "core/pyfmt.hpp"

namespace si::rx::rc_link::protocols {

namespace {
using Bands = std::vector<std::pair<double, double>>;
struct Proto { const char* name; double R, dev; Bands bands; double per; const char* note; };  // dev/per 0 = None
const Bands kS24 = {{2400.0, 2484.0}};
const Bands kSub = {{420.0, 480.0}, {860.0, 930.0}};

const std::vector<Proto>& table() {
  static const std::vector<Proto> p = {
      {"FrSky D8 (ACCST, CC2500)", 31.0e3, 31.7e3, kS24, 9.0, "47 ch on a 299.927 kHz grid"},
      {"RadioLink (CC2500)", 38.0e3, 38.0e3, kS24, 20.0, "16 ch, 400 kHz grid from 2402 MHz"},
      {"FrSky D16 / ACCST v1 (FCC)", 70.0e3, 57.0e3, kS24, 9.0, "47 ch, 299.927 kHz grid"},
      {"FrSky ACCST v2 (D16 v2)", 77.0e3, 57.0e3, kS24, 9.0, "hop table from the 16-bit TX ID"},
      {"FrSky D16 LBT (EU)", 100.0e3, 70.0e3, kS24, 9.0, "2-FSK, base 2405 MHz"},
      {"TBS Crossfire 150 Hz (FSK)", 85.1e3, 42.48e3, {{860.0, 930.0}}, 6.666, "260 kHz grid, downlink +50 channels"},
      {"mLRS 50 Hz FSK", 100.0e3, 50.0e3, kSub, 20.0, "91-byte frames, whitening"},
      {"FLARM-compatible (2-FSK, Manchester)", 100.0e3, 50.0e3, {{860.0, 870.0}, {902.0, 928.0}}, 0,
       "India 866.0 MHz; ~6 ms every 0.6-1.4 s; aircraft type 13 = UAV"},
      {"Futaba S-FHSS (CC2500)", 128.1e3, 38.1e3, kS24, 6.8, "30 ch, 1.5 MHz hop grid, packets in pairs"},
      {"Graupner HoTT (MSK, CC2500)", 250.0e3, 62.5e3, kS24, 10.0, "MSK + FEC, 338 kHz grid"},
      {"Hitec AFHSS / Minima (MSK, CC2500)", 250.0e3, 62.5e3, kS24, 22.5, "sync word = TX ID, 300 kHz grid"},
      {"ExpressLRS LR1121 GFSK", 300.0e3, 100.0e3, {{2400.0, 2484.0}, {860.0, 930.0}}, 1.0, "1000 Hz / DVDA"},
      {"FlySky AFHDS 2A (A7105)", 500.0e3, 0, kS24, 3.85, "16 ch from the TX ID, 0.5 MHz grid"},
      {"ExpressLRS / mLRS / TBS Tracer FLRC (SX1280)", 650.0e3, 0, kS24, 0, "ELRS FLRC sync word = UID[2..5]"},
      {"nRF24 / XN297 toy RC (Bayang, V2x2, ...)", 1.0e6, 160.0e3, kS24, 0, "Bayang 2 ms, V2x2 4 ms"},
      {"OpenLRS-ng 57.6 kb/s", 57.6e3, 28.75e3, kSub, 0, "Si443x, rf_magic bind ID"},
      {"OpenLRS-ng 125 kb/s", 125.0e3, 125.0e3, kSub, 0, "Si443x"},
      {"XBee-PRO 900HP 200 kb/s", 200.0e3, 0, {{902.0, 928.0}}, 0, "FSK FHSS, 64 ch"},
  };
  return p;
}
const std::vector<std::pair<int, double>> kSik = {{2, 1.875},  {4, 3.75},   {8, 8.125},   {16, 16.25},   {19, 18.75},
                                                  {24, 23.75}, {32, 31.875}, {48, 48.125}, {64, 63.75}, {96, 96.25},
                                                  {128, 128.125}, {192, 158.75}, {250, 158.75}};

bool in_bands(double f, const Bands& b) {
  for (const auto& [lo, hi] : b)
    if (lo <= f && f <= hi) return true;
  return false;
}
}  // namespace

Json candidates(double sym_rate, std::optional<double> dev_hz, std::optional<double> f_mhz,
                std::optional<double> interval_ms, size_t top) {
  struct C { Json j; int score; size_t nc; };
  std::vector<C> out;
  if (!(sym_rate > 0)) return Json::array();
  auto add = [&](const std::string& name, const Json& m, const Json& c, const char* note) {
    const int score = static_cast<int>(m.size()) - static_cast<int>(c.size());
    out.push_back({Json{{"protocol", name}, {"score", score}, {"matched", m}, {"conflicts", c}, {"note", note}}, score, c.size()});
  };
  const bool have_dev = dev_hz && *dev_hz != 0;
  for (const auto& p : table()) {
    if (std::fabs(sym_rate - p.R) > 0.05 * p.R) continue;
    Json m = Json::array(), c = Json::array();
    m.push_back(fmt("rate %g kb/s", p.R / 1e3));
    if (p.dev != 0 && have_dev) (std::fabs(*dev_hz - p.dev) <= 0.25 * p.dev ? m : c).push_back(fmt("deviation %g kHz", p.dev / 1e3));
    if (f_mhz) (in_bands(*f_mhz, p.bands) ? m : c).push_back("band");
    if (p.per != 0 && interval_ms)
      (std::fabs(*interval_ms - p.per) <= 0.05 * p.per ? m : c).push_back(fmt("packet period %g ms", p.per));
    add(p.name, m, c, p.note);
  }
  static const Bands sik_bands = {{433.0, 435.0}, {470.0, 471.0}, {868.0, 870.0}, {902.0, 928.0}};
  for (const auto& [kb, dk] : kSik) {
    if (std::fabs(sym_rate - kb * 1e3) > 0.05 * kb * 1e3) continue;
    Json m = Json::array(), c = Json::array();
    m.push_back(fmt("rate %d kb/s", kb));
    if (have_dev) (std::fabs(*dev_hz - dk * 1e3) <= 0.25 * dk * 1e3 ? m : c).push_back(fmt("deviation %g kHz (h ~ 2)", dk));
    if (f_mhz) (in_bands(*f_mhz, sik_bands) ? m : c).push_back("band");
    add(fmt("SiK telemetry radio %d kb/s (3DR / Holybro / mRo / RFD900 compatible)", kb), m, c,
        "GFSK, PN9 whitening, NETID header; TDM window ~114 ms at 64 kb/s");
  }
  std::stable_sort(out.begin(), out.end(), [](const C& a, const C& b) { return a.score != b.score ? a.score > b.score : a.nc < b.nc; });
  Json res = Json::array();
  for (size_t i = 0; i < out.size() && i < top; ++i) res.push_back(out[i].j);
  return res;
}

}  // namespace si::rx::rc_link::protocols
