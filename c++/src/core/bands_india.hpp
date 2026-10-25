// India frequency-allocation prior, 300 MHz - 6 GHz (python: core/bands_india.py):
// NFAP-2022 allocation table, channel rasters, Wi-Fi / BLE channel centres, FPV video
// ranges. Operator names are a prior only; reports never assert an operator from
// frequency alone.
#pragma once

#include <cstdint>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "core/types.hpp"

namespace si::bands {

constexpr double MHz = 1e6;

struct Band {
  double lo, hi;                      // Hz, [lo, hi)
  std::string key;
  std::string desc;
  std::vector<std::string> classes;   // expected signal classes
};

const std::vector<Band>& bands();                       // BANDS, in table order
std::vector<const Band*> lookup(double f_hz);           // every allocation containing f
std::vector<std::string> band_keys(double f_lo, double f_hi);  // allocations overlapping (lo, hi)
std::string describe(double f_hz);                      // "; "-joined descriptions (identical strings)

std::optional<int> gsm_arfcn(double f_dl_hz);           // GSM-900 / DCS-1800 downlink ARFCN
std::optional<std::pair<int, int>> lte_earfcn_dl(double f_hz);  // (band, EARFCN)
int64_t nr_arfcn(double f_hz);                          // NR-ARFCN, global raster

const std::vector<std::pair<int, double>>& wifi_24();   // WIFI_24 {ch: centre Hz} in dict order
const std::vector<std::pair<int, double>>& wifi_5();    // WIFI_5
const std::vector<std::pair<int, double>>& ble_adv();   // BLE_ADV {37, 38, 39}
// PULSED_BANDS: allocations where radar / DME pulses are expected (classify r_pulsed; features pulse_train)
const std::set<std::string>& pulsed_bands();
bool pulsed_band(double f_lo, double f_hi, double fc);   // python pulsed_band()
std::optional<std::pair<int, double>> nearest_wifi_channel(double f_hz);  // (channel, |f - centre|)

const std::vector<std::pair<double, double>>& fpv_video_ranges();  // FPV_VIDEO_RANGES
const std::vector<std::pair<double, double>>& analog_vtx_ranges();    // ANALOG_VTX_RANGES (+ 3.3 GHz, PLAN G9)
const std::vector<std::pair<double, double>>& digital_link_ranges();  // DIGITAL_LINK_RANGES (PLAN G8)
std::optional<std::string> uav_link_band(double f_hz);  // python uav_link_band: band context outside ISM / FPV
bool cellular_at(double f_hz);                          // python cellular_at: inside an LTE / NR band key

}  // namespace si::bands
