#include "core/bands_india.hpp"

#include <cctype>
#include <cmath>

namespace si::bands {

// Generated from python/spectrum_intel/core/bands_india.py BANDS (doubles printed with
// repr, so every edge is the exact value Python computed, e.g. 328.6 * MHz).
const std::vector<Band>& bands() {
  static const std::vector<Band> t = {
    {70000000.0, 74800000.0, "govt_mobile_70", "Fixed/mobile (government, defence)", {"narrowband_fm", "pmr"}},
    {74800000.0, 75200000.0, "marker_75", "Aeronautical marker beacons 75 MHz (ILS outer / middle / inner)", {"marker_beacon"}},
    {75200000.0, 87500000.0, "govt_mobile_75", "Fixed/mobile (government)", {"narrowband_fm", "pmr"}},
    {87500000.0, 108000000.0, "fm_bcast", "FM broadcasting (AIR / private FM), RDS", {"fm_broadcast", "rds"}},
    {108000000.0, 111975000.0, "vor_ils_loc", "VOR (even tenths) and ILS localizer (odd tenths)", {"vor", "ils"}},
    {111975000.0, 117975000.0, "vor", "VOR (and GBAS VDB)", {"vor"}},
    {117975000.0, 137000000.0, "air_vhf", "Aeronautical mobile: ATC AM voice 118-136.975, 121.5 emergency, ACARS, VDL Mode 2 136.725-136.975", {"am_voice", "acars", "vdl2"}},
    {137000000.0, 138000000.0, "metsat_137", "Meteorological satellite downlinks (Meteor-M LRPT)", {"metsat"}},
    {138000000.0, 144000000.0, "govt_mobile_138", "Fixed/mobile (defence)", {"narrowband_fm", "pmr"}},
    {144000000.0, 146000000.0, "ham_2m", "Amateur 2 m: FM repeaters, APRS, amateur satellites", {"aprs", "narrowband_fm"}},
    {146000000.0, 156000000.0, "govt_mobile_146", "Land mobile: police, state government, railways (FM, DMR)", {"narrowband_fm", "dmr", "pmr"}},
    {156000000.0, 162050000.0, "marine_vhf", "Maritime VHF: ch 16 156.8, DSC ch 70 156.525, AIS 161.975 / 162.025", {"marine_fm", "dsc", "ais"}},
    {162050000.0, 174000000.0, "govt_mobile_162", "Land mobile (police, railways, forest, PSU)", {"narrowband_fm", "dmr", "pmr"}},
    {174000000.0, 225000000.0, "vhf_band3", "Broadcasting Band III (VHF TV largely switched off; no DAB)", {"analog_tv"}},
    {225000000.0, 300000000.0, "mil_air_uhf", "Military aeronautical UHF (AM voice 225-400, 243.0 guard); fixed / mobile (defence)", {"am_voice"}},
    {300000000.0, 328600000.0, "govt_mobile_300", "Fixed/mobile (government, defence)", {"narrowband_fm", "pmr"}},
    {328600000.0, 335400000.0, "ils_glidepath", "ILS glide path (aeronautical)", {"ils"}},
    {335400000.0, 380000000.0, "govt_mobile_335", "Fixed/mobile (government, defence)", {"narrowband_fm", "pmr"}},
    {380000000.0, 400000000.0, "tetra_ppdr", "TETRA / PPDR public-safety trunking", {"tetra", "pmr"}},
    {400050000.0, 406000000.0, "metaids", "Meteorological aids (radiosondes) / Met-sat", {"radiosonde"}},
    {406000000.0, 406100000.0, "cospas_sarsat", "COSPAS-SARSAT distress beacons (EPIRB/PLB)", {"epirb"}},
    {406100000.0, 433050000.0, "govt_mobile_406", "Fixed/mobile (government), PMR", {"pmr", "narrowband_fm"}},
    {433050000.0, 434790000.0, "ism_433", "433 MHz SRD (delicensed, GSR 1048(E))", {"ism_srd", "lora"}},
    {434790000.0, 450000000.0, "govt_mobile_435", "Fixed/mobile (government), amateur 430-440", {"pmr", "narrowband_fm"}},
    {450000000.0, 470000000.0, "pmr_450", "PMR / land-mobile (police, utilities), DMR/analog FM", {"pmr", "dmr", "narrowband_fm"}},
    {470000000.0, 585000000.0, "uhf_tv", "UHF TV broadcasting (Doordarshan DVB-T2 / analog)", {"dvb_t2", "analog_tv"}},
    {585000000.0, 698000000.0, "uhf_bcast_imt", "Broadcasting / IMT (600 MHz, unassigned)", {"dvb_t2"}},
    {703000000.0, 748000000.0, "b28_ul", "LTE B28 / NR n28 uplink (FDD)", {"lte_ul", "nr_ul"}},
    {748000000.0, 758000000.0, "b28_gap", "B28 duplex gap", {}},
    {758000000.0, 803000000.0, "b28_dl", "LTE B28 / NR n28 downlink (FDD) — Jio 5G SA, Airtel, Railways LTE-R", {"lte_dl", "nr_dl"}},
    {806000000.0, 824000000.0, "trunk_800", "Trunked radio / IMT 800 (B26 lower)", {"pmr", "lte_ul"}},
    {824000000.0, 849000000.0, "b5_ul", "LTE B5 / NR n5 uplink (850 MHz, FDD)", {"lte_ul"}},
    {851000000.0, 869000000.0, "trunk_850", "Trunked radio / B26 downlink lower", {"pmr", "lte_dl"}},
    {869000000.0, 894000000.0, "b5_dl", "LTE B5 / NR n5 downlink (850 MHz) — Jio, Airtel", {"lte_dl", "nr_dl"}},
    {865000000.0, 867000000.0, "ism_865", "865-867 MHz SRD (delicensed): LoRaWAN IN865, RFID", {"lora", "ism_srd", "rfid"}},
    {880000000.0, 915000000.0, "b8_ul", "GSM-900 / LTE B8 uplink (FDD)", {"gsm_ul", "lte_ul"}},
    {915000000.0, 925000000.0, "b8_gap", "900 MHz duplex gap / RFID (India 920-925? not delicensed)", {}},
    {925000000.0, 960000000.0, "b8_dl", "GSM-900 / LTE B8 / NR n8 downlink — Airtel, Vi, BSNL", {"gsm_dl", "lte_dl", "nr_dl"}},
    {960000000.0, 1164000000.0, "arns_l", "Aeronautical radionavigation: DME/TACAN 962-1213, SSR 1030/1090, ADS-B", {"dme", "ssr", "adsb"}},
    {1164000000.0, 1215000000.0, "rnss_l5", "GNSS L5/E5a/E5b, NavIC L5 1176.45 (below noise) + DME", {"gnss", "dme"}},
    {1215000000.0, 1300000000.0, "rnss_l2_radar", "GNSS L2/G2/E6 + L-band surveillance radar", {"gnss", "radar"}},
    {1300000000.0, 1350000000.0, "arns_radar", "Aeronautical / L-band radar", {"radar"}},
    {1350000000.0, 1400000000.0, "fixed_radar_1350", "Fixed / mobile / radiolocation", {"radar"}},
    {1400000000.0, 1427000000.0, "passive_1400", "Passive (radio astronomy / EESS) — no emissions allowed", {}},
    {1427000000.0, 1518000000.0, "l_band_imt", "L-band IMT / fixed links (not deployed for mobile)", {"fixed_link"}},
    {1518000000.0, 1559000000.0, "mss_dl", "MSS downlink (Inmarsat 1525-1559, Thuraya)", {"mss"}},
    {1559000000.0, 1610000000.0, "rnss_l1", "GNSS L1/E1/B1/G1 (below noise)", {"gnss"}},
    {1610000000.0, 1626500000.0, "mss_iridium", "MSS: Iridium 1616-1626.5 (TDMA bursts), Globalstar uplink", {"iridium"}},
    {1626500000.0, 1660500000.0, "mss_ul", "MSS uplink (Inmarsat terminals)", {"mss"}},
    {1660500000.0, 1710000000.0, "met_1675", "Met-sat / MetAids 1668-1710", {"metsat"}},
    {1710000000.0, 1785000000.0, "b3_ul", "LTE B3 / NR n3 uplink (1800 MHz FDD), GSM-1800 UL", {"lte_ul", "gsm_ul"}},
    {1785000000.0, 1805000000.0, "b3_gap", "1800 duplex gap", {}},
    {1805000000.0, 1880000000.0, "b3_dl", "LTE B3 / NR n3 / GSM-1800 downlink — all operators", {"lte_dl", "nr_dl", "gsm_dl"}},
    {1880000000.0, 1920000000.0, "b39", "IMT 1900 (unassigned) / DECT-like cordless", {}},
    {1920000000.0, 1980000000.0, "b1_ul", "LTE B1 / UMTS 2100 uplink", {"lte_ul", "umts_ul"}},
    {1980000000.0, 2010000000.0, "mss_s_ul", "MSS S-band uplink", {}},
    {2010000000.0, 2110000000.0, "space_2025", "Space operation / research (uplink TT&C)", {}},
    {2110000000.0, 2170000000.0, "b1_dl", "LTE B1 / NR n1 / UMTS 2100 downlink — Airtel, Vi, BSNL", {"lte_dl", "nr_dl", "umts_dl"}},
    {2170000000.0, 2200000000.0, "mss_s_dl", "MSS S-band downlink", {}},
    {2200000000.0, 2300000000.0, "space_tlm", "Space research / telemetry (ISRO S-band TLM)", {"telemetry"}},
    {2300000000.0, 2400000000.0, "b40", "LTE B40 / NR n40 TDD (2300 MHz) — Jio, Airtel", {"lte_tdd", "nr_tdd"}},
    {2400000000.0, 2483500000.0, "ism_24", "2.4 GHz ISM (delicensed): Wi-Fi, BLE/BT, Zigbee, drones, microwave ovens", {"wifi", "wifi_b", "ble", "bt_classic", "zigbee", "dji", "fhss", "lora", "microwave_oven", "analog_video"}},
    {2483500000.0, 2500000000.0, "mss_24", "MSS downlink (Globalstar) / NavIC S-band 2492.028", {}},
    {2500000000.0, 2690000000.0, "b41", "LTE B41 / NR n41 TDD (2500 MHz) — Vi, BSNL / BSS", {"lte_tdd", "nr_tdd"}},
    {2700000000.0, 2900000000.0, "radar_s_atc", "S-band ATC / weather radar", {"radar"}},
    {2900000000.0, 3100000000.0, "radar_s_marine", "Maritime / radiolocation radar", {"radar"}},
    {3100000000.0, 3300000000.0, "radar_s_def", "Radiolocation (defence radar)", {"radar"}},
    {3300000000.0, 3670000000.0, "n78", "NR n78 TDD (3.3-3.67 GHz) — Jio, Airtel, Vi 5G", {"nr_tdd"}},
    {3670000000.0, 3700000000.0, "c_band_low", "C-band (guard / FSS)", {}},
    {3700000000.0, 4200000000.0, "fss_c", "C-band FSS satellite downlink (satellite TV; dish needed)", {"fss"}},
    {4200000000.0, 4400000000.0, "radalt", "Aircraft radio altimeters (FMCW)", {"radar_altimeter"}},
    {4400000000.0, 4990000000.0, "govt_4400", "Fixed/mobile (defence); 4.8-4.99 IMT-identified", {"fixed_link"}},
    {4990000000.0, 5000000000.0, "passive_4990", "Radio astronomy", {}},
    {5000000000.0, 5150000000.0, "arns_mls", "ARNS / MLS 5030-5091 / FSS", {}},
    {5150000000.0, 5350000000.0, "unii_12", "5 GHz Wi-Fi (U-NII-1/2, delicensed)", {"wifi", "dji", "analog_video"}},
    {5350000000.0, 5470000000.0, "radar_5400", "Radiolocation (no Wi-Fi)", {"radar"}},
    {5470000000.0, 5725000000.0, "unii_2e", "5 GHz Wi-Fi (U-NII-2e, delicensed) / weather radar 5600-5650", {"wifi", "radar"}},
    {5725000000.0, 5875000000.0, "ism_58", "5.8 GHz ISM / Wi-Fi U-NII-3 (delicensed): Wi-Fi, FPV video, DJI", {"wifi", "dji", "analog_video", "digital_vtx"}},
    {5875000000.0, 5925000000.0, "its_59", "ITS / 5.9 GHz", {"wifi"}},
    {5925000000.0, 6425000000.0, "wifi6e", "Lower 6 GHz (Wi-Fi 6E delicensing in progress in India)", {"wifi"}},
  };
  return t;
}

std::vector<const Band*> lookup(double f_hz) {
  std::vector<const Band*> out;
  for (const Band& b : bands())
    if (b.lo <= f_hz && f_hz < b.hi) out.push_back(&b);
  return out;
}

std::vector<std::string> band_keys(double f_lo, double f_hi) {
  std::vector<std::string> out;
  for (const Band& b : bands())
    if (b.lo < f_hi && b.hi > f_lo) out.push_back(b.key);
  return out;
}

std::string describe(double f_hz) {
  const auto m = lookup(f_hz);
  if (m.empty()) return "no specific allocation in table";
  std::string s;
  for (size_t i = 0; i < m.size(); ++i) {
    if (i) s += "; ";
    s += m[i]->desc;
  }
  return s;
}

std::optional<int> gsm_arfcn(double f_dl_hz) {
  const double f = f_dl_hz / MHz;
  if (935 <= f && f <= 960) return static_cast<int>(py_round_int((f - 935) / 0.2));
  if (1805.2 <= f && f <= 1879.8) return static_cast<int>(py_round_int((f - 1805.2) / 0.2)) + 512;
  return std::nullopt;
}

std::optional<std::pair<int, int>> lte_earfcn_dl(double f_hz) {
  struct E { int band; double flo; int noff; double span; };
  static const E tab[] = {{1, 2110, 0, 60}, {3, 1805, 1200, 75}, {5, 869, 2400, 25}, {8, 925, 3450, 35},
                          {28, 758, 9210, 45}, {40, 2300, 38650, 100}, {41, 2496, 39650, 194}};
  const double f = f_hz / MHz;
  for (const E& e : tab)
    if (e.flo <= f && f <= e.flo + e.span)
      return std::make_pair(e.band, static_cast<int>(py_round_int((f - e.flo) * 10)) + e.noff);
  return std::nullopt;
}

int64_t nr_arfcn(double f_hz) {
  const double f = f_hz / 1e3;
  if (f < 3e6) return py_round_int(f / 5);
  return py_round_int(600000 + (f - 3e6) / 15);
}

const std::vector<std::pair<int, double>>& wifi_24() {
  static const std::vector<std::pair<int, double>> t = [] {
    std::vector<std::pair<int, double>> v;
    for (int ch = 1; ch < 14; ++ch) v.emplace_back(ch, 2407e6 + 5e6 * ch);
    v.emplace_back(14, 2484e6);
    return v;
  }();
  return t;
}

const std::vector<std::pair<int, double>>& wifi_5() {
  static const std::vector<std::pair<int, double>> t = [] {
    std::vector<std::pair<int, double>> v;
    for (int ch = 32; ch < 178; ch += 2) v.emplace_back(ch, 5000e6 + 5e6 * ch);
    return v;
  }();
  return t;
}

const std::vector<std::pair<int, double>>& ble_adv() {
  static const std::vector<std::pair<int, double>> t = {{37, 2402e6}, {38, 2426e6}, {39, 2480e6}};
  return t;
}

const std::set<std::string>& pulsed_bands() {
  static const std::set<std::string> t = {"arns_l",       "rnss_l5",       "rnss_l2_radar",  "arns_radar",
                                          "fixed_radar_1350", "radar_s_atc", "radar_s_marine", "radar_s_def",
                                          "radar_5400",   "unii_2e",       "radalt"};
  return t;
}

bool pulsed_band(double f_lo, double f_hi, double fc) {
  const auto& pb = pulsed_bands();
  for (const auto& k : band_keys(f_lo, f_hi))
    if (pb.count(k)) return true;
  for (const auto* b : lookup(fc))
    if (pb.count(b->key)) return true;
  return false;
}

std::optional<std::pair<int, double>> nearest_wifi_channel(double f_hz) {
  std::optional<std::pair<int, double>> best;
  for (const auto* tab : {&wifi_24(), &wifi_5()})
    for (const auto& [ch, fc] : *tab) {
      const double d = std::fabs(f_hz - fc);
      if (!best || d < best->second) best = std::make_pair(ch, d);
    }
  return best;
}

const std::vector<std::pair<double, double>>& fpv_video_ranges() {
  static const std::vector<std::pair<double, double>> t = {{1060000000.0, 1380000000.0}, {2330000000.0, 2530000000.0}, {4900000000.0, 6000000000.0}};
  return t;
}

const std::vector<std::pair<double, double>>& analog_vtx_ranges() {
  static const std::vector<std::pair<double, double>> t = {
      {1060000000.0, 1380000000.0}, {2330000000.0, 2530000000.0}, {3150000000.0, 3490000000.0}, {4900000000.0, 6000000000.0}};
  return t;
}

const std::vector<std::pair<double, double>>& digital_link_ranges() {
  static const std::vector<std::pair<double, double>> t = {
      {1350000000.0, 1444000000.0}, {1600000000.0, 2530000000.0}, {4400000000.0, 6000000000.0}};
  return t;
}

std::optional<std::string> uav_link_band(double f_hz) {
  struct B { double lo, hi; const char* text; };
  static const B tab[] = {
      {1430, 1444, "1.4 GHz UAV data-link band (China 1430-1444 MHz; Taisync-like)"},
      {1350, 1444, "L-band 1.35-1.44 GHz (government / professional UAV links)"},
      {1600, 2330, "1.6-2.33 GHz government / professional band (e.g. Skydio Connect MH, Microhard 1.8 GHz, "
                   "Silvus / Doodle Labs)"},
      {4400, 4900, "4.4-4.9 GHz government / professional band (e.g. Silvus, Doodle Labs, Persistent MPU5)"},
  };
  const double f = f_hz / 1e6;
  for (const auto& b : tab)
    if (b.lo <= f && f < b.hi) return std::string(b.text);
  return std::nullopt;
}

bool cellular_at(double f_hz) {
  for (const auto& k : band_keys(f_hz, f_hz))
    if (k.size() >= 2 && (k[0] == 'b' || k[0] == 'n') && std::isdigit(static_cast<unsigned char>(k[1]))) return true;
  return false;
}

}  // namespace si::bands
