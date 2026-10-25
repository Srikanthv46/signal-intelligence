// GFSK / FSK narrowband bursts, BLE advertising and 2 Mb/s ESB-like packets (python:
// receivers/rc_link/detector.py, "FSK / GFSK bursts": fsk_scan, _symbol_rate, BLE helpers, fsk_refine,
// _fsk_measure, _gsm_burst, _fsk_classify, _esb_try, detect_fsk; esb.py).
#pragma once

#include <optional>
#include <string>
#include <vector>

#include "core/model.hpp"
#include "core/types.hpp"

namespace si::rx::rc_link {

// python detect_fsk(x, X, fs, fc) (X is not used by Python either)
std::vector<Emission> detect_fsk(const cfvec& x, double fs, double fc);

namespace detail {

struct FskCand {
  double t0_us = 0, t1_us = 0, f_lo = 0, f_hi = 0, fc_b = 0, zin = 0, peak_e = 0;
};
// python fsk_scan(x, fs, fc, nfft=512, k_sigma=8.0, widths=(1, 3, 7), t_frames=12, max_bw=3.0e6,
// min_us=30.0, edge=0.45)
std::vector<FskCand> fsk_scan(const cfvec& x, double fs, double fc);

// python _symbol_rate(ifr, ofs) -> (rate or None, concentration, phase)
struct SymRate {
  std::optional<double> rate;
  double conc = 0, ph = 0;
};
SymRate symbol_rate(const fvec& ifr, double ofs, double r_lo = 20e3, double r_hi = 1.3e6);

int ble_channel_index(double f_hz);  // -1 = None
// python ble_try_decode(y, ofs, f_hz): the dict, or null
Json ble_try_decode(const cfvec& y, double ofs, double f_hz);
// python ble_coherent_decode(y, ofs, ch, aa_bits): the low-SNR path alone (ch = whitening channel, -1 = unknown)
Json ble_coherent_decode(const cfvec& y, double ofs, int ch);
// python ble_coded_decode(y, ofs, ch, aa_bits): LE Coded PHY (BT5 Long Range) path
Json ble_coded_decode(const cfvec& y, double ofs, int ch);

// one measured burst (python dict returned by _fsk_measure)
struct FskMeas {
  double t0_us = 0, t1_us = 0, f_center = 0, f_lo = 0, f_hi = 0, rho = 0, ofs = 0, n0 = 0, bw99 = 0, snr_db = 0;
  std::optional<bool> const_env;
  bool has_rate = false;            // sym_rate / sym_conc keys present
  std::optional<double> sym_rate;
  double sym_conc = 0;
  bool has_dev = false;             // dev_hz / fsk_center / mod_index keys present
  double dev_hz = 0, fsk_center = 0, mod_index = 0;
  Json ble = nullptr;
  Json sik = nullptr;               // python r["sik"] (SiK telemetry decode, research)
  Json crossfire = nullptr;         // python r["crossfire"] (TBS Crossfire sync + packet, research)
  Json frsky = nullptr;             // python r["frsky"] (FrSky D8 / D16 identification, research)
  Json flarm = nullptr;             // python r["flarm"] (FLARM clear header, research)
  Json ieee802154 = nullptr;        // python r["ieee802154"] (802.15.4 frame, research)
  Json xn297 = nullptr;             // python r["xn297"] (XN297 / Bayang, research)
  Json bredr = nullptr;             // python r["bredr"] (Bluetooth BR/EDR access code, research S7)
  Json afhds2a = nullptr;           // python r["afhds2a"] (FlySky AFHDS 2A identification, research)
  Json frame = nullptr;             // python r["frame"] (blind FSK frame structure, research)
  Json radiolink = nullptr;         // python r["radiolink"] (RadioLink CC2500, research, PLAN G4)
  Json hubsan = nullptr;            // python r["hubsan"] (Hubsan A7105, research, PLAN G5)
  Json elrs_gfsk = nullptr;         // python r["elrs_gfsk"] (ExpressLRS LR1121 GFSK SYNC, research, PLAN G3)
  Json cc2500_rc = nullptr;         // python r["cc2500_rc"] (HoTT / Hitec / S-FHSS, research, PLAN G12)
  Json adsl = nullptr;              // python r["adsl"] (ADS-L / OGN tracker on the FLARM PHY, research, PLAN G10)
  bool esb = false;                 // python r["esb"] set (a 2 Mb/s ESB-like packet was sliced)
};
// python _refine_local(x, fs, fc, cand): fsk_refine on a local FFT, capture-relative times
std::vector<FskMeas> refine_local(const cfvec& x, double fs, double fc, const FskCand& cand);

struct FskLabel {
  std::optional<std::string> family, label;
  double conf = 0;
  std::optional<double> nominal;
};
FskLabel fsk_classify(const FskMeas& r, double fc);

}  // namespace detail
}  // namespace si::rx::rc_link
