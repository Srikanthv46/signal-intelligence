// 802.11 OFDM at a non-standard clock, at any frequency, in either spectrum orientation
// (python receivers/wifi/ofdm_scaled.py): bursty wideband bands are found on a ~60 kHz
// spectrogram, their OFDM numerology (useful-symbol time from the cyclic-prefix correlation,
// symbol period from its periodicity) gives the 802.11 clock, and the capture is decoded by
// ofdm::WifiDetector as if sampled at fs * 20 MHz / clock, plain and mirrored.
#pragma once

#include <optional>
#include <string>
#include <vector>

#include "core/types.hpp"
#include "receivers/wifi/ofdm_rx.hpp"

namespace si::rx::wifi::scaled {

// python find_bands() band dict
struct Band {
  double f_lo = 0, f_hi = 0;
  int64_t n_frames = 0;
  std::vector<uint8_t> frames;
  size_t frame_len = 0;
  bool merged = false;
};
// python find_bands(x, fs, fc) -> (bands, unions)
std::pair<std::vector<Band>, std::vector<Band>> find_bands(const cfvec& x, double fs, double fc, double thresh_db = 10.0);

// python numerology() dict; `json` holds it in Python's key order
struct Numerology {
  double center_hz = 0, occupied_hz = 0, tu_s = 0, ts_s = 0, guard = 0, cp_corr = 0, cp_corr_over_median = 0,
         symbol_line_over_median = 0;
  std::optional<double> clock_hz, clock_measured_hz, occupied_over_clock;
  std::optional<std::string> reject;
  Json json() const;
};
std::optional<Numerology> numerology(const cfvec& x, double fs, double fc, const Band& band);

// python Link
struct Link {
  double center_hz = 0, clock_hz = 0;
  bool inverted = false;
  Numerology numerology;
  std::vector<ofdm::Packet> packets;  // real units
  size_t n_rejects = 0;
};
// python scan(x, fs, fc)
std::vector<Link> scan(const cfvec& x, double fs, double fc, size_t max_bands = 12, size_t max_unions = 4);

}  // namespace si::rx::wifi::scaled
