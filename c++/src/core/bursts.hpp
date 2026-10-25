// Occupied band and burst extent of a KNOWN link from its in-band power, whatever the waveform
// (python: core/bursts.py). A receiver finds a link by its structure (CP, sync sequence); in-band
// power over many subcarriers is far more sensitive per symbol but blind to identity, so power gives
// the extent and the receiver's own evidence (seed runs) gives the membership.
#pragma once

#include <optional>
#include <utility>
#include <vector>

#include "core/types.hpp"

namespace si::bursts {

using Runs = std::vector<std::pair<double, double>>;  // [(t0, t1)] seconds

constexpr double kFrameS = 16.67e-6;
constexpr int kMaxSb = 8;
constexpr double kSbHz = 2e6;
constexpr double kGuardGapHz = 0.3e6;
constexpr double kGuardHz = 1.0e6;
constexpr double kMaxBridgeS = 0.6e-3;
constexpr size_t kSmooth = 3;

// python TFPower: |fftshift(fft(Hann-windowed frames))|^2, L = 2^round(log2(fs * frame_s)) (>= 16),
// no overlap; f: bin offsets from fc (Hz).
class TFPower {
 public:
  TFPower(const cfvec& x, double fs, double frame_s = kFrameS);
  std::vector<size_t> cols(double lo, double hi) const;
  std::vector<uint8_t> frames(const Runs& runs) const;            // frame entirely inside a run
  std::vector<uint8_t> touched(const Runs& runs, int64_t pad = 0) const;  // frame overlaps a run (+-pad)

  double fs = 0, dt = 0;
  size_t L = 0, n = 0;
  dmat P;
  dvec f;
};

struct RectFit {
  size_t lo = 0, hi = 0;
  double inside = 0, outside = 0;
};
// python rect_fit(d)
RectFit rect_fit(const dvec& d);

struct Band {
  double lo = 0, hi = 0, contrast_db = 0, z = 0;
};
// python fit_band(tf, runs, lo, hi, search_lo, search_hi, smooth_hz=0.25e6); nullopt = None
std::optional<Band> fit_band(const TFPower& tf, const Runs& runs, double lo, double hi, double search_lo,
                             double search_hi, double smooth_hz = 0.25e6);

struct EnergyInfo {
  bool used = false;
  bool has_level = false;
  double level_db = 0, sigma = 0, thr = 0;
  size_t n_energy = 0;
};
// python energy_bursts(tf, seeds, lo, hi, usable_half, min_sym_s=0.0)
std::pair<Runs, EnergyInfo> energy_bursts(const TFPower& tf, Runs seeds, double lo, double hi, double usable_half,
                                          double min_sym_s = 0.0);

}  // namespace si::bursts
