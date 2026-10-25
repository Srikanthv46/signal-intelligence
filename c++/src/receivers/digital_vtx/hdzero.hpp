// HDZero (Divimath DM6300) digital FPV video link, 27 MHz mode (python: receivers/digital_vtx/hdzero.py).
// Measured on an HDZero Freestyle V2 (testdata/live/LIVE_HDZERO_5658MHz): CP-OFDM, 1024 FFT at 27 MSps, CP 32,
// 889 used subcarriers, 26 continual pilots, a per-frame preamble block. A link is accepted on the continual-pilot
// coherence against the data carriers beside them.
#pragma once

#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "core/types.hpp"

namespace si::rx::digital_vtx::hdzero {

constexpr double FS = 27e6;
constexpr int64_t N = 1024, CP = 32, TS = N + CP, HALF = 444;
constexpr double SCS = FS / N;
constexpr int N_CONTINUAL = 26;
extern const int CONTINUAL[N_CONTINUAL];

struct Link {
  std::string channel;
  double channel_mhz = 0, off = 0, cfo_hz = 0, ppm = 0, lo = 0, hi = 0;
  double pilot_coherence = 0, data_coherence = 0, fold_contrast = 0, p_band = 0, half_bw = 0;
  int64_t n_symbols = 0;
  dvec preamble_ms;
  std::vector<int64_t> preamble_symbols;
  std::optional<int64_t> frame_symbols;
  std::optional<double> frame_ms;
  double t0 = 0, t1 = 0;
};

struct LogEntry {
  std::string channel;
  double fold_contrast = 0, pilot_coherence = 0, data_coherence = 0;
  int64_t shift = 0;
  bool accepted = false;
};

// python slice_plan(n_samples, fs) -> (slice length, per-symbol FFT size)
std::pair<int64_t, int64_t> slice_plan(int64_t n_samples, double fs);
bool smooth(int64_t n);

// sl(f_off, half_bw, out_fs, &rate) -> slice; python hdzero_links(sl, fs, fc, usable, n_samples, log)
using Slicer = std::function<cdvec(double, double, double, double*)>;
std::vector<Link> hdzero_links(const Slicer& sl, double fs, double fc, double usable, int64_t n_samples,
                               std::vector<LogEntry>* log = nullptr);

}  // namespace si::rx::digital_vtx::hdzero
