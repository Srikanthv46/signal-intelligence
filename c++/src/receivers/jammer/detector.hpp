// Swept-carrier jammer detector (python: receivers/jammer/detector.py).
//
//   detect(x, fs, fc) -> Emission family "jammer", label "swept_carrier" (at most one per block)
//   measure(x, fs, fc) -> Sweep (found = false when none), with the per-candidate diagnostics
//
// One carrier swept periodically across a band lights every STFT bin it crosses once or twice per sweep: a
// cyclic line in the per-bin occupancy gives the rate, a fold of the STFT on the period shows one carrier
// moving continuously. Depends on core/ only.
#pragma once

#include <string>
#include <utility>
#include <vector>

#include "core/model.hpp"
#include "core/types.hpp"

namespace si::rx::jammer {

constexpr size_t kNfft = 64;            // STFT size: fs/64 bins, hop kNfft/2
constexpr double kTrigDb = 10.0;        // cyclic line over its local median
constexpr double kHarmDb = 10.0;        // a harmonic counts for the rate refinement
constexpr size_t kNCand = 3;            // cyclic-line candidates tried on the fold
constexpr int kMaxDiv = 8;              // the true period is up to this many line periods
constexpr double kContrastFrac = 0.8;   // the period: shortest fold with this fraction of the best structure
constexpr double kMaxPhases = 128;
constexpr double kNSeg = 4;             // block segments for the rate (drift) fit, at most ...
constexpr double kSegPeriods = 20;      // ... and at least this many periods each (>= 2 segments)
constexpr double kActiveDb = 3.0;       // folded phase holds the carrier: peak bin over the bin floor ...
constexpr double kActiveRelDb = 10.0;   // ... and >= (median peak of the strongest half) - this
constexpr double kMinActive = 0.3;      // carrier in the captured band for >= 30 % of the period
constexpr double kMinSpanHz = 10e6;     // narrower periodic FM (video lines, FSK) is not a sweep
constexpr int kMaxFlyback = 1;          // frequency jumps per period (a sawtooth's flyback)
constexpr int kMaxExits = 2;            // band exits per period
constexpr double kMinDomDb = 4.0;       // folded peak over the swept span's median (75th pct of phases)
constexpr size_t kSmooth = 5;           // circular median over phases (none below 40 phases per period)
constexpr double kBridgeDb = 10.0;      // a jump is bridged when the bins between stay within this of the peaks
constexpr size_t kMaxOvershoot = 3;     // occupied band: at most this many bins beyond the trajectory's range
constexpr double kMinRateHz = 300.0;
constexpr size_t kTrajPoints = 32;      // trajectory points reported

// python occupancy(): Hann STFT power (hop 32) summed over frame pairs, (n/64 - 1) rows x 64 bins (FFT order)
fmat occupancy(const cfvec& x);

struct Sweep {
  bool found = false;
  double period = 0;                    // s
  size_t nph = 0;                       // fold phases
  dvec fk;                              // fractional peak bin per phase
  std::vector<uint8_t> active;          // phase holds the carrier
  dvec pk_db;                           // peak over the bin floor per phase
  double span = 0, dom75 = 0;
  int flyback = 0, exits = 0;
  std::vector<std::string> beyond;      // "above" / "below" the captured band, sorted
  dvec freqs;                           // absolute frequency per phase (Hz)
  std::pair<double, double> occupied_hz;
  double top_line_db = 0;
  Json tried = Json::array();           // per-candidate diagnostics (python measure()[1]["tried"])
};

Sweep measure(const cfvec& x, double fs, double fc);
std::vector<Emission> detect(const cfvec& x, double fs, double fc);

}  // namespace si::rx::jammer
