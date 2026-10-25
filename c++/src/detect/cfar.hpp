// Energy detection on the noise-normalised spectrogram -> time-frequency boxes
// (python: detect/cfar.py).
//
// continuous  per-bin 20th-percentile-over-time profile vs. the noise floor; each run of
//             bins >= 3 dB above it (grown to 1.5 dB) is one always-on emitter.
// bursts      hysteresis CFAR on P / max(floor, continuous level), two passes:
//             3x3 cells seed 6 dB / grow 3 dB, and 15x5 cells seed 3 dB / grow 1.5 dB;
//             8-connected components that hold a seed.
#pragma once

#include <array>
#include <optional>
#include <vector>

#include "core/spectrogram.hpp"

namespace si::detect {

struct BurstPass {
  size_t nt, nf;
  double seed_db, grow_db;
};
extern const std::vector<BurstPass> kBurstPasses;  // ((3,3,6,3), (15,5,3,1.5))
extern const double kContSeedDb, kContGrowDb;       // 3.0, 1.5

std::vector<Box> detect_continuous(const Spectrogram& S, const dvec& prof_db, const dvec& floor_db,
                                   const std::string& res, size_t min_bins = 2,
                                   std::optional<double> flank_db = std::nullopt);

std::vector<Box> detect_bursts(const Spectrogram& S, const dvec& floor_db, const dvec& baseline_db,
                               const std::string& res, const std::vector<BurstPass>& passes = kBurstPasses,
                               size_t min_cells = 6, double strong_db = 15.0,
                               std::optional<double> flank_db = std::nullopt);

// numpy mean of a 2-D float64 slice [r0, r1) x [c0, c1): rows reduced pairwise, rows added in order
double mean2d(const dmat& m, size_t r0, size_t r1, size_t c0, size_t c1);

}  // namespace si::detect
