// Noise-floor estimation (python: detect/noise.py), per block.
//
// Per bin, the 20th percentile over time of the coarse spectrogram is the floor wherever a
// signal is off at least 80 % of the time; for exponential (complex-Gaussian) noise that
// percentile sits -ln(0.8) below the mean, which is corrected. A morphological opening
// (3 MHz) removes everything narrower than that from the profile. The Python library can
// also use a sweep template (neighbouring tunings); a production block has no sweep, so
// that path is not ported.
#pragma once

#include "core/spectrogram.hpp"

namespace si::detect {

extern const double kQ;          // percentile over time (20)
extern const double kQBiasDb;    // +6.52 dB: percentile -> mean for exponential noise

dvec profile_db(const fmat& P);                                   // per-bin noise-mean estimate, dB
dvec opening_floor(const dvec& prof_db, double df, double win_hz = 3e6);
std::pair<dvec, Json> estimate_floor(const dvec& prof_db, double fc, double df);
// per-frame noise-level multiplier g(t) >= 1 (float32)
fvec time_gain(const fmat& P, const dvec& floor_db, const dvec& prof_db, double q = 30.0, size_t smooth = 5,
               double max_excess_db = 1.0, size_t min_bins = 32);

}  // namespace si::detect
