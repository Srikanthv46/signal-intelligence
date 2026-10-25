// Receiver-generated lines and LO-harmonic images (python: core/spurs.py), shared by the pipeline and receivers.
//
// Receiver spurs: a line made by the receiver's own clocks is displayed at its exact nominal frequency (the LO
// and the sample clock share the reference): k x fs, or an integer MHz locked to the reference (that test needs
// the receiver ppm: an external line there would be displayed off by the receiver error).
// LO-harmonic images: a switching mixer also converts RF near n x LO (n odd); harmonic_image_sources() lists the
// strong-emitter bands an image would come from (a caution only).
#pragma once

#include <optional>
#include <utility>
#include <vector>

#include "core/bands_india.hpp"
#include "core/types.hpp"

namespace si::spurs {

constexpr double kLineTolHz = 150.0;
constexpr double kRefGridHz = 1e6;
constexpr double kLineMinDb = 10.0;

// python refine_line(X, fs, fc, f_abs, half_hz) -> (frequency Hz (numpy float64), peak over window median dB)
std::pair<double, double> refine_line(const cfvec& X, double fs, double fc, double f_abs, double half_hz);
std::optional<int64_t> sample_clock_harmonic(double f_abs, double fs, double tol_hz = kLineTolHz);
bool reference_locked(double f_abs, std::optional<double> ppm, double tol_hz = kLineTolHz);
// python receiver_lines(X, fs, fc, ppm, band): [{"freq_hz", "kind", "harmonic", "peak_db"}, ...]
Json receiver_lines(const cfvec& X, double fs, double fc, std::optional<double> ppm,
                    std::pair<double, double> band, double win_hz = 5e3);
// python harmonic_image_sources(f0, f1, lo, BANDS) / in_strong_band(f0, f1, BANDS)
Json harmonic_image_sources(double f0, double f1, double lo);
bool in_strong_band(double f0, double f1);

}  // namespace si::spurs
