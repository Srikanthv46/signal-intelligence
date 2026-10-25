// Receiver artifacts that look like signals (python: detect/artifacts.py): LO leakage at
// DC, IQ-imbalance images, ADC clipping.
#pragma once

#include <optional>
#include <vector>

#include "core/model.hpp"

namespace si::detect {

extern const double kAdcFullScale;  // AD936x 12-bit, as written by signal-vault (raw counts)

Json clipping(const cfvec& x, double full_scale = kAdcFullScale);
void flag_lo_leakage(std::vector<Signal>& sigs, double fc, double df_fine);
// python flag_receiver_spurs(sigs, lines, fs, ppm): a narrow detection on a receiver line (core::spurs) is flagged
// receiver_spur; a wider one that contains a line only records it
void flag_receiver_spurs(std::vector<Signal>& sigs, const Json& lines, double fs, std::optional<double> ppm,
                         double max_bw = 25e3);
void flag_iq_images(std::vector<Signal>& sigs, double fc, double min_rej_db = 15.0);
std::vector<Signal> group_images(std::vector<Signal> sigs);

}  // namespace si::detect
