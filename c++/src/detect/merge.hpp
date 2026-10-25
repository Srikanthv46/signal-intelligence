// Box clean-up (python: detect/merge.py): fuse the two resolutions, merge fragments,
// split adjacent carriers at spectral valleys, group bursts into per-emitter signals.
#pragma once

#include <utility>
#include <vector>

#include "core/spectrogram.hpp"

namespace si::detect {

extern const double kNarrowHz;  // 300 kHz: below it the fine (3.75 kHz) spectrogram owns the box

std::vector<Box> fuse_resolutions(const std::vector<Box>& coarse, const std::vector<Box>& fine);
std::vector<Box> merge_fragments(std::vector<Box> boxes, double f_gap = 130e3, double t_gap = 40e-6);
std::vector<Box> split_valleys(const Box& b, const Spectrogram& S, const dvec& floor_db, double depth_db = 6.0,
                               double min_part_hz = 600e3);
std::vector<Signal> group_signals(std::vector<Box> boxes, double iou_min = 0.6);
std::vector<Signal> absorb_into_continuous(std::vector<Signal> sigs, double frac = 0.8);
// python merge_pulse_trains(): pieces of one pulse train (same PRI and pulse times, from features "pri_s" /
// "_pulse_t_s") and the burst-only slivers on its pulses -> one signal. Returns the signals and the indices (into
// them) of the merged ones, whose features must be measured again.
std::pair<std::vector<Signal>, std::vector<size_t>> merge_pulse_trains(std::vector<Signal> sigs);

}  // namespace si::detect
