// Standalone burst source for emitter association (python: receivers/emitters/bursts.py
// extract_bursts): 15 kHz-SCS OFDM bursts (DJI OcuSync-class RC uplink hops, DroneID, video)
// from a capture via core/ofdm15k, as Emission records. The synthetic-scene helpers of the
// Python module (synth_scene, inject_copy) are test tools and are not ported.
#pragma once

#include <vector>

#include "core/model.hpp"
#include "core/types.hpp"

namespace si::rx::emitters {

// python extract_bursts(x, fs, fc, usable_half=None): sorted by start
std::vector<Emission> extract_bursts(const cfvec& x, double fs, double fc, double usable_half = 0);

}  // namespace si::rx::emitters
