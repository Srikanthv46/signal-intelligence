// Detection stage: IQ -> energy signals with provenance (python: detect/pipeline.py).
#pragma once

#include "core/spectrogram.hpp"
#include "io/sigmf.hpp"

namespace si::detect {

extern const size_t kCoarseNfft;  // 1024: 60 kHz x 16.7 us at 61.44 MSps
extern const size_t kFineNfft;    // 16384: 3.75 kHz x 267 us
extern const size_t kFineHop;     // 16384: no overlap (CFAR statistics assume independent cells)

DetectResult detect(const Capture& cap, double usable_frac = 0.95);

}  // namespace si::detect
