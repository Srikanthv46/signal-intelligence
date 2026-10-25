// Physical-layer features of one detected signal, measured on its own IQ
// (python: features/measure.py). Each feature records how it was measured, so the
// classifier evidence and the report can cite it.
#pragma once

#include <string>
#include <utility>
#include <vector>

#include "core/dsp.hpp"
#include "core/model.hpp"

namespace si::features {

extern const size_t kMaxSamples;  // 2^19

Json timing(const Signal& sig, double capture_s);
Json pulse_train(const Signal& sig, const Channelizer& ch);
std::pair<Json, std::vector<std::string>> iq_features(const Signal& sig, const Channelizer& ch, size_t n);
std::pair<std::string, std::vector<std::string>> infer_modulation(const Signal& sig, const Json& f, double df_fine);
void extract(Signal& sig, const Channelizer& ch, size_t n, double capture_s, double df_fine);

}  // namespace si::features
