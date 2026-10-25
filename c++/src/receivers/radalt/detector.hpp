// FMCW radio altimeters, 4200-4400 MHz (research, signal gaps S8; python: receivers/radalt/detector.py).
#pragma once

#include <vector>

#include "core/model.hpp"
#include "core/types.hpp"

namespace si::rx::radalt {

// python detect(x, fs, fc): at most one emission, family "radalt", label "radio_altimeter"
std::vector<Emission> detect(const cfvec& x, double fs, double fc);

}  // namespace si::rx::radalt
