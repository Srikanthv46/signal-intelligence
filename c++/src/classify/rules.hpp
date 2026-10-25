// Rule-based classifier: India band prior + measured features -> label (python: classify/rules.py).
//
// Every rule returns (label, score, reasons) or nothing. Scores add up the independent
// evidence a rule found (band, bandwidth, numerology, timing); the best rule wins; below
// kMinConf the signal is "unknown" with the candidates listed. Inferred labels say "-like".
// A decoder result always overrides the rules.
#pragma once

#include "core/model.hpp"

namespace si::classify {

extern const double kMinConf;  // 0.5

void classify(Signal& s);

}  // namespace si::classify
