// Pulse deinterleaving (research, signal gaps S10; python: core/deinterleave.py): one signal's pulse list (TOA, PW)
// split into emitters by PRI (PW groups, CDIF, sequence search with a chance test, stagger, jitter mode).
#pragma once

#include <vector>

#include "core/types.hpp"

namespace si::deint {

// python deinterleave(toa_s, pw_s) -> {"emitters": [...], "unassigned": n}
Json deinterleave(const std::vector<double>& toa_s, const std::vector<double>& pw_s);

}  // namespace si::deint
