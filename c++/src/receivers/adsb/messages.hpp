// Per-frame Mode S / ADS-B field decoding (python: receivers/adsb/messages.py).
#pragma once

#include <vector>

#include "core/types.hpp"

namespace si::rx::adsb {

// python fields(msgs): adds the decoded fields of every DF / typecode / Comm-B register to each
// frame dict from demod() (in Python's key order). Multi-frame fields are left to tracks().
void fields(std::vector<Json>& msgs);

}  // namespace si::rx::adsb
