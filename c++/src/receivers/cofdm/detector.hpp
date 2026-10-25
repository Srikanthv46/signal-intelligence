// DVB-T-framed COFDM links (python receivers/cofdm/detector.py; research, synthetic validation only).
//   detect(x, fs, fc) -> Emission family "cofdm", label dvbt | dvbt_like
#pragma once

#include <vector>

#include "core/model.hpp"
#include "core/types.hpp"

namespace si::rx::cofdm {

std::vector<Emission> detect(const cfvec& x, double fs, double fc);
// python tps_decode on DBPSK bits b[l] (b[0] unused): the dict, or null
Json tps_decode(const std::vector<int>& bits);

}  // namespace si::rx::cofdm
