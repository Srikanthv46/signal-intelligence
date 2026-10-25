// Which transmitter sent each 802.11 PPDU (python receivers/wifi/links.py), for the OFDM
// Wi-Fi plugins: MAC header, control responses after SIFS, sequential scrambler-seed
// counters, then a weak received-power hint.
#pragma once

#include <vector>

#include "core/types.hpp"

namespace si::rx::wifi::links {

// recs: the PPDU records (Json objects, any order) with t0/t1 (s), fmt, mac (dissected MPDU
// when the FCS passed), scrambler_seed / seed_ok, power_db, cfo_khz. Like Python, every
// record is updated in place (ta, ta_evidence, ta_likely removed, then set where attributed;
// set keys are appended at the end of the record). Returns one row per transmitter.
Json attribute(std::vector<Json>& recs, double sifs_us = 16.0);
// same, on records that live elsewhere (e.g. inside a larger Json); none may be null
Json attribute(const std::vector<Json*>& recs, double sifs_us = 16.0);

}  // namespace si::rx::wifi::links
