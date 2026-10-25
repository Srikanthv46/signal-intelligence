// Per-aircraft summary from decoded Mode S frames (python: receivers/adsb/tracks.py): identity, CPR
// positions, kinematics, intent, integrity, capabilities, Comm-B air data, derived wind /
// temperature, the interrogation environment and plausibility flags.
#pragma once

#include <optional>
#include <vector>

#include "core/types.hpp"

namespace si::rx::adsb {

// python tracks(msgs, ref): msgs from fields() with "t" (seconds); ref = [lat, lon(, alt_m)] of the
// receiver (a Json array, numbers kept as given) or null. Sets frame["commb"]["register"].
// Throws std::runtime_error("KeyError: ...") where the Python raises (an ambiguous BDS50/60 reply
// that is50or60 cannot resolve), so the pipeline reports the same plugin error.
Json tracks(std::vector<Json>& msgs, const Json& ref);

}  // namespace si::rx::adsb
