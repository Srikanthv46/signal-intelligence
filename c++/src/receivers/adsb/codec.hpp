// Mode S / ADS-B field encodings, implemented from the public standards (no pyModeS code):
//   ICAO Annex 10 Vol IV (3.1.2.6.5.4 altitude code AC, 3.1.2.6.7.1 identity code ID, CRC-24),
//   DO-260B / ICAO Doc 9871 (ADS-B ME formats, CPR A.1.7, Comm-B MB registers BDS 1,0-6,0),
//   ICAO Doc 7488 (ISA atmosphere) for the airspeed conversions of the Comm-B plausibility checks.
// The Python receiver uses pyModeS for these encodings; the presentation conventions its output
// carries (int vs float values, '_' for a space in callsigns, register names "BDS50", the
// plausibility limits of the register inference) are reproduced and verified black-box against it
// (tests/vectors/adsb, test_rx_adsb.cpp).
#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <utility>

#include "core/types.hpp"

namespace si::rx::adsb {

// "8D4840D6..." -> "10001101..." (python bin(int(h, 16))[2:].zfill(4 * len(h)))
std::string hex_bits(const std::string& hex);
// int(b[i:j], 2)
int64_t ubits(const std::string& b, size_t i, size_t j);

// 13-bit altitude code (AC, Annex 10 Vol IV 3.1.2.6.5.4): ft, or nullopt when not available
std::optional<int64_t> altitude_ac13(const std::string& b13);
// 13-bit identity code (ID: C1 A1 C2 A2 C4 A4 X B1 D1 B2 D2 B4 D4) -> "ABCD" octal digits
std::string squawk_id13(const std::string& b13);
// 6-bit character set (IA-5 subset): A-Z, '_' for space, 0-9, '#' for unassigned codes
char callsign_char(int c);

// ---- DF17/18 extended squitter (h = 28 hex digits)
int es_typecode(const std::string& bits);
std::string es_callsign(const std::string& bits);                 // TC1-4, unassigned codes dropped
Json es_altitude(const std::string& bits, int tc);                  // TC9-18 baro (AC12), TC20-22 GNSS
std::pair<Json, Json> es_surface_movement(const std::string& bits); // TC5-8: (speed kt, track deg)
struct Velocity {                                                   // TC19 (DO-260B 2.2.3.2.6)
  Json spd, ang, vr;
  std::string kind, ang_src, vr_src;
};
Velocity es_velocity(const std::string& bits);
Json es_gnss_baro_diff(const std::string& bits);                    // TC19 ME 49-56
std::string es_emergency_squawk(const std::string& bits);           // TC28/1 Mode A code

// ---- CPR (DO-260B A.1.7)
int cpr_nl(double lat);
std::optional<std::pair<double, double>> cpr_airborne_global(const std::string& even_hex, const std::string& odd_hex,
                                                             double t_even, double t_odd);
std::pair<double, double> cpr_local(const std::string& hex, double lat_ref, double lon_ref, bool surface);

// ---- Comm-B (Doc 9871), mb = the 56 MB bits of a DF20/21 reply
// register inference over BDS 1,0 1,7 2,0 3,0 4,0 4,4 4,5 5,0 6,0 ("BDS40,BDS50"), "EMPTY" or null
Json commb_infer(const std::string& bits);
// the field values of one register: name of the field function -> value (null when not available)
Json commb_field(const std::string& fn, const std::string& bits);
Json commb_cap17(const std::string& bits);
// BDS50 or BDS60 for an MB that is both, from the aircraft's ground speed, track and altitude
Json commb_is50or60(const std::string& bits, double spd_ref, double trk_ref, double alt_ref);

}  // namespace si::rx::adsb
