// Lookup tables for Mode S / ADS-B fields (python: receivers/adsb/tables.py): ICAO address blocks,
// operators, code meanings.
#pragma once

#include <optional>
#include <string>
#include <utility>

namespace si::rx::adsb {

// State of registration of a 24-bit ICAO address (Annex 10 Vol III Part I Ch. 9, a subset)
std::optional<std::string> country(const std::string& icao_hex);
// ICAO operator designator (Doc 8585) -> (operator, IATA code)
std::optional<std::pair<std::string, std::string>> operator_of(const std::string& icao3);
// emitter category of (typecode, ME CA field), DO-260B 2.2.3.2.5.2; nullopt -> "reserved"
std::optional<std::string> category(int tc, int ca);
std::optional<std::string> capability(int ca);        // CA field of DF11/DF17
const char* control_field(int cf);                    // DF18 CF, 0..7
const char* flight_status(int fs);                    // 0..7
const char* surveillance_status(int ss);              // 0..3
const char* emergency(int e);                         // 0..7
std::optional<std::string> acas_ri(int ri);           // DF0/DF16 reply information
const char* acas_sl(int sl);                          // 0..7
std::optional<int> nacp_epu_m(int nacp);

}  // namespace si::rx::adsb
