// python: receivers/adsb/tables.py (same entries, same order).
#include "receivers/adsb/tables.hpp"

#include <cstdint>
#include <cstdlib>
#include <map>

namespace si::rx::adsb {

namespace {

struct Block {
  uint32_t lo, hi;
  const char* name;
};

const Block kBlocks[] = {
    {0x004000, 0x0043FF, "Zimbabwe"}, {0x008000, 0x00FFFF, "South Africa"}, {0x010000, 0x017FFF, "Egypt"},
    {0x018000, 0x01FFFF, "Libya"}, {0x020000, 0x027FFF, "Morocco"}, {0x040000, 0x040FFF, "Ethiopia"},
    {0x04C000, 0x04CFFF, "Kenya"}, {0x05A000, 0x05A3FF, "Maldives"}, {0x064000, 0x064FFF, "Nigeria"},
    {0x06A000, 0x06A3FF, "Qatar"}, {0x0A0000, 0x0A7FFF, "Algeria"}, {0x0D0000, 0x0D7FFF, "Mexico"},
    {0x100000, 0x1FFFFF, "Russia"}, {0x300000, 0x33FFFF, "Italy"}, {0x340000, 0x37FFFF, "Spain"},
    {0x380000, 0x3BFFFF, "France"}, {0x3C0000, 0x3FFFFF, "Germany"}, {0x400000, 0x43FFFF, "United Kingdom"},
    {0x440000, 0x447FFF, "Austria"}, {0x448000, 0x44FFFF, "Belgium"}, {0x450000, 0x457FFF, "Bulgaria"},
    {0x458000, 0x45FFFF, "Denmark"}, {0x460000, 0x467FFF, "Finland"}, {0x468000, 0x46FFFF, "Greece"},
    {0x470000, 0x477FFF, "Hungary"}, {0x478000, 0x47FFFF, "Norway"}, {0x480000, 0x487FFF, "Netherlands"},
    {0x488000, 0x48FFFF, "Poland"}, {0x490000, 0x497FFF, "Portugal"}, {0x498000, 0x49FFFF, "Czechia"},
    {0x4A0000, 0x4A7FFF, "Romania"}, {0x4A8000, 0x4AFFFF, "Sweden"}, {0x4B0000, 0x4B7FFF, "Switzerland"},
    {0x4B8000, 0x4BFFFF, "Turkey"}, {0x4CA000, 0x4CAFFF, "Ireland"}, {0x508000, 0x50FFFF, "Ukraine"},
    {0x680000, 0x6803FF, "Bhutan"}, {0x700000, 0x700FFF, "Afghanistan"}, {0x702000, 0x702FFF, "Bangladesh"},
    {0x704000, 0x704FFF, "Myanmar"}, {0x706000, 0x706FFF, "Kuwait"}, {0x708000, 0x708FFF, "Laos"},
    {0x70A000, 0x70AFFF, "Nepal"}, {0x70C000, 0x70C3FF, "Oman"}, {0x70E000, 0x70EFFF, "Cambodia"},
    {0x710000, 0x717FFF, "Saudi Arabia"}, {0x718000, 0x71FFFF, "South Korea"}, {0x720000, 0x727FFF, "North Korea"},
    {0x728000, 0x72FFFF, "Iraq"}, {0x730000, 0x737FFF, "Iran"}, {0x738000, 0x73FFFF, "Israel"},
    {0x740000, 0x747FFF, "Jordan"}, {0x748000, 0x74FFFF, "Lebanon"}, {0x750000, 0x757FFF, "Malaysia"},
    {0x758000, 0x75FFFF, "Philippines"}, {0x760000, 0x767FFF, "Pakistan"}, {0x768000, 0x76FFFF, "Singapore"},
    {0x770000, 0x777FFF, "Sri Lanka"}, {0x778000, 0x77FFFF, "Syria"}, {0x780000, 0x7BFFFF, "China"},
    {0x7C0000, 0x7FFFFF, "Australia"}, {0x800000, 0x83FFFF, "India"}, {0x840000, 0x87FFFF, "Japan"},
    {0x880000, 0x887FFF, "Thailand"}, {0x888000, 0x88FFFF, "Viet Nam"}, {0x890000, 0x890FFF, "Yemen"},
    {0x894000, 0x894FFF, "Bahrain"}, {0x895000, 0x8953FF, "Brunei"}, {0x896000, 0x896FFF, "United Arab Emirates"},
    {0x899000, 0x8997FF, "Taiwan"}, {0x8A0000, 0x8A7FFF, "Indonesia"}, {0xA00000, 0xAFFFFF, "United States"},
    {0xC00000, 0xC3FFFF, "Canada"}, {0xC80000, 0xC87FFF, "New Zealand"}, {0xE00000, 0xE3FFFF, "Argentina"},
    {0xE40000, 0xE7FFFF, "Brazil"}, {0xE80000, 0xE80FFF, "Chile"},
};

const std::map<std::string, std::pair<std::string, std::string>>& operators() {
  static const std::map<std::string, std::pair<std::string, std::string>> m = {
      {"IGO", {"IndiGo", "6E"}}, {"AIC", {"Air India", "AI"}}, {"AXB", {"Air India Express", "IX"}},
      {"SEJ", {"SpiceJet", "SG"}}, {"AKJ", {"Akasa Air", "QP"}}, {"VTI", {"Vistara", "UK"}},
      {"IAD", {"AirAsia India", "I5"}}, {"LLR", {"Alliance Air", "9I"}}, {"FLG", {"Fly91", "IC"}},
      {"SDG", {"Star Air", "S5"}}, {"BDA", {"Blue Dart Aviation", "BZ"}}, {"QTR", {"Qatar Airways", "QR"}},
      {"UAE", {"Emirates", "EK"}}, {"ETD", {"Etihad", "EY"}}, {"FDB", {"flydubai", "FZ"}}, {"ABY", {"Air Arabia", "G9"}},
      {"OMA", {"Oman Air", "WY"}}, {"SVA", {"Saudia", "SV"}}, {"GFA", {"Gulf Air", "GF"}}, {"KAC", {"Kuwait Airways", "KU"}},
      {"SIA", {"Singapore Airlines", "SQ"}}, {"MAS", {"Malaysia Airlines", "MH"}}, {"THA", {"Thai Airways", "TG"}},
      {"ALK", {"SriLankan Airlines", "UL"}}, {"BBC", {"Biman Bangladesh", "BG"}}, {"RNA", {"Nepal Airlines", "RA"}},
      {"CPA", {"Cathay Pacific", "CX"}}, {"BAW", {"British Airways", "BA"}}, {"DLH", {"Lufthansa", "LH"}},
      {"AFR", {"Air France", "AF"}}, {"KLM", {"KLM", "KL"}}, {"THY", {"Turkish Airlines", "TK"}},
      {"UAL", {"United Airlines", "UA"}}, {"AAL", {"American Airlines", "AA"}}, {"FDX", {"FedEx", "FX"}},
      {"UPS", {"UPS", "5X"}}, {"SWR", {"Swiss", "LX"}}, {"VIR", {"Virgin Atlantic", "VS"}},
      {"AXM", {"AirAsia", "AK"}}, {"TGW", {"Scoot", "TR"}}, {"JAI", {"Jet Airways", "9W"}},
  };
  return m;
}

}  // namespace

std::optional<std::string> country(const std::string& icao_hex) {
  const unsigned long a = std::strtoul(icao_hex.c_str(), nullptr, 16);
  for (const auto& b : kBlocks)
    if (b.lo <= a && a <= b.hi) return std::string(b.name);
  return std::nullopt;
}

std::optional<std::pair<std::string, std::string>> operator_of(const std::string& icao3) {
  auto it = operators().find(icao3);
  if (it == operators().end()) return std::nullopt;
  return it->second;
}

std::optional<std::string> category(int tc, int ca) {
  static const char* const k4[] = {"no category info", "light (< 15,500 lb)", "small (15,500-75,000 lb)",
                                   "large (75,000-300,000 lb)", "high-vortex large (B757)", "heavy (> 300,000 lb)",
                                   "high performance (> 5 g, > 400 kt)", "rotorcraft"};
  static const char* const k3[] = {"no category info", "glider / sailplane", "lighter-than-air",
                                   "parachutist / skydiver", "ultralight / hang-glider", "reserved",
                                   "unmanned aerial vehicle", "space / trans-atmospheric vehicle"};
  static const char* const k2[] = {"no category info", "surface emergency vehicle", "surface service vehicle",
                                   "point obstacle", "cluster obstacle", "line obstacle", "reserved", "reserved"};
  if (ca < 0 || ca > 7) return std::nullopt;
  if (tc == 4) return std::string(k4[ca]);
  if (tc == 3) return std::string(k3[ca]);
  if (tc == 2) return std::string(k2[ca]);
  return std::nullopt;
}

std::optional<std::string> capability(int ca) {
  switch (ca) {
    case 0: return std::string("level 1 transponder");
    case 4: return std::string("level 2+ transponder, on ground");
    case 5: return std::string("level 2+ transponder, airborne");
    case 6: return std::string("level 2+ transponder, ground or airborne");
    case 7: return std::string("DR != 0 or FS = 2-5 (alert/SPI)");
    default: return std::nullopt;
  }
}

const char* control_field(int cf) {
  static const char* const k[] = {"ADS-B ES/NT, ICAO address", "ADS-B ES/NT, non-ICAO address", "TIS-B fine format",
                                  "TIS-B coarse format", "TIS-B / ADS-R management", "TIS-B non-ICAO address",
                                  "ADS-R (rebroadcast of UAT)", "reserved"};
  return k[cf & 7];
}

const char* flight_status(int fs) {
  static const char* const k[] = {"airborne", "on ground", "alert, airborne", "alert, on ground",
                                  "alert + SPI", "SPI (ident)", "reserved", "not assigned"};
  return k[fs & 7];
}

const char* surveillance_status(int ss) {
  static const char* const k[] = {"no condition", "permanent alert (emergency squawk)",
                                  "temporary alert (squawk changed)", "SPI (ident)"};
  return k[ss & 3];
}

const char* emergency(int e) {
  static const char* const k[] = {"no emergency", "general emergency", "lifeguard / medical", "minimum fuel",
                                  "no communications", "unlawful interference", "downed aircraft", "reserved"};
  return k[e & 7];
}

std::optional<std::string> acas_ri(int ri) {
  switch (ri) {
    case 0: return std::string("no operating ACAS");
    case 2: return std::string("ACAS, resolution capability inhibited");
    case 3: return std::string("ACAS, vertical-only RAs");
    case 4: return std::string("ACAS, vertical and horizontal RAs");
    case 8: return std::string("no max airspeed data");
    case 9: return std::string("max airspeed <= 75 kt");
    case 10: return std::string("max airspeed 75-150 kt");
    case 11: return std::string("max airspeed 150-300 kt");
    case 12: return std::string("max airspeed 300-600 kt");
    case 13: return std::string("max airspeed 600-1200 kt");
    case 14: return std::string("max airspeed > 1200 kt");
    default: return std::nullopt;
  }
}

const char* acas_sl(int sl) {
  static const char* const k[] = {"ACAS inoperative", "ACAS standby", "TA only", "TA/RA",
                                  "TA/RA", "TA/RA", "TA/RA", "TA/RA"};
  return k[sl & 7];
}

std::optional<int> nacp_epu_m(int nacp) {
  static const int k[] = {0, 18520, 7408, 3704, 1852, 926, 556, 185, 93, 30, 10, 3};
  if (nacp >= 1 && nacp <= 11) return k[nacp];
  return std::nullopt;
}

}  // namespace si::rx::adsb
