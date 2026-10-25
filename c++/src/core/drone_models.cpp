// python core/drone_models.py (see there for the sources; facts re-typed, no code copied). The tables below are
// generated from the python module (DJI_PRODUCT_TYPES, MFR_CODES, SERIAL_MODELS) so both stay identical.
#include "core/drone_models.hpp"

#include <map>
#include <set>

namespace si::drone_models {
namespace {
const std::map<int64_t, const char*> kDji = {
    {1, "Inspire 1"},
    {2, "Phantom 3 series"},
    {3, "Phantom 3 series"},
    {4, "Phantom 3 Standard"},
    {5, "Matrice 100"},
    {6, "ACEONE"},
    {7, "WooKong-M"},
    {8, "NAZA"},
    {9, "A2"},
    {10, "A3"},
    {11, "Phantom 4"},
    {12, "MG-1"},
    {14, "Matrice 600"},
    {15, "Phantom 3 4K"},
    {16, "Mavic Pro"},
    {17, "Inspire 2"},
    {18, "Phantom 4 Pro"},
    {20, "N2"},
    {21, "Spark"},
    {23, "Matrice 600 Pro"},
    {24, "Mavic Air"},
    {25, "Matrice 200"},
    {26, "Phantom 4 series"},
    {27, "Phantom 4 Advanced"},
    {28, "Matrice 210"},
    {30, "Matrice 210 RTK"},
    {31, "A3 AG"},
    {32, "MG-2"},
    {34, "MG-1A"},
    {35, "Phantom 4 RTK"},
    {36, "Phantom 4 Pro V2.0"},
    {38, "MG-1P"},
    {40, "MG-1P RTK"},
    {41, "Mavic 2"},
    {44, "Matrice 200 V2 series"},
    {51, "Mavic 2 Enterprise"},
    {53, "Mavic Mini"},
    {58, "Mavic Air 2"},
    {59, "Phantom 4 Multispectral"},
    {60, "Matrice 300 RTK"},
    {61, "DJI FPV"},
    {63, "Mini 2"},
    {64, "Agras T10"},
    {65, "Agras T30"},
    {66, "Air 2S"},
    {67, "Matrice 30"},
    {68, "Mavic 3"},
    {69, "Mavic 2 Enterprise Advanced"},
    {70, "Mini SE"},
    {72, "Agras T40"},
    {73, "Mini 3 Pro"},
    {75, "Avata"},
    {76, "Inspire 3"},
    {77, "Mavic 3 Enterprise"},
    {78, "FlyCart 30"},
    {82, "Agras T25"},
    {83, "Agras T50"},
    {84, "Mavic 3 Pro"},
    {86, "Mavic 3 Classic"},
    {87, "Mini 3"},
    {88, "Mini 2 SE"},
    {89, "Matrice 350 RTK"},
    {90, "Air 3"},
    {91, "Matrice 3D / 3TD"},
    {93, "Mini 4 Pro"},
    {94, "Avata 2"},
    {95, "Agras T60"},
    {96, "Agras T25P"},
    {99, "Matrice 4T / 4E"},
    {100, "Matrice 4D / 4TD"},
    {102, "Flip"},
    {103, "Matrice 400"},
    {104, "Neo"},
    {105, "Air 3S"},
    {108, "Mavic 4 Pro"},
    {109, "Agras T100"},
    {112, "Mini 5 Pro"},
    {124, "Neo 2"},
    {127, "Avata 360"},
    {128, "Agras T100S / T70S"},
    {130, "Matrice 5TD"},
    {132, "Agras T55"},
    {133, "Lito X1"},
    {140, "Lito 1"},
    {166, "Agras T200"},
    {167, "FlyCart 200"},
    {173, "Mini 4K"}};
constexpr int64_t kDjiConfirmed = 70;
struct Mfr { const char* maker; const char* kind; };
const std::map<std::string, Mfr> kMfr = {
    {"1581", {"DJI", "aircraft"}},
    {"1748", {"Autel Robotics", "aircraft"}},
    {"1588", {"Parrot", "aircraft"}},
    {"1668", {"Skydio", "aircraft"}},
    {"1596", {"Dronetag", "module"}},
    {"1863", {"XAG", "aircraft"}},
    {"2104", {"Talos", "aircraft"}},
    {"1865", {"Holy Stone", "aircraft"}},
    {"1910", {"Potensic", "aircraft"}},
    {"1587", {"AgEagle / senseFly", "aircraft"}},
    {"1792", {"uAvionix", "module"}},
    {"1817", {"Freefly", "aircraft"}},
    {"1849", {"Spektrum", "module"}},
    {"1930", {"SwellPro", "aircraft"}},
    {"1972", {"Lumenier", "aircraft"}},
    {"1973", {"FrSky", "module"}},
    {"2071", {"NewBeeDrone", "aircraft"}},
    {"1892", {"SJRC", "aircraft"}},
    {"1875", {"Bwine", "aircraft"}},
    {"1710", {"Atom UAV", "aircraft"}},
    {"1902", {"Elsight", "module"}},
    {"1595", {"Sunflower Labs", "aircraft"}},
    {"1605", {"Pyka", "aircraft"}},
    {"1705", {"Pyka", "aircraft"}},
    {"1607", {"Wingcopter", "aircraft"}},
    {"1620", {"Aero Systems West", "aircraft"}},
    {"1633", {"Wing Aviation", "aircraft"}},
    {"1663", {"Wing Aviation", "aircraft"}},
    {"1634", {"ACSL", "aircraft"}},
    {"1641", {"Applied Aeronautics", "aircraft"}},
    {"1644", {"Microdrones", "aircraft"}},
    {"1656", {"Quantum-Systems", "aircraft"}},
    {"1698", {"RigiTech", "aircraft"}},
    {"1701", {"Percepto", "aircraft"}},
    {"1718", {"Censys", "aircraft"}},
    {"1720", {"Leica Geosystems", "aircraft"}},
    {"1728", {"Sony", "aircraft"}},
    {"1744", {"Wingtra", "aircraft"}},
    {"1759", {"ABZ Innovation", "aircraft"}},
    {"1767", {"Teledyne FLIR", "aircraft"}},
    {"1774", {"unknown (declared make 'UAS')", "module"}},
    {"1785", {"Watts Innovations", "aircraft"}},
    {"1786", {"Amazon Prime Air", "aircraft"}},
    {"1787", {"BlueMark", "module"}},
    {"1788", {"Vision Aerial", "aircraft"}},
    {"1790", {"Hylio", "aircraft"}},
    {"1791", {"Inspired Flight", "aircraft"}},
    {"1795", {"Zephyr Systems", "aircraft"}},
    {"1796", {"Drone Defence", "module"}},
    {"1798", {"Sentera", "module"}},
    {"1801", {"Drone Volt", "aircraft"}},
    {"1807", {"Event38", "aircraft"}},
    {"1808", {"Draganfly", "aircraft"}},
    {"1809", {"Aerobits (idME)", "module"}},
    {"1810", {"Matternet", "aircraft"}},
    {"1812", {"Ascent AeroSystems", "aircraft"}},
    {"1813", {"Harris Aerial", "aircraft"}},
    {"1816", {"Field / Explorer", "aircraft"}},
    {"1836", {"Coretronic", "aircraft"}},
    {"1839", {"Teal Drones", "aircraft"}},
    {"1840", {"Voliro", "aircraft"}},
    {"1843", {"blueflite", "aircraft"}},
    {"1860", {"XiL", "aircraft"}},
    {"1862", {"DeltaQuad", "aircraft"}},
    {"1866", {"SmartDrone", "aircraft"}},
    {"1867", {"Hoverfly", "aircraft"}},
    {"1869", {"Ruko", "aircraft"}},
    {"1871", {"Camflite", "aircraft"}},
    {"1872", {"Veeniix", "module"}},
    {"1873", {"Zipline", "aircraft"}},
    {"1877", {"Lucid Bots", "aircraft"}},
    {"1879", {"Flightwave", "aircraft"}},
    {"1884", {"Workhorse Aero", "aircraft"}},
    {"1886", {"Flytrex", "aircraft"}},
    {"1889", {"WISPR", "aircraft"}},
    {"1891", {"Zing Drone Solutions", "module"}},
    {"1893", {"Pierce Aerospace", "module"}},
    {"1894", {"Drone-Clone Xperts", "aircraft"}},
    {"1895", {"EXO", "aircraft"}},
    {"1898", {"Fotokite", "aircraft"}},
    {"1900", {"Tucok", "aircraft"}},
    {"1901", {"Meteodrone", "aircraft"}},
    {"1904", {"ShenZhen C-Fly", "aircraft"}},
    {"1906", {"Aurelia", "module"}},
    {"1907", {"Anzu Robotics", "aircraft"}},
    {"1908", {"MicroFlight", "module"}},
    {"1911", {"Cider Copter", "aircraft"}},
    {"1914", {"BRINC", "aircraft"}},
    {"1922", {"TRIPLEFINE", "aircraft"}},
    {"1924", {"Arcsky", "aircraft"}},
    {"1925", {"DIY-D", "module"}},
    {"1926", {"Skyfront", "aircraft"}},
    {"1928", {"XM2 Labs", "aircraft"}},
    {"1933", {"Triad Drones", "aircraft"}},
    {"1937", {"PrecisionVision", "aircraft"}},
    {"1939", {"SPECTA", "aircraft"}},
    {"1943", {"LISHI", "aircraft"}},
    {"1947", {"Swoop Aero", "aircraft"}},
    {"1949", {"SKYFISH", "aircraft"}},
    {"1961", {"Phoenix UAS", "module"}},
    {"1968", {"Light Quadcopter UAV", "aircraft"}},
    {"1969", {"SkySpecs", "aircraft"}},
    {"1971", {"FIMI", "aircraft"}},
    {"1984", {"Nokia", "aircraft"}},
    {"1990", {"FACEGLE", "aircraft"}},
    {"1995", {"Blackbird", "aircraft"}},
    {"2003", {"Holyton", "module"}},
    {"2009", {"Cinnamon Drones", "module"}},
    {"2021", {"SHANGHENG", "aircraft"}},
    {"2030", {"Telegrid", "aircraft"}},
    {"2046", {"Aurora Flight Sciences", "aircraft"}},
    {"2047", {"Airobotics", "aircraft"}},
    {"2049", {"Aeroo", "aircraft"}},
    {"2054", {"Global Drone", "aircraft"}},
    {"2055", {"Gleesfun", "module"}},
    {"2061", {"Aerostar Dynamics", "aircraft"}},
    {"2068", {"Black Swift", "aircraft"}},
    {"2078", {"CHCNAV", "aircraft"}},
    {"2095", {"Zero Zero Robotics", "aircraft"}},
    {"2098", {"XPECE", "aircraft"}},
    {"2105", {"Flyby Robotics", "aircraft"}},
    {"2108", {"Trace / Vesper UAS", "aircraft"}},
    {"2110", {"DMR Technologies", "aircraft"}},
    {"2111", {"TopXGun", "aircraft"}},
    {"2123", {"Dunlevy Consulting", "aircraft"}},
    {"2124", {"Avol", "aircraft"}},
    {"2125", {"Exedy Drones", "aircraft"}},
    {"2126", {"Caelone", "aircraft"}},
    {"2129", {"Heilongjiang Huida", "aircraft"}},
    {"2132", {"Angel Aerial Systems", "aircraft"}},
    {"2137", {"Antigravity", "aircraft"}},
    {"2146", {"Contixo", "aircraft"}},
    {"2148", {"Sustainable Skylines", "aircraft"}},
    {"2157", {"Ceres Air", "aircraft"}},
    {"2160", {"Rainmaker", "aircraft"}},
    {"2161", {"Tofect", "aircraft"}},
    {"2164", {"Terra Drone", "aircraft"}},
    {"2165", {"A2Z", "aircraft"}},
    {"2183", {"SubUAS", "aircraft"}},
    {"2187", {"Amor Fati", "aircraft"}},
    {"2199", {"Voltair", "aircraft"}},
    {"2218", {"Manna", "aircraft"}}};
struct SerialRow { const char* prefix; std::vector<const char*> models; bool verified; const char* doc; };
const std::vector<SerialRow> kSerial = {
    {"1581F0UY", {"Phantom 4 RTK"}, true, "RID000000752"},
    {"1581F11V", {"Phantom 4 Pro V2.0"}, true, "RID000000759"},
    {"1581F0M6", {"Mavic 2 Zoom"}, true, "RID000000758"},
    {"1581F163", {"Mavic 2 Pro"}, true, "RID000000757"},
    {"1581F4GC", {"Mavic 2 Enterprise Advanced"}, true, "RID000000753"},
    {"1581F1WN", {"Mavic Air 2"}, true, "RID000000754"},
    {"1581F3N3", {"Mavic Air 2"}, true, "RID000000755"},
    {"1581F37Q", {"DJI FPV"}, true, "RID000000756"},
    {"1581F3YT", {"Air 2S"}, true, "RID000000007"},
    {"1581F4XF", {"Mini 3 Pro"}, true, "RID000000001"},
    {"1581F6L5", {"Mini 3 Pro"}, true, "RID000000001"},
    {"1581F5YH", {"Mini 3"}, true, "RID000000062"},
    {"1581F6Z9", {"Mini 4 Pro"}, true, "RID000000713"},
    {"1581F9DE", {"Mini 5 Pro"}, true, "RID000002581"},
    {"1581F6N8", {"Air 3"}, true, "RID000000714"},
    {"1581F895", {"Air 3S"}, true, "RID000002046"},
    {"1581F4QW", {"Avata"}, true, "RID000000008"},
    {"1581F6W8", {"Avata 2"}, true, "RID000001376"},
    {"1581FA8J", {"Avata 360"}, true, "RID000002651"},
    {"1581FBLK", {"Avata 360"}, true, "RID000002650"},
    {"1581F45T", {"Mavic 3"}, true, "RID000000009"},
    {"1581F4QZ", {"Mavic 3 Cine"}, true, "RID000000010"},
    {"1581F67P", {"Mavic 3 Classic"}, true, "RID000000041"},
    {"1581F67Q", {"Mavic 3 Pro"}, true, "RID000000204"},
    {"1581F6MK", {"Mavic 3 Pro Cine"}, true, "RID000000205"},
    {"1581F8LQ", {"Mavic 4 Pro"}, true, "RID000002375"},
    {"1581F986", {"Mavic 4 Pro"}, true, "RID000002376"},
    {"1581F4Z4", {"Inspire 3"}, true, "RID000000192"},
    {"1581F5FH", {"Mavic 3 Enterprise M3E"}, true, "RID000000019"},
    {"1581F5FJ", {"Mavic 3 Enterprise M3T"}, true, "RID000000020"},
    {"1581F5FK", {"Mavic 3 Enterprise M3M"}, true, "RID000000057"},
    {"1581FBEX", {"M3TA"}, true, "RID000002584"},
    {"1581F6QA", {"Matrice 3D"}, true, "RID000001166"},
    {"1581F6Q8", {"Matrice 3TD"}, true, "RID000001167"},
    {"1581F4BN", {"Matrice 30"}, true, "RID000000011"},
    {"1581F5BK", {"Matrice 30T"}, true, "RID000000012"},
    {"1581F5BL", {"Matrice 30 Dock"}, true, "RID000000273"},
    {"1581F5BM", {"Matrice 30T Dock"}, true, "RID000000274"},
    {"1581F62H", {"Matrice 30T (EU / KR)"}, true, "RID000001136"},
    {"1581F1ZN", {"Matrice 300 RTK"}, true, "RID000000021"},
    {"1581F6GK", {"Matrice 350 RTK"}, true, "RID000000256"},
    {"1581F7FV", {"Matrice 4E"}, true, "RID000002198"},
    {"1581F7K3", {"Matrice 4T"}, true, "RID000002199"},
    {"1581F8HH", {"Matrice 4D"}, true, "RID000002246"},
    {"1581F8HG", {"Matrice 4TD"}, true, "RID000002247"},
    {"1581F8DB", {"Matrice 400"}, true, "RID000002486"},
    {"1581F7C6", {"FlyCart 30"}, true, "RID000000760"},
    {"1581FAN4", {"FlyCart 100"}, true, "RID000002620"},
    {"1581F446", {"Agras T10", "Agras T30"}, true, "RID000001163 / 038"},
    {"1581F574", {"Agras T20P", "Agras T40"}, true, "RID000001181 / 039"},
    {"1581F6BU", {"Agras T50 / T25"}, true, "RID000001594"},
    {"1581F836", {"Agras T25P"}, true, "RID000002587"},
    {"1581F8ZX", {"Agras T70P"}, true, "RID000002579"},
    {"1581F8ZL", {"Agras T100"}, true, "RID000002576"},
    {"1748CH", {"Autel EVO II"}, true, "RID000000030"},
    {"1748FEV2H", {"Autel EVO II V3"}, true, "RID000001969"},
    {"1748CLTC", {"Autel EVO Lite", "Autel EVO Lite+"}, true, "RID000000032"},
    {"1748FEC1LTC", {"Autel EVO Lite"}, true, "RID000001086"},
    {"1748FEC2LTC", {"Autel EVO Lite"}, true, "RID000001086"},
    {"1748FEV5LTE", {"Autel EVO Lite"}, true, "RID000001086"},
    {"1748CNA", {"Autel EVO Nano"}, false, "RID000000320"},
    {"1748FEC3NAC", {"Autel EVO Nano"}, true, "RID000001038"},
    {"1748FEC4NAC", {"Autel EVO Nano"}, true, "RID000001038"},
    {"1748FEV3HMA", {"Autel EVO Max"}, true, "RID000000151"},
    {"1748FEV3HMM", {"Autel MDX-1"}, true, "RID000002344"},
    {"1748FEV3HMP", {"Autel MDX-1"}, true, "RID000002344"},
    {"1748FEV3HMQ", {"Autel MDX-1"}, true, "RID000002344"},
    {"1748FEV4AL", {"Autel Alpha"}, true, "RID000001332"},
    {"1748CJD1", {"Autel Dragonfish Standard"}, true, "RID000000048"},
    {"1748CJD2", {"Autel Dragonfish Lite"}, true, "RID000000046"},
    {"1748CJD3", {"Autel Dragonfish Pro"}, true, "RID000000047"},
    {"1748FDF1JD1", {"Autel Dragonfish Standard"}, true, "RID000000048"},
    {"1748FDF2JD2", {"Autel Dragonfish Lite"}, true, "RID000000046"},
    {"1748FDF3JD3", {"Autel Dragonfish Pro"}, true, "RID000000047"},
    {"1588E040452", {"Parrot ANAFI Ai"}, true, "RID000000067"},
    {"1588E040505", {"Parrot ANAFI USA"}, true, "RID000000068"},
    {"1588E040465", {"Parrot ANAFI UKR GOV"}, true, "RID000002663"},
    {"1668BE10JA00", {"Skydio X2E"}, true, "RID000000097"},
    {"1668BS20", {"Skydio 2+"}, true, "RID000000099"},
    {"1668BR40", {"Skydio X10 family"}, false, "RID000000749"},
    {"1668BR70", {"Skydio R10"}, false, "RID000002670"},
    {"1907F7USE23B", {"Anzu Raptor"}, false, "RID000001035"},
    {"1907F7UEE23B", {"Anzu Raptor T"}, false, "RID000001036"},
    {"1910F14", {"Potensic ATOM 2"}, true, "RID000002237"},
    {"1910F916", {"Potensic RID-916 module"}, true, "RID000001125"},
    {"1865AHS600", {"Holy Stone HS600"}, false, "RID000000327"},
    {"1865CHS600D", {"Holy Stone HS600D"}, false, "RID000001175"},
    {"1865CHS720G0", {"Holy Stone HS720G"}, false, "RID000000411"},
    {"1865CHS920", {"Holy Stone HS920"}, false, "RID000001150"},
    {"18699F11", {"Ruko F11"}, false, "RID000000287"},
    {"1869AF11P", {"Ruko F11PRO"}, false, "RID000000288"},
    {"1869BF11P2", {"Ruko F11PRO 2"}, false, "RID000000683"},
    {"1869BF11GM", {"Ruko F11GIM"}, false, "RID000000521"},
    {"1869CF11GM2", {"Ruko F11GIM2"}, false, "RID000000251"},
    {"1869AU11P", {"Ruko U11PRO"}, false, "RID000000286"},
    {"1869AU11S", {"Ruko U11S"}, false, "RID000000285"},
    {"1892AF110", {"SJRC F11 4K PRO"}, false, "RID000000414"},
    {"1892CF11S23", {"SJRC F11s 4K PRO"}, false, "RID000000394"},
    {"1892CF22A23", {"SJRC F22 4K PRO"}, false, "RID000000492"},
    {"1892CF22B23", {"SJRC F22s 4K PRO"}, false, "RID000000379"},
    {"1892CF5SA00", {"SJRC F5s PRO"}, false, "RID000000473"},
    {"1892BF7A238", {"SJRC F7 4K PRO"}, false, "RID000000413"},
    {"1875BF7GB2", {"Bwine F7GB2"}, false, "RID000000418"},
    {"1875CF7GB3", {"Bwine F7GB3"}, true, "RID000002596"},
    {"1875CF7GM2", {"Bwine F7GIM2"}, false, "RID000001711"},
    {"2095FH05", {"Zero Zero Robotics ZZ-H-1-003"}, true, "RID000002593"},
    {"2095FH04", {"Zero Zero Robotics ZZ-H-1-004"}, true, "RID000002594"},
    {"17287200", {"Sony Airpeak S1"}, true, "RID000000017"},
    {"1767CSRH", {"Teledyne FLIR SkyRanger / SkyRaider"}, false, "RID000000191"},
    {"1930FD02", {"SwellPro FD2"}, true, "RID000001501"},
    {"1930FD03", {"SwellPro FD3"}, true, "RID000001762"},
    {"1930FARD", {"SwellPro A-RID module"}, true, "RID000001321"},
    {"1863FPEG8218", {"XAG P100"}, true, "RID000000162"},
    {"1863FPEG8219", {"XAG P100 Pro"}, true, "RID000000163"},
    {"1863FPEG8223", {"XAG P150"}, true, "RID000002240"},
    {"1863F251", {"XAG P150 Max"}, true, "RID000002602"},
    {"1863F242", {"XAG P60"}, true, "RID000002588"}};
const std::string kSerialChars = "0123456789ABCDEFGHJKLMNPQRSTUVWXYZ";   // CTA-2063-A excludes O and I

bool ends_with(const std::string& s, const std::string& t) {
  return s.size() >= t.size() && s.compare(s.size() - t.size(), t.size(), t) == 0;
}
}  // namespace

std::string dji_product(int64_t code) {
  const auto it = kDji.find(code);
  return it == kDji.end() ? std::string() : it->second;
}

std::string dji_product_source(int64_t code) {
  return code <= kDjiConfirmed ? "public decoder tables (RUB-SysSec / Kismet), real decodes"
                               : "DJI AeroScope type map (not yet confirmed by a real decode)";
}

Json serial_model(const std::string& serial) {
  const SerialRow* best = nullptr;   // python max(..., key=len): the first of the longest
  for (const auto& r : kSerial) {
    const std::string p = r.prefix;
    if (serial.compare(0, p.size(), p) == 0 && serial.size() >= p.size() &&
        (!best || p.size() > std::string(best->prefix).size()))
      best = &r;
  }
  if (!best) return nullptr;
  Json models = Json::array();
  for (const char* m : best->models) models.push_back(m);
  return {{"models", models}, {"verified", best->verified}, {"doc", best->doc}, {"prefix", best->prefix}};
}

Json serial_manufacturer(const std::string& serial) {
  if (serial.size() < 6) return nullptr;
  for (char c : serial)
    if (kSerialChars.find(c) == std::string::npos) return nullptr;
  const std::string lens = "123456789ABCDEF";
  const size_t p = lens.find(serial[4]);
  const size_t n = p == std::string::npos ? 0 : p + 1;
  if (n == 0 || serial.size() != 5 + n) return nullptr;
  const auto it = kMfr.find(serial.substr(0, 4));
  if (it == kMfr.end()) return nullptr;
  return {{"code", serial.substr(0, 4)}, {"maker", it->second.maker}, {"kind", it->second.kind}};
}

Json model_block(const std::string& family, const std::vector<std::string>& candidates,
                 const std::vector<std::string>& evidence, const std::string& basis,
                 const std::string& unresolved_reason) {
  std::vector<std::string> cands;                      // python list(dict.fromkeys(candidates)): first-seen order
  std::set<std::string> seen;
  for (const auto& c : candidates)
    if (seen.insert(c).second) cands.push_back(c);
  Json model = nullptr;
  if (cands.size() == 1) {
    std::string m = cands[0];
    if (basis != "decoded" && !ends_with(m, "-like")) m += "-like";
    model = m;
  }
  Json out = Json::object();
  out["family"] = family;
  out["model"] = model;
  out["candidates"] = cands;
  out["basis"] = basis;
  out["evidence"] = evidence;
  if (model.is_null())
    out["unresolved_reason"] = !unresolved_reason.empty() ? unresolved_reason
                               : cands.empty() ? "no model candidates" : "several models share these on-air features";
  return out;
}

Json remote_id_model(const Json& rid) {
  const Json mf = rid.contains("manufacturer") ? rid["manufacturer"] : Json(nullptr);
  const std::string uas = rid.contains("uas_id") && rid["uas_id"].is_string() ? rid["uas_id"].get<std::string>() : "";
  std::vector<std::string> ev;
  if (rid.contains("ua_type") && !rid["ua_type"].is_null())
    ev.push_back("Remote ID UA type: " + (rid["ua_type"].is_string() ? rid["ua_type"].get<std::string>()
                                                                     : rid["ua_type"].dump()));
  if (mf.is_object() && mf["kind"] == "aircraft") {
    ev.insert(ev.begin(), "Remote ID serial " + uas + ": CTA-2063-A manufacturer code " + mf["code"].get<std::string>() +
                              " = " + mf["maker"].get<std::string>());
    const Json sm = serial_model(uas);
    if (sm.is_object()) {
      std::vector<std::string> models;
      std::string joined;
      for (const auto& m : sm["models"]) models.push_back(m.get<std::string>()), joined += (joined.empty() ? "" : " / ") + m.get<std::string>();
      const bool ver = sm["verified"].get<bool>();
      ev.insert(ev.begin() + 1, "serial prefix " + sm["prefix"].get<std::string>() + " -> " + joined + " (FAA DOC " +
                                    sm["doc"].get<std::string>() +
                                    (ver ? "" : "; prefix inferred from the listed serials, unverified") + ")");
      return model_block(mf["maker"].get<std::string>(), models, ev, ver ? "decoded" : "inferred",
                         "the DOC prefix is shared by several models");
    }
    return model_block(mf["maker"].get<std::string>(), {}, ev, "decoded", "no FAA DOC serial prefix for this serial");
  }
  if (mf.is_object()) {
    ev.insert(ev.begin(), "Remote ID serial " + uas + ": add-on broadcast module by " + mf["maker"].get<std::string>() +
                              " (code " + mf["code"].get<std::string>() + ")");
    return model_block("unknown airframe (Remote ID module)", {}, ev, "decoded",
                       "an add-on module names the module maker, not the airframe");
  }
  if (!uas.empty()) ev.insert(ev.begin(), "Remote ID ID " + uas + ": no known manufacturer code");
  return model_block("unknown (Remote ID)", {}, ev, "decoded", "no manufacturer code in the Remote ID");
}

}  // namespace si::drone_models
