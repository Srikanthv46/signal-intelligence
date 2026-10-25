// python receivers/digital_vtx/catalog.py (same table, same order, same strings).
#include "receivers/digital_vtx/catalog.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "core/drone_models.hpp"
#include "core/pyfmt.hpp"

namespace si::rx::digital_vtx::catalog {

namespace {

struct Band { double lo, hi; };
struct System {
  const char* name;
  double tu;  // 0 = None
  std::vector<double> cp, bw;
  std::vector<Band> bands;
  int raster;  // 0 none, 1 dji_x5, 2 1mhz, 3 list
  std::vector<double> channels;
  int antennas;  // 0 = None
  double frame_ms;  // 0 = None
  const char* src;
};

std::vector<double> cat(std::vector<double> a, const std::vector<double>& b) { a.insert(a.end(), b.begin(), b.end()); return a; }

const std::vector<System>& systems() {
  static const std::vector<double> race = {5658, 5695, 5732, 5769, 5806, 5843, 5880, 5917};
  static const std::vector<double> djifpv = {5660, 5695, 5735, 5770, 5805, 5839, 5878, 5914};
  static std::vector<double> ar_extra = [] {
    std::vector<double> v = {5758, 5788, 5828};
    for (int k = 0; k < 8; ++k) v.push_back(5340 + 40 * k);
    return v;
  }();
  static const std::vector<double> hdz_l = [] {   // HDZero L1-L8
    std::vector<double> v;
    for (int k = 0; k < 8; ++k) v.push_back(5362 + 37 * k);
    return v;
  }();
  const Band ism24{2400, 2484}, ism58{5650, 5950}, u51{5150, 5350}, u58l{5340, 5650};
  static const std::vector<System> s = {
      {"DJI O4 (Mini 4 Pro, Air 3/3S, Avata 2, Neo/Flip, O4 Air Unit Pro)", 33.33, {0.0703}, {9.0, 18.0, 36.0, 54.0},
       {ism24, u51, ism58}, 1, {}, 2, 2.5, "M-lib + DJI specs (2T4R drones, O4 Air Unit Pro 2T2R)"},
      {"DJI O4 Air Unit (Lite, 1T1R)", 33.33, {0.0703}, {9.0, 18.0, 36.0, 54.0}, {u51, ism58}, 1, {}, 1, 2.5,
       "DJI O4 Air Unit specs (1T1R)"},
      {"DJI O3 / O3+ (Air 2S, Mini 3/3 Pro, Mavic 3, Avata, O3 Air Unit)", 66.67, {0.0703}, {9.0, 18.0, 36.0},
       {ism24, ism58}, 1, {}, 2, 0,
       "M-3p (Air 2S: 2048 @ 30.72 MSps, CP 144/160, SFBC), O3 Air Unit channels 5768.5 / 5804.5 / 5839.5"},
      {"DJI OcuSync 2 (Mini 2, Mini SE 2, Air 2, Mavic 2)", 66.67, {0.0703, 0.0781}, {9.0, 18.0}, {ism24, ism58}, 0, {},
       0, 0, "M-lib (Mini 2)"},
      {"DJI OcuSync 1 (Mavic Pro)", 66.67, {0.0703}, {1.3, 9.0, 18.0}, {ism24}, 2, {}, 0, 0,
       "DJI SDK (1.4 / 10 / 20 MHz), 1 MHz raster"},
      {"Herelink (CubePilot)", 66.67, {0.0703}, {18.0, 9.0}, {ism24}, 0, {}, 0, 5.0,
       "M-lib (LTE PSS/SSS; the herelink path names it when its PSS is found)"},
      {"Autel SkyLink (EVO II / Lite+ / Nano / Max)", 66.67, {0.0703}, {9.0, 18.0}, {ism24, ism58, {900, 930}}, 0, {}, 0,
       0, "vendor C-UAS white paper (600 sc); RFUAV (Nano ~20 MHz)"},
      {"Walksnail Avatar HD (dual-antenna air unit: Pro / V2 / Moonlight / GT)", 91.43, {0.125}, {18.9}, {ism58, u51}, 3,
       cat(djifpv, race), 2, 7.163, "waveshark #226 (802.16e 20 MHz, 1729 sc, STBC matrix A, frame 7.163 ms)"},
      {"Walksnail Avatar HD (single-antenna air unit: Mini 1S / Mini 1S Lite / Nano)", 91.43, {0.125}, {18.9},
       {ism58, u51}, 3, cat(djifpv, race), 1, 7.163, "waveshark #226; Caddx: Mini 1S / Lite / Nano V3 single antenna"},
      {"Artosyn AR8030 link (Walksnail Ascent, BETAFPV ArtLynk, HGLRC Draco, Skyzone SteadyDigital, StartRC VT5)", 91.43,
       {0.125}, {4.7, 9.4, 18.9}, {ism58}, 3, cat(race, ar_extra), 0, 0,
       "waveshark #228 / fpvwiki (5 / 10 / 20 MHz, 19 centres)"},
      {"DJI Lightbridge 2 (Phantom 4 / Inspire 2 / M600)", 91.43, {0.125}, {9.2, 9.6}, {ism24, ism58}, 0, {}, 0, 0,
       "M-lib (one Phantom 4 capture)"},
      {"HDZero / Shark Byte (Divimath)", 37.926, {0.03125}, {23.44, 17.0}, {ism58, u58l}, 3,
       cat(cat(race, {5760, 5800, 5705, 5740}), hdz_l), 1, 16.67,
       "measured on an HDZero Freestyle V2 (27 MHz mode: 1024 FFT at 27 MSps, CP 32, 889 sc, 60 fps frame); "
       "17 MHz mode and Shark Byte not measured"},
      {"DJI O1 (DJI FPV system, Air Unit, Caddx Vista, Polar, Runcam Link)", 0, {}, {18.0, 36.0}, {ism58}, 3, djifpv, 0, 0,
       "DJI FCC filings (20 / 40 MHz)"},
      {"Cambium PMP 450 (fixed-wireless infrastructure, not a drone)", 22.86, {0.0625}, {18.5}, {ism58, u51}, 0, {}, 0,
       2.5, "M-lib"},
  };
  return s;
}

bool raster_ok(double f, const System& s) {
  if (s.raster == 1) return std::fabs((f - 0.5) - std::nearbyint(f - 0.5)) <= 0.15;
  if (s.raster == 2) return std::fabs(f - std::nearbyint(f)) <= 0.15;
  double d = 1e300;
  for (double c : s.channels) d = std::min(d, std::fabs(f - c));
  return d <= 1.0;
}

std::string num(double v) { return py_str(Json(v)); }  // python str(float)

}  // namespace

Json candidates(double tu_us, std::optional<double> cp_frac, std::optional<double> bw_mhz, std::optional<double> carrier_mhz,
                std::optional<int> tx_antennas, std::optional<double> frame_ms, size_t top) {
  struct C { Json j; int score; size_t nconf; };
  std::vector<C> out;
  for (const auto& s : systems()) {
    Json m = Json::array(), c = Json::array();
    if (s.tu != 0) {
      if (std::fabs(tu_us - s.tu) > 0.015 * s.tu + 0.1) continue;
      m.push_back("Tu " + num(s.tu) + " us");
    }
    if (!s.cp.empty() && cp_frac) {
      double d = 1e300;
      for (double x : s.cp) d = std::min(d, std::fabs(*cp_frac - x));
      (d <= 0.008 ? m : c).push_back(fmt("CP %.4f Tu", *cp_frac));
    }
    if (bw_mhz) {
      double near = s.bw[0];
      for (double b : s.bw)
        if (std::fabs(*bw_mhz - b) < std::fabs(*bw_mhz - near)) near = b;
      (std::fabs(*bw_mhz - near) <= 0.12 * near ? m : c).push_back(fmt("BW %.1f MHz (mode ", *bw_mhz) + num(near) + ")");
    }
    if (carrier_mhz) {
      bool inb = false;
      for (const auto& b : s.bands) inb = inb || (b.lo <= *carrier_mhz && *carrier_mhz <= b.hi);
      (inb ? m : c).push_back("band");
      if (s.raster != 0) (raster_ok(*carrier_mhz, s) ? m : c).push_back("channel raster");
    }
    if (tx_antennas && s.antennas != 0) {
      if (*tx_antennas == 2 && s.antennas == 1) c.push_back("2 TX antennas measured, system has 1");
      else if (*tx_antennas == 2 && s.antennas == 2) m.push_back("2 TX antennas");
    }
    if (frame_ms && s.frame_ms != 0 && std::fabs(*frame_ms - s.frame_ms) <= 0.05 * s.frame_ms)
      m.push_back("frame " + num(s.frame_ms) + " ms");
    if (s.tu == 0 && m.size() < 2) continue;
    const int score = static_cast<int>(m.size()) - static_cast<int>(c.size());
    out.push_back({Json{{"system", s.name}, {"score", score}, {"matched", m}, {"conflicts", c}, {"source", s.src}}, score, c.size()});
  }
  std::stable_sort(out.begin(), out.end(), [](const C& a, const C& b) {
    return a.score != b.score ? a.score > b.score : a.nconf < b.nconf;
  });
  Json res = Json::array();
  for (size_t i = 0; i < out.size() && i < top; ++i) res.push_back(out[i].j);
  return res;
}

namespace {
struct ModelEntry { const char* system; const char* family; std::vector<const char*> models; };
// python MODELS (file 11 decision tables)
const std::vector<ModelEntry>& model_table() {
  static const std::vector<ModelEntry> t = {
      {"DJI O4 (Mini 4 Pro, Air 3/3S, Avata 2, Neo/Flip, O4 Air Unit Pro)", "DJI O4",
       {"Mini 4 Pro", "Air 3", "Air 3S", "Avata 2", "Neo", "Flip", "Matrice 4", "O4 Air Unit Pro"}},
      {"DJI O4 Air Unit (Lite, 1T1R)", "DJI O4", {"O4 Air Unit (Lite)"}},
      {"DJI O3 / O3+ (Air 2S, Mini 3/3 Pro, Mavic 3, Avata, O3 Air Unit)", "DJI O3",
       {"Air 2S", "Mini 3 Pro", "Mavic 3 / 3 Classic / 3 Pro", "Avata", "Matrice 30", "O3 Air Unit"}},
      {"DJI OcuSync 2 (Mini 2, Mini SE 2, Air 2, Mavic 2)", "DJI OcuSync 2",
       {"Mini 2 / Mini 4K / Mini 2 SE", "Mini 3", "Mavic Air 2", "Mavic 2 Pro / Zoom"}},
      {"DJI OcuSync 1 (Mavic Pro)", "DJI OcuSync 1", {"Mavic Pro"}},
      {"Herelink (CubePilot)", "Herelink", {"Herelink (v1.1 / v2: identical on air)"}},
      {"Autel SkyLink (EVO II / Lite+ / Nano / Max)", "Autel SkyLink", {"EVO II", "EVO Lite+", "EVO Nano", "EVO Max"}},
      {"Walksnail Avatar HD (dual-antenna air unit: Pro / V2 / Moonlight / GT)", "Walksnail Avatar HD",
       {"Avatar HD V1 / Pro", "Avatar HD V2 (dual antenna)", "Avatar GT", "Moonlight"}},
      {"Walksnail Avatar HD (single-antenna air unit: Mini 1S / Mini 1S Lite / Nano)", "Walksnail Avatar HD",
       {"Mini 1S", "Mini 1S Lite", "Nano V3", "Avatar HD V2 (single antenna)", "Avatar GT2"}},
      {"Artosyn AR8030 link (Walksnail Ascent, BETAFPV ArtLynk, HGLRC Draco, Skyzone SteadyDigital, StartRC VT5)",
       "Artosyn AR8030", {"Walksnail Ascent", "BETAFPV ArtLynk", "HGLRC Draco", "Skyzone SteadyDigital", "StartRC VT5"}},
      {"DJI Lightbridge 2 (Phantom 4 / Inspire 2 / M600)", "DJI Lightbridge 2", {"Phantom 4 Pro", "Inspire 2", "Matrice 600"}},
      {"HDZero / Shark Byte (Divimath)", "HDZero / Shark Byte",
       {"HDZero VTX (any: Race / Freestyle / Whoop / ECO / AIO)", "Shark Byte TX5S.1 / TX5M.1"}},
      {"DJI O1 (DJI FPV system, Air Unit, Caddx Vista, Polar, Runcam Link)", "DJI O1 (DJI FPV system)",
       {"DJI FPV Air Unit", "DJI FPV drone", "Caddx Vista", "RunCam Link"}}};
  return t;
}
const ModelEntry* find_model(const std::string& sys) {
  for (const auto& m : model_table())
    if (sys == m.system) return &m;
  return nullptr;
}
}  // namespace

Json model_candidates(const Json& system_candidates, std::optional<double> carrier_mhz, std::optional<double> bw_mhz,
                      std::optional<int> tx_antennas) {
  std::vector<const Json*> ok;
  if (system_candidates.is_array())
    for (const auto& c : system_candidates)
      if (c["conflicts"].empty() && find_model(c["system"].get<std::string>())) ok.push_back(&c);
  if (ok.empty()) return nullptr;
  int64_t best = (*ok[0])["score"].get<int64_t>();
  for (const Json* c : ok) best = std::max(best, (*c)["score"].get<int64_t>());
  std::vector<const ModelEntry*> syst;
  std::vector<std::string> ev, fams, models;
  for (const Json* c : ok)
    if ((*c)["score"].get<int64_t>() == best) {
      const ModelEntry* m = find_model((*c)["system"].get<std::string>());
      syst.push_back(m);
      ev.push_back(std::string("system ") + m->system + " (research catalog)");
      fams.push_back(m->family);
      for (const char* x : m->models) models.push_back(x);
    }
  const bool dji = std::any_of(fams.begin(), fams.end(), [](const std::string& f) { return f == "DJI O3" || f == "DJI O4"; });
  if (carrier_mhz && dji) {
    auto near = [&](std::initializer_list<double> cs) {
      double d = 1e300;
      for (double c : cs) d = std::min(d, std::fabs(*carrier_mhz - c));
      return d <= 0.3;
    };
    std::string want;
    if (near({5789.5, 5814.5})) {
      want = "DJI O4";
      ev.push_back(fmt("channel %.1f MHz is an O4 centre (not used by O3)", *carrier_mhz));
    } else if (near({5804.5, 5839.5})) {
      want = "DJI O3";
      ev.push_back(fmt("channel %.1f MHz is an O3 centre (not used by O4)", *carrier_mhz));
    }
    if (!want.empty() && std::find(fams.begin(), fams.end(), want) != fams.end()) {
      std::vector<std::string> f2, m2;
      for (size_t i = 0; i < syst.size(); ++i)
        if (fams[i] == want) {
          f2.push_back(fams[i]);
          for (const char* x : syst[i]->models) m2.push_back(x);
        }
      fams = f2;
      models = m2;
    }
  }
  if (std::find(fams.begin(), fams.end(), "HDZero / Shark Byte") != fams.end() && bw_mhz && std::fabs(*bw_mhz - 17.0) <= 2.0) {
    for (auto& m : models)
      if (m.rfind("HDZero", 0) == 0) m = "HDZero VTX with the RunCam Nano 90 camera (Whoop / Lite / AIO)";
    ev.push_back("17 MHz mode: only with the RunCam Nano 90 camera");
  }
  if (tx_antennas && *tx_antennas == 2) ev.push_back("2 transmit antennas (STBC / SFBC)");
  std::string fam;
  std::vector<std::string> seen;
  for (const auto& f : fams)
    if (std::find(seen.begin(), seen.end(), f) == seen.end()) {
      fam += (seen.empty() ? "" : " / ") + f;
      seen.push_back(f);
    }
  return drone_models::model_block(fam, models, ev, "inferred");
}

Json hdzero_model(const std::string& bandwidth_mode) {
  static const char* kVtx[] = {"Race V1", "Race V2", "Race V3", "Freestyle V1", "Freestyle V2",
                               "Whoop", "Whoop Lite", "ECO", "AIO5", "AIO15"};
  std::vector<std::string> models;
  for (const char* m : kVtx) models.push_back(std::string("HDZero ") + m + " VTX");
  const std::vector<std::string> ev = {
      "HDZero (Divimath DM6300) pilot structure, " + bandwidth_mode + " mode, measured on an HDZero Freestyle V2",
      "every HDZero VTX transmits the same PHY; only the power ceiling differs (research file 11 section 2.3)",
      "Shark Byte (the Divimath predecessor) is not measured and not listed"};
  return drone_models::model_block("HDZero", models, ev, "inferred", "every HDZero VTX is identical on air (power only)");
}

Json hop_model(std::optional<double> hop_bw_mhz, std::optional<double> hop_us, std::optional<double> period_ms) {
  // python HOP_MODELS / HOP_RANGES (research, PLAN M6 / D3); p <= 0 = interval not used
  struct H { const char* model; const char* brand; double w, t, p; };
  static const H kHop[] = {
      {"DJI FPV", "DJI", 5.0, 0.64, 4.0}, {"DJI Avata 2", "DJI", 7.0, 0.41, 2.0},
      {"DJI Mavic 3 Pro", "DJI", 4.78, 1.7, 6.0}, {"DJI Mini 3", "DJI", 3.5, 0.56, 5.96},
      {"DJI Mini 4 Pro", "DJI", 6.54, 0.404, 2.5},
      {"DJI OcuSync 2 aircraft (Mini 2-class; model unresolved)", "DJI", 1.1, 0.51, 0.0},
      {"Herelink", "Herelink", 2.96, 0.52, 5.16}, {"Skydroid H12", "Skydroid", 6.06, 0.25, 2.97},
      {"Skydroid T10", "Skydroid", 5.27, 0.27, 2.91}, {"Siyi MK15", "Siyi", 3.32, 1.55, 4.99},
      {"Siyi MK32", "Siyi", 7.23, 0.47, 5.01}, {"Siyi FT24", "Siyi", 4.59, 1.61, 0.0}};
  struct R { const char* model; const char* brand; double w0, w1, t, p0, p1; };
  static const R kRange[] = {{"DJI O4 aircraft (model unresolved)", "DJI", 4.4, 7.0, 0.41, 2.0, 2.5}};
  if (!hop_bw_mhz || !hop_us || *hop_bw_mhz == 0 || *hop_us == 0) return nullptr;
  const double w = *hop_bw_mhz, on = *hop_us / 1e3;
  const bool hp = period_ms && *period_ms != 0;
  std::vector<std::string> names, brands;
  for (const auto& h : kHop)
    if (std::fabs(w - h.w) <= 0.15 * h.w && std::fabs(on - h.t) <= 0.2 * h.t &&
        (h.p <= 0 || (hp && std::fabs(*period_ms - h.p) <= 0.1 * h.p)))
      names.push_back(h.model), brands.push_back(h.brand);
  for (const auto& r : kRange)
    if (0.85 * r.w0 <= w && w <= 1.15 * r.w1 && std::fabs(on - r.t) <= 0.2 * r.t && hp && 0.9 * r.p0 <= *period_ms &&
        *period_ms <= 1.1 * r.p1)
      names.push_back(r.model), brands.push_back(r.brand);
  if (names.empty()) return nullptr;
  std::string list, fam;
  std::vector<std::string> seen;
  for (const auto& m : names) list += (list.empty() ? "" : ", ") + m;
  for (const auto& b : brands)
    if (std::find(seen.begin(), seen.end(), b) == seen.end()) fam += (seen.empty() ? "" : " / ") + b, seen.push_back(b);
  const std::string ev = "hops " + py_str(Json(w)) + " MHz x " + py_str(Json(*hop_us)) + " us every " +
                         (period_ms ? py_str(Json(*period_ms)) : std::string("None")) + " ms match the published timing of " +
                         list + " (one measured unit per model: low confidence)";
  return drone_models::model_block(fam + " (hop timing)", names, {ev}, "inferred");
}

}  // namespace si::rx::digital_vtx::catalog
