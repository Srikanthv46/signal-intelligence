// Pipeline plugin: LTE cells as capture-level signals (python: receivers/lte/plugin.py).
// Publishes the receiver ppm it uses in ctx.cache["ppm"] (Json, null if unknown) for nr/gsm.
#include <algorithm>
#include <set>

#include "core/bands_india.hpp"
#include "core/pyfmt.hpp"
#include "decode/plugin.hpp"
#include "receivers/lte/cellsearch.hpp"

namespace si::rx::lte {

namespace {

bool in_lte_band(const std::vector<std::string>& keys) {
  for (const auto& [k, v] : lte_bands())
    if (std::find(keys.begin(), keys.end(), k) != keys.end()) return true;
  return false;
}

// energy SNR of a candidate: the PSS correlates over many symbols, so a carrier only a few dB over the floor in
// energy still gives a decisive PSS/SSS (python CAND_SNR_DB)
constexpr double kCandSnrDb = 3.0;

// a Python dict is truthy when not None and not empty
bool truthy(const Json* j) { return j && j->is_object() && !j->empty(); }

class LteDecoder : public Decoder {
 public:
  std::string name() const override { return "lte"; }
  std::string package() const override { return "lte"; }
  std::string protocol() const override { return "3GPP LTE (E-UTRA) downlink"; }
  bool capture_level() const override { return true; }

  bool can_decode(const Signal& sig, DecodeContext&) const override {
    return in_lte_band(bands::band_keys(sig.f0, sig.f1));
  }

  // python _ppm(ctx, cand)
  Json ppm_cal(DecodeContext& ctx, const std::vector<const Signal*>& cand) const {
    if (ctx.cache.contains("ppm")) {
      const Json* p = ctx.cache.get<Json>("ppm");
      return p ? *p : Json(nullptr);
    }
    const Json* sw = ctx.cache.get<Json>("sweep_ppm");
    Json cal = sw ? *sw : Json(nullptr);
    if (cal.is_null() && !cand.empty()) {
      std::vector<const Signal*> top = cand;
      std::stable_sort(top.begin(), top.end(), [](const Signal* a, const Signal* b) { return -a->snr_db < -b->snr_db; });
      if (top.size() > 3) top.resize(3);
      dvec centres;
      for (const Signal* s : top) centres.push_back(static_cast<double>(py_round_int(s->fc() / 100e3)) * 100e3);
      auto w = wide_ppm_search(ctx.ch, centres);
      if (w) {
        cal = *w;
        cal["source"] = "this capture: strongest PSS, +-500 kHz search";
      }
    }
    ctx.cache.set<Json>("ppm", cal);
    return cal;
  }

  // python capture_signals: continuous carriers >= 6 dB (energy) are searched as always; bursty or weaker energy
  // (>= kCandSnrDb) too, widest first, skipping rasters already searched, and there only SSS-resolved cells count.
  std::vector<Signal> capture_signals(DecodeContext& ctx) const override {
    std::vector<const Signal*> inb, base, extra;
    for (const Signal& s : ctx.det->signals) {
      if (!s.flags.empty() || !(s.snr_db >= kCandSnrDb) || !(s.bw() >= 1.0e6) ||
          !in_lte_band(bands::band_keys(s.f0, s.f1)))
        continue;
      inb.push_back(&s);
      bool cont = false;
      for (const Box& b : s.bursts) cont = cont || b.kind == "continuous";
      (s.snr_db >= 6 && cont ? base : extra).push_back(&s);
    }
    if (inb.empty()) return {};
    std::stable_sort(extra.begin(), extra.end(), [](const Signal* a, const Signal* b) { return -a->bw() < -b->bw(); });
    const Json cal = ppm_cal(ctx, base.empty() ? inb : base);
    const bool have = truthy(&cal);
    const double ppm = have ? cal["ppm"].get<double>() : 0.0;
    std::vector<Signal> out;
    std::set<std::pair<int64_t, int>> seen;
    std::set<int64_t> done;
    std::vector<std::pair<const Signal*, bool>> order;
    for (const Signal* s : base) order.emplace_back(s, false);
    for (const Signal* s : extra) order.emplace_back(s, true);
    for (const auto& [sig, strict] : order) {
      const std::vector<Cell> cells = search_cells(ctx.ch, *sig, ppm, strict ? &done : nullptr, have, strict);
      if (!strict)
        for (double r : raster_grid(*sig)) done.insert(py_round_int(r / 100e3));
      for (const Cell& c : cells) {
        if (band_of(c.raster_hz) == "?") continue;  // an LTE downlink carrier lies in an LTE band
        const auto key = std::make_pair(py_round_int(c.raster_hz / 100e3), c.nid2);
        if (seen.count(key)) continue;
        seen.insert(key);
        out.push_back(make_signal(ctx, *sig, c, have ? &cal : nullptr));
      }
    }
    return out;
  }

  // python _band(raster)
  static std::string band_of(double raster) {
    const auto ea = bands::lte_earfcn_dl(raster);
    if (ea) return "B" + std::to_string(ea->first);
    const auto keys = bands::band_keys(raster, raster);
    for (const auto& [k, v] : lte_bands())
      if (std::find(keys.begin(), keys.end(), k) != keys.end()) return v;
    return "?";
  }

  Signal make_signal(DecodeContext& ctx, const Signal& sig, const Cell& c, const Json* cal) const {
    std::vector<std::string> steps;
    if (cal) {
      steps.push_back(fmt("carrier hypotheses: 100 kHz EARFCN raster inside %.2f-%.2f MHz, receiver error %+.2f ppm "
                          "applied (",
                          sig.f0 / 1e6, sig.f1 / 1e6, (*cal)["ppm"].get<double>()) +
                      py_str((*cal)["source"]) + ")");
    } else {
      steps.push_back("carrier hypotheses: 100 kHz EARFCN raster (receiver ppm unknown, 0 assumed)");
    }
    steps.push_back(
        "PSS: 3.84 MSps slice, 3 roots, 0 Hz (receiver error known, < 3 kHz residual) or +-5 kHz, normalised "
        "correlation folded on 5 ms; z >= 15");
    steps.push_back(fmt("PSS N_ID2 %d z %.1f over %zu half-frames; CFO refined +-3 kHz", c.nid2, c.pss_z,
                        c.inst.size()));
    steps.push_back(
        "SSS: 168 N_ID1 x {SF0,SF5} at FDD/TDD x normal/extended CP, maximum-ratio combined with the PSS channel "
        "estimate (smoothed over 225 kHz), accumulated over half-frames");
    const double raster = c.raster_hz;
    const auto ea = bands::lte_earfcn_dl(raster);
    const std::string band = band_of(raster);
    const auto bw = bandwidth_from_energy(*ctx.det, c.carrier_hz);
    const double half = (bw ? occ(*bw) : 1.08e6) / 2;
    Json fields = Json::object();
    fields["carrier_hz"] = py_round_int(c.carrier_hz);
    fields["raster_carrier_hz"] = py_round_int(raster);
    fields["earfcn"] = ea ? Json(ea->second) : Json(nullptr);
    fields["band"] = band;
    fields["nid2"] = c.nid2;
    fields["pss_fold_z"] = py_round(c.pss_z, 1);
    fields["n_half_frames"] = c.inst.size();
    fields["pss_isolation"] = py_round(c.pss_isolation, 2);
    fields["receiver_freq_error_hz"] = py_round_int(c.err_hz);
    fields["receiver_freq_error_ppm"] = py_round(c.err_hz / raster * 1e6, 2);
    fields["channel_bandwidth_hz_from_energy"] = bw ? Json(*bw) : Json(nullptr);
    fields["mib"] = "not decoded (PBCH decoder not implemented)";
    const auto& s = c.sss;
    const bool ok = sss_decisive(s);
    const std::pair<double, double> span{c.carrier_hz - half, c.carrier_hz + half};
    DecodeResult res;
    if (ok) {
      const int pci = 3 * s->nid1 + c.nid2;
      // subframe 0 among the PSS instances found (none when they all fall on subframe-5 half-frames)
      Json sf0 = nullptr;
      for (const auto& [k, p] : c.inst)
        if ((k + s->parity) % 2 == 0) {
          sf0 = py_round(static_cast<double>(p) / kOfs, 6);
          break;
        }
      fields["pci"] = pci;
      fields["nid1"] = s->nid1;
      fields["duplex"] = s->duplex;
      fields["cyclic_prefix"] = s->cp;
      fields["sss_z"] = py_round(s->z, 1);
      fields["sss_z_runner_up"] = py_round(s->z2, 1);
      fields["subframe0_time_s"] = sf0;
      steps.push_back(fmt("SSS N_ID1 %d (%s, %s CP): z %.1f vs runner-up %.1f -> PCI %d", s->nid1, s->duplex.c_str(),
                          s->cp.c_str(), s->z, s->z2, pci));
      std::string lab = fmt("LTE %s %s downlink PCI %d", band.c_str(), s->duplex.c_str(), pci) +
                        (ea ? fmt(" EARFCN %d", ea->second) : std::string()) +
                        (bw ? fmt(" %g MHz", *bw / 1e6) : std::string());
      std::string ev = fmt("LTE PSS (root %d) every 5 ms z %.1f; SSS N_ID1 %d z %.1f vs %.1f; carrier on EARFCN raster "
                           "within %+.1f kHz",
                           kRoots[c.nid2], c.pss_z, s->nid1, s->z, s->z2, c.err_hz / 1e3);
      res = result("partial", "LTE PSS/SSS cell search (this package)", steps, fields,
                   "physical cell identity decoded; MIB/SIB (PBCH/PDSCH) decoding not implemented", lab, 0.95,
                   std::string("OFDMA 15 kHz SCS (PSS/SSS sync)"), ev, span);
    } else {
      steps.push_back("SSS not decisive" + (s ? fmt(" (best z %.1f vs %.1f)", s->z, s->z2) : std::string()));
      res = result("partial", "LTE PSS/SSS cell search (this package)", steps, fields,
                   "LTE PSS found but SSS ambiguous: N_ID1 not resolved",
                   fmt("LTE %s downlink (PSS N_ID2 %d)", band.c_str(), c.nid2), 0.85, std::string("OFDMA 15 kHz SCS"),
                   fmt("LTE PSS root %d every 5 ms (z %.1f)", kRoots[c.nid2], c.pss_z), span);
    }
    const bool tdd = ok && s->duplex == "TDD";
    std::vector<std::pair<double, double>> bursts;
    if (tdd) {
      for (const Box& b : sig.bursts) bursts.emplace_back(b.t0, b.t1);
    } else {
      bursts.emplace_back(0.0, ctx.cap.duration());
    }
    Signal o = receiver_signal(bursts, span.first, span.second, res, sig.snr_db, name());
    for (Box& b : o.bursts) b.kind = tdd ? "burst" : "continuous";
    o.features["owns_band"] = true;
    return o;
  }
};

}  // namespace

void register_plugins(PluginList& out) { out.push_back(std::make_unique<LteDecoder>()); }

}  // namespace si::rx::lte
