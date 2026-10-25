// Generic half of the digital_vtx detector (python: receivers/digital_vtx/detector.py,
// KNOWN_NUMEROLOGIES ... generic_detect): every CP-OFDM link with Tu >= 10 us, found blind
// (blind.cpp), one Emission per burst; the numerology table names what it resembles. Then the
// hopping narrowband CP-OFDM links (blind_hops), one Emission per hop.
#include <algorithm>
#include <mutex>
#include <cmath>

#include "core/bands_india.hpp"
#include "core/bursts.hpp"
#include "core/npcompat.hpp"
#include "core/pyfmt.hpp"
#include "core/spectral.hpp"
#include "receivers/digital_vtx/blind.hpp"
#include "receivers/digital_vtx/common.hpp"
#include "receivers/digital_vtx/detector.hpp"
#include "receivers/digital_vtx/catalog.hpp"
#include "receivers/digital_vtx/streams.hpp"

namespace si::rx::digital_vtx {

namespace {
struct KnownNum {
  const char* name;
  double tu_us;  // 0 = None
  std::vector<double> bw_mhz;
  const char* role;
};
// KNOWN_NUMEROLOGIES (name, tu_us, bw_mhz, role) in table order
const std::vector<KnownNum>& known() {
  static const std::vector<KnownNum> k = {
      {"wifi_he", 12.8, {4.4, 19.1, 38.3, 78.1}, "wifi"},
      {"herelink", 66.67, {18.0, 9.0}, "drone"},
      {"lte_like", 66.67, {1.08, 4.5, 9.0, 18.0}, "drone"},
      {"dji_like", 66.67, {1.4, 10.0, 20.0, 40.0}, "drone"},
      {"autel_skylink_like", 66.67, {9.0, 1.1}, "drone"},
      {"dji_o4_like", 33.33, {9.0, 18.0, 36.0}, "drone"},
      {"lightbridge_like", 91.43, {9.2}, "drone"},
      {"walksnail_like", 91.43, {18.9}, "drone"},
      {"pmp450_like", 22.86, {18.5}, "infrastructure"},
      {"hdzero_like", 37.926, {23.44}, "drone"},
  };
  return k;
}
}  // namespace

std::pair<std::string, std::vector<std::string>> classify_numerology(double tu_us, std::optional<double> bw_mhz) {
  std::vector<std::string> hits;
  for (const auto& n : known())
    if (n.tu_us != 0 && std::fabs(tu_us - n.tu_us) < 0.015 * n.tu_us + 0.1) hits.push_back(n.name);
  if (hits.empty()) return {"unknown", {}};
  if (std::find(hits.begin(), hits.end(), "lte_like") != hits.end()) {
    std::vector<std::string> h;
    for (const auto& x : hits)
      if (x != "herelink") h.push_back(x);
    return {"lte_like", h};
  }
  if (bw_mhz && hits.size() > 1) {
    size_t best = 0;
    double bd = 0;
    for (size_t i = 0; i < hits.size(); ++i) {
      const auto& n = *std::find_if(known().begin(), known().end(), [&](const KnownNum& k) { return hits[i] == k.name; });
      double d = std::fabs(*bw_mhz - n.bw_mhz[0]);
      for (double b : n.bw_mhz) d = std::min(d, std::fabs(*bw_mhz - b));
      if (i == 0 || d < bd) bd = d, best = i;
    }
    return {hits[best], hits};
  }
  return {hits[0], hits};
}

std::string numerology_role(const std::string& name) {
  for (const auto& n : known())
    if (name == n.name) return n.role;
  return "unknown";
}

namespace detail {

std::vector<Emission> generic_detect(const cfvec& x, double fs, double fc, const cdvec& X,
                                     const std::vector<Emission>& exclude) {
  const double usable = std::min(fs, 56e6) / 2;
  const double noise = noise_floor_psd(x).floor_power_per_sample;
  // python _SliceCache: the last 6 distinct slices are reused (scan -> confirm -> segment, links -> hops; PLAN P1)
  struct Entry { double f, h, r, ofs; cdvec y; };
  std::vector<Entry> cache;
  std::mutex mu;   // ofdm15k's scan calls the slicer from several threads: lookup / insert only under the lock
  const blind::Slicer sl = [&](double f, double h, double r, double* o) {
    {
      std::lock_guard<std::mutex> lk(mu);
      for (size_t i = 0; i < cache.size(); ++i)
        if (cache[i].f == f && cache[i].h == h && cache[i].r == r) {
          Entry e = std::move(cache[i]);
          cache.erase(cache.begin() + static_cast<std::ptrdiff_t>(i));
          cache.push_back(std::move(e));
          if (o) *o = cache.back().ofs;
          return cache.back().y;
        }
    }
    double ofs = 0;
    cdvec y = slice(X, fs, f, h, r, &ofs);
    if (o) *o = ofs;
    std::lock_guard<std::mutex> lk(mu);
    if (cache.size() >= 6) cache.erase(cache.begin());
    cache.push_back(Entry{f, h, r, ofs, y});
    return y;
  };
  const bursts::TFPower tf(x, fs);
  // a candidate inside a licensed cellular band, and whose IQ mirror is too, is the operator's carrier: the plugin
  // drops it, so it is not confirmed here (PLAN G8)
  const blind::SkipFn cellular = [fc](double off) { return bands::cellular_at(fc + off) && bands::cellular_at(fc - off); };
  auto links = blind::blind_links(sl, fs, fc, usable, noise, &tf, streams::SNR_MIN, cellular);
  std::vector<Emission> out;
  for (auto& lk : links) {
    const double lo = (fc + lk.lo) / 1e6, hi = (fc + lk.hi) / 1e6;
    const double tu_us = lk.tu * 1e6, ts_us = lk.ts * 1e6;
    {  // a link of the same Tu over more than half of whose band an exclude emission lies: Herelink (PSS, 66.67 us),
       // HDZero (pilots, 37.93 us)
      bool hit = false;
      for (const auto& e : exclude) {
        const double etu = e.features.contains("tu_us") && e.features["tu_us"].is_number() ? e.features["tu_us"].get<double>() : 66.667;
        hit = hit || (std::fabs(tu_us / etu - 1) < 0.03 && std::min(hi, e.f_hi_mhz) - std::max(lo, e.f_lo_mhz) > 0.5 * (hi - lo));
      }
      if (hit) continue;
    }
    const double cp_us = ts_us - tu_us;
    const double scs = 1e3 / tu_us;
    const auto [name, matches] = classify_numerology(tu_us, lk.bw / 1e6);
    const std::string role = numerology_role(name);
    const bool tn = lk.ts_np;
    Json feats = Json::object();
    feats["numerology"] = name;
    feats["matches"] = matches;
    feats["role"] = role;
    feats["tu_us"] = py_round(tu_us, 3);
    feats["symbol_us"] = rnd(ts_us, 3, tn);
    feats["cp_us"] = rnd(cp_us, 3, tn);
    feats["cp_fraction"] = rnd(cp_us / tu_us, 4, tn);
    feats["scs_khz"] = py_round(scs, 3);
    feats["bw_mhz"] = py_round(lk.bw / 1e6, 2);
    feats["carrier_mhz"] = py_round((fc + lk.cen) / 1e6, 3);
    feats["n_bursts"] = lk.runs.size();
    feats["duty_cycle"] = py_round(lk.duty, 3);
    feats["frame_period_ms"] = lk.frame_s ? Json(py_round(*lk.frame_s * 1e3, 3)) : Json(nullptr);
    feats["frame_contrast"] = py_round(lk.frame_contrast, 1);
    feats["symbol_grids"] = lk.n_grids;
    feats["line_z"] = py_round(lk.line_z, 1);
    feats["strong_windows"] = lk.win_z;
    feats["acf_z"] = py_round(lk.acf_z, 1);
    feats["link_snr_db"] = py_round(lk.snr_db, 1);
    feats["cfo_frac_hz"] = py_round(lk.cfo_frac_hz, 1);
    feats["cfo_note"] = "CP phase: carrier offset modulo the subcarrier spacing, vs band centre";
    {
      const auto bc = bands::uav_link_band(fc + lk.cen);
      feats["band_context"] = bc ? Json(*bc) : Json(nullptr);
    }
    // transmit streams per subcarrier (1 TX vs 2-TX STBC/SFBC), from the link's slice at N / Tu samples/s
    feats["tx_antennas"] = streams::tx_streams(lk.slice, lk.slice_rate, lk.tu, lk.ts, lk.bw, lk.runs, lk.snr_db,
                                               static_cast<double>(x.size()) / fs);
    cdvec().swap(lk.slice);
    // research catalog: which system / model the measured PHY most resembles
    {
      const Json& tx = feats["tx_antennas"];
      std::optional<int> na;
      if (tx.contains("tx_antennas") && tx["tx_antennas"].is_number()) na = tx["tx_antennas"].get<int>();
      feats["system_candidates"] = catalog::candidates(tu_us, cp_us / tu_us, lk.bw / 1e6, (fc + lk.cen) / 1e6, na,
                                                       lk.frame_s ? std::optional<double>(*lk.frame_s * 1e3) : std::nullopt);
    }
    const double conf = role == "unknown" ? 0.45 : 0.6;  // inferred from the PHY alone: "-like" / unknown
    const std::string lab = std::fabs(scs - static_cast<double>(py_round_int(scs))) < 0.05 ? fmt("ofdm_link_%.0fkHz", scs)
                                                                                           : fmt("ofdm_link_%.2fkHz", scs);
    for (const auto& r : lk.runs) {
      Json f = feats;
      f["burst_us"] = py_round((r.t1 - r.t0) * 1e6, 1);
      f["n_symbols"] = r.n;
      Emission e;
      e.family = "digital_vtx";
      e.label = lab;
      e.start_us = py_round(r.t0 * 1e6, 1);
      e.end_us = py_round(r.t1 * 1e6, 1);
      e.f_lo_mhz = py_round(lo, 3);
      e.f_hi_mhz = py_round(hi, 3);
      e.center_mhz = py_round((lo + hi) / 2, 3);
      e.snr_db = py_round(lk.snr_db, 1);
      e.confidence = conf;
      e.features = std::move(f);
      out.push_back(std::move(e));
    }
  }
  // hopping narrowband CP-OFDM links (RC uplinks): one Emission per hop, label ofdm_hop_<SCS>kHz
  for (const auto& hk : blind::blind_hops(sl, fs, fc, usable, links, cellular)) {
    const double tu_us = hk.tu * 1e6, ts_us = hk.ts * 1e6;
    const double cp_us = ts_us - tu_us;
    const double scs = 1e3 / tu_us;
    const auto [name, matches] = classify_numerology(tu_us, hk.bw / 1e6);
    const std::string role = numerology_role(name);
    Json feats = Json::object();
    feats["numerology"] = name;
    feats["matches"] = matches;
    feats["role"] = role;
    feats["tu_us"] = py_round(tu_us, 3);
    feats["symbol_us"] = py_round(ts_us, 3);
    feats["cp_us"] = py_round(cp_us, 3);
    feats["cp_fraction"] = py_round(cp_us / tu_us, 4);
    feats["scs_khz"] = py_round(scs, 3);
    feats["hop_bw_mhz"] = py_round(hk.bw / 1e6, 2);
    feats["hop_us"] = py_round(hk.dur_s * 1e6, 1);
    feats["n_hops"] = hk.hops.size();
    feats["n_channels"] = hk.n_channels;
    feats["hop_period_ms"] = hk.frame_s ? Json(py_round(*hk.frame_s * 1e3, 3)) : Json(nullptr);
    Json pl = Json::array();
    for (double c : hk.paired) pl.push_back(py_round((fc + c) / 1e6, 3));
    feats["paired_link_mhz"] = pl;
    {
      const auto bc = hk.hops.empty() ? std::nullopt : bands::uav_link_band(fc + hk.hops[0].channel);
      feats["band_context"] = bc ? Json(*bc) : Json(nullptr);
    }
    const double conf = role == "drone" || role == "infrastructure" ? 0.6 : 0.45;
    const std::string lab = std::fabs(scs - static_cast<double>(py_round_int(scs))) < 0.05 ? fmt("ofdm_hop_%.0fkHz", scs)
                                                                                           : fmt("ofdm_hop_%.2fkHz", scs);
    for (const auto& h : hk.hops) {
      const double lo = (fc + h.lo) / 1e6, hi = (fc + h.hi) / 1e6;
      const double snr = 10 * std::log10(std::max(std::pow(10.0, h.contrast / 10) - 1, 1e-3));
      Json f = feats;
      f["channel_mhz"] = py_round((fc + h.channel) / 1e6, 2);
      f["bw_mhz"] = py_round(hi - lo, 2);
      f["burst_us"] = py_round((h.t1 - h.t0) * 1e6, 1);
      f["n_symbols"] = py_round_int((h.t1 - h.t0) / hk.ts);
      f["fold_z"] = py_round(h.z, 1);
      f["contrast_db"] = py_round(h.contrast, 1);
      f["tracked"] = h.tracked;
      Emission e;
      e.family = "digital_vtx";
      e.label = lab;
      e.start_us = py_round(h.t0 * 1e6, 1);
      e.end_us = py_round(h.t1 * 1e6, 1);
      e.f_lo_mhz = py_round(lo, 3);
      e.f_hi_mhz = py_round(hi, 3);
      e.center_mhz = py_round((lo + hi) / 2, 3);
      e.snr_db = np::round(snr, 1);
      e.confidence = conf;
      e.features = std::move(f);
      out.push_back(std::move(e));
    }
  }
  return out;
}

}  // namespace detail
}  // namespace si::rx::digital_vtx
