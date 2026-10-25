// python: receivers/emitters/bursts.py extract_bursts
#include "receivers/emitters/bursts.hpp"

#include <algorithm>
#include <map>
#include <string>
#include <utility>

#include "core/ofdm15k.hpp"
#include "core/spectral.hpp"

namespace si::rx::emitters {

std::vector<Emission> extract_bursts(const cfvec& x, double fs, double fc, double usable_half) {
  static const std::map<std::string, std::pair<std::string, std::string>> LABELS = {
      {"ocusync_uplink_like", {"dji_ocusync", "ocusync_uplink"}},
      {"droneid_candidate", {"dji_droneid", "droneid"}},
      {"ocusync_video_like", {"dji_ocusync", "ocusync_video"}},
      {"ofdm15k_other", {"unknown_ofdm", "ofdm15k_other"}}};
  if (usable_half == 0) usable_half = std::min(fs, 56e6) / 2;
  std::vector<std::pair<double, double>> notch;
  for (const auto& s : find_spurs(x, fs, fc)) notch.emplace_back(s.offset_hz, 30e3);
  cdvec xd(x.size());
  for (size_t i = 0; i < x.size(); ++i) xd[i] = cd(x[i].real(), x[i].imag());
  const SpectralChannelizer ch(xd, fs, notch);
  const ofdm15k::Slicer slicer = [&](double off, double half, double ofs) { return ch.slice(off, half, ofs); };
  const auto scan = ofdm15k::ofdm15k_scan(slicer, ofdm15k::profiles(), 0.5e-3, 7.0, usable_half);
  std::vector<ofdm15k::Burst> refined;
  for (const auto& h : scan.hits)
    if (auto b = ofdm15k::refine_ofdm15k(slicer, h)) refined.push_back(*b);
  auto ob = ofdm15k::merge_ofdm15k(refined);
  auto tracked = ofdm15k::track_hops(slicer, ob, scan.stats, static_cast<double>(x.size()) / fs, 4.0, usable_half);
  ob.insert(ob.end(), tracked.begin(), tracked.end());
  ob = ofdm15k::merge_ofdm15k(ob);
  std::stable_sort(ob.begin(), ob.end(), [](const auto& a, const auto& b) { return a.start_s < b.start_s; });
  std::vector<Emission> out;
  for (const auto& b : ob) {
    const auto& [fam, label] = LABELS.at(ofdm15k::label_ofdm15k(b).first);
    const double cen = (fc + b.f_off) / 1e6;
    Emission e;
    e.family = fam;
    e.label = label;
    e.start_us = b.start_s * 1e6;
    e.end_us = b.end_s * 1e6;
    e.f_lo_mhz = cen - b.bw_hz / 2e6;
    e.f_hi_mhz = cen + b.bw_hz / 2e6;
    e.center_mhz = cen;
    e.confidence = std::min(0.95, 0.5 + 0.5 * b.score);
    e.features = Json{{"scs_khz", 15.0},          {"bw_mhz", b.bw_hz / 1e6},
                      {"symbol_us", b.period_us}, {"burst_us", (b.end_s - b.start_s) * 1e6},
                      {"nsym", b.nsym},           {"profile", b.profile},
                      {"cp_score", b.score},      {"tracked", b.tracked},
                      {"source", "assoc_bursts"}};
    out.push_back(std::move(e));
  }
  return out;
}

}  // namespace si::rx::emitters
