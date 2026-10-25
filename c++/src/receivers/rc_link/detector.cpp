#include "receivers/rc_link/detector.hpp"
#include "receivers/rc_link/elrs.hpp"
#include "receivers/rc_link/tables.hpp"

#include <algorithm>
#include <cmath>

#include "core/fft.hpp"
#include "core/npcompat.hpp"
#include "core/parallel.hpp"
#include "receivers/rc_link/dsp_util.hpp"
#include "receivers/rc_link/fsk.hpp"
#include "receivers/rc_link/lora.hpp"

namespace si::rx::rc_link {
namespace detail {

namespace {
// python (e.snr_db or -99): None and 0.0 are falsy
double snr_or(const Emission& e) { return (e.snr_db && *e.snr_db != 0.0) ? *e.snr_db : -99.0; }
}  // namespace

std::vector<Emission> suppress_images(const std::vector<Emission>& ems, double fc, double min_db) {
  std::vector<Emission> keep;
  for (size_t i = 0; i < ems.size(); ++i) {
    const Emission& e = ems[i];
    const double fm = 2 * fc / 1e6 - e.center_mhz;
    const double bw = std::max(e.f_hi_mhz - e.f_lo_mhz, 0.2);
    bool img = false;
    for (size_t j = 0; j < ems.size() && !img; ++j) {
      if (!(std::fabs(ems[j].center_mhz - fm) <= bw / 2 + 0.1) || j == i) continue;
      const Emission& o = ems[j];
      const double ov = std::min(o.end_us, e.end_us) - std::max(o.start_us, e.start_us);
      if (ov > 0.5 * (e.end_us - e.start_us) && snr_or(o) > snr_or(e) + min_db &&
          (e.end_us - e.start_us) <= 1.3 * (o.end_us - o.start_us))
        img = true;
    }
    if (!img) keep.push_back(e);
  }
  return keep;
}

bool is_raw_image(const cfvec& x, double fs, double fc, const Emission& e, double min_db) {
  const double bw = std::max((e.f_hi_mhz - e.f_lo_mhz) * 1e6, 200e3);
  const double ofs = std::min({std::max(2 * bw, 1e6), 8e6, fs});
  const auto [X, off] = local_fft(x, fs, e.start_us, e.end_us, 50.0);
  if (X.empty()) return false;
  const double a = (e.start_us - off) * 1e-6, b = (e.end_us - off) * 1e-6;
  float p[2];
  int k = 0;
  for (double f : {e.center_mhz * 1e6, 2 * fc - e.center_mhz * 1e6}) {
    const SliceOut s = slice(X, fs, FOff{f - fc, false}, bw / 2, ofs);
    const int64_t n = static_cast<int64_t>(s.y.size());
    const int64_t i0 = static_cast<int64_t>(a * s.ofs);
    const int64_t i1 = std::max(static_cast<int64_t>(b * s.ofs), i0 + 1);
    const int64_t lo = std::clamp<int64_t>(i0, 0, n), hi = std::clamp<int64_t>(i1, 0, n);
    fvec v;
    for (int64_t i = lo; i < hi; ++i) v.push_back(np::cabs2(s.y[i]));
    p[k++] = v.empty() ? std::nanf("") : mean_f(v);
  }
  const float den = (1e-30 > static_cast<double>(p[0])) ? static_cast<float>(1e-30) : p[0];
  return 10.0f * std::log10(p[1] / den) > static_cast<float>(min_db);
}

// python _elrs_hops(ems): ExpressLRS hop recovery (research, PLAN R2d) from a CRC-verified SYNC packet and the link's
// other packets in the block; every domain whose grid holds the SYNC packet is tried with the packets on its grid
void elrs_hops(std::vector<Emission>& ems) {
  for (auto& s : ems) {
    Json& d = s.decoded;
    if (!d.is_object() || d.value("protocol", Json(nullptr)) != "ExpressLRS" || d.value("packet_type", Json(nullptr)) != "SYNC" ||
        !(d.contains("crc_ok") && d["crc_ok"].is_boolean() && d["crc_ok"].get<bool>()))
      continue;
    const Json modes = s.features.value("elrs_modes", Json(nullptr));
    if (!modes.is_array() || modes.empty()) continue;
    const ElrsMode* mode = nullptr;
    for (const auto& m : elrs_lora())
      if (modes[0] == m.name) mode = &m;
    if (!mode) continue;
    std::vector<const Emission*> link;
    for (const auto& e : ems)
      if (&e != &s && e.features.value("modulation", Json(nullptr)) == "lora" &&
          e.features.value("sf", Json(nullptr)) == s.features.value("sf", Json(nullptr)) &&
          e.features.value("bw_hz", Json(nullptr)) == s.features.value("bw_hz", Json(nullptr)))
        link.push_back(&e);
    const Json fsj = s.features.value("f_carrier_hz", Json(nullptr));
    if (!fsj.is_number()) continue;
    std::vector<Json> fits;
    for (const auto& dom : elrs_domains()) {
      const double sp = (dom.f1 - dom.f0) / (dom.n - 1);
      auto idx = [&](double f) { return py_round_int((f - dom.f0) / sp); };
      auto ok = [&](double f) {
        const int64_t i = idx(f);
        return 0 <= i && i < dom.n && std::fabs(f - (dom.f0 + static_cast<double>(i) * sp)) < 150e3;
      };
      if (!ok(fsj.get<double>())) continue;
      std::vector<std::pair<double, int>> pk;
      for (const Emission* e : link) {
        const double f = e->features["f_carrier_hz"].get<double>();
        if (ok(f)) pk.push_back({e->start_us, static_cast<int>(idx(f))});
      }
      Json r = elrs::recover_hops(d, s.start_us, static_cast<int>(idx(fsj.get<double>())), mode->interval_us,
                                  mode->hop_every, dom.n, pk);
      if (r.is_object() && r["uid23_candidates"].get<int64_t>() > 0) {
        r["mode"] = mode->name;
        r["domain"] = dom.name;
        fits.push_back(r);
      }
    }
    if (!fits.empty()) {
      Json h = fits[0];
      Json others = Json::array();
      for (size_t i = 1; i < fits.size(); ++i) others.push_back(fits[i]["domain"]);
      h["other_domains"] = others;
      d["hop_recovery"] = h;
    }
  }
}

}  // namespace detail

std::vector<Emission> detect(const cfvec& x, double fs, double fc, const cfvec* X,
                             const std::pair<double, double>* lora_only) {
  using namespace detail;
  cfvec own;
  if (!X || X->size() != x.size()) {
    own = fft::fft(x);
    X = &own;
  }
  std::vector<Emission> ems = detect_lora(*X, fs, fc, lora_only);
  std::vector<Emission> fsk = lora_only ? std::vector<Emission>{} : detect_fsk(x, fs, fc);
  // an FSK candidate covering a LoRa packet is the same burst seen by energy
  std::vector<Emission> all = ems;
  for (auto& e : fsk) {
    bool cov = false;
    for (const auto& l : ems)
      if (std::fabs(e.center_mhz - l.center_mhz) < 0.6 * (l.f_hi_mhz - l.f_lo_mhz) &&
          std::min(e.end_us, l.end_us) - std::max(e.start_us, l.start_us) > 0) {
        cov = true;
        break;
      }
    if (!cov) all.push_back(std::move(e));
  }
  all = suppress_images(all, fc);
  const auto img = pmap(all, [&](const Emission& e) { return is_raw_image(x, fs, fc, e) ? 1 : 0; });
  std::vector<Emission> out;
  for (size_t i = 0; i < all.size(); ++i)
    if (!img[i]) out.push_back(std::move(all[i]));
  std::stable_sort(out.begin(), out.end(), [](const Emission& a, const Emission& b) { return a.start_us < b.start_us; });
  elrs_hops(out);
  return out;
}

}  // namespace si::rx::rc_link
