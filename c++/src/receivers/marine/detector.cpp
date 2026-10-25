// python receivers/marine/detector.py (see there for the method and the sources).
#include "receivers/marine/detector.hpp"

#include <algorithm>
#include <cmath>
#include <map>

#include "core/npcompat.hpp"
#include "receivers/marine/ais.hpp"

namespace si::rx::marine {

namespace {
constexpr double kRate = 9600.0, kFs = 48e3, kHalfBw = 10e3;

std::string hex(const std::vector<uint8_t>& b) {
  static const char* d = "0123456789abcdef";
  std::string s;
  for (uint8_t v : b) s += d[v >> 4], s += d[v & 15];
  return s;
}
}  // namespace

std::vector<Json> channel_frames(const Channelizer& ch, double f_meas) {
  const auto sl = ch.slice(f_meas, kHalfBw, kFs);
  const cfvec& y = sl.y;
  const double fs = sl.fs;
  if (y.size() < 512) return {};
  dvec d(y.size() - 1);
  // numpy temporary elision (see voice_fm): a conj() temporary >= 256 KiB receives the product in place
  const bool elide = (y.size() - 1) * sizeof(cf) >= 256 * 1024;
  for (size_t i = 1; i < y.size(); ++i) {
    const cf z = elide ? np::cmul(std::conj(y[i - 1]), y[i]) : np::cmul(y[i], std::conj(y[i - 1]));
    d[i - 1] = static_cast<double>(std::atan2(z.imag(), z.real()));
  }
  const double sps = fs / kRate;
  const size_t w = static_cast<size_t>(py_round_int(32 * sps));
  if (d.size() > w) {
    const dvec ma = np::convolve(d, dvec(w, 1.0 / static_cast<double>(w)), np::ConvMode::Same);
    for (size_t i = 0; i < d.size(); ++i) d[i] = d[i] - ma[i];
  }
  const int64_t r = py_round_int(sps);
  d = np::convolve(d, dvec(static_cast<size_t>(r), 1.0 / static_cast<double>(r)), np::ConvMode::Same);
  const int64_t nb = static_cast<int64_t>((static_cast<double>(d.size()) - 1) / sps);
  std::map<std::string, double> seen;
  std::vector<Json> out;
  for (int64_t ph = 0; ph < r; ++ph) {
    std::vector<int64_t> idx;
    for (int64_t k = 0; k < nb; ++k) {
      const int64_t v = static_cast<int64_t>(static_cast<double>(ph) + static_cast<double>(k) * sps);
      if (v < static_cast<int64_t>(d.size())) idx.push_back(v);
    }
    std::vector<int> lv(idx.size());
    for (size_t k = 0; k < idx.size(); ++k) lv[k] = d[idx[k]] > 0 ? 1 : 0;
    for (const auto& [pos, octs] : ais::frames(ais::unnrzi(lv))) {
      const double t = static_cast<double>(idx[pos]) / fs;
      const std::string key = hex(octs);
      const auto it = seen.find(key);
      if (it != seen.end() && std::fabs(it->second - t) < 0.01) continue;
      seen[key] = t;
      const int64_t nbits = 8 * static_cast<int64_t>(octs.size());
      Json o = Json::object();
      o["t_s"] = py_round(t, 5);
      o["octets"] = key;
      o["bits"] = nbits;
      o["end_s"] = py_round(t + static_cast<double>(nbits + 40) / kRate, 5);
      out.push_back(o);
    }
  }
  std::stable_sort(out.begin(), out.end(),
                   [](const Json& a, const Json& b) { return a["t_s"].get<double>() < b["t_s"].get<double>(); });
  return out;
}

std::vector<Emission> detect(const Channelizer& ch, double ppm) {
  const double lo = ch.fc() - 0.45 * ch.fs(), hi = ch.fc() + 0.45 * ch.fs();
  std::vector<Emission> out;
  for (const auto& [name, f] : std::vector<std::pair<std::string, double>>{{"AIS 1 (87B)", 161.975e6},
                                                                           {"AIS 2 (88B)", 162.025e6}}) {
    const double f_meas = f / (1 + ppm * 1e-6);
    if (!(lo + 15e3 < f_meas && f_meas < hi - 15e3)) continue;
    for (const Json& r : channel_frames(ch, f_meas)) {
      const std::string h = r["octets"].get<std::string>();
      std::vector<uint8_t> octs;
      for (size_t i = 0; i + 1 < h.size(); i += 2) octs.push_back(static_cast<uint8_t>(std::stoul(h.substr(i, 2), nullptr, 16)));
      Json dec = Json::object();
      dec["crc_ok"] = true;
      const Json msg = ais::fields(octs);
      for (auto it = msg.begin(); it != msg.end(); ++it) dec[it.key()] = it.value();
      dec["raw_hex"] = h;
      Emission e;
      e.family = "marine";
      e.label = "ais";
      e.start_us = r["t_s"].get<double>() * 1e6;
      e.end_us = r["end_s"].get<double>() * 1e6;
      e.f_lo_mhz = (f - 12.5e3) / 1e6;
      e.f_hi_mhz = (f + 12.5e3) / 1e6;
      e.center_mhz = f / 1e6;
      e.confidence = 0.97;
      e.decoded = dec;
      e.features = Json{{"channel", name}, {"modulation", "GMSK 9600 bit/s"}};
      out.push_back(std::move(e));
    }
  }
  return out;
}

}  // namespace si::rx::marine
