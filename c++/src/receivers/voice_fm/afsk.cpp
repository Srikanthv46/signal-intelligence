// python receivers/voice_fm/afsk.py (see there for the method and the sources).
#include "receivers/voice_fm/afsk.hpp"

#include <algorithm>
#include <cmath>
#include <set>

#include "core/hdlc.hpp"
#include "core/npcompat.hpp"

namespace si::rx::voice_fm::afsk {

namespace {

constexpr double kBaud = 1200.0, kMark = 1200.0, kSpace = 2200.0;
constexpr int kPhases = 8;
constexpr double kMinS = 0.12;

// python _tone_energy(d, fs, f, w): |sum of the last w samples of d e^-j2pi f n / fs|^2
dvec tone_energy(const fvec& d, double fs, double f, size_t w) {
  const size_t n = d.size();
  std::vector<cd> c(n + 1, cd(0, 0));
  const double k = -2 * kPi * f;
  cd acc(0, 0);
  for (size_t i = 0; i < n; ++i) {
    const double th = k * (static_cast<double>(i) / fs);
    const double v = d[i];
    acc += cd(v * std::cos(th), v * std::sin(th));
    c[i + 1] = acc;
  }
  dvec out(n + 1 - w);
  for (size_t i = 0; i + w <= n; ++i) {
    const double a = np::cabs(c[i + w] - c[i]);
    out[i] = a * a;
  }
  return out;
}

std::string hex(const std::vector<uint8_t>& b) {
  static const char* d = "0123456789abcdef";
  std::string s;
  for (uint8_t v : b) s += d[v >> 4], s += d[v & 15];
  return s;
}

}  // namespace

std::vector<Json> decode(const fvec& d, double fs) {
  if (static_cast<double>(d.size()) / fs < kMinS) return {};
  const double spb = fs / kBaud;
  const size_t w = static_cast<size_t>(py_round_int(spb));
  const dvec m = tone_energy(d, fs, kMark, w), s = tone_energy(d, fs, kSpace, w);
  dvec soft(m.size());
  for (size_t i = 0; i < m.size(); ++i) soft[i] = m[i] - s[i];
  const int64_t nb = static_cast<int64_t>(static_cast<double>(soft.size()) / spb) - 1;
  std::vector<Json> out;
  std::set<std::vector<uint8_t>> seen;
  for (int ph = 0; ph < kPhases; ++ph) {
    std::vector<int64_t> idx;
    for (int64_t k = 0; k < nb; ++k) {
      const int64_t v = static_cast<int64_t>(static_cast<double>(ph) * spb / kPhases + static_cast<double>(k) * spb);
      if (v < static_cast<int64_t>(soft.size())) idx.push_back(v);
    }
    std::vector<int> lv(idx.size());
    for (size_t k = 0; k < idx.size(); ++k) lv[k] = soft[idx[k]] > 0 ? 1 : 0;
    for (const auto& [pos, octs] : hdlc::deframe(hdlc::unnrzi(lv))) {
      if (seen.count(octs)) continue;
      const auto f = hdlc::ax25(octs);
      if (!f) continue;
      seen.insert(octs);
      Json r = Json::object();
      r["t_s"] = py_round(static_cast<double>(idx[pos] + static_cast<int64_t>(w)) / fs, 4);
      r["raw_hex"] = hex(octs);
      r["src"] = f->src;
      r["dest"] = f->dest;
      r["path"] = Json(f->path);
      r["info"] = hdlc::latin1_utf8(f->info);
      const Json a = hdlc::aprs(f->dest, f->info);
      if (!a.empty()) r["aprs"] = a;
      out.push_back(r);
    }
  }
  std::stable_sort(out.begin(), out.end(),
                   [](const Json& a, const Json& b) { return a["t_s"].get<double>() < b["t_s"].get<double>(); });
  return out;
}

}  // namespace si::rx::voice_fm::afsk
