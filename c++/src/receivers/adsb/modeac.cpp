// python receivers/adsb/modeac.py (see there for the format and the sources).
#include "receivers/adsb/modeac.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <set>
#include <tuple>

#include "core/npcompat.hpp"

namespace si::rx::adsb::modeac {

namespace {

const char* const kPositions[13] = {"C1", "A1", "C2", "A2", "C4", "A4", "X", "B1", "D1", "B2", "D2", "B4", "D4"};
constexpr double kTolUs = 0.1;
constexpr double kAmpDb = 3.0;
constexpr double kWeakDb = 9.0;
constexpr double kPulseDb = 12.0;

int gray_to_bin(int g) {
  int b = 0;
  while (g) {
    b ^= g;
    g >>= 1;
  }
  return b;
}

struct Pulse {
  double t, w, a;
  bool operator<(const Pulse& o) const { return std::tie(t, w, a) < std::tie(o.t, o.w, o.a); }
};

// python _pulses(env, fs, thr)
std::vector<Pulse> pulses(const fvec& env, double fs, double thr) {
  const float th = static_cast<float>(thr);
  std::vector<Pulse> out;
  const size_t n = env.size();
  size_t a = 0;
  while (a < n) {
    if (!(env[a] > th)) {
      ++a;
      continue;
    }
    size_t b = a;
    while (b < n && env[b] > th) ++b;
    size_t c = a;
    for (size_t k = a; k < b; ++k)
      if (env[k] > env[c]) c = k;
    const double pk = env[c], h = 0.5 * pk;
    size_t i = c;
    while (i > 0 && static_cast<double>(env[i - 1]) > h) --i;
    size_t j = c + 1;
    while (j < n && static_cast<double>(env[j]) > h) ++j;
    double e0;
    if (i > 0) {
      const double lo = env[i - 1], hi = env[i];
      e0 = static_cast<double>(i) - 1 + (h - lo) / std::max(hi - lo, 1e-12);
    } else {
      e0 = static_cast<double>(i);
    }
    out.push_back({e0 / fs, static_cast<double>(j - i) / fs, pk});
    a = b;
  }
  return out;
}

Json opt_int(const std::optional<int>& v) { return v ? Json(*v) : Json(nullptr); }

}  // namespace

std::optional<int> gillham_ft(const Bits& p) {
  if (p.at("D1") || p.at("X")) return std::nullopt;
  int g500 = 0;
  for (const char* k : {"D2", "D4", "A1", "A2", "A4", "B1", "B2", "B4"}) g500 = (g500 << 1) | p.at(k);
  const int g100 = (p.at("C1") << 2) | (p.at("C2") << 1) | p.at("C4");
  const int n500 = gray_to_bin(g500);
  int n100 = gray_to_bin(g100);
  if (n100 == 0 || n100 == 5 || n100 == 6) return std::nullopt;
  if (n100 == 7) n100 = 5;
  if (n500 % 2) n100 = 6 - n100;
  return n500 * 500 + n100 * 100 - 1300;
}

std::string squawk(const Bits& p) {
  std::string s;
  for (char d : std::string("ABCD")) {
    const std::string k(1, d);
    s += std::to_string(p.at(k + "4") * 4 + p.at(k + "2") * 2 + p.at(k + "1"));
  }
  return s;
}

std::vector<Json> decode(const fvec& env, double fs, const std::vector<std::pair<double, double>>& busy, double noise) {
  if (noise < 0) noise = static_cast<double>(np::median(env));
  const double thr = std::max(noise, 1e-12) * std::pow(10.0, kPulseDb / 20);
  std::vector<Pulse> P;
  for (const Pulse& p : pulses(env, fs, thr)) {
    if (!(0.25e-6 <= p.w && p.w <= 0.8e-6)) continue;
    bool b = false;
    for (const auto& [a0, a1] : busy) b = b || (a0 - 1e-6 <= p.t && p.t <= a1);
    if (!b) P.push_back(p);
  }
  std::sort(P.begin(), P.end());
  std::vector<Json> out;
  if (P.size() < 2) return out;
  std::vector<double> t(P.size());
  for (size_t i = 0; i < P.size(); ++i) t[i] = P[i].t;
  const double tol = kTolUs * 1e-6 + 1 / fs;
  auto search = [&](double v) { return static_cast<size_t>(std::lower_bound(t.begin(), t.end(), v) - t.begin()); };
  std::set<size_t> used;
  for (size_t i = 0; i < P.size(); ++i) {
    if (used.count(i)) continue;
    const double t1 = P[i].t, a1 = P[i].a;
    const size_t j = search(t1 + kFrameUs * 1e-6 - tol);
    if (j >= P.size() || std::abs(t[j] - t1 - kFrameUs * 1e-6) > tol) continue;
    const double a2 = P[j].a;
    if (std::abs(20 * std::log10(a1 / a2)) > kAmpDb) continue;
    Bits bits;
    for (const char* k : kPositions) bits[k] = 0;
    bool ok = true;
    const double weak = std::min(a1, a2) * std::pow(10.0, -kWeakDb / 20);
    for (size_t k = i + 1; k < j; ++k) {
      if (P[k].a < weak) continue;
      const double q = (t[k] - t1) / (kStepUs * 1e-6);
      const int64_t n = py_round_int(q);
      if (!(1 <= n && n <= 13) || std::abs(t[k] - t1 - static_cast<double>(n) * kStepUs * 1e-6) > tol ||
          bits[kPositions[n - 1]]) {
        ok = false;
        break;
      }
      bits[kPositions[n - 1]] = 1;
    }
    if (!ok || bits["X"]) continue;
    const size_t ks = search(t1 + kSpiUs * 1e-6 - tol);
    const bool spi = ks < P.size() && std::abs(t[ks] - t1 - kSpiUs * 1e-6) <= tol;
    std::string word;
    for (const char* k : kPositions) word += bits[k] ? '1' : '0';
    Json r = Json::object();
    r["t_s"] = py_round(t1, 7);
    r["code"] = squawk(bits);
    r["altitude_ft_if_mode_c"] = opt_int(gillham_ft(bits));
    r["spi"] = spi;
    r["pulses"] = word;
    r["snr_db"] = py_round(20 * std::log10(0.5 * (a1 + a2) / std::max(noise, 1e-12)), 1);
    out.push_back(std::move(r));
    for (size_t k = i; k <= j; ++k) used.insert(k);
  }
  return out;
}

Json summarise(const std::vector<Json>& replies) {
  static const std::map<std::string, std::string> special = {{"7500", "unlawful interference (hijack)"},
                                                             {"7600", "radio failure"},
                                                             {"7700", "emergency"}};
  std::vector<std::string> order;
  std::map<std::string, Json> by;
  for (const Json& r : replies) {
    const std::string w = r["pulses"].get<std::string>();
    if (!by.count(w)) {
      Json g = Json::object();
      g["code"] = r["code"];
      g["altitude_ft_if_mode_c"] = r["altitude_ft_if_mode_c"];
      g["n"] = 0;
      g["spi"] = false;
      g["t_first_s"] = r["t_s"];
      by[w] = g;
      order.push_back(w);
    }
    Json& g = by[w];
    g["n"] = g["n"].get<int64_t>() + 1;
    g["spi"] = g["spi"].get<bool>() || r["spi"].get<bool>();
  }
  std::vector<Json> words;
  for (const auto& w : order) words.push_back(by[w]);
  std::stable_sort(words.begin(), words.end(),
                   [](const Json& a, const Json& b) { return a["n"].get<int64_t>() > b["n"].get<int64_t>(); });
  int64_t rep = 0;
  for (Json& g : words) {
    g["repeated"] = g["n"].get<int64_t>() >= 2;
    rep += g["repeated"].get<bool>();
    const auto it = special.find(g["code"].get<std::string>());
    if (it != special.end()) g["special_code"] = it->second;
  }
  Json o = Json::object();
  o["n_replies"] = static_cast<int64_t>(replies.size());
  o["n_distinct"] = static_cast<int64_t>(words.size());
  o["n_repeated"] = rep;
  o["replies"] = Json::array();
  for (size_t i = 0; i < words.size() && i < 50; ++i) o["replies"].push_back(words[i]);
  return o;
}

}  // namespace si::rx::adsb::modeac
