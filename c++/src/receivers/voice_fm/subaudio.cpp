// python: receivers/voice_fm/subaudio.py
#include "receivers/voice_fm/subaudio.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <set>
#include <stdexcept>

#include "core/fft.hpp"
#include "core/npcompat.hpp"
#include "core/signal.hpp"
#include "receivers/voice_fm/channels.hpp"

namespace si::rx::voice_fm {

namespace {
std::string oct3(uint32_t v) {
  char b[8];
  std::snprintf(b, sizeof b, "%03o", v);
  return b;
}
inline uint32_t rot(uint32_t w, int r) { return ((w >> r) | (w << (23 - r))) & kMask23; }
}  // namespace

uint32_t golay_syndrome(uint32_t w, uint32_t g) {
  for (int i = 22; i > 10; --i)
    if (w >> i & 1) w ^= g << (i - 11);
  return w;
}

uint32_t dcs_encode(const std::string& code_octal) {
  const uint32_t d = static_cast<uint32_t>(std::stoul(code_octal, nullptr, 8)) | 0x800;
  for (uint32_t p = 0; p < 2048; ++p)
    if (golay_syndrome(d | p << 12) == 0) return d | p << 12;
  throw std::invalid_argument(code_octal);
}

std::vector<std::pair<std::string, std::string>> dcs_aliases(uint32_t word) {
  std::set<std::pair<std::string, std::string>> out;
  for (int inv = 0; inv < 2; ++inv) {
    const uint32_t v = word ^ (inv ? kMask23 : 0);
    for (int r = 0; r < 23; ++r) {
      const uint32_t x = rot(v, r);
      if (((x >> 9) & 7) == 4 && golay_syndrome(x) == 0) out.insert({oct3(x & 0x1FF), inv ? "I" : "N"});
    }
  }
  return {out.begin(), out.end()};
}

dvec to_subaudio(const fvec& disc, double fs) {
  const dmat sos = sig::butter_sos(6, 300.0, sig::BType::Lowpass, fs);
  const float m = np::mean(disc);
  fvec x(disc.size());
  for (size_t i = 0; i < x.size(); ++i) x[i] = disc[i] - m;
  const dvec a = sig::sosfiltfilt(sos, x);
  const dvec t = np::arange<double>(0.0, static_cast<double>(a.size()) / fs, 1 / kSubFs);
  dvec xp(a.size());
  for (size_t i = 0; i < a.size(); ++i) xp[i] = static_cast<double>(i) / fs;
  return np::interp(t, xp, a);
}

std::optional<Json> dcs(const dvec& a, int min_words) {
  if (static_cast<double>(a.size()) / kSubFs * kDcsRate < 23 * (min_words + 1)) return std::nullopt;
  dvec t(a.size());
  for (size_t i = 0; i < a.size(); ++i) t[i] = static_cast<double>(i) / kSubFs;
  double best_eye = 0;
  dvec v;
  bool have = false;
  for (double ph : np::linspace<double>(0.0, 1 / kDcsRate, 16, false)) {
    const dvec vv = np::interp(np::arange<double>(ph + 0.5 / kDcsRate, t.back(), 1 / kDcsRate), t, a);
    dvec m(vv.size());
    for (size_t i = 0; i < vv.size(); ++i) m[i] = std::fabs(vv[i]);
    const double eye = np::mean(m) / (np::stddev(m) + 1e-9);
    if (!have || eye > best_eye) {
      best_eye = eye;
      v = vv;
      have = true;
    }
  }
  std::vector<uint8_t> bits(v.size());
  for (size_t i = 0; i < v.size(); ++i) bits[i] = v[i] > 0;
  using Key = std::pair<std::string, std::string>;
  std::vector<std::pair<Key, int>> votes;   // Counter, in insertion order
  std::vector<std::pair<Key, uint32_t>> words;
  auto find = [](auto& vec, const Key& k) {
    for (auto& e : vec)
      if (e.first == k) return &e;
    return static_cast<decltype(&vec[0])>(nullptr);
  };
  for (size_t s = 0; s + 22 < bits.size(); ++s) {
    uint32_t w = 0;
    for (int i = 0; i < 23; ++i) w |= static_cast<uint32_t>(bits[s + i]) << i;
    for (int inv = 0; inv < 2; ++inv) {
      const uint32_t x = w ^ (inv ? kMask23 : 0);
      if (((x >> 9) & 7) == 4 && golay_syndrome(x) == 0) {
        const Key k{oct3(x & 0x1FF), inv ? "I" : "N"};
        if (auto* e = find(votes, k)) ++e->second;
        else votes.push_back({k, 1});
        if (!find(words, k)) words.push_back({k, w});
      }
    }
  }
  if (votes.empty()) return std::nullopt;
  size_t ti = 0;
  for (size_t i = 1; i < votes.size(); ++i)
    if (votes[i].second > votes[ti].second) ti = i;
  const Key top = votes[ti].first;
  if (votes[ti].second < min_words) return std::nullopt;
  const auto family = dcs_aliases(find(words, top)->second);
  int n_family = 0;
  for (const auto& k : family)
    if (auto* e = find(votes, k)) n_family += e->second;
  std::vector<Key> stdk;
  for (const auto& k : family)
    if (is_standard_dcs(k.first)) stdk.push_back(k);
  const std::vector<Key>& pool = stdk.empty() ? family : stdk;
  Key primary = pool[0];
  for (const auto& k : pool) {
    const auto kk = std::make_pair(k.second != "N", k.first), kp = std::make_pair(primary.second != "N", primary.first);
    if (kk < kp) primary = k;
  }
  size_t agree = 0, nrep = bits.size() > 23 ? bits.size() - 23 : 0;
  for (size_t i = 0; i < nrep; ++i) agree += bits[i + 23] == bits[i];
  const double rep = nrep ? static_cast<double>(agree) / static_cast<double>(nrep) : std::nan("");
  dvec av(v.size());
  for (size_t i = 0; i < v.size(); ++i) av[i] = std::fabs(v[i]);
  Json aliases = Json::array();
  for (const auto& k : family)
    if (k != primary) aliases.push_back("D" + k.first + k.second);
  Json r = Json::object();
  r["type"] = "DCS";
  r["dcs_code"] = "D" + primary.first + primary.second;
  r["dcs_aliases"] = aliases;
  r["dcs_standard"] = is_standard_dcs(primary.first);
  r["dcs_words_verified"] = n_family;
  r["dcs_bits"] = static_cast<int64_t>(bits.size());
  r["dcs_bit_repeat_agreement"] = py_round(rep, 3);
  r["dcs_deviation_hz"] = py_round(np::median(av), 0);
  r["dcs_note"] = "Golay(23,12)-verified; N/I assumes a 1 bit = positive deviation";
  return r;
}

std::optional<Json> ctcss(const dvec& a, double min_s) {
  if (static_cast<double>(a.size()) / kSubFs < min_s) return std::nullopt;
  const size_t n = std::max(a.size(), static_cast<size_t>(kSubFs * 10));
  const dvec w = np::hanning(a.size());
  dvec aw(a.size());
  for (size_t i = 0; i < a.size(); ++i) aw[i] = a[i] * w[i];
  const cdvec X = fft::rfft(aw, n);
  dvec A(X.size());
  for (size_t i = 0; i < X.size(); ++i) A[i] = np::cabs2(X[i]);
  const dvec F = fft::rfftfreq(n, 1 / kSubFs);
  size_t k = 0;
  bool have = false;
  dvec band;
  for (size_t i = 0; i < F.size(); ++i)
    if (F[i] >= 60 && F[i] <= 260) {
      band.push_back(A[i]);
      if (!have || A[i] > A[k]) {
        k = i;
        have = true;
      }
    }
  if (!have) return std::nullopt;
  const double snr = 10 * std::log10(A[k] / (np::median(band) + 1e-30));
  if (snr < 15) return std::nullopt;
  const double y0 = std::log(A[k - 1] + 1e-30), y1 = std::log(A[k] + 1e-30), y2 = std::log(A[k + 1] + 1e-30);
  const double f = F[k] + 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2 + 1e-30) * (F[1] - F[0]);
  const auto& tones = ctcss_tones();
  double stdt = tones[0];
  for (double c : tones)
    if (std::fabs(c - f) < std::fabs(stdt - f)) stdt = c;
  const double res = kSubFs / static_cast<double>(a.size());
  const bool ok = std::fabs(stdt - f) <= std::max(0.6, 0.3 * res);
  Json r = Json::object();
  r["type"] = "CTCSS";
  r["ctcss_hz"] = ok ? Json(stdt) : Json(nullptr);
  r["ctcss_measured_hz"] = py_round(f, 2);
  r["ctcss_snr_db"] = py_round(snr, 1);
  r["ctcss_standard"] = ok;
  r["ctcss_deviation_hz"] = py_round(std::sqrt(2.0) * np::stddev(a), 0);
  return r;
}

std::optional<Json> subaudible(const fvec& disc, double fs) {
  const dvec a = to_subaudio(disc, fs);
  if (auto d = dcs(a)) return d;
  return ctcss(a);
}

}  // namespace si::rx::voice_fm
