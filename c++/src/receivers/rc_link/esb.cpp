#include "receivers/rc_link/esb.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <string>

#include "core/ndimage.hpp"
#include "core/npcompat.hpp"
#include "core/pyfmt.hpp"
#include "receivers/rc_link/dsp_util.hpp"

namespace si::rx::rc_link::esb {
namespace {

constexpr double kPi = 3.14159265358979323846;

// python _crc0_prefix(bits): out[k] = CRC-16/CCITT register (init 0) after bits[:k]
std::vector<uint32_t> crc0_prefix(const std::vector<uint8_t>& bits) {
  std::vector<uint32_t> out(bits.size() + 1, 0);
  uint32_t c = 0;
  for (size_t k = 0; k < bits.size(); ++k) {
    const uint32_t fb = ((c >> 15) & 1u) ^ bits[k];
    c = ((c << 1) & 0xFFFFu) ^ (fb ? kCrcPoly : 0u);
    out[k + 1] = c;
  }
  return out;
}

// python _field16(bits): out[k] = the 16 bits starting at k, MSB first
std::vector<uint32_t> field16(const std::vector<uint8_t>& bits) {
  const size_t n = bits.size() >= 16 ? bits.size() - 15 : 0;
  std::vector<uint32_t> out(n, 0);
  for (size_t k = 0; k < n; ++k)
    for (size_t j = 0; j < 16; ++j) out[k] |= static_cast<uint32_t>(bits[k + j]) << (15 - j);
  return out;
}

// python _unstep(c, n): n inverse zero-input register steps
uint32_t unstep(uint32_t c, int n) {
  for (int i = 0; i < n; ++i) c = (c & 1u) ? (((c ^ kCrcPoly) >> 1) | 0x8000u) : c >> 1;
  return c;
}

// python _hex(bits): MSB-first bytes (zero-padded), upper-case hex
std::string hexbits(const uint8_t* b, size_t n) {
  std::string s;
  for (size_t i = 0; i < n; i += 8) {
    int v = 0;
    for (size_t j = 0; j < 8; ++j) v = (v << 1) | (i + j < n ? b[i + j] : 0);
    s += fmt("%02X", v);
  }
  return s;
}

}  // namespace

std::optional<Packet> slice_packet(const cf* y, size_t len, double ofs) {
  using namespace si::rx::rc_link::detail;
  const double sps = ofs / kRate;
  if (static_cast<double>(len) < 48 * sps || len < 2) return std::nullopt;
  const size_t n = len - 1;
  fvec ifr(n), pw(n);
  const bool el = elides(n, sizeof(cf));
  for (size_t i = 0; i < n; ++i) {
    const cf p = mul_conj(y[i + 1], y[i], el);
    ifr[i] = std::atan2(p.imag(), p.real());
    pw[i] = np::cabs2(y[i + 1]);
  }
  const float lim = 0.25f * np::percentile(pw, 90);
  std::vector<int64_t> on;
  for (size_t i = 0; i < n; ++i)
    if (pw[i] > lim) on.push_back(static_cast<int64_t>(i));
  if (static_cast<double>(on.size()) < 48 * sps) return std::nullopt;
  const int64_t a = on.front();
  const int64_t b = std::min<int64_t>(static_cast<int64_t>(n) - 1,
                                      static_cast<int64_t>(static_cast<double>(on.back()) + 4 * sps));
  const fvec ifs = ndi::uniform_filter1d(ifr, static_cast<size_t>(std::max<int64_t>(1, static_cast<int64_t>(sps / 2))));
  const auto [lo, hi] = kmeans2(fvec(ifs.begin() + a, ifs.begin() + b + 1));   // python _two_means == _kmeans2
  const float thr = (lo + hi) / 2.0f;
  float best_sc = -1.0f;
  fvec best_v, v;
  const int64_t nph = static_cast<int64_t>(std::ceil(sps));
  for (int64_t ph = 0; ph < nph; ++ph) {   // bit timing: the phase with the widest eye
    const int64_t nb = static_cast<int64_t>(static_cast<double>(b - a - ph) / sps);
    v.assign(static_cast<size_t>(std::max<int64_t>(nb, 0)), 0.0f);
    fvec av(v.size());
    for (int64_t i = 0; i < nb; ++i) {
      const int64_t idx = static_cast<int64_t>(static_cast<double>(a + ph) + static_cast<double>(i) * sps);
      v[i] = ifs[idx] - thr;
      av[i] = std::fabs(v[i]);
    }
    const float sc = av.empty() ? std::nanf("") : mean_f(av);
    if (best_sc < 0 || sc > best_sc) {
      best_sc = sc;
      best_v = v;
    }
  }
  std::vector<uint8_t> bits(best_v.size());
  for (size_t i = 0; i < bits.size(); ++i) bits[i] = best_v[i] > 0;
  const int64_t nalt = static_cast<int64_t>(bits.size()) - 1;   // alt[j]: bits j and j + 1 differ
  auto alt = [&](int64_t j) { return bits[j] != bits[j + 1]; };
  int64_t start = -1;
  for (int64_t s0 = 0; s0 < std::min<int64_t>(kPreambleStart, nalt); ++s0) {
    if (s0 > 0 && alt(s0 - 1)) continue;   // inside a run, not its first bit
    int64_t k = s0;
    while (k < nalt && alt(k)) ++k;
    if (k - s0 + 1 >= kMinPreamble) {      // bits s0..k alternate
      start = k;
      break;
    }
  }
  if (start < 0 || static_cast<int64_t>(bits.size()) - start < kMinSpan + 16) return std::nullopt;
  Packet p;
  p.bits.assign(bits.begin() + start, bits.end());
  p.preamble_bits = start;
  // float32 arithmetic as numpy (NEP 50: the Python floats are weak)
  p.dev_hz = static_cast<double>((hi - lo) / 2.0f * static_cast<float>(ofs) / static_cast<float>(2 * kPi));
  p.f_off_hz = static_cast<double>(thr * static_cast<float>(ofs) / static_cast<float>(2 * kPi));
  return p;
}

std::optional<double> frame_period_us(std::vector<double> t) {
  std::sort(t.begin(), t.end());
  if (t.size() < 3) return std::nullopt;
  dvec d;
  for (size_t i = 1; i < t.size(); ++i)
    if (t[i] - t[i - 1] > 20.0) d.push_back(t[i] - t[i - 1]);
  if (d.empty()) return std::nullopt;
  const double T = np::percentile(d, 10);
  dvec q(d.size());
  for (size_t i = 0; i < d.size(); ++i) q[i] = d[i] / std::max(std::nearbyint(d[i] / T), 1.0);
  return np::median(q);
}

std::vector<Json> verify_links(const std::vector<std::vector<uint8_t>>& frames, const std::vector<double>& starts_us) {
  const size_t nf = frames.size();
  std::vector<Json> out(nf);
  std::map<std::string, std::vector<double>> nodes;
  for (size_t i = 0; i < nf; ++i) {
    out[i] = Json::object();
    out[i]["crc_ok"] = false;
    out[i]["sync_hex"] = hexbits(frames[i].data(), std::min<size_t>(32, frames[i].size()));
    nodes[out[i]["sync_hex"].get<std::string>()].push_back(starts_us[i]);
  }
  for (size_t i = 0; i < nf; ++i)
    if (auto T = frame_period_us(nodes[out[i]["sync_hex"].get<std::string>()])) out[i]["node_frame_period_us"] = *T;
  std::vector<std::string> order;
  std::map<std::string, std::vector<size_t>> groups;
  for (size_t i = 0; i < nf; ++i) {
    const std::string key(frames[i].begin(), frames[i].begin() + std::min<size_t>(kSyncBits, frames[i].size()));
    if (!groups.count(key)) order.push_back(key);
    groups[key].push_back(i);
  }
  for (const auto& gk : order) {
    const auto& idx = groups[gk];
    if (idx.size() < static_cast<size_t>(kMinVerified)) continue;
    std::map<size_t, std::vector<uint32_t>> pre, fld;
    size_t n_max = 0;
    for (size_t i : idx) {
      pre[i] = crc0_prefix(frames[i]);
      fld[i] = field16(frames[i]);
      n_max = std::max(n_max, frames[i].size());
    }
    bool found = false;
    size_t best_n = 0;
    int bestE = 0;
    uint32_t bestK = 0;
    for (int E = kMinSpan; E < static_cast<int>(n_max) - 15; ++E) {
      std::map<uint32_t, std::set<std::string>> seen;
      for (size_t i : idx) {
        const auto& f = frames[i];
        if (f.size() < static_cast<size_t>(E) + 16) continue;
        const uint32_t K = pre[i][E] ^ fld[i][E];
        seen[K].insert(std::string(f.begin(), f.begin() + E + 16));
      }
      if (seen.empty()) continue;
      // python max(seen.items(), key=(len(set), -K)): most distinct packets, ties to the smallest K
      uint32_t K = 0;
      size_t ns = 0;
      for (const auto& [k, s] : seen)
        if (s.size() > ns) {
          ns = s.size();
          K = k;
        }
      auto twin = seen.find(K ^ kShortByOne);
      if (twin != seen.end() && 4 * twin->second.size() >= ns) continue;   // one bit short of the true span
      if (ns >= static_cast<size_t>(kMinVerified) && 2 * ns >= idx.size() && (!found || ns > best_n)) {
        found = true;
        best_n = ns;
        bestE = E;
        bestK = K;
      }
    }
    if (!found) continue;
    const uint32_t init = unstep(bestK, bestE);
    int n_ok = 0;
    for (size_t i : idx)
      if (frames[i].size() >= static_cast<size_t>(bestE) + 16 && (pre[i][bestE] ^ fld[i][bestE]) == bestK) ++n_ok;
    for (size_t i : idx) {
      const auto& f = frames[i];
      const bool ok = f.size() >= static_cast<size_t>(bestE) + 16 && (pre[i][bestE] ^ fld[i][bestE]) == bestK;
      Json& o = out[i];
      o["crc_ok"] = ok;
      o["crc"] = "CRC-16/CCITT (0x1021)";
      o["crc_span_bits"] = bestE;
      o["crc_init_hex"] = fmt("%04X", init);
      o["link_packets_crc_ok"] = n_ok;
      if (ok) {
        o["frame_hex"] = hexbits(f.data(), static_cast<size_t>(bestE));
        o["crc_hex"] = fmt("%04X", fld[i][bestE]);
      }
    }
  }
  return out;
}

}  // namespace si::rx::rc_link::esb
