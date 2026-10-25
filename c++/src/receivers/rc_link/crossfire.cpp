// python receivers/rc_link/crossfire.py (see there for the sources and the packet).
#include "receivers/rc_link/crossfire.hpp"

#include <cmath>
#include <cstdio>

#include "core/types.hpp"
#include "receivers/rc_link/sik.hpp"

namespace si::rx::rc_link::crossfire {
namespace {
constexpr uint32_t kSync = 0x3E6B3E;
struct Grid { const char* band; double f0, step; int nmax; };
constexpr Grid kGrids[] = {{"868", 860.165e6, 260e3, 100}, {"915", 902.165e6, 260e3, 100}, {"868 CE", 863.093e6, 112e3, 50}};

std::vector<std::pair<int64_t, bool>> find_sync(const std::vector<uint8_t>& bits, int max_err = 2) {
  std::vector<std::pair<int64_t, bool>> hits;
  const int64_t n = static_cast<int64_t>(bits.size());
  if (n < 48) return hits;
  int pat[24];
  for (int i = 0; i < 24; ++i) pat[i] = ((kSync >> (23 - i)) & 1) ? 1 : -1;
  std::vector<int> c(static_cast<size_t>(n - 23));
  for (int64_t p = 0; p + 24 <= n; ++p) {
    int v = 0;
    for (int i = 0; i < 24; ++i) v += (bits[p + i] ? 1 : -1) * pat[i];
    c[static_cast<size_t>(p)] = v;
  }
  for (bool inv : {false, true})
    for (int64_t p = 0; p < static_cast<int64_t>(c.size()); ++p) {
      const int cc = inv ? -c[static_cast<size_t>(p)] : c[static_cast<size_t>(p)];
      if (cc < 24 - 2 * max_err || p < 8) continue;
      int tr = 0;
      for (int64_t j = p - 8; j + 1 < p; ++j) tr += bits[j] != bits[j + 1];
      if (tr >= 6) hits.emplace_back(p + 24, inv);
    }
  return hits;
}

std::string hexs(const std::vector<uint8_t>& v) {
  std::string s;
  char b[3];
  for (uint8_t x : v) std::snprintf(b, sizeof b, "%02x", x), s += b;
  return s;
}
}  // namespace

Json grid_channel(double f_hz, double tol_hz) {
  double best = -1;
  Json out = nullptr;
  for (const auto& g : kGrids) {
    const int64_t n = py_round_int((f_hz - g.f0) / g.step);
    if (n < 0 || n >= g.nmax) continue;
    const double err = std::fabs(f_hz - (g.f0 + static_cast<double>(n) * g.step));
    if (err > tol_hz || (best >= 0 && err >= best)) continue;
    best = err;
    const bool ce = std::string(g.band) == "868 CE";
    const std::string link = ce ? "unknown" : (n < 50 ? "uplink" : "downlink");
    out = {{"band", g.band}, {"channel", (n < 50 || ce) ? n : n - 50}, {"link", link}};
  }
  return out;
}

Json try_decode(const cfvec& y, double ofs, double f_hz) {
  const Json ch = grid_channel(f_hz);
  const size_t nbytes = (ch.is_object() && ch["link"] == "downlink") ? 13 : 23;
  for (const auto& [ph, bits] : sik::demod_bits(y, ofs, kRate))
    for (const auto& [pos, inv] : find_sync(bits)) {
      if (static_cast<size_t>(pos) + 8 * nbytes > bits.size()) continue;
      std::vector<uint8_t> raw(nbytes), dw(nbytes);
      for (size_t b = 0; b < nbytes; ++b) {
        unsigned v = 0;
        for (int i = 0; i < 8; ++i) v = (v << 1) | static_cast<unsigned>(bits[static_cast<size_t>(pos) + 8 * b + i] ^ (inv ? 1 : 0));
        raw[b] = static_cast<uint8_t>(v);
        dw[b] = raw[b] ^ sik::whiten_pn9_msb()[b];
      }
      Json out = Json::object();
      out["sync"] = "0x3E6B3E";
      out["polarity_inverted"] = inv;
      out["packet_bytes"] = nbytes;
      out["packet_hex"] = hexs(raw);
      out["dewhitened_hex"] = hexs(dw);
      out["crc_ok"] = nullptr;
      out["crc_note"] = "per-hop CRC init unknown: packet not verifiable (research)";
      if (ch.is_object()) out["grid"] = ch;
      return out;
    }
  return nullptr;
}

}  // namespace si::rx::rc_link::crossfire
