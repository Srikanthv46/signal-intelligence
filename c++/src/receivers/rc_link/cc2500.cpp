// python receivers/rc_link/cc2500.py (see there for the sources and the packets).
#include "receivers/rc_link/cc2500.hpp"

#include <algorithm>
#include <array>
#include <cstdio>
#include <optional>

#include "receivers/rc_link/sik.hpp"

namespace si::rx::rc_link::cc2500 {
namespace {
constexpr uint8_t kD8Sync[] = {0xD3, 0x91};
constexpr uint8_t kD16Sync[] = {0xD3, 0x91, 0xD3, 0x91};
constexpr uint8_t kD8Len = 0x11;

const std::array<uint16_t, 256>& table() {
  static const std::array<uint16_t, 256> t = [] {
    std::array<uint16_t, 256> a{};
    for (unsigned v = 0; v < 256; ++v) {
      unsigned c = v;
      for (int i = 0; i < 8; ++i) c = (c & 1) ? (c >> 1) ^ 0x8408u : c >> 1;
      a[v] = static_cast<uint16_t>(c);   // KERMIT table: T[1] = 0x1189
    }
    return a;
  }();
  return t;
}

using Hit = SyncHit;
}  // namespace

std::vector<SyncHit> find_sync(const std::vector<uint8_t>& bits, const uint8_t* sync, int nsync) {
  std::vector<Hit> hits;
  const int L = 8 * nsync;
  const int64_t n = static_cast<int64_t>(bits.size());
  if (n < L + 8) return hits;
  std::vector<int> pat(static_cast<size_t>(L));
  for (int i = 0; i < L; ++i) pat[static_cast<size_t>(i)] = ((sync[i / 8] >> (7 - i % 8)) & 1) ? 1 : -1;
  std::vector<int> c(static_cast<size_t>(n - L + 1));
  for (int64_t p = 0; p + L <= n; ++p) {
    int v = 0;
    for (int i = 0; i < L; ++i) v += (bits[static_cast<size_t>(p + i)] ? 1 : -1) * pat[static_cast<size_t>(i)];
    c[static_cast<size_t>(p)] = v;
  }
  for (bool inv : {false, true})
    for (int64_t p = 0; p < static_cast<int64_t>(c.size()); ++p) {
      const int cc = inv ? -c[static_cast<size_t>(p)] : c[static_cast<size_t>(p)];
      if (cc < L - 4 || p < 8) continue;
      int tr = 0;
      for (int64_t j = p - 8; j + 1 < p; ++j) tr += bits[static_cast<size_t>(j)] != bits[static_cast<size_t>(j + 1)];
      if (tr >= 6) hits.push_back({p + L, inv, (L - cc) / 2});
    }
  return hits;
}

namespace {

std::string hexs(const uint8_t* p, size_t n) {
  std::string s;
  char b[3];
  for (size_t i = 0; i < n; ++i) std::snprintf(b, sizeof b, "%02x", p[i]), s += b;
  return s;
}

std::string hex_sync(const uint8_t* s, int n) {
  std::string out = "0x";
  char b[3];
  for (int i = 0; i < n; ++i) std::snprintf(b, sizeof b, "%02X", s[i]), out += b;
  return out;
}

struct Check { std::string variant; bool crc_ok; Json crc; size_t len; };

std::optional<Check> check_d8(const std::vector<uint8_t>& raw) {
  if (raw.size() < kD8Len + 3u || raw[0] != kD8Len) return std::nullopt;
  const size_t n = kD8Len + 1;
  const bool ok = crc16_cc2500(raw.data(), n) == ((raw[n] << 8) | raw[n + 1]);
  return Check{"D8", ok, ok ? Json("hardware 0x8005") : Json(nullptr), n + 2};
}

std::optional<Check> check_d16(const std::vector<uint8_t>& raw) {
  if (raw.empty() || (raw[0] != 0x1D && raw[0] != 0x20)) return std::nullopt;
  const size_t L = raw[0] + 1u;
  if (raw.size() < L) return std::nullopt;
  const bool ok = crc16_d16(raw.data() + 3, L - 5) == ((raw[L - 2] << 8) | raw[L - 1]);
  const bool hw = raw.size() >= L + 2 && crc16_cc2500(raw.data(), L) == ((raw[L] << 8) | raw[L + 1]);
  const std::string var = ok && hw ? "D16 v2" : ok ? (raw[0] == 0x1D ? "D16 v1" : "D16 LBT v1") : "D16";
  return Check{var, ok, ok ? Json(hw ? "software + hardware" : "software") : Json(nullptr), hw ? L + 2 : L};
}
}  // namespace

double nominal_rate(const std::string& variant) {
  std::string v = variant;
  if (v == "D16 LBT v1") v = "D16 LBT";
  for (const auto& r : kRates)
    if (v == r.variant) return r.bps;
  return kRates[1].bps;
}

uint16_t crc16_cc2500(const uint8_t* p, size_t n) {
  unsigned c = 0xFFFF;
  for (size_t k = 0; k < n; ++k) {
    c ^= static_cast<unsigned>(p[k]) << 8;
    for (int i = 0; i < 8; ++i) c = (c & 0x8000) ? ((c << 1) ^ 0x8005u) & 0xFFFF : (c << 1) & 0xFFFF;
  }
  return static_cast<uint16_t>(c);
}

uint16_t crc16_d16(const uint8_t* p, size_t n) {
  unsigned c = 0;
  for (size_t k = 0; k < n; ++k) c = ((c << 8) & 0xFFFF) ^ table()[((c >> 8) ^ p[k]) & 0xFF];
  return static_cast<uint16_t>(c);
}

Json try_identify(const cfvec& y, double ofs, double rate) {
  struct Var { const char* name; const uint8_t* sync; int nsync; bool d16; };
  const Var vars[] = {{"D16", kD16Sync, 4, true}, {"D8", kD8Sync, 2, false}};
  const auto phases = sik::demod_bits(y, ofs, rate);
  for (const auto& v : vars) {
    Json best = nullptr, sync_only = nullptr;
    int best_err = 1 << 30;
    for (const auto& [ph, bits] : phases)
      for (const Hit& h : find_sync(bits, v.sync, v.nsync)) {
        const size_t nb = (bits.size() - static_cast<size_t>(h.pos)) / 8;
        std::vector<uint8_t> raw(nb);
        for (size_t b = 0; b < nb; ++b) {
          unsigned x = 0;
          for (int i = 0; i < 8; ++i)
            x = (x << 1) | static_cast<unsigned>(bits[static_cast<size_t>(h.pos) + 8 * b + i] ^ (h.inv ? 1 : 0));
          raw[b] = static_cast<uint8_t>(x);
        }
        const auto r = v.d16 ? check_d16(raw) : check_d8(raw);
        if (!r) {
          if (v.d16 && h.err == 0 && !raw.empty() && sync_only.is_null()) {
            sync_only = Json::object();
            sync_only["protocol"] = "FrSky D16-like (ACCST)";
            sync_only["variant"] = "D16";
            sync_only["sync"] = "0xD391D391";
            sync_only["length"] = raw[0];
            sync_only["tx_id"] = nullptr;
            sync_only["crc_ok"] = nullptr;
            sync_only["crc"] = nullptr;
            sync_only["packet_hex"] = hexs(raw.data(), std::min<size_t>(raw.size(), raw[0] + 3u));
            sync_only["polarity_inverted"] = h.inv;
            sync_only["note"] = "32-bit sync only; length not a data frame";
          }
          continue;
        }
        char id[8];
        std::snprintf(id, sizeof id, "0x%02X%02X", raw[2], raw[1]);
        Json out = Json::object();
        out["protocol"] = r->crc_ok ? "FrSky " + r->variant + " (ACCST)" : std::string("FrSky ") + v.name + "-like (ACCST)";
        out["variant"] = r->variant;
        out["sync"] = hex_sync(v.sync, v.nsync);
        out["length"] = raw[0];
        out["tx_id"] = id;
        out["crc_ok"] = r->crc_ok;
        out["crc"] = r->crc;
        out["packet_hex"] = hexs(raw.data(), r->len);
        out["polarity_inverted"] = h.inv;
        if (r->crc_ok) return out;
        if (h.err < best_err) best_err = h.err, best = std::move(out);
      }
    if (!best.is_null()) return best;   // D16 sync contains D8's: do not fall through
    if (!sync_only.is_null()) return sync_only;
  }
  return nullptr;
}

}  // namespace si::rx::rc_link::cc2500
