// python receivers/rc_link/elrs.py (see there for the packet layouts and sources).
#include "receivers/rc_link/elrs.hpp"

#include <algorithm>
#include <array>
#include <cstdio>
#include <map>
#include <set>
#include <string>

#include "receivers/rc_link/sik.hpp"

namespace si::rx::rc_link::elrs {
namespace {

std::array<int, 256> table(int bits, int poly) {
  std::array<int, 256> t{};
  for (int i = 0; i < 256; ++i) {
    int64_t c = static_cast<int64_t>(i) << (bits - 8);
    for (int k = 0; k < 8; ++k) c = (c & (int64_t(1) << (bits - 1))) ? ((c << 1) ^ poly) : (c << 1);
    t[i] = static_cast<int>(c & ((int64_t(1) << bits) - 1));
  }
  return t;
}
const std::array<int, 256>& t14() { static const auto t = table(14, 0x2E57); return t; }
const std::array<int, 256>& t16() { static const auto t = table(16, 0x3D65); return t; }

Json sync_fields(const uint8_t* b, int ver) {
  Json f = Json::object();
  f["fhss_index"] = b[0];
  f["nonce"] = b[1];
  if (ver >= 4) {
    const int v = b[3];
    f["rate_enum"] = b[2];
    f["switch_enc_mode"] = v & 1;
    f["tlm_ratio_code"] = (v >> 1) & 7;
    f["gemini_mode"] = (v >> 4) & 1;
    f["ota_protocol"] = (v >> 5) & 3;
    f["uid4"] = b[4];
    f["uid5"] = b[5];
  } else {
    const int v = b[2];
    f["switch_enc_mode"] = v & 1;
    f["tlm_ratio_code"] = (v >> 1) & 7;
    f["rate_index"] = v >> 4;
    f["uid3"] = b[3];
    f["uid4"] = b[4];
    f["uid5"] = b[5];
  }
  return f;
}

std::string hexs(const std::vector<uint8_t>& v) {
  std::string s;
  char b[3];
  for (uint8_t x : v) std::snprintf(b, sizeof b, "%02x", x), s += b;
  return s;
}
}  // namespace

int crc(const std::vector<uint8_t>& data, int init, int bits) {
  const auto& t = bits == 14 ? t14() : t16();
  int64_t c = init;
  for (uint8_t b : data) c = (c << 8) ^ t[static_cast<size_t>(((c >> (bits - 8)) ^ b) & 0xFF)];
  return static_cast<int>(c & ((int64_t(1) << bits) - 1));
}

Json parse(const std::vector<uint8_t>& p) {
  if (p.size() != 8 && p.size() != 13) return nullptr;
  static const char* kTypes[4] = {"RCDATA / LINKSTATS", "DATA", "SYNC", "?"};
  const int typ = p[0] & 3;
  Json out = Json::object();
  out["protocol"] = "ExpressLRS";
  out["ota"] = p.size() == 8 ? "OTA4" : "OTA8";
  out["packet_type"] = kTypes[typ];
  out["payload_hex"] = hexs(p);
  out["crc_ok"] = nullptr;
  if (typ != 2) return out;
  int rx, bits;
  std::vector<uint8_t> data;
  if (p.size() == 8) {
    rx = ((p[0] >> 2) << 8) | p[7];
    data.push_back(static_cast<uint8_t>(p[0] & 3));
    data.insert(data.end(), p.begin() + 1, p.begin() + 7);
    bits = 14;
  } else {
    rx = p[11] | (p[12] << 8);
    data.assign(p.begin(), p.begin() + 11);
    bits = 16;
  }
  for (int ver : {4, 3}) {
    const Json f = sync_fields(p.data() + 1, ver);
    const int u = (f["uid4"].get<int>() << 8) | f["uid5"].get<int>();
    const std::pair<int, bool> inits[2] = {{u ^ (ver << 8), false}, {ver >= 4 ? ver : 0, true}};
    for (const auto& [init, binding] : inits)
      if (crc(data, init, bits) == rx) {
        out["crc_ok"] = true;
        out["ota_version"] = ver;
        out["binding"] = binding;
        out["sync"] = f;
        return out;
      }
  }
  out["crc_ok"] = false;
  out["note"] = "SYNC CRC not verified with its own UID bytes (model match on, other OTA version, or bit errors)";
  return out;
}

uint32_t fhss_seed(int uid2, int uid3, int uid4, int uid5, int ver) {
  return (static_cast<uint32_t>(uid2) << 24) | (static_cast<uint32_t>(uid3) << 16) | (static_cast<uint32_t>(uid4) << 8) |
         static_cast<uint32_t>(uid5 ^ ver);
}

std::vector<int> fhss_sequence(uint32_t seed, int count) {
  const int sync = count / 2, n = (256 / count) * count;
  std::vector<int> seq(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i) seq[static_cast<size_t>(i)] = i % count == 0 ? sync : i % count == sync ? 0 : i % count;
  uint64_t st = seed;
  for (int i = 0; i < n; ++i) {
    if (i % count == 0) continue;
    st = (214013ULL * st + 2531011ULL) % (1ULL << 31);
    const int r = static_cast<int>((st >> 16) % static_cast<uint64_t>(count - 1)) + 1;
    std::swap(seq[static_cast<size_t>(i)], seq[static_cast<size_t>((i / count) * count + r)]);
  }
  return seq;
}

std::vector<std::pair<int, int>> uid23_candidates(int uid4, int uid5, int count,
                                                  const std::vector<std::pair<int64_t, int>>& observations, int ver,
                                                  int uid3) {
  const int n_seq = (256 / count) * count, sync = count / 2;
  std::map<int, std::vector<std::pair<int, int>>> blocks;
  for (const auto& [i, ch] : observations) {
    const int ii = static_cast<int>(((i % n_seq) + n_seq) % n_seq);
    blocks[ii / count].push_back({ii % count, ch});
  }
  // LCG jump to each block start: st_k = a st + c (mod 2^31)
  auto jump = [](uint64_t k) {
    uint64_t a = 1, c = 0, A = 214013, C = 2531011;
    const uint64_t M = 1ULL << 31;
    while (k) {
      if (k & 1) a = (A * a) % M, c = (A * c + C) % M;
      C = (A * C + C) % M;
      A = (A * A) % M;
      k >>= 1;
    }
    return std::make_pair(a, c);
  };
  std::vector<std::pair<int, int>> out;
  std::vector<int> blk(static_cast<size_t>(count));
  for (int u3 = uid3 < 0 ? 0 : uid3; u3 <= (uid3 < 0 ? 255 : uid3); ++u3)
    for (int u2 = 0; u2 < 128; ++u2) {   // the LCG drops seed bit 31: UID2 bit 7 never changes the sequence
      const uint64_t seed = fhss_seed(u2, u3, uid4, uid5, ver) % (1ULL << 31);
      bool ok = true;
      for (const auto& [b, obs] : blocks) {
        const auto [a, c] = jump(static_cast<uint64_t>(b) * static_cast<uint64_t>(count - 1));
        uint64_t st = (a * seed + c) % (1ULL << 31);
        for (int j = 0; j < count; ++j) blk[static_cast<size_t>(j)] = j == 0 ? sync : j == sync ? 0 : j;
        for (int j = 1; j < count; ++j) {
          st = (214013ULL * st + 2531011ULL) % (1ULL << 31);
          const int r = static_cast<int>((st >> 16) % static_cast<uint64_t>(count - 1)) + 1;
          std::swap(blk[static_cast<size_t>(j)], blk[static_cast<size_t>(r)]);
        }
        for (const auto& [j, ch] : obs) ok = ok && blk[static_cast<size_t>(j)] == ch;
        if (!ok) break;
      }
      if (ok) out.push_back({u2, u3}), out.push_back({u2 | 0x80, u3});
    }
  std::sort(out.begin(), out.end());
  return out;
}

Json recover_hops(const Json& sync, double t_sync_us, int sync_ch, double interval_us, int hop_every, int count,
                  const std::vector<std::pair<double, int>>& packets) {
  const Json& f = sync.contains("sync") && sync["sync"].is_object() ? sync["sync"] : sync;
  const int ver = sync.contains("ota_version") ? sync["ota_version"].get<int>() : 4;
  const int64_t n0 = f["nonce"].get<int64_t>(), i0 = f["fhss_index"].get<int64_t>();
  const int64_t seq_len = (256 / count) * count;
  auto index_at = [&](int64_t dn) {
    int64_t hops = 0;
    if (dn >= 0) {
      for (int64_t n = n0; n < n0 + dn; ++n) hops += ((n + 1) % hop_every == 0);
    } else {
      for (int64_t n = n0 + dn; n < n0; ++n) hops -= (((n + 1) % hop_every + hop_every) % hop_every == 0);
    }
    return (((i0 + hops) % seq_len) + seq_len) % seq_len;
  };
  std::vector<std::pair<int64_t, int>> obs = {{i0, sync_ch}};
  for (const auto& [t, ch] : packets) {
    const int64_t dn = py_round_int((t - t_sync_us) / interval_us);
    if (dn != 0) obs.push_back({index_at(dn), ch});
  }
  if (obs.size() < 2) return nullptr;
  const int uid3 = f.contains("uid3") && f["uid3"].is_number_integer() ? f["uid3"].get<int>() : -1;
  const auto c = uid23_candidates(f["uid4"].get<int>(), f["uid5"].get<int>(), count, obs, ver, uid3);
  Json out = Json::object();
  out["observations"] = obs.size();
  out["uid23_candidates"] = c.size();
  Json ids = Json::array();
  char b[8];
  for (size_t k = 0; k < c.size() && k < 8; ++k) std::snprintf(b, sizeof b, "%02X%02X", c[k].first, c[k].second), ids.push_back(b);
  out["uid2_uid3"] = ids;
  std::set<std::vector<int>> seqs;
  for (const auto& [a, u3] : c) seqs.insert(fhss_sequence(fhss_seed(a, u3, f["uid4"].get<int>(), f["uid5"].get<int>(), ver), count));
  if (seqs.size() == 1) {
    const auto& s = *seqs.begin();
    Json nx = Json::array();
    for (int k = 1; k <= 8; ++k) nx.push_back(s[static_cast<size_t>((i0 + k) % seq_len)]);
    out["next_channels"] = nx;
  }
  return out;
}

// ---------------------------------------------------------------- LR1121 GFSK modes (python gfsk_*; PLAN G3)
std::vector<uint8_t> gfsk_whiten(const std::vector<uint8_t>& b, const char* variant) {
  const std::string v = variant;
  if (v == "none") return b;
  static const std::vector<uint8_t> pn = [] {
    std::vector<uint8_t> out;
    unsigned st = 0x1FF;
    for (int i = 0; i < 16; ++i) {
      out.push_back(static_cast<uint8_t>(st & 0xFF));
      for (int k = 0; k < 8; ++k) st = (st >> 1) | (((st ^ (st >> 5)) & 1u) << 8);
    }
    return out;
  }();
  std::vector<uint8_t> out(b.size());
  for (size_t i = 0; i < b.size(); ++i) {
    uint8_t w = pn[i];
    if (v == "pn9_lsb") {
      uint8_t r = 0;
      for (int k = 0; k < 8; ++k) r = static_cast<uint8_t>((r << 1) | ((w >> k) & 1));
      w = r;
    }
    out[i] = b[i] ^ w;
  }
  return out;
}

Json gfsk_try_decode(const cfvec& y, double ofs, double rate) {
  static const char* kVariants[] = {"pn9_msb", "pn9_lsb", "none"};
  for (const auto& [ph, bits] : sik::demod_bits(y, ofs, rate))
    for (bool inv : {false, true}) {
      std::vector<uint8_t> b(bits.size());
      for (size_t i = 0; i < b.size(); ++i) b[i] = inv ? 1 - bits[i] : bits[i];
      std::vector<int64_t> edges = {-1};
      for (size_t i = 0; i + 1 < b.size(); ++i)
        if (b[i] == b[i + 1]) edges.push_back(static_cast<int64_t>(i));
      edges.push_back(static_cast<int64_t>(b.size()) - 1);
      for (size_t e = 0; e + 1 < edges.size(); ++e) {
        const int64_t a1 = edges[e + 1];
        if (a1 - edges[e] < 12) continue;
        for (int64_t s0 = a1 - 1; s0 < a1 + 3; ++s0)
          for (int L : {8, 13}) {
            if (s0 < 0 || static_cast<size_t>(s0) + 16 + 8 * L > b.size()) continue;
            std::vector<uint8_t> raw(static_cast<size_t>(2 + L));
            for (size_t k = 0; k < raw.size(); ++k) {
              unsigned v = 0;
              for (int i = 0; i < 8; ++i) v = (v << 1) | b[static_cast<size_t>(s0) + 8 * k + i];
              raw[k] = static_cast<uint8_t>(v);
            }
            const std::vector<uint8_t> pay(raw.begin() + 2, raw.end());
            for (const char* wv : kVariants) {
              Json p = parse(gfsk_whiten(pay, wv));
              if (p.is_null() || p["packet_type"] != "SYNC" || !(p["crc_ok"].is_boolean() && p["crc_ok"].get<bool>()))
                continue;
              const Json& f = p["sync"];
              if (raw[0] != f["uid5"].get<int>() || raw[1] != f["uid4"].get<int>()) continue;
              char sw[8];
              std::snprintf(sw, sizeof sw, "0x%02X%02X", raw[0], raw[1]);
              Json out = Json::object();
              out["protocol"] = "ExpressLRS GFSK (LR1121)";
              out["sync_word"] = sw;
              out["uid4"] = f["uid4"];
              out["uid5"] = f["uid5"];
              out["whitening"] = wv;
              out["ota"] = p["ota"];
              out["ota_version"] = p["ota_version"];
              out["binding"] = p["binding"];
              out["sync"] = f;
              out["polarity_inverted"] = inv;
              out["crc_ok"] = true;
              return out;
            }
          }
      }
    }
  return nullptr;
}

}  // namespace si::rx::rc_link::elrs
