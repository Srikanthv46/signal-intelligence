// python receivers/rc_link/hott.py (see there for the sources, the PHY, the FEC and the packets).
#include "receivers/rc_link/hott.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <string>

#include "core/pyfmt.hpp"
#include "receivers/rc_link/cc2500.hpp"
#include "receivers/rc_link/radiolink.hpp"
#include "receivers/rc_link/sik.hpp"

namespace si::rx::rc_link::hott {
namespace {
constexpr unsigned kSyncTx = 0xD391, kSyncTlm = 0x2C6E;
constexpr uint16_t kTags[] = {0xC06B, 0xC34A, 0xDB24, 0x8E09, 0x272E, 0x217F, 0x155B, 0xEDE8,
                              0x1D31, 0x0986, 0x56F7, 0x6454, 0xC42D, 0x01D2, 0xC253, 0x1180};
constexpr int kT[16] = {0, 3, 1, 2, 3, 0, 2, 1, 3, 0, 2, 1, 0, 3, 1, 2};

std::string hexs(const uint8_t* p, size_t n) {
  std::string s;
  char b[3];
  for (size_t i = 0; i < n; ++i) std::snprintf(b, sizeof b, "%02x", p[i]), s += b;
  return s;
}

std::string hx(unsigned v, int w) {
  char b[12];
  std::snprintf(b, sizeof b, "0x%0*X", w, v);
  return b;
}

Json parse_tx(const uint8_t* p) {
  const unsigned tag = p[0] | (p[1] << 8);
  const bool bind = std::all_of(p + 40, p + 45, [](uint8_t v) { return v == 0xFA; });
  const unsigned grp = p[3];
  Json ch = Json::array();
  for (int i = 0; i < 12; ++i) ch.push_back(py_round((p[4 + 2 * i] | (p[5 + 2 * i] << 8)) / 2.0, 1));
  int table = -1;
  for (int k = 0; k < 16; ++k)
    if (kTags[k] == tag) table = k;
  Json out = Json::object();
  out["protocol"] = "Graupner HoTT";
  out["packet_name"] = bind ? "bind" : ((grp & 0x40) ? "failsafe" : "data");
  out["hop_table_tag"] = hx(tag, 4);
  out["hop_table"] = table >= 0 ? Json(table) : Json(nullptr);
  out["hop_index"] = p[2];
  out["channel_group"] = (grp & 1) ? "CH1-8 + CH13-16" : "CH1-8 + CH9-12";
  out["channels_us"] = ch;
  out["sensor_byte"] = hx(p[28], 2);
  out["mode_byte"] = hx(p[29], 2);
  out["tx_id"] = bind ? hexs(p + 45, 5) : hexs(p + 40, 5);
  out["rx_id"] = bind ? Json(nullptr) : Json(hexs(p + 45, 5));
  out["crc_ok"] = true;
  if (!std::all_of(p + 30, p + 39, [](uint8_t v) { return v == 0xFF; }) || p[39] != 0x07)
    out["note"] = "constant bytes 30-39 differ from the documented FF x 9, 07";
  return out;
}

Json parse_tlm(const uint8_t* p) {
  Json out = Json::object();
  out["protocol"] = "Graupner HoTT receiver telemetry";
  out["packet_name"] = "telemetry";
  out["tx_id"] = hexs(p, 5);
  out["rx_id"] = hexs(p + 5, 5);
  out["warnings"] = p[10];
  out["page"] = p[11];
  out["crc_ok"] = true;
  if (p[11] == 0) {
    out["rx_voltage_v"] = py_round(p[13] / 10.0, 1);
    out["temperature_c"] = static_cast<int>(p[14]) - 20;
    out["rx_rssi_raw"] = p[15];
    out["lqi_pct"] = p[16];
  }
  return out;
}
}  // namespace

int coded_len(int n_data) {
  const int n = n_data + 2;
  return 2 * 8 * (n + (n % 2 ? 1 : 2));
}

std::vector<uint8_t> fec_encode(const std::vector<uint8_t>& data) {
  std::vector<uint8_t> b = data;
  b.push_back(0x0B);
  if (data.size() % 2 == 0) b.push_back(0x0B);
  unsigned st = 0;
  std::vector<int> sym;
  for (uint8_t v : b)
    for (int i = 0; i < 8; ++i) {
      const unsigned x = (v >> (7 - i)) & 1u;
      sym.push_back(kT[(st << 1) | x]);
      st = ((st << 1) | x) & 7u;
    }
  std::vector<uint8_t> out;
  for (size_t k = 0; k + 16 <= sym.size(); k += 16)
    for (int j = 0; j < 16; ++j) {
      const int v = sym[k + static_cast<size_t>(4 * (3 - (j % 4)) + (3 - j / 4))];
      out.push_back(static_cast<uint8_t>((v >> 1) & 1));
      out.push_back(static_cast<uint8_t>(v & 1));
    }
  return out;
}

std::pair<std::vector<uint8_t>, int> fec_decode(const std::vector<uint8_t>& bits, int n_bytes) {
  const size_t nb = bits.size() / 32;
  std::vector<int> sym(16 * nb);
  for (size_t k = 0; k < nb; ++k)
    for (int j = 0; j < 16; ++j)
      sym[16 * k + static_cast<size_t>(4 * (3 - (j % 4)) + (3 - j / 4))] = bits[32 * k + 2 * j] * 2 + bits[32 * k + 2 * j + 1];
  const size_t n = std::min(sym.size(), static_cast<size_t>(8 * n_bytes + 3));
  constexpr int64_t kInf = int64_t(1) << 30;
  std::array<int64_t, 8> pm;
  pm.fill(kInf);
  pm[0] = 0;
  std::vector<std::array<int8_t, 8>> prev(n);
  static const int pc[4] = {0, 1, 1, 2};
  for (size_t t = 0; t < n; ++t) {
    std::array<int64_t, 8> nw;
    nw.fill(kInf);
    std::array<int8_t, 8> arg{};
    for (int k = 0; k < 16; ++k) {   // from state k >> 1 with input k & 1 -> state k & 7
      const int64_t c = pm[static_cast<size_t>(k >> 1)] + pc[kT[k] ^ sym[t]];
      if (c < nw[static_cast<size_t>(k & 7)]) {
        nw[static_cast<size_t>(k & 7)] = c;
        arg[static_cast<size_t>(k & 7)] = static_cast<int8_t>(k);
      }
    }
    pm = nw;
    prev[t] = arg;
  }
  size_t st = static_cast<size_t>(std::min_element(pm.begin(), pm.end()) - pm.begin());
  const int metric = static_cast<int>(pm[st]);
  std::vector<uint8_t> ob(n);
  for (size_t t = n; t-- > 0;) {
    const int k = prev[t][st];
    ob[t] = static_cast<uint8_t>(k & 1);
    st = static_cast<size_t>(k >> 1);
  }
  std::vector<uint8_t> out(static_cast<size_t>(n_bytes), 0);
  for (size_t i = 0; i < static_cast<size_t>(8 * n_bytes) && i < n; ++i) out[i / 8] |= static_cast<uint8_t>(ob[i] << (7 - i % 8));
  return {out, metric};
}

Json try_decode(const cfvec& y, double ofs, std::optional<double> f_hz, double rate) {
  static const std::vector<uint8_t> w = radiolink::pn9(64);
  struct Kind { unsigned sync; int n; bool tx; };
  const Kind kinds[] = {{kSyncTx, kLenTx, true}, {kSyncTlm, kLenTlm, false}};
  for (const auto& [ph, bits] : sik::demod_bits_mid(y, ofs, rate))
    for (const auto& kd : kinds) {
      const int nc = coded_len(kd.n);
      const uint8_t s4[4] = {static_cast<uint8_t>(kd.sync >> 8), static_cast<uint8_t>(kd.sync & 0xFF),
                             static_cast<uint8_t>(kd.sync >> 8), static_cast<uint8_t>(kd.sync & 0xFF)};
      for (const auto& h : cc2500::find_sync(bits, s4, 4)) {
        if (static_cast<int64_t>(bits.size()) - h.pos < nc) continue;
        std::vector<uint8_t> b(bits.begin() + h.pos, bits.begin() + h.pos + nc);
        if (h.inv)
          for (auto& v : b) v ^= 1;
        auto [dec, metric] = fec_decode(b, kd.n + 2);
        for (size_t i = 0; i < dec.size(); ++i) dec[i] ^= w[i];
        if (cc2500::crc16_cc2500(dec.data(), static_cast<size_t>(kd.n)) !=
            ((dec[static_cast<size_t>(kd.n)] << 8) | dec[static_cast<size_t>(kd.n) + 1]))
          continue;
        Json d = kd.tx ? parse_tx(dec.data()) : parse_tlm(dec.data());
        d["sync"] = hx(kd.sync, 4);
        d["sync_bit_errors"] = h.err;
        d["polarity_inverted"] = h.inv;
        d["fec_corrected_bits"] = metric;
        if (f_hz) d["channel"] = static_cast<int64_t>(std::nearbyint((*f_hz - kF0) / kStep));
        return d;
      }
    }
  return nullptr;
}

}  // namespace si::rx::rc_link::hott
