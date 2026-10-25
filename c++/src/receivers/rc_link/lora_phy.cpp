#include "receivers/rc_link/lora_phy.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <utility>

#include "core/mavlink.hpp"
#include "core/pybytes.hpp"
#include "core/pyfmt.hpp"

namespace si::rx::rc_link::lora_phy {
namespace {

int64_t pymod(int64_t a, int64_t b) {
  int64_t m = a % b;
  return m < 0 ? m + b : m;
}

using Bits = std::vector<int>;

Bits bits_msb(int64_t v, int n) {
  Bits b(n);
  for (int i = 0; i < n; ++i) b[i] = static_cast<int>((v >> (n - 1 - i)) & 1);
  return b;
}

// python _deinterleave: cw_len symbols of sf_app bits -> sf_app codewords of cw_len bits
std::vector<Bits> deinterleave(const std::vector<int64_t>& block, int sf_app, int cw_len) {
  std::vector<Bits> inter;
  for (auto s : block) inter.push_back(bits_msb(s, sf_app));
  std::vector<Bits> de(sf_app, Bits(cw_len, 0));
  for (int i = 0; i < cw_len; ++i)
    for (int j = 0; j < sf_app; ++j) de[pymod(i - j - 1, sf_app)][i] = inter[i][j];
  return de;
}

struct Nib {
  int nib;
  bool corrected, bad;
};

// python _hamming
Nib hamming(const Bits& cw, int cw_len) {
  int d[4] = {cw[3], cw[2], cw[1], cw[0]};
  bool corrected = false, bad = false;
  if (cw_len >= 7) {
    const int s0 = cw[0] ^ cw[1] ^ cw[2] ^ cw[4];
    const int s1 = cw[1] ^ cw[2] ^ cw[3] ^ cw[5];
    const int s2 = cw[0] ^ cw[1] ^ cw[3] ^ cw[6];
    const int syn = s0 | (s1 << 1) | (s2 << 2);
    int sum = 0;
    for (int b : cw) sum += b;
    const bool even = sum % 2 == 0;
    if (syn && (cw_len == 7 || !even)) {
      int flip = -1;
      if (syn == 5) flip = 3;
      else if (syn == 7) flip = 2;
      else if (syn == 3) flip = 1;
      else if (syn == 6) flip = 0;
      if (flip >= 0) {
        d[flip] ^= 1;
        corrected = true;
      } else {
        bad = true;
      }
    } else if (syn) {
      bad = true;
    }
  } else if (cw_len == 6) {
    bad = (cw[4] ^ cw[0] ^ cw[1] ^ cw[2]) || (cw[5] ^ cw[1] ^ cw[2] ^ cw[3]);
  } else if (cw_len == 5) {
    int sum = 0;
    for (int b : cw) sum += b;
    bad = sum % 2;
  }
  return {(d[0] << 3) | (d[1] << 2) | (d[2] << 1) | d[3], corrected, bad};
}

struct Block {
  std::vector<int> nib;
  int corr = 0, bad = 0;
};

// python _block_nibbles
Block block_nibbles(const std::vector<int64_t>& syms, int sf, int sf_app, int cw_len) {
  const int64_t N = int64_t(1) << sf;
  std::vector<int64_t> g(syms.size());
  for (size_t i = 0; i < syms.size(); ++i) {
    int64_t v = pymod(syms[i] - 1, N);
    if (sf_app < sf) v >>= 2;
    g[i] = v ^ (v >> 1);
  }
  Block b;
  for (const auto& cw : deinterleave(g, sf_app, cw_len)) {
    Nib n = hamming(cw, cw_len);
    b.nib.push_back(n.nib);
    b.corr += n.corrected;
    b.bad += n.bad;
  }
  return b;
}

std::string hex(const uint8_t* p, size_t n) {
  static const char* d = "0123456789abcdef";
  std::string s;
  for (size_t i = 0; i < n; ++i) {
    s += d[p[i] >> 4];
    s += d[p[i] & 15];
  }
  return s;
}

uint64_t le(const std::vector<uint8_t>& p, size_t a, size_t b) {
  uint64_t v = 0;
  for (size_t i = b; i-- > a;) v = (v << 8) | p[i];
  return v;
}

}  // namespace

std::optional<Implicit> decode_implicit(const std::vector<int64_t>& raw, int sf, int cr, int length, bool ldro) {
  if (sf < 6 || raw.size() < 8 || cr < 1 || cr > 4) return std::nullopt;
  Block first = block_nibbles(std::vector<int64_t>(raw.begin(), raw.begin() + 8), sf, sf - 2, 8);
  std::vector<int> nib = first.nib;
  int corr = first.corr, bad = first.bad;
  const int sf_app = ldro ? sf - 2 : sf, cw_len = 4 + cr;
  size_t j = 8;
  while (static_cast<int>(nib.size()) < 2 * length && j + static_cast<size_t>(cw_len) <= raw.size()) {
    Block b = block_nibbles(std::vector<int64_t>(raw.begin() + static_cast<ptrdiff_t>(j),
                                                 raw.begin() + static_cast<ptrdiff_t>(j) + cw_len), sf, sf_app, cw_len);
    nib.insert(nib.end(), b.nib.begin(), b.nib.end());
    corr += b.corr;
    bad += b.bad;
    j += static_cast<size_t>(cw_len);
  }
  if (static_cast<int>(nib.size()) < 2 * length) return std::nullopt;
  const auto wh = whitening_sequence(static_cast<size_t>(length));
  Implicit out;
  out.payload.resize(static_cast<size_t>(length));
  for (int i = 0; i < length; ++i)
    out.payload[i] = static_cast<uint8_t>(((nib[2 * i + 1] << 4) | nib[2 * i]) ^ wh[static_cast<size_t>(i)]);
  out.cr = "4/" + std::to_string(cw_len);
  out.fec_corrected = corr;
  out.fec_detected_errors = bad;
  out.symbols_used = static_cast<int64_t>(j);
  return out;
}

std::vector<uint8_t> whitening_sequence(size_t n) {
  std::vector<uint8_t> out;
  unsigned r = 0xFF;
  for (size_t i = 0; i < n; ++i) {
    out.push_back(static_cast<uint8_t>(r));
    const unsigned fb = ((r >> 7) ^ (r >> 5) ^ (r >> 4) ^ (r >> 3)) & 1;
    r = ((r << 1) | fb) & 0xFF;
  }
  return out;
}

uint16_t crc16(const uint8_t* data, size_t n) {
  unsigned crc = 0;
  for (size_t i = 0; i < n; ++i)
    for (int j = 7; j >= 0; --j) {
      const unsigned fb = ((crc >> 15) & 1) ^ ((data[i] >> j) & 1);
      crc = (crc << 1) & 0xFFFF;
      if (fb) crc ^= 0x1021;
    }
  return static_cast<uint16_t>(crc);
}

int header_checksum(int n0, int n1, int n2) {
  auto b = [](int v, int k) { return (v >> k) & 1; };
  const int c4 = b(n0, 3) ^ b(n0, 2) ^ b(n0, 1) ^ b(n0, 0);
  const int c3 = b(n0, 3) ^ b(n1, 3) ^ b(n1, 2) ^ b(n1, 1) ^ b(n2, 0);
  const int c2 = b(n0, 2) ^ b(n1, 3) ^ b(n1, 0) ^ b(n2, 3) ^ b(n2, 1);
  const int c1 = b(n0, 1) ^ b(n1, 2) ^ b(n1, 0) ^ b(n2, 2) ^ b(n2, 1) ^ b(n2, 0);
  const int c0 = b(n0, 0) ^ b(n1, 1) ^ b(n2, 3) ^ b(n2, 2) ^ b(n2, 1) ^ b(n2, 0);
  return (c4 << 4) | (c3 << 3) | (c2 << 2) | (c1 << 1) | c0;
}

int64_t n_symbols(int sf, int length, int cr, int crc, bool ldro) {
  const int num = 8 * length - 4 * sf + 28 + 16 * crc;
  const double q = static_cast<double>(num) / static_cast<double>(4 * (sf - 2 * (ldro ? 1 : 0)));
  return 8 + std::max<int64_t>(static_cast<int64_t>(std::ceil(q)) * (cr + 4), 0);
}

bool default_ldro(int sf, double bw) { return static_cast<double>(int64_t(1) << sf) / bw > 16e-3; }

std::optional<Json> decode(const std::vector<int64_t>& raw, int sf, double bw, std::optional<bool> ldro_opt) {
  const bool ldro = ldro_opt ? *ldro_opt : default_ldro(sf, bw);
  if (sf < 7 || raw.size() < 8) return std::nullopt;
  Block hdr = block_nibbles(std::vector<int64_t>(raw.begin(), raw.begin() + 8), sf, sf - 2, 8);
  const int n0 = hdr.nib[0], n1 = hdr.nib[1], n2 = hdr.nib[2], n3 = hdr.nib[3], n4 = hdr.nib[4];
  if ((((n3 & 1) << 4) | n4) != header_checksum(n0, n1, n2)) return std::nullopt;
  const int length = (n0 << 4) | n1, crc_on = n2 & 1, cr = n2 >> 1;
  if (!(1 <= cr && cr <= 4) || length == 0) return std::nullopt;
  std::vector<int> nib(hdr.nib.begin() + 5, hdr.nib.end());
  int corr = hdr.corr, bad = hdr.bad;
  const size_t need = 2 * static_cast<size_t>(length) + 4 * static_cast<size_t>(crc_on);
  const int sf_app = ldro ? sf - 2 : sf;
  const int cw_len = 4 + cr;
  size_t j = 8;
  while (nib.size() < need && j + cw_len <= raw.size()) {
    Block b = block_nibbles(std::vector<int64_t>(raw.begin() + j, raw.begin() + j + cw_len), sf, sf_app, cw_len);
    nib.insert(nib.end(), b.nib.begin(), b.nib.end());
    corr += b.corr;
    bad += b.bad;
    j += cw_len;
  }
  Json out = Json::object();
  out["header_ok"] = true;
  out["length"] = length;
  out["cr"] = "4/" + std::to_string(4 + cr);
  out["crc_present"] = static_cast<bool>(crc_on);
  out["ldro"] = ldro;
  out["n_symbols"] = n_symbols(sf, length, cr, crc_on, ldro);
  out["symbols_used"] = j;
  if (nib.size() < 2 * static_cast<size_t>(length)) {
    out["complete"] = false;
    out["reason"] = "packet ends after " + std::to_string(nib.size() / 2) + " of " + std::to_string(length) + " bytes";
    return out;
  }
  std::vector<uint8_t> pay(length);
  const auto wh = whitening_sequence(length);
  for (int i = 0; i < length; ++i) pay[i] = static_cast<uint8_t>(((nib[2 * i + 1] << 4) | nib[2 * i]) ^ wh[i]);
  out["complete"] = true;
  out["payload_hex"] = hex(pay.data(), pay.size());
  out["fec_corrected"] = corr;
  out["fec_detected_errors"] = bad;
  if (crc_on) {
    if (nib.size() < need) {
      out["crc_ok"] = nullptr;
      out["reason"] = "CRC bytes missing";
      return out;
    }
    const size_t L2 = 2 * static_cast<size_t>(length);
    const unsigned rx = ((nib[L2 + 1] << 4) | nib[L2]) | (((nib[L2 + 3] << 4) | nib[L2 + 2]) << 8);
    unsigned c;
    if (length >= 2)
      c = crc16(pay.data(), length - 2) ^ pay[length - 1] ^ (static_cast<unsigned>(pay[length - 2]) << 8);
    else
      c = crc16(nullptr, 0) ^ pay[length - 1];
    out["crc_ok"] = c == rx;
  } else {
    out["crc_ok"] = nullptr;
  }
  return out;
}

int64_t n_symbols_sf5(int length, int cr, bool crc, int sf) {
  const double q = static_cast<double>(8 * length + 16 * (crc ? 1 : 0) - 4 * sf) / static_cast<double>(4 * sf);
  return SF5_EXTRA_SYMS + 8 + std::max<int64_t>(static_cast<int64_t>(std::ceil(q)), 0) * (cr + 4);
}

std::optional<Json> decode_sf5_implicit(const std::vector<int64_t>& raw, int64_t n_pay, int sf) {
  if (sf != 5 || raw.size() < static_cast<size_t>(SF5_EXTRA_SYMS + 8)) return std::nullopt;
  const Block first =
      block_nibbles(std::vector<int64_t>(raw.begin() + SF5_EXTRA_SYMS, raw.begin() + SF5_EXTRA_SYMS + 8), sf, sf, 8);
  if (first.bad) return std::nullopt;
  std::vector<Json> passes;
  for (int cr = 1; cr <= 4; ++cr) {
    const int cw_len = 4 + cr;
    std::vector<int> nib = first.nib;
    int corr = first.corr, bad = 0;
    int64_t j = SF5_EXTRA_SYMS + 8;
    for (int length = 1; length < 256; ++length) {
      const int64_t ns = n_symbols_sf5(length, cr, true, sf);
      if (ns > std::min<int64_t>(n_pay + 8, static_cast<int64_t>(raw.size()))) break;
      if (ns < n_pay - 8) continue;
      while (j < ns) {  // extend the nibble stream to ns symbols
        Block b = block_nibbles(std::vector<int64_t>(raw.begin() + j, raw.begin() + j + cw_len), sf, sf, cw_len);
        nib.insert(nib.end(), b.nib.begin(), b.nib.end());
        corr += b.corr;
        bad += b.bad;
        j += cw_len;
      }
      if (bad || nib.size() < 2 * static_cast<size_t>(length) + 4) continue;
      std::vector<uint8_t> pay(length);
      const auto wh = whitening_sequence(length);
      for (int i = 0; i < length; ++i) pay[i] = static_cast<uint8_t>(((nib[2 * i + 1] << 4) | nib[2 * i]) ^ wh[i]);
      const size_t L2 = 2 * static_cast<size_t>(length);
      const unsigned rx = ((nib[L2 + 1] << 4) | nib[L2]) | (((nib[L2 + 3] << 4) | nib[L2 + 2]) << 8);
      const unsigned c = length >= 2
                             ? crc16(pay.data(), length - 2) ^ pay[length - 1] ^ (static_cast<unsigned>(pay[length - 2]) << 8)
                             : crc16(nullptr, 0) ^ pay[length - 1];
      if (c != rx) continue;
      Json out = Json::object();
      out["header_ok"] = false;
      out["implicit_header"] = true;
      out["length"] = length;
      out["cr"] = "4/" + std::to_string(cw_len);
      out["crc_present"] = true;
      out["ldro"] = false;
      out["n_symbols"] = ns;
      out["symbols_used"] = ns;
      out["complete"] = true;
      out["payload_hex"] = hex(pay.data(), pay.size());
      out["fec_corrected"] = corr;
      out["fec_detected_errors"] = 0;
      out["crc_ok"] = true;
      passes.push_back(out);
    }
  }
  if (passes.empty()) return std::nullopt;
  for (const auto& p : passes)
    if (p["payload_hex"] != passes[0]["payload_hex"]) return std::nullopt;  // ambiguous payload
  size_t bi = 0;
  auto key = [&](const Json& p) {
    return std::make_pair(std::llabs(p["n_symbols"].get<int64_t>() - n_pay), p["fec_corrected"].get<int>());
  };
  for (size_t i = 1; i < passes.size(); ++i)
    if (key(passes[i]) < key(passes[bi])) bi = i;
  Json best = passes[bi];
  bool amb = false;
  for (const auto& p : passes) amb = amb || p["cr"] != best["cr"];
  best["cr_ambiguous"] = amb;
  return best;
}

std::optional<int> sync_word(const std::vector<int64_t>& sync_vals) {
  if (sync_vals.size() != 2) return std::nullopt;
  int nib[2];
  for (int i = 0; i < 2; ++i) {
    const int64_t v = sync_vals[i];
    const int64_t q = py_round_int(static_cast<double>(v) / 8.0);
    if (std::llabs(v - 8 * q) > 1 || !(0 <= q && q < 16)) return std::nullopt;
    nib[i] = static_cast<int>(q);
  }
  return (nib[0] << 4) | nib[1];
}

Json sync_word_network(std::optional<int> sw) {
  if (!sw) return nullptr;
  if (*sw == 0x12) return "private (Semtech default, e.g. RadioLib / RadioHead / point-to-point)";
  if (*sw == 0x34) return "LoRaWAN public network";
  if (*sw == 0x2B) return "Meshtastic";
  if (*sw == 0xF1) return "FANET";
  return nullptr;
}

namespace {
// python FANET_MFR / FANET_TYPES / FANET_AIRCRAFT (research, PLAN G10)
const char* fanet_mfr(int m) {
  switch (m) {
    case 0x01: return "Skytraxx"; case 0x03: return "BitBroker"; case 0x04: return "AirWhere";
    case 0x05: return "Windline"; case 0x06: return "Burnair"; case 0x07: return "SoftRF";
    case 0x08: return "GXAirCom"; case 0x09: return "Airtribune"; case 0x0A: return "FLARM";
    case 0x0B: return "FlyBeeper"; case 0x0C: return "Leaf Vario"; case 0x0D: return "LYNKX";
    case 0x0E: return "Syride"; case 0x0F: return "FlyBrief"; case 0x10: return "alfapilot";
    case 0x11: return "FANET+ (FLARM + FANET device)"; case 0x20: return "XC Tracer"; case 0xBA: return "Burnair";
    case 0xCB: return "Cloudbuddy"; case 0xE0: return "OGN Tracker"; case 0xE4: return "4aviation";
    case 0xFA: return "various"; case 0xFB: return "Espressif base station"; case 0xFC: return "unregistered";
    case 0xFD: return "unregistered";
    default: return nullptr;
  }
}
const char* kFanetTypes[] = {"ACK", "tracking", "name", "message", "service", "landmarks", "remote config",
                             "ground tracking", "HW info", "thermal", "HW info v2"};
const char* kFanetAircraft[] = {"other", "paraglider", "hang glider", "balloon", "glider", "powered aircraft",
                                "helicopter", "UAV"};
int64_t fs24(const uint8_t* b) {
  const int64_t v = b[0] | (b[1] << 8) | (b[2] << 16);
  return (v & 0x800000) ? v - (int64_t(1) << 24) : v;
}
}  // namespace

std::optional<Json> fanet_parse(const std::vector<uint8_t>& p, bool strict) {
  if (p.size() < 4) return std::nullopt;
  const bool ext = (p[0] & 0x80) != 0, fwd = (p[0] & 0x40) != 0;
  const int typ = p[0] & 0x3F, mfr = p[1];
  const int dev = p[2] | (p[3] << 8);
  if (typ > 10 || mfr == 0x00 || mfr == 0xFF || dev == 0 || dev == 0xFFFF) return std::nullopt;
  size_t k = 4;
  Json out = Json::object();
  out["protocol"] = "FANET";
  out["message_type"] = kFanetTypes[typ];
  out["forward"] = fwd;
  out["manufacturer_id"] = fmt("0x%02X", mfr);
  out["manufacturer"] = fanet_mfr(mfr) ? fanet_mfr(mfr) : "unknown";
  out["device_id"] = fmt("%02X:%04X", mfr, dev);
  if (ext) {
    if (p.size() < 5) return std::nullopt;
    const int e = p[4];
    k = 5;
    out["ack"] = (e >> 6) & 3;
    out["unicast"] = (e & 0x20) != 0;
    out["signed"] = (e & 0x10) != 0;
    if (e & 0x20) {
      if (p.size() < k + 3) return std::nullopt;
      out["destination"] = fmt("%02X:%04X", p[k], p[k + 1] | (p[k + 2] << 8));
      k += 3;
    }
    if (e & 0x10) k += 4;
  }
  const std::vector<uint8_t> body(p.begin() + static_cast<std::ptrdiff_t>(std::min(k, p.size())), p.end());
  if (strict && (typ != 1 || !fanet_mfr(mfr) || ext || body.size() < 11 || body.size() > 13)) return std::nullopt;
  if (typ == 1 || typ == 7 || typ == 4) {
    const size_t off = typ == 4 ? 1 : 0;
    if (body.size() < off + 6) return std::nullopt;
    const double lat = static_cast<double>(fs24(body.data() + off)) / 93206, lon = static_cast<double>(fs24(body.data() + off + 3)) / 46603;
    if (std::fabs(lat) > 90 || std::fabs(lon) > 180) return std::nullopt;
    out["lat"] = py_round(lat, 5);
    out["lon"] = py_round(lon, 5);
  }
  if (typ == 1) {
    if (body.size() < 11) return std::nullopt;
    const int w = body[6] | (body[7] << 8);
    const int alt = (w & 0x7FF) * ((w & 0x800) ? 4 : 1);
    int clb = body[9] & 0x7F;
    if (clb & 0x40) clb -= 128;
    out["aircraft_type"] = kFanetAircraft[(w >> 12) & 7];
    out["uav"] = ((w >> 12) & 7) == 7;
    out["online_tracking"] = (w & 0x8000) != 0;
    out["alt_m"] = alt;
    out["speed_kmh"] = ((body[8] & 0x7F) * ((body[8] & 0x80) ? 5 : 1)) / 2.0;
    out["climb_ms"] = py_round(clb * ((body[9] & 0x80) ? 5 : 1) / 10.0, 1);
    out["heading_deg"] = py_round(body[10] * 360.0 / 256, 1);
  } else if (typ == 7 && body.size() >= 7) {
    out["ground_type"] = body[6] >> 4;
    out["online_tracking"] = (body[6] & 1) != 0;
  } else if (typ == 2) {
    out["name"] = pybytes::decode_utf8(pybytes::Bytes(body));
  }
  return out;
}

std::optional<Json> identify(std::optional<int> sw, const std::vector<uint8_t>& p) {
  if (sw && *sw == 0x34 && p.size() >= 12) {
    const int mtype = p[0] >> 5, major = p[0] & 3;
    static const char* names[] = {"join_request",        "join_accept",         "unconfirmed_data_up",
                                  "unconfirmed_data_down", "confirmed_data_up",   "confirmed_data_down",
                                  "rejoin_request",      "proprietary"};
    if (major != 0) return std::nullopt;
    Json out = Json::object();
    out["protocol"] = "LoRaWAN";
    out["mtype"] = names[mtype];
    if (mtype == 0 && p.size() == 23) {
      std::vector<uint8_t> je(p.begin() + 1, p.begin() + 9), de(p.begin() + 9, p.begin() + 17);
      std::reverse(je.begin(), je.end());
      std::reverse(de.begin(), de.end());
      out["join_eui"] = hex(je.data(), je.size());
      out["dev_eui"] = hex(de.data(), de.size());
      out["dev_nonce"] = le(p, 17, 19);
    } else if (mtype >= 2 && mtype <= 5) {
      const int fopts = p[5] & 0x0F;
      out["dev_addr"] = fmt("%08llx", static_cast<unsigned long long>(le(p, 1, 5)));
      out["adr"] = static_cast<bool>(p[5] & 0x80);
      out["ack"] = static_cast<bool>(p[5] & 0x20);
      out["fcnt"] = le(p, 6, 8);
      out["fopts_len"] = fopts;
      if (p.size() > static_cast<size_t>(8 + fopts + 4)) out["fport"] = p[8 + fopts];
    }
    return out;
  }
  if (sw && (*sw == 0xF1 || *sw == 0x12)) {
    auto f = fanet_parse(p, *sw == 0x12);
    if (f) return f;
  }
  if (sw && *sw == 0x2B && p.size() >= 16) {
    const uint64_t to = le(p, 0, 4), frm = le(p, 4, 8), pid = le(p, 8, 12);
    const int fl = p[12];
    Json out = Json::object();
    out["protocol"] = "Meshtastic";
    out["to"] = to == 0xFFFFFFFFull ? std::string("broadcast") : fmt("!%08llx", static_cast<unsigned long long>(to));
    out["sender"] = fmt("!%08llx", static_cast<unsigned long long>(frm));
    out["packet_id"] = pid;
    out["hop_limit"] = fl & 7;
    out["hop_start"] = fl >> 5;
    out["want_ack"] = static_cast<bool>(fl & 8);
    out["via_mqtt"] = static_cast<bool>(fl & 16);
    out["channel_hash"] = p[13];
    return out;
  }
  return std::nullopt;
}

Json mavlink_in(const std::vector<uint8_t>& payload) {
  Json ok = Json::array();
  for (const auto& f : mavlink::scan(payload))
    if (f.contains("crc_ok") && f["crc_ok"].is_boolean() && f["crc_ok"].get<bool>()) ok.push_back(f);
  if (ok.empty()) return nullptr;
  Json o = Json::object();
  o["mavlink_frames"] = ok.size();
  const Json sm = mavlink::summary(ok);
  for (auto it = sm.begin(); it != sm.end(); ++it) o[it.key()] = it.value();
  return o;
}

std::string ascii_view(const std::vector<uint8_t>& p) {
  std::string s;
  for (uint8_t b : p) s += (b >= 32 && b < 127) ? static_cast<char>(b) : '.';
  return s;
}

}  // namespace si::rx::rc_link::lora_phy
