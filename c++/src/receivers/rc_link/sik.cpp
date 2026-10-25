// python receivers/rc_link/sik.py (same variants, order and keys).
#include "receivers/rc_link/sik.hpp"

#include <array>
#include <cmath>
#include <cstdlib>
#include <string>

#include "core/mavlink.hpp"
#include "core/npcompat.hpp"
#include "core/pyfmt.hpp"

namespace si::rx::rc_link::sik {

const std::vector<double> kRates = {2e3, 4e3, 8e3, 16e3, 19e3, 24e3, 32e3, 48e3, 64e3, 96e3, 128e3, 192e3, 250e3};

namespace {
constexpr int kSync = 0x2DD4;
constexpr size_t kWBytes = 260;

std::vector<uint8_t> pn9_bits(size_t n) {
  std::vector<uint8_t> out(n);
  unsigned reg = 0x1FF;
  for (size_t i = 0; i < n; ++i) {
    out[i] = reg & 1;
    const unsigned fb = (reg ^ (reg >> 5)) & 1;
    reg = (reg >> 1) | (fb << 8);
  }
  return out;
}

std::vector<uint8_t> whiten_bytes(bool msb_first) {
  const auto w = pn9_bits(8 * kWBytes);
  std::vector<uint8_t> out(kWBytes);
  for (size_t b = 0; b < kWBytes; ++b) {
    unsigned v = 0;
    for (int i = 0; i < 8; ++i) v |= w[8 * b + i] << (msb_first ? 7 - i : i);
    out[b] = static_cast<uint8_t>(v);
  }
  return out;
}

const std::array<std::pair<const char*, std::vector<uint8_t>>, 2>& whiten() {
  static const std::array<std::pair<const char*, std::vector<uint8_t>>, 2> w = {
      std::make_pair("pn9_lsb_first", whiten_bytes(false)), std::make_pair("pn9_msb_first", whiten_bytes(true))};
  return w;
}

// python crc16_variants: [(name, value)] in python's order
std::vector<std::pair<std::string, int>> crc16_variants(const std::vector<uint8_t>& d) {
  std::vector<std::pair<std::string, int>> out;
  for (unsigned init : {0x0000u, 0xFFFFu}) {
    unsigned c = init;
    for (uint8_t b : d) {
      c ^= static_cast<unsigned>(b) << 8;
      for (int k = 0; k < 8; ++k) c = (c & 0x8000) ? ((c << 1) ^ 0x8005) & 0xFFFF : (c << 1) & 0xFFFF;
    }
    out.emplace_back(fmt("ibm_msb_init%04x", init), static_cast<int>(c));
    c = init;
    for (uint8_t b : d) {
      c ^= b;
      for (int k = 0; k < 8; ++k) c = (c & 1) ? (c >> 1) ^ 0xA001 : c >> 1;
    }
    out.emplace_back(fmt("ibm_refl_init%04x", init), static_cast<int>(c));
  }
  return out;
}

std::string hexs(const uint8_t* p, size_t n) {
  static const char* d = "0123456789abcdef";
  std::string s;
  for (size_t i = 0; i < n; ++i) s += d[p[i] >> 4], s += d[p[i] & 15];
  return s;
}
}  // namespace

std::vector<std::pair<int64_t, bool>> find_sync(const std::vector<uint8_t>& bits) {
  std::vector<std::pair<int64_t, bool>> hits;
  const int64_t n = static_cast<int64_t>(bits.size());
  if (n < 40) return hits;
  int pat[16], pre[8];
  for (int i = 0; i < 16; ++i) pat[i] = ((kSync >> (15 - i)) & 1) * 2 - 1;
  for (int i = 0; i < 8; ++i) pre[i] = ((i + 1) % 2) * 2 - 1;
  std::vector<int> c(static_cast<size_t>(n - 15));
  for (int64_t p = 0; p + 16 <= n; ++p) {
    int v = 0;
    for (int i = 0; i < 16; ++i) v += (bits[p + i] * 2 - 1) * pat[i];
    c[p] = v;
  }
  for (bool inv : {false, true})
    for (int64_t p = 0; p < static_cast<int64_t>(c.size()); ++p) {
      if ((inv ? -c[p] : c[p]) < 14) continue;
      if (p < 8) continue;
      int q = 0;
      for (int i = 0; i < 8; ++i) q += (bits[p - 8 + i] * 2 - 1) * (inv ? -1 : 1) * pre[i];
      if (std::abs(q) >= 6) hits.emplace_back(p + 16, inv);
    }
  return hits;
}

Json decode_packet(const std::vector<uint8_t>& bits, bool inverted) {
  const size_t nby = bits.size() / 8;
  if (nby < 6) return nullptr;
  std::vector<uint8_t> raw(nby);
  for (size_t b = 0; b < nby; ++b) {
    unsigned v = 0;
    for (int i = 0; i < 8; ++i) v = (v << 1) | ((bits[8 * b + i] ^ (inverted ? 1 : 0)) & 1);
    raw[b] = static_cast<uint8_t>(v);
  }
  Json best = nullptr;
  int best_score = -1;
  for (const auto& [wname, w] : whiten()) {
    std::vector<uint8_t> d(nby);
    for (size_t i = 0; i < nby; ++i) d[i] = raw[i] ^ (i < w.size() ? w[i] : 0);
    const size_t L = d[2];
    if (3 + L + 2 > nby || L < 3) continue;
    const std::vector<uint8_t> hld(d.begin(), d.begin() + static_cast<ptrdiff_t>(3 + L));
    const int be = (d[3 + L] << 8) | d[4 + L], le = (d[4 + L] << 8) | d[3 + L];
    Json match = Json::array();
    for (const auto& [k, v] : crc16_variants(hld)) {
      if (v == be) match.push_back(k + "_hi_first");
      if (v == le) match.push_back(k + "_lo_first");
    }
    const std::vector<uint8_t> data(d.begin() + 3, d.begin() + static_cast<ptrdiff_t>(3 + L));
    const std::vector<uint8_t> mav = L > 2 ? std::vector<uint8_t>(data.begin(), data.end() - 2) : data;
    const Json frames = mavlink::scan(mav);
    int n_ok = 0;
    for (const auto& f : frames)
      if (f["crc_ok"].is_boolean() && f["crc_ok"].get<bool>()) ++n_ok;
    const int score = 2 * n_ok + (match.empty() ? 0 : 3);
    if (score == 0) continue;
    Json cand = {{"netid", (d[0] << 8) | d[1]},
                 {"length", L},
                 {"whitening", wname},
                 {"hw_crc_variant", match.empty() ? Json(nullptr) : match},
                 {"hw_crc_ok", !match.empty()},
                 {"mavlink_frames", frames},
                 {"mavlink_crc_ok", n_ok},
                 {"tdm_trailer_hex", L >= 2 ? Json(hexs(d.data() + 3 + L - 2, 2)) : Json(nullptr)},
                 {"data_hex", hexs(data.data(), data.size())}};
    if (score > best_score) best = cand, best_score = score;
  }
  return best;
}

double two_level_threshold(const std::vector<double>& v, int iters) {
  double thr = 0.5 * (np::percentile(v, 5.0) + np::percentile(v, 95.0));
  std::vector<double> lo, hi;
  for (int it = 0; it < iters; ++it) {
    lo.clear();
    hi.clear();
    for (double x : v) (x < thr ? lo : hi).push_back(x);
    if (lo.empty() || hi.empty()) break;
    thr = 0.5 * (np::mean(lo) + np::mean(hi));
  }
  return thr;
}

namespace {
std::vector<std::pair<double, std::vector<uint8_t>>> demod_impl(const cfvec& y, double ofs, double rate, bool mid);
}

std::vector<std::pair<double, std::vector<uint8_t>>> demod_bits_mid(const cfvec& y, double ofs, double rate) {
  return demod_impl(y, ofs, rate, true);
}

std::vector<std::pair<double, std::vector<uint8_t>>> demod_bits(const cfvec& y, double ofs, double rate) {
  return demod_impl(y, ofs, rate, false);
}

namespace {
std::vector<std::pair<double, std::vector<uint8_t>>> demod_impl(const cfvec& y, double ofs, double rate, bool mid) {
  std::vector<std::pair<double, std::vector<uint8_t>>> out;
  if (y.size() < 3) return out;
  const double sps = ofs / rate;
  std::vector<double> ifr(y.size() - 1);
  for (size_t i = 0; i + 1 < y.size(); ++i) {
    const cf z = y[i + 1] * std::conj(y[i]);
    ifr[i] = static_cast<double>(std::atan2(z.imag(), z.real()));
  }
  const int64_t k = std::max<int64_t>(1, static_cast<int64_t>(sps / 2));
  const int64_t n = static_cast<int64_t>(ifr.size()), off = (k - 1) / 2;
  std::vector<double> ifs(static_cast<size_t>(n), 0.0);
  for (int64_t i = 0; i < n; ++i) {
    double s = 0;
    for (int64_t j = 0; j < k; ++j) {
      const int64_t m = i + off - j;
      if (m >= 0 && m < n) s += ifr[m] * (1.0 / static_cast<double>(k));
    }
    ifs[i] = s;
  }
  const double thr = mid ? two_level_threshold(ifs) : np::median(ifs);
  const double step = std::max(sps / 8, 1.0);
  for (int64_t pi = 0; static_cast<double>(pi) * step < sps; ++pi) {
    const double ph = static_cast<double>(pi) * step;
    const int64_t nb = static_cast<int64_t>((static_cast<double>(n) - ph) / sps);
    std::vector<uint8_t> bits(static_cast<size_t>(std::max<int64_t>(nb, 0)));
    for (int64_t i = 0; i < nb; ++i) bits[i] = ifs[static_cast<int64_t>(ph + static_cast<double>(i) * sps)] > thr;
    out.emplace_back(ph, std::move(bits));
  }
  return out;
}
}  // namespace

const std::vector<uint8_t>& whiten_pn9_msb() { return whiten()[1].second; }

Json sik_try_decode(const cfvec& y, double ofs, double rate) {
  if (y.size() < 3) return nullptr;
  const double sps = ofs / rate;
  std::vector<double> ifr(y.size() - 1);
  for (size_t i = 0; i + 1 < y.size(); ++i) {
    const cf z = y[i + 1] * std::conj(y[i]);
    ifr[i] = static_cast<double>(std::atan2(z.imag(), z.real()));
  }
  const int64_t k = std::max<int64_t>(1, static_cast<int64_t>(sps / 2));
  const int64_t n = static_cast<int64_t>(ifr.size()), off = (k - 1) / 2;
  std::vector<double> ifs(static_cast<size_t>(n), 0.0);  // np.convolve(ifr, ones(k)/k, "same")
  for (int64_t i = 0; i < n; ++i) {
    double s = 0;
    for (int64_t j = 0; j < k; ++j) {
      const int64_t m = i + off - j;
      if (m >= 0 && m < n) s += ifr[m] * (1.0 / static_cast<double>(k));
    }
    ifs[i] = s;
  }
  const double thr = np::median(ifs);
  Json phy = nullptr;
  const double step = std::max(sps / 8, 1.0);
  for (int64_t pi = 0; static_cast<double>(pi) * step < sps; ++pi) {  // np.arange(0, sps, step)
    const double ph = static_cast<double>(pi) * step;
    const int64_t nb = static_cast<int64_t>((static_cast<double>(n) - ph) / sps);
    std::vector<uint8_t> bits(static_cast<size_t>(std::max<int64_t>(nb, 0)));
    for (int64_t i = 0; i < nb; ++i) bits[i] = ifs[static_cast<int64_t>(ph + static_cast<double>(i) * sps)] > thr;
    for (const auto& [pos, inv] : find_sync(bits)) {
      Json r = decode_packet(std::vector<uint8_t>(bits.begin() + pos, bits.end()), inv);
      if (r.is_object()) {
        r["sync"] = "0x2DD4";
        r["polarity_inverted"] = inv;
        r["air_rate_bps"] = rate;
        return r;
      }
      if (phy.is_null()) phy = {{"sync", "0x2DD4"}, {"polarity_inverted", inv}, {"air_rate_bps", rate}, {"decoded", false}};
    }
  }
  return phy;
}

}  // namespace si::rx::rc_link::sik
