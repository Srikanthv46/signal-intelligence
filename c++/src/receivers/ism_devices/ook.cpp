// OOK remotes (research, signal gaps S6; python: receivers/ism_devices/ook.py): EV1527 / PT2262 fixed codes and
// KeeLoq HCS clear fields.
#include "receivers/ism_devices/ook.hpp"

#include <algorithm>

#include "core/pyfmt.hpp"

namespace si::rx::ism_devices::ook {

namespace {

double median(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  const size_t n = v.size();
  return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

struct Word {
  double t0;
  int64_t code;
  double te, t_end;
};

}  // namespace

std::optional<Json> fixed_code(const std::vector<Pulse>& p) {
  std::vector<Word> words;
  size_t i = 0;
  while (i < p.size()) {
    const double t0 = p[i].t_us, m = p[i].mark_us, s = p[i].space_us;
    if (!(100.0 <= m && m <= 1500.0 && 20 * m <= s && s <= 45 * m)) {
      ++i;
      continue;
    }
    std::vector<int> bits;
    std::vector<double> per;
    for (int k = 0; k < 24; ++k) {
      const size_t j = i + 1 + static_cast<size_t>(k);
      if (j >= p.size()) break;
      const double mk = p[j].mark_us, sp = p[j].space_us;
      const double r = mk / m;
      int b;
      if (0.4 <= r && r <= 1.8)
        b = 0;
      else if (2.2 <= r && r <= 4.2)
        b = 1;
      else
        break;
      if (k < 23) {
        if (!(2.6 * m <= mk + sp && mk + sp <= 5.4 * m)) break;
        per.push_back(mk + sp);
      }
      bits.push_back(b);
    }
    if (bits.size() == 24) {
      int64_t code = 0;
      for (int q = 0; q < 24; ++q) code |= static_cast<int64_t>(bits[q]) << (23 - q);
      if (0 < code && code < 0xFFFFFF) {
        words.push_back({t0, code, median(per) / 4, p[i + 24].t_us + p[i + 24].mark_us});
        i += 25;
        continue;
      }
    }
    ++i;
  }
  if (words.empty()) return std::nullopt;
  int64_t best = words[0].code;
  size_t best_n = 0;
  for (const Word& w : words) {    // most repeated; ties: the first seen
    const size_t n = static_cast<size_t>(
        std::count_if(words.begin(), words.end(), [&](const Word& x) { return x.code == w.code; }));
    if (n > best_n) {
      best = w.code;
      best_n = n;
    }
  }
  std::vector<Word> sel;
  for (const Word& w : words)
    if (w.code == best) sel.push_back(w);
  std::vector<double> tes;
  for (const Word& w : sel) tes.push_back(w.te);
  Json out = Json::object();
  out["code"] = fmt("%06llX", static_cast<long long>(best));
  out["id"] = fmt("%05llX", static_cast<long long>(best >> 4));
  out["data"] = best & 0xF;
  out["repeats"] = static_cast<int64_t>(sel.size());
  out["te_us"] = py_round(median(tes), 1);
  out["t_first_us"] = py_round(sel.front().t0, 1);
  out["t_last_us"] = py_round(sel.back().t0, 1);
  out["t_end_us"] = py_round(sel.back().t_end, 1);
  out["words"] = static_cast<int64_t>(words.size());
  std::string tri;
  for (int q = 0; q < 12; ++q) {
    const int a = (best >> (23 - 2 * q)) & 1, b = (best >> (22 - 2 * q)) & 1;
    if (a == 0 && b == 0)
      tri += '0';
    else if (a == 1 && b == 1)
      tri += '1';
    else if (a == 0 && b == 1)
      tri += 'F';
    else {
      tri.clear();
      break;
    }
  }
  if (tri.size() == 12) out["tristate"] = tri;
  return out;
}

std::optional<Json> keeloq(const std::vector<Pulse>& p) {
  std::vector<Json> pkts;
  size_t i = 0;
  while (i + 12 + 66 <= p.size()) {
    const double te = p[i].mark_us;
    bool ok = 200.0 <= te && te <= 800.0;
    for (int k = 0; k < 12 && ok; ++k) {
      const double mk = p[i + k].mark_us, sp = p[i + k].space_us;
      ok = 0.6 * te <= mk && mk <= 1.4 * te &&
           (k < 11 ? (0.6 * te <= sp && sp <= 1.4 * te) : (7 * te <= sp && sp <= 13 * te));
    }
    if (!ok) {
      ++i;
      continue;
    }
    std::vector<int> bits;
    for (int k = 0; k < 66; ++k) {
      const double mk = p[i + 12 + k].mark_us, sp = p[i + 12 + k].space_us;
      const double r = mk / te;
      if (0.5 <= r && r <= 1.5)
        bits.push_back(1);
      else if (1.5 < r && r <= 2.6)
        bits.push_back(0);
      else
        break;
      if (k < 65 && !(2.2 * te <= mk + sp && mk + sp <= 3.8 * te)) break;
    }
    const bool all_ones = bits.size() == 66 && std::all_of(bits.begin(), bits.begin() + 64, [](int b) { return b == 1; });
    if (bits.size() == 66 && !all_ones) {
      int64_t serial = 0;
      for (int q = 0; q < 28; ++q) serial |= static_cast<int64_t>(bits[32 + q]) << q;
      const int btn = (bits[60] << 3) | (bits[61] << 2) | (bits[62] << 1) | bits[63];
      Json k = Json::object();
      k["t_us"] = py_round(p[i].t_us, 1);
      k["serial"] = fmt("%07llX", static_cast<long long>(serial));
      k["button"] = static_cast<int64_t>((btn & 0x08) | ((btn & 0x01) << 2) | (btn & 0x02) | ((btn & 0x04) >> 2));
      k["learn"] = btn == 0x0F;
      k["battery_low"] = bits[64] == 1;
      k["repeat"] = bits[65] == 1;
      k["te_us"] = py_round(te, 1);
      pkts.push_back(k);
      i += 78;
      continue;
    }
    ++i;
  }
  if (pkts.empty()) return std::nullopt;
  const Json& first = pkts[0];
  int64_t same = 0;
  for (const Json& q : pkts) same += q["serial"] == first["serial"] && q["button"] == first["button"];
  Json out = Json::object();
  for (const char* k : {"serial", "button", "learn", "battery_low", "repeat", "te_us"}) out[k] = first[k];
  out["packets"] = static_cast<int64_t>(pkts.size());
  out["matching_packets"] = same;
  out["t_first_us"] = first["t_us"];
  return out;
}

}  // namespace si::rx::ism_devices::ook
