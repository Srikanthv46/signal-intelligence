// python receivers/broadcast_fm/rds.py (see there for the code and the sources).
#include "receivers/broadcast_fm/rds.hpp"

#include <algorithm>
#include <map>
#include <set>

#include "core/pyfmt.hpp"

namespace si::rx::broadcast_fm::rds {

namespace {

constexpr uint32_t kG = 0x5B9;
const std::map<std::string, int>& offsets() {
  static const std::map<std::string, int> m = {{"A", 0x0FC}, {"B", 0x198}, {"C", 0x168}, {"Cp", 0x350}, {"D", 0x1B4}};
  return m;
}
std::string by_syndrome(int s) {
  for (const auto& [k, v] : offsets())
    if (v == s) return k;
  return "";
}
bool follows(const std::string& cur, const std::string& o) {
  if (cur == "A") return o == "B";
  if (cur == "B") return o == "C" || o == "Cp";
  if (cur == "C" || cur == "Cp") return o == "D";
  if (cur == "D") return o == "A";
  return false;
}
const char* kPty[32] = {"none", "News", "Current affairs", "Information", "Sport", "Education", "Drama", "Culture",
                        "Science", "Varied", "Pop music", "Rock music", "Easy listening", "Light classical",
                        "Serious classical", "Other music", "Weather", "Finance", "Children's programmes",
                        "Social affairs", "Religion", "Phone-in", "Travel", "Leisure", "Jazz music", "Country music",
                        "National music", "Oldies music", "Folk music", "Documentary", "Alarm test", "Alarm"};

uint32_t mod26(uint32_t v, int nbits) {
  for (int i = nbits - 1; i > 9; --i)
    if ((v >> i) & 1) v ^= kG << (i - 10);
  return v;
}

std::string chars(int v) {
  std::string s;
  for (int c : {(v >> 8) & 0xFF, v & 0xFF}) s += (32 <= c && c < 127) ? static_cast<char>(c) : (c == 13 ? '\r' : '.');
  return s;
}

}  // namespace

int syndrome(uint32_t w) { return static_cast<int>(mod26(w, 26)); }

uint32_t encode_block(uint16_t d, const std::string& off) {
  return (static_cast<uint32_t>(d) << 10) | (mod26(static_cast<uint32_t>(d) << 10, 26) ^ offsets().at(off));
}

std::vector<int> encode_group(uint16_t a, uint16_t b, uint16_t c, uint16_t d, bool vb) {
  std::vector<int> out;
  const uint16_t data[4] = {a, b, c, d};
  const std::string offs[4] = {"A", "B", vb ? "Cp" : "C", "D"};
  for (int k = 0; k < 4; ++k) {
    const uint32_t w = encode_block(data[k], offs[k]);
    for (int i = 0; i < 26; ++i) out.push_back((w >> (25 - i)) & 1);
  }
  return out;
}

std::vector<Block> blocks(const std::vector<int>& bits) {
  const int64_t n = static_cast<int64_t>(bits.size());
  std::vector<uint32_t> words;
  for (int64_t i = 0; i + 25 < n; ++i) {
    uint32_t w = 0;
    for (int64_t k = i; k < i + 26; ++k) w = (w << 1) | static_cast<uint32_t>(bits[k]);
    words.push_back(w);
  }
  const int64_t nw = static_cast<int64_t>(words.size());
  std::vector<std::string> syn(nw);
  for (int64_t i = 0; i < nw; ++i) syn[i] = by_syndrome(syndrome(words[i]));
  std::vector<Block> out;
  std::set<int64_t> used;
  int64_t i = 0;
  while (i + 26 < nw) {
    const std::string& a = syn[i];
    const std::string& b = syn[i + 26];
    if (!a.empty() && !b.empty() && follows(a, b) && !used.count(i)) {
      int64_t j = i;
      std::string cur = a;
      while (j < nw) {
        const std::string& o = syn[j];
        if (o.empty() || (j != i && !follows(cur, o))) break;
        out.push_back({j, o, static_cast<int>(words[j] >> 10)});
        used.insert(j);
        cur = o;
        j += 26;
      }
      i = j;
      continue;
    }
    ++i;
  }
  return out;
}

std::vector<Group> groups(const std::vector<Block>& blks) {
  std::vector<Group> out;
  for (size_t k = 0; k < blks.size(); ++k) {
    if (blks[k].offset != "A") continue;
    Group g;
    g.pos = blks[k].pos;
    g.a = blks[k].data;
    for (size_t q = k + 1; q < blks.size() && q < k + 4; ++q) {
      const int64_t dp = blks[q].pos - g.pos;
      const int64_t step = dp >= 0 ? dp / 26 : -((-dp + 25) / 26);   // python floor division
      if (dp != 26 * step || step > 3) break;
      const std::string& o = blks[q].offset;
      if (o == "B") g.b = blks[q].data;
      else if (o == "C" || o == "Cp") g.c = blks[q].data;
      else if (o == "D") g.d = blks[q].data;
      if (o == "Cp") g.version_b = true;
    }
    out.push_back(g);
  }
  return out;
}

Json decode(const std::vector<Group>& grps) {
  Json out = Json::object();
  if (grps.empty()) return out;
  std::vector<int> pis;
  for (const auto& g : grps) pis.push_back(g.a);
  int pi = pis[0];
  int64_t best_n = -1, best_i = 0;
  for (size_t i = 0; i < pis.size(); ++i) {
    const int p = pis[i];
    const int64_t n = std::count(pis.begin(), pis.end(), p);
    const int64_t first = std::find(pis.begin(), pis.end(), p) - pis.begin();
    if (n > best_n || (n == best_n && first < best_i)) best_n = n, best_i = first, pi = p;
  }
  out["pi"] = fmt("%04X", pi);
  out["pi_country_nibble"] = pi >> 12;
  std::string ps(8, '_'), rt(64, '_');
  std::set<std::string> types;
  for (const auto& g : grps) {
    if (g.a != pi || !g.b) continue;
    const int b = *g.b;
    const int gt = b >> 12;
    const bool vb = (b >> 11) & 1;
    types.insert(std::to_string(gt) + (vb ? "B" : "A"));
    out["pty"] = (b >> 5) & 31;
    out["pty_name"] = kPty[(b >> 5) & 31];
    out["tp"] = static_cast<bool>((b >> 10) & 1);
    if (gt == 0 && g.d) {
      const int seg = b & 3;
      const std::string c = chars(*g.d);
      ps[2 * seg] = c[0];
      ps[2 * seg + 1] = c[1];
      out["ta"] = static_cast<bool>((b >> 4) & 1);
      out["music"] = static_cast<bool>((b >> 3) & 1);
    } else if (gt == 2) {
      const int seg = b & 15;
      if (!vb && g.c && g.d) {
        const std::string c = chars(*g.c) + chars(*g.d);
        for (int k = 0; k < 4; ++k) rt[4 * seg + k] = c[k];
      } else if (vb && g.d) {
        const std::string c = chars(*g.d);
        for (int k = 0; k < 2; ++k) rt[2 * seg + k] = c[k];
      }
    } else if (gt == 4 && !vb && g.c && g.d) {
      const int64_t mjd = (static_cast<int64_t>(b & 3) << 15) | (*g.c >> 1);
      const int hour = ((*g.c & 1) << 4) | (*g.d >> 12);
      const int minute = (*g.d >> 6) & 63;
      const int off = (*g.d & 31) * (((*g.d >> 5) & 1) ? -1 : 1);
      const int64_t y1 = static_cast<int64_t>((static_cast<double>(mjd) - 15078.2) / 365.25);
      const int64_t m1 = static_cast<int64_t>((static_cast<double>(mjd) - 14956.1 -
                                               static_cast<double>(static_cast<int64_t>(static_cast<double>(y1) * 365.25))) /
                                              30.6001);
      const int64_t day = mjd - 14956 - static_cast<int64_t>(static_cast<double>(y1) * 365.25) -
                          static_cast<int64_t>(static_cast<double>(m1) * 30.6001);
      const int64_t k = (m1 == 14 || m1 == 15) ? 1 : 0;
      Json ct = Json::object();
      ct["utc"] = fmt("%04lld-%02lld-%02lldT%02d:%02dZ", static_cast<long long>(1900 + y1 + k),
                      static_cast<long long>(m1 - 1 - 12 * k), static_cast<long long>(day), hour, minute);
      ct["local_offset_h"] = off / 2.0;
      out["clock_time"] = ct;
    }
  }
  if (ps != std::string(8, '_') && ps.find_first_not_of('_') != std::string::npos) out["ps"] = ps;
  std::string txt = rt.substr(0, rt.find('\r'));
  while (!txt.empty() && txt.back() == '_') txt.pop_back();
  if (txt.find_first_not_of('_') != std::string::npos) out["radiotext"] = txt;
  out["group_types"] = Json::array();
  for (const auto& t : types) out["group_types"].push_back(t);
  return out;
}

}  // namespace si::rx::broadcast_fm::rds
