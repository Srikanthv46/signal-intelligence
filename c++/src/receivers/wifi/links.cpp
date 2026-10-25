// Port of python receivers/wifi/links.py (see links.hpp).
#include "receivers/wifi/links.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <string>

#include "core/npcompat.hpp"
#include "receivers/wifi/pyutil.hpp"

namespace si::rx::wifi::links {

using detail::jget;
using detail::truthy;

namespace {

int delta(int64_t a, int64_t b) {  // counter steps from seed a to seed b (1..127)
  int64_t d = ((b - a) % 127 + 127) % 127;
  return d ? static_cast<int>(d) : 127;
}

bool set_ta(Json& r, const Json& ta, const char* how) {
  if (jget(r, "ta").is_null() && truthy(ta)) {
    r["ta"] = ta;
    r["ta_evidence"] = how;
    return true;
  }
  return false;
}

// collections.Counter -> dict as json.dumps writes it (first-seen order)
void count(Json& d, const Json& key) {
  const std::string k = detail::json_key(key);
  if (d.contains(k)) d[k] = d[k].get<int64_t>() + 1;
  else d[k] = 1;
}

bool seeded(const Json& r) { return truthy(jget(r, "seed_ok")) && truthy(jget(r, "scrambler_seed")); }
int64_t seed(const Json& r) { return r["scrambler_seed"].get<int64_t>(); }

}  // namespace

Json attribute(std::vector<Json>& recs, double sifs_us) {
  std::vector<Json*> p;
  for (auto& r : recs) p.push_back(&r);
  return attribute(p, sifs_us);
}

Json attribute(const std::vector<Json*>& in, double sifs_us) {
  std::vector<Json*> recs = in;
  std::stable_sort(recs.begin(), recs.end(),
                   [](const Json* a, const Json* b) { return (*a)["t0"].get<double>() < (*b)["t0"].get<double>(); });
  const size_t n = recs.size();
  std::vector<Json> sub(n), ra(n);
  for (size_t i = 0; i < n; ++i) {
    Json& r = *recs[i];
    r.erase("ta");
    r.erase("ta_evidence");
    r.erase("ta_likely");
    const Json& m = jget(r, "mac");
    sub[i] = jget(m, "subtype");
    ra[i] = jget(m, "ra");
    if (truthy(jget(m, "ta"))) set_ta(r, m["ta"], "mac header");
  }
  auto is_resp = [](const Json& s) { return s == "ack" || s == "cts" || s == "block-ack"; };
  // ---- 2. control responses
  const double lo = -4e-6, hi = 84e-6 * sifs_us / 16.0;
  for (size_t i = 0; i + 1 < n; ++i) {
    Json& a = *recs[i];
    Json& b = *recs[i + 1];
    const double gap = b["t0"].get<double>() - a["t1"].get<double>();
    if (!(lo <= gap && gap <= hi)) continue;
    if (is_resp(sub[i + 1]) && truthy(ra[i + 1]) && !is_resp(sub[i])) set_ta(a, ra[i + 1], "response");
    if ((sub[i + 1] == "ack" || sub[i + 1] == "cts") && truthy(ra[i])) set_ta(b, ra[i], "response");
  }
  // ---- 3. scrambler counters, until nothing changes
  std::map<std::string, bool> sequential;
  auto ta_of = [&](size_t k) -> const Json& { return jget(*recs[k], "ta"); };
  for (int it = 0; it < 4; ++it) {
    bool changed = false;
    std::vector<std::pair<std::string, std::vector<size_t>>> by;
    for (size_t i = 0; i < n; ++i) {
      if (truthy(ta_of(i)) && seeded(*recs[i])) {
        const std::string t = ta_of(i).get<std::string>();
        auto f = std::find_if(by.begin(), by.end(), [&](const auto& e) { return e.first == t; });
        if (f == by.end()) by.push_back({t, {i}});
        else f->second.push_back(i);
      }
    }
    for (const auto& entry : by) {
      const std::string& ta = entry.first;
      const std::vector<size_t>& idx = entry.second;
      auto between = [&](size_t i, size_t j) {
        std::vector<size_t> out;
        for (size_t k = i + 1; k < j; ++k) {
          const Json& t = ta_of(k);
          if (t.is_null() || t == ta) out.push_back(k);
        }
        return out;
      };
      size_t ngood = 0;
      for (size_t q = 0; q + 1 < idx.size(); ++q)
        ngood += delta(seed(*recs[idx[q]]), seed(*recs[idx[q + 1]])) <=
                 static_cast<int>(between(idx[q], idx[q + 1]).size()) + 1;
      const size_t npairs = idx.size() > 0 ? idx.size() - 1 : 0;
      sequential[ta] = npairs >= 2 && static_cast<double>(ngood) / static_cast<double>(npairs) >= 0.8;
      if (!sequential[ta]) continue;
      for (size_t q = 0; q + 1 < idx.size(); ++q) {  // exact gap fill
        auto bt = between(idx[q], idx[q + 1]);
        if (!bt.empty() &&
            delta(seed(*recs[idx[q]]), seed(*recs[idx[q + 1]])) == static_cast<int>(bt.size()) + 1)
          for (size_t k : bt) changed |= set_ta(*recs[k], Json(ta), "scrambler counter");
      }
    }
    // single PPDUs whose own seed fits one sequential counter on both sides
    for (size_t k = 0; k < n; ++k) {
      Json& r = *recs[k];
      if (truthy(jget(r, "ta")) || !seeded(r)) continue;
      std::vector<std::string> fits;
      for (const auto& entry : by) {
        const std::string& ta = entry.first;
        const std::vector<size_t>& idx = entry.second;
        if (!sequential[ta]) continue;
        std::vector<size_t> before, after;
        for (size_t i : idx) (i < k ? before : after).push_back(i);
        after.erase(std::remove_if(after.begin(), after.end(), [&](size_t j) { return j <= k; }), after.end());
        const bool both = !before.empty() && !after.empty();
        bool ok = !before.empty() || !after.empty();
        auto count_open = [&](size_t a, size_t b) {
          int c = 0;
          for (size_t q = a; q < b; ++q) {
            const Json& t = ta_of(q);
            c += t.is_null() || t == ta;
          }
          return c;
        };
        if (!before.empty()) {
          const size_t i = before.back();
          const int m = count_open(i + 1, k);
          ok = ok && delta(seed(*recs[i]), seed(r)) <= (both ? std::min(m + 3, 5) : std::min(m + 1, 2));
        }
        if (!after.empty()) {
          const size_t j = after.front();
          const int m = count_open(k + 1, j);
          ok = ok && delta(seed(r), seed(*recs[j])) <= (both ? std::min(m + 3, 5) : std::min(m + 1, 2));
        }
        ok = ok && (both || idx.size() >= 4);
        if (ok) fits.push_back(ta);
      }
      if (fits.size() == 1) changed |= set_ta(r, Json(fits[0]), "scrambler counter");
    }
    if (!changed) break;
  }
  // ---- 4. per-transmitter physical summary, and a weak power hint for what is left
  std::vector<std::pair<std::string, std::vector<size_t>>> tx;
  for (size_t i = 0; i < n; ++i) {
    if (!truthy(ta_of(i))) continue;
    const std::string t = ta_of(i).get<std::string>();
    auto f = std::find_if(tx.begin(), tx.end(), [&](const auto& e) { return e.first == t; });
    if (f == tx.end()) tx.push_back({t, {i}});
    else f->second.push_back(i);
  }
  std::stable_sort(tx.begin(), tx.end(), [](const auto& a, const auto& b) { return a.second.size() > b.second.size(); });
  Json rows = Json::array();
  std::vector<std::pair<std::string, double>> levels;
  for (const auto& [ta, rs] : tx) {
    dvec pw, cfo;
    Json seeds = Json::array();
    Json vendor;
    bool have_vendor = false;
    Json evidence = Json::object(), formats = Json::object(), frame_types = Json::object();
    for (size_t i : rs) {
      const Json& r = *recs[i];
      if (!jget(r, "power_db").is_null()) pw.push_back(r["power_db"].get<double>());
      if (!jget(r, "cfo_khz").is_null()) cfo.push_back(r["cfo_khz"].get<double>());
      if (seeded(r)) seeds.push_back(r["scrambler_seed"]);
      const Json& m = jget(r, "mac");
      if (!have_vendor && truthy(m) && jget(m, "ta") == ta) {
        const Json* v = &jget(m, "ta_vendor");
        if (!truthy(*v)) v = &jget(m, "addr2_ta_vendor");
        if (!truthy(*v)) v = &jget(m, "addr2_vendor");
        vendor = *v;
        have_vendor = true;
      }
      count(evidence, r["ta_evidence"]);
      count(formats, jget(r, "fmt"));
      if (truthy(sub[i])) count(frame_types, sub[i]);
    }
    Json row = Json::object();
    row["ta"] = ta;
    row["vendor"] = vendor;
    row["n_ppdus"] = rs.size();
    row["evidence"] = evidence;
    row["formats"] = formats;
    row["frame_types"] = frame_types;
    row["power_db_median"] = pw.empty() ? Json(nullptr) : Json(py_round(np::median(pw), 1));
    if (pw.size() > 1) {
      const auto [mn, mx] = std::minmax_element(pw.begin(), pw.end());
      row["power_db_spread"] = py_round(*mx - *mn, 1);
    } else {
      row["power_db_spread"] = nullptr;
    }
    row["cfo_khz_median"] = cfo.empty() ? Json(nullptr) : Json(py_round(np::median(cfo), 2));
    row["scrambler_counter"] = sequential.count(ta) && sequential[ta] ? Json("sequential")
                               : seeds.size() >= 3                      ? Json("not sequential")
                                                                        : Json(nullptr);
    Json s64 = Json::array();
    for (size_t q = 0; q < seeds.size() && q < 64; ++q) s64.push_back(seeds[q]);
    row["seeds"] = s64;
    row["t_first_ms"] = py_round((*recs[rs.front()])["t0"].get<double>() * 1e3, 3);
    row["t_last_ms"] = py_round((*recs[rs.back()])["t0"].get<double>() * 1e3, 3);
    if (!row["power_db_median"].is_null()) levels.push_back({ta, row["power_db_median"].get<double>()});
    rows.push_back(row);
  }
  for (size_t i = 0; i < n; ++i) {
    Json& r = *recs[i];
    if (truthy(jget(r, "ta")) || jget(r, "power_db").is_null() || levels.empty()) continue;
    const double p0 = r["power_db"].get<double>();
    std::vector<std::string> near;
    for (const auto& [t, p] : levels)
      if (std::fabs(p - p0) <= 3) near.push_back(t);
    bool far = true;
    for (const auto& [t, p] : levels)
      if (std::find(near.begin(), near.end(), t) == near.end() && !(std::fabs(p - p0) > 6)) far = false;
    if (near.size() == 1 && far) r["ta_likely"] = near[0];
  }
  for (auto& row : rows) {
    int c = 0;
    for (size_t i = 0; i < n; ++i)
      if (!truthy(ta_of(i)) && jget(*recs[i], "ta_likely") == row["ta"]) ++c;
    row["n_likely_by_power"] = c;
  }
  return rows;
}

}  // namespace si::rx::wifi::links
