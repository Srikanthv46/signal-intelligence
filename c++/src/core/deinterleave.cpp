// Pulse deinterleaving (research, signal gaps S10; python: core/deinterleave.py). Same order of operations as the
// Python lists (sequential sums, stable sorts, value-based removal of extracted TOAs).
#include "core/deinterleave.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <optional>
#include <set>
#include <string>

namespace si::deint {

namespace {

constexpr double kMinD = 20e-6, kMaxD = 20e-3;
constexpr double kTolAbs = 1.5e-6, kTolRel = 0.02;
constexpr double kPwAbs = 1.5e-6, kPwRel = 0.30;
constexpr int kMaxLevel = 8;
constexpr double kCdifX = 0.25;
constexpr int kMaxMiss = 2;
constexpr size_t kMinChain = 4;
constexpr double kMinFill = 0.6;
constexpr size_t kMaxEmitters = 4;
constexpr size_t kMaxCand = 12;
constexpr double kPChance = 1e-4;
constexpr int64_t kMinReport = 6;
constexpr double kChainCost = 6.0;
constexpr double kTolJit = 0.10;
constexpr double kJitGain = 1.25;

double median(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  const size_t n = v.size();
  return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

struct Cand {
  double pri;
  int level;
};

std::vector<Cand> cdif(const std::vector<double>& r) {
  const size_t e = r.size();
  std::vector<Cand> out;
  std::vector<double> acc;
  const int top = static_cast<int>(std::min<size_t>(kMaxLevel, e - 1));
  for (int c = 1; c <= top; ++c) {
    for (size_t i = 0; i + static_cast<size_t>(c) < e; ++i) {
      const double v = r[i + c] - r[i];
      if (kMinD <= v && v <= kMaxD) acc.push_back(v);
    }
    std::vector<double> d = acc;
    std::sort(d.begin(), d.end());
    const double need = std::max(3.0, kCdifX * static_cast<double>(e - c));
    size_t i = 0;
    while (i < d.size()) {
      size_t j = i + 1;
      while (j < d.size() && d[j] - d[j - 1] <= std::max(kTolAbs, kTolRel * d[j - 1])) ++j;
      if (static_cast<double>(j - i) >= need) {
        const double m = median(std::vector<double>(d.begin() + static_cast<ptrdiff_t>(i), d.begin() + static_cast<ptrdiff_t>(j)));
        bool fresh = true;
        for (const Cand& q : out) fresh = fresh && std::fabs(m - q.pri) > std::max(kTolAbs, kTolRel * q.pri);
        if (fresh) out.push_back({m, c});
      }
      i = j;
    }
    if (out.size() >= kMaxCand) break;
  }
  if (out.size() > kMaxCand) out.resize(kMaxCand);
  return out;
}

// python chain(r, pri, tol, rate, n_all) -> (indices, missed slots)
std::pair<std::vector<size_t>, int> chain(const std::vector<double>& r, double pri, double tol, double rate, size_t n_all) {
  std::vector<size_t> best;
  std::vector<int> best_k;
  std::vector<bool> seen(r.size(), false);
  for (size_t s = 0; s < r.size(); ++s) {
    if (seen[s]) continue;
    std::vector<size_t> cur{s};
    std::vector<int> ks;
    size_t j = s;
    while (true) {
      std::optional<size_t> nxt;
      int kk = 0;
      for (int k = 1; k <= kMaxMiss + 1; ++k) {
        const double want = r[cur.back()] + k * pri;
        size_t q = j + 1;
        while (q < r.size() && r[q] < want - tol) ++q;
        if (q < r.size() && std::fabs(r[q] - want) <= tol) {
          nxt = q;
          kk = k;
          break;
        }
      }
      if (!nxt) break;
      cur.push_back(*nxt);
      ks.push_back(kk);
      seen[*nxt] = true;
      j = *nxt;
    }
    if (cur.size() > best.size()) {
      best = cur;
      best_k = ks;
    }
  }
  if (best.size() < kMinChain) return {{}, 0};
  const double span = r[best.back()] - r[best.front()];
  if (static_cast<double>(best.size()) < kMinFill * (span / pri + 1)) return {{}, 0};
  double s = 0.0;
  for (int k : best_k) s += std::log(std::min(1.0, std::max(k * rate * 2 * tol, 1e-300)));
  const double logp = std::log(static_cast<double>(n_all * kMaxLevel)) + s;
  if (logp >= std::log(kPChance)) return {{}, 0};
  int miss = 0;
  for (int k : best_k) miss += k - 1;
  return {best, miss};
}

struct Chain {
  std::vector<double> t, w;
  int miss;
};

struct Try {
  double score;
  std::vector<Chain> used;
  double raw;
};

std::optional<Try> try_pri(const std::vector<double>& r, const std::vector<double>& w, double pri, double tol, bool single) {
  const double rate = r.back() > r.front() ? static_cast<double>(r.size() - 1) / (r.back() - r.front()) : 0.0;
  std::vector<Chain> chains;
  std::vector<double> rt = r, rw = w;
  while (rt.size() >= kMinChain && !(single && !chains.empty())) {
    const auto [c, miss] = chain(rt, pri, tol, rate, r.size());
    if (c.empty()) break;
    Chain ch{{}, {}, miss};
    std::vector<bool> drop(rt.size(), false);
    for (size_t i : c) {
      ch.t.push_back(rt[i]);
      ch.w.push_back(rw[i]);
      drop[i] = true;
    }
    chains.push_back(ch);
    std::vector<double> nt, nw;
    for (size_t i = 0; i < rt.size(); ++i)
      if (!drop[i]) {
        nt.push_back(rt[i]);
        nw.push_back(rw[i]);
      }
    rt = nt;
    rw = nw;
  }
  if (chains.empty()) return std::nullopt;
  const double lo = chains[0].t.front(), hi = chains[0].t.back();
  Try t{0.0, {}, 0.0};
  for (const Chain& ch : chains)
    if (ch.t.front() < hi && ch.t.back() > lo) t.used.push_back(ch);
  double raw = 0;
  for (const Chain& ch : t.used) raw += static_cast<double>(static_cast<int64_t>(ch.t.size()) - ch.miss);
  t.raw = raw;
  t.score = raw - kChainCost * static_cast<double>(t.used.size() - 1);
  return t;
}

// python _same_phase(a, b, pri): b continues the train through a on the PRI grid (a beam-scan gap broke the chain)
bool same_phase(double a, double b, double pri) {
  const double n = static_cast<double>(py_round_int((b - a) / pri));
  return std::fabs(b - a - n * pri) <= std::max(kTolAbs, kTolRel * pri);
}

Json describe(double pri, const std::vector<Chain>& used) {
  std::vector<double> ts, ws;
  for (const Chain& ch : used) {
    ts.insert(ts.end(), ch.t.begin(), ch.t.end());
    ws.insert(ws.end(), ch.w.begin(), ch.w.end());
  }
  std::sort(ts.begin(), ts.end());
  double lo = used[0].t.front();
  for (const Chain& ch : used) lo = std::min(lo, ch.t.front());
  const double tol = std::max(kTolAbs, kTolRel * pri / static_cast<double>(used.size()));
  std::vector<double> raw;
  for (const Chain& ch : used) raw.push_back(std::fmod(ch.t.front() - lo, pri));
  std::sort(raw.begin(), raw.end());
  std::vector<double> ph;   // distinct phases: a continuation chain repeats one
  for (double x : raw)
    if (ph.empty() || x - ph.back() > tol) ph.push_back(x);
  if (ph.size() >= 2 && ph[0] + pri - ph.back() <= tol) ph.pop_back();   // the circular wrap
  const size_t K = ph.size();
  std::vector<double> sub;
  for (size_t i = 0; i < K; ++i) sub.push_back((i + 1 < K ? ph[i + 1] : ph[0] + pri) - ph[i]);
  Json stagger = nullptr;
  const double smax = *std::max_element(sub.begin(), sub.end()), smin = *std::min_element(sub.begin(), sub.end());
  if (K >= 2 && smax - smin > 2 * tol) {
    stagger = Json::array();
    for (double x : sub) stagger.push_back(py_round(x * 1e6, 2));
  } else if (K >= 2) {
    pri = pri / static_cast<double>(K);   // equal steps: a constant PRI seen as a frame of K chains
  }
  const double q = stagger.is_null() && K >= 2 ? pri * static_cast<double>(K) : pri;
  std::vector<double> steps;
  for (const Chain& ch : used)
    for (size_t i = 0; i + 1 < ch.t.size(); ++i) {
      const double d = ch.t[i + 1] - ch.t[i];
      steps.push_back(d / static_cast<double>(py_round_int(d / q)));
    }
  double jit = 0.0;
  if (!steps.empty()) {
    double s = 0.0;
    for (double x : steps) s += std::pow(x - q, 2.0);
    jit = std::sqrt(s / static_cast<double>(steps.size())) / q;
  }
  Json e = Json::object();
  e["pri_us"] = py_round(pri * 1e6, 2);
  e["n_pulses"] = static_cast<int64_t>(ts.size());
  e["pw_us"] = py_round(median(ws) * 1e6, 2);
  e["stagger_us"] = stagger;
  e["jitter"] = py_round(jit, 4);
  e["t_first_s"] = py_round(ts.front(), 7);
  return e;
}

struct Found {
  double pri;
  std::vector<Chain> used;
  bool strict;
};

std::pair<std::vector<Json>, std::vector<double>> emitters(std::vector<double> r, std::vector<double> w) {
  std::vector<Found> found;
  while (r.size() >= kMinChain && found.size() < kMaxEmitters) {
    std::optional<Try> best, jit;
    double best_pri = 0, jit_pri = 0;
    const auto cand = cdif(r);
    for (const Cand& c : cand) {   // strict: a frame period sums `level` intervals
      auto t = try_pri(r, w, c.pri, std::max(kTolAbs, kTolRel * c.pri / std::sqrt(static_cast<double>(c.level))), false);
      if (t && (!best || t->score > best->score)) {
        best = t;
        best_pri = c.pri;
      }
    }
    for (const Cand& c : cand) {   // jitter mode: one chain at 10 %
      auto t = try_pri(r, w, c.pri, std::max(kTolAbs, kTolJit * c.pri), true);
      if (t && (!jit || t->score > jit->score)) {
        jit = t;
        jit_pri = c.pri;
      }
    }
    bool strict = true;
    if (jit && (!best || jit->raw > kJitGain * best->raw)) {
      best = jit;
      best_pri = jit_pri;
      strict = false;
    }
    if (!best) break;
    const double pri = best_pri;
    const auto& used = best->used;
    double ulo = 1e300, uhi = -1e300;
    for (const Chain& ch : used) {
      ulo = std::min(ulo, ch.t.front());
      uhi = std::max(uhi, ch.t.back());
    }
    Found* prev = nullptr;
    for (Found& f : found) {
      if (!(strict && f.strict && std::fabs(f.pri - pri) <= std::max(kTolAbs, kTolRel * pri))) continue;
      double flo = 1e300, fhi = -1e300;
      for (const Chain& ch : f.used) {
        flo = std::min(flo, ch.t.front());
        fhi = std::max(fhi, ch.t.back());
      }
      if ((flo < uhi && fhi > ulo) || same_phase(f.used[0].t.front(), used[0].t.front(), pri)) {
        prev = &f;
        break;
      }
    }
    if (prev) {   // more phases of a stagger, or the same train after a gap
      prev->used.insert(prev->used.end(), used.begin(), used.end());
    } else {
      found.push_back({pri, used, strict});
    }
    std::set<double> gone;
    for (const Chain& ch : used) gone.insert(ch.t.begin(), ch.t.end());
    std::vector<double> nr, nw;
    for (size_t i = 0; i < r.size(); ++i)
      if (!gone.count(r[i])) {
        nr.push_back(r[i]);
        nw.push_back(w[i]);
      }
    r = nr;
    w = nw;
  }
  std::vector<Json> out;
  for (const Found& f : found) out.push_back(describe(f.pri, f.used));
  return {out, r};
}

}  // namespace

Json deinterleave(const std::vector<double>& toa_s, const std::vector<double>& pw_s) {
  std::vector<size_t> order(toa_s.size());
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    return pw_s[a] < pw_s[b] || (pw_s[a] == pw_s[b] && toa_s[a] < toa_s[b]);
  });
  std::vector<std::vector<size_t>> groups;
  for (size_t i : order) {
    if (!groups.empty() && pw_s[i] <= pw_s[groups.back()[0]] + std::max(kPwAbs, kPwRel * pw_s[groups.back()[0]]))
      groups.back().push_back(i);
    else
      groups.push_back({i});
  }
  std::stable_sort(groups.begin(), groups.end(), [](const auto& a, const auto& b) { return a.size() > b.size(); });
  std::vector<Json> ems;
  int64_t left = 0;
  for (auto g : groups) {
    std::stable_sort(g.begin(), g.end(), [&](size_t a, size_t b) { return toa_s[a] < toa_s[b]; });
    std::vector<double> t, w;
    for (size_t i : g) {
      t.push_back(toa_s[i]);
      w.push_back(pw_s[i]);
    }
    auto [em, rest] = emitters(t, w);
    for (const Json& e : em) {
      if (e["n_pulses"].get<int64_t>() >= kMinReport)
        ems.push_back(e);
      else
        left += e["n_pulses"].get<int64_t>();
    }
    left += static_cast<int64_t>(rest.size());
  }
  std::stable_sort(ems.begin(), ems.end(),
                   [](const Json& a, const Json& b) { return a["n_pulses"].get<int64_t>() > b["n_pulses"].get<int64_t>(); });
  Json out = Json::object();
  out["emitters"] = ems;
  out["unassigned"] = left;
  return out;
}

}  // namespace si::deint
