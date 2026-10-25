// One block end to end (python: analyze.py): detect -> features -> receivers -> classify.
#include "analyze.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <future>
#include <map>
#include <cmath>
#include <limits>
#include <set>
#include <thread>

#include "classify/rules.hpp"
#include "core/parallel.hpp"
#include "core/pyfmt.hpp"
#include "core/spurs.hpp"
#include "detect/artifacts.hpp"
#include "detect/merge.hpp"
#include "detect/pipeline.hpp"
#include "features/measure.hpp"

namespace si {
namespace {

using Clock = std::chrono::steady_clock;

double since(Clock::time_point t) { return py_round(std::chrono::duration<double>(Clock::now() - t).count(), 2); }

std::string err_text(const std::string& prefix, const std::exception& e) {
  std::string s = prefix + e.what();
  return s.substr(0, 200);
}

bool continuous(const Signal& s) {
  return std::any_of(s.bursts.begin(), s.bursts.end(), [](const Box& b) { return b.kind == "continuous"; });
}

constexpr double kCoincTS = 20e-6;         // pieces of one burst start and end together
constexpr size_t kCoincMinN = 3;           // a wideband burst split at fading notches leaves >= 3 pieces
constexpr double kCoincMinSpanHz = 8e6;    // spanning >= 8 MHz together (a narrow hop / RC packet is one piece)
constexpr double kCoincFov = 0.8;          // each piece >= 80 % inside the receiver signal's span

// python _absorb_coincident: >= 3 burst-only pieces inside one receiver signal's span whose bursts start and end
// together, spanning >= 8 MHz, are one wideband burst of that receiver's emitter
void absorb_coincident(const std::vector<Signal>& receivers, const std::vector<Signal>& energy,
                       const std::vector<size_t>& alive, std::map<size_t, std::vector<size_t>>& taken) {
  std::set<size_t> used;
  for (const auto& [k, v] : taken) used.insert(v.begin(), v.end());
  auto match = [&](const Signal& a, const Signal& b) {
    for (const auto& x : a.bursts) {
      bool hit = false;
      for (const auto& y : b.bursts)
        if (std::fabs(x.t0 - y.t0) <= kCoincTS && std::fabs(x.t1 - y.t1) <= kCoincTS) {
          hit = true;
          break;
        }
      if (!hit) return false;
    }
    return true;
  };
  for (size_t k = 0; k < receivers.size(); ++k) {
    const Signal& r = receivers[k];
    std::vector<size_t> inside;
    for (size_t ei : alive) {
      const Signal& e = energy[ei];
      if (used.count(ei) || !e.flags.empty() || e.bursts.empty() || continuous(e)) continue;
      if (std::max(0.0, std::min(e.f1, r.f1) - std::max(e.f0, r.f0)) >= kCoincFov * std::max(e.bw(), 1.0))
        inside.push_back(ei);
    }
    std::stable_sort(inside.begin(), inside.end(),
                     [&](size_t a, size_t b) { return energy[a].bursts.size() > energy[b].bursts.size(); });
    std::vector<std::vector<size_t>> groups;
    for (size_t ei : inside) {
      auto g = std::find_if(groups.begin(), groups.end(), [&](const auto& g) { return match(energy[ei], energy[g[0]]); });
      if (g == groups.end()) groups.push_back({ei});
      else g->push_back(ei);
    }
    for (const auto& g : groups) {
      double lo = std::numeric_limits<double>::infinity(), hi = -lo;
      for (size_t ei : g) {
        lo = std::min(lo, energy[ei].f0);
        hi = std::max(hi, energy[ei].f1);
      }
      if (g.size() >= kCoincMinN && hi - lo >= kCoincMinSpanHz) {
        auto& t = taken[k];
        t.insert(t.end(), g.begin(), g.end());
        used.insert(g.begin(), g.end());
      }
    }
  }
}

// Assign energy detections to the receiver-defined signals that explain them (python _absorb):
//   burst energy:      >= 50 % of its band inside a receiver signal's span and >= 30 % of its
//                      on-time covered by that receiver's bursts; or the energy box contains the
//                      receiver's whole band (centre inside it) and >= 60 % of its on-time is covered
//   continuous energy: explained when receiver spans cover >= 60 % of its band and those receivers
//                      are active >= 5 % of the capture (channel occupancy).
std::map<size_t, std::vector<size_t>> absorb(const std::vector<Signal>& receivers, const std::vector<Signal>& energy,
                                             const std::vector<size_t>& alive, double cap_s, double pad = 60e-6) {
  std::map<size_t, std::vector<size_t>> taken;
  for (size_t ei : alive) {
    const Signal& e = energy[ei];
    // a box cut at the capture edge is still an emission, but only a receiver that owns its band explains it;
    // artifact boxes are never absorbed
    const bool edge_cut = e.flags.size() == 1 && e.flags[0] == "truncated_at_capture_edge";
    if (!e.flags.empty() && !edge_cut) continue;
    std::optional<size_t> best;
    double best_score = 0.0;
    const bool cont = continuous(e);
    double cov_f = 0.0;
    for (size_t k = 0; k < receivers.size(); ++k) {
      const Signal& r = receivers[k];
      if (edge_cut && !r.features.value("owns_band", false)) continue;
      const double fov = std::max(0.0, std::min(e.f1, r.f1) - std::max(e.f0, r.f0)) / std::max(e.bw(), 1.0);
      if (fov <= 0) continue;
      double on_r = 0.0;
      for (const auto& b : r.bursts) on_r += b.t1 - b.t0;
      on_r /= cap_s;
      double score;
      const bool owns = r.features.value("owns_band", false);
      if (owns && fov >= 0.8) {
        score = fov;  // licensed carrier: everything inside its band is it
      } else if (cont) {
        cov_f += on_r >= 0.05 ? fov : 0.0;
        score = fov * (on_r >= 0.05 ? 1.0 : 0.0);
      } else {
        double on = 0.0, cov = 0.0;
        for (const auto& b : e.bursts) on += b.t1 - b.t0;
        for (const auto& b : e.bursts)
          for (const auto& x : r.bursts) cov += std::max(0.0, std::min(b.t1, x.t1 + pad) - std::max(b.t0, x.t0 - pad));
        score = (fov >= 0.5 && on > 0 && cov / on >= 0.3) ? fov : 0.0;
        if (score == 0.0 && on > 0 && cov / on >= 0.6 && e.f0 <= r.f0 && e.f1 >= r.f1 && r.f0 <= e.fc() &&
            e.fc() <= r.f1)
          score = 0.5;  // receiver band inside a wider, time-aligned box
      }
      if (score > best_score) {
        best = k;
        best_score = score;
      }
    }
    const bool owned = best && receivers[*best].features.value("owns_band", false);
    if (best && (owned || !cont || std::min(cov_f, 1.0) >= 0.6)) taken[*best].push_back(ei);
  }
  absorb_coincident(receivers, energy, alive, taken);
  return taken;
}

// python _image_caution(s, lo): single-block LO-harmonic hint (core/spurs). An unidentified signal outside the
// cellular bands that a 3rd / 5th LO-harmonic image of a cellular downlink would explain gets a caution; label and
// confidence stay.
void image_caution(Signal& s, double lo) {
  if (s.decoded.status != "not_possible") return;
  for (const char* f : kArtifactFlags)
    if (s.has_flag(f)) return;
  if (spurs::in_strong_band(s.f0, s.f1)) return;
  const Json src = spurs::harmonic_image_sources(s.f0, s.f1, lo);
  if (src.empty()) return;
  s.flags.push_back("possible_lo_harmonic_image");
  s.features["lo_harmonic_image_sources"] = src;
  std::string parts;
  for (size_t i = 0; i < src.size(); ++i) {
    const Json& c = src[i];
    if (i) parts += ", ";
    parts += fmt("%.1f-%.1f MHz (", c["source_low_hz"].get<int64_t>() / 1e6, c["source_high_hz"].get<int64_t>() / 1e6) +
             c["band"].get<std::string>() + fmt(", %d x LO", c["n"].get<int>()) +
             (c["inverted"].get<bool>() ? ", inverted" : "") + ")";
  }
  s.evidence = (s.evidence.empty() ? std::string() : s.evidence + "; ") + "CAUTION possible LO-harmonic image of " +
               parts + " (no preselector: confirm with a second tuning)";
}

std::string why_not(const Signal& s, const std::vector<std::string>& attempts) {
  if (s.has_flag("lo_leakage") || s.has_flag("iq_image")) return "receiver artifact, not an emission";
  const std::string& lab = s.label;
  if (lab.rfind("CW", 0) == 0) return "unmodulated carrier: no information content";
  if (s.snr_db < 6) return "insufficient SNR (" + py_float_repr(s.snr_db) + " dB) for demodulation";
  if (s.has_flag("truncated_at_capture_edge")) return "signal extends beyond the captured bandwidth (truncated)";
  static const std::vector<std::pair<const char*, const char*>> known = {
      {"LTE", "no LTE/NR decoder matched (PSS/SSS not found in this span)"},
      {"GSM", "GSM FCCH/SCH not present in this 32.5 ms window"},
      {"UMTS", "WCDMA decoding (scrambling-code search, P-CCPCH) not implemented"},
      {"DVB-T2", "DVB-T2 L1 decoding not implemented; 32 ms covers < 1 T2 frame"},
      {"PMR", "32.5 ms is too short for voice; digital PMR (DMR/TETRA) vocoder/encryption not handled"},
      {"TETRA", "TETRA air interface is typically encrypted (TEA); 32.5 ms covers < 1 multiframe"},
      {"radar", "pulsed radar carries no decodable payload; PRI/pulse width reported"},
      {"DME", "DME pulse pairs carry no payload; X/Y channel from pair spacing only"},
      {"Bluetooth", "classic BT / non-advertising BLE uses hopping + unknown access address"},
      {"ISM", "proprietary SRD protocol (OOK/FSK) — no decoder for the specific device"},
  };
  for (const auto& [k, v] : known)
    if (lab.find(k) != std::string::npos) return v;
  if (lab == "unknown") return "unknown protocol: no decoder applicable";
  std::string out = "no decoder implemented for this signal class";
  if (!attempts.empty()) {
    out += "; tried: ";
    for (size_t i = 0; i < attempts.size(); ++i) out += (i ? ", " : "") + attempts[i];
  }
  return out;
}

// capture_signals() of every plugin, grouped by receiver package. Plugins of one package run in
// order; packages run concurrently, each on its own copy of the capture cache, which sees only
// what the packages named in a plugin's after() published (python analyze._capture_signals).
std::vector<Signal> capture_signals(DecodeContext& ctx, const PluginList& plugins, Json& T) {
  std::vector<std::string> names;
  std::map<std::string, std::vector<const Decoder*>> groups;
  for (const auto& p : plugins) {
    if (!groups.count(p->package())) names.push_back(p->package());
    groups[p->package()].push_back(p.get());
  }
  const size_t G = names.size();
  const Cache base = ctx.cache;
  std::vector<Cache> caches(G, base);
  std::vector<std::vector<Signal>> results(G);
  std::map<std::string, double> timings;
  std::map<std::string, std::string> errors;
  std::mutex mu;
  std::vector<std::promise<void>> done(G);
  std::vector<std::shared_future<void>> ready;
  for (auto& d : done) ready.push_back(d.get_future().share());
  auto index_of = [&](const std::string& n) -> std::optional<size_t> {
    auto it = std::find(names.begin(), names.end(), n);
    return it == names.end() ? std::nullopt : std::optional<size_t>(it - names.begin());
  };
  auto run = [&](size_t g) {
    Cache cache = base;
    for (const auto* dec : groups[names[g]])
      for (const auto& d : dec->after()) {
        auto j = index_of(d);
        if (!j || *j == g) continue;
        ready[*j].wait();  // packages earlier in registry order: already running or done
        for (const auto& [k, v] : caches[*j].items())
          if (!base.contains(k)) cache.items()[k] = v;
      }
    DecodeContext gctx{ctx.cap, ctx.ch, ctx.det, std::move(cache)};
    std::vector<Signal> got;
    for (const auto* dec : groups[names[g]]) {
      const auto t = Clock::now();
      try {
        auto s = dec->capture_signals(gctx);
        got.insert(got.end(), std::make_move_iterator(s.begin()), std::make_move_iterator(s.end()));
      } catch (const std::exception& e) {
        std::lock_guard<std::mutex> lk(mu);
        errors[dec->name()] = err_text("", e);
      }
      std::lock_guard<std::mutex> lk(mu);
      timings[dec->name()] = since(t);
    }
    caches[g] = std::move(gctx.cache);
    results[g] = std::move(got);
    done[g].set_value();
  };
  const size_t n = std::min<size_t>(static_cast<size_t>(threads()), G);
  if (n <= 1) {
    for (size_t g = 0; g < G; ++g) run(g);
  } else {
    // FIFO in registry order: whatever a package waits on started before it
    std::atomic<size_t> next{0};
    std::vector<std::thread> pool;
    for (size_t i = 0; i < n; ++i)
      pool.emplace_back([&] {
        for (size_t g = next.fetch_add(1); g < G; g = next.fetch_add(1)) run(g);
      });
    for (auto& th : pool) th.join();
  }
  std::vector<Signal> out;
  for (size_t g = 0; g < G; ++g) {
    out.insert(out.end(), std::make_move_iterator(results[g].begin()), std::make_move_iterator(results[g].end()));
    for (const auto& [k, v] : caches[g].items())
      if (!base.contains(k)) ctx.cache.items()[k] = v;
  }
  for (const auto& p : plugins) {
    if (errors.count(p->name())) T[p->name() + "_error"] = errors[p->name()];
    T["receiver_s"][p->name()] = timings.count(p->name()) ? timings[p->name()] : 0.0;
  }
  return out;
}

Analysis analyze_once(const Capture& in, const PluginList& plugins, std::optional<double> fc_override,
                      const std::optional<Json>& ppm) {
  Json T = Json::object();
  auto t = Clock::now();
  Analysis A;
  A.cap = in;  // shares the samples
  if (fc_override) A.cap.fc = *fc_override;
  const Capture& cap = A.cap;
  A.quality = detect::clipping(cap.iq());
  A.det = detect::detect(cap);
  T["detect_s"] = since(t);

  t = Clock::now();
  Channelizer ch(cap.iq(), cap.fs, cap.fc);
  std::vector<Signal>& energy = A.det.signals;
  parallel_for(energy.size(),
               [&](size_t i) { features::extract(energy[i], ch, cap.n(), cap.duration(), A.det.fine.df()); });
  {
    auto [sigs, merged] = detect::merge_pulse_trains(std::move(energy));
    energy = std::move(sigs);
    for (size_t i : merged) features::extract(energy[i], ch, cap.n(), cap.duration(), A.det.fine.df());
  }
  T["features_s"] = since(t);

  // 1) receiver-defined signals absorb the energy signals they explain
  t = Clock::now();
  DecodeContext ctx{cap, ch, &A.det, Cache()};
  if (ppm) ctx.cache.set<Json>("sweep_ppm", *ppm);  // receiver frequency error (caller / sweep calibration)
  // lines from the receiver's own clocks (the caller's / sweep's ppm only, so the result does not depend on which
  // receivers ran)
  std::optional<double> ppm_v;
  if (ppm && ppm->is_object() && ppm->contains("ppm") && (*ppm)["ppm"].is_number()) ppm_v = (*ppm)["ppm"].get<double>();
  {
    Json lines = spurs::receiver_lines(ch.X(), cap.fs, cap.fc, ppm_v, A.det.usable_hz);
    detect::flag_receiver_spurs(energy, lines, cap.fs, ppm_v);
    ctx.cache.set<Json>("receiver_lines", std::move(lines));  // receivers skip these (e.g. voice_fm)
  }
  std::vector<Signal> out = capture_signals(ctx, plugins, T);
  for (auto& r : out)
    for (auto& b : r.bursts) b.t1 = std::min(b.t1, cap.duration());
  for (const auto& dec : plugins) {
    try {
      out = dec->arbitrate(std::move(out), ctx);
    } catch (const std::exception& e) {
      T[dec->name() + "_arbitrate_error"] = err_text("", e);
    }
  }
  std::vector<size_t> alive(energy.size());
  for (size_t i = 0; i < alive.size(); ++i) alive[i] = i;
  const auto absorbed = absorb(out, energy, alive, cap.duration());
  std::set<size_t> eaten_all;
  for (size_t k = 0; k < out.size(); ++k) {
    Signal& r = out[k];
    auto it = absorbed.find(k);
    if (it != absorbed.end() && !it->second.empty()) {
      const auto& eaten = it->second;
      Json lst = Json::array();
      double mx = -1e300;
      for (size_t ei : eaten) {
        const Signal& e = energy[ei];
        lst.push_back({{"f_low_hz", py_round_int(e.f0)},
                       {"f_high_hz", py_round_int(e.f1)},
                       {"t_start_s", py_round(e.t0(), 6)},
                       {"t_end_s", py_round(e.t1(), 6)},
                       {"snr_db", e.snr_db},
                       {"kind", continuous(e) ? "channel occupancy (aggregate traffic)" : "burst"}});
        mx = std::max(mx, e.snr_db);
      }
      r.detection["absorbed_energy_detections"] = lst;
      r.detection["description"] = r.detection.value("description", std::string()) + "; explains " +
                                   std::to_string(eaten.size()) +
                                   " energy detection(s) (STFT/CFAR) in the same band and time";
      if (r.snr_db == 0) r.snr_db = mx;
      for (const char* key : {"envelope_cv", "cyclic_line_hz", "ofdm_scs_hz"}) {
        for (size_t ei : eaten) {
          const Json& ef = energy[ei].features;
          if (ef.contains(key) && !ef[key].is_null()) {
            if (!r.features.contains(key)) r.features[key] = ef[key];
            break;
          }
        }
      }
      eaten_all.insert(eaten.begin(), eaten.end());
    }
    const Json tm = features::timing(r, cap.duration());
    for (auto jt = tm.begin(); jt != tm.end(); ++jt)
      if (!r.features.contains(jt.key()) || jt.key() == "n_bursts" || jt.key() == "burst_period_s")
        r.features[jt.key()] = jt.value();
  }
  T["receivers_s"] = since(t);

  t = Clock::now();
  Json extras = Json::object();
  for (const auto& dec : plugins) {
    try {
      Json x = dec->postprocess(ctx, out);
      for (auto it = x.begin(); it != x.end(); ++it) extras[it.key()] = it.value();
    } catch (const std::exception& e) {
      T[dec->name() + "_postprocess_error"] = err_text("", e);
    }
  }
  T["postprocess_s"] = since(t);

  // 2) per-signal decoders on what is left, then rules
  t = Clock::now();
  for (size_t ei = 0; ei < energy.size(); ++ei) {
    if (eaten_all.count(ei)) continue;
    Signal s = std::move(energy[ei]);
    std::vector<std::string> attempts;
    if (s.flags.empty()) {
      for (const auto& dec : plugins) {
        DecodeResult r;
        try {
          if (dec->capture_level() || !dec->can_decode(s, ctx)) continue;
          r = dec->decode(s, ctx);
        } catch (const std::exception& e) {
          attempts.push_back(err_text(dec->name() + ": error ", e));
          continue;
        }
        attempts.push_back(dec->name() + ": " + r.status + (r.reason.empty() ? "" : " (" + r.reason + ")"));
        if (r.status == "success" || r.status == "partial") {
          s.decoded.status = r.status;
          s.decoded.reason = r.reason;
          s.decoded.method = r.method;
          s.decoded.steps = r.steps;
          s.decoded.fields = r.fields;
          if (r.label && !r.label->empty()) {
            s.label = *r.label;
            s.confidence = r.confidence.value_or(0.0);
            if (r.modulation && !r.modulation->empty()) s.modulation = *r.modulation;
            s.evidence = r.evidence.value_or("");
          }
          if (r.span) s.features["decoder_span_hz"] = Json::array({py_round_int(r.span->first), py_round_int(r.span->second)});
          break;
        }
      }
    }
    if (s.decoded.status == "not_possible" || !s.decoded.method || s.decoded.method->empty()) {
      classify::classify(s);
      s.decoded.status = "not_possible";
      s.decoded.reason = why_not(s, attempts);
      s.decoded.method.reset();
      s.decoded.steps = attempts;
      s.decoded.fields = Json::object();
    } else {
      std::vector<std::string> steps = attempts;
      steps.insert(steps.end(), s.decoded.steps.begin(), s.decoded.steps.end());
      s.decoded.steps = steps;
    }
    out.push_back(std::move(s));
  }
  T["decode_classify_s"] = since(t);

  for (auto& s : out) image_caution(s, cap.fc);
  if (A.quality.value("clipping", false))
    for (auto& s : out) s.flags.push_back("adc_clipping_in_capture");
  std::stable_sort(out.begin(), out.end(),
                   [](const Signal& a, const Signal& b) { return a.f0 != b.f0 ? a.f0 < b.f0 : a.t0() < b.t0(); });
  A.det.signals.clear();  // moved into out; the report does not need them again
  A.signals = std::move(out);
  A.timings = T;
  A.extras = extras;
  for (const auto& p : plugins) A.plugins.push_back(p.get());
  return A;
}

std::optional<Json> label_check(const Analysis& a, const PluginList& plugins) {
  for (const auto& p : plugins) {
    std::optional<Json> chk;
    try {
      chk = p->frequency_label_check(a.signals);
    } catch (...) {
      chk.reset();
    }
    if (chk && !chk->is_null() && !chk->empty()) {
      Json c = *chk;
      c["plugin"] = p->name();
      return c;
    }
  }
  return std::nullopt;
}

}  // namespace

Analysis analyze(const Capture& cap, const PluginList& plugins, const std::optional<Json>& ppm) {
  Analysis a = analyze_once(cap, plugins, std::nullopt, ppm);
  auto chk = label_check(a, plugins);
  if (chk && chk->value("suspect", false)) {
    Analysis b = analyze_once(cap, plugins, a.cap.fc - (*chk)["offset_mhz"].get<double>() * 1e6, ppm);
    Json lc = *chk;
    lc["status"] = "corrected_by_" + (*chk)["plugin"].get<std::string>();
    lc["fc_metadata_hz"] = a.cap.meta["captures"][0]["core:frequency"].get<double>();
    lc["fc_used_hz"] = b.cap.fc;
    b.label_check = lc;
    double first = 0.0;
    for (auto it = a.timings.begin(); it != a.timings.end(); ++it)
      if (it.value().is_number_float()) first += it.value().get<double>();
    b.timings["first_pass_s"] = py_round(first, 2);
    return b;
  }
  if (chk) {
    Json lc = *chk;
    lc["status"] = "consistent";
    a.label_check = lc;
  }
  return a;
}

}  // namespace si
