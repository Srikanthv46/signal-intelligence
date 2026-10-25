// Production entry point (python: api/engine.py and api/stream.py).
#include <spectrum_intel/spectrum_intel.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>

#include "analyze.hpp"
#include "api/contract.hpp"
#include "api/response.hpp"
#include "core/parallel.hpp"
#include "decode/registry.hpp"
#include "io/sigmf.hpp"
#include "report/writer.hpp"

#ifndef SI_VERSION
#define SI_VERSION "0.0.0"
#endif

namespace si {
namespace {

using Clock = std::chrono::steady_clock;

double seconds_since(Clock::time_point t) {
  return std::chrono::duration<double>(Clock::now() - t).count();
}

const char* datatype(IqFormat f) {
  switch (f) {
    case IqFormat::CF32: return "cf32_le";
    case IqFormat::CF64: return "cf64_le";
    case IqFormat::CI16: return "ci16_le";
    case IqFormat::CI8: return "ci8";
    case IqFormat::CU8: return "cu8";
  }
  return "cf32_le";
}

api::BlockMeta meta_of(const IQBlock& b) {
  api::BlockMeta m;
  m.center_frequency_hz = b.center_frequency_hz;
  m.sample_rate_hz = b.sample_rate_hz;
  m.bandwidth_hz = b.bandwidth_hz;
  m.gain_db = b.gain_db;
  m.receiver_ppm = b.receiver_ppm;
  if (b.start_time_utc) m.start = UtcTime::parse(*b.start_time_utc);
  m.source_id = b.source_id;
  m.request_id = b.request_id;
  m.n_samples = b.n_samples;
  return m;
}

// Why this block cannot be analysed ("" if it can). python api.engine.validate
std::string validate(const IQBlock& b, const cfvec* x) {
  if (!std::isfinite(b.sample_rate_hz)) return "sample_rate_hz must be a finite number";
  if (!std::isfinite(b.center_frequency_hz)) return "center_frequency_hz must be a finite number";
  if (b.sample_rate_hz <= 0) return "sample_rate_hz must be > 0, got " + std::to_string(b.sample_rate_hz);
  if (b.bandwidth_hz && !(std::isfinite(*b.bandwidth_hz) && *b.bandwidth_hz > 0))
    return "bandwidth_hz must be > 0 or unset";
  if (b.receiver_ppm && !(std::isfinite(*b.receiver_ppm) && std::fabs(*b.receiver_ppm) < 1000))
    return "receiver_ppm must be a finite number of ppm (|x| < 1000)";
  if (b.start_time_utc && !UtcTime::parse(*b.start_time_utc))
    return "start_time_utc is not an ISO 8601 time: " + *b.start_time_utc;
  if (!b.data && b.n_samples) return "iq: data is null";
  const Json& p = contract::requirements().at("pipeline");
  const size_t need = p.at("min_samples").get<size_t>();
  if (b.n_samples < need)
    return std::to_string(b.n_samples) + " samples is below the minimum of " + std::to_string(need) +
           " for analysis (recommended " + std::to_string(p.at("recommended_samples").get<int64_t>()) +
           "); accumulate blocks with StreamAnalyzer";
  if (x)
    for (const auto& v : *x)
      if (!std::isfinite(v.real()) || !std::isfinite(v.imag())) return "iq contains NaN or Inf";
  return "";
}

}  // namespace

const char* version() { return SI_VERSION; }

std::vector<std::string> available_receivers() { return receiver_names(); }

const Json& requirements() { return contract::requirements(); }

struct Engine::Impl {
  EngineConfig cfg;
  PluginList plugins;
  api::Packages packages;
  std::vector<std::string> names;
};

Engine::Engine(EngineConfig cfg) : impl_(std::make_unique<Impl>()) {
  impl_->cfg = std::move(cfg);
  if (impl_->cfg.threads > 0) set_threads(impl_->cfg.threads);
  impl_->plugins = plugins(impl_->cfg.enable, impl_->cfg.disable);
  for (const auto& p : impl_->plugins) {
    const std::string pkg = p->package();
    auto it = std::find_if(impl_->packages.begin(), impl_->packages.end(),
                           [&](const auto& e) { return e.first == pkg; });
    if (it == impl_->packages.end()) {
      impl_->packages.push_back({pkg, {p->name()}});
      impl_->names.push_back(pkg);
    } else {
      it->second.push_back(p->name());
    }
  }
}

Engine::~Engine() = default;

std::vector<std::string> Engine::receivers() const { return impl_->names; }

Json Engine::info() const {
  return {{"name", "spectrum_intel"},
          {"implementation", "cpp"},
          {"version", SI_VERSION},
          {"receivers_enabled", impl_->names},
          {"receivers_unavailable", Json::object()},
          {"threads", threads()}};
}

Json Engine::analyze(const IQBlock& b) const {
  const auto t = Clock::now();
  const api::BlockMeta meta = meta_of(b);
  std::string why = validate(b, nullptr);
  if (!why.empty()) return api::failure_response(meta, info(), "insufficient_input", "INVALID_INPUT", why, seconds_since(t));
  try {
    cfvec x = convert_iq(b.data, b.n_samples, datatype(b.format));
    why = validate(b, &x);
    if (!why.empty())
      return api::failure_response(meta, info(), "insufficient_input", "INVALID_INPUT", why, seconds_since(t));
    Capture cap = Capture::from_iq(std::move(x), b.sample_rate_hz, b.center_frequency_hz, meta.start,
                                   b.source_id.value_or("block"), b.gain_db, b.bandwidth_hz);
    std::optional<Json> ppm;
    if (b.receiver_ppm) ppm = Json{{"ppm", *b.receiver_ppm}, {"source", "caller (IQBlock.receiver_ppm)"}};
    Analysis a = si::analyze(cap, impl_->plugins, ppm);
    for (auto& s : a.signals) {  // binary products are for files, not the response
      s.video_frame.reset();
      s.audio.reset();
    }
    Json rep = build_report(a);
    return api::build_response(rep, a, meta, info(), impl_->packages, seconds_since(t), impl_->cfg.compact);
  } catch (const std::exception& e) {
    return api::failure_response(meta, info(), "error", "ANALYSIS_FAILED", e.what(), seconds_since(t));
  } catch (...) {
    return api::failure_response(meta, info(), "error", "ANALYSIS_FAILED", "unknown exception", seconds_since(t));
  }
}

std::string Engine::analyze_json(const IQBlock& b, int indent) const {
  return analyze(b).dump(indent, ' ', false, Json::error_handler_t::replace);
}

// ---------------------------------------------------------------- StreamAnalyzer
struct StreamAnalyzer::Impl {
  const Engine* engine;
  size_t target, max_gap;
  bool have_key = false;
  double fs = 0, fc = 0;
  std::optional<double> bw;
  std::optional<UtcTime> t0;           // time of the first buffered sample (or of the next one, when empty)
  cfvec buf;
  IQBlock meta;                        // gain / source / request ids of the current block

  bool same_key(const IQBlock& b) const {
    return have_key && b.sample_rate_hz == fs && b.center_frequency_hz == fc && b.bandwidth_hz == bw;
  }
  bool gap(const std::optional<UtcTime>& t) const {
    if (!t || !t0 || buf.empty()) return false;
    UtcTime expected = t0->plus(static_cast<double>(buf.size()) / fs);
    double d = static_cast<double>(t->epoch_s - expected.epoch_s) + (t->micro - expected.micro) * 1e-6;
    return std::fabs(d) * fs > static_cast<double>(max_gap);
  }
  Json analyse(bool keep_going) {
    IQBlock b = meta;
    b.data = buf.data();
    b.n_samples = buf.size();
    b.format = IqFormat::CF32;
    b.sample_rate_hz = fs;
    b.center_frequency_hz = fc;
    b.bandwidth_hz = bw;
    b.start_time_utc = t0 ? std::optional<std::string>(t0->iso_z()) : std::nullopt;
    Json r = engine->analyze(b);
    if (keep_going) {
      if (t0) t0 = t0->plus(static_cast<double>(buf.size()) / fs);
      buf.clear();
    }
    return r;
  }
};

StreamAnalyzer::StreamAnalyzer(const Engine& engine, size_t target_samples, size_t max_gap_samples)
    : impl_(std::make_unique<Impl>()) {
  const Json& p = requirements().at("pipeline");
  impl_->engine = &engine;
  impl_->target = target_samples ? target_samples : p.at("recommended_samples").get<size_t>();
  if (impl_->target < p.at("min_samples").get<size_t>())
    throw std::invalid_argument("target_samples must be >= " + p.at("min_samples").dump());
  impl_->max_gap = max_gap_samples;
  impl_->buf.reserve(impl_->target);
}

StreamAnalyzer::~StreamAnalyzer() = default;

size_t StreamAnalyzer::buffered() const { return impl_->buf.size(); }

std::vector<Json> StreamAnalyzer::push(const IQBlock& b) {
  auto& s = *impl_;
  std::vector<Json> out;
  cfvec x = convert_iq(b.data, b.n_samples, datatype(b.format));
  std::optional<UtcTime> t;
  if (b.start_time_utc) t = UtcTime::parse(*b.start_time_utc);
  if (s.have_key && (!s.same_key(b) || s.gap(t))) {
    if (auto r = flush()) out.push_back(std::move(*r));
  }
  if (!s.have_key || s.buf.empty()) {
    const bool continues = s.same_key(b);
    s.have_key = true;
    s.fs = b.sample_rate_hz;
    s.fc = b.center_frequency_hz;
    s.bw = b.bandwidth_hz;
    s.meta = b;
    s.t0 = t ? t : (continues ? s.t0 : std::nullopt);
  }
  size_t pos = 0;
  while (pos < x.size()) {
    const size_t take = std::min(x.size() - pos, s.target - s.buf.size());
    s.buf.insert(s.buf.end(), x.begin() + static_cast<ptrdiff_t>(pos), x.begin() + static_cast<ptrdiff_t>(pos + take));
    pos += take;
    if (s.buf.size() >= s.target) out.push_back(s.analyse(true));
  }
  return out;
}

std::optional<Json> StreamAnalyzer::flush() {
  auto& s = *impl_;
  std::optional<Json> r;
  if (s.buf.size() >= requirements().at("pipeline").at("min_samples").get<size_t>()) r = s.analyse(false);
  s.have_key = false;
  s.t0.reset();
  s.buf.clear();
  return r;
}

}  // namespace si
