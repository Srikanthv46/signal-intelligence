#include <spectrum_intel/c_api.h>

#include <cstdlib>
#include <cstring>
#include <string>

#include <spectrum_intel/spectrum_intel.hpp>

struct si_engine {
  si::Engine engine;
  explicit si_engine(si::EngineConfig c) : engine(std::move(c)) {}
};

namespace {

char* dup(const std::string& s) {
  char* p = static_cast<char*>(std::malloc(s.size() + 1));
  if (p) std::memcpy(p, s.c_str(), s.size() + 1);
  return p;
}

si::IqFormat format_of(int f) {
  switch (f) {
    case SI_CF64: return si::IqFormat::CF64;
    case SI_CI16: return si::IqFormat::CI16;
    case SI_CI8: return si::IqFormat::CI8;
    case SI_CU8: return si::IqFormat::CU8;
    default: return si::IqFormat::CF32;
  }
}

}  // namespace

extern "C" {

si_engine* si_engine_create(const char* config_json, char** error) {
  try {
    si::EngineConfig cfg;
    if (config_json && *config_json) {
      si::Json j = si::Json::parse(config_json);
      if (j.contains("enable")) cfg.enable = j["enable"].get<std::vector<std::string>>();
      if (j.contains("disable")) cfg.disable = j["disable"].get<std::vector<std::string>>();
      cfg.threads = j.value("threads", 0);
      cfg.compact = j.value("compact", false);
    }
    return new si_engine(std::move(cfg));
  } catch (const std::exception& e) {
    if (error) *error = dup(e.what());
    return nullptr;
  }
}

void si_engine_destroy(si_engine* e) { delete e; }

char* si_analyze(const si_engine* e, const void* iq, size_t n, int format, double fs, double fc, double bw,
                 const char* meta_json) {
  si::IQBlock b;
  b.data = iq;
  b.n_samples = n;
  b.format = format_of(format);
  b.sample_rate_hz = fs;
  b.center_frequency_hz = fc;
  if (bw > 0) b.bandwidth_hz = bw;
  try {
    if (meta_json && *meta_json) {
      si::Json m = si::Json::parse(meta_json);
      if (m.contains("start_time_utc") && m["start_time_utc"].is_string()) b.start_time_utc = m["start_time_utc"].get<std::string>();
      if (m.contains("gain_db") && m["gain_db"].is_number()) b.gain_db = m["gain_db"].get<double>();
      if (m.contains("receiver_ppm") && m["receiver_ppm"].is_number()) b.receiver_ppm = m["receiver_ppm"].get<double>();
      if (m.contains("source_id") && m["source_id"].is_string()) b.source_id = m["source_id"].get<std::string>();
      if (m.contains("request_id") && m["request_id"].is_string()) b.request_id = m["request_id"].get<std::string>();
    }
  } catch (const std::exception& ex) {
    si::Json r = {{"status", "insufficient_input"},
                  {"errors", si::Json::array({{{"code", "INVALID_META"}, {"receiver", nullptr}, {"message", ex.what()}}})}};
    return dup(r.dump());
  }
  if (!e) return dup(R"({"status": "error", "errors": [{"code": "NO_ENGINE", "receiver": null, "message": "engine is NULL"}]})");
  return dup(e->engine.analyze_json(b));
}

char* si_requirements(void) { return dup(si::requirements().dump()); }

const char* si_version(void) { return si::version(); }

void si_free(char* p) { std::free(p); }

}  // extern "C"
