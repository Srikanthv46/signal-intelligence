// spectrum_intel — detect, localise, identify and decode every signal in an IQ block.
//
// Production C++ API. One IQ block (samples + sample rate + centre frequency + analog
// bandwidth) in, one JSON response out; the response follows schema/response.schema.json
// and is the same as the Python implementation's (python/spectrum_intel/api).
//
//   si::Engine engine;                                   // load receivers once
//   si::IQBlock b;
//   b.data = samples; b.n_samples = n; b.format = si::IqFormat::CI16;
//   b.sample_rate_hz = 61.44e6; b.center_frequency_hz = 2.437e9; b.bandwidth_hz = 56e6;
//   b.start_time_utc = "2026-10-09T10:00:00.000000Z";
//   std::string json = engine.analyze_json(b);           // never throws for bad input
//
// Blocks of a few thousand samples are too short to analyse: accumulate them with
// StreamAnalyzer (see requirements() for what each receiver needs).
#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

namespace si {

using Json = nlohmann::ordered_json;

const char* version();

// Sample formats an SDR pipeline delivers (SigMF names). Values are not rescaled: raw
// ADC counts keep the clipping check meaningful.
enum class IqFormat {
  CF32,  // interleaved float32 I/Q (std::complex<float>)
  CF64,  // interleaved float64 I/Q
  CI16,  // interleaved int16 I/Q
  CI8,   // interleaved int8 I/Q
  CU8,   // interleaved uint8 I/Q (offset binary, centred on 127.5)
};

// One block from the acquisition pipeline. The samples are borrowed (not copied until
// analysis) and must stay valid during the call.
struct IQBlock {
  const void* data = nullptr;
  size_t n_samples = 0;                       // complex samples (pairs)
  IqFormat format = IqFormat::CF32;
  double sample_rate_hz = 0;                  // complex sample rate
  double center_frequency_hz = 0;             // RF of baseband 0 Hz (the LO); trusted as given
  std::optional<double> bandwidth_hz;         // analog front-end bandwidth; unset = whole sampled band
  std::optional<std::string> start_time_utc;  // ISO 8601 time of the first sample
  std::optional<double> receiver_ppm;         // the SDR's frequency error (ppm), if known; else estimated per block
  std::optional<double> gain_db;
  std::optional<std::string> source_id;
  std::optional<std::string> request_id;
};

struct EngineConfig {
  std::vector<std::string> enable;   // receivers to run (empty = all built in)
  std::vector<std::string> disable;  // receivers to skip
  int threads = 0;                   // 0 = SPECTRUM_INTEL_THREADS or min(4, cores)
  bool compact = false;              // drop provenance detail (features, detection, steps, bursts)
};

class Engine {
 public:
  explicit Engine(EngineConfig cfg = {});  // throws std::invalid_argument for unknown receiver names
  ~Engine();
  Engine(const Engine&) = delete;
  Engine& operator=(const Engine&) = delete;

  // Analyse one block. Never throws for bad input or a failing receiver: the response's
  // status is "ok", "insufficient_input" or "error", with reasons in errors[].
  Json analyze(const IQBlock& block) const;
  std::string analyze_json(const IQBlock& block, int indent = -1) const;

  std::vector<std::string> receivers() const;  // enabled receivers, in priority order
  Json info() const;                            // the response's "engine" block

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

// Every receiver built into this binary, in priority order.
std::vector<std::string> available_receivers();

// Sample rate and duration each receiver needs, and the minimum / recommended block size
// (schema/requirements.json, compiled in).
const Json& requirements();

// Accumulates contiguous blocks of one tuning (2048, 4096, ... samples) into analysable
// blocks of target_samples. A change of centre frequency, sample rate or bandwidth, or a
// gap in start times, closes the current block.
class StreamAnalyzer {
 public:
  explicit StreamAnalyzer(const Engine& engine, size_t target_samples = 0 /* 0 = recommended */,
                          size_t max_gap_samples = 1);
  ~StreamAnalyzer();
  // Add one block; returns the responses of the blocks it completed (usually none or one).
  std::vector<Json> push(const IQBlock& block);
  // Analyse what is buffered if it is long enough (else drop it) and forget the tuning.
  std::optional<Json> flush();
  size_t buffered() const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace si
