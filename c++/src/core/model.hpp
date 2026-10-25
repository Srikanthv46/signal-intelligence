// Records passed between pipeline stages (python: spectrum_intel/model.py and
// core/emission.py). Free-form dicts of the Python version (features, detection,
// decoded fields) are Json objects; Python's private objects kept under "_" keys of
// features (_emissions, _video_frame, _audio) are typed members here.
#pragma once

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "core/types.hpp"

namespace si {

// One connected time-frequency region above the noise floor (absolute units).
struct Box {
  double t0 = 0, t1 = 0;          // s from capture start
  double f0 = 0, f1 = 0;          // Hz, absolute RF
  double peak_snr_db = 0;
  double mean_snr_db = 0;         // power-weighted mean SNR inside the region
  std::string res = "coarse";     // "coarse" | "fine" | "receiver"
  int64_t n_cells = 0;
  std::string kind = "burst";     // "burst" | "continuous"

  double fc() const { return 0.5 * (f0 + f1); }
  double bw() const { return f1 - f0; }
  double dur() const { return t1 - t0; }
};

// One burst / continuous segment found by a receiver (python core.emission.Emission).
struct Emission {
  std::string family;             // dji_ocusync | herelink | analog_vtx | lora | bluetooth | ...
  std::string label;              // finer label
  double start_us = 0, end_us = 0;
  double f_lo_mhz = 0, f_hi_mhz = 0, center_mhz = 0;
  std::optional<double> snr_db;
  double confidence = 0.5;
  Json decoded = nullptr;         // payload if demodulated (null = none)
  Json features = Json::object(); // physical-layer fingerprint (cfo_hz, cfo_ppm, symbol_us, ...)
};

// Decode outcome attached to a signal (python: Signal.decoded dict).
struct Decoded {
  std::string status = "not_possible";      // success | partial | not_possible
  std::string reason = "no decoder run";
  std::optional<std::string> method;        // None until a decoder claims it
  std::vector<std::string> steps;
  Json fields = Json::object();
};

// A grayscale picture (analog video field), row-major.
struct Image {
  size_t rows = 0, cols = 0;
  std::vector<float> pixels;
};

// One emitter: a set of bursts sharing a frequency span.
// flags that make a signal a receiver artifact, not an emission (response category "artifact")
inline constexpr const char* kArtifactFlags[4] = {"lo_leakage", "iq_image", "receiver_spur", "lo_harmonic_response"};

struct Signal {
  std::vector<Box> bursts;
  double f0 = 0, f1 = 0;
  double snr_db = 0;
  std::string label = "unknown";
  double confidence = 0;
  std::string modulation = "unknown";
  std::string evidence;
  Json features = Json::object();     // measured features (serialised in the report)
  Json detection = Json::object();    // provenance: how it was found
  Decoded decoded;
  std::vector<std::string> flags;     // artifact / quality flags
  std::vector<std::string> classify_trace;

  // typed side channels (python: features["_emissions"], ["_video_frame"], ["_audio"])
  std::vector<Emission> emissions;
  std::shared_ptr<Image> video_frame;
  std::shared_ptr<std::vector<float>> audio;
  double audio_fs = 0;

  double t0() const;
  double t1() const;
  double fc() const { return 0.5 * (f0 + f1); }
  double bw() const { return f1 - f0; }
  bool has_flag(const std::string& f) const;
};

}  // namespace si
