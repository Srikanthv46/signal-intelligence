// Emitter association (python: receivers/emitters/association.py): which bursts (Emission
// records from any receiver) come from one radio, which radios form one drone system.
//
//   measure_fingerprints(x, fs, fc, emissions)   refines 15 kHz-OFDM timing/carrier features in place
//   same_emitter_score(a, b, models)             pairwise log-likelihood ratio + evidence cues
//   associate(emissions, models)                 -> emitters (average-linkage on the LLR matrix)
//
// Depends on core/ only.
#pragma once

#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "core/model.hpp"
#include "core/spectral.hpp"
#include "core/types.hpp"

namespace si::rx::emitters {

// burst starts of ONE radio lie on t0 + n * frame_us + offset (its own frame clock)
struct Lattice {
  double frame_us = 0;
  dvec offsets_us;
  double weight = 1.0;
};

// Per (family, label) association model (python MODELS / DEFAULT_MODEL dicts; defaults = DEFAULT_MODEL).
struct Model {
  bool variable_length = false;
  std::vector<Lattice> lattices;
  double jitter_us = 2.0;
  double chan_delay_us_per_mhz = 0.05;
  std::optional<std::pair<double, double>> raster;   // (offset_hz, step_hz)
  double ppm_h1 = 0.1;
  double ppm_drift_per_s = 0.1;
  double ppm_pop = 5.0;
  double snr_h1_db = 5.0;
  bool single_tx = true;
};
using ModelKey = std::pair<std::string, std::string>;   // (family, label)
using Models = std::map<ModelKey, Model>;

const Models& MODELS();

// One cue of same_emitter_score: {"cue", "llr", "detail", extra keys...} as python's dict.
struct Cue {
  std::string cue;
  double llr = 0;
  std::string detail;
  Json extra = Json::object();     // residual_us / sigma_us (timing), dppm (oscillator), p_same (total)
  Json to_json() const;
};

struct Emitter {
  std::string emitter_id, family, label;
  std::vector<int64_t> members;          // indices into the input emission list (time order)
  dvec frequencies;                      // distinct channel centres (MHz)
  Json hop_pattern = Json::object();
  std::optional<double> frame_period_ms;
  Json fingerprint = Json::object();
  std::optional<double> confidence;
  std::vector<std::string> evidence;
  Json separation = Json::array();
  Json linked = Json::array();
  Json to_dict() const;                  // python dataclasses.asdict (field order)
};

// python measure_fingerprints(x, fs, fc, emissions, ch=None, notch_hz=()) -> info
Json measure_fingerprints(const cfvec& x, double fs, double fc, std::vector<Emission>& emissions,
                          const SpectralChannelizer* ch = nullptr,
                          const std::vector<std::pair<double, double>>& notch_hz = {});

// python same_emitter_score(a, b, models=None, prior_same=0.5) -> (llr, evidence)
std::pair<double, std::vector<Cue>> same_emitter_score(const Emission& a, const Emission& b,
                                                       const Models* models = nullptr, double prior_same = 0.5);

struct Association {
  std::vector<Emitter> emitters;
  dmat M;                                                // pairwise LLR matrix
  std::map<std::pair<int64_t, int64_t>, std::vector<Cue>> EV;
};
// python associate(emissions, fc=None, models=None, prior_same=0.5, return_matrix=True)
Association associate_full(const std::vector<Emission>& emissions, const Models* models = nullptr,
                           double prior_same = 0.5);
inline std::vector<Emitter> associate(const std::vector<Emission>& emissions, const Models* models = nullptr,
                                      double prior_same = 0.5) {
  return associate_full(emissions, models, prior_same).emitters;
}

}  // namespace si::rx::emitters
