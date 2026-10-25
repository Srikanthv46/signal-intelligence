// Data records of the detection stage that receivers may read (python:
// detect/stft.Spectrogram and detect/pipeline.DetectResult). The algorithms that fill
// them live in detect/; receivers only see these plain records.
#pragma once

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include "core/model.hpp"
#include "core/types.hpp"

namespace si {

struct Spectrogram {
  fmat P;              // (T, F) power per bin, float32; white noise of variance s2 -> E[P] = s2
  dvec f;              // (F,) absolute bin-centre frequency, Hz
  dvec t;              // (T,) frame start time, s
  size_t nfft = 0, hop = 0;
  double fs = 0, fc = 0;

  double df() const { return fs / static_cast<double>(nfft); }
  double dt() const { return static_cast<double>(hop) / fs; }
  double frame_len() const { return static_cast<double>(nfft) / fs; }
  // nearest bin / frame index, clipped (python Spectrogram.f_index / t_index)
  int64_t f_index(double f_hz) const {
    double v = std::nearbyint((f_hz - f[0]) / df());
    return std::clamp<int64_t>(static_cast<int64_t>(v), 0, static_cast<int64_t>(f.size()) - 1);
  }
  int64_t t_index(double t_s) const {
    double v = std::nearbyint(t_s / dt());
    return std::clamp<int64_t>(static_cast<int64_t>(v), 0, static_cast<int64_t>(t.size()) - 1);
  }
};

struct DetectResult {
  std::vector<Signal> signals;     // energy detections (with features once measured)
  Spectrogram coarse;
  dvec floor_coarse_db;
  Spectrogram fine;
  dvec floor_fine_db;
  std::pair<double, double> usable_hz{0, 0};
  Json info = Json::object();
};

}  // namespace si
