// A capture (one IQ block plus its metadata) and the SigMF reader
// (python: spectrum_intel/io/sigmf.py).
#pragma once

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "core/timeutil.hpp"
#include "core/types.hpp"

namespace si {

struct Capture {
  std::string name = "live";               // file stem or a caller-chosen name
  std::string meta_path;                   // .sigmf-meta, "" for an in-memory block
  double fc = 0;                           // centre frequency, Hz (trusted as given)
  double fs = 0;                           // sample rate, Hz
  std::optional<UtcTime> t0;               // time of the first sample
  std::optional<double> bandwidth_hz;      // analog front-end bandwidth, if known
  std::optional<double> gain_db;
  std::string datatype = "cf32_le";
  Json meta = Json::object();
  std::shared_ptr<const cfvec> samples;    // complex64, DC kept (LO leakage is reported)

  const cfvec& iq() const { return *samples; }
  size_t n() const { return samples ? samples->size() : 0; }
  double duration() const { return fs > 0 ? static_cast<double>(n()) / fs : 0.0; }

  static Capture from_iq(cfvec iq, double fs, double fc, std::optional<UtcTime> t0 = std::nullopt,
                         std::string name = "live", std::optional<double> gain_db = std::nullopt,
                         std::optional<double> bandwidth_hz = std::nullopt);
};

// Raw samples in a SigMF datatype -> complex64 (cf32_le, cf64_le, ci16_le, ci8, cu8; values
// not rescaled, cu8 centred on 127.5).
cfvec convert_iq(const void* data, size_t n_samples, const std::string& datatype);

// A .sigmf-meta, its .sigmf-data, or a folder holding exactly one .sigmf-meta.
Capture load_sigmf(const std::string& path);

// Every .sigmf-meta under folder, sorted by centre frequency then time.
std::vector<std::string> find_captures(const std::string& folder);

}  // namespace si
