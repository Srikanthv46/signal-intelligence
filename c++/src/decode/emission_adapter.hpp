// Adapter for receivers whose detector returns Emission records (one per burst)
// (python: decode/emission_adapter.py). EmissionReceiver groups the bursts into
// pipeline signals (one per link / channel), summarises them into a decode result and
// keeps the Emission records in Signal::emissions for later plugins (emitter
// association). A plugin subclasses it and implements run_detector(); keep(),
// group_key(), summarise(), finish() and in_range() refine the behaviour.
#pragma once

#include <set>
#include <string>
#include <vector>

#include "decode/plugin.hpp"

namespace si {

bool in_ranges(double lo, double hi, const std::vector<std::pair<double, double>>& ranges);

// Emission features as the summary shows them: scalars only (lists/dicts dropped, except sync_syms).
Json scalar_features(const Json& features);

class EmissionReceiver : public Decoder {
 public:
  bool capture_level() const override { return true; }

  virtual std::string key() const = 0;                     // ctx.cache key of the capture-wide run
  virtual std::set<std::string> bands() const { return {}; }
  virtual std::string method() const = 0;
  virtual std::vector<std::string> steps() const = 0;
  virtual bool owns_band() const { return false; }         // everything inside the signal's band belongs to it

  virtual bool in_range(double lo, double hi) const;
  virtual std::vector<Emission> run_detector(DecodeContext& ctx) const = 0;
  // the capture-wide run, cached in ctx.cache under key() (empty when out of range)
  const std::vector<Emission>& emissions(DecodeContext& ctx) const;

  virtual bool keep(const Emission& e) const;
  // grouping key (python: a tuple; here any string that is equal for the same tuple)
  virtual std::string group_key(const Emission& e) const;
  virtual DecodeResult summarise(const std::vector<Emission>& em, const Signal* sig) const;
  virtual void finish(Signal& sig, const std::vector<Emission>& em, DecodeContext& ctx) const {}
  virtual std::string pretty(const std::string& top, const std::vector<Emission>& em) const { return top; }
  virtual std::string modulation(const std::vector<Emission>& em) const;

  bool can_decode(const Signal& sig, DecodeContext& ctx) const override;
  std::vector<Signal> capture_signals(DecodeContext& ctx) const override;
  DecodeResult decode(const Signal& sig, DecodeContext& ctx) const override;
};

}  // namespace si
