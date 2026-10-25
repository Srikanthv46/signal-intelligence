// Plugin interface between the pipeline and the receivers (python: decode/base.py).
//
// A receiver (src/receivers/<name>/) keeps its DSP in files that include only core/
// and io/ headers, plus a plugin.cpp that adapts it to the pipeline:
//
//   class LteDecoder : public si::Decoder { ... };
//   namespace si::rx::lte {
//   void register_plugins(si::PluginList& out) { out.push_back(std::make_unique<LteDecoder>()); }
//   }
//
// Every hook is optional (the defaults do nothing):
//   capture_signals(ctx)            receiver-defined signals; energy detections they explain are absorbed
//   arbitrate(signals, ctx)         after all capture_signals: drop/merge other receivers' signals
//   frequency_label_check(signals)  {"offset_mhz", "suspect", ...} if the plugin's decodes prove the
//                                   capture's centre-frequency label wrong (the pipeline re-analyses)
//   postprocess(ctx, signals)       report-level extras, e.g. {"emitters": [...]}
//   finalize_report(rep)            edit the report after signal ids are final
//   can_decode / decode             per-signal decoding of energy detections
//
// capture_signals() of different receiver packages run concurrently, each package on its
// own copy of the cache; a plugin that reads what another package's capture_signals() put
// in the cache names that package in after() (nr and gsm: after lte, for the receiver
// ppm). Every other hook runs sequentially in registry order.
#pragma once

#include <any>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <vector>

#include "core/dsp.hpp"
#include "core/model.hpp"
#include "core/spectrogram.hpp"
#include "io/sigmf.hpp"

namespace si {

// Capture-wide work shared between hooks of one run (python ctx.cache dict).
class Cache {
 public:
  bool contains(const std::string& k) const { return m_.count(k) > 0; }
  template <class T>
  void set(const std::string& k, T v) { m_[k] = std::move(v); }
  // nullptr if absent or of another type
  template <class T>
  const T* get(const std::string& k) const {
    auto it = m_.find(k);
    return it == m_.end() ? nullptr : std::any_cast<T>(&it->second);
  }
  template <class T>
  T* get_mut(const std::string& k) {
    auto it = m_.find(k);
    return it == m_.end() ? nullptr : std::any_cast<T>(&it->second);
  }
  const std::map<std::string, std::any>& items() const { return m_; }
  std::map<std::string, std::any>& items() { return m_; }

 private:
  std::map<std::string, std::any> m_;
};

struct DecodeContext {
  const Capture& cap;
  const Channelizer& ch;          // one FFT of the whole capture (core/dsp.hpp)
  const DetectResult* det;        // energy detections, spectrograms, floors
  Cache cache;
};

// What a decoder returns (python decode.base.result()).
struct DecodeResult {
  std::string status = "not_possible";      // success | partial | not_possible
  std::string method;
  std::vector<std::string> steps;
  Json fields = Json::object();
  std::string reason;
  std::optional<std::string> label;          // overrides the classifier when set
  std::optional<double> confidence;
  std::optional<std::string> modulation;
  std::optional<std::string> evidence;
  std::optional<std::pair<double, double>> span;  // refined occupied band (f_low, f_high)
};

DecodeResult result(std::string status, std::string method, std::vector<std::string> steps,
                    Json fields = Json::object(), std::string reason = "",
                    std::optional<std::string> label = std::nullopt, std::optional<double> confidence = std::nullopt,
                    std::optional<std::string> modulation = std::nullopt,
                    std::optional<std::string> evidence = std::nullopt,
                    std::optional<std::pair<double, double>> span = std::nullopt);

class Decoder {
 public:
  virtual ~Decoder() = default;
  virtual std::string name() const = 0;          // plugin name (e.g. "wifi_dsss", "nr_ssb")
  virtual std::string package() const = 0;       // receiver folder (e.g. "wifi", "nr")
  virtual std::string protocol() const { return ""; }
  virtual bool capture_level() const { return false; }
  virtual std::vector<std::string> after() const { return {}; }

  virtual bool can_decode(const Signal&, DecodeContext&) const { return false; }
  virtual DecodeResult decode(const Signal&, DecodeContext&) const { return {}; }
  virtual std::vector<Signal> capture_signals(DecodeContext&) const { return {}; }
  virtual std::vector<Signal> arbitrate(std::vector<Signal> signals, DecodeContext&) const { return signals; }
  virtual std::optional<Json> frequency_label_check(const std::vector<Signal>&) const { return std::nullopt; }
  virtual Json postprocess(DecodeContext&, std::vector<Signal>&) const { return Json::object(); }
  virtual void finalize_report(Json& /*rep*/) const {}
};

using PluginList = std::vector<std::unique_ptr<Decoder>>;

// A signal defined by a receiver's own decodes (python decode.base.receiver_signal).
Signal receiver_signal(const std::vector<std::pair<double, double>>& bursts, double f0, double f1,
                       const DecodeResult& res, std::optional<double> snr_db, const std::string& decoder);

// Does sig overlap the (t0, t1, f0, f1) box (python decode.base.overlaps)?
bool overlaps(const Signal& sig, double t0, double t1, double f0, double f1, double f_frac = 0.3);

}  // namespace si
