// One block end to end: detect -> features -> receivers -> classify (python: analyze.py).
#pragma once

#include <optional>
#include <vector>

#include "core/model.hpp"
#include "core/spectrogram.hpp"
#include "decode/plugin.hpp"
#include "io/sigmf.hpp"

namespace si {

struct Analysis {
  Capture cap;                         // as analysed (fc corrected when a plugin proved the label wrong)
  DetectResult det;
  std::vector<Signal> signals;         // final signals, sorted by (f0, t0)
  Json quality = Json::object();       // detect/artifacts clipping() result
  std::optional<Json> label_check;     // centre-frequency label check by a plugin (BLE adverts)
  Json timings = Json::object();       // python Analysis.timings: detect_s, features_s, receiver_s{}, *_error, ...
  Json extras = Json::object();        // report-level additions from plugins (emitters, devices, ...)
  std::vector<const Decoder*> plugins; // the plugins that ran (finalize_report hooks)
};

// Analyse one capture with the given plugins (python analyze.analyze(cap, None, plugins, ppm)):
// if a plugin's decodes prove the centre-frequency label wrong, analyse again at the
// corrected frequency and record both in label_check. ppm: the receiver's frequency error
// {"ppm": x, "source": "..."} when the caller knows it; it goes into the capture cache as
// "sweep_ppm" (python: ctx.cache["sweep_ppm"]), where the LTE/NR/GSM/voice receivers read it.
Analysis analyze(const Capture& cap, const PluginList& plugins, const std::optional<Json>& ppm = std::nullopt);

}  // namespace si
