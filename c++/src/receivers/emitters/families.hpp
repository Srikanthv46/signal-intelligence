// Which emission families are drone links, and the association model for each detector
// label (python: receivers/emitters/families.py).
#pragma once

#include <map>
#include <string>

#include "receivers/emitters/association.hpp"

namespace si::rx::emitters {

// what an emission family means for the operator
inline const std::map<std::string, std::string>& DRONE_FAMILIES() {
  static const std::map<std::string, std::string> m = {
      {"dji_ocusync", "DJI drone link"},       {"dji_lightbridge", "DJI drone link"},
      {"dji_droneid", "DJI DroneID"},          {"herelink", "Herelink drone link"},
      {"analog_vtx", "analog FPV video"},      {"digital_vtx", "unknown digital video/datalink"},
      {"rc_fhss", "RC control link"},          {"telemetry", "telemetry link"}};
  return m;
}

// detector labels -> keys of MODELS
inline const std::map<ModelKey, ModelKey>& LABEL_ALIASES() {
  static const std::map<ModelKey, ModelKey> m = {{{"dji_ocusync", "uplink"}, {"dji_ocusync", "ocusync_uplink"}}};
  return m;
}

// python _assoc_models(): MODELS plus the aliases (setdefault)
inline const Models& assoc_models() {
  static const Models m = [] {
    Models out = MODELS();
    for (const auto& [k, v] : LABEL_ALIASES()) out.emplace(k, MODELS().at(v));
    return out;
  }();
  return m;
}

}  // namespace si::rx::emitters
