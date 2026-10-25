// Device tracking over every decoded 802.11 frame of a capture (python: receivers/devices/plugin.py,
// the part that is not pipeline glue).
//
//   frames_of(fields, source)  one signal's decoded.fields.ppdus -> normalised frames
//   device_report(signals)     report extras {devices, findings, device_summary[, device_groups]};
//                              sets decoded.fields.device_macs on the signals the devices were heard in
//   finalize(rep)              points devices and findings at the report's signal ids
#pragma once

#include <string>
#include <vector>

#include "core/model.hpp"
#include "receivers/devices/threats.hpp"
#include "receivers/devices/tracker.hpp"

namespace si::rx::devices {

std::vector<Frame> frames_of(const Json& fields, const std::string& source);

Json device_report(std::vector<Signal>& signals);

void finalize(Json& rep);

}  // namespace si::rx::devices
