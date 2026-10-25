// Decoded 802.11 frames -> one record per transmitting device (python: receivers/devices/tracker.py).
//
//   tracks(frames)           -> device dicts, strongest first
//   link_randomized(devices) -> randomised MACs that are probably one device
//
// Devices and groups are Json objects with the keys, order and value types Python writes.
#pragma once

#include <string>
#include <vector>

#include "core/types.hpp"

namespace si::rx::devices {

// One decoded frame, normalised (python tracker.frame()). `mac` is the dissected 802.11 frame.
struct Frame {
  Json mac = Json::object();
  double t0 = 0.0;
  double dur_s = 0.0;
  Json snr_db = nullptr;
  Json cfo_hz = nullptr;
  Json channel = nullptr;
  Json phy = nullptr;
  Json rate_mbps = nullptr;
  Json length_bytes = nullptr;
  std::string source;
};

Frame frame(const Json& mac, double t0, double dur_s, Json snr_db, Json cfo_hz, Json channel, Json phy,
            Json rate_mbps, std::string source, Json length_bytes = nullptr);

Json tracks(const std::vector<Frame>& frames);

Json link_randomized(const Json& devices, double cfo_tol_hz = 2500.0);

}  // namespace si::rx::devices
