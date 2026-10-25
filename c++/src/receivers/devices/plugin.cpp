// Pipeline plugin: device tracking and findings across every decoded 802.11 frame in a capture
// (python: receivers/devices/plugin.py). Runs after all receivers, reads their decoded frames off
// the generic report contract (decoded.fields.ppdus) and imports no receiver.
#include <memory>

#include "decode/plugin.hpp"
#include "receivers/devices/devices.hpp"

namespace si::rx::devices {
namespace {

class DeviceTracking : public Decoder {
 public:
  std::string name() const override { return "devices"; }
  std::string package() const override { return "devices"; }
  std::string protocol() const override { return "device tracking and findings (802.11)"; }
  bool capture_level() const override { return true; }
  std::vector<Signal> capture_signals(DecodeContext&) const override { return {}; }
  Json postprocess(DecodeContext&, std::vector<Signal>& signals) const override { return device_report(signals); }
  void finalize_report(Json& rep) const override { finalize(rep); }
};

}  // namespace

void register_plugins(PluginList& out) { out.push_back(std::make_unique<DeviceTracking>()); }

}  // namespace si::rx::devices
