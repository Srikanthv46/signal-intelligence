#include "decode/registry.hpp"

#include <algorithm>
#include <stdexcept>

namespace si {

const std::vector<std::string> kReceiverOrder = {"wifi",        "rc_link", "dji", "analog_vtx", "voice_fm", "pmr_digital",
                                                 "digital_vtx", "lte",     "nr",  "gsm",        "adsb",     "emitters"};

std::vector<std::string> receiver_names() {
  std::vector<std::string> found;
  for (const auto& e : built_receivers()) found.emplace_back(e.name);
  std::sort(found.begin(), found.end());
  std::vector<std::string> out;
  for (const auto& n : kReceiverOrder)
    if (std::find(found.begin(), found.end(), n) != found.end()) out.push_back(n);
  for (const auto& n : found)
    if (std::find(kReceiverOrder.begin(), kReceiverOrder.end(), n) == kReceiverOrder.end()) out.push_back(n);
  return out;
}

PluginList plugins(const std::vector<std::string>& enable, const std::vector<std::string>& disable) {
  const auto names = receiver_names();
  std::vector<std::string> unknown;
  for (const auto* lst : {&enable, &disable})
    for (const auto& n : *lst)
      if (std::find(names.begin(), names.end(), n) == names.end()) unknown.push_back(n);
  if (!unknown.empty()) {
    std::string msg = "unknown receiver(s):";
    for (const auto& u : unknown) msg += " " + u;
    msg += "; available:";
    for (const auto& n : names) msg += " " + n;
    throw std::invalid_argument(msg);
  }
  PluginList out;
  for (const auto& n : names) {
    if (!enable.empty() && std::find(enable.begin(), enable.end(), n) == enable.end()) continue;
    if (std::find(disable.begin(), disable.end(), n) != disable.end()) continue;
    for (const auto& e : built_receivers())
      if (n == e.name) e.register_plugins(out);
  }
  return out;
}

}  // namespace si
