// Capture centre-frequency label check from CRC-verified BLE advertising packets
// (python: receivers/rc_link/label_check.py). A CRC-valid advert only passes when
// de-whitened with its true channel index (37/38/39 -> 2402/2426/2480 MHz), so its
// labelled frequency minus the true advertising-channel frequency is the label error.
#pragma once

#include <optional>
#include <vector>

#include "core/model.hpp"
#include "core/types.hpp"

namespace si::rx::rc_link {

// python check_frequency_label(emissions): {"offset_mhz", "n_adverts", "spread_mhz", "suspect",
// "note"} or nullopt without CRC-verified adverts
std::optional<Json> check_frequency_label(const std::vector<Emission>& emissions);

}  // namespace si::rx::rc_link
