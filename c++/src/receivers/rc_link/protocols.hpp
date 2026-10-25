// Research catalog of RC / telemetry / identification FSK protocols (python receivers/rc_link/protocols.py).
#pragma once

#include <optional>

#include "core/types.hpp"

namespace si::rx::rc_link::protocols {

Json candidates(double sym_rate, std::optional<double> dev_hz, std::optional<double> f_mhz,
                std::optional<double> interval_ms, size_t top = 3);

}  // namespace si::rx::rc_link::protocols
