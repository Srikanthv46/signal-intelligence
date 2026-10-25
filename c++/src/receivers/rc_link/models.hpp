// Model level for RC and telemetry links (python receivers/rc_link/models.py; research, PLAN M5).
#pragma once

#include <vector>

#include "core/model.hpp"

namespace si::rx::rc_link::models {

// python model(em): the model block or null
Json model(const std::vector<Emission>& em);

}  // namespace si::rx::rc_link::models
