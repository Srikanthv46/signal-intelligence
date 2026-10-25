// Model-level identification tables and the common `model` block (python core/drone_models.py; research, PLAN M0-M2;
// research/drone_transmitters/11_model_level_identification.md).
#pragma once

#include <string>
#include <vector>

#include "core/types.hpp"

namespace si::drone_models {

// python dji_product(code): model name or "" (unknown codes are never mapped to the nearest name)
std::string dji_product(int64_t code);
// python dji_product_source(code)
std::string dji_product_source(int64_t code);
// python serial_model(serial): {models, verified, doc, prefix} or null
Json serial_model(const std::string& serial);
// python serial_manufacturer(serial): {code, maker, kind} or null
Json serial_manufacturer(const std::string& serial);
// python model_block(family, candidates, evidence, basis, unresolved_reason); empty reason = python None
Json model_block(const std::string& family, const std::vector<std::string>& candidates,
                 const std::vector<std::string>& evidence, const std::string& basis = "inferred",
                 const std::string& unresolved_reason = "");
// python remote_id_model(rid)
Json remote_id_model(const Json& rid);

}  // namespace si::drone_models
