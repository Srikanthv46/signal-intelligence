// The shared contract in <repo>/schema, compiled in (python: api/contract.py):
// technology codes and families, identifier keys, receiver requirements.
#pragma once

#include <string>
#include <utility>
#include <vector>

#include "core/model.hpp"
#include "core/types.hpp"

namespace si::contract {

extern const char* const kSchemaVersion;          // response schema version, e.g. "1.0.0"
using si::kArtifactFlags;                         // core/model.hpp (python model.ARTIFACT_FLAGS, re-exported)

const Json& technologies();                       // schema/technologies.json
const Json& requirements();                       // schema/requirements.json
const Json& response_schema();                    // schema/response.schema.json
const Json& capabilities();                       // schema/capabilities.json
const Json& descriptions();                       // schema/descriptions.json

// proven / research validation tag of a signal (python validation()): the technology's level from
// schema/capabilities.json ('research' when missing) plus the research-level decoded fields present.
Json validation(const std::string& tech, const std::string& decoder, const Json& fields);

// (technology code, family) of a signal from its label and the receiver that defined it.
std::pair<std::string, std::string> technology(const std::string& label, const std::string& decoder);

const std::vector<std::string>& identifier_keys();

// python description(label, decoder, flags): plain-language description (schema/descriptions.json); a signal
// flagged as a receiver artifact is described as that artifact; else the technology family's description
std::string description(const std::string& label, const std::string& decoder, const Json& flags = Json::array());

}  // namespace si::contract
