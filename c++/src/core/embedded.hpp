// Text files compiled into the library (cmake/EmbedFiles.cmake): the shared contract in
// <repo>/schema (technologies.json, requirements.json, response.schema.json).
#pragma once

#include <string>

namespace si {
// File contents by base name, nullptr if not embedded.
const char* embedded_file(const std::string& name);
}  // namespace si
