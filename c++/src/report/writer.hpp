// Per-block report (python: report/writer.py): the user schema + provenance blocks.
// The production response (api/response) is built from it.
#pragma once

#include "analyze.hpp"

namespace si {

// python report.writer.build(a): signals -> dicts, weak detections split off,
// plugins' finalize_report hooks applied.
Json build_report(const Analysis& a);

}  // namespace si
