// Python-compatible text formatting, for report strings (evidence, steps, reasons) that
// must read exactly as the Python reference writes them.
#pragma once

#include <string>

#include "core/types.hpp"

namespace si {

// repr(float): shortest round-trip digits; fixed notation for 1e-4 <= |x| < 1e16 (always
// with a fractional part, "100.0"), otherwise scientific "1e-05", "1.5e+16"; nan, inf.
std::string py_float_repr(double x);

// repr() of a JSON value as Python would print the equivalent object:
// {'a': 1, 'b': [1.5, 'x']}, True/False/None, strings in single quotes unless they
// contain a single quote and no double quote.
std::string py_repr(const Json& v);

// str() of a JSON value: strings unquoted, everything else as py_repr.
std::string py_str(const Json& v);

// printf-style formatting into a std::string (for f"{x:.1f}", f"{x:+.2f}", f"{x:g}" ...).
// Python's 'g' and printf's %g agree for the values used here; f"{x:,}" etc. are not provided.
std::string fmt(const char* format, ...)
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((format(printf, 1, 2)))
#endif
    ;

}  // namespace si
