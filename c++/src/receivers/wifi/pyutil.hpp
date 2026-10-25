// Python semantics the Wi-Fi MAC dissector relies on, now in core/pybytes.hpp (shared with core/remote_id);
// re-exported into the wifi receiver's detail namespace so its code is unchanged.
#pragma once

#include "core/pybytes.hpp"

namespace si::rx::wifi::detail {
using si::pybytes::IndexError;
using si::pybytes::Bytes;
using si::pybytes::u16le;
using si::pybytes::be16;
using si::pybytes::le_uint;
using si::pybytes::be_uint;
using si::pybytes::hex;
using si::pybytes::mac_str;
using si::pybytes::decode_utf8;
using si::pybytes::decode_ascii;
using si::pybytes::py_strip;
using si::pybytes::py_strip_chars;
using si::pybytes::py_split_ws;
using si::pybytes::ascii_lower;
using si::pybytes::ascii_upper;
using si::pybytes::head_codepoints;
using si::pybytes::len_codepoints;
using si::pybytes::truthy;
using si::pybytes::jget;
using si::pybytes::pstr;
using si::pybytes::json_key;
}  // namespace si::rx::wifi::detail
