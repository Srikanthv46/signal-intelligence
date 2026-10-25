#include "receivers/rc_link/label_check.hpp"

#include <algorithm>
#include <cmath>

#include "core/npcompat.hpp"
#include "core/pyfmt.hpp"

namespace si::rx::rc_link {

std::optional<Json> check_frequency_label(const std::vector<Emission>& emissions) {
  dvec offs;
  for (const auto& e : emissions) {
    if (e.label != "ble_adv" || !e.decoded.is_object()) continue;
    const Json& d = e.decoded;
    if (!(d.contains("crc_ok") && d["crc_ok"].is_boolean() && d["crc_ok"].get<bool>())) continue;
    if (!d.contains("whitening_channel") || !d["whitening_channel"].is_number_integer()) continue;
    const int w = d["whitening_channel"].get<int>();
    const double adv = w == 37 ? 2402.0 : w == 38 ? 2426.0 : w == 39 ? 2480.0 : 0.0;
    if (adv == 0.0) continue;
    offs.push_back(e.center_mhz - adv);
  }
  if (offs.empty()) return std::nullopt;
  const double off = np::median(offs);
  const double ptp = *std::max_element(offs.begin(), offs.end()) - *std::min_element(offs.begin(), offs.end());
  const bool suspect = std::fabs(off) > 0.5;
  Json out = Json::object();
  out["offset_mhz"] = py_round(off, 3);
  out["n_adverts"] = offs.size();
  out["spread_mhz"] = py_round(ptp, 3);
  out["suspect"] = suspect;
  out["note"] = suspect ? fmt("capture frequency label looks wrong by %+.1f MHz: rerun with --fc <labelled - %.1f MHz>", off, off)
                        : std::string("frequency label consistent with BLE adverts");
  return out;
}

}  // namespace si::rx::rc_link
