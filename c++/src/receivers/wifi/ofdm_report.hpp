// The report side of the 802.11 OFDM plugins (python plugin_ofdm.py minus the pipeline glue in
// plugin.cpp): PPDU records with decoded payloads, VHT Partial-AID attribution, and the decode
// results of wifi_ofdm (per channel / per signal) and wifi_ofdm_scaled (per link).
#pragma once

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "core/types.hpp"
#include "receivers/wifi/dsss_report.hpp"
#include "receivers/wifi/ofdm_rx.hpp"
#include "receivers/wifi/ofdm_scaled.hpp"

namespace si::rx::wifi::ofdm_report {

using Summary = dsss::Summary;

// python _records(pkts, cap, clock_hz): payloads decoded on up to nthreads threads.
// scaled: the packets come from ofdm_scaled (bandwidths are floats there).
std::vector<Json> records(const std::vector<ofdm::Packet>& pkts, double cap_fs, double clock_hz, bool scaled,
                          int nthreads = 0);
// python _attribute_vht(recs)
void attribute_vht(std::vector<Json>& recs);

// python WifiOfdmDecoder.summarise(mine, n_all, n_rej) (mutates the records: links.attribute)
// all: every record (BSS colour -> BSSID); block_s > 0 adds the field.py link profiles (PLAN W1)
Summary summarise_ofdm(const std::vector<Json*>& mine, size_t n_all, size_t n_rej,
                       const std::vector<Json>* all = nullptr, double block_s = 0);
// python WifiOfdmScaledDecoder.summarise(ln, recs)
Summary summarise_scaled(const scaled::Link& ln, const std::vector<Json*>& recs, double block_s = 0);
// python plugin_ofdm.field_summary(frames, block_s): wifi_links (top 10) + drone_links
Json field_summary(const std::vector<Json>& frames, double block_s);
// python plugin_ofdm.drone_label(fields): (label suffix, evidence suffix)
std::pair<std::string, std::string> drone_label(const Json& fields);

}  // namespace si::rx::wifi::ofdm_report
