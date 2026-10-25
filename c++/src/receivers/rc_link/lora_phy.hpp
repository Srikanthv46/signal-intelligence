// LoRa payload decoding from hard symbol values, explicit-header packets
// (python: receivers/rc_link/lora_phy.py; the Semtech SX127x chain as gr-lora_sdr
// implements it: symbol - 1, reduced rate in the header block / with LDRO, Gray,
// diagonal deinterleaving, Hamming syndrome correction, 5-bit header checksum,
// de-whitening, CRC-16/CCITT). Implicit-header packets (ExpressLRS) are not decoded.
#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "core/types.hpp"

namespace si::rx::rc_link::lora_phy {

// LoRa whitening bytes 0xFF, 0xFE, 0xFC, ... (python whitening_sequence)
std::vector<uint8_t> whitening_sequence(size_t n);
// CRC-16/CCITT, polynomial 0x1021, init 0, MSB first (python crc16)
uint16_t crc16(const uint8_t* data, size_t n);
int header_checksum(int n0, int n1, int n2);
// payload symbols after the SFD of an explicit-header packet (python n_symbols)
int64_t n_symbols(int sf, int length, int cr, int crc, bool ldro);
// low-data-rate optimisation on when a symbol lasts more than 16 ms (python default_ldro)
bool default_ldro(int sf, double bw);

// python decode(raw_syms, sf, bw, ldro=None): the result dict (keys in Python's order), or
// nullopt when the header block is missing or its checksum fails.
std::optional<Json> decode(const std::vector<int64_t>& raw_syms, int sf, double bw = 125e3,
                           std::optional<bool> ldro = std::nullopt);

// SF5 (SX126x, measured on the bench): 2 extra symbols after the SFD, first block CR 4/8 at full SF.
constexpr int SF5_EXTRA_SYMS = 2;
// python n_symbols_sf5(length, cr, crc=True, sf=5)
int64_t n_symbols_sf5(int length, int cr, bool crc = true, int sf = 5);
// python decode_sf5_implicit(raw_syms, n_pay, sf=5): implicit-header SF5 packet, length and CR
// found by the payload CRC-16 within n_pay -8..+8 symbols; nullopt when nothing (or more than one
// payload) passes.
std::optional<Json> decode_sf5_implicit(const std::vector<int64_t>& raw_syms, int64_t n_pay, int sf = 5);

// python decode_implicit(raw_syms, sf, cr, length, ldro): implicit-header payload of known length / CR (SF >= 7)
struct Implicit {
  std::vector<uint8_t> payload;
  std::string cr;
  int fec_corrected = 0, fec_detected_errors = 0;
  int64_t symbols_used = 0;
};
std::optional<Implicit> decode_implicit(const std::vector<int64_t>& raw_syms, int sf, int cr, int length, bool ldro = false);

// python sync_word(sync_vals): the sync byte from the two sync symbols, or nullopt
std::optional<int> sync_word(const std::vector<int64_t>& sync_vals);
// python SYNC_WORDS.get(sw) (null when unknown)
Json sync_word_network(std::optional<int> sw);
// python identify(sw, payload): plain header fields of LoRaWAN / Meshtastic, protocol first
std::optional<Json> identify(std::optional<int> sw, const std::vector<uint8_t>& payload);
// python fanet_parse(payload, strict): FANET MAC header + payload (research, PLAN G10)
std::optional<Json> fanet_parse(const std::vector<uint8_t>& payload, bool strict = false);
// python mavlink_in(payload): CRC-valid MAVLink frames in a LoRa payload -> Json object or null (PLAN G2)
Json mavlink_in(const std::vector<uint8_t>& payload);
// python ascii_view(payload)
std::string ascii_view(const std::vector<uint8_t>& payload);

}  // namespace si::rx::rc_link::lora_phy
