// RC / telemetry control links and BLE (python: receivers/rc_link/detector.py detect()).
//
//   detect(x, fs, fc[, X]) -> one Emission per burst: LoRa / ExpressLRS CSS packets
//   (blind dechirp scan, SFD verification, symbol walk, explicit-header payload decode),
//   GFSK/FSK bursts (matched-area energy, envelope / deviation / symbol rate) and BLE
//   advertising (access address, de-whitening, CRC-24); IQ-imbalance mirror images dropped.
//
// group_hops() (link grouping) is not used by the pipeline and is not ported.
#pragma once

#include <utility>
#include <vector>

#include "core/model.hpp"
#include "core/types.hpp"

namespace si::rx::rc_link {

// x: complex64 capture at fs, centre fc. X: fft(x) (scipy, complex64) if already computed.
// lora_only: (lo, hi) Hz offsets from fc -> only the LoRa scan, on sub-bands centred there (python ctx
// 'lora_only'; no FSK / BLE / ESB).
std::vector<Emission> detect(const cfvec& x, double fs, double fc, const cfvec* X = nullptr,
                             const std::pair<double, double>* lora_only = nullptr);

namespace detail {
// python suppress_images(ems, fc, min_db=3.0)
std::vector<Emission> suppress_images(const std::vector<Emission>& ems, double fc, double min_db = 3.0);
// python _is_raw_image(x, fs, fc, e, min_db=3.0)
bool is_raw_image(const cfvec& x, double fs, double fc, const Emission& e, double min_db = 3.0);
}  // namespace detail

}  // namespace si::rx::rc_link
