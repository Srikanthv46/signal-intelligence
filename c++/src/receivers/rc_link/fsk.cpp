#include "receivers/rc_link/fsk.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <mutex>

#include "core/remote_id.hpp"
#include "receivers/rc_link/bredr.hpp"
#include "receivers/rc_link/sik.hpp"
#include "core/bands_india.hpp"
#include "core/fft.hpp"
#include "core/ndimage.hpp"
#include "core/npcompat.hpp"
#include "core/parallel.hpp"
#include "core/pyfmt.hpp"
#include "core/spectral.hpp"
#include "receivers/rc_link/dsp_util.hpp"
#include "receivers/rc_link/esb.hpp"
#include "receivers/rc_link/ble_coded.hpp"
#include "receivers/rc_link/afhds2a.hpp"
#include "receivers/rc_link/cc2500.hpp"
#include "receivers/rc_link/crossfire.hpp"
#include "receivers/rc_link/adsl.hpp"
#include "receivers/rc_link/flarm.hpp"
#include "receivers/rc_link/fsk_frame.hpp"
#include "receivers/rc_link/ieee802154.hpp"
#include "receivers/rc_link/xn297.hpp"
#include "receivers/rc_link/elrs.hpp"
#include "receivers/rc_link/hitec.hpp"
#include "receivers/rc_link/hott.hpp"
#include "receivers/rc_link/hubsan.hpp"
#include "receivers/rc_link/radiolink.hpp"
#include "receivers/rc_link/toys.hpp"
#include "receivers/rc_link/gfsk_coherent.hpp"
#include "receivers/rc_link/tables.hpp"

namespace si::rx::rc_link {
namespace detail {
namespace {

constexpr uint32_t kBleAdvAA = 0x8E89BED6u;
// GSM TDMA bursts (GMSK, constant envelope) look like a sub-GHz FSK packet: one time slot (or k adjacent slots)
// long, on the 200 kHz channel raster, inside a GSM-900 / 1800 band (traffic carriers have no FCCH)
constexpr double kGsmSlotUs = 576.9;
constexpr double kGsmSlotTolUs = 60.0;      // a single burst without its guard / ramps: 546-577 us
constexpr double kGsmRasterTolHz = 15e3;    // receiver ppm (~10 ppm at 960 MHz) + carrier offset
// 2 Mb/s ESB-like packets: 2.4 GHz, ~50 us .. 1.1 ms, 99 % bandwidth 0.5-3 MHz (python ESB_BAND / ESB_US / ESB_BW)
constexpr double kEsbBandLo = 2370e6, kEsbBandHi = 2520e6;
constexpr double kEsbUsLo = 40.0, kEsbUsHi = 1100.0;
constexpr double kEsbBwLo = 0.5e6, kEsbBwHi = 3.0e6;
constexpr double kEsbMarginUs = 50.0;

// python _gsm_burst(): one to eight GSM time slots on the 200 kHz raster in a GSM band
bool gsm_burst(const FskMeas& r, double f) {
  const double dur = r.t1_us - r.t0_us;
  const int64_t k = py_round_int(dur / kGsmSlotUs);
  if (k < 1 || k > 8 || std::fabs(dur - static_cast<double>(k) * kGsmSlotUs) > kGsmSlotTolUs) return false;
  bool gsm = false;
  for (const bands::Band* b : bands::lookup(f))
    for (const auto& c : b->classes)
      if (c == "gsm_ul" || c == "gsm_dl") gsm = true;
  if (!gsm) return false;
  return std::fabs(f - 200e3 * static_cast<double>(py_round_int(f / 200e3))) <= kGsmRasterTolHz;
}

// float32 angle of a*conj(b) (np.angle of `a_arr * np.conj(b_arr)`, complex64)
inline float angle_f(const cf& a, const cf& b, bool elided) {
  const cf p = mul_conj(a, b, elided);
  return std::atan2(p.imag(), p.real());
}

// numpy np.maximum (NaN propagates)
inline float npmax(float a, float b) {
  if (std::isnan(a) || std::isnan(b)) return std::numeric_limits<float>::quiet_NaN();
  return a > b ? a : b;
}

// first index i with double(a[i]) >= v (np.searchsorted(float32 array, python float), left)
size_t searchsorted_d(const fvec& a, double v) {
  size_t lo = 0, hi = a.size();
  while (lo < hi) {
    const size_t mid = (lo + hi) / 2;
    if (static_cast<double>(a[mid]) < v) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

// ---- BLE (Core spec Vol 6 Part B): LSB-first bits, whitening x^7+x^4+1, CRC-24
// python _ble_dewhiten(bits, ch), first n bits only
std::vector<uint8_t> ble_dewhiten(const uint8_t* bits, size_t n, int ch) {
  uint8_t reg[7];
  reg[0] = 1;
  for (int i = 0; i < 6; ++i) reg[1 + i] = (ch >> (5 - i)) & 1;
  std::vector<uint8_t> out(n);
  for (size_t i = 0; i < n; ++i) {
    const uint8_t o = reg[6];
    out[i] = bits[i] ^ o;
    for (int k = 6; k > 0; --k) reg[k] = reg[k - 1];
    reg[0] = o;
    reg[4] ^= o;
  }
  return out;
}

// python _ble_crc_ok(bits, n_pdu_bits, init=0x555555)
bool ble_crc_ok(const std::vector<uint8_t>& bits, size_t n_pdu, uint32_t init = 0x555555) {
  uint8_t reg[24];
  for (int i = 0; i < 24; ++i) reg[i] = (init >> i) & 1;
  for (size_t j = 0; j < n_pdu; ++j) {
    const uint8_t fb = reg[23] ^ bits[j];
    for (int k = 23; k > 0; --k) reg[k] = reg[k - 1];
    reg[0] = fb;
    for (int p : {1, 3, 4, 6, 9, 10}) reg[p] ^= fb;
  }
  if (bits.size() < n_pdu + 24) return false;
  for (int i = 0; i < 24; ++i)
    if (bits[n_pdu + i] != reg[23 - i]) return false;
  return true;
}

const char* pdu_name(int t) {
  switch (t) {
    case 0: return "ADV_IND";
    case 1: return "ADV_DIRECT_IND";
    case 2: return "ADV_NONCONN_IND";
    case 3: return "SCAN_REQ";
    case 4: return "SCAN_RSP";
    case 5: return "CONNECT_IND";
    case 6: return "ADV_SCAN_IND";
    case 7: return "ADV_EXT_IND/AUX_*";
    default: return "?";
  }
}

// np.convolve(pw (float32 -> float64), np.ones(k) / k, "same")[lo:hi]: each output is numpy's
// _pyarray_correlate dot over the in-range part of its window, in ascending order (as
// np::convolve; numpy's BLAS ddot may differ in the last bit)
dvec conv_same_ones(const fvec& pw, int64_t k, int64_t lo, int64_t hi) {
  const int64_t n = static_cast<int64_t>(pw.size());
  const double w = 1.0 / static_cast<double>(k);
  const int64_t n_left = k / 2;
  dvec out(static_cast<size_t>(hi - lo));
  if (k > n) {  // numpy swaps the operands when v is longer: fall back to the generic routine
    const dvec full = np::convolve(dvec(pw.begin(), pw.end()), dvec(static_cast<size_t>(k), w), np::ConvMode::Same);
    for (int64_t i = lo; i < hi; ++i) out[i - lo] = full[i];
    return out;
  }
  for (int64_t i = lo; i < hi; ++i) {
    const int64_t m0 = std::max<int64_t>(0, n_left - i), m1 = std::min<int64_t>(k, n - i + n_left);
    double s = 0.0;
    for (int64_t m = m0; m < m1; ++m) s += static_cast<double>(pw[i - n_left + m]) * w;
    out[i - lo] = s;
  }
  return out;
}

// python fsk_refine's psd(v): mean |fftshift(fft(segments * hanning))|^2 (complex128)
std::optional<dvec> psd256(const cf* v, size_t len, const dvec& win) {
  const size_t nseg = 256, m = len / nseg;
  if (m < 1) return std::nullopt;
  dvec acc(nseg, 0.0);
  cdvec seg(nseg);
  for (size_t r = 0; r < m; ++r) {
    for (size_t i = 0; i < nseg; ++i) {
      const cf x = v[r * nseg + i];
      seg[i] = cd(static_cast<double>(x.real()) * win[i], static_cast<double>(x.imag()) * win[i]);
    }
    fft::c2c(seg.data(), seg.data(), nseg, true);
    for (size_t i = 0; i < nseg; ++i) acc[(i + nseg / 2) % nseg] += np::cabs2(seg[i]);
  }
  for (auto& a : acc) a /= static_cast<double>(m);
  return acc;
}

// python _fsk_measure(X, fs, fc, y, ofs, f_act, half, n0, i0, i1)
std::optional<FskMeas> fsk_measure(const cfvec& X, double fs, double fc, const cfvec& y, double ofs, double f_act,
                                   double half, double n0, int64_t i0, int64_t i1) {
  if (i1 - i0 < 16) return std::nullopt;
  const cf* seg = y.data() + i0;
  const size_t n = static_cast<size_t>(i1 - i0);
  fvec a2(n), a4(n);
  for (size_t i = 0; i < n; ++i) {
    const float a = np::cabs(seg[i]);
    a2[i] = a * a;
    a4[i] = std::pow(a, 4.0f);
  }
  const double m2 = mean_f(a2);
  const double s_pow = std::max(m2 - n0, 1e-30);
  const double snr_slice = s_pow / n0;
  const double m4 = mean_f(a4);
  const double rho = m4 / ((s_pow * s_pow + 4 * s_pow * n0) + 2 * (n0 * n0));
  const size_t nf = size_t(1) << static_cast<int>(std::ceil(std::log2(static_cast<double>(std::max<size_t>(n, 256)))));
  const cfvec F = fft::np_fft(seg, n, nf);
  fvec Sx(nf);
  const float sub = static_cast<float>(n0 * static_cast<double>(n));
  for (size_t i = 0; i < nf; ++i) {
    float v = np::cabs(F[i]);
    v = v * v - sub;
    Sx[(i + nf / 2) % nf] = v < 0 ? 0.0f : v;
  }
  const dvec fr = fft::fftshift(fft::fftfreq(nf, 1.0 / ofs));
  const float ssum = np::pairwise_sum(Sx.data(), nf);
  const bool floor = 1e-30 > static_cast<double>(ssum);
  const float den_f = floor ? static_cast<float>(1e-30) : ssum;
  fvec cs(nf);
  float acc = 0;
  for (size_t i = 0; i < nf; ++i) {
    acc = i == 0 ? Sx[0] : acc + Sx[i];
    cs[i] = acc / den_f;
  }
  const size_t il = std::min(searchsorted_d(cs, 0.005), nf - 1);  // python: IndexError past the end
  const size_t ih = std::min(searchsorted_d(cs, 0.995), nf - 1);
  const double f_l = fr[il], f_h = fr[ih];
  dvec fw(nf);
  for (size_t i = 0; i < nf; ++i) fw[i] = fr[i] * static_cast<double>(Sx[i]);
  const double centroid = np::pairwise_sum(fw.data(), nf) / (floor ? 1e-30 : static_cast<double>(ssum));
  FskMeas r;
  r.t0_us = static_cast<double>(i0) / ofs * 1e6;
  r.t1_us = static_cast<double>(i1) / ofs * 1e6;
  r.f_center = fc + f_act + centroid;
  r.f_lo = fc + f_act + f_l;
  r.f_hi = fc + f_act + f_h;
  r.rho = rho;
  r.ofs = ofs;
  r.n0 = n0;
  r.bw99 = f_h - f_l;
  r.snr_db = 10 * std::log10(std::max(snr_slice * (2 * half) / std::max(f_h - f_l, 50e3), 1e-6));
  if (snr_slice > 2) r.const_env = rho < 1.25;
  if (snr_slice > 4 && n > 64) {
    fvec ifr(n - 1);
    const float of = static_cast<float>(ofs), tp = static_cast<float>(2 * kPi);
    const bool el = elides(n - 1, sizeof(cf));
    for (size_t i = 0; i + 1 < n; ++i) ifr[i] = (angle_f(seg[i + 1], seg[i], el) * of) / tp;
    const SymRate sr = symbol_rate(ifr, ofs);
    r.has_rate = true;
    r.sym_rate = sr.rate;
    r.sym_conc = sr.conc;
    if (sr.rate) {
      const double R = *sr.rate, sps = ofs / R;
      const fvec v = ndi::uniform_filter1d(ifr, static_cast<size_t>(std::max<int64_t>(1, static_cast<int64_t>(sps / 2))));
      const double t0 = (sr.ph / (2 * kPi)) / R;
      const int64_t ni = static_cast<int64_t>(static_cast<double>(v.size()) / sps);
      if (ni < 8) return r;
      fvec vi(static_cast<size_t>(ni));
      const double off = (t0 + 0.5 / R) * ofs;
      for (int64_t i = 0; i < ni; ++i)
        vi[i] = v[static_cast<size_t>(pymod(static_cast<int64_t>(off + static_cast<double>(i) * sps),
                                            static_cast<int64_t>(v.size())))];
      const auto [lo, hi] = kmeans2(vi);
      r.has_dev = true;
      r.dev_hz = static_cast<double>((hi - lo) / 2.0f);
      r.fsk_center = static_cast<double>(static_cast<float>(fc + f_act) + (hi + lo) / 2.0f);
      r.mod_index = static_cast<double>((hi - lo) / static_cast<float>(R));
      if (0.85e6 < R && R < 1.15e6 && 2370e6 <= r.fsk_center && r.fsk_center <= 2520e6 &&
          static_cast<double>(i1 - i0) / ofs >= 120e-6) {
        // BLE: re-slice +-0.75 MHz at exactly 8 samples/bit around the FSK centre
        const SliceOut sb = slice(X, fs, FOff{r.fsk_center - fc, false}, 0.75e6, 8e6);
        const double ob = sb.ofs;
        const int64_t ny = static_cast<int64_t>(sb.y.size());
        const int64_t a = std::min(ny, std::max<int64_t>(0, static_cast<int64_t>(static_cast<double>(i0) / ofs * ob) - 80));
        int64_t b = static_cast<int64_t>(static_cast<double>(i1) / ofs * ob) + 120;
        b = std::min(ny, b);
        const cfvec yb(sb.y.begin() + a, sb.y.begin() + std::max(a, b));
        Json d = ble_try_decode(yb, ob, r.fsk_center);
        if (d.is_object() && !d.empty()) {
          r.ble = std::move(d);
        } else {
          // XN297 toy drones / Bayang (research, PLAN R8): prefix 0x710F55 + Bayang checksum
          Json x = xn297::try_decode(yb, ob);
          if (x.is_object()) {
            r.xn297 = std::move(x);
          } else {
            // plain nRF24 toy drones (V2x2, SymaX; research, PLAN G6): fixed address + nRF CRC-16
            Json t = toys::nrf24_try_decode(yb, ob);
            if (t.is_object()) {
              r.xn297 = std::move(t);
            } else {
              // Bluetooth BR/EDR (research, signal gaps S7): (64,30) BCH access code -> LAP
              Json br = bredr::try_decode(yb, ob);
              if (br.is_object()) r.bredr = std::move(br);
            }
          }
        }
      } else if (0.45e6 < R && R < 0.55e6 && 2370e6 <= r.fsk_center && r.fsk_center <= 2520e6 && 0.6e6 <= r.bw99 &&
                 r.bw99 <= 2.0e6 && static_cast<double>(i1 - i0) / ofs >= 400e-6) {
        // BLE LE Coded S=8 (BT5 Long Range): the 0011 / 1100 patterns look like 500 kBd (research, PLAN T2b)
        const SliceOut sb = slice(X, fs, FOff{r.fsk_center - fc, false}, 0.75e6, 8e6);
        const double ob = sb.ofs;
        const int64_t ny = static_cast<int64_t>(sb.y.size());
        const int64_t a = std::min(ny, std::max<int64_t>(0, static_cast<int64_t>(static_cast<double>(i0) / ofs * ob) - 80));
        const int64_t b = std::min(ny, static_cast<int64_t>(static_cast<double>(i1) / ofs * ob) + 120);
        const cfvec yb(sb.y.begin() + a, sb.y.begin() + std::max(a, b));
        Json d = ble_coded_decode(yb, ob, ble_channel_index(r.fsk_center));
        if (d.is_object() && !d.empty()) r.ble = std::move(d);
      } else if (2370e6 <= r.fsk_center && r.fsk_center <= 2520e6 && 27e3 <= R && R <= 110e3 && 15e3 <= r.dev_hz &&
                 r.dev_hz <= 95e3) {
        // FrSky D8 / D16 (CC2500; research, PLAN R6): a nominal rate within 8 % -> re-slice at 8 samples/bit, sync +
        // length + CRC -> variant and TX ID; read from this piece's start over the longest packet
        if (std::fabs(R - radiolink::kRate) <= 0.08 * radiolink::kRate &&
            static_cast<double>(i1 - i0) / ofs * radiolink::kRate >= 32) {
          // RadioLink (CC2500; research, PLAN G4): doubled sync, PN9 dewhitening, CRC + checksum, hop-table check
          const SliceOut sr = slice(X, fs, FOff{r.fsk_center - fc, false},
                                    std::min(r.dev_hz + radiolink::kRate, 0.45 * 8 * radiolink::kRate), 8 * radiolink::kRate);
          const double orr = sr.ofs;
          const int64_t ny = static_cast<int64_t>(sr.y.size());
          const int64_t j0 = std::max<int64_t>(0, static_cast<int64_t>(static_cast<double>(i0) / ofs * orr) - 40);
          const int64_t j1 = std::max<int64_t>(static_cast<int64_t>(static_cast<double>(i1) / ofs * orr) + 40,
                                               j0 + (8 + radiolink::kMaxLen + 3) * 8 * 8 + 80);
          const int64_t a = std::min(ny, j0), b = std::min(ny, j1);
          Json d = radiolink::try_decode(cfvec(sr.y.begin() + a, sr.y.begin() + std::max(a, b)), orr, r.fsk_center);
          if (d.is_object()) r.radiolink = std::move(d);
        }
        double kb = cc2500::kRates[0].bps;
        for (const auto& v : cc2500::kRates)
          if (std::fabs(v.bps - R) < std::fabs(kb - R)) kb = v.bps;
        if (std::fabs(R - kb) <= 0.08 * kb && static_cast<double>(i1 - i0) / ofs * kb >= 64) {
          const SliceOut sr = slice(X, fs, FOff{r.fsk_center - fc, false}, std::min(r.dev_hz + kb, 0.45 * 8 * kb), 8 * kb);
          const double orr = sr.ofs;
          const int64_t ny = static_cast<int64_t>(sr.y.size());
          const int64_t j0 = std::max<int64_t>(0, static_cast<int64_t>(static_cast<double>(i0) / ofs * orr) - 40);
          const int64_t j1 = std::max<int64_t>(static_cast<int64_t>(static_cast<double>(i1) / ofs * orr) + 40,
                                               j0 + (8 + 36) * 8 * 8 + 80);
          const int64_t a = std::min(ny, j0), b = std::min(ny, j1);
          const cfvec yr(sr.y.begin() + a, sr.y.begin() + std::max(a, b));
          Json d = cc2500::try_identify(yr, orr, kb);
          if (d.is_object()) r.frsky = std::move(d);
        }
      } else if (r.fsk_center < 1.5e9 && r.dev_hz > 0) {
        // SiK telemetry (research): air rate within 5 %, deviation ~ rate (h ~ 2); decode from this piece's start over
        // the longest SiK packet (python _fsk_measure)
        double kb = sik::kRates[0];
        for (double v : sik::kRates)
          if (std::fabs(v - R) < std::fabs(kb - R)) kb = v;
        if (std::fabs(R - kb) <= 0.05 * kb && r.dev_hz / kb >= 0.6 && r.dev_hz / kb <= 1.4 &&
            static_cast<double>(i1 - i0) / ofs * kb >= 32) {
          const double half_s = std::min(r.dev_hz + kb, 0.45 * 8 * kb);
          const SliceOut ss = slice(X, fs, FOff{r.fsk_center - fc, false}, half_s, 8 * kb);
          const double os = ss.ofs;
          const int64_t ny = static_cast<int64_t>(ss.y.size());
          const int64_t j0 = std::max<int64_t>(0, static_cast<int64_t>(static_cast<double>(i0) / ofs * os) - 40);
          const int64_t j1 = std::max<int64_t>(static_cast<int64_t>(static_cast<double>(i1) / ofs * os) + 40,
                                               j0 + static_cast<int64_t>((64 + 16 + 8 * 260) * 8 + 80));
          const int64_t a = std::min(ny, j0), b = std::min(ny, j1);
          const cfvec ys(ss.y.begin() + a, ss.y.begin() + std::max(a, b));
          Json d = sik::sik_try_decode(ys, os, kb);
          if (d.is_object() && !d.empty()) {
            d["crc_ok"] = (d.contains("mavlink_crc_ok") && d["mavlink_crc_ok"].get<int>() > 0) ||
                          (d.contains("hw_crc_ok") && d["hw_crc_ok"].get<bool>());
            r.sik = std::move(d);
          }
        }
        if (std::fabs(R - flarm::kChipRate) <= 0.06 * flarm::kChipRate && 30e3 <= r.dev_hz && r.dev_hz <= 70e3 &&
            flarm::in_band(r.fsk_center) && static_cast<double>(i1 - i0) / ofs * flarm::kChipRate >= 40) {
          // FLARM beacon (research, PLAN T5): Manchester 100 kchip/s, sync F5 31 FA B6 -> ID, address type;
          // >= 40 chips: a packet split into pieces is read from the piece holding preamble + sync
          const SliceOut sl = slice(X, fs, FOff{r.fsk_center - fc, false}, flarm::kDev + flarm::kChipRate, 8 * flarm::kChipRate);
          const double ol = sl.ofs;
          const int64_t ny = static_cast<int64_t>(sl.y.size());
          const int64_t j0 = std::max<int64_t>(0, static_cast<int64_t>(static_cast<double>(i0) / ofs * ol) - 40);
          const int64_t j1 = std::max<int64_t>(static_cast<int64_t>(static_cast<double>(i1) / ofs * ol) + 40,
                                               j0 + (16 + 2 * 8 * 30) * 8 + 80);
          auto seg = [&](int64_t e) {
            const int64_t a = std::min(ny, j0), b = std::min(ny, e);
            return cfvec(sl.y.begin() + a, sl.y.begin() + std::max(a, b));
          };
          // ADS-L SRD-860 M-band, then the OGN tracker (research, PLAN G10): same PHY, CRC-24 / Manchester-clean
          Json da = adsl::try_decode(seg(std::max(j1, j0 + (16 + 2 * 8 * 32) * 8 + 80)), ol);
          if (!da.is_object()) da = ogn::try_decode(seg(std::max(j1, j0 + (16 + 2 * 8 * 30) * 8 + 80)), ol);
          if (da.is_object()) {
            r.adsl = std::move(da);
          } else {
            Json d = flarm::try_identify(seg(j1), ol);
            if (d.is_object()) r.flarm = std::move(d);
          }
        }
        if (std::fabs(R - crossfire::kRate) <= 0.05 * crossfire::kRate && 30e3 <= r.dev_hz && r.dev_hz <= 55e3 &&
            850e6 <= r.fsk_center && r.fsk_center <= 935e6 && static_cast<double>(i1 - i0) / ofs * crossfire::kRate >= 48) {
          // TBS Crossfire 150 Hz (research, PLAN R9): read from this piece's start over the longest packet
          const SliceOut sc = slice(X, fs, FOff{r.fsk_center - fc, false}, crossfire::kDev + crossfire::kRate,
                                    8 * crossfire::kRate);
          const double oc = sc.ofs;
          const int64_t ny = static_cast<int64_t>(sc.y.size());
          const int64_t j0 = std::max<int64_t>(0, static_cast<int64_t>(static_cast<double>(i0) / ofs * oc) - 40);
          const int64_t j1 = std::max<int64_t>(static_cast<int64_t>(static_cast<double>(i1) / ofs * oc) + 40,
                                               j0 + (8 * 8 + 24 + 8 * 23) * 8 + 80);
          const int64_t a = std::min(ny, j0), b = std::min(ny, j1);
          const cfvec yc(sc.y.begin() + a, sc.y.begin() + std::max(a, b));
          Json d = crossfire::try_decode(yc, oc, r.fsk_center);
          if (d.is_object()) r.crossfire = std::move(d);
        }
      }
    }
  }
  {
    // XN297 toy drones at h ~0.3: the symbol-rate estimate is unreliable, so a packet-sized 2.4 GHz burst gets one
    // XN297 attempt at 1 Mb/s (research, PLAN R8)
    const double R0 = r.sym_rate ? *r.sym_rate : 0.0;
    const double fcen = r.has_dev ? r.fsk_center : 0.0;
    const double dur = static_cast<double>(i1 - i0) / ofs;
    const bool ble_branch = 0.85e6 < R0 && R0 < 1.15e6 && 2370e6 <= fcen && fcen <= 2520e6 && dur >= 120e-6;
    if (!r.ble.is_object() && !r.xn297.is_object() && !ble_branch && !(r.const_env && !*r.const_env) &&
        2370e6 <= r.f_center && r.f_center <= 2520e6 && 100e-6 <= dur && dur <= 600e-6 && 0.3e6 <= r.bw99 &&
        r.bw99 <= 1.6e6 && !(0.45e6 < R0 && R0 < 0.55e6)) {
      const double fx = (r.has_dev && r.fsk_center != 0.0) ? r.fsk_center : r.f_center;
      const SliceOut sx = slice(X, fs, FOff{fx - fc, false}, 0.75e6, 8e6);
      const double ox = sx.ofs;
      const int64_t ny = static_cast<int64_t>(sx.y.size());
      const int64_t a = std::min(ny, std::max<int64_t>(0, static_cast<int64_t>(static_cast<double>(i0) / ofs * ox) - 80));
      const int64_t b = std::min(ny, static_cast<int64_t>(static_cast<double>(i1) / ofs * ox) + 120);
      const cfvec yxs(sx.y.begin() + a, sx.y.begin() + std::max(a, b));
      Json d = xn297::try_decode(yxs, ox);
      if (!d.is_object()) d = toys::nrf24_try_decode(yxs, ox);   // + plain nRF24 toys (research, PLAN G6)
      if (d.is_object()) r.xn297 = std::move(d);
    }
  }
  if (r.sym_rate && *r.sym_rate != 0.0 && r.has_dev && r.dev_hz != 0.0 && !r.ble.is_object() && 0.45e6 < *r.sym_rate &&
      *r.sym_rate < 0.55e6 && 2370e6 <= r.fsk_center && r.fsk_center <= 2520e6 && 120e3 <= r.dev_hz && r.dev_hz <= 260e3 &&
      static_cast<double>(i1 - i0) / ofs * *r.sym_rate >= 64) {
    // FlySky AFHDS 2A (A7105; research, PLAN R7): ID code -> identification; FEC + CRC under candidate conventions
    const SliceOut sa = slice(X, fs, FOff{r.fsk_center - fc, false}, afhds2a::kDev + afhds2a::kRate, 8 * afhds2a::kRate);
    const double oa = sa.ofs;
    const int64_t ny = static_cast<int64_t>(sa.y.size());
    const int64_t j0 = std::max<int64_t>(0, static_cast<int64_t>(static_cast<double>(i0) / ofs * oa) - 40);
    const int64_t j1 = std::max<int64_t>(static_cast<int64_t>(static_cast<double>(i1) / ofs * oa) + 40,
                                         j0 + (64 + (afhds2a::kPayload + 2) * 14) * 8 + 80);
    const int64_t a = std::min(ny, j0), b = std::min(ny, j1);
    Json d = afhds2a::try_identify(cfvec(sa.y.begin() + a, sa.y.begin() + std::max(a, b)), oa, afhds2a::kRate,
                                   r.fsk_center);
    if (d.is_object()) r.afhds2a = std::move(d);
  }
  if (r.sym_rate && *r.sym_rate != 0.0 && r.has_dev && r.dev_hz != 0.0 && !r.ble.is_object() && 250e3 <= *r.sym_rate &&
      *r.sym_rate <= 350e3 && 60e3 <= r.dev_hz && r.dev_hz <= 140e3 &&
      ((2400e6 <= r.fsk_center && r.fsk_center <= 2480e6) || (860e6 <= r.fsk_center && r.fsk_center <= 930e6)) &&
      static_cast<double>(i1 - i0) / ofs <= 700e-6) {
    // ExpressLRS LR1121 GFSK modes (research, PLAN G3): a self-verifying SYNC packet whose UID bytes = the sync word
    const SliceOut se = slice(X, fs, FOff{r.fsk_center - fc, false}, elrs::kGfskDev + elrs::kGfskRate, 8 * elrs::kGfskRate);
    const double oe = se.ofs;
    const int64_t ny = static_cast<int64_t>(se.y.size());
    const int64_t a = std::min(ny, std::max<int64_t>(0, static_cast<int64_t>(static_cast<double>(i0) / ofs * oe) - 80));
    const int64_t b = std::min(ny, static_cast<int64_t>(static_cast<double>(i1) / ofs * oe) + 120);
    Json d = elrs::gfsk_try_decode(cfvec(se.y.begin() + a, se.y.begin() + std::max(a, b)), oe);
    if (d.is_object()) r.elrs_gfsk = std::move(d);
  }
  if (r.sym_rate && *r.sym_rate != 0.0 && r.has_dev && r.dev_hz != 0.0 && !r.ble.is_object() && 2400e6 <= r.fsk_center &&
      r.fsk_center <= 2490e6 && 120e3 <= r.dev_hz && r.dev_hz <= 260e3 && 1.5e-3 <= static_cast<double>(i1 - i0) / ofs &&
      static_cast<double>(i1 - i0) / ofs <= 2.6e-3) {
    // Hubsan (A7105 100 kb/s; research, PLAN G5): gated on deviation and packet length; the decode verifies itself
    const SliceOut sh = slice(X, fs, FOff{r.fsk_center - fc, false}, hubsan::kDev + hubsan::kRate, 8 * hubsan::kRate);
    const double oh = sh.ofs;
    const int64_t ny = static_cast<int64_t>(sh.y.size());
    const int64_t j0 = std::max<int64_t>(0, static_cast<int64_t>(static_cast<double>(i0) / ofs * oh) - 40);
    const int64_t j1 = std::max<int64_t>(static_cast<int64_t>(static_cast<double>(i1) / ofs * oh) + 40,
                                         j0 + (64 + 8 * (hubsan::kLen + 2)) * 8 + 80);
    const int64_t a = std::min(ny, j0), b = std::min(ny, j1);
    Json d = hubsan::try_decode(cfvec(sh.y.begin() + a, sh.y.begin() + std::max(a, b)), oh, r.fsk_center);
    if (d.is_object()) r.hubsan = std::move(d);
  }
  const double dur = static_cast<double>(i1 - i0) / ofs;
  if (r.sym_rate && *r.sym_rate != 0.0 && r.has_dev && r.dev_hz != 0.0 && !r.ble.is_object() && 2400e6 <= r.fsk_center &&
      r.fsk_center <= 2490e6 && 200e3 <= *r.sym_rate && *r.sym_rate <= 310e3 && 30e3 <= r.dev_hz && r.dev_hz <= 110e3 &&
      0.6e-3 <= dur && dur <= 4.4e-3) {
    // Graupner HoTT (MSK + CC2500 FEC) and Hitec Optima / Minima (MSK, sync = TX ID; research, PLAN G12)
    const SliceOut sk = slice(X, fs, FOff{r.fsk_center - fc, false}, hott::kDev + hott::kRate, 8 * hott::kRate);
    const double ok_ = sk.ofs;
    const int64_t ny = static_cast<int64_t>(sk.y.size());
    const int64_t j0 = std::max<int64_t>(0, static_cast<int64_t>(static_cast<double>(i0) / ofs * ok_) - 80);
    const int64_t e1 = static_cast<int64_t>(static_cast<double>(i1) / ofs * ok_) + 80;
    auto seg = [&](int64_t e) {
      const int64_t a = std::min(ny, j0), b = std::min(ny, e);
      return cfvec(sk.y.begin() + a, sk.y.begin() + std::max(a, b));
    };
    Json d = hott::try_decode(seg(std::max<int64_t>(e1, j0 + (64 + hott::coded_len(hott::kLenTx)) * 8 + 160)), ok_, r.fsk_center);
    if (!d.is_object()) d = hitec::hitec_try_decode(seg(std::max<int64_t>(e1, j0 + (8 * 60) * 8 + 160)), ok_, r.fsk_center);
    if (d.is_object()) r.cc2500_rc = std::move(d);
  }
  if (r.sym_rate && *r.sym_rate != 0.0 && r.has_dev && r.dev_hz != 0.0 && !r.ble.is_object() && 2400e6 <= r.fsk_center &&
      r.fsk_center <= 2450e6 && 100e3 <= *r.sym_rate && *r.sym_rate <= 160e3 && 20e3 <= r.dev_hz && r.dev_hz <= 60e3 &&
      0.4e-3 <= dur && dur <= 1.9e-3) {
    // Futaba S-FHSS (CC2500 2-FSK 128 kBd, sync D391, CRC in CC2400 mode; research, PLAN G12)
    const SliceOut ss = slice(X, fs, FOff{r.fsk_center - fc, false}, hitec::kSfDev + hitec::kSfRate, 8 * hitec::kSfRate);
    const double os_ = ss.ofs;
    const int64_t ny = static_cast<int64_t>(ss.y.size());
    const int64_t j0 = std::max<int64_t>(0, static_cast<int64_t>(static_cast<double>(i0) / ofs * os_) - 80);
    const int64_t e1 = std::max<int64_t>(static_cast<int64_t>(static_cast<double>(i1) / ofs * os_) + 80, j0 + 8 * 8 * 24 + 160);
    const int64_t a = std::min(ny, j0), b = std::min(ny, e1);
    Json d = hitec::sfhss_try_decode(cfvec(ss.y.begin() + a, ss.y.begin() + std::max(a, b)), os_, r.fsk_center);
    if (d.is_object()) r.cc2500_rc = std::move(d);
  }
  if (r.sym_rate && *r.sym_rate != 0.0 && r.has_dev && r.dev_hz != 0.0 && !r.ble.is_object() && !r.sik.is_object() &&
      !r.crossfire.is_object() && !r.frsky.is_object() && !r.afhds2a.is_object() && !r.xn297.is_object() && !r.flarm.is_object() &&
      !r.radiolink.is_object() && !r.hubsan.is_object() && !r.elrs_gfsk.is_object() && !r.cc2500_rc.is_object() &&
      !r.adsl.is_object() && 20e3 <= *r.sym_rate && *r.sym_rate <= 1.2e6 && r.snr_db >= 10 &&
      static_cast<double>(i1 - i0) / ofs * *r.sym_rate >= fsk_frame::kMinPre + fsk_frame::kSyncBits + 16) {
    // any other FSK packet: blind frame structure (preamble, 32-bit sync-word candidate, length; research)
    const double R = *r.sym_rate;
    const SliceOut sf = slice(X, fs, FOff{r.fsk_center - fc, false}, std::min(r.dev_hz + R, 0.45 * 8 * R), 8 * R);
    const double of = sf.ofs;
    const int64_t ny = static_cast<int64_t>(sf.y.size());
    const int64_t a = std::min(ny, std::max<int64_t>(0, static_cast<int64_t>(static_cast<double>(i0) / ofs * of) - 16));
    const int64_t b = std::min(ny, static_cast<int64_t>(static_cast<double>(i1) / ofs * of) + 16);
    Json fr = fsk_frame::frame(cfvec(sf.y.begin() + a, sf.y.begin() + std::max(a, b)), of, R);
    if (fr.is_object()) r.frame = std::move(fr);
  }
  if (!(r.sym_rate && *r.sym_rate != 0.0) && 2370e6 <= r.f_center && r.f_center <= 2520e6 &&
      static_cast<double>(i1 - i0) / ofs >= 120e-6 && 0.5e6 <= r.bw99 && r.bw99 <= 2.5e6) {
    // low SNR (PLAN L1): no symbol rate below ~6 dB slice SNR, but a BLE-wide, advert-long 2.4 GHz burst is worth
    // one coherent access-address search + decode (nothing is reported without the CRC)
    const SliceOut sb = slice(X, fs, FOff{r.f_center - fc, false}, 0.75e6, 8e6);
    const double ob = sb.ofs;
    const int64_t ny = static_cast<int64_t>(sb.y.size());
    const int64_t a = std::min(ny, std::max<int64_t>(0, static_cast<int64_t>(static_cast<double>(i0) / ofs * ob) - 80));
    const int64_t b = std::min(ny, static_cast<int64_t>(static_cast<double>(i1) / ofs * ob) + 120);
    const cfvec yb(sb.y.begin() + a, sb.y.begin() + std::max(a, b));
    const int chx = ble_channel_index(r.f_center);
    Json d = ble_coherent_decode(yb, ob, chx);
    if (!d.is_object()) d = ble_coded_decode(yb, ob, chx);
    if (d.is_object() && !d.empty()) r.ble = std::move(d);
  }
  return r;
}

// python fsk_refine(X, fs, fc, cand, noise_win_us=2000.0) on the local FFT X
std::vector<FskMeas> fsk_refine(const cfvec& X, double fs, double fc, const FskCand& cand, double noise_win_us = 2000.0) {
  const double bw0 = cand.f_hi - cand.f_lo;
  double half = std::max(0.75 * bw0 + 60e3, 150e3);
  double ofs = std::min({std::max(8 * half, 2e6), 16e6, fs});
  SliceOut sl = slice(X, fs, FOff{cand.fc_b - fc, false}, half, ofs);
  cfvec y = std::move(sl.y);
  ofs = sl.ofs;
  double f_act = sl.f_act;
  int64_t ny = static_cast<int64_t>(y.size());
  int64_t a = std::max<int64_t>(static_cast<int64_t>(cand.t0_us * 1e-6 * ofs), 0);
  int64_t b = std::min<int64_t>(static_cast<int64_t>(cand.t1_us * 1e-6 * ofs), ny);
  if (b - a < 16) return {};
  // pass 1: occupied bandwidth (PSD minus the local noise PSD), then re-slice tightly
  const size_t nseg = 256;
  const dvec win = np::hanning(nseg);
  const int64_t w0 = static_cast<int64_t>(noise_win_us * 1e-6 * ofs);
  auto ps = psd256(y.data() + a, static_cast<size_t>(b - a), win);
  cfvec noise;
  {
    const int64_t s1 = std::max<int64_t>(0, a - w0), e1 = std::max<int64_t>(0, a - 32);
    const int64_t s2 = std::min(ny, b + 32), e2 = std::min(ny, b + w0);
    if (e1 > s1) noise.insert(noise.end(), y.begin() + s1, y.begin() + e1);
    if (e2 > s2) noise.insert(noise.end(), y.begin() + s2, y.begin() + e2);
  }
  auto pn = psd256(noise.data(), noise.size(), win);
  if (ps && pn) {
    dvec ex(nseg);
    for (size_t i = 0; i < nseg; ++i) {
      const double d = (*ps)[i] - (*pn)[i];
      ex[i] = d < 0 ? 0.0 : d;
    }
    ex = ndi::uniform_filter1d(ex, 5);
    const dvec frq = fft::fftshift(fft::fftfreq(nseg, 1.0 / ofs));
    const double mx = *std::max_element(ex.begin(), ex.end());
    int64_t first = -1, last = -1;
    for (size_t i = 0; i < nseg; ++i)
      if (ex[i] > 0.1 * mx) {
        if (first < 0) first = static_cast<int64_t>(i);
        last = static_cast<int64_t>(i);
      }
    if (first >= 0) {
      const double f_lo_b = frq[first], f_hi_b = frq[last];
      half = std::max(0.6 * (f_hi_b - f_lo_b) + 2 * ofs / nseg, 60e3);
      const double ofs2_req = std::min({std::max(8 * half, 1e6), 16e6, fs});
      SliceOut s2 = slice(X, fs, FOff{f_act + (f_lo_b + f_hi_b) / 2, false}, half, ofs2_req);
      y = std::move(s2.y);
      ny = static_cast<int64_t>(y.size());
      f_act = s2.f_act;
      const double ofs2 = s2.ofs;
      a = static_cast<int64_t>(static_cast<double>(a) / ofs * ofs2);
      b = static_cast<int64_t>(static_cast<double>(b) / ofs * ofs2);
      ofs = ofs2;
    }
  }
  fvec pw(y.size());
  for (size_t i = 0; i < y.size(); ++i) pw[i] = np::cabs2(y[i]);
  const int64_t w = static_cast<int64_t>(noise_win_us * 1e-6 * ofs);
  const int64_t i10 = static_cast<int64_t>(10e-6 * ofs);
  fvec outv;
  {
    const int64_t s1 = std::max<int64_t>(0, a - w), e1 = std::max<int64_t>(0, a - i10);
    const int64_t s2 = std::min(ny, b + i10), e2 = std::min(ny, b + w);
    for (int64_t i = s1; i < e1; ++i) outv.push_back(pw[i]);
    for (int64_t i = s2; i < e2; ++i) outv.push_back(pw[i]);
  }
  const double n0 = outv.size() > 32 ? static_cast<double>(np::median(outv)) / std::log(2.0)
                                     : static_cast<double>(np::percentile(pw, 20) / 0.223f);
  // precise edges: smoothed power crossing halfway between noise and the burst level
  const int64_t k = std::max<int64_t>(1, i10);
  // env = np.convolve(pw, ones(k)/k, "same"), evaluated only on [base, end) (the samples used)
  const int64_t base = std::max<int64_t>(0, a - 5 * k), end = std::min(ny, b + 5 * k);
  const dvec env = conv_same_ones(pw, k, base, std::max(base, end));
  const int64_t ea = std::min(a, ny), eb = std::min(b, ny);
  if (eb <= ea) return {};  // python: percentile of an empty slice (NaN); never seen
  const dvec env_ab(env.begin() + (ea - base), env.begin() + (eb - base));
  double lvl = np::percentile(env_ab, 75);
  if (lvl < 1.15 * n0) {
    // a short strong packet in a long box (fsk_scan's 800 us integration smears a 1 kHz train of ~70-100 us
    // packets into ~870 us boxes): its level is the top of the box; noise-only stays below 1.3 n0 there
    lvl = np::percentile(env_ab, 98);
    if (lvl < 4 * n0) return {};
  }
  const double thr = n0 + 0.5 * (lvl - n0);
  std::vector<uint8_t> on(static_cast<size_t>(std::max<int64_t>(0, end - base)));
  for (int64_t i = base; i < end; ++i) on[i - base] = env[i - base] > thr;
  on = ndi::binary_closing(on, ndi::ones1(static_cast<size_t>(std::max<int64_t>(1, static_cast<int64_t>(8e-6 * ofs)))));
  const auto lab = ndi::label(on);
  std::vector<FskMeas> out;
  int nruns = 0;
  for (const auto& [r0, r1] : ndi::find_objects(lab.lab, lab.n)) {
    if (static_cast<double>(r1 - r0) / ofs * 1e6 < 20) continue;
    ++nruns;
    if (auto m = fsk_measure(X, fs, fc, y, ofs, f_act, half, n0, base + static_cast<int64_t>(r0),
                             base + static_cast<int64_t>(r1)))
      out.push_back(std::move(*m));
  }
  auto crc = [](const FskMeas& m) {
    return m.ble.is_object() && m.ble.contains("crc_ok") && m.ble["crc_ok"].is_boolean() && m.ble["crc_ok"].get<bool>();
  };
  const double f_box = fc + f_act;
  const bool ble_box = 2370e6 <= f_box && f_box <= 2520e6 && std::none_of(out.begin(), out.end(), crc);
  const bool cf_box = 850e6 <= f_box && f_box <= 935e6 &&
                      std::none_of(out.begin(), out.end(), [](const FskMeas& m) { return m.crossfire.is_object(); });
  if (nruns >= 2 && static_cast<double>(b - a) / ofs >= 400e-6 && (ble_box || cf_box)) {
    // a long packet splits into pieces (BT5 Long Range near the noise, or an S=2 block after the S=8 one): one
    // decode attempt over the whole box (python fsk_refine; research, PLAN T2b)
    auto m = fsk_measure(X, fs, fc, y, ofs, f_act, half, n0, std::max<int64_t>(a, 0),
                         std::min<int64_t>(b, static_cast<int64_t>(y.size())));
    if (m && (crc(*m) || m->crossfire.is_object())) return {std::move(*m)};
  }
  return out;
}

}  // namespace

std::vector<FskCand> fsk_scan(const cfvec& x, double fs, double fc) {
  const size_t nfft = 512, hop = nfft / 2;
  const double k_sigma = 8.0, max_bw = 3.0e6, min_us = 30.0, edge = 0.45;
  const int t_frames = 12;
  const SpectralSpectrogram sg = spectrogram(x, fs, fc, nfft, hop);
  const dmat& P = sg.P;
  const size_t R = P.rows, Cn = P.cols;
  if (R == 0) return {};
  dvec fl = np::percentile_axis0(P, 10);
  for (auto& v : fl) v /= 0.10536;
  fl = ndi::median_filter(fl, 15);
  dmat Z(R, Cn);
  for (size_t r = 0; r < R; ++r)
    for (size_t c = 0; c < Cn; ++c) Z(r, c) = P(r, c) / fl[c];
  std::vector<uint8_t> keep(Cn);
  for (size_t c = 0; c < Cn; ++c) keep[c] = std::fabs(sg.f_mhz[c] - fc / 1e6) <= edge * fs / 1e6;
  // power over the local (+-3 MHz, per frame) median: rows are independent
  dmat loc(R, Cn);
  {
    const size_t chunk = 256, nch = (R + chunk - 1) / chunk;
    parallel_for(nch, [&](size_t ci) {
      const size_t r0 = ci * chunk, r1 = std::min(R, r0 + chunk);
      dmat part(r1 - r0, Cn);
      std::copy(P.row(r0), P.row(r0) + (r1 - r0) * Cn, part.data.begin());
      const dmat m = ndi::median_filter(part, 1, 49, ndi::Mode::Nearest);
      std::copy(m.data.begin(), m.data.end(), loc.row(r0));
    });
  }
  dmat Zl(R, Cn);
  for (size_t r = 0; r < R; ++r)
    for (size_t c = 0; c < Cn; ++c) Zl(r, c) = P(r, c) / std::max(loc(r, c), fl[c]);
  loc = dmat();
  struct Job {
    const dmat* Z;
    size_t wf, tf;
  };
  std::vector<Job> jobs;
  for (size_t wf : {1, 3, 7})
    for (size_t tf : {t_frames, 4 * t_frames, 16 * t_frames}) jobs.push_back({&Z, wf, tf});
  for (size_t wf : {1, 3, 7})
    for (size_t tf : {t_frames, 4 * t_frames}) jobs.push_back({&Zl, wf, tf});
  fmat E(R, Cn, -std::numeric_limits<float>::infinity());
  std::mutex mu;
  parallel_for(jobs.size(), [&](size_t j) {
    const Job& jb = jobs[j];
    const dmat Zf = ndi::uniform_filter(*jb.Z, jb.tf, jb.wf, ndi::Mode::Nearest);
    const double sq = std::sqrt(static_cast<double>(jb.wf * jb.tf) / 1.5);
    fvec v(Zf.size());
    for (size_t i = 0; i < v.size(); ++i) v[i] = static_cast<float>((Zf.data[i] - 1) * sq);
    std::lock_guard<std::mutex> lk(mu);
    for (size_t i = 0; i < v.size(); ++i) E.data[i] = npmax(E.data[i], v[i]);
  });
  for (size_t r = 0; r < R; ++r)
    for (size_t c = 0; c < Cn; ++c)
      if (!keep[c]) E(r, c) = -std::numeric_limits<float>::infinity();
  const double dt = static_cast<double>(hop) / fs * 1e6, df = fs / static_cast<double>(nfft) / 1e6;
  std::vector<FskCand> cands;
  bmat claimed(R, Cn, 0);
  const int64_t gap = std::max<int64_t>(1, py_round_int(60.0 / (static_cast<double>(hop) / fs * 1e6)));
  const bmat st_gap = ndi::ones(static_cast<size_t>(gap), 1);
  const bmat st3 = ndi::ones(3, 3);
  const double t_last = sg.t_us.back();
  const auto& f_mhz = sg.f_mhz;
  for (double thr : {k_sigma, 10 * k_sigma, 100 * k_sigma, 1000 * k_sigma}) {
    bmat mk(R, Cn, 0);
    const float thf = static_cast<float>(thr);
    for (size_t i = 0; i < mk.size(); ++i) mk.data[i] = E.data[i] > thf;
    mk = ndi::binary_closing(mk, st_gap);
    for (size_t i = 0; i < mk.size(); ++i) mk.data[i] = mk.data[i] && std::isfinite(E.data[i]);
    const ndi::Labels lab = ndi::label(mk, st3);
    const auto objs = ndi::find_objects(lab.lab, lab.n);
    for (size_t i = 0; i < objs.size(); ++i) {
      const auto& o = objs[i];
      if (o.empty()) continue;
      const int32_t id = static_cast<int32_t>(i + 1);
      bool cl = false;
      for (size_t r = o.r0; r < o.r1 && !cl; ++r)
        for (size_t c = o.c0; c < o.c1; ++c)
          if (lab.lab(r, c) == id && claimed(r, c)) {
            cl = true;
            break;
          }
      if (cl) continue;
      const double dur = static_cast<double>(o.r1 - o.r0) * dt;
      const double bw = static_cast<double>(o.c1 - o.c0) * df;
      if (dur < min_us || bw > max_bw / 1e6 || dur > 0.9 * t_last) continue;
      const size_t nc = o.c1 - o.c0;
      dvec spec(nc, 0.0);
      for (size_t r = o.r0; r < o.r1; ++r)
        for (size_t c = 0; c < nc; ++c) spec[c] += Z(r, o.c0 + c) * (lab.lab(r, o.c0 + c) == id ? 1.0 : 0.0);
      dvec fw(nc);
      for (size_t c = 0; c < nc; ++c) fw[c] = f_mhz[o.c0 + c] * spec[c];
      const double fc_b = np::pairwise_sum(fw.data(), nc) / std::max(np::pairwise_sum(spec.data(), nc), 1e-12);
      // flanks while the burst is on: in-band power must stand >= 2x above BOTH flanks
      double zacc = 0;
      for (size_t r = o.r0; r < o.r1; ++r) zacc += np::pairwise_sum(Z.row(r) + o.c0, nc);
      const double zin = zacc / static_cast<double>((o.r1 - o.r0) * nc);
      std::vector<double> fl_z;
      const double flo = f_mhz[o.c0], fhi = f_mhz[o.c1 - 1];
      for (int side = 0; side < 2; ++side) {
        std::vector<size_t> cols;
        for (size_t c = 0; c < Cn; ++c) {
          const bool in = side == 0 ? (f_mhz[c] >= flo - 2.0 && f_mhz[c] < flo - 0.4)
                                    : (f_mhz[c] > fhi + 0.4 && f_mhz[c] <= fhi + 2.0);
          if (in && keep[c]) cols.push_back(c);
        }
        if (cols.size() < 3) continue;
        dvec m(cols.size(), 0.0);
        for (size_t r = o.r0; r < o.r1; ++r)
          for (size_t q = 0; q < cols.size(); ++q) m[q] += Z(r, cols[q]);
        for (auto& v : m) v /= static_cast<double>(o.r1 - o.r0);
        fl_z.push_back(np::median(m));
      }
      if (fl_z.empty()) continue;
      const double fmax = *std::max_element(fl_z.begin(), fl_z.end());
      if (fmax - 1 > 0.3 * (zin - 1) && zin < 2.0 * fmax) continue;
      float pe = -std::numeric_limits<float>::infinity();
      for (size_t r = o.r0; r < o.r1; ++r)
        for (size_t c = o.c0; c < o.c1; ++c) {
          if (lab.lab(r, c) == id) claimed(r, c) = 1;
          pe = std::max(pe, E(r, c));
        }
      FskCand cd;
      cd.t0_us = sg.t_us[o.r0] - dt;
      cd.t1_us = sg.t_us[o.r1 - 1] + dt;
      cd.f_lo = (f_mhz[o.c0] - df / 2) * 1e6;
      cd.f_hi = (f_mhz[o.c1 - 1] + df / 2) * 1e6;
      cd.fc_b = fc_b * 1e6;
      cd.zin = zin;
      cd.peak_e = pe;
      cands.push_back(cd);
    }
  }
  return cands;
}

SymRate symbol_rate(const fvec& ifr, double ofs, double r_lo, double r_hi) {
  std::optional<SymRate> best;
  const int64_t sms[3] = {std::max<int64_t>(1, static_cast<int64_t>(ofs / 4e6)),
                          std::max<int64_t>(1, static_cast<int64_t>(ofs / 1e6)),
                          std::max<int64_t>(1, static_cast<int64_t>(ofs / 250e3))};
  for (int64_t sm : sms) {
    const fvec v = ndi::uniform_filter1d(ifr, static_cast<size_t>(sm));
    const float med = np::median(v);
    std::vector<int64_t> tr;
    float prev = 0;
    for (size_t i = 0; i < v.size(); ++i) {
      const float d = v[i] - med;
      const float s = d > 0 ? 1.0f : (d < 0 ? -1.0f : 0.0f);
      if (i > 0 && s - prev != 0) tr.push_back(static_cast<int64_t>(i - 1));
      prev = s;
    }
    if (tr.size() < 8) continue;
    const size_t n = v.size();
    dvec ind(n, 0.0);
    for (auto t : tr) ind[t] = 1.0;
    const double mu = np::pairwise_sum(ind.data(), n) / static_cast<double>(n);
    for (auto& q : ind) q -= mu;
    const cdvec F = fft::rfft(ind.data(), n, 4 * n);
    const dvec ff = fft::rfftfreq(4 * n, 1.0 / ofs);
    const double hi = std::min(r_hi, ofs / (2.5 * static_cast<double>(sm)));
    int64_t kbest = -1;
    double fbest = 0;
    for (size_t i = 0; i < ff.size(); ++i)
      if (ff[i] > r_lo && ff[i] < hi) {
        const double p = np::cabs2(F[i]);
        if (kbest < 0 || p > fbest) {
          fbest = p;
          kbest = static_cast<int64_t>(i);
        }
      }
    if (kbest < 0) continue;
    const double Rr = ff[kbest];
    dvec t(tr.size());
    for (size_t i = 0; i < tr.size(); ++i) t[i] = static_cast<double>(tr[i]) / ofs;
    const dvec lin = np::linspace<double>(-0.03, 0.03, 241);
    dvec conc(241);
    cdvec e(t.size());
    for (size_t g = 0; g < 241; ++g) {
      const double gr = Rr * (1 + lin[g]);
      for (size_t i = 0; i < t.size(); ++i) {
        const double th = (2 * kPi) * (gr * t[i]);
        e[i] = cd(std::cos(th), std::sin(th));
      }
      conc[g] = np::cabs(np::mean(e.data(), e.size()));
    }
    const size_t j = np::argmax(conc);
    if (!best || conc[j] > best->conc) {
      const double gj = Rr * (1 + lin[j]);
      const double c = (2 * kPi) * gj;
      for (size_t i = 0; i < t.size(); ++i) {
        const double th = c * t[i];
        e[i] = cd(std::cos(th), std::sin(th));
      }
      const cd m = np::mean(e.data(), e.size());
      best = SymRate{gj, conc[j], std::atan2(m.imag(), m.real())};
    }
  }
  return best ? *best : SymRate{};
}

int ble_channel_index(double f_hz) {
  const double k = (f_hz - 2402e6) / 2e6;
  const int64_t ki = py_round_int(k);
  if (std::fabs(k - static_cast<double>(ki)) > 0.1 || !(0 <= ki && ki <= 39)) return -1;
  if (ki == 0) return 37;
  if (ki == 12) return 38;
  if (ki == 39) return 39;
  return static_cast<int>(ki < 12 ? ki - 1 : ki - 2);
}

// python _ble_pdu(d, ln, ch, wch): dewhitened, CRC-verified PDU bits -> the decoded dict
Json ble_pdu(const std::vector<uint8_t>& d, int ln, const Json& ch_j, int wch, int ch) {
  int ptype = 0;
  for (int i = 0; i < 4; ++i) ptype |= d[i] << i;
  std::vector<uint8_t> by(ln);
  for (int j = 0; j < ln; ++j) {
    int v = 0;
    for (int i = 0; i < 8; ++i) v |= d[16 + 8 * j + i] << i;
    by[j] = static_cast<uint8_t>(v);
  }
  Json out = Json::object();
  out["aa"] = fmt("%08X", kBleAdvAA);
  out["crc_ok"] = true;
  out["channel_index_rf"] = ch_j;
  out["whitening_channel"] = wch;
  out["channel_mismatch"] = wch != ch;
  out["pdu_type"] = ptype;
  out["pdu_len"] = ln;
  out["pdu_name"] = pdu_name(ptype);
  out["tx_addr_random"] = static_cast<bool>(d[6]);
  if ((ptype == 0 || ptype == 2 || ptype == 4 || ptype == 6) && ln >= 6) {
    std::string addr;
    for (int j = 5; j >= 0; --j) addr += (j < 5 ? ":" : "") + fmt("%02X", by[j]);
    out["adv_addr"] = addr;
    std::string hx;
    for (int j = 6; j < ln; ++j) hx += fmt("%02x", by[j]);
    out["adv_data_hex"] = hx;
    // AD structures; Remote ID (ASTM F3411) if present (python RID.parse_ble_ad)
    const Json ad = si::remote_id::parse_ble_ad(std::vector<uint8_t>(by.begin() + 6, by.begin() + ln));
    if (!ad.empty()) out["ad"] = ad;
  } else if (ptype == 7 && ln >= 1) {
    // extended advertising (ADV_EXT_IND -> AUX_ADV_IND; BT5 Long Range Remote ID); research, PLAN T2b
    std::vector<uint8_t> data;
    Json ext = coded::parse_extended(by, data);
    out["pdu_name"] = (ch == 37 || ch == 38 || ch == 39) ? "ADV_EXT_IND" : "AUX_ADV_IND";
    if (ext.contains("adv_addr")) {
      out["adv_addr"] = ext["adv_addr"];
      ext.erase("adv_addr");
    }
    out["extended"] = ext;
    if (!data.empty()) {
      std::string hx;
      for (uint8_t v : data) hx += fmt("%02x", v);
      out["adv_data_hex"] = hx;
      Json ad = si::remote_id::parse_ble_ad(data);
      if (!ad.empty()) {
        if (ad.contains("remote_id") && ad["remote_id"].is_object())
          ad["remote_id"]["transport"] = "Bluetooth LE extended advertising (BT5)";
        out["ad"] = ad;
      }
    }
  }
  return out;
}

// python ble_coded_decode(y, ofs, ch, aa_bits): LE Coded PHY (research, PLAN T2b)
Json ble_coded_decode(const cfvec& y, double ofs, int ch) {
  if (ch < 0 || std::fabs(ofs / (gfsk::kSps * 1e6) - 1) > 1e-3) return nullptr;
  std::vector<uint8_t> aa_bits(32);
  for (int i = 0; i < 32; ++i) aa_bits[i] = (kBleAdvAA >> i) & 1;
  const cdvec p = gfsk::diff_products(y);
  for (const auto& [i0, nc] : coded::find_coded(p, aa_bits)) {
    const auto r = coded::decode(y, i0, aa_bits);
    if (!r || r->aa_bit_errors > 4 || !r->has_bits) continue;
    const auto d = ble_dewhiten(r->bits.data(), r->bits.size(), ch);
    if (d.size() < 16 + 24) continue;
    int ln = 0;
    for (int i = 0; i < 8; ++i) ln |= d[8 + i] << i;
    const size_t n_pdu = 16 + 8 * static_cast<size_t>(ln);
    if (d.size() < n_pdu + 24 || !ble_crc_ok(d, n_pdu)) continue;
    Json out = ble_pdu(d, ln, Json(ch), ch, ch);
    out["phy"] = r->ci == 0 ? "LE Coded S=8" : "LE Coded S=2";
    out["demod"] = "coherent";
    out["aa_bit_errors"] = r->aa_bit_errors;
    out["mod_index"] = py_round(r->h, 2);
    out["aa_ncorr"] = py_round(nc, 3);
    return out;
  }
  return nullptr;
}

// python ble_coherent_decode(y, ofs, ch, aa_bits): access address by soft correlation, coherent detection
// (gfsk_coherent), CRC with Chase correction on the whitening channel of the RF frequency (research, PLAN L1)
Json ble_coherent_decode(const cfvec& y, double ofs, int ch) {
  if (ch < 0 || std::fabs(ofs / (gfsk::kSps * 1e6) - 1) > 1e-3) return nullptr;
  std::vector<uint8_t> aa_bits(32);
  for (int i = 0; i < 32; ++i) aa_bits[i] = (kBleAdvAA >> i) & 1;
  const cdvec p = gfsk::diff_products(y);
  for (const auto& [i0, nc] : gfsk::find_access_address(p, aa_bits)) {
    const int64_t nb_pre = std::min<int64_t>(8, i0 / gfsk::kSps);
    const int64_t s0 = i0 - nb_pre * gfsk::kSps;
    const auto res = gfsk::demodulate(y, s0, nb_pre + 32 + 16 + 8 * 255 + 24);
    if (!res) continue;
    const auto& bits = res->bits;
    auto aa_errs = [&](int64_t q) -> int64_t {
      if (nb_pre + q < 0) return 99;
      int64_t e = 0;
      for (int64_t i = 0; i < 32; ++i) {
        const int64_t j = nb_pre + q + i;
        if (j >= static_cast<int64_t>(bits.size()) || bits[j] != aa_bits[i]) ++e;
      }
      return e;
    };
    int64_t sh = -1;
    for (int64_t q : {0, 1})
      if (aa_errs(q) < aa_errs(sh)) sh = q;
    const int64_t aa_err = aa_errs(sh);
    if (aa_err > 4) continue;
    const int64_t off = nb_pre + sh + 32;
    if (static_cast<int64_t>(bits.size()) - off < 16 + 24) continue;
    const size_t npdu_avail = bits.size() - static_cast<size_t>(off);
    const auto d = ble_dewhiten(bits.data() + off, npdu_avail, ch);
    int ln = 0;
    for (int i = 0; i < 8; ++i) ln |= d[8 + i] << i;
    const size_t n_pdu = 16 + 8 * static_cast<size_t>(ln);
    if (d.size() < n_pdu + 24) continue;
    const auto r = gfsk::chase(d, n_pdu, res->rel, off);
    if (!r) continue;
    Json out = ble_pdu(r->first, ln, Json(ch), ch, ch);
    out["demod"] = "coherent";
    out["bit_corrections"] = r->second;
    out["aa_bit_errors"] = aa_err;
    out["mod_index"] = py_round(res->h, 2);
    out["aa_ncorr"] = py_round(nc, 3);
    return out;
  }
  return nullptr;
}

Json ble_try_decode(const cfvec& y, double ofs, double f_hz) {
  const int ch = ble_channel_index(f_hz);
  const Json ch_j = ch < 0 ? Json(nullptr) : Json(ch);
  // coherent detector + Chase first (~4-5 dB more sensitive); the hard-decision slicer is the fallback
  if (Json coh = ble_coherent_decode(y, ofs, ch); coh.is_object()) return coh;
  if (Json cod = ble_coded_decode(y, ofs, ch); cod.is_object()) return cod;
  const double sps = ofs / 1e6;
  if (y.size() < 2) return nullptr;
  const size_t n = y.size() - 1;
  fvec ifr(n), pw(n);
  const bool el = elides(n, sizeof(cf));
  for (size_t i = 0; i < n; ++i) {
    ifr[i] = angle_f(y[i + 1], y[i], el);
    pw[i] = np::cabs2(y[i + 1]);
  }
  const float lim = 0.25f * np::percentile(pw, 90);
  std::vector<uint8_t> on(n);
  size_t n_on = 0;
  for (size_t i = 0; i < n; ++i) n_on += (on[i] = pw[i] > lim);
  int aa[32];
  for (int i = 0; i < 32; ++i) aa[i] = ((kBleAdvAA >> i) & 1) ? 1 : -1;
  Json fallback = nullptr;
  const int64_t sms[2] = {std::max<int64_t>(1, static_cast<int64_t>(sps / 2)),
                          std::max<int64_t>(1, static_cast<int64_t>(sps * 0.7))};
  const dvec phs = np::arange<double>(0.0, sps, std::max(sps / 8, 1.0));
  for (int64_t sm : sms) {
    const fvec ifs = ndi::uniform_filter1d(ifr, static_cast<size_t>(sm));
    float thr = 0.0f;
    if (n_on > 16) {
      fvec v;
      for (size_t i = 0; i < n; ++i)
        if (on[i]) v.push_back(ifs[i]);
      const auto [lo, hi] = kmeans2(v);
      thr = (lo + hi) / 2.0f;
    }
    for (double ph : phs) {
      const int64_t nb = static_cast<int64_t>((static_cast<double>(n) - ph) / sps);
      if (nb < 80) return fallback;
      std::vector<uint8_t> bits(static_cast<size_t>(nb));
      for (int64_t i = 0; i < nb; ++i) bits[i] = ifs[static_cast<size_t>(ph + static_cast<double>(i) * sps)] > thr;
      std::vector<int64_t> pos;
      for (int64_t p = 0; p + 32 <= nb; ++p) {
        int c = 0;
        for (int i = 0; i < 32; ++i) c += (bits[p + i] ? 1 : -1) * aa[i];
        if (c >= 30) pos.push_back(p);
      }
      for (size_t pi = 0; pi < pos.size() && pi < 3; ++pi) {
        const int64_t p = pos[pi];
        const uint8_t* pdu = bits.data() + p + 32;
        const size_t np_ = static_cast<size_t>(nb - p - 32);
        if (np_ < 16 + 24) continue;
        std::vector<int> order;
        if (ch >= 0) order.push_back(ch);
        for (int c = 0; c < 40; ++c)
          if (c != ch) order.push_back(c);
        for (int wch : order) {
          const auto h = ble_dewhiten(pdu, 16, wch);
          int ln = 0;
          for (int i = 0; i < 8; ++i) ln |= h[8 + i] << i;
          const size_t n_pdu = 16 + 8 * static_cast<size_t>(ln);
          if (np_ < n_pdu + 24) continue;
          const auto d = ble_dewhiten(pdu, n_pdu + 24, wch);
          if (!ble_crc_ok(d, n_pdu)) continue;
          return ble_pdu(d, ln, ch_j, wch, ch);
        }
      }
      if (!pos.empty() && fallback.is_null()) {
        fallback = Json::object();
        fallback["aa"] = fmt("%08X", kBleAdvAA);
        fallback["crc_ok"] = false;
        fallback["channel_index_rf"] = ch_j;
      }
    }
  }
  return fallback;
}

std::vector<FskMeas> refine_local(const cfvec& x, double fs, double fc, const FskCand& cand) {
  const double noise_win_us = 2000.0;
  auto [Xs, off] = local_fft(x, fs, cand.t0_us, cand.t1_us, noise_win_us + 500.0);
  FskCand c = cand;
  c.t0_us = cand.t0_us - off;
  c.t1_us = cand.t1_us - off;
  auto rs = Xs.empty() ? std::vector<FskMeas>() : fsk_refine(Xs, fs, fc, c, noise_win_us);
  for (auto& r : rs) {
    r.t0_us += off;
    r.t1_us += off;
  }
  return rs;
}

FskLabel fsk_classify(const FskMeas& r, double fc) {
  if (r.esb && !(r.ble.is_object() && !r.ble.empty())) return {"esb_link", "esb_like_2m", 0.6, std::nullopt};
  if (r.ieee802154.is_object()) {
    const bool ok = r.ieee802154["crc_ok"].get<bool>();
    return {"telemetry", ok ? "ieee802154" : "ieee802154_like", ok ? 0.9 : 0.6, std::nullopt};
  }
  const double dur = r.t1_us - r.t0_us;
  if (r.const_env && !*r.const_env) {
    if (r.snr_db >= 12 && dur >= 60) return {"unknown", "nb_nonconst_envelope", 0.4, std::nullopt};
    return {std::nullopt, std::nullopt, 0, std::nullopt};
  }
  const std::optional<double> R = r.sym_rate;
  const double dev = r.has_dev ? r.dev_hz : 0.0;
  if (r.sik.is_object() && r.sik.value("crc_ok", false))
    return {"telemetry", "sik_telemetry", 0.9, r.sik["air_rate_bps"].get<double>()};
  if (r.crossfire.is_object()) return {"rc_fhss", "crossfire_sync", 0.8, crossfire::kRate};
  if (r.adsl.is_object()) {
    if (r.adsl["protocol"].get<std::string>().rfind("OGN", 0) == 0) return {"telemetry", "ogn_like", 0.75, flarm::kChipRate};
    return {"telemetry", "adsl", 0.9, flarm::kChipRate};
  }
  if (r.flarm.is_object()) return {"telemetry", "flarm_like", 0.8, flarm::kChipRate};
  if (r.xn297.is_object()) {
    const bool ok = r.xn297["checksum_ok"].is_boolean() && r.xn297["checksum_ok"].get<bool>();
    if (r.xn297.contains("payload_len")) return {"rc_fhss", "toy_drone", ok ? 0.85 : 0.6, xn297::kRate};
    return {"rc_fhss", ok ? "bayang" : "xn297_like", ok ? 0.85 : 0.6, xn297::kRate};
  }
  if (r.elrs_gfsk.is_object()) return {"rc_fhss", "elrs_gfsk", 0.9, elrs::kGfskRate};
  if (r.radiolink.is_object()) return {"rc_fhss", "radiolink", 0.9, radiolink::kRate};
  if (r.hubsan.is_object()) return {"rc_fhss", "hubsan", 0.85, hubsan::kRate};
  if (r.cc2500_rc.is_object()) {
    const std::string p = r.cc2500_rc["protocol"].get<std::string>();
    if (p.rfind("Graupner", 0) == 0) return {"rc_fhss", "hott", 0.9, hott::kRate};
    if (p.rfind("Hitec", 0) == 0) return {"rc_fhss", "hitec", 0.9, hitec::kRate};
    return {"rc_fhss", "sfhss", 0.85, hitec::kSfRate};
  }
  if (r.afhds2a.is_object()) {
    const bool ok = r.afhds2a["crc_ok"].is_boolean() && r.afhds2a["crc_ok"].get<bool>();
    if (r.afhds2a["protocol"].get<std::string>().rfind("FlySky AFHDS 1", 0) == 0)
      return {"rc_fhss", ok ? "afhds1" : "afhds1_like", ok ? 0.85 : 0.7, afhds2a::kRate};
    return {"rc_fhss", ok ? "afhds2a" : "afhds2a_id", ok ? 0.85 : 0.7, afhds2a::kRate};
  }
  if (r.frsky.is_object()) {
    const std::string var = r.frsky["variant"].get<std::string>();
    const bool ok = r.frsky["crc_ok"].is_boolean() && r.frsky["crc_ok"].get<bool>();
    return {"rc_fhss", "frsky_" + std::string(var.rfind("D16", 0) == 0 ? "d16" : "d8"), ok ? 0.9 : 0.7,
            cc2500::nominal_rate(var)};
  }
  if (r.ble.is_object() && !r.ble.empty()) {
    const bool ok = r.ble.value("crc_ok", false);
    return {"bluetooth", ok ? "ble_adv" : "bluetooth_like", ok ? 0.95 : 0.8, 1e6};
  }
  if (r.bredr.is_object()) return {"bluetooth", "bluetooth_bredr", 0.9, 1e6};   // BCH-verified access code (S7)
  const double conc = r.has_rate ? r.sym_conc : 0.0;
  if (!R || *R == 0 || dev == 0 || conc < 0.5 || r.snr_db < 6 || !(r.const_env && *r.const_env))
    return {"narrowband", "nb_burst", 0.3, std::nullopt};
  if (dur * 1e-6 * *R < 64) return {"narrowband", "nb_burst", 0.3, std::nullopt};
  if (0.85e6 < *R && *R < 1.15e6 && 100e3 < dev && dev < 400e3) return {"bluetooth", "bluetooth_like", 0.6, 1e6};
  const double carson = 2 * (dev + *R / 2);
  if (r.bw99 > 3 * carson + 100e3) return {"narrowband", "fsk_unknown", 0.3, std::nullopt};
  for (const auto& c : fsk_classes()) {
    const std::string lab = c.label;
    const bool sub_ghz = lab == "crossfire_fsk_like" || lab == "sik_like";
    if (sub_ghz != (fc < 1.5e9)) continue;
    if (c.r0 <= *R && *R <= c.r1 && c.d0 <= dev && dev <= c.d1)
      return {"rc_fhss", lab, 0.6, c.nominal > 0 ? std::optional<double>(c.nominal) : std::nullopt};
  }
  return {"narrowband", "fsk_unknown", 0.4, std::nullopt};
}

// python _esb_try(x, fs, fc, r): a 2 Mb/s ESB-like packet in a measured 2.4 GHz burst 0.5-3 MHz wide, sliced
// +-1.25 MHz at 8 samples/bit from an FFT of the burst +-50 us, cropped to the burst -2 / +6 us
std::optional<esb::Packet> esb_try(const cfvec& x, double fs, double fc, const FskMeas& r) {
  const double dur = r.t1_us - r.t0_us;
  if (!(kEsbBandLo <= r.f_center && r.f_center <= kEsbBandHi && kEsbUsLo <= dur && dur <= kEsbUsHi &&
        kEsbBwLo <= r.bw99 && r.bw99 <= kEsbBwHi && r.snr_db >= 6.0))
    return std::nullopt;
  auto [Xs, off] = local_fft(x, fs, r.t0_us, r.t1_us, kEsbMarginUs);
  if (Xs.empty()) return std::nullopt;
  const SliceOut s = slice(Xs, fs, FOff{r.f_center - fc, false}, esb::kHalfBw, esb::kOutFs);
  const int64_t ny = static_cast<int64_t>(s.y.size());
  const int64_t i0 = std::min(ny, std::max<int64_t>(0, static_cast<int64_t>((r.t0_us - off - 2.0) * 1e-6 * s.ofs)));
  const int64_t i1 = std::min(ny, static_cast<int64_t>((r.t1_us - off + 6.0) * 1e-6 * s.ofs));
  auto d = esb::slice_packet(s.y.data() + i0, static_cast<size_t>(std::max<int64_t>(0, i1 - i0)), s.ofs);
  if (d) d->f_carrier = fc + s.f_act + d->f_off_hz;
  return d;
}

// python _zb_try(x, fs, fc, r): an IEEE 802.15.4 frame (research, PLAN T4) in a 2.4 GHz burst on an 802.15.4 channel
Json zb_try(const cfvec& x, double fs, double fc, const FskMeas& r) {
  const auto ch = ieee802154::channel(r.f_center, 300e3);
  if (!ch || !(1.2e6 <= r.bw99 && r.bw99 <= 4e6) || r.t1_us - r.t0_us < 250 || r.snr_db < 3.0) return nullptr;
  const double f0 = 2405e6 + 5e6 * (*ch - 11);
  auto [Xs, off] = local_fft(x, fs, r.t0_us, r.t1_us, kEsbMarginUs);
  if (Xs.empty()) return nullptr;
  const SliceOut s = slice(Xs, fs, FOff{f0 - fc, false}, 1.6e6, 8e6);
  const int64_t ny = static_cast<int64_t>(s.y.size());
  const int64_t i0 = std::min(ny, std::max<int64_t>(0, static_cast<int64_t>((r.t0_us - off - 20.0) * 1e-6 * s.ofs)));
  const int64_t i1 = std::min(ny, static_cast<int64_t>((r.t1_us - off + 20.0) * 1e-6 * s.ofs));
  Json d = ieee802154::demod(cfvec(s.y.begin() + i0, s.y.begin() + std::max(i0, i1)), s.ofs);
  if (d.is_object()) d["channel"] = *ch;
  return d;
}

}  // namespace detail

std::vector<Emission> detect_fsk(const cfvec& x, double fs, double fc) {
  using namespace detail;
  const std::vector<FskCand> cands = fsk_scan(x, fs, fc);
  std::vector<FskMeas> refined;
  for (auto& rs : pmap(cands, [&](const FskCand& c) { return refine_local(x, fs, fc, c); }))
    for (auto& r : rs) refined.push_back(std::move(r));
  // fragments of one long low-SNR burst: merge same-frequency pieces < 100 us apart, re-measure
  std::stable_sort(refined.begin(), refined.end(), [](const FskMeas& a, const FskMeas& b) { return a.t0_us < b.t0_us; });
  struct Merged {
    double f_center, bw99, t0, t1;
    std::vector<FskMeas> parts;
  };
  std::vector<Merged> merged;
  for (auto& r : refined) {
    bool done = false;
    for (auto& m : merged)
      if (std::fabs(m.f_center - r.f_center) < 0.5 * std::max({m.bw99, r.bw99, 100e3}) && r.t0_us - m.t1 < 100 &&
          r.t1_us - m.t0 > -100) {
        m.parts.push_back(r);
        m.t1 = std::max(m.t1, r.t1_us);
        m.t0 = std::min(m.t0, r.t0_us);
        done = true;
        break;
      }
    if (!done) merged.push_back({r.f_center, r.bw99, r.t0_us, r.t1_us, {r}});
  }
  const auto finals = pmap(merged, [&](const Merged& m) {
    if (m.parts.size() == 1) return m.parts;
    double lo = m.parts[0].f_lo, hi = m.parts[0].f_hi;
    for (const auto& p : m.parts) {
      lo = std::min(lo, p.f_lo);
      hi = std::max(hi, p.f_hi);
    }
    FskCand c;
    c.t0_us = m.t0;
    c.t1_us = m.t1;
    c.f_lo = lo;
    c.f_hi = hi;
    c.fc_b = (lo + hi) / 2;
    auto rr = refine_local(x, fs, fc, c);
    // a decode found on one piece (the new PLAN G3-G6 decoders) stays with the merged burst that covers it
    auto carried = [](const FskMeas& q) {
      return q.xn297.is_object() || q.radiolink.is_object() || q.hubsan.is_object() || q.elrs_gfsk.is_object() ||
             q.afhds2a.is_object() || q.cc2500_rc.is_object();
    };
    for (const auto& part : m.parts) {
      auto give = [&](Json FskMeas::*key, bool ok) {
        if (!ok) return;
        for (auto& q : rr)
          if (q.t0_us - 50 <= part.t0_us && part.t0_us <= q.t1_us && !carried(q)) {
            q.*key = part.*key;
            return;
          }
      };
      give(&FskMeas::xn297, part.xn297.is_object() && part.xn297.contains("payload_len"));
      give(&FskMeas::radiolink, part.radiolink.is_object());
      give(&FskMeas::hubsan, part.hubsan.is_object());
      give(&FskMeas::elrs_gfsk, part.elrs_gfsk.is_object());
      give(&FskMeas::afhds2a, part.afhds2a.is_object() &&
                                  part.afhds2a["protocol"].get<std::string>().rfind("FlySky AFHDS 1", 0) == 0);
      give(&FskMeas::cc2500_rc, part.cc2500_rc.is_object());
    }
    return rr.empty() ? m.parts : rr;
  });
  std::vector<Emission> out;
  std::vector<size_t> esb_idx;
  std::vector<std::vector<uint8_t>> esb_bits;
  struct Covered { double f, t0, t1; size_t idx; };
  std::vector<Covered> covered;   // decoded long packets (RadioLink 9.3 ms): later pieces inside their air time join
  for (const auto& group : finals)
    for (const auto& r0 : group) {
      FskMeas r = r0;
      {
        const double f0 = r.has_dev ? r.fsk_center : r.f_center;
        const Covered* host = nullptr;
        for (const auto& c : covered)
          if (std::fabs(c.f - f0) < 100e3 && c.t0 <= r.t0_us && r.t0_us < c.t1) {
            host = &c;
            break;
          }
        if (host) {
          out[host->idx].end_us = std::max(out[host->idx].end_us, r.t1_us);
          continue;
        }
      }
      // energy significance of the refined burst: snr * sqrt(B*T)
      const double sig = std::pow(10.0, r.snr_db / 10) * std::sqrt(std::max(r.bw99, 50e3) * (r.t1_us - r.t0_us) * 1e-6);
      if (r.t1_us - r.t0_us < 20 || sig < 8.0) continue;
      double f = r.has_dev ? r.fsk_center : r.f_center;
      if (gsm_burst(r, f)) continue;
      r.ieee802154 = zb_try(x, fs, fc, r);
      const std::optional<esb::Packet> ep = r.ieee802154.is_object() ? std::nullopt : esb_try(x, fs, fc, r);
      r.esb = ep.has_value();
      const FskLabel cl = fsk_classify(r, fc);
      if (!cl.family) continue;  // weak non-constant-envelope piece: not ours
      Json feats = Json::object();
      feats["modulation"] = (r.has_dev && r.dev_hz != 0) ? Json("gfsk") : Json(nullptr);
      feats["bw_hz"] = r.bw99;
      feats["bw_mhz"] = r.bw99 / 1e6;
      feats["burst_us"] = r.t1_us - r.t0_us;
      feats["symbol_rate"] = r.sym_rate ? Json(*r.sym_rate) : Json(nullptr);
      feats["dev_hz"] = r.has_dev ? Json(r.dev_hz) : Json(nullptr);
      feats["mod_index"] = r.has_dev ? Json(r.mod_index) : Json(nullptr);
      feats["const_env_rho"] = r.rho;
      feats["sym_conc"] = r.has_rate ? Json(r.sym_conc) : Json(nullptr);
      if (cl.nominal && r.sym_rate && *r.sym_rate != 0) {
        feats["clock_ppm"] = (*r.sym_rate / *cl.nominal - 1) * 1e6;
        feats["symbol_us"] = 1e6 / *r.sym_rate;
      }
      if (*cl.family == "bluetooth" && r.ble.is_object() && !r.ble.empty()) {
        const double fn = 2402e6 + 2e6 * static_cast<double>(py_round_int((f - 2402e6) / 2e6));
        feats["cfo_hz"] = f - fn;
        feats["cfo_ppm"] = (f - fn) / fn * 1e6;
      }
      if (*cl.family == "esb_link") {
        f = ep->f_carrier;
        const double fn = 1e6 * static_cast<double>(py_round_int(f / 1e6));   // ESB channels: 2400 + k MHz
        feats["modulation"] = "gfsk";
        feats["symbol_rate"] = esb::kRate;
        feats["symbol_us"] = 1e6 / esb::kRate;
        feats["dev_hz"] = ep->dev_hz;
        feats["mod_index"] = 2 * ep->dev_hz / esb::kRate;
        feats["sym_conc"] = nullptr;
        feats["preamble_bits"] = ep->preamble_bits;
        feats["cfo_hz"] = f - fn;
        feats["cfo_ppm"] = (f - fn) / fn * 1e6;
        esb_idx.push_back(out.size());
        esb_bits.push_back(ep->bits);
      }
      if (r.frame.is_object() && *cl.family != "esb_link") feats["frame"] = r.frame;
      Emission e;
      e.family = *cl.family;
      e.label = *cl.label;
      e.start_us = r.t0_us;
      e.end_us = r.t1_us;
      e.f_lo_mhz = std::min(r.f_lo, f - r.bw99 / 2) / 1e6;
      e.f_hi_mhz = std::max(r.f_hi, f + r.bw99 / 2) / 1e6;
      e.center_mhz = f / 1e6;
      e.snr_db = r.snr_db;
      e.confidence = cl.conf;
      e.decoded = (r.ble.is_object() && !r.ble.empty()) ? r.ble : r.bredr.is_object() ? r.bredr : (r.sik.is_object() && !r.sik.empty()) ? r.sik : r.crossfire.is_object() ? r.crossfire : r.frsky.is_object() ? r.frsky : r.afhds2a.is_object() ? r.afhds2a : r.xn297.is_object() ? r.xn297 : r.ieee802154.is_object() ? r.ieee802154 : r.flarm.is_object() ? r.flarm
                  : r.adsl.is_object() ? r.adsl : r.radiolink.is_object() ? r.radiolink : r.hubsan.is_object() ? r.hubsan
                  : r.elrs_gfsk.is_object() ? r.elrs_gfsk : r.cc2500_rc;
      e.features = std::move(feats);
      out.push_back(std::move(e));
      if (r.radiolink.is_object()) {
        const double n_bits = 8.0 * (4 + 4 + r.radiolink["length"].get<int>() + 3);
        covered.push_back({f, r.t0_us, r.t0_us + n_bits / radiolink::kRate * 1e6 + 200.0, out.size() - 1});
      }
      if (r.adsl.is_object()) {   // ADS-L / OGN (FLARM PHY): the packet's chips after the sync
        const bool og = r.adsl["protocol"].get<std::string>().rfind("OGN", 0) == 0;
        const double n_chips = og ? 2.0 * 8 * (4 + 26) : 2.0 * 8 * (r.adsl.value("length", 26) + 5);
        covered.push_back({f, r.t0_us, r.t0_us + n_chips / flarm::kChipRate * 1e6 + 200.0, out.size() - 1});
      }
    }
  if (!esb_idx.empty()) {   // one CRC search per link (packets sharing a sync word)
    std::vector<double> starts;
    for (size_t i : esb_idx) starts.push_back(out[i].start_us);
    const std::vector<Json> dec = esb::verify_links(esb_bits, starts);
    for (size_t j = 0; j < esb_idx.size(); ++j) {
      out[esb_idx[j]].decoded = dec[j];
      if (dec[j]["crc_ok"].get<bool>()) out[esb_idx[j]].confidence = 0.85;
    }
  }
  return out;
}

}  // namespace si::rx::rc_link
