// Analog FM video (PAL / NTSC composite) -> a picture, one field per block
// (python: receivers/analog_vtx/video.py, FieldDemodulator / VideoFrame).
//
// A streaming demodulator for a VTX that detect() already found: build it from the
// detection, then feed it consecutive (or non-contiguous) IQ blocks; every block that holds
// a vertical sync gives one field: luma, and colour (rgb) when the VTX sends a colour burst.
//
//   auto demod = FieldDemodulator::from_fields(signal["decoded"]["fields"]);   // report fields
//   auto demod = FieldDemodulator::from_emission(e);                           // detect() output
//   if (auto fr = demod.process(iq, fs, fc)) use(fr->image, fr->rows, fr->cols, fr->meta);
//
// Same steps as Python: slice + power-weighted FM discriminator (the low-passed lag-1 product,
// complex64, as core.dsp.Channelizer), 31-tap low-pass; the line clock as an H-PLL over the block
// (sync phase from the fold at the detector's rate, every line's sync edge placed by the 4.7 us
// box and timed by a short step filter, robust straight-line fit, weighted against the detector's
// rate); raster on the line clock, black = back porch, white = 99.5th percentile, vertical-sync
// field selection; sync-tip noise, and a 3x3 median above NOISY. Colour (python
// FieldDemodulator.chroma): the plain discriminator low-passed to 6 MHz on a 1024-per-line grid,
// mixed with the line-locked subcarrier, 95-tap 1 MHz chroma filter, burst gate tried at three
// offsets, burst phase averaged over same-parity lines (PAL: two-line reference, V switch, delay
// line), sync-depth levels (from the picture's white when the sync is compressed) and burst ACC,
// field-level and noise colour killers. The FIR filters are direct convolutions in float64 where
// scipy uses oaconvolve; outputs agree to rounding.
#pragma once

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "core/dsp.hpp"
#include "core/model.hpp"
#include "core/types.hpp"

namespace si::rx::analog_vtx {

constexpr int kVideoWidth = 720;

struct VideoStandard {
  double active_us0, active_len_us, porch_us0, porch_us1;
  int skip, lines;
  double fl;
};
const VideoStandard& video_standard(const std::string& std_name);   // "pal" | "ntsc"

struct VideoFrame {
  size_t rows = 0, cols = 0;
  std::vector<uint8_t> image;   // row-major, 0 = black (luma, chroma included)
  std::vector<uint8_t> rgb;     // row-major [rows, cols, 3]; empty: no colour burst (or color off)
  Json meta = Json::object();   // line_rate_hz, line_rate_ppm, vsync, lines, noise, denoise,
                                // sync_jitter_us, standard, center_mhz, color, color_killer,
                                // burst_ratio, burst_field, burst_gate_us, sync_compressed
  std::shared_ptr<Image> to_image() const;
};

class FieldDemodulator {
 public:
  FieldDemodulator(double center_hz, double line_rate_hz, const std::string& standard = "ntsc",
                   std::optional<double> sync_hz = std::nullopt, std::optional<double> bw_hz = std::nullopt,
                   const std::string& polarity = "sync_low", bool color = true);
  static FieldDemodulator from_fields(const Json& fields);     // report decoded.fields (analog_vtx plugin)
  static FieldDemodulator from_emission(const Emission& e);    // detect() output

  bool in_band(double fs, double fc) const;
  // one block -> one luma field, or nullopt (out of band, too short, no picture). ch: an
  // existing Channelizer of this block (optional).
  std::optional<VideoFrame> process(const cfvec& iq, double fs, double fc, const Channelizer* ch = nullptr);

  // the stages, public for tests: the weighted discriminator low-passed to 4.8 MHz (python
  // process() before raster()); raw: the plain discriminator (colour), before any low-pass
  fvec discriminate(const cfvec& iq, double fs, double fc, double& ofs, const Channelizer* ch = nullptr) const;
  fvec discriminate_raw(const cfvec& iq, double fs, double fc, double& ofs, const Channelizer* ch = nullptr) const;
  std::optional<VideoFrame> raster(const fvec& d, double ofs);
  // colour of the field raster() just selected, from the discriminator low-passed to 6 MHz:
  // fills fr.rgb (or leaves it empty) and fr.meta color / burst_ratio
  void chroma(const fvec& d6, double ofs, VideoFrame& fr) const;

  bool color() const { return color_; }
  void set_color(bool on) { color_ = on; }

  const std::string& standard() const { return std_; }
  double rate() const { return rate_; }
  double center_hz() const { return center_hz_; }

 private:
  std::string std_;
  VideoStandard p_;
  double center_hz_, sync_hz_, half_bw_, rate0_, rate_, polarity_;
  bool color_ = true;
  struct Timing {
    double t0, T, sync, black, white;
    std::vector<int64_t> rows;
    double noise;                  // sync-tip noise / (white - black)
  };
  std::optional<Timing> timing_;   // set by raster(), used by chroma()
  struct EdgeFit {
    double t0, T, se, jitter_us;
  };
  // the slice of the VTX band (python process(): y, ofs)
  Channelizer::Slice slice(const cfvec& iq, double fs, double fc, const Channelizer* ch) const;
  std::optional<EdgeFit> edge_fit(const fvec& d, double ofs, double t0, double T) const;
  void prior_blend(double& t0, double& T, double se, double dur) const;
  double sync_noise(const fvec& d, double ofs, double t0, double T, int64_t n_lines) const;
};

}  // namespace si::rx::analog_vtx
