// Pipeline plugin: analog FM video transmitters (PAL/NTSC), demodulated to a picture
// (python: receivers/analog_vtx/plugin.py).
#include <algorithm>
#include <cmath>
#include <exception>
#include <memory>
#include <set>

#include "core/drone_models.hpp"
#include "core/bands_india.hpp"
#include "core/pyfmt.hpp"
#include "decode/emission_adapter.hpp"
#include "receivers/analog_vtx/detector.hpp"

namespace si::rx::analog_vtx {

namespace {

class AnalogVtxDecoder : public EmissionReceiver {
 public:
  std::string name() const override { return "analog_vtx"; }
  std::string package() const override { return "analog_vtx"; }
  std::string protocol() const override { return "analog FM video (PAL/NTSC)"; }
  std::string key() const override { return "analog_vtx"; }
  std::string method() const override { return "receivers.analog_vtx (research receiver) + FM discriminator"; }
  bool owns_band() const override { return true; }
  std::vector<std::string> steps() const override {
    return {"~1 us x ~1 MHz time-frequency power; fold every column on the PAL (15625 Hz) and NTSC "
            "(15734.26 Hz) line period over a +-500 ppm grid: the sync tip parks the carrier at one "
            "frequency for 4.7 us per line -> line-locked z-score (works to ~-12 dB in-band SNR)",
            "score >= 10 -> emitter: line rate refined (fold sharpness + drift fit), sync-tip frequency, "
            "occupied band, field timing (lines per field from the vertical sync), standard channel table",
            "FM discriminator on the occupied band: deviation, polarity, audio subcarrier",
            "picture: discriminator output sampled on the measured line clock, active lines stacked "
            "-> <report>_sig<id>_video.png"};
  }

  std::vector<Emission> run_detector(DecodeContext& ctx) const override {
    // python: schema/requirements.json analog_vtx.min_sample_rate_hz (20 MHz): a 6+ MHz FM video
    const double rate = std::min(ctx.cap.fs, ctx.cap.bandwidth_hz.value_or(ctx.cap.fs));
    if (rate < 20e6) return {};
    return detect(ctx.cap.iq(), ctx.cap.fs, ctx.cap.fc);
  }
  bool in_range(double lo, double hi) const override { return in_ranges(lo, hi, bands::analog_vtx_ranges()); }
  std::string group_key(const Emission& e) const override { return py_float_repr(e.center_mhz); }

  DecodeResult summarise(const std::vector<Emission>& em, const Signal*) const override {
    const Emission& e = em[0];
    const Json& F = e.features;
    const std::string std_name = F.contains("standard") ? F["standard"].get<std::string>() : e.label;
    static const char* keys[] = {"standard", "line_rate_hz", "line_rate_ppm", "line_rate_ppm_sigma",
                                 "lines_per_field", "field_rate_hz", "scan", "sync_tip_mhz", "center_mhz",
                                 "deviation_mhz", "deviation_source", "polarity", "audio_subcarrier_mhz",
                                 "channel", "channel_candidates", "nominal_mhz", "cfo_hz", "bw_mhz", "score",
                                 "color", "color_standard", "chroma_fsc_hz", "burst_swing_deg", "chroma_coherence",
                                 "color_note", "chroma_candidates", "chroma_gate_note", "channel_note",
                                 "secam_db_hz", "secam_dr_hz", "audio_subcarrier_hz", "audio_dev_rms_khz",
                                 "audio_peak_hz", "audio_peak_db", "audio_flatness"};
    Json fields = Json::object();
    for (const char* k : keys)
      if (F.contains(k) && !F[k].is_null()) fields[k] = F[k];
    {  // model level (research, PLAN M8; python _model)
      std::vector<std::string> cands, mev;
      if (F.contains("channel") && F["channel"].is_string() && !F["channel"].get<std::string>().empty())
        mev.push_back("channel " + F["channel"].get<std::string>());
      if (std::fabs(e.center_mhz - 5584.0) <= 2.0) {
        cands.push_back("TBS Unify Pro in pit mode (very low power)");
        mev.push_back("5584 MHz is the TBS Unify Pro pit-mode frequency");
      }
      if (std::fabs(e.center_mhz - 1258.0) <= 3.0 || std::fabs(e.center_mhz - 1280.0) <= 3.0) {
        cands.push_back("2-channel 1.2 / 1.3 GHz VTX (VTX-1G3SE class)");
        cands.push_back("9-channel 1.2 / 1.3 GHz VTX");
        mev.push_back("1258 / 1280 MHz: the two channels of locked 1.2 / 1.3 GHz VTXs");
      }
      const bool has_a = F.contains("audio_subcarrier_mhz") && F["audio_subcarrier_mhz"].is_number() &&
                         F["audio_subcarrier_mhz"].get<double>() != 0.0;
      mev.push_back(has_a ? "audio subcarrier " + py_str(F["audio_subcarrier_mhz"]) + " MHz: a VTX with a microphone input"
                          : std::string("no audio subcarrier"));
      fields["model"] = drone_models::model_block(
          "analog FPV VTX", cands, mev, "inferred",
          cands.empty() ? "analog VTX models are not separable on air (power, control and camera standard are not visible)" : "");
    }
    auto num = [&](const char* k) { return F.contains(k) && F[k].is_number() ? F[k].get<double>() : 0.0; };
    std::string ch;
    if (F.contains("channel") && F["channel"].is_string() && !F["channel"].get<std::string>().empty())
      ch = ", channel " + F["channel"].get<std::string>();
    const std::string lpf = F.contains("lines_per_field") ? py_str(F["lines_per_field"]) : "?";
    // python STD_NAMES.get(e.label, STD_NAMES.get(std, "FM video"))
    auto std_label = [](const std::string& s) -> std::string {
      return s == "pal" ? "PAL" : s == "ntsc" ? "NTSC" : s == "pal_m" ? "PAL-M" : s == "pal_n" ? "PAL-N"
             : s == "secam" ? "SECAM" : "";
    };
    std::string sl = std_label(e.label);
    if (sl.empty()) sl = std_label(std_name);
    if (sl.empty()) sl = "FM video";
    const std::string label = "Analog FPV video (" + sl + ")";
    const std::string ev = fmt("line-locked sync tips at %.2f Hz (%+.0f ppm), fold score %.1f (threshold 10)",
                               num("line_rate_hz"), num("line_rate_ppm"), num("score")) +
                           ch + "; " + lpf + " lines/field";
    return result("success", method(), steps(), fields, "", label, e.confidence, "FM (analog composite video)", ev,
                  std::make_pair(e.f_lo_mhz * 1e6, e.f_hi_mhz * 1e6));
  }

  void finish(Signal& sig, const std::vector<Emission>& em, DecodeContext& ctx) const override {
    std::shared_ptr<Image> img;
    try {
      img = video_frame(ctx.cap.iq(), ctx.cap.fs, ctx.cap.fc, em[0]);
    } catch (const std::exception& e) {   // the picture is a bonus, the detection stands
      sig.decoded.fields["video_frame_error"] = std::string("Error: ") + std::string(e.what()).substr(0, 190);
      return;
    }
    if (img) {
      sig.video_frame = img;
      sig.decoded.fields["video_frame_lines"] = static_cast<int64_t>(img->rows);
    }
    sig.decoded.reason =
        "analog video demodulated to a picture (no CRC exists for analog video; "
        "line/field timing locks are the verification)";
  }

  // python AnalogVtxDecoder.arbitrate: weak receiver detections inside the VTX's FM spectrum
  // (max(6 MHz, 0.2 x its bandwidth) beyond its band) are pieces of it.
  std::vector<Signal> arbitrate(std::vector<Signal> signals, DecodeContext&) const override {
    static const std::set<std::string> weak = {"rc_link", "digital_vtx"};
    auto decoder = [](const Signal& s) {
      return s.features.contains("decoder") && s.features["decoder"].is_string() ? s.features["decoder"].get<std::string>()
                                                                                 : std::string();
    };
    std::vector<size_t> vtx;
    for (size_t i = 0; i < signals.size(); ++i)
      if (decoder(signals[i]) == name()) vtx.push_back(i);
    if (vtx.empty()) return signals;
    std::vector<size_t> keep;
    for (size_t i = 0; i < signals.size(); ++i) {
      const Signal& s = signals[i];
      std::optional<size_t> host;
      if (weak.count(decoder(s)) && s.decoded.status != "success" && s.label.find("Herelink") == std::string::npos) {
        for (size_t v : vtx) {
          const Signal& h = signals[v];
          const double m = std::max(6e6, 0.2 * h.bw());
          if (h.f0 - m <= s.fc() && s.fc() <= h.f1 + m) {
            host = v;
            break;
          }
        }
      }
      if (!host) {
        keep.push_back(i);
        continue;
      }
      Json& det = signals[*host].detection;
      if (!det.contains("explained_receiver_detections")) det["explained_receiver_detections"] = Json::array();
      det["explained_receiver_detections"].push_back(
          {{"decoder", decoder(s)},
           {"label", s.label},
           {"f_low_hz", py_round_int(s.f0)},
           {"f_high_hz", py_round_int(s.f1)},
           {"note", "inside the analog VTX's FM spectrum: attributed to it"}});
    }
    std::vector<Signal> out;
    for (size_t i : keep) out.push_back(std::move(signals[i]));
    return out;
  }
};

}  // namespace

void register_plugins(si::PluginList& out) { out.push_back(std::make_unique<AnalogVtxDecoder>()); }

}  // namespace si::rx::analog_vtx
