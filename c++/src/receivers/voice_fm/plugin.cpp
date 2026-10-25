// Pipeline plugin: narrowband FM voice (walkie-talkies / PMR), demodulated to audio
// (python: receivers/voice_fm/plugin.py).
#include <algorithm>
#include <exception>
#include <memory>

#include "core/pyfmt.hpp"
#include "decode/emission_adapter.hpp"
#include "receivers/voice_fm/audio.hpp"
#include "receivers/voice_fm/detector.hpp"

namespace si::rx::voice_fm {

namespace {

bool truthy(const Json& j, const char* k) {
  if (!j.contains(k)) return false;
  const Json& v = j[k];
  if (v.is_null()) return false;
  if (v.is_boolean()) return v.get<bool>();
  if (v.is_number()) return v.get<double>() != 0;
  if (v.is_string()) return !v.get<std::string>().empty();
  return !v.empty();
}

class VoiceFmDecoder : public EmissionReceiver {
 public:
  std::string name() const override { return "voice_fm"; }
  std::string package() const override { return "voice_fm"; }
  std::string protocol() const override { return "analog narrowband FM voice"; }
  std::string key() const override { return "voice_fm"; }
  std::string method() const override {
    return "receivers.voice_fm: FM discriminator + CTCSS / DCS (Golay) + channel plan";
  }
  bool owns_band() const override { return true; }
  std::vector<std::string> steps() const override {
    return {"1 kHz-bin spectrum of the capture: narrow peaks >= 10 dB over the floor and their "
            "surroundings in 118-174 / 225-520 / 806-870 MHz",
            "per peak, the +-10 kHz channel at 48 kSps: keyed time (5 ms frames >= 10 dB over the floor), "
            "constant envelope, FM discriminator -> carrier, peak deviation (<= 6 kHz), voice-band activity; "
            "2/4-level FSK data and unmodulated carriers are rejected; without a sub-audible code the envelope "
            "must vary no more than noise explains (two-tone beats / clock lines are not FM)",
            "the receiver's own clock lines (pipeline: k x fs, reference-locked) are skipped",
            "sub-audible signalling on the longest keyed segment: DCS 134.4 bit/s NRZ, every 23-bit window "
            "checked against Golay(23,12) and the 100 marker (code + aliases), else a CTCSS tone matched to "
            "the EIA list",
            "clipping products / IQ image of a stronger carrier (k x its LO offset, |k| x its deviation) "
            "dropped and listed on it",
            "channel plan (PMR446, LPD433, FRS) with the receiver ppm when the pipeline has one",
            "audio: channel -> discriminator -> 300-3400 Hz -> 750 us de-emphasis, carrier squelch "
            "-> <name>_voice<MHz>_audio.wav"};
  }

  std::vector<Emission> run_detector(DecodeContext& ctx) const override {
    // python: cal = ctx.cache.get("ppm") or ctx.cache.get("sweep_ppm")
    const Json* cal = ctx.cache.get<Json>("ppm");
    if (!cal || cal->is_null() || cal->empty()) cal = ctx.cache.get<Json>("sweep_ppm");
    double ppm = 0.0;
    if (cal && cal->is_object() && cal->contains("ppm") && (*cal)["ppm"].is_number()) ppm = (*cal)["ppm"].get<double>();
    dvec excl;
    if (const Json* ln = ctx.cache.get<Json>("receiver_lines"))
      for (const auto& l : *ln) excl.push_back(l["freq_hz"].get<double>());
    return detect(ctx.ch, ppm, excl);
  }
  bool in_range(double lo, double hi) const override { return in_ranges(lo, hi, voice_ranges()); }
  std::string group_key(const Emission& e) const override { return py_float_repr(py_round(e.center_mhz, 3)); }

  DecodeResult summarise(const std::vector<Emission>& em, const Signal*) const override {
    size_t bi = 0;
    for (size_t i = 1; i < em.size(); ++i)
      if (em[i].confidence > em[bi].confidence) bi = i;
    const Emission& e = em[bi];
    const Json& F = e.features;
    if (e.label == "am_voice") return summarise_am(e, em);
    if (e.label == "aprs" || e.label == "ax25_packet") return summarise_packets(e, em);
    static const char* keys[] = {
        "center_mhz", "deviation_hz", "channel_bw_khz", "service", "channel", "nominal_mhz", "offset_hz",
        "offset_ppm", "receiver_ppm_applied", "type", "dcs_code", "dcs_aliases", "dcs_standard",
        "dcs_words_verified", "dcs_bit_repeat_agreement", "dcs_deviation_hz", "dcs_note", "ctcss_hz",
        "ctcss_measured_hz", "ctcss_snr_db", "ctcss_standard", "ctcss_deviation_hz", "keyed_s",
        "keyed_fraction", "voice_activity", "voice_rms_hz", "voice_kurtosis", "envelope_cv", "snr_db",
        "segments_s", "suppressed_artifacts"};
    Json fields = Json::object();
    for (const char* k : keys)
      if (F.contains(k) && !F[k].is_null()) fields[k] = F[k];
    fields["audio"] = "nbfm";   // the viewer offers 'Listen' on this
    std::string code;
    if (truthy(F, "dcs_code")) code = F["dcs_code"].get<std::string>();
    else if (truthy(F, "ctcss_hz")) code = "CTCSS " + py_str(F["ctcss_hz"]) + " Hz";
    const std::string where = truthy(F, "service")
                                  ? F["service"].get<std::string>() + " ch " + py_str(F["channel"])
                                  : fmt("%.4f MHz", e.center_mhz);
    const bool verified = (truthy(F, "dcs_code") && F.contains("dcs_words_verified") &&
                           F["dcs_words_verified"].get<double>() >= 2) ||
                          truthy(F, "ctcss_hz");
    const bool speech = F.contains("voice_activity") && F["voice_activity"].get<double>() >= 0.1;
    const std::string label = (!code.empty() || speech)
                                  ? "Analog FM voice (" + where + (code.empty() ? "" : ", " + code) + ")"
                                  : "Analog FM voice-like (" + where + ")";
    const std::string reason =
        std::string("FM voice demodulated to audio; ") +
        (verified ? "sub-audible squelch code verified (DCS Golay / standard CTCSS tone)" : "no sub-audible code to verify") +
        (speech ? "" : "; little speech-band activity in this capture (carrier keyed, quiet mic?)");
    std::string ev;
    if (F.contains("evidence"))
      for (size_t i = 0; i < F["evidence"].size(); ++i) ev += (i ? "; " : "") + F["evidence"][i].get<std::string>();
    double lo = em[0].f_lo_mhz, hi = em[0].f_hi_mhz;
    for (const auto& x : em) {
      lo = std::min(lo, x.f_lo_mhz);
      hi = std::max(hi, x.f_hi_mhz);
    }
    return result(verified ? "success" : "partial", method(), steps(), fields, reason, label, e.confidence,
                  "FM (analog narrowband voice)", ev, std::make_pair(lo * 1e6, hi * 1e6));
  }

  // AM voice (research, signal gaps S13; python _summarise_am)
  DecodeResult summarise_am(const Emission& e, const std::vector<Emission>& em) const {
    const Json& F = e.features;
    static const char* keys[] = {"center_mhz", "modulation_rms", "channel_bw_khz", "service", "channel_mhz",
                                 "offset_hz", "emergency", "keyed_s", "keyed_fraction", "voice_activity",
                                 "phase_rms_hz", "envelope_cv", "snr_db", "segments_s", "suppressed_artifacts"};
    Json fields = Json::object();
    for (const char* k : keys)
      if (F.contains(k) && !F[k].is_null()) fields[k] = F[k];
    fields["audio"] = "am";
    std::string where;
    if (truthy(F, "emergency")) where = F["emergency"].get<std::string>();
    else if (truthy(F, "service")) where = F["service"].get<std::string>() + ", " + py_str(F["channel_mhz"]) + " MHz";
    else where = fmt("%.4f MHz", e.center_mhz);
    std::string ev;
    if (F.contains("evidence"))
      for (size_t i = 0; i < F["evidence"].size(); ++i) ev += (i ? "; " : "") + F["evidence"][i].get<std::string>();
    double lo = em[0].f_lo_mhz, hi = em[0].f_hi_mhz;
    for (const auto& x : em) {
      lo = std::min(lo, x.f_lo_mhz);
      hi = std::max(hi, x.f_hi_mhz);
    }
    const std::vector<std::string> st = {
        steps()[0],
        "per peak, the +-10 kHz channel at 48 kSps: keyed time; AM: the envelope's "
        "300-3400 Hz rms relative to the carrier (>= 4 % and 3x what noise gives) in >= 10 % of 20 ms "
        "frames, while the discriminator's speech-band rms stays at the noise level (phase quiet)",
        "aeronautical band and 25 kHz channel (VHF civil 118-137 MHz, UHF military 225-400 MHz); "
        "121.5 / 243.0 MHz emergency",
        "audio: channel -> envelope -> 300-3400 Hz, gain from the carrier level, carrier squelch "
        "-> <name>_voice<MHz>_audio.wav"};
    return result("partial", "receivers.voice_fm: AM envelope detector (research)", st, fields,
                  "AM voice demodulated to audio; AM carries no code to verify", "Analog AM voice (" + where + ")",
                  e.confidence, "AM (analog voice, full carrier)", ev, std::make_pair(lo * 1e6, hi * 1e6));
  }

  // AX.25 / APRS packet radio (research, signal gaps S19; python _summarise_packets)
  DecodeResult summarise_packets(const Emission& e, const std::vector<Emission>& em) const {
    const Json& F = e.features;
    static const char* keys[] = {"center_mhz", "deviation_hz", "channel_bw_khz", "keyed_s", "snr_db", "stations",
                                 "packets", "segments_s"};
    Json fields = Json::object();
    for (const char* k : keys)
      if (F.contains(k) && !F[k].is_null()) fields[k] = F[k];
    const std::string what = e.label == "aprs" ? "APRS" : "AX.25 packet radio";
    const Json& st = F["stations"];
    std::string who;
    for (size_t i = 0; i < st.size() && i < 3; ++i) who += (i ? ", " : "") + st[i].get<std::string>();
    if (st.size() > 3) who += ", ...";
    std::string ev;
    for (size_t i = 0; i < F["evidence"].size(); ++i) ev += (i ? "; " : "") + F["evidence"][i].get<std::string>();
    double lo = em[0].f_lo_mhz, hi = em[0].f_hi_mhz;
    for (const auto& x : em) {
      lo = std::min(lo, x.f_lo_mhz);
      hi = std::max(hi, x.f_hi_mhz);
    }
    const std::vector<std::string> st2 = {
        steps()[0], "per peak, the +-10 kHz channel at 48 kSps: FM discriminator",
        "Bell 202 AFSK: mark 1200 / space 2200 Hz tone energy over one bit, 8 bit phases, NRZI decoding, HDLC "
        "frames with a valid CRC-16/X.25, AX.25 UI addresses, APRS position / status / message fields"};
    return result("success", "receivers.voice_fm: Bell 202 AFSK + AX.25 HDLC (research)", st2, fields,
                  "every packet verified by its CRC-16/X.25", what + " (" + who + ")", e.confidence,
                  "AFSK 1200 bit/s (Bell 202) on narrowband FM", ev, std::make_pair(lo * 1e6, hi * 1e6));
  }

  void finish(Signal& sig, const std::vector<Emission>&, DecodeContext& ctx) const override {
    if (sig.decoded.fields.contains("packets")) return;   // packet radio: decoded data, no audio
    std::optional<AudioChunk> a;
    try {
      a = demodulate(ctx.cap.iq(), ctx.cap.fs, ctx.cap.fc, sig.decoded.fields, 16000, &ctx.ch);
    } catch (const std::exception& e) {   // audio is a bonus, the detection stands
      sig.decoded.fields["audio_error"] = std::string("Error: ") + std::string(e.what()).substr(0, 190);
      return;
    }
    if (a && !a->samples.empty()) {
      sig.audio = std::make_shared<fvec>(std::move(a->samples));
      sig.audio_fs = a->fs;
      sig.decoded.fields["audio_s"] = py_round(static_cast<double>(sig.audio->size()) / a->fs, 3);
      sig.decoded.fields["audio_squelch_open"] = py_round(a->meta["squelch_open"].get<double>(), 3);
    }
  }
};

}  // namespace

void register_plugins(si::PluginList& out) { out.push_back(std::make_unique<VoiceFmDecoder>()); }

}  // namespace si::rx::voice_fm
