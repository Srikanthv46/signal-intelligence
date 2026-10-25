// receivers/voice_fm against the Python receiver (tests/vectors/voice_fm/gen_voice_fm.py):
// DCS/Golay tables, candidates, per-carrier measurements, Emission features, one-shot audio and
// the streaming VoiceDemodulator, on the real PMR446 capture (DCS D053N) and a synthetic
// speech + CTCSS transmitter. WAV bytes are checked against Python's wave module separately.
#include <doctest/doctest.h>

#include <cmath>
#include <cstring>
#include <fstream>
#include <string>

#include "core/dsp.hpp"
#include "receivers/voice_fm/audio.hpp"
#include "receivers/voice_fm/channels.hpp"
#include "receivers/voice_fm/detector.hpp"
#include "receivers/voice_fm/subaudio.hpp"

using namespace si;
namespace vf = si::rx::voice_fm;

namespace {

const std::string kDir = std::string(SI_TESTVEC) + "/voice_fm/";

template <class T>
std::vector<T> load_bin(const std::string& path) {
  std::ifstream f(path, std::ios::binary);
  REQUIRE_MESSAGE(f.good(), "missing " << path << " (run tests/vectors/voice_fm/gen_voice_fm.py)");
  std::vector<char> b((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
  std::vector<T> v(b.size() / sizeof(T));
  std::memcpy(v.data(), b.data(), v.size() * sizeof(T));
  return v;
}

Json load_json(const std::string& name) {
  std::ifstream f(kDir + name);
  REQUIRE(f.good());
  return Json::parse(f);
}

bool num_ok(const Json& a, const Json& b, double rel) {
  const double x = a.get<double>(), y = b.get<double>();
  const std::string s = a.dump();
  int dec = 15;
  if (s.find('e') == std::string::npos && s.find('E') == std::string::npos) {
    const auto p = s.find('.');
    dec = p == std::string::npos ? 0 : static_cast<int>(s.size() - p - 1);
  }
  const double tol = std::max(rel * std::fabs(x), a.is_number_float() ? std::pow(10.0, -dec) : 0.0);
  return std::fabs(x - y) <= tol * (1 + 1e-9);
}

void same(const Json& a, const Json& b, const std::string& path, double rel = 1e-6) {
  CAPTURE(path);
  if (a.is_object()) {
    REQUIRE(b.is_object());
    CHECK(a.size() == b.size());
    for (auto it = a.begin(); it != a.end(); ++it) {
      CHECK_MESSAGE(b.contains(it.key()), "missing " << it.key());
      if (b.contains(it.key())) same(it.value(), b[it.key()], path + "." + it.key(), rel);
    }
  } else if (a.is_array()) {
    REQUIRE(b.is_array());
    REQUIRE(a.size() == b.size());
    for (size_t i = 0; i < a.size(); ++i) same(a[i], b[i], path + "[" + std::to_string(i) + "]", rel);
  } else if (a.is_number() && !a.is_boolean()) {
    REQUIRE(b.is_number());
    CHECK_MESSAGE(num_ok(a, b, rel), a.dump() << " != " << b.dump());
  } else {
    CHECK_MESSAGE(a == b, a.dump() << " != " << b.dump());
  }
}

double max_abs_diff(const fvec& a, const std::vector<float>& b) {
  REQUIRE(a.size() == b.size());
  double m = 0;
  for (size_t i = 0; i < a.size(); ++i) m = std::max(m, static_cast<double>(std::fabs(a[i] - b[i])));
  return m;
}

void check(const std::string& tag, const cfvec& x) {
  CAPTURE(tag);
  const Json ref = load_json(tag + ".json");
  const double fs = ref["fs"].get<double>(), fc = ref["fc"].get<double>();
  const Channelizer ch(x, fs, fc);
  const auto cands = vf::candidates(ch);
  REQUIRE(cands.size() == ref["candidates"].size());
  for (size_t i = 0; i < cands.size(); ++i) {
    const Json& r = ref["candidates"][i];
    CHECK(cands[i].f_hz == doctest::Approx(r["f_hz"].get<double>()).epsilon(1e-14));
    CHECK(cands[i].peak_db == r["peak_db"].get<double>());
    CHECK(cands[i].width_hz == doctest::Approx(r["width_hz"].get<double>()).epsilon(1e-12));
    CHECK(cands[i].noise_density == doctest::Approx(r["noise_density"].get<double>()).epsilon(1e-14));
    const auto m = vf::measure(ch, cands[i]);
    const Json& rm = ref["measure"][i];
    REQUIRE(m.has_value() == !rm.is_null());
    if (!m) continue;
    Json cm = {{"f_hz", m->f_hz}, {"keyed_s", m->keyed_s}, {"duration_s", m->duration_s},
               {"envelope_cv", m->envelope_cv}, {"deviation_hz", m->deviation_hz}, {"voice_rms_hz", m->voice_rms_hz},
               {"voice_kurtosis", m->voice_kurtosis}, {"voice_activity", m->voice_activity},
               {"subaudible", m->sub ? *m->sub : Json(nullptr)},
               {"snr_db", m->snr_db ? Json(*m->snr_db) : Json(nullptr)}, {"peak_db", m->peak_db}};
    same(rm, cm, "measure[" + std::to_string(i) + "]", 1e-9);
  }
  const auto em = vf::detect(ch, 0.0);
  REQUIRE(em.size() == ref["emissions"].size());
  for (size_t i = 0; i < em.size(); ++i) {
    const Json& r = ref["emissions"][i];
    CHECK(em[i].label == r["label"].get<std::string>());
    CHECK(em[i].confidence == r["confidence"].get<double>());
    CHECK(em[i].start_us == doctest::Approx(r["start_us"].get<double>()));
    CHECK(em[i].end_us == doctest::Approx(r["end_us"].get<double>()));
    CHECK(em[i].f_lo_mhz == doctest::Approx(r["f_lo_mhz"].get<double>()).epsilon(1e-14));
    CHECK(em[i].f_hi_mhz == doctest::Approx(r["f_hi_mhz"].get<double>()).epsilon(1e-14));
    same(r["features"], em[i].features, "features");
    same(r["decoded"], em[i].decoded, "decoded");
  }
  if (em.empty()) return;

  // one-shot audio
  const auto a = vf::demodulate(x, fs, fc, ref["fields"]);
  REQUIRE(a.has_value());
  CHECK(a->fs == ref["audio"]["fs"].get<double>());
  same(ref["audio"]["meta"], a->meta, "audio.meta");
  CHECK(max_abs_diff(a->samples, load_bin<float>(kDir + tag + "_audio.f32")) < 1e-5);

  // streaming: 50000-sample blocks, filter states carried across blocks
  vf::VoiceDemodulator demod = vf::VoiceDemodulator::from_fields(ref["fields"]);
  fvec all;
  for (size_t k = 0, b = 0; k < x.size(); k += 50000, ++b) {
    const cfvec blk(x.begin() + static_cast<std::ptrdiff_t>(k),
                    x.begin() + static_cast<std::ptrdiff_t>(std::min(x.size(), k + 50000)));
    const auto c = demod.process(blk, fs, fc);
    const Json& rm = ref["stream"][b];
    REQUIRE(c.has_value() == !rm.is_null());
    if (!c) continue;
    same(rm, c->meta, "stream.meta");
    all.insert(all.end(), c->samples.begin(), c->samples.end());
  }
  CHECK(max_abs_diff(all, load_bin<float>(kDir + tag + "_stream.f32")) < 1e-5);
}

}  // namespace

TEST_CASE("voice_fm: DCS Golay tables match python") {
  const Json ref = load_json("dcs.json");
  const auto& codes = vf::dcs_codes();
  REQUIRE(codes.size() == ref["codes"].size());
  for (size_t i = 0; i < codes.size(); ++i) {
    CHECK(codes[i] == ref["codes"][i].get<std::string>());
    CHECK(vf::dcs_encode(codes[i]) == ref["words"][i].get<uint32_t>());
  }
  const auto al = vf::dcs_aliases(vf::dcs_encode("023"));
  REQUIRE(al.size() == ref["aliases_023"].size());
  for (size_t i = 0; i < al.size(); ++i) {
    CHECK(al[i].first == ref["aliases_023"][i][0].get<std::string>());
    CHECK(al[i].second == ref["aliases_023"][i][1].get<std::string>());
  }
  const uint32_t w[] = {0, 1, 0x7FFFFF, 0x123456, 0x400C75};
  for (size_t i = 0; i < 5; ++i) CHECK(vf::golay_syndrome(w[i]) == ref["syndromes"][i].get<uint32_t>());
  const auto ch = vf::match_channel(446.091001e6);
  REQUIRE(ch.has_value());
  CHECK((*ch)["channel"] == 8);
  CHECK((*ch)["offset_hz"] == -2749);
}

TEST_CASE("voice_fm: real PMR446 capture (DCS D053N) matches python") {
  check("real", load_bin<cf>(std::string(SI_TESTDATA) + "/voice/real_pmr446_ch8_dcs.cf32"));
}

TEST_CASE("voice_fm: synthetic speech + CTCSS 88.5 Hz matches python") { check("synth", load_bin<cf>(kDir + "synth.cf32")); }

TEST_CASE("voice_fm: WAV header") {
  const auto b = vf::wav_bytes(fvec{0.f, 1.f, -1.f}, 16000.0);
  REQUIRE(b.size() == 44 + 6);
  CHECK(std::string(b.begin(), b.begin() + 4) == "RIFF");
  CHECK(b[44 + 2] == 0xFF);
  CHECK(b[44 + 3] == 0x7F);   // 32767
}
