// receivers/analog_vtx against the Python receiver (tests/vectors/analog_vtx/gen_analog_vtx.py):
// STFT power and fold bit-identical, scan profile, every Emission feature, the picture and the
// streaming FieldDemodulator, on testdata/Analog_20260804-161815 and a DJI capture (no VTX).
#include <doctest/doctest.h>

#include <cmath>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>

#include "core/npcompat.hpp"
#include "io/sigmf.hpp"
#include "receivers/analog_vtx/detector.hpp"
#include "receivers/analog_vtx/png.hpp"
#include "receivers/analog_vtx/video.hpp"

using namespace si;
namespace av = si::rx::analog_vtx;

namespace {

const std::string kDir = std::string(SI_TESTVEC) + "/analog_vtx/";

template <class T>
std::vector<T> load(const std::string& name) {
  std::ifstream f(kDir + name, std::ios::binary);
  REQUIRE_MESSAGE(f.good(), "missing vector " << name << " (run tests/vectors/analog_vtx/gen_analog_vtx.py)");
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

// a number as Python printed it may differ by one unit in its last decimal (parity rule)
bool num_ok(const Json& a, const Json& b) {
  const double x = a.get<double>(), y = b.get<double>();
  const std::string s = a.dump();
  int dec = 15;
  if (s.find('e') == std::string::npos && s.find('E') == std::string::npos) {
    const auto p = s.find('.');
    dec = p == std::string::npos ? 0 : static_cast<int>(s.size() - p - 1);
  }
  const double tol = std::max(1e-6 * std::fabs(x), a.is_number_float() ? std::pow(10.0, -dec) : 0.0);
  return std::fabs(x - y) <= tol * (1 + 1e-9);
}

void same(const Json& a, const Json& b, const std::string& path) {
  CAPTURE(path);
  if (a.is_object()) {
    REQUIRE(b.is_object());
    CHECK(a.size() == b.size());
    for (auto it = a.begin(); it != a.end(); ++it) {
      CHECK_MESSAGE(b.contains(it.key()), "missing " << it.key());
      if (b.contains(it.key())) same(it.value(), b[it.key()], path + "." + it.key());
    }
  } else if (a.is_array()) {
    REQUIRE(b.is_array());
    REQUIRE(a.size() == b.size());
    for (size_t i = 0; i < a.size(); ++i) same(a[i], b[i], path + "[" + std::to_string(i) + "]");
  } else if (a.is_number() && !a.is_boolean()) {
    REQUIRE(b.is_number());
    CHECK_MESSAGE(num_ok(a, b), a.dump() << " != " << b.dump());
  } else {
    CHECK_MESSAGE(a == b, a.dump() << " != " << b.dump());
  }
}

Capture capture(const std::string& tag) {
  const std::string td = SI_TESTDATA;
  return tag == "analog" ? load_sigmf(td + "/Analog_20260804-161815")
                         : load_sigmf(td + "/DJI_Mini_2_20251119-165139");
}

void check_capture(const std::string& tag) {
  CAPTURE(tag);
  const Json ref = load_json(tag + ".json");
  const Capture cap = capture(tag);
  const cfvec& x = cap.iq();
  const size_t n = std::min(x.size(), av::kMaxSamples);
  const av::TfPower tf = av::tf_power(x.data(), n, cap.fs);
  REQUIRE(tf.t.size() == ref["n_frames"].get<size_t>());
  const auto P = load<float>(tag + "_P1000.f32");
  const auto q = load<float>(tag + "_q1000.f32");
  size_t bad = 0;
  for (size_t i = 0; i < P.size(); ++i) bad += P[i] != tf.P.data[i] || q[i] != tf.q.data[i];
  CHECK(bad == 0);

  const av::Scanner sc(tf);
  const dvec prof = sc.profile();
  const auto pref = load<double>(tag + "_profile.f64");
  REQUIRE(prof.size() == pref.size());
  double worst = 0;
  for (size_t i = 0; i < prof.size(); ++i) worst = std::max(worst, std::fabs(prof[i] - pref[i]) / std::max(1.0, std::fabs(pref[i])));
  CHECK(worst < 1e-9);
  const av::ScanBest best = sc.best(nullptr);
  CHECK(best.std_name == ref["best"]["std"].get<std::string>());
  CHECK(best.b0 == ref["best"]["b0"].get<int>());
  CHECK(best.k0 == ref["best"]["k0"].get<int>());
  CHECK(best.ppm == doctest::Approx(ref["best"]["ppm"].get<double>()).epsilon(1e-12));
  CHECK(best.score == doctest::Approx(ref["best"]["score"].get<double>()).epsilon(1e-9));
  const av::Fold F = av::fold(tf.q, tf.t, best.rate);
  const auto Fref = load<double>(tag + "_fold_best.f64");
  REQUIRE(F.F.data.size() == Fref.size());
  bad = 0;
  for (size_t i = 0; i < Fref.size(); ++i) bad += F.F.data[i] != Fref[i];
  CHECK(bad == 0);

  const std::vector<Emission> em = av::detect(x, cap.fs, cap.fc);
  REQUIRE(em.size() == ref["emissions"].size());
  for (size_t i = 0; i < em.size(); ++i) {
    const Json& r = ref["emissions"][i];
    CHECK(em[i].label == r["label"].get<std::string>());
    CHECK(em[i].f_lo_mhz == r["f_lo_mhz"].get<double>());
    CHECK(em[i].f_hi_mhz == r["f_hi_mhz"].get<double>());
    CHECK(em[i].center_mhz == r["center_mhz"].get<double>());
    CHECK(*em[i].snr_db == r["snr_db"].get<double>());
    CHECK(em[i].confidence == doctest::Approx(r["confidence"].get<double>()).epsilon(1e-12));
    same(r["features"], em[i].features, "features");
  }
  if (em.empty()) return;

  // the picture (report): identical pixels
  const auto img = av::video_frame(x, cap.fs, cap.fc, em[0]);
  REQUIRE(img);
  CHECK(img->rows == ref["video_shape"][0].get<size_t>());
  CHECK(img->cols == ref["video_shape"][1].get<size_t>());
  const auto vref = load<uint8_t>(tag + "_video.u8");
  REQUIRE(vref.size() == img->pixels.size());
  size_t diff = 0;
  for (size_t i = 0; i < vref.size(); ++i) diff += static_cast<float>(vref[i]) != img->pixels[i];
  CHECK(diff == 0);
  const auto png = av::png_gray8(*img);
  CHECK(png.size() > img->pixels.size());

  // streaming fields: same field timing; pixels within 1 grey level (31-tap low-pass in float64)
  av::FieldDemodulator demod = av::FieldDemodulator::from_emission(em[0]);
  const size_t starts[2] = {0, 500000};
  for (int k = 0; k < 2; ++k) {
    CAPTURE(k);
    const cfvec blk(x.begin() + static_cast<std::ptrdiff_t>(starts[k]), x.end());
    const auto fr = demod.process(blk, cap.fs, cap.fc);
    const Json& rf = ref["fields"][k];
    REQUIRE(fr.has_value() == !rf.is_null());
    if (!fr) continue;
    same(rf["meta"], fr->meta, "meta");
    const auto fref = load<uint8_t>(tag + "_field" + std::to_string(k) + ".u8");
    REQUIRE(fref.size() == fr->image.size());
    size_t off = 0, far = 0;
    for (size_t i = 0; i < fref.size(); ++i) {
      const int d = std::abs(static_cast<int>(fref[i]) - static_cast<int>(fr->image[i]));
      off += d != 0;
      far += d > 1;
    }
    CHECK(far == 0);
    CHECK(off < fref.size() / 1000);
    // colour: present exactly when python has it; within 2 levels (FIR filters in float64)
    REQUIRE(fr->rgb.empty() == !rf["rgb"].get<bool>());
    if (fr->rgb.empty()) continue;
    const auto cref = load<uint8_t>(tag + "_field" + std::to_string(k) + "_rgb.u8");
    REQUIRE(cref.size() == fr->rgb.size());
    REQUIRE(cref.size() == fr->rows * fr->cols * 3);
    size_t coff = 0, cfar = 0;
    for (size_t i = 0; i < cref.size(); ++i) {
      const int d = std::abs(static_cast<int>(cref[i]) - static_cast<int>(fr->rgb[i]));
      coff += d != 0;
      cfar += d > 2;
    }
    CHECK(cfar == 0);
    CHECK(coff < cref.size() / 100);
  }
  // colour off: no rgb, luma unchanged
  av::FieldDemodulator mono = av::FieldDemodulator::from_emission(em[0]);
  mono.set_color(false);
  const auto fm = mono.process(x, cap.fs, cap.fc);
  if (fm) CHECK(fm->rgb.empty());
}

}  // namespace

TEST_CASE("analog_vtx: NTSC VTX capture matches python") { check_capture("analog"); }
TEST_CASE("analog_vtx: DJI capture, no false VTX") { check_capture("dji"); }

TEST_CASE("analog_vtx: channel table") {
  const auto c = av::match_channel(5806.4);
  REQUIRE(!c.empty());
  CHECK(c[0].first == "R5");   // RaceBand first among equals
  CHECK(av::match_channel(5015.7).empty());
  CHECK(av::match_channel(3409.0)[0].first == "3G3A5");  // 3.3 GHz groups (PLAN G9)
  CHECK(av::match_channel(3171.0)[0].first == "3G3B1");
}

TEST_CASE("analog_vtx: colour standard on synthetic composite (NTSC / PAL-M / PAL / PAL-N / monochrome / SECAM)") {
  std::ifstream jf(std::string(SI_TESTVEC) + "/analog_vtx/chroma.json");
  REQUIRE(jf);
  const Json ref = Json::parse(jf);
  for (const auto& c : ref["cases"]) {
    std::ifstream f(std::string(SI_TESTVEC) + "/analog_vtx/" + c["file"].get<std::string>(), std::ios::binary | std::ios::ate);
    REQUIRE(f);
    const size_t n = static_cast<size_t>(f.tellg()) / sizeof(float);
    std::vector<float> v(n);
    f.seekg(0);
    f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(n * sizeof(float)));
    const dvec d(v.begin(), v.end());
    const Json got = si::rx::analog_vtx::chroma_standard(d, ref["fs"], c["rate"], c["t0"], c["line_std"]);
    const Json& want = c["out"];
    INFO(c["std"].get<std::string>());
    CHECK(got["color"] == want["color"]);
    CHECK(got.value("color_standard", Json(nullptr)) == want.value("color_standard", Json(nullptr)));
    for (auto it = want["chroma_candidates"].begin(); it != want["chroma_candidates"].end(); ++it)
      for (const char* k : {"coherence", "prominence_db", "swing_deg", "fsc_hz"})
        CHECK(std::fabs(got["chroma_candidates"][it.key()][k].get<double>() - it.value()[k].get<double>()) <= 0.11);
    if (want.contains("chroma_fsc_hz")) CHECK(std::fabs(got["chroma_fsc_hz"].get<double>() - want["chroma_fsc_hz"].get<double>()) <= 0.11);
    CHECK(got.contains("secam_test") == want.contains("secam_test"));
    if (want.contains("secam_test") && got.contains("secam_test")) {
      const Json &g = got["secam_test"], &w = want["secam_test"];
      CHECK(std::fabs(g["alternation"].get<double>() - w["alternation"].get<double>()) <= 0.0011);
      CHECK(std::fabs(g["t_stat"].get<double>() - w["t_stat"].get<double>()) <= 0.11);
      CHECK(std::fabs(g["prominence_db"].get<double>() - w["prominence_db"].get<double>()) <= 0.11);
      CHECK(g.value("db_hz", Json(nullptr)) == w.value("db_hz", Json(nullptr)));
      CHECK(g.value("dr_hz", Json(nullptr)) == w.value("dr_hz", Json(nullptr)));
    }
    for (const char* k : {"secam_db_hz", "secam_dr_hz"})
      CHECK(got.value(k, Json(nullptr)) == want.value(k, Json(nullptr)));
  }
}

TEST_CASE("analog_vtx: audio subcarrier features == python (research A5)") {
  std::ifstream jf(std::string(SI_TESTVEC) + "/analog_vtx/audio.json");
  REQUIRE(jf);
  const Json ref = Json::parse(jf);
  std::ifstream f(std::string(SI_TESTVEC) + "/analog_vtx/audio.f32", std::ios::binary | std::ios::ate);
  REQUIRE(f);
  const size_t n = static_cast<size_t>(f.tellg()) / sizeof(float);
  std::vector<float> v(n);
  f.seekg(0);
  f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(n * sizeof(float)));
  const dvec d(v.begin(), v.end());
  const Json got = si::rx::analog_vtx::audio_features(d, ref["fs"], ref["fa"]);
  const Json& want = ref["out"];
  REQUIRE(got.size() == want.size());
  CHECK(std::fabs(got["audio_subcarrier_hz"].get<double>() - want["audio_subcarrier_hz"].get<double>()) <= 1.0);
  CHECK(std::fabs(got["audio_dev_rms_khz"].get<double>() - want["audio_dev_rms_khz"].get<double>()) <= 0.011);
  CHECK(std::fabs(got["audio_peak_hz"].get<double>() - want["audio_peak_hz"].get<double>()) <= 0.11);
  CHECK(std::fabs(got["audio_peak_db"].get<double>() - want["audio_peak_db"].get<double>()) <= 0.11);
  CHECK(std::fabs(got["audio_flatness"].get<double>() - want["audio_flatness"].get<double>()) <= 0.0011);
}
