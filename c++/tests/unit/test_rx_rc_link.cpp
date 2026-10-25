// rc_link receiver against the Python reference (tests/vectors/rc_link/gen_rc_link.py).
#include <doctest/doctest.h>

#include <cmath>
#include <fstream>
#include <string>

#include "io/sigmf.hpp"
#include "receivers/rc_link/detector.hpp"
#include "receivers/rc_link/dsp_util.hpp"
#include "receivers/rc_link/fsk.hpp"
#include "receivers/rc_link/label_check.hpp"
#include "receivers/rc_link/lora.hpp"
#include "receivers/rc_link/lora_phy.hpp"
#include "core/fft.hpp"

using namespace si;
using namespace si::rx::rc_link;

namespace {

const std::string kVec = std::string(SI_TESTVEC) + "/rc_link/";

template <class T>
std::vector<T> load_bin(const std::string& path) {
  std::ifstream f(path, std::ios::binary | std::ios::ate);
  REQUIRE_MESSAGE(f.good(), "missing " << path << " (run tests/vectors/rc_link/gen_rc_link.py)");
  const auto n = static_cast<size_t>(f.tellg()) / sizeof(T);
  std::vector<T> v(n);
  f.seekg(0);
  f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(n * sizeof(T)));
  return v;
}

const Json& ref() {
  static const Json j = [] {
    std::ifstream f(kVec + "ref.json");
    REQUIRE_MESSAGE(f.good(), "missing ref.json (run tests/vectors/rc_link/gen_rc_link.py)");
    return Json::parse(f);
  }();
  return j;
}

const cfvec& lora_x() {
  static const cfvec x = load_bin<cf>(std::string(SI_TESTDATA) + "/lora/real_434MHz_sf9.cf32");
  return x;
}
const cfvec& lora_X() {
  static const cfvec X = fft::fft(lora_x());
  return X;
}

bool close(double a, double b, double rel = 1e-6) {
  if (std::isnan(a) || std::isnan(b)) return std::isnan(a) && std::isnan(b);
  return std::fabs(a - b) <= rel * std::max(1.0, std::fabs(a));
}

size_t mismatches(const cfvec& a, const cfvec& b) {
  if (a.size() != b.size()) return std::max(a.size(), b.size());
  size_t bad = 0;
  for (size_t i = 0; i < a.size(); ++i) bad += !(a[i] == b[i]);
  return bad;
}

// numbers / strings / lists in a Python walk or emission dict vs the C++ value
void check_num(const Json& py, double v, const char* what, double rel = 1e-9) {
  CAPTURE(what);
  if (py.is_null()) {
    CHECK(std::isnan(v));
  } else {
    CHECK(close(py.get<double>(), v, rel));
  }
}

}  // namespace

TEST_CASE("rc_link: LoRa building blocks bit-identical to Python") {
  const Json& R = ref();
  const double fs = 500e3;
  CHECK(mismatches(detail::upchirp(512), load_bin<cf>(kVec + "upchirp_512.bin")) == 0);
  auto s = detail::slice(lora_X(), fs, {87500.0, false}, 62500.0, 125e3);
  CHECK(mismatches(s.y, load_bin<cf>(kVec + "slice_plain.bin")) == 0);
  CHECK(s.ofs == R["slice"]["ofs"].get<double>());
  CHECK(s.f_act == R["slice"]["f_act"].get<double>());
  auto s2 = detail::slice(lora_X(), fs, {static_cast<float>(87512.34), true}, 62500.0, 125e3, 0.37);
  CHECK(mismatches(s2.y, load_bin<cf>(kVec + "slice_f32_tadv.bin")) == 0);
  auto st = detail::streams(lora_X(), fs, {87500.0, false}, 125e3, 512, 1, 0.37);
  CHECK(mismatches(st.z_up, load_bin<cf>(kVec + "streams_zup.bin")) == 0);
  for (const auto& q : R["sym_spec"]) {
    auto p = detail::sym_spec(st.z_up, q["start"].get<int64_t>(), 512);
    REQUIRE(p);
    CHECK(p->pk == static_cast<float>(q["pk"].get<double>()));
    CHECK(p->r == static_cast<float>(q["r"].get<double>()));
    CHECK(p->s2 == static_cast<float>(q["s2"].get<double>()));
  }
  const Json& cm = R["circ_mean"];
  CHECK(detail::circ_mean(cm["v"].get<dvec>(), 512, cm["w"].get<dvec>()) == cm["out"].get<double>());
}

TEST_CASE("rc_link: LoRa scan, refine and walk on the real SF9 packet") {
  const Json& R = ref();
  const double fs = 500e3, fc = 433.92e6;
  auto raw = detail::lora_scan(lora_X(), fs, 125e3, {7, 8, 9, 10, 11, 12});
  CHECK(raw.size() == R["lora_scan_n"].get<size_t>());
  auto cands = detail::nms_candidates(raw, 125e3);
  REQUIRE(cands.size() == R["lora_nms"].size());
  for (size_t i = 0; i < cands.size(); ++i) {
    const Json& c = R["lora_nms"][i];
    CAPTURE(i);
    CHECK(cands[i].f_sb == c["f_sb"].get<double>());
    CHECK(cands[i].sf == c["sf"].get<int>());
    CHECK(cands[i].pol == c["pol"].get<int>());
    CHECK(cands[i].win == c["win"].get<int64_t>());
    CHECK(cands[i].ratio == c["ratio"].get<double>());
    CHECK(cands[i].bin == c["bin"].get<int64_t>());
    auto w = detail::lora_refine(lora_X(), fs, fc, cands[i]);
    const Json& rw = R["lora_refine"][i];
    REQUIRE(w.has_value() == !rw.is_null());
    if (!w) continue;
    CHECK(w->partial == rw["partial"].get<bool>());
    CHECK(w->prem.n_pre_est == rw["n_pre_est"].get<int64_t>());
    CHECK(w->prem.pre_start == rw["pre_start"].get<int64_t>());
    CHECK(w->prem.k_u == static_cast<float>(rw["k_u"].get<double>()));
    if (w->partial) continue;
    const Json& pw = rw["walk"];
    const auto& cw = w->walk;
    CHECK(cw.n_pre == pw["n_pre"].get<int64_t>());
    CHECK(cw.n_pay == pw["n_pay"].get<int64_t>());
    CHECK(cw.sync_vals == pw["sync_vals"].get<std::vector<int64_t>>());
    CHECK(cw.pay_vals == pw["pay_vals"].get<std::vector<int64_t>>());
    CHECK(cw.all_vals == pw["all_vals"].get<std::vector<int64_t>>());
    check_num(pw["start_samp"], cw.start_samp, "start_samp");
    check_num(pw["f_carrier"], cw.f_carrier, "f_carrier", 1e-12);
    check_num(pw["snr_db"], cw.snr_db, "snr_db");
    check_num(pw["z"], cw.z, "z");
    check_num(pw["z_sfd"], cw.z_sfd, "z_sfd");
    check_num(pw["pre_ratio"], cw.pre_ratio, "pre_ratio");
    check_num(pw["sfd_ratio"], cw.sfd_ratio, "sfd_ratio");
    // weighted LS through BLAS in numpy: agrees to rounding only
    check_num(pw["slope"], cw.slope, "slope", 1e-6);
  }
}

TEST_CASE("rc_link: walk_packet on the Python walk input") {
  const Json& W = ref()["walk_direct"];
  const cfvec zu = load_bin<cf>(kVec + "walk_zup.bin"), zd = load_bin<cf>(kVec + "walk_zdn.bin");
  auto w = detail::walk_packet(zu, zd, W["b"].get<double>(), W["N"].get<int64_t>(), W["cfo"].get<double>(),
                               W["pol"].get<int>());
  REQUIRE(w);
  const Json& o = W["out"];
  CHECK(w->n_pay == o["n_pay"].get<int64_t>());
  CHECK(w->all_vals == o["all_vals"].get<std::vector<int64_t>>());
  check_num(o["snr_db"], w->snr_db, "snr_db");
  check_num(o["z"], w->z, "z");
  check_num(o["slope"], w->slope, "slope", 1e-6);
  check_num(o["slope_err"], w->slope_err, "slope_err", 1e-6);
}

TEST_CASE("rc_link: SF5 implicit-header decode (python decode_sf5_implicit)") {
  // on-air symbols of a live SX126x packet at 890.03 MHz (testdata/live/LIVE_LORA_SF5_890MHz)
  const std::vector<int64_t> real = {1,  1,  11, 25, 11, 29, 27, 1,  15, 31, 6,  27, 11, 27, 0,  31, 14, 28,
                                     15, 11, 19, 26, 4,  9,  27, 18, 23, 12, 7,  7,  13, 22, 21, 19, 2,  30,
                                     27, 23, 13, 30, 25, 3,  6,  22, 3,  15, 31, 29, 7,  29, 9,  29, 1,  5,
                                     31, 2,  5,  8,  1,  19, 17, 1,  1,  16, 26, 26, 31, 2,  1,  7,  26, 16,
                                     0,  9,  16, 3,  7,  28, 2,  27, 1,  4,  17, 21, 2,  2,  2,  2,  2};
  auto d = lora_phy::decode_sf5_implicit(real, 89);
  REQUIRE(d);
  CHECK((*d)["payload_hex"] == "e883000b3e000032000000aedbd71e14e000e0e80300007001002000000000421f0008000200");
  CHECK((*d)["length"] == 38);
  CHECK((*d)["cr"] == "4/5");
  CHECK((*d)["n_symbols"] == 85);
  CHECK((*d)["crc_ok"] == true);
  CHECK((*d)["cr_ambiguous"] == false);
  auto bad = real;  // two symbol errors in one block: uncorrectable -> nothing
  bad[40] = (bad[40] + 7) % 32;
  bad[41] = (bad[41] + 3) % 32;
  CHECK(!lora_phy::decode_sf5_implicit(bad, 89));
  // 1-byte payload 0x57 at CR 4/6 (synth.lora_sf5_implicit_symbols): its last block is zero nibbles,
  // a codeword at every CR; the walk ended it at 10 symbols
  const std::vector<int64_t> one = {1, 1, 9, 3, 3, 9, 9, 15, 25, 17, 1, 1, 1, 1, 1, 1, 28, 11, 2, 31, 17, 14, 9, 14};
  auto o = lora_phy::decode_sf5_implicit(one, 10);
  REQUIRE(o);
  CHECK((*o)["payload_hex"] == "57");
  CHECK((*o)["cr"] == "4/5");
  CHECK((*o)["n_symbols"] == 15);
  CHECK((*o)["cr_ambiguous"] == true);
}

TEST_CASE("rc_link: LoRa emissions and explicit-header payload") {
  const Json& R = ref();
  auto ems = detect_lora(lora_X(), 500e3, 433.92e6);
  REQUIRE(ems.size() == R["lora_emissions"].size());
  for (size_t i = 0; i < ems.size(); ++i) {
    const Json& p = R["lora_emissions"][i];
    CHECK(ems[i].label == p["label"].get<std::string>());
    CHECK(ems[i].decoded == p["decoded"]);
    check_num(p["center_mhz"], ems[i].center_mhz, "center_mhz", 1e-12);
  }
  CHECK(lora_phy::whitening_sequence(16) == R["whitening"].get<std::vector<uint8_t>>());
  const std::string s = "123456789";
  CHECK(lora_phy::crc16(reinterpret_cast<const uint8_t*>(s.data()), s.size()) == R["crc16"].get<int>());
  for (const auto& d : R["lora_phy"]) CHECK(d["decoded"]["payload_hex"] == "a7f1d92a82c8d8fe000000e1000287013439");
}

TEST_CASE("rc_link: FSK scan, refinement, BLE and the frequency-label check (Phantom 4)") {
  const Json& R = ref();
  const Capture cap = load_sigmf(std::string(SI_TESTDATA) + "/DJI_Phantom_4_20251119-182450/D_20251119-182450.sigmf-meta");
  const cfvec& x = cap.iq();
  auto cands = detail::fsk_scan(x, cap.fs, cap.fc);
  REQUIRE(cands.size() == R["fsk_cands"].size());
  for (size_t i = 0; i < cands.size(); ++i) {
    const Json& c = R["fsk_cands"][i];
    CAPTURE(i);
    CHECK(cands[i].t0_us == c["t0_us"].get<double>());
    CHECK(cands[i].f_lo == c["f_lo"].get<double>());
    CHECK(cands[i].fc_b == doctest::Approx(c["fc_b"].get<double>()).epsilon(1e-12));
    CHECK(cands[i].zin == doctest::Approx(c["zin"].get<double>()).epsilon(1e-12));
    CHECK(cands[i].peak_e == c["peak_e"].get<double>());
  }
  size_t k = 0;
  for (size_t i = 0; i < 40 && i < cands.size(); ++i)
    for (const auto& r : detail::refine_local(x, cap.fs, cap.fc, cands[i])) {
      REQUIRE(k < R["fsk_meas"].size());
      const Json& p = R["fsk_meas"][k++];
      CAPTURE(k);
      check_num(p["t0_us"], r.t0_us, "t0_us", 1e-12);
      check_num(p["t1_us"], r.t1_us, "t1_us", 1e-12);
      check_num(p["bw99"], r.bw99, "bw99", 1e-12);
      check_num(p["snr_db"], r.snr_db, "snr_db");
      check_num(p["rho"], r.rho, "rho");
      if (p.contains("sym_rate")) {
        check_num(p["sym_rate"], r.sym_rate ? *r.sym_rate : NAN, "sym_rate");
        check_num(p["sym_conc"], r.sym_conc, "sym_conc");
      }
      if (p.contains("dev_hz")) {
        check_num(p["dev_hz"], r.dev_hz, "dev_hz");
        check_num(p["fsk_center"], r.fsk_center, "fsk_center", 1e-12);
      }
      CHECK(r.ble == (p.contains("ble") ? p["ble"] : Json(nullptr)));
    }
  CHECK(k == R["fsk_meas"].size());
  auto fe = detect_fsk(x, cap.fs, cap.fc);
  REQUIRE(fe.size() == R["fsk_emissions"].size());
  for (size_t i = 0; i < fe.size(); ++i) {
    CHECK(fe[i].label == R["fsk_emissions"][i]["label"].get<std::string>());
    CHECK(fe[i].decoded == R["fsk_emissions"][i]["decoded"]);
  }
  auto lc = check_frequency_label(fe);
  REQUIRE(lc);
  CHECK(*lc == R["label_check"]);
}
