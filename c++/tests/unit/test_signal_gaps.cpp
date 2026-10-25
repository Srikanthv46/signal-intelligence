// Batch S (signal gaps, research/signal_gaps/PLAN.md) vs python (tests/vectors/signal_gaps/gen_signal_gaps.py):
// one TEST_CASE group per item, each on its <item>.json reference and binary inputs.
#include <doctest/doctest.h>

#include <algorithm>
#include <fstream>
#include <limits>
#include <string>

#include "analyze.hpp"
#include "core/deinterleave.hpp"
#include "core/hdlc.hpp"
#include "core/pyfmt.hpp"
#include "core/remote_id.hpp"
#include "decode/emission_adapter.hpp"
#include "decode/registry.hpp"
#include "io/sigmf.hpp"
#include "receivers/adsb/modeac.hpp"
#include "receivers/broadcast_fm/detector.hpp"
#include "receivers/broadcast_fm/rds.hpp"
#include "receivers/gnss/acquire.hpp"
#include "receivers/ism_devices/detector.hpp"
#include "receivers/jammer/barrage.hpp"
#include "receivers/ism_devices/ook.hpp"
#include "receivers/ism_devices/tpms.hpp"
#include "receivers/marine/ais.hpp"
#include "receivers/marine/detector.hpp"
#include "receivers/navaids/detector.hpp"
#include "receivers/radalt/detector.hpp"
#include "receivers/rc_link/bredr.hpp"
#include "receivers/rc_link/detector.hpp"
#include "receivers/voice_fm/afsk.hpp"
#include "receivers/voice_fm/detector.hpp"

using namespace si;

namespace {

Json ref(const std::string& item) {
  std::ifstream jf(std::string(SI_TESTVEC) + "/signal_gaps/" + item + ".json");
  REQUIRE(jf);
  return Json::parse(jf);
}

template <class T>
std::vector<T> load(const std::string& name) {
  std::ifstream f(std::string(SI_TESTVEC) + "/signal_gaps/" + name, std::ios::binary | std::ios::ate);
  REQUIRE(f);
  const size_t n = static_cast<size_t>(f.tellg()) / sizeof(T);
  std::vector<T> x(n);
  f.seekg(0);
  f.read(reinterpret_cast<char*>(x.data()), static_cast<std::streamsize>(n * sizeof(T)));
  return x;
}

// the pipeline's signals of these receivers on a capture == python signals_of()
void check_scene(const Json& sc, const std::string& file, const std::vector<std::string>& receivers) {
  const Capture cap = Capture::from_iq(load<cf>(file), sc["fs"].get<double>(), sc["fc"].get<double>());
  const Analysis a = analyze(cap, plugins(receivers));
  Json got = Json::array();
  for (const auto& s : a.signals) {
    const std::string dec = s.features.contains("decoder") && s.features["decoder"].is_string()
                                ? s.features["decoder"].get<std::string>() : "";
    const bool mine = std::find(receivers.begin(), receivers.end(), dec) != receivers.end();
    if (!mine && s.decoded.status == "not_possible") continue;
    Json j = Json::object();
    j["label"] = s.label;
    j["status"] = s.decoded.status;
    j["confidence"] = s.confidence;
    j["fields"] = s.decoded.fields;
    j["f0"] = s.f0;
    j["f1"] = s.f1;
    got.push_back(j);
  }
  REQUIRE(got.size() == sc["signals"].size());
  for (size_t i = 0; i < got.size(); ++i) {
    CAPTURE(i);
    CHECK(got[i]["label"] == sc["signals"][i]["label"]);
    CHECK(got[i]["status"] == sc["signals"][i]["status"]);
    CHECK(got[i]["fields"].dump() == sc["signals"][i]["fields"].dump());
  }
}

}  // namespace

TEST_CASE("S11 Mode A/C: Gillham + squawk table == python") {
  const Json r = ref("modeac");
  const char* keys[11] = {"D2", "D4", "A1", "A2", "A4", "B1", "B2", "B4", "C1", "C2", "C4"};
  for (const auto& c : r["gillham"]) {
    rx::adsb::modeac::Bits p;
    for (const char* k : {"C1", "A1", "C2", "A2", "C4", "A4", "X", "B1", "D1", "B2", "D2", "B4", "D4"}) p[k] = 0;
    const int v = c[0].get<int>();
    for (int i = 0; i < 11; ++i) p[keys[i]] = (v >> (10 - i)) & 1;
    const auto alt = rx::adsb::modeac::gillham_ft(p);
    CAPTURE(v);
    CHECK((alt ? Json(*alt) : Json(nullptr)) == c[1]);
    CHECK(rx::adsb::modeac::squawk(p) == c[2].get<std::string>());
  }
}

TEST_CASE("S11 Mode A/C: decode / summarise on the envelope == python") {
  const Json r = ref("modeac");
  const fvec env = load<float>("modeac_env.f32");
  const auto dec = rx::adsb::modeac::decode(env, r["env_fs"].get<double>(), {});
  CHECK(Json(dec).dump() == r["decode"].dump());
  CHECK(rx::adsb::modeac::summarise(dec).dump() == r["summary"].dump());
  const auto busy = rx::adsb::modeac::decode(env, r["env_fs"].get<double>(), {{6e-3, 6e-3 + 120e-6}});
  CHECK(Json(busy).dump() == r["decode_busy"].dump());
}

TEST_CASE("S11 Mode A/C: adsb plugin on a capture == python") {
  check_scene(ref("modeac")["scene"], "modeac_scene.cf32", {"adsb"});
}

TEST_CASE("S13 AM voice: air_channel + voice_fm detect == python") {
  const Json r = ref("am");
  for (const auto& c : r["air_channel"]) {
    CAPTURE(c[0].get<double>());
    CHECK(rx::voice_fm::air_channel(c[0].get<double>(), c[1].get<double>()).dump() == c[2].dump());
  }
  const auto em = rx::voice_fm::detect(load<cf>("am_scene.cf32"), r["fs"].get<double>(), r["fc"].get<double>());
  REQUIRE(em.size() == r["emissions"].size());
  for (size_t i = 0; i < em.size(); ++i) {
    CAPTURE(i);
    CHECK(em[i].label == r["emissions"][i]["label"].get<std::string>());
    CHECK(em[i].features.dump() == r["emissions"][i]["features"].dump());
    CHECK(em[i].confidence == r["emissions"][i]["confidence"].get<double>());
  }
}

TEST_CASE("S13 AM voice: voice_fm plugin on a capture == python") {
  check_scene(ref("am")["scene"], "am_scene.cf32", {"voice_fm"});
}

TEST_CASE("S9 / S21 navaids: channel table + detect == python") {
  const Json r = ref("navaids");
  const auto ch = rx::navaids::channels();
  REQUIRE(ch.size() == r["channels"].size());
  for (size_t i = 0; i < ch.size(); ++i) {
    CHECK(ch[i].kind == r["channels"][i][0].get<std::string>());
    CHECK(ch[i].mhz == r["channels"][i][1].get<double>());
    CHECK((ch[i].paired ? Json(*ch[i].paired) : Json(nullptr)) == r["channels"][i][2]);
  }
  for (const char* name : {"nav_vhf", "nav_uhf"}) {
    const Json& c = r[name];
    const cfvec x = load<cf>(std::string(name) + ".cf32");
    const auto em = rx::navaids::detect(Channelizer(x, c["fs"].get<double>(), c["fc"].get<double>()),
                                        c["ppm"].get<double>());
    REQUIRE(em.size() == c["emissions"].size());
    for (size_t i = 0; i < em.size(); ++i) {
      CAPTURE(name);
      CAPTURE(i);
      CHECK(em[i].label == c["emissions"][i]["label"].get<std::string>());
      CHECK(em[i].features.dump() == c["emissions"][i]["features"].dump());
      CHECK(em[i].confidence == c["emissions"][i]["confidence"].get<double>());
    }
  }
}

TEST_CASE("S9 / S21 navaids: plugin on captures == python") {
  const Json r = ref("navaids");
  check_scene(r["nav_vhf"], "nav_vhf.cf32", {"navaids"});
  check_scene(r["nav_uhf"], "nav_uhf.cf32", {"navaids"});
}

TEST_CASE("S3 GNSS: C/A codes + chi2 thresholds == python") {
  const Json r = ref("gnss");
  for (int p = 1; p <= 32; ++p) {
    std::string s;
    for (uint8_t c : rx::gnss::ca_code(p)) s += static_cast<char>('0' + c);
    CHECK(s == r["codes"][std::to_string(p)].get<std::string>());
  }
  for (const auto& c : r["threshold"])
    CHECK(rx::gnss::threshold(c[0].get<int>(), c[1].get<double>()) == c[2].get<double>());
}

TEST_CASE("S3 GNSS: analyse on synthetic skies == python") {
  const Json r = ref("gnss");
  const double fs = r["fs"].get<double>(), fc = r["fc"].get<double>();
  for (const char* name : {"gnss_sky", "gnss_sky_noppm", "gnss_spoof", "gnss_cw"}) {
    const std::string file = std::string(name) == "gnss_sky_noppm" ? "gnss_sky" : name;
    const Json& c = r[name];
    std::optional<double> ppm;
    if (!c["ppm"].is_null()) ppm = c["ppm"].get<double>();
    const auto got = rx::gnss::analyse(Channelizer(load<cf>(file + ".cf32"), fs, fc), ppm);
    REQUIRE(got);
    CAPTURE(name);
    CHECK(got->dump() == c["result"].dump());
  }
}

TEST_CASE("S3 GNSS: plugin on a capture == python") { check_scene(ref("gnss")["scene"], "gnss_scene.cf32", {"gnss"}); }

TEST_CASE("S7 BLE trackers: parse_ble_ad (with device recognition) == python") {
  const Json r = ref("ble_trackers");
  for (const auto& c : r["ad"]) {
    const std::string h = c[0].get<std::string>();
    std::vector<uint8_t> b;
    for (size_t i = 0; i + 1 < h.size(); i += 2) b.push_back(static_cast<uint8_t>(std::stoul(h.substr(i, 2), nullptr, 16)));
    CAPTURE(h);
    CHECK(remote_id::parse_ble_ad(b).dump() == c[1].dump());
  }
}

TEST_CASE("S7 BLE trackers: AirTag advert through rc_link == python") {
  const Json r = ref("ble_trackers");
  std::vector<Emission> em;
  for (auto& e : rx::rc_link::detect(load<cf>("ble_tracker_scene.cf32"), r["fs"].get<double>(), r["fc"].get<double>()))
    if (e.decoded.is_object() && e.decoded.contains("crc_ok") && e.decoded["crc_ok"] == true) em.push_back(e);
  REQUIRE(em.size() == r["decoded"].size());
  for (size_t i = 0; i < em.size(); ++i) CHECK(em[i].decoded.dump() == r["decoded"][i].dump());
  const PluginList pl = plugins({"rc_link"});
  const EmissionReceiver* rc = nullptr;
  for (const auto& p : pl)
    if (p->name() == "rc_link") rc = dynamic_cast<const EmissionReceiver*>(p.get());
  REQUIRE(rc);
  CHECK(rc->pretty("ble", em) == r["label"].get<std::string>());
  const auto out = rc->summarise(em, nullptr);
  CHECK(out.fields["ble_devices"].dump() == r["ble_devices"].dump());
  CHECK(out.evidence.value_or("") == r["evidence"].get<std::string>());
}

TEST_CASE("S20 FM + RDS: block sync, groups and fields == python") {
  const Json r = ref("broadcast_fm");
  std::vector<int> bits;
  for (char c : r["bits"].get<std::string>()) bits.push_back(c - '0');
  const auto blks = rx::broadcast_fm::rds::blocks(bits);
  Json b = Json::array();
  for (const auto& x : blks) b.push_back(Json::array({x.pos, x.offset, x.data}));
  CHECK(b.dump() == r["blocks"].dump());
  CHECK(rx::broadcast_fm::rds::decode(rx::broadcast_fm::rds::groups(blks)).dump() == r["decode"].dump());
}

TEST_CASE("S20 FM + RDS: detect and plugin on a capture == python") {
  const Json r = ref("broadcast_fm");
  const auto em = rx::broadcast_fm::detect(Channelizer(load<cf>("fm_scene.cf32"), r["fs"].get<double>(),
                                                       r["fc"].get<double>()));
  REQUIRE(em.size() == r["emissions"].size());
  for (size_t i = 0; i < em.size(); ++i) {
    CAPTURE(i);
    CHECK(em[i].features.dump() == r["emissions"][i]["features"].dump());
    CHECK(em[i].confidence == r["emissions"][i]["confidence"].get<double>());
  }
  check_scene(r["scene"], "fm_scene.cf32", {"broadcast_fm"});
}

TEST_CASE("S15 AIS: NRZI / HDLC frames / CRC / fields == python") {
  const Json r = ref("ais");
  for (auto& [k, c] : r["codec"].items()) {
    std::vector<int> lv;
    for (char ch : c["levels"].get<std::string>()) lv.push_back(ch - '0');
    std::vector<int> bits = {0, 1, 0, 1, 0};
    for (int b : rx::marine::ais::unnrzi(lv)) bits.push_back(b);
    const auto fr = rx::marine::ais::frames(bits);
    CAPTURE(k);
    REQUIRE(fr.size() == c["frames"].size());
    for (size_t i = 0; i < fr.size(); ++i) {
      CHECK(fr[i].first == c["frames"][i][0].get<int64_t>());
      CHECK(rx::marine::ais::fields(fr[i].second).dump() == c["fields"][i].dump());
    }
  }
}

TEST_CASE("S15 AIS: detect and plugin on a capture == python") {
  const Json r = ref("ais");
  const auto em = rx::marine::detect(Channelizer(load<cf>("ais_scene.cf32"), r["fs"].get<double>(), r["fc"].get<double>()));
  REQUIRE(em.size() == r["emissions"].size());
  for (size_t i = 0; i < em.size(); ++i) {
    CHECK(em[i].decoded.dump() == r["emissions"][i]["decoded"].dump());
    CHECK(em[i].start_us == r["emissions"][i]["start_us"].get<double>());
  }
  check_scene(r["scene"], "ais_scene.cf32", {"marine"});
}

TEST_CASE("S19 APRS: AX.25 + APRS parsing and HDLC deframing == python") {
  const Json r = ref("aprs");
  for (const auto& c : r["frames"]) {
    const std::string h = c[0].get<std::string>();
    std::vector<uint8_t> o;
    for (size_t i = 0; i + 1 < h.size(); i += 2) o.push_back(static_cast<uint8_t>(std::stoul(h.substr(i, 2), nullptr, 16)));
    const auto f = hdlc::ax25(o);
    CAPTURE(h);
    REQUIRE(static_cast<bool>(f) == !c[1].is_null());
    if (!f) continue;
    CHECK(f->dest == c[1]["dest"].get<std::string>());
    CHECK(f->src == c[1]["src"].get<std::string>());
    CHECK(Json(f->path) == c[1]["path"]);
    CHECK(hdlc::latin1_utf8(f->info) == c[1]["info"].get<std::string>());
    CHECK(hdlc::aprs(f->dest, f->info).dump() == c[2].dump());
  }
  std::vector<int> bits;
  for (char ch : r["stream"].get<std::string>()) bits.push_back(ch - '0');
  Json got = Json::array();
  for (const auto& [p, o] : hdlc::deframe(bits)) {
    std::string h;
    for (uint8_t v : o) h += fmt("%02x", v);
    got.push_back(Json::array({p, h}));
  }
  CHECK(got.dump() == r["deframe"].dump());
}

TEST_CASE("S19 APRS: AFSK demodulator and voice_fm plugin == python") {
  const Json r = ref("aprs");
  CHECK(Json(rx::voice_fm::afsk::decode(load<float>("aprs_disc.f32"), 48e3)).dump() == r["afsk"].dump());
  check_scene(r["scene"], "aprs_scene.cf32", {"voice_fm"});
}

TEST_CASE("S12 GNSS swept jammer (PPD) == python") {
  const Json r = ref("gnss")["ppd"];
  for (size_t k = 0; k < r["results"].size(); ++k) {
    const auto got = rx::gnss::swept_ppd(
        Channelizer(load<cf>("gnss_ppd" + std::to_string(k) + ".cf32"), r["fs"].get<double>(), r["fc"].get<double>()),
        1575.42e6);
    CAPTURE(k);
    CHECK((got ? *got : Json(nullptr)).dump() == r["results"][k].dump());
  }
}

TEST_CASE("S8 radio altimeter: detect == python") {
  const Json r = ref("radalt");
  for (size_t i = 0; i < r["cases"].size(); ++i) {
    CAPTURE(i);
    const Json& c = r["cases"][i];
    const auto em = rx::radalt::detect(load<cf>("radalt_" + std::to_string(i) + ".cf32"), c["fs"].get<double>(),
                                       c["fc"].get<double>());
    REQUIRE(em.size() == c["emissions"].size());
    for (size_t k = 0; k < em.size(); ++k) {
      const Json& w = c["emissions"][k];
      CHECK(em[k].label == w["label"].get<std::string>());
      CHECK(em[k].confidence == w["confidence"].get<double>());
      CHECK(em[k].features.dump() == w["features"].dump());
      CHECK(em[k].start_us == w["start_us"].get<double>());
      CHECK(em[k].f_lo_mhz == w["f_lo_mhz"].get<double>());
      CHECK(em[k].f_hi_mhz == w["f_hi_mhz"].get<double>());
    }
  }
}

TEST_CASE("S8 radio altimeter: radalt plugin on a capture == python") {
  check_scene(ref("radalt")["scene"], "radalt_0.cf32", {"radalt"});
}

TEST_CASE("S6 315 / 433 MHz: TPMS symbols, OOK pulses == python") {
  const Json r = ref("ism_devices");
  for (const auto& c : r["symbols"]) {
    std::vector<int> st;
    for (char ch : c["s"].get<std::string>()) st.push_back(ch - '0');
    const auto d = rx::ism_devices::tpms::decode(st);
    CHECK((d ? *d : Json(nullptr)).dump() == c["decode"].dump());
  }
  for (const auto& c : r["pulses"]) {
    std::vector<rx::ism_devices::ook::Pulse> p;
    for (const auto& q : c["p"]) {
      const double sp = q[2].get<double>();
      p.push_back({q[0].get<double>(), q[1].get<double>(), sp < 0 ? std::numeric_limits<double>::infinity() : sp});
    }
    const auto fc = rx::ism_devices::ook::fixed_code(p);
    const auto kq = rx::ism_devices::ook::keeloq(p);
    CHECK((fc ? *fc : Json(nullptr)).dump() == c["fixed_code"].dump());
    CHECK((kq ? *kq : Json(nullptr)).dump() == c["keeloq"].dump());
  }
}

TEST_CASE("S6 315 / 433 MHz: detect on a scene == python") {
  const Json r = ref("ism_devices");
  const auto em = rx::ism_devices::detect(
      Channelizer(load<cf>("ism_scene.cf32"), r["fs"].get<double>(), r["fc"].get<double>()));
  REQUIRE(em.size() == r["emissions"].size());
  for (size_t i = 0; i < em.size(); ++i) {
    CAPTURE(i);
    const Json& w = r["emissions"][i];
    CHECK(em[i].label == w["label"].get<std::string>());
    CHECK(em[i].confidence == w["confidence"].get<double>());
    CHECK(em[i].start_us == w["start_us"].get<double>());
    CHECK(em[i].end_us == w["end_us"].get<double>());
    CHECK(em[i].center_mhz == w["center_mhz"].get<double>());
    CHECK(em[i].features.dump() == w["features"].dump());
    CHECK(em[i].decoded.dump() == w["decoded"].dump());
  }
}

TEST_CASE("S6 315 / 433 MHz: ism_devices + rc_link plugins on the scene == python") {
  check_scene(ref("ism_devices")["scene"], "ism_scene.cf32", {"ism_devices", "rc_link"});
}

TEST_CASE("S12 noise jammers: barrage analyse / covers_band == python") {
  const Json r = ref("barrage");
  for (size_t i = 0; i < r["cases"].size(); ++i) {
    CAPTURE(i);
    const Json& c = r["cases"][i];
    const auto a = rx::jammer::barrage::analyse(
        Channelizer(load<cf>("barrage_" + std::to_string(i) + ".cf32"), r["fs"].get<double>(), r["fc"].get<double>()),
        c["f0"].get<double>(), c["f1"].get<double>());
    CHECK((a ? *a : Json(nullptr)).dump() == c["analyse"].dump());
  }
  for (const auto& c : r["covers"])
    CHECK(Json(rx::jammer::barrage::covers_band(c[0].get<double>(), c[1].get<double>())) == c[2]);
}

TEST_CASE("S12 noise jammers: jammer plugins on a capture == python") {
  check_scene(ref("barrage")["scene"], "barrage_0.cf32", {"jammer"});
}

TEST_CASE("S10 pulse deinterleaving: pulse lists == python") {
  const Json r = ref("deinterleave");
  for (size_t i = 0; i < r["cases"].size(); ++i) {
    CAPTURE(i);
    const Json& c = r["cases"][i];
    const auto out = deint::deinterleave(c["toa"].get<std::vector<double>>(), c["pw"].get<std::vector<double>>());
    CHECK(out.dump() == c["out"].dump());
  }
}

TEST_CASE("S10 pulse deinterleaving: pulse features and radar evidence on a scene == python") {
  const Json r = ref("deinterleave")["scene"];
  const Capture cap = Capture::from_iq(load<cf>("radar_scene.cf32"), r["fs"].get<double>(), r["fc"].get<double>());
  const PluginList none;
  const Analysis a = analyze(cap, none);
  Json got = Json::array();
  for (const auto& s : a.signals) {
    if (!s.features.contains("pulse_n") || s.features["pulse_n"].is_null() || s.features["pulse_n"].get<int64_t>() == 0)
      continue;
    Json f = Json::object();
    for (const char* k : {"pulse_n", "pri_s", "pulse_emitters", "pulse_unassigned"})
      f[k] = s.features.contains(k) ? s.features[k] : Json(nullptr);
    got.push_back({{"label", s.label}, {"confidence", s.confidence}, {"evidence", s.evidence}, {"features", f}});
  }
  REQUIRE(got.size() == r["signals"].size());
  for (size_t i = 0; i < got.size(); ++i) {
    CAPTURE(i);
    CHECK(got[i]["label"] == r["signals"][i]["label"]);
    CHECK(got[i]["evidence"] == r["signals"][i]["evidence"]);
    CHECK(got[i]["features"].dump() == r["signals"][i]["features"].dump());
  }
}

TEST_CASE("S7 Bluetooth BR/EDR: sync words, access-code search, decode == python") {
  const Json r = ref("bredr");
  for (const auto& c : r["syncwords"])
    CHECK(fmt("%016llx", static_cast<unsigned long long>(rx::rc_link::bredr::syncword(c[0].get<uint32_t>()))) ==
          c[1].get<std::string>());
  for (const auto& c : r["finds"]) {
    std::vector<int64_t> b;
    for (char ch : c["bits"].get<std::string>()) b.push_back(ch - '0');
    const auto f = rx::rc_link::bredr::find(b);
    const Json got = f ? Json::array({std::get<0>(*f), std::get<1>(*f), std::get<2>(*f)}) : Json(nullptr);
    CHECK(got.dump() == c["find"].dump());
  }
  CHECK(rx::rc_link::bredr::try_decode(load<cf>("bredr_seg.cf32"), r["seg_ofs"].get<double>()).dump() ==
        r["try_decode"].dump());
}

TEST_CASE("S7 Bluetooth BR/EDR: rc_link on a capture == python") {
  check_scene(ref("bredr")["scene"], "bredr_scene.cf32", {"rc_link"});
}
