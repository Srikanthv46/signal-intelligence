// receivers/adsb against python receivers.adsb (tests/vectors/adsb/gen_adsb.py; pyModeS only as the
// oracle behind the Python package): demod on the real_1090_commb envelope, per-frame fields of
// every DF / typecode / Comm-B register, per-aircraft tracks (CPR global/local, BDS50/60, met, flags).
#include <doctest/doctest.h>

#include <fstream>
#include <string>

#include "core/pyfmt.hpp"
#include "receivers/adsb/codec.hpp"
#include "receivers/adsb/messages.hpp"
#include "receivers/adsb/modes.hpp"
#include "receivers/adsb/tracks.hpp"

using namespace si;
using namespace si::rx::adsb;

namespace {

Json load_json(const std::string& name) {
  std::ifstream f(std::string(SI_TESTVEC) + "/adsb/" + name);
  REQUIRE(f.good());
  return Json::parse(f);
}

fvec load_f32(const std::string& name) {
  std::ifstream f(std::string(SI_TESTVEC) + "/adsb/" + name, std::ios::binary | std::ios::ate);
  REQUIRE(f.good());
  const auto n = static_cast<size_t>(f.tellg()) / sizeof(float);
  fvec v(n);
  f.seekg(0);
  f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(n * sizeof(float)));
  return v;
}

// identical serialisation (types, key order, float repr)
bool same(const Json& a, const Json& b, const std::string& what, int& reported) {
  if (a.dump() == b.dump()) return true;
  if (reported++ < 8) {
    INFO(what);
    const std::string x = a.dump(), y = b.dump();
    size_t i = 0;
    while (i < x.size() && i < y.size() && x[i] == y[i]) ++i;
    const size_t s = i > 150 ? i - 150 : 0;
    INFO("c++:    ..." << x.substr(s, 400));
    INFO("python: ..." << y.substr(s, 400));
    CHECK(false);
  }
  return false;
}

Json frames_json(const std::vector<Json>& v) {
  Json a = Json::array();
  for (const auto& x : v) a.push_back(x);
  return a;
}

}  // namespace

TEST_CASE("adsb: CRC-24 and Annex 10 codes") {
  // DF17 from the real capture: residual 0; flipping one bit breaks it
  std::string b = hex_bits("8D801648583B4760C66D9D704ECB");
  std::vector<uint8_t> bits(b.size());
  for (size_t i = 0; i < b.size(); ++i) bits[i] = b[i] == '1';
  CHECK(parity_residual(bits.data(), bits.size()) == 0);
  bits[40] ^= 1;
  CHECK(parity_residual(bits.data(), bits.size()) != 0);
  CHECK_THROWS(squawk_id13("0000000000000"));   // all zeros / all ones are rejected like the reference
  CHECK_THROWS(squawk_id13("1111111111111"));
  CHECK(squawk_id13("1111111011111") == "7677");
  CHECK(!altitude_ac13("0000000000000"));
  CHECK(*altitude_ac13("0100000010000") == 11800);   // Q = 1, N = 512
  CHECK(cpr_nl(0.0) == 59);
  CHECK(cpr_nl(87.0) == 2);
  CHECK(cpr_nl(-88.0) == 1);
  CHECK(cpr_nl(10.0) == 59);
  CHECK(cpr_nl(10.48) == 58);
}

TEST_CASE("adsb: demod on the real_1090_commb envelope") {
  const Json ref = load_json("demod.json");
  const fvec env = load_f32("env_full.f32");
  REQUIRE(env.size() == ref["n_env"].get<size_t>());
  const float thr = ref["thr"].get<float>();
  int rep = 0;
  same(frames_json(demod(env, thr)), ref["full"], "full envelope", rep);
  CHECK(ref["full"].size() == 6);
  for (const auto& w : ref["windows"]) {
    const size_t s0 = w["start"].get<size_t>(), n = w["n"].get<size_t>();
    const fvec part(env.begin() + static_cast<std::ptrdiff_t>(s0), env.begin() + static_cast<std::ptrdiff_t>(s0 + n));
    same(frames_json(demod(part, thr)), w["msgs"], "window " + std::to_string(s0), rep);
  }
}

TEST_CASE("adsb: fields of every DF / typecode / Comm-B register") {
  const Json ref = load_json("fields.json");
  std::vector<Json> in;
  for (const auto& f : ref) {   // the demod() keys come first in each frame dict
    Json m = Json::object();
    for (const char* k : {"t_s", "df", "bits", "hex", "parity_residual", "snr_db", "aa", "confirmed"}) m[k] = f[k];
    in.push_back(m);
  }
  fields(in);
  int rep = 0, bad = 0;
  for (size_t i = 0; i < in.size(); ++i) bad += !same(in[i], ref[i], "frame " + std::to_string(i), rep);
  CHECK(bad == 0);
  MESSAGE(ref.size() << " frames, " << bad << " differ");
}

TEST_CASE("adsb: tracks (CPR, kinematics, Comm-B resolution, met, interrogation, flags)") {
  const Json cases = load_json("tracks.json");
  int rep = 0, bad = 0, n_err = 0, n_pos = 0;
  for (size_t c = 0; c < cases.size(); ++c) {
    const Json& k = cases[c];
    std::vector<Json> frames;
    for (const auto& f : k["frames"]) {
      Json m = f;
      if (m.contains("commb") && m["commb"].contains("register")) m["commb"].erase("register");
      frames.push_back(m);
    }
    const std::string what = "case " + std::to_string(c);
    if (k.contains("error")) {
      ++n_err;
      std::string err;
      try {
        tracks(frames, k["ref"]);
      } catch (const std::exception& e) {
        err = e.what();
      }
      if (err != k["error"].get<std::string>()) {
        ++bad;
        if (rep++ < 8) {
          INFO(what << " c++ error '" << err << "' python '" << k["error"].get<std::string>() << "'");
          CHECK(false);
        }
      }
      continue;
    }
    const Json ac = tracks(frames, k["ref"]);
    for (const auto& a : ac) n_pos += a.contains("position");
    bool ok = same(ac, k["aircraft"], what + " aircraft", rep);
    ok = same(frames_json(frames), k["frames"], what + " frames", rep) && ok;
    bad += !ok;
  }
  CHECK(bad == 0);
  MESSAGE(cases.size() << " cases (" << n_err << " raising in python), " << bad << " differ, " << n_pos
                       << " aircraft with a position");
}

TEST_CASE("adsb: CPR local / global decoding and NL (DO-260B A.1.7)") {
  const Json ref = load_json("cpr.json");
  int rep = 0, bad = 0;
  for (const auto& c : ref["nl"]) {
    const double lat = c[0].get<double>();
    if (cpr_nl(lat) != c[1].get<int>()) {
      ++bad;
      if (rep++ < 8) {
        INFO("NL(" << py_float_repr(lat) << ") c++ " << cpr_nl(lat) << " python " << c[1].get<int>());
        CHECK(false);
      }
    }
  }
  for (const auto& c : ref["local"]) {
    const auto p = cpr_local(c["hex"].get<std::string>(), c["lat_ref"].get<double>(), c["lon_ref"].get<double>(),
                             c["surface"].get<bool>());
    bad += !same(Json::array({p.first, p.second}), c["pos"], "local " + c.dump(), rep);
  }
  for (const auto& c : ref["global"]) {
    const auto p = cpr_airborne_global(c["a"].get<std::string>(), c["b"].get<std::string>(), c["ta"].get<double>(),
                                       c["tb"].get<double>());
    bad += !same(p ? Json::array({p->first, p->second}) : Json(), c["pos"], "global " + c.dump(), rep);
  }
  CHECK(bad == 0);
  MESSAGE(ref["nl"].size() << " NL, " << ref["local"].size() << " local, " << ref["global"].size() << " global; "
                           << bad << " differ");
}

TEST_CASE("adsb: BDS50 / BDS60 resolution against the aircraft's ADS-B state") {
  const Json ref = load_json("is50or60.json");
  int rep = 0, bad = 0;
  for (const auto& c : ref) {
    const Json r = commb_is50or60(hex_bits(c["hex"].get<std::string>()), c["spd"].get<double>(), c["trk"].get<double>(),
                                  c["alt"].get<double>());
    bad += !same(r, c["reg"], c.dump(), rep);
  }
  CHECK(bad == 0);
  MESSAGE(ref.size() << " cases, " << bad << " differ");
}

TEST_CASE("adsb: emitter category B6 -> UAV flag (python test_uav_emitter_category)") {
  // identification ES of 801648 / IGO311 with ME byte 0 = TC 3, CA 6; parity from the Mode S CRC-24
  std::string h = "8D8016481E2473F3C71820";
  uint64_t hi = std::stoull(h.substr(0, 11), nullptr, 16), lo = std::stoull(h.substr(11), nullptr, 16);
  std::vector<uint8_t> bits;
  for (int i = 43; i >= 0; --i) bits.push_back((hi >> i) & 1);
  for (int i = 43; i >= 0; --i) bits.push_back((lo >> i) & 1);
  for (int i = 0; i < 24; ++i) bits.push_back(0);
  for (size_t i = 0; i < 88; ++i)
    if (bits[i]) {
      const uint32_t g = 0x1FFF409;
      for (int j = 0; j < 25; ++j) bits[i + j] ^= (g >> (24 - j)) & 1;
    }
  uint32_t par = 0;
  for (size_t i = 88; i < 112; ++i) par = (par << 1) | bits[i];
  h += fmt("%06X", par);
  Json m = Json::object();
  m["t_s"] = 0.0; m["df"] = 17; m["bits"] = 112; m["hex"] = h; m["parity_residual"] = "000000";
  m["snr_db"] = 20.0; m["aa"] = nullptr; m["confirmed"] = true;
  std::vector<Json> in{m};
  fields(in);
  CHECK(in[0]["category"] == "unmanned aerial vehicle");
  in[0]["t"] = 0.0;
  const Json ac = tracks(in, nullptr);
  REQUIRE(ac.size() == 1);
  CHECK(ac[0]["uav"] == true);
  bool flagged = false;
  for (const auto& f : ac[0]["flags"]) flagged |= f == "UAV (ADS-B emitter category B6)";
  CHECK(flagged);
}
