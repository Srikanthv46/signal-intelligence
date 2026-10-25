// receivers/lte against the Python receiver on a synthetic FDD downlink (tests/vectors/lte/gen_lte.py).
#include <doctest/doctest.h>

#include <cmath>
#include <fstream>

#include "core/dsp.hpp"
#include "core/spectrogram.hpp"
#include "receivers/lte/cellsearch.hpp"

using namespace si;

namespace {

const Json& vec() {
  static const Json j = [] {
    std::ifstream f(std::string(SI_TESTVEC) + "/lte/lte.json");
    REQUIRE_MESSAGE(f.good(), "missing tests/vectors/lte/lte.json (run gen_lte.py)");
    return Json::parse(f);
  }();
  return j;
}

cfvec load_iq() {
  std::ifstream f(std::string(SI_TESTVEC) + "/lte/lte_iq.cf32", std::ios::binary | std::ios::ate);
  REQUIRE(f.good());
  const auto n = static_cast<size_t>(f.tellg()) / sizeof(cf);
  cfvec v(n);
  f.seekg(0);
  f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(n * sizeof(cf)));
  return v;
}

bool close(double a, double b, double rel = 1e-9) { return std::fabs(a - b) <= rel * std::max(1.0, std::fabs(b)); }

}  // namespace

TEST_CASE("lte: PSS template spectra equal Python") {
  const auto bank = rx::lte::template_bank(512);
  const auto& t = vec()["tpl_512"];
  double worst = 0;
  for (size_t i = 0; i < t.size(); ++i)
    worst = std::max(worst, std::abs(bank->T.data[i] - cd(t[i][0].get<double>(), t[i][1].get<double>())));
  CHECK(worst < 1e-12);
}

TEST_CASE("lte: search_cells / wide_ppm_search on a synthetic FDD cell equal Python") {
  const cfvec iq = load_iq();
  const Channelizer ch(iq, vec()["fs"].get<double>(), vec()["fc"].get<double>());
  Signal sig;
  sig.f0 = vec()["sig"][0];
  sig.f1 = vec()["sig"][1];
  const auto cells = rx::lte::search_cells(ch, sig, 0.0);
  const auto& ref = vec()["cells"];
  REQUIRE(cells.size() == ref.size());
  REQUIRE(!cells.empty());
  for (size_t i = 0; i < cells.size(); ++i) {
    const auto& c = cells[i];
    const auto& r = ref[i];
    CHECK(c.raster_hz == r["raster_hz"].get<double>());
    CHECK(c.carrier_hz == r["carrier_hz"].get<double>());
    CHECK(close(c.pss_z, r["pss_z"]));
    CHECK(close(c.pss_isolation, r["pss_isolation"]));
    CHECK(c.nid2 == r["nid2"].get<int>());
    CHECK(c.phase == r["phase"].get<int64_t>());
    CHECK(c.sss_ok == r["sss_ok"].get<bool>());
    REQUIRE(c.inst.size() == r["inst"].size());
    for (size_t k = 0; k < c.inst.size(); ++k) {
      CHECK(c.inst[k].first == r["inst"][k][0].get<int64_t>());
      CHECK(c.inst[k].second == r["inst"][k][1].get<int64_t>());
    }
    REQUIRE(c.sss.has_value() == !r["sss"].is_null());
    if (c.sss) {
      CHECK(c.sss->duplex == r["sss"]["duplex"].get<std::string>());
      CHECK(c.sss->cp == r["sss"]["cp"].get<std::string>());
      CHECK(c.sss->nid1 == r["sss"]["nid1"].get<int>());
      CHECK(c.sss->parity == r["sss"]["parity"].get<int>());
      CHECK(c.sss->n == r["sss"]["n"].get<int>());
      CHECK(close(c.sss->z, r["sss"]["z"]));
      CHECK(close(c.sss->z2, r["sss"]["z2"]));
    }
  }
  const auto w = rx::lte::wide_ppm_search(ch, {945.0e6}, 20e3, 5e3);
  REQUIRE(w.has_value() == !vec()["wide"].is_null());
  if (w) {
    for (const char* k : {"ppm", "z", "carrier_hz", "raster_hz"})
      CHECK(close((*w)[k].get<double>(), vec()["wide"][k].get<double>()));
  }
}

TEST_CASE("lte: bandwidth_from_energy equals Python") {
  DetectResult det;
  const auto P = vec()["bw_P"].get<std::vector<float>>();
  det.coarse.f = vec()["bw_f"].get<dvec>();
  det.coarse.P = fmat(P.size() / det.coarse.f.size(), det.coarse.f.size());
  det.coarse.P.data = P;
  det.floor_coarse_db = dvec(det.coarse.f.size(), 8.0);
  for (const auto& x : vec()["bw"]) {
    const auto bw = rx::lte::bandwidth_from_energy(det, x["fcar"].get<double>());
    REQUIRE(bw.has_value() == !x["bw"].is_null());
    if (bw) CHECK(*bw == x["bw"].get<double>());
  }
}
