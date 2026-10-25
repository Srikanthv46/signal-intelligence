// receivers/nr and core/seq3gpp against the Python receiver and py3gpp (black-box oracle):
// tests/vectors/nr/gen_nr.py.
#include <doctest/doctest.h>

#include <cmath>
#include <fstream>

#include "core/seq3gpp.hpp"
#include "core/types.hpp"
#include "receivers/nr/bch.hpp"
#include "receivers/nr/ssb.hpp"

using namespace si;

namespace {

const Json& vec() {
  static const Json j = [] {
    std::ifstream f(std::string(SI_TESTVEC) + "/nr/nr.json");
    REQUIRE_MESSAGE(f.good(), "missing tests/vectors/nr/nr.json (run gen_nr.py)");
    return Json::parse(f);
  }();
  return j;
}

cdvec cvec(const Json& a) {
  cdvec v;
  for (const auto& p : a) v.emplace_back(p[0].get<double>(), p[1].get<double>());
  return v;
}

double rel_err(double a, double b) { return std::fabs(a - b) / std::max(1e-300, std::max(std::fabs(a), std::fabs(b))); }

}  // namespace

TEST_CASE("seq3gpp gold / qpsk: bit-exact against Python and py3gpp nrPRBS") {
  for (const auto& g : vec()["gold"]) {
    const auto c = seq3gpp::gold(g["c_init"].get<uint32_t>(), 3000);
    const std::string bits = g["bits"];
    size_t bad = 0;
    for (size_t i = 0; i < 3000; ++i) bad += c[i] != bits[i] - '0';
    CHECK(bad == 0);
  }
  const cdvec q = seq3gpp::qpsk(seq3gpp::gold(1234, 288));
  const cdvec r = cvec(vec()["qpsk_1234"]);
  REQUIRE(q.size() == r.size());
  for (size_t i = 0; i < q.size(); ++i) CHECK(q[i] == r[i]);
}

TEST_CASE("nr: PSS / SSS / PBCH DMRS sequences (TS 38.211) equal py3gpp") {
  for (int n = 0; n < 3; ++n) {
    const dvec d = rx::nr::pss_seq(n);
    for (int i = 0; i < 127; ++i) CHECK(d[i] == vec()["pss"][n][i].get<int>());
  }
  for (const auto& s : vec()["sss"]) {
    const dvec d = rx::nr::sss_seq(s["pci"].get<int>());
    for (int i = 0; i < 127; ++i) CHECK(d[i] == s["d"][i].get<int>());
  }
  for (const auto& s : vec()["dmrs"]) {
    const cdvec d = rx::nr::pbch_dmrs(s["pci"].get<int>(), s["ib"].get<int>());
    const cdvec r = cvec(s["d"]);
    REQUIRE(d.size() == r.size());
    for (size_t i = 0; i < d.size(); ++i) CHECK(std::abs(d[i] - r[i]) < 1e-15);
  }
}

TEST_CASE("nr: polar code construction (TS 38.212 5.3.1) and rate recovery (5.4.1)") {
  const auto& q = rx::nr::polar_reliability();
  REQUIRE(q.size() == 1024);
  for (size_t i = 0; i < 1024; ++i) CHECK(q[i] == vec()["reliability"][i].get<int>());
  CHECK(rx::nr::polar_info_positions(56, 512) == vec()["info_56_512"].get<std::vector<int>>());
  CHECK(rx::nr::polar_input_interleaver(56) == vec()["interleaver_56"].get<std::vector<int>>());
  const dvec rr = rx::nr::rate_recover_polar(vec()["rr_in"].get<dvec>(), 512);
  CHECK(rr == vec()["rr_out"].get<dvec>());
}

TEST_CASE("nr: BCH decode (SC polar, CRC-24C, descrambling) identical to the Python receiver") {
  int ok = 0;
  for (const auto& b : vec()["bch"]) {
    const dvec soft = b["soft"].get<dvec>();
    const int pci = b["pci"].get<int>();
    CHECK(rx::nr::polar_decode(rx::nr::rate_recover_polar(soft, 512), 56) == b["polar"].get<std::vector<int>>());
    const auto r = rx::nr::bch_decode(soft, pci);
    CHECK(r.crc_ok == b["crc_ok"].get<bool>());
    CHECK(r.scrblk == b["scrblk"].get<std::vector<int>>());
    CHECK(r.payload == b["payload"].get<std::vector<int>>());
    CHECK(r.sfn_lsb == b["sfn_lsb"].get<std::vector<int>>());
    CHECK(r.hrf == b["hrf"].get<int>());
    ok += r.crc_ok;
  }
  CHECK(ok > 0);
}

TEST_CASE("nr: SSB pieces (SSS scores, PBCH LLRs, peaks, GSCN, PSS templates)") {
  const cdvec g = cvec(vec()["grid"]);
  Mat<cd> G(4, 240);
  for (size_t i = 0; i < g.size(); ++i) G.data[i] = g[i];
  for (int n2 = 0; n2 < 3; ++n2) {
    const dvec s = rx::nr::sss_scores(G, n2);
    const dvec r = vec()["sss_scores"][std::to_string(n2)].get<dvec>();
    for (size_t i = 0; i < s.size(); ++i) CHECK(rel_err(s[i], r[i]) < 1e-12);
  }
  for (const auto& p : vec()["pbch_llr"]) {
    const auto [ib, llr] = rx::nr::pbch_llr(G, p["pci"].get<int>());
    CHECK(ib == p["ib"].get<int>());
    const dvec r = p["llr"].get<dvec>();
    REQUIRE(llr.size() == r.size());
    double worst = 0;
    for (size_t i = 0; i < llr.size(); ++i) worst = std::max(worst, std::fabs(llr[i] - r[i]));
    CHECK(worst < 1e-12);
  }
  const auto pk = rx::nr::peaks(vec()["peaks_in"].get<dvec>());
  CHECK(pk == vec()["peaks"].get<std::vector<int64_t>>());
  for (const auto& x : vec()["gscn"]) CHECK(rx::nr::gscn(x["lo"], x["hi"]) == x["f"].get<dvec>());
  const auto bank = rx::nr::template_bank(30e3, 256 * 30e3, 1024);
  const cdvec t = cvec(vec()["tpl_30k_1024"]);
  double worst = 0;
  for (size_t i = 0; i < t.size(); ++i) worst = std::max(worst, std::abs(bank->T.data[i] - t[i]));
  CHECK(worst < 1e-12);
}
