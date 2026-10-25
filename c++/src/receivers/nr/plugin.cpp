// Pipeline plugin: NR cells (SSB + MIB) as capture-level signals (python: receivers/nr/plugin.py).
// Runs after lte: reads the receiver ppm it publishes (ctx.cache["ppm"]), else the caller's
// ctx.cache["sweep_ppm"].
#include <algorithm>
#include <map>
#include <set>

#include "core/bands_india.hpp"
#include "core/npcompat.hpp"
#include "core/parallel.hpp"
#include "core/pyfmt.hpp"
#include "decode/plugin.hpp"
#include "receivers/nr/ssb.hpp"

namespace si::rx::nr {

namespace {

bool has(const std::vector<std::string>& v, const std::string& k) { return std::find(v.begin(), v.end(), k) != v.end(); }
bool truthy(const Json* j) { return j && j->is_object() && !j->empty(); }

class NrDecoder : public Decoder {
 public:
  std::string name() const override { return "nr_ssb"; }
  std::string package() const override { return "nr"; }
  std::string protocol() const override { return "3GPP 5G NR SS/PBCH block"; }
  bool capture_level() const override { return true; }
  std::vector<std::string> after() const override { return {"lte"}; }

  bool can_decode(const Signal& sig, DecodeContext&) const override {
    const auto keys = bands::band_keys(sig.f0, sig.f1);
    for (const auto& b : nr_bands())
      if (has(keys, b.key)) return true;
    return false;
  }

  std::vector<Signal> capture_signals(DecodeContext& ctx) const override {
    std::vector<const Signal*> cand;
    for (const Signal& s : ctx.det->signals)
      if (s.flags.empty() && s.snr_db >= 3 && s.bw() >= 3e6 && can_decode(s, ctx)) cand.push_back(&s);
    if (cand.empty()) return {};
    const Json* p = ctx.cache.get<Json>("ppm");
    const Json* sw = ctx.cache.get<Json>("sweep_ppm");
    const Json* cal = truthy(p) ? p : (truthy(sw) ? sw : nullptr);
    const double ppm = cal ? (*cal)["ppm"].get<double>() : 0.0;
    struct Todo {
      const Signal* sig;
      double f_ss;
      std::string bname;
      double scs;
    };
    std::set<int64_t> seen;
    std::vector<Todo> todo;
    const Capture& cap = ctx.cap;
    const double lo_cap = cap.fc - 0.49 * cap.fs, hi_cap = cap.fc + 0.49 * cap.fs;  // AD936x passband
    for (const Signal* sig : cand) {
      const auto keys = bands::band_keys(sig->f0, sig->f1);
      const NrBand* nb = nullptr;
      for (const auto& b : nr_bands())
        if (has(keys, b.key)) {
          nb = &b;
          break;
        }
      const auto bl = bands_for(nb->key);
      double mn = bl[0].first, mx = bl[0].second;
      for (const auto& x : bl) mn = std::min(mn, x.first), mx = std::max(mx, x.second);
      const double blo = std::max(lo_cap, mn) + 120 * nb->scs;
      const double bhi = std::min(hi_cap, mx) - 120 * nb->scs;
      for (double f_ss : gscn(blo, bhi)) {
        const int64_t k = py_round_int(f_ss / 10e3);
        if (seen.count(k)) continue;
        seen.insert(k);
        todo.push_back({sig, f_ss, nb->name, nb->scs});
      }
    }
    if (todo.empty()) return {};
    // one template bank per SSB SCS (every slice of a capture has the same length and rate)
    std::map<double, std::shared_ptr<const TemplateBank>> banks;
    for (const auto& t : todo)
      if (!banks.count(t.scs)) {
        const SliceDims d = slice_dims(ctx.ch, kNfft * t.scs);
        banks[t.scs] = template_bank(t.scs, d.ofs, d.nf);
      }
    const dvec coarse = np::arange<double>(-12e3, 12.01e3, 3e3);
    auto slice_of = [&](const Todo& t) { return SSBSlice(ctx.ch, t.f_ss * (1 + ppm * 1e-6), t.scs, banks.at(t.scs)); };
    // GSCN hypotheses are independent
    const std::vector<PssBest> scan = pmap(todo, [&](const Todo& t) { return slice_of(t).best(coarse); });
    dvec zs;
    for (const auto& b : scan) zs.push_back(b.z);
    // z of a loaded TDD band has a heavy tail: an SSB must stand well clear of the scan
    const double z_thr = std::max(40.0, 2.5 * np::median(zs));
    std::vector<size_t> order(scan.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return -scan[a].z < -scan[b].z; });
    std::vector<Signal> out;
    for (size_t i : order) {
      PssBest b = scan[i];
      const Todo& t = todo[i];
      if (b.z < z_thr) break;
      bool dup = false;
      for (const Signal& o : out)
        if (std::fabs(o.features["ssb_center_hz"].get<double>() - t.f_ss) < 240 * t.scs) dup = true;
      if (dup) continue;
      const SSBSlice sl = slice_of(t);
      const PssBest fine = sl.best(np::arange<double>(b.df - 2e3, b.df + 2.01e3, 500.0));
      if (fine.z >= b.z) b = fine;
      auto [pci, info] = decode_ssb_multi(sl, b);
      if (!pci) continue;
      if (!info.crc_ok && (info.sss_z < 4.5 || info.sss_z / std::max(info.sss_z2, 1e-3) < 1.6)) continue;
      b.nid2 = info.nid2;
      Signal o = make_signal(ctx, *t.sig, t.f_ss, t.scs, t.bname, b, *pci, info, cal);
      o.features["ssb_center_hz"] = py_round_int(t.f_ss);
      out.push_back(std::move(o));
    }
    return out;
  }

  Signal make_signal(DecodeContext& ctx, const Signal& sig, double f_ss, double scs, const std::string& bname,
                     const PssBest& b, int pci, const SsbInfo& info, const Json* cal) const {
    std::vector<std::string> steps;
    steps.push_back(cal ? fmt("SSB centre hypotheses on the GSCN raster inside %.2f-%.2f MHz (%g kHz SCS), receiver "
                              "error %+.2f ppm applied",
                              sig.f0 / 1e6, sig.f1 / 1e6, scs / 1e3, (*cal)["ppm"].get<double>())
                        : std::string("SSB centre hypotheses on the GSCN raster (ppm unknown)"));
    steps.push_back(fmt("PSS N_ID2 %d z %.1f at %.3f ms", b.nid2, b.z, static_cast<double>(b.pos) / (256 * scs) * 1e3));
    steps.push_back(fmt("SSS (differential, channel-free; summed over SSB instances) N_ID1 %d (z %.1f vs %.1f) -> PCI %d",
                        *info.nid1, info.sss_z, info.sss_z2, pci));
    Json idx = Json::array();
    for (int v : info.ssb_indices) idx.push_back(v);
    steps.push_back(fmt("%zu SSB instance(s) (beam sweep / bursts) combined: SSS scores summed; per-SSB PBCH DMRS -> "
                        "i_SSB ",
                        info.n_ssb_instances) +
                    py_repr(idx) + ", DMRS channel estimate, MRC");
    steps.push_back(
        "PBCH LLRs descrambled per i_SSB and summed per half-frame; polar decode + CRC-24C (py3gpp.nrBCHDecode)");
    Json others = Json::array();
    for (const auto& c : info.other_cells) others.push_back(Json{{"pci", c.pci}, {"sss_z", c.sss_z}});
    Json fields = Json::object();
    fields["pci"] = pci;
    fields["nid1"] = *info.nid1;
    fields["nid2"] = b.nid2;
    fields["co_channel_cells"] = others;
    fields["ssb_index"] = info.ibar_ssb ? Json(*info.ibar_ssb) : Json(nullptr);
    fields["ssb_center_hz"] = py_round_int(f_ss);
    fields["gscn_freq_hz"] = py_round_int(f_ss);
    fields["nr_arfcn"] = bands::nr_arfcn(f_ss);
    fields["ssb_scs_hz"] = scs;
    fields["pss_z"] = py_round(b.z, 1);
    fields["sss_z"] = py_round(info.sss_z, 1);
    fields["ssb_time_s"] = py_round(static_cast<double>(b.pos) / (256 * scs), 6);
    fields["crc_ok"] = info.crc_ok;
    fields["receiver_cfo_residual_hz"] = py_round_int(b.df);
    std::string status, reason;
    double conf;
    if (info.crc_ok) {
      fields["mib"] = info.mib;
      steps.push_back(fmt("CRC-24C OK (L_max %d) -> MIB decoded", info.lmax));
      status = "success";
      conf = 0.97;
    } else {
      steps.push_back("CRC-24C failed for all scrambling/L_max hypotheses");
      status = "partial";
      reason = "PCI from PSS/SSS; PBCH CRC failed (SNR/clipping/channel)";
      conf = 0.9;
    }
    const double ssb_half = 120 * scs;
    // the cell's carrier: the energy signal holding the SSB, else the SSB itself
    const Signal* host = &sig;
    for (const Signal& e : ctx.det->signals)
      if (e.f0 - 1e6 <= f_ss && f_ss <= e.f1 + 1e6 && e.bw() >= 3e6) {
        host = &e;
        break;
      }
    const double f0 = std::min(host->f0, f_ss - ssb_half), f1 = std::max(host->f1, f_ss + ssb_half);
    DecodeResult res =
        result(status, "NR SSB receiver (this package; py3gpp sequences + polar/CRC)", steps, fields, reason,
               fmt("5G NR %s PCI %d (SSB %.2f MHz)", bname.c_str(), pci, f_ss / 1e6), conf,
               fmt("CP-OFDM %g kHz SCS (SSB: PSS/SSS BPSK, PBCH QPSK)", scs / 1e3),
               fmt("NR PSS N_ID2 %d (z %.1f) on the GSCN raster, SSS N_ID1 %d", b.nid2, b.z, *info.nid1) +
                   (info.crc_ok ? ", PBCH CRC-24C passed" : ""),
               std::make_pair(f0, f1));
    std::vector<std::pair<double, double>> bursts;
    for (const Box& x : host->bursts) bursts.emplace_back(x.t0, x.t1);
    Signal o = receiver_signal(bursts, f0, f1, res, host->snr_db, name());
    o.features["ssb_span_hz"] = Json::array({py_round_int(f_ss - ssb_half), py_round_int(f_ss + ssb_half)});
    o.features["owns_band"] = true;
    for (size_t i = 0; i < std::min(o.bursts.size(), host->bursts.size()); ++i) o.bursts[i].kind = host->bursts[i].kind;
    return o;
  }
};

}  // namespace

void register_plugins(PluginList& out) { out.push_back(std::make_unique<NrDecoder>()); }

}  // namespace si::rx::nr
