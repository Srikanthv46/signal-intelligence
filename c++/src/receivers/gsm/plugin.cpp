// Pipeline plugin: GSM BCCH carriers as capture-level signals (python: receivers/gsm/plugin.py).
// Runs after lte: reads the receiver ppm it publishes (ctx.cache["ppm"]), else ctx.cache["sweep_ppm"].
#include <algorithm>
#include <cmath>
#include <set>

#include "core/bands_india.hpp"
#include "core/npcompat.hpp"
#include "core/parallel.hpp"
#include "core/pyfmt.hpp"
#include "decode/plugin.hpp"
#include "receivers/gsm/sch.hpp"

namespace si::rx::gsm {

namespace {

bool truthy(const Json* j) { return j && j->is_object() && !j->empty(); }

constexpr double kFcchPpmKnown = 5.0;     // python FCCH_PPM_KNOWN
constexpr double kFcchPpmUnknown = 30.0;  // python FCCH_PPM_UNKNOWN
constexpr double kFcchPpmAgree = 1.5;     // python FCCH_PPM_AGREE
constexpr double kFcchClean = 0.8;        // python FCCH_CLEAN

class GsmDecoder : public Decoder {
 public:
  std::string name() const override { return "gsm"; }
  std::string package() const override { return "gsm"; }
  std::string protocol() const override { return "GSM BCCH (FCCH/SCH)"; }
  bool capture_level() const override { return true; }
  std::vector<std::string> after() const override { return {"lte"}; }

  bool can_decode(const Signal& sig, DecodeContext&) const override {
    const auto keys = bands::band_keys(sig.f0, sig.f1);
    for (const auto& k : gsm_bands())
      if (std::find(keys.begin(), keys.end(), k) != keys.end()) return true;
    return false;
  }

  std::vector<Signal> capture_signals(DecodeContext& ctx) const override {
    // adjacent GSM carriers merge into one energy region, so every 200 kHz ARFCN inside a
    // GSM-band region up to 6 MHz wide is tried (LTE carriers are wider and excluded)
    std::vector<const Signal*> cand;
    for (const Signal& s : ctx.det->signals)
      if (s.flags.empty() && 100e3 <= s.bw() && s.bw() <= 6e6 && s.snr_db >= 8 && can_decode(s, ctx))
        cand.push_back(&s);
    const Json* p = ctx.cache.get<Json>("ppm");
    const Json* sw = ctx.cache.get<Json>("sweep_ppm");
    const Json* cal = truthy(p) ? p : (truthy(sw) ? sw : nullptr);
    const double ppm = cal ? (*cal)["ppm"].get<double>() : 0.0;
    struct Item {
      const Signal* sig;
      int a;
      double fr;
    };
    std::vector<Item> items;  // python _arfcns, deduplicated in order
    std::set<int> seen;
    for (const Signal* sig : cand) {
      double lo = sig->f0 + 60e3, hi = sig->f1 - 60e3;
      if (hi < lo) lo = hi = sig->fc();
      for (double f : np::arange<double>(std::ceil(lo / 200e3) * 200e3, hi + 1, 200e3)) {
        const auto af = arfcn_freq(f);
        if (!af || seen.count(af->first)) continue;
        seen.insert(af->first);
        items.push_back({sig, af->first, af->second});
      }
    }
    // python: FCCH per ARFCN (tone searched within +-5 ppm of a known receiver ppm, else +-30 ppm), in parallel
    const double span_ppm = cal ? kFcchPpmKnown : kFcchPpmUnknown;
    struct Hit {
      cfvec y;
      double ofs = 0;
      std::optional<Fcch> fc;
      double fppm = 0;
    };
    std::vector<Hit> hits = pmap(items, [&](const Item& it) {
      Hit h;
      auto sl = ctx.ch.slice(it.fr * (1 + ppm * 1e-6), 110e3, kOfs);
      h.y = std::move(sl.y);
      h.ofs = sl.fs;
      h.fc = find_fcch_cfo(h.y, span_ppm * 1e-6 * it.fr + 2e3);
      if (h.fc) h.fppm = ppm + h.fc->cfo / it.fr * 1e6;
      return h;
    });
    // a GSM BTS clock is within 0.05 ppm: every FCCH tone of one capture shows the same receiver error; the median
    // of the clean ones (purity >= 0.8, else all) is the reference, tones off it are not FCCH
    std::optional<double> ref;
    {
      dvec clean, all;
      for (const Hit& h : hits)
        if (h.fc) {
          all.push_back(h.fppm);
          if (h.fc->purity >= kFcchClean) clean.push_back(h.fppm);
        }
      dvec& v = clean.empty() ? all : clean;
      if (!v.empty()) {
        std::sort(v.begin(), v.end());
        const size_t n = v.size();
        ref = n % 2 ? v[n / 2] : (v[n / 2 - 1] + v[n / 2]) / 2.0;
      }
    }
    std::vector<size_t> idx(items.size());
    for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
    auto results = pmap(idx, [&](size_t i) {
      return carrier(*items[i].sig, items[i].a, items[i].fr, ppm, span_ppm, hits[i], ref);
    });
    std::vector<Signal> out;
    for (auto& r : results)
      if (r) out.push_back(std::move(*r));
    return out;
  }

  template <class H>
  std::optional<Signal> carrier(const Signal& sig, int a, double fr, double ppm, double span_ppm, const H& h,
                                std::optional<double> ref) const {
    const cfvec& y = h.y;
    const double ofs = h.ofs;
    std::vector<std::string> steps = {
        fmt("ARFCN %d (%.1f MHz) channelised at 4 samples/symbol, receiver %+.2f ppm applied", a, fr / 1e6, ppm),
        "FCCH: strongest +67.708 kHz tone within +-" + fmt("%g", span_ppm) +
            " ppm on 0.4 ms STFT frames, then the sliding fraction of power in that tone (142 zero bits, 546 us)"};
    std::optional<Fcch> fc = h.fc;
    Json fields = Json::object();
    fields["arfcn"] = a;
    fields["carrier_hz"] = py_round_int(fr);
    const std::pair<double, double> span{fr - 100e3, fr + 100e3};
    if (fc && std::fabs(h.fppm - *ref) > kFcchPpmAgree) {
      steps.push_back(fmt("tone at %+.2f ppm rejected: the capture's FCCH tones agree on %+.2f ppm (+-", h.fppm, *ref) +
                      fmt("%g", kFcchPpmAgree) + ")");
      fc.reset();
    }
    if (!fc && sig.bw() > 400e3) return std::nullopt;  // inside a multi-carrier region only FCCH-confirmed ARFCNs
    DecodeResult res;
    if (!fc) {
      steps.push_back("no FCCH in this 32.5 ms window (it repeats every 46 ms)");
      res = result("not_possible", "GSM FCCH/SCH receiver (this package)", steps, fields,
                   "no FCCH burst inside the capture (FCCH every 10 TDMA frames = 46 ms)",
                   fmt("GSM-like downlink ARFCN %d", a), 0.6, std::string("GMSK 270.833 ksym/s"),
                   fmt("always-on ~%.0f kHz carrier on the GSM 200 kHz raster", sig.bw() / 1e3), span);
    } else {
      const int64_t s0 = fc->start;
      const double cfo = fc->cfo, peak = fc->purity;
      // cfo = tone frequency minus 67.708 kHz = carrier offset after the ppm correction
      const double ar = -(2.0 * kPi) * cfo, inv = 1.0 / ofs;
      cdvec y2(y.size());
      for (size_t n = 0; n < y.size(); ++n) {
        const double ph = (ar * static_cast<double>(n)) * inv;
        y2[n] = np::cmul(cd(y[n].real(), y[n].imag()), cd(std::cos(ph), std::sin(ph)));
      }
      steps.push_back(fmt("FCCH at %.3f ms, tone purity %.2f, residual CFO %+.0f Hz (receiver %+.2f ppm)",
                          static_cast<double>(s0) / ofs * 1e3, peak, cfo, h.fppm));
      const int64_t sch_start = s0 + py_round_int(1250 * kSps);  // 8 timeslots x 156.25 symbols later
      const auto d = demod_sch(y2, sch_start);
      fields["fcch_time_s"] = py_round(static_cast<double>(s0) / ofs, 6);
      fields["residual_cfo_hz"] = py_round_int(cfo);
      fields["receiver_ppm_fcch"] = py_round(h.fppm, 2);
      bool ok = false;
      if (d) {
        const std::vector<int> bits = viterbi(d->soft);
        const std::vector<int> b25(bits.begin(), bits.begin() + 25), p10(bits.begin() + 25, bits.begin() + 35);
        ok = crc10_ok(b25, p10);
        steps.push_back(fmt("SCH one TDMA frame later: training sequence %d/64, Viterbi, CRC-10 %s", d->ts_match,
                            ok ? "OK" : "failed"));
        if (ok) {
          const Json f = sch_fields(b25);
          for (auto it = f.begin(); it != f.end(); ++it) fields[it.key()] = it.value();
        }
      } else {
        steps.push_back("SCH training sequence not found one frame after the FCCH");
      }
      if (ok) {
        res = result("success", "GSM FCCH/SCH receiver (this package)", steps, fields, "",
                     fmt("GSM BCCH ARFCN %d BSIC %d%d", a, fields["ncc"].get<int>(), fields["bcc"].get<int>()), 0.95,
                     std::string("GMSK 270.833 ksym/s"),
                     fmt("FCCH tone (purity %.2f) + SCH with valid CRC-10", peak), span);
      } else {
        res = result("partial", "GSM FCCH/SCH receiver (this package)", steps, fields,
                     "FCCH found (GSM BCCH carrier confirmed); SCH not decoded", fmt("GSM BCCH ARFCN %d", a), 0.9,
                     std::string("GMSK 270.833 ksym/s"),
                     fmt("FCCH frequency-correction burst (pure +67.7 kHz tone, purity %.2f)", peak), span);
      }
    }
    std::vector<std::pair<double, double>> bursts;
    for (const Box& b : sig.bursts) bursts.emplace_back(b.t0, b.t1);
    Signal o = receiver_signal(bursts, span.first, span.second, res, sig.snr_db, name());
    o.features["owns_band"] = true;
    for (size_t i = 0; i < std::min(o.bursts.size(), sig.bursts.size()); ++i) o.bursts[i].kind = sig.bursts[i].kind;
    return o;
  }
};

}  // namespace

void register_plugins(PluginList& out) { out.push_back(std::make_unique<GsmDecoder>()); }

}  // namespace si::rx::gsm
