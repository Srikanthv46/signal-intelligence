"""Pipeline plugin: swept-carrier jammers (wraps .detector)."""
from __future__ import annotations

from spectrum_intel.core.bands_india import band_keys
from spectrum_intel.decode.base import Decoder, result
from spectrum_intel.decode.emission_adapter import EmissionReceiver, _clean

# a continuous swept carrier here is the band's own FMCW radar, not a jammer
FMCW_BANDS = {"radalt"}


class JammerDecoder(EmissionReceiver):
    name = "jammer"
    protocol = "swept-carrier jammer"
    key = "jammer"
    method = "receivers.jammer (swept-carrier detector)"
    owns_band = True                 # the swept carrier's clipping products / band crossings would otherwise
                                     # become 'unknown' boxes all over its band
    steps = ["occupancy: Hann STFT power, fs/64 bins, ~1 us rows; cyclic spectrum of per-bin log occupancy, "
             "whitened per bin and summed -> up to 3 sweep-rate candidates (lines >= 10 dB over the local median)",
             "period: the shortest multiple (<= 8) of the line period whose fold has >= 80 % of the best phase "
             "structure, refused if 2-5x that period has more (pulsed); rate corrected from the fold drift "
             "across 2-4 block segments",
             "fold on the period: one dominant carrier per phase moving continuously (<= 1 flyback, <= 2 band "
             "exits), in band >= 30 % of the period, spanning >= 10 MHz -> sweep rate, extent, shape, trajectory"]

    MIN_RATE_HZ = 20e6               # schema/requirements.json jammer.min_sample_rate_hz: a >= 10 MHz sweep

    def run_detector(self, ctx):
        from .detector import detect
        rate = min(ctx.cap.fs, ctx.cap.bandwidth_hz or ctx.cap.fs)
        if rate < self.MIN_RATE_HZ:
            return []
        return detect(ctx.cap.iq, ctx.cap.fs, ctx.cap.fc)

    def in_range(self, lo, hi):
        return True                  # jammers are not tied to an allocation

    def group_key(self, e):
        return e.center_mhz

    def finish(self, sig, em, ctx):
        F = em[0].features
        sig.detection["description"] = (
            f"defined by the jammer receiver: a cyclic line in the per-bin STFT occupancy at the sweep rate, "
            f"confirmed by folding the STFT on the {F['sweep_period_us']:.2f} us period (one carrier moving "
            f"continuously through {F['sweep_span_mhz']:.1f} MHz), not by energy alone")

    def summarise(self, em, sig):
        e = em[0]
        F = e.features
        keys = ("sweep_rate_hz", "sweep_period_us", "sweep_low_mhz", "sweep_high_mhz", "sweep_span_mhz",
                "beyond_capture", "shape", "in_band_fraction", "flybacks", "band_exits", "dominance_db",
                "fold_peak_db", "trajectory_us_mhz")
        fields = _clean({k: F[k] for k in keys if k in F})
        beyond = " and ".join(F.get("beyond_capture") or [])
        extent = (f"{F['sweep_low_mhz']:.1f}-{F['sweep_high_mhz']:.1f} MHz seen"
                  + (f", continuing {beyond} the captured band" if beyond else ""))
        fmcw = bool(FMCW_BANDS & set(band_keys(e.f_lo_mhz * 1e6, e.f_hi_mhz * 1e6)))
        label = "Swept carrier (FMCW-like)" if fmcw else "Jammer (swept carrier)"
        conf = 0.6 if fmcw else e.confidence
        return result("success", self.method, self.steps, fields,
                      "sweep measured (a jammer carries no payload; the repeating sweep is the verification)",
                      label=label, confidence=conf, modulation=F.get("modulation", "swept CW"),
                      evidence=f"one carrier swept {F['shape']} every {F['sweep_period_us']:.2f} us "
                               f"({F['sweep_rate_hz']:.1f} Hz), {extent}; in band "
                               f"{F['in_band_fraction']:.0%} of the period, {F['dominance_db']:.1f} dB over the "
                               f"swept span" + ("; in the radio-altimeter band: FMCW radar, not called a jammer"
                                                if fmcw else ""),
                      span=(e.f_lo_mhz * 1e6, e.f_hi_mhz * 1e6))


class BarrageDecoder(Decoder):
    """Barrage-noise / noise-FM jammers (research, signal gaps S12): runs on the energy detector's wide continuous
    boxes that no receiver explained (barrage.py)."""
    name = "jammer_barrage"
    protocol = "barrage-noise / noise-FM jammer"
    method = "receivers.jammer.barrage: plateau, Gaussianity, no CP / cyclic structure, stationarity (research)"
    steps = ["gate: an unexplained energy box >= 1 MHz wide, continuous over the block, >= 6 dB",
             "plateau (block FFT, 64 sub-bands of the inner 90 %): >= 90 % within +-3 dB, outside strips >= 6 dB "
             "lower, no narrow line",
             "statistics at 1.25 x the width: envelope k = E|x|^4 / E^2|x|^2, I / Q kurtosis",
             "no structure: Hann-tapered autocorrelation (cyclic prefix) < 0.025 over 1-1000 us, no cyclic line in "
             "|x|^2, x^2 or the instantaneous frequency, power per 1 ms within +-1 dB",
             "barrage: Gaussian + no structure + flat plateau, >= 10 MHz or edges on an allocation (covers_band); "
             "noise-FM: constant envelope, wide Gaussian-like instantaneous frequency, no structure"]

    def can_decode(self, sig, ctx):
        return sig.bw >= 1e6 and sig.snr_db is not None and sig.snr_db >= 6.0 \
            and any(b.kind == "continuous" for b in sig.bursts)

    def decode(self, sig, ctx):
        from .barrage import analyse
        r = analyse(ctx.ch, sig.f0, sig.f1)
        if r is None or r["class"] is None:
            why = "not a noise jammer"
            if r is not None:
                why += (f" (envelope k {r['envelope_k']}, cyclic-prefix correlation {r['cp_corr']}, plateau "
                        f"{'flat' if r['plateau_ok'] else 'not flat'})")
            return result("not_possible", self.method, self.steps, reason=why)
        cb = r["covers_band"]
        if r["class"] == "barrage":
            label, conf = "Jammer (barrage noise)-like", 0.75 if cb else 0.6
            ev = (f"flat Gaussian plateau {r['plateau_bw_mhz']} MHz (envelope k {r['envelope_k']}, I / Q kurtosis "
                  f"{r['kurtosis_i']} / {r['kurtosis_q']}), edges {r['edge_db']} dB; no cyclic prefix (max "
                  f"{r['cp_corr']}), no cyclic line, stationary within {r['stationarity_db']} dB")
            mod = "noise"
        else:
            label, conf = "Jammer (noise-FM)-like", 0.7
            ev = (f"constant envelope (k {r['envelope_k']}) over {r['plateau_bw_mhz']} MHz, instantaneous frequency "
                  f"spread {r['if_std_mhz']} MHz with kurtosis {r['if_kurtosis']}; no cyclic line, stationary")
            mod = "noise-FM"
        if cb:
            ev += f"; edges on the {', '.join(cb)} allocation"
        return result("partial", self.method, self.steps, _clean(r),
                      "a noise jammer carries no payload: the statistics are the evidence",
                      label=label, confidence=conf, modulation=mod, evidence=ev, span=(sig.f0, sig.f1))


PLUGINS = [JammerDecoder(), BarrageDecoder()]
