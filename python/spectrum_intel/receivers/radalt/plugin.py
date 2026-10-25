"""Pipeline plugin: FMCW radio altimeters, 4200-4400 MHz (wraps .detector; research, signal gaps S8)."""
from __future__ import annotations

from spectrum_intel.decode.base import result
from spectrum_intel.decode.emission_adapter import EmissionReceiver, _clean


class RadaltDecoder(EmissionReceiver):
    name = "radalt"
    protocol = "FMCW radio altimeter"
    key = "radalt"
    bands = {"radalt"}
    method = "receivers.radalt: STFT ridge crossings of the FMCW sweep (research)"
    owns_band = True                 # the sweep crosses the whole capture: its pieces would otherwise be 'unknown'
    steps = ["gate: the capture overlaps 4200-4400 MHz, fs >= 20 MSps",
             "ridge: 128-point Hann FFT every 64 samples, the peak bin (parabolic) >= 10 dB over the frame median",
             "crossings: ridge frames joined into tracks (<= 4 bins per elapsed frame, gaps <= 8 frames), >= 5 frames "
             "on in >= half the frames spanned, linear fit R^2 >= 0.99 (a V split at its extreme), same-direction "
             "pieces within 20 us merged, >= 15 MHz covered",
             "sweep: rate = median |slope| (2-1500 MHz/ms); same-direction spacing in whole periods -> period; "
             "up / down pairs -> turnarounds, sweep bandwidth and the triangle period 2 B / r",
             "confidence 0.7, 0.8 when the bandwidth is 95-190 MHz or the period matches an ITU-R M.2059 type "
             "(20, 6.67, 7, 1 ms)"]

    MIN_RATE_HZ = 20e6

    def run_detector(self, ctx):
        from .detector import detect
        if min(ctx.cap.fs, ctx.cap.bandwidth_hz or ctx.cap.fs) < self.MIN_RATE_HZ:
            return []
        return detect(ctx.cap.iq, ctx.cap.fs, ctx.cap.fc)

    def group_key(self, e):
        return e.label

    def finish(self, sig, em, ctx):
        F = em[0].features
        sig.detection["description"] = (
            f"defined by the radalt receiver: {F['n_crossings']} straight frequency ramp(s) at "
            f"{F['rate_mhz_per_ms']} MHz/ms through the capture, not by energy alone")

    def summarise(self, em, sig):
        e = em[0]
        F = e.features
        fields = _clean(dict(F))
        what = [F["shape"], f"{F['rate_mhz_per_ms']} MHz/ms"]
        if F.get("sweep_bw_mhz") is not None:
            what.append(f"{F['sweep_bw_mhz']} MHz sweep")
        if F.get("period_ms") is not None:
            what.append(f"period {F['period_ms']} ms")
        ev = (f"{F['n_crossings']} linear crossing(s) (R^2 >= {F['r2_min']}) at {F['rate_mhz_per_ms']} MHz/ms, "
              f"seen {F['seen_mhz'][0]}-{F['seen_mhz'][1]} MHz")
        if F.get("f_top_mhz") is not None:
            ev += f"; top turnaround {F['f_top_mhz']} MHz"
        if F.get("f_bottom_mhz") is not None:
            ev += f"; bottom turnaround {F['f_bottom_mhz']} MHz"
        return result("partial", self.method, self.steps, fields,
                      "sweep measured (an FMCW altimeter carries no payload); the platform is not identified",
                      label=f"Radio altimeter (FMCW)-like ({', '.join(what)})", confidence=e.confidence,
                      modulation="FMCW", evidence=ev, span=(e.f_lo_mhz * 1e6, e.f_hi_mhz * 1e6))


PLUGINS = [RadaltDecoder()]
