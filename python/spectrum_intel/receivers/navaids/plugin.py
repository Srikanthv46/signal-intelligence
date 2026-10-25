"""Pipeline plugin: ILS glide slope / localizer and VOR (wraps .detector; research, signal gaps S9 / S21)."""
from __future__ import annotations

from spectrum_intel.decode.base import result
from spectrum_intel.decode.emission_adapter import EmissionReceiver, _clean

KIND = {"ils_gs": "ILS glide slope", "ils_loc": "ILS localizer", "vor": "VOR beacon"}


class NavaidsDecoder(EmissionReceiver):
    name = "navaids"
    protocol = "ILS / VOR navigation aids"
    key = "navaids"
    bands = {"ils_glidepath", "vor_ils_loc", "vor"}
    method = "receivers.navaids: carrier on the ICAO channel grid + tone fits (research)"
    owns_band = True                 # the tone sidebands would otherwise become 'unknown' boxes
    steps = ["channels inside the capture: ILS glide path 329.15-335.00 MHz and localizer 108.10-111.95 MHz "
             "(ICAO pairing), VOR 108-117.95 MHz; a carrier line >= 15 dB over the +-35 kHz median, receiver ppm "
             "applied; a second line 4-32 kHz away = two-frequency system",
             "ILS: envelope of the strongest carrier (+-1.5 kHz, 6 kSps) over whole 1/30 s windows, least squares "
             "on [1, 90 s, 150 s Hz] with the common tone scale s searched over 0.97-1.03 -> m90, m150, DDM, SDM",
             "VOR: envelope (+-12.5 kHz, 32 kSps) over whole 30 Hz periods: 30 Hz AM by least squares; 9960 Hz "
             "subcarrier (+-700 Hz) from the envelope FFT, its instantaneous frequency fitted at 30 Hz -> "
             "deviation; radial = AM phase - FM phase",
             "validity: ILS tones >= 5 sigma and SDM of the kind (glide path 0.6-1.0, localizer 0.25-0.65); VOR "
             "subcarrier >= 10 dB, AM depth 0.15-0.45, deviation 300-700 Hz; else '-like'"]

    def run_detector(self, ctx):
        from .detector import detect
        cal = ctx.cache.get("ppm") or ctx.cache.get("sweep_ppm")
        return detect(ctx.ch, cal["ppm"] if cal else 0.0)

    def group_key(self, e):
        return (e.label, e.features["channel_mhz"])

    def summarise(self, em, sig):
        e = max(em, key=lambda e: e.confidence)
        F = e.features
        fields = _clean({k: v for k, v in F.items()})
        kind = KIND[e.label]
        ok = F["verified_structure"]
        if e.label == "vor":
            what = f"{F['channel_mhz']:.2f} MHz"
            ev = (f"carrier {F['line_db']} dB over its surroundings; 30 Hz AM {F['am_depth'] * 100:.0f} %, 9960 Hz "
                  f"subcarrier +{F['subcarrier_db']} dB with 30 Hz FM +-{F['deviation_hz']} Hz; radial "
                  f"{F['radial_deg']} deg from the station")
            mod = "AM 30 Hz + 9960 Hz subcarrier FM (VOR)"
        else:
            what = f"{F['channel_mhz']:.2f} MHz, paired {'LOC' if e.label == 'ils_gs' else 'GS'} {F['paired_mhz']:.2f} MHz"
            ev = (f"carrier {F['line_db']} dB over its surroundings; 90 Hz {F['m90'] * 100:.1f} % and 150 Hz "
                  f"{F['m150'] * 100:.1f} % (tone scale {F['tone_scale']}); DDM {F['ddm']:+.4f} ({F['side']}), "
                  f"SDM {F['sdm']:.3f}")
            mod = "AM 90 / 150 Hz (ILS)"
        if F.get("two_frequency_hz"):
            what += ", two-frequency"
            ev += f"; second carrier {F['two_frequency_hz']} Hz away (two-frequency system)"
        label = f"{kind} ({what})" if ok else f"{kind}-like ({what})"
        return result("partial", self.method, self.steps, fields,
                      "navigation tones measured (an ILS / VOR carries no payload); ident (Morse 1020 Hz) needs a "
                      "dwell of several seconds",
                      label=label, confidence=e.confidence, modulation=mod, evidence=ev,
                      span=(min(x.f_lo_mhz for x in em) * 1e6, max(x.f_hi_mhz for x in em) * 1e6))


PLUGINS = [NavaidsDecoder()]
