"""Pipeline plugin: GPS L1 C/A satellites, GNSS jamming / spoofing indicators (wraps .acquire; research, S3)."""
from __future__ import annotations

from spectrum_intel.decode.base import Decoder, receiver_signal, result

from .codes import F_L1


class GnssDecoder(Decoder):
    name = "gnss"
    protocol = "GPS L1 C/A"
    capture_level = True
    after = ("lte",)                 # reads the receiver ppm the LTE receiver measured
    steps = ["gate: 1575.42 MHz +-1.1 MHz in the usable band, fs >= 2.1 MSps, >= 16 ms",
             "block FFT cropped to +-1.024 MHz around L1 / (1 + ppm) -> 2.048 MSps; gate at 1 sample / chip",
             "32 PRNs (IS-GPS-200 Gold codes, ICD self-test octals), Doppler 125 Hz bins over +-8.2 kHz with the "
             "receiver ppm (+-25 kHz without), 4 ms coherent folds x 4 non-coherent, one FFT per fold shared by "
             "all PRNs; threshold chi2(2M) at P_fa 1e-3 per PRN per block",
             "each hit refined at 2 samples / chip over the whole block: code phase, Doppler, C/N0 (parabolic "
             "straddle correction + measured implementation loss)",
             "indicators: in-band J/N over the capture's floor outside 1559-1610 MHz, narrow lines near L1 (receiver "
             "lines excluded); spoofing S1 C/N0 > 55, S2 uniform C/N0 (>= 4 PRNs > 45 dB-Hz within 1.5 dB), S3 J/N "
             ">= 1.5 dB with >= 4 strong PRNs, S4 second correlation peaks in >= 2 PRNs",
             "swept jammer (personal privacy device): 32-point FFT ridge over 400 us, sweep period from its "
             "autocorrelation (0.5-50 us), fold coherence >= 0.6 and span >= 1 MHz (fs >= 30 MSps)"]

    def can_decode(self, sig, ctx):
        return False

    def capture_signals(self, ctx):
        from .acquire import analyse
        cal = ctx.cache.get("ppm") or ctx.cache.get("sweep_ppm")
        r = analyse(ctx.ch, cal["ppm"] if cal else None,
                    [ln["freq_hz"] for ln in ctx.cache.get("receiver_lines", [])])
        if r is None:
            return []
        sats, ind = r["satellites"], r["indicators"]
        n = len(sats)
        names = [s["prn"] for s in sats]
        span = (F_L1 - 1.023e6, F_L1 + 1.023e6)
        if ind["spoofing"]["verdict"]:
            k = ind["spoofing"]["n"]
            label, conf = f"GNSS spoofing-like (GPS L1 C/A, {n} satellites)", 0.75 if k >= 3 else 0.6
            ev = (f"{n} GPS satellites {names}; spoofing indicators "
                  f"{[k for k, v in ind['spoofing'].items() if v is True and k != 'verdict']}")
        elif ind["jamming"]["verdict"]:
            itf = r["interference"]
            sw = itf.get("swept")
            label, conf = ("GNSS jamming-like (GPS L1, swept PPD)", 0.85) if sw else ("GNSS jamming-like (GPS L1)", 0.6)
            ev = (f"in-band J/N {itf.get('jn_db')} dB, {len(itf.get('lines') or [])} narrow line(s) near L1; "
                  f"{n} satellite(s) acquired")
            if sw:
                ev += (f"; swept carrier ({sw['shape']}) every {sw['period_us']} us over {sw['span_mhz']} MHz, centre "
                       f"{sw['centre_offset_mhz']:+} MHz from L1{'' if sw['covers_l1'] else ' (does not cross L1)'}")
        elif n >= 2 or (n == 1 and sats[0]["cn0_dbhz"] >= 40):
            label, conf = f"GNSS: GPS L1 C/A ({n} satellite{'s' if n > 1 else ''})", 0.9
            ev = (f"{n} PRN(s) above the chi2 threshold (P_fa 1e-3 per PRN): "
                  + ", ".join(f"PRN {s['prn']} {s['cn0_dbhz']} dB-Hz" for s in sats))
        else:
            ctx.cache["gnss_note"] = (f"L1 searched: {n} satellite(s) above threshold, no jamming / spoofing "
                                      f"indicator")
            return []
        res = result("partial", "receivers.gnss: GPS L1 C/A acquisition (research)", self.steps, r,
                     "satellites acquired by code correlation; the navigation message is not decoded",
                     label=label, confidence=conf, modulation="BPSK(1) DSSS 1.023 Mcps (GPS C/A)", evidence=ev,
                     span=span)
        dur = ctx.ch.L / ctx.ch.fs
        return [receiver_signal([(0.0, dur)], span[0], span[1], res, None, self.name)]


PLUGINS = [GnssDecoder()]
