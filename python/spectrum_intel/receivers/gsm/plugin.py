"""Pipeline plugin: GSM BCCH carriers as capture-level signals (owns its band)."""
from __future__ import annotations

import numpy as np
from spectrum_intel.core.bands_india import band_keys

from spectrum_intel.decode.base import Decoder, receiver_signal, result
from .sch import GSM_BANDS, OFS, SPS, _arfcn_freq, crc10_ok, demod_sch, find_fcch_cfo, sch_fields, viterbi

FCCH_PPM_KNOWN = 5.0                # tone search around a known receiver ppm (TCXO drift: 2-3 ppm seen)
FCCH_PPM_UNKNOWN = 30.0             # ... and without one
FCCH_PPM_AGREE = 1.5                # an FCCH tone off the capture's consensus by more is not an FCCH
FCCH_CLEAN = 0.8                    # purity of the tones that set the consensus


class GsmDecoder(Decoder):
    name = "gsm"
    after = ("lte",)                 # reads the receiver ppm the lte plugin measures (ctx.cache["ppm"])
    protocol = "GSM BCCH (FCCH/SCH)"
    capture_level = True

    def can_decode(self, sig, ctx):
        return bool(set(band_keys(sig.f0, sig.f1)) & GSM_BANDS)

    @staticmethod
    def _arfcns(cand):
        for sig in cand:
            lo, hi = sig.f0 + 60e3, sig.f1 - 60e3
            if hi < lo:
                lo = hi = sig.fc
            for f in np.arange(np.ceil(lo / 200e3) * 200e3, hi + 1, 200e3):
                a, fr = _arfcn_freq(f)
                if a is not None:
                    yield sig, a, fr

    def capture_signals(self, ctx):
        # adjacent GSM carriers merge into one energy region, so every 200 kHz ARFCN inside
        # a GSM-band region up to 6 MHz wide is tried (LTE carriers are wider and excluded)
        cand = [s for s in ctx.det.signals if not s.flags and 100e3 <= s.bw <= 6e6 and s.snr_db >= 8
                and self.can_decode(s, ctx)]
        cal = ctx.cache.get("ppm") or ctx.cache.get("sweep_ppm")
        ppm = cal["ppm"] if cal else 0.0
        span_ppm = FCCH_PPM_KNOWN if cal else FCCH_PPM_UNKNOWN
        items, seen = [], set()
        for sig, a, fr in self._arfcns(cand):
            if a in seen:
                continue
            seen.add(a)
            y, ofs = ctx.ch.slice(fr * (1 + ppm * 1e-6), 110e3, OFS)
            fc = find_fcch_cfo(y, span_ppm * 1e-6 * fr + 2e3)
            items.append((sig, a, fr, y, ofs, fc, None if fc is None else ppm + fc[1] / fr * 1e6))
        # a GSM BTS clock is within 0.05 ppm, so every FCCH tone of one capture shows the same receiver error:
        # the median of the clean ones (purity >= 0.8, else all) is the reference, outliers are not FCCH
        hits = [(it[6], it[5][2]) for it in items if it[5] is not None]
        ref = None
        if hits:
            clean = [p for p, q in hits if q >= FCCH_CLEAN] or [p for p, _ in hits]
            ref = float(np.median(clean))
        out = []
        for sig, a, fr, y, ofs, fc, fppm in items:
            steps = [f"ARFCN {a} ({fr / 1e6:.1f} MHz) channelised at 4 samples/symbol, receiver "
                     f"{ppm:+.2f} ppm applied",
                     f"FCCH: strongest +67.708 kHz tone within +-{span_ppm:g} ppm on 0.4 ms STFT frames, then the "
                     "sliding fraction of power in that tone (142 zero bits, 546 us)"]
            fields = {"arfcn": a, "carrier_hz": round(fr)}
            if fc is not None and abs(fppm - ref) > FCCH_PPM_AGREE:
                steps.append(f"tone at {fppm:+.2f} ppm rejected: the capture's FCCH tones agree on {ref:+.2f} ppm "
                             f"(+-{FCCH_PPM_AGREE:g})")
                fc = None
            if fc is None and sig.bw > 400e3:
                continue                       # inside a multi-carrier region only FCCH-confirmed ARFCNs count
            if fc is None:
                steps.append("no FCCH in this 32.5 ms window (it repeats every 46 ms)")
                res = result("not_possible", "GSM FCCH/SCH receiver (this package)", steps, fields,
                             reason="no FCCH burst inside the capture (FCCH every 10 TDMA frames = 46 ms)",
                             label=f"GSM-like downlink ARFCN {a}", confidence=0.6,
                             modulation="GMSK 270.833 ksym/s", span=(fr - 100e3, fr + 100e3),
                             evidence=f"always-on ~{sig.bw / 1e3:.0f} kHz carrier on the GSM 200 kHz raster")
            else:
                s0, cfo, peak = fc
                # cfo = tone frequency minus 67.708 kHz = carrier offset after the ppm correction
                y2 = y * np.exp(-2j * np.pi * cfo * np.arange(len(y)) / ofs)
                steps.append(f"FCCH at {s0 / ofs * 1e3:.3f} ms, tone purity {peak:.2f}, residual CFO {cfo:+.0f} Hz "
                             f"(receiver {fppm:+.2f} ppm)")
                sch_start = s0 + int(round(1250 * SPS))         # 8 timeslots x 156.25 symbols later
                d = demod_sch(y2, sch_start)
                fields.update(fcch_time_s=round(s0 / ofs, 6), residual_cfo_hz=round(cfo),
                              receiver_ppm_fcch=round(fppm, 2))
                ok = False
                if d is not None:
                    soft, inf = d
                    bits = viterbi(soft)
                    ok = crc10_ok(bits[:25], bits[25:35])
                    steps.append(f"SCH one TDMA frame later: training sequence {inf['ts_match']}/64, Viterbi, "
                                 f"CRC-10 {'OK' if ok else 'failed'}")
                    if ok:
                        fields.update(sch_fields(bits[:25]))
                else:
                    steps.append("SCH training sequence not found one frame after the FCCH")
                if ok:
                    res = result("success", "GSM FCCH/SCH receiver (this package)", steps, fields,
                                 label=f"GSM BCCH ARFCN {a} BSIC {fields['ncc']}{fields['bcc']}", confidence=0.95,
                                 modulation="GMSK 270.833 ksym/s", span=(fr - 100e3, fr + 100e3),
                                 evidence=f"FCCH tone (purity {peak:.2f}) + SCH with valid CRC-10")
                else:
                    res = result("partial", "GSM FCCH/SCH receiver (this package)", steps, fields,
                                 reason="FCCH found (GSM BCCH carrier confirmed); SCH not decoded",
                                 label=f"GSM BCCH ARFCN {a}", confidence=0.9, modulation="GMSK 270.833 ksym/s",
                                 span=(fr - 100e3, fr + 100e3),
                                 evidence=f"FCCH frequency-correction burst (pure +67.7 kHz tone, purity {peak:.2f})")
            o = receiver_signal([(b.t0, b.t1) for b in sig.bursts], fr - 100e3, fr + 100e3, res, sig.snr_db,
                                self.name)
            o.features["owns_band"] = True
            for x, yb in zip(o.bursts, sig.bursts):
                x.kind = yb.kind
            out.append(o)
        return out


PLUGINS = [GsmDecoder()]
