"""Pipeline plugin: LTE cells as capture-level signals (owns its band)."""
from __future__ import annotations

from spectrum_intel.core.bands_india import band_keys, lte_earfcn_dl

from spectrum_intel.decode.base import Decoder, receiver_signal, result
from .cellsearch import (LTE_BANDS, OCC, OFS, ROOTS, bandwidth_from_energy, raster_grid, search_cells, sss_decisive,
                         wide_ppm_search)


# energy SNR of a candidate: the PSS correlates over many symbols, so a carrier only a few dB over the floor in
# energy still gives a decisive PSS/SSS (B28 at 760.5 MHz: 3.6 dB energy, PSS z 77)
CAND_SNR_DB = 3.0


class LteDecoder(Decoder):
    name = "lte"
    protocol = "3GPP LTE (E-UTRA) downlink"
    capture_level = True

    def can_decode(self, sig, ctx):
        return bool(set(band_keys(sig.f0, sig.f1)) & set(LTE_BANDS))

    def _ppm(self, ctx, cand):
        if "ppm" in ctx.cache:
            return ctx.cache["ppm"]
        cal = ctx.cache.get("sweep_ppm")
        if cal is None and cand:
            top = sorted(cand, key=lambda s: -s.snr_db)[:3]
            cal = wide_ppm_search(ctx.ch, [round(s.fc / 100e3) * 100e3 for s in top])
            if cal:
                cal["source"] = "this capture: strongest PSS, +-500 kHz search"
        ctx.cache["ppm"] = cal
        return cal

    def capture_signals(self, ctx):
        # Continuous carriers >= 6 dB (energy) are searched as always. Bursty or weaker energy is searched too
        # (a lightly loaded downlink carries little more than its reference signals, so its energy detection
        # is split into bursts; the PSS gains over many symbols), but there only an SSS-resolved cell counts.
        inb = [s for s in ctx.det.signals
               if not s.flags and s.snr_db >= CAND_SNR_DB and s.bw >= 1.0e6 and self.can_decode(s, ctx)]
        base = [s for s in inb if s.snr_db >= 6 and any(b.kind == "continuous" for b in s.bursts)]
        extra = sorted((s for s in inb if s not in base), key=lambda s: -s.bw)
        if not inb:
            return []
        cal = self._ppm(ctx, base or inb)
        ppm = cal["ppm"] if cal else 0.0
        out, seen, done = [], set(), set()
        for sig, strict in [(s, False) for s in base] + [(s, True) for s in extra]:
            cells = search_cells(ctx, sig, ppm, done if strict else None, ppm_known=cal is not None,
                                 require_sss=strict)
            if not strict:
                done.update(int(round(c / 100e3)) for c in raster_grid(sig))
            for c in cells:
                s = c["sss"]
                if self._band(c["raster_hz"]) == "?":
                    continue                     # an LTE downlink carrier lies in an LTE band (Herelink at 2.41 GHz)
                key = (round(c["raster_hz"] / 100e3), c["nid2"])
                if key in seen:
                    continue
                seen.add(key)
                out.append(self._signal(ctx, sig, c, s, cal))
        return out

    @staticmethod
    def _band(raster):
        ea = lte_earfcn_dl(raster)
        return f"B{ea[0]}" if ea else next((v for k, v in LTE_BANDS.items() if k in band_keys(raster, raster)), "?")

    def _signal(self, ctx, sig, c, s, cal):
        steps = [f"carrier hypotheses: 100 kHz EARFCN raster inside {sig.f0 / 1e6:.2f}-{sig.f1 / 1e6:.2f} MHz, "
                 f"receiver error {cal['ppm']:+.2f} ppm applied ({cal['source']})" if cal else
                 "carrier hypotheses: 100 kHz EARFCN raster (receiver ppm unknown, 0 assumed)",
                 "PSS: 3.84 MSps slice, 3 roots, 0 Hz (receiver error known, < 3 kHz residual) or +-5 kHz, "
                 "normalised correlation folded on 5 ms; z >= 15",
                 f"PSS N_ID2 {c['nid2']} z {c['pss_z']:.1f} over {len(c['inst'])} half-frames; CFO refined +-3 kHz",
                 "SSS: 168 N_ID1 x {SF0,SF5} at FDD/TDD x normal/extended CP, maximum-ratio combined with the PSS "
                 "channel estimate (smoothed over 225 kHz), accumulated over half-frames"]
        raster = c["raster_hz"]
        ea = lte_earfcn_dl(raster)
        band = self._band(raster)
        bw = bandwidth_from_energy(ctx.det, c["carrier_hz"])
        half = (OCC[bw] if bw else 1.08e6) / 2
        fields = {"carrier_hz": round(c["carrier_hz"]), "raster_carrier_hz": round(raster),
                  "earfcn": ea[1] if ea else None, "band": band, "nid2": c["nid2"], "pss_fold_z": round(c["pss_z"], 1),
                  "n_half_frames": len(c["inst"]), "pss_isolation": round(c["pss_isolation"] or 0, 2),
                  "receiver_freq_error_hz": round(c["err_hz"]),
                  "receiver_freq_error_ppm": round(c["err_hz"] / raster * 1e6, 2),
                  "channel_bandwidth_hz_from_energy": bw,
                  "mib": "not decoded (PBCH decoder not implemented)"}
        ok = sss_decisive(s)
        if ok:
            pci = 3 * s["nid1"] + c["nid2"]
            # subframe 0 among the PSS instances found (none when they all fall on subframe-5 half-frames)
            p0 = next((p for k, p in c["inst"] if (k + s["parity"]) % 2 == 0), None)
            sf0 = round(p0 / OFS, 6) if p0 is not None else None
            fields.update(pci=pci, nid1=s["nid1"], duplex=s["duplex"], cyclic_prefix=s["cp"],
                          sss_z=round(s["z"], 1), sss_z_runner_up=round(s["z2"], 1),
                          subframe0_time_s=sf0)
            steps.append(f"SSS N_ID1 {s['nid1']} ({s['duplex']}, {s['cp']} CP): z {s['z']:.1f} vs runner-up "
                         f"{s['z2']:.1f} -> PCI {pci}")
            lab = (f"LTE {band} {s['duplex']} downlink PCI {pci}" + (f" EARFCN {ea[1]}" if ea else "") +
                   (f" {bw / 1e6:g} MHz" if bw else ""))
            ev = (f"LTE PSS (root {ROOTS[c['nid2']]}) every 5 ms z {c['pss_z']:.1f}; SSS N_ID1 {s['nid1']} "
                  f"z {s['z']:.1f} vs {s['z2']:.1f}; carrier on EARFCN raster within {c['err_hz'] / 1e3:+.1f} kHz")
            res = result("partial", "LTE PSS/SSS cell search (this package)", steps, fields,
                         reason="physical cell identity decoded; MIB/SIB (PBCH/PDSCH) decoding not implemented",
                         label=lab, confidence=0.95, modulation="OFDMA 15 kHz SCS (PSS/SSS sync)", evidence=ev,
                         span=(c["carrier_hz"] - half, c["carrier_hz"] + half))
        else:
            steps.append("SSS not decisive" + (f" (best z {s['z']:.1f} vs {s['z2']:.1f})" if s else ""))
            res = result("partial", "LTE PSS/SSS cell search (this package)", steps, fields,
                         reason="LTE PSS found but SSS ambiguous: N_ID1 not resolved",
                         label=f"LTE {band} downlink (PSS N_ID2 {c['nid2']})", confidence=0.85,
                         modulation="OFDMA 15 kHz SCS", span=(c["carrier_hz"] - half, c["carrier_hz"] + half),
                         evidence=f"LTE PSS root {ROOTS[c['nid2']]} every 5 ms (z {c['pss_z']:.1f})")
        tdd = ok and s["duplex"] == "TDD"
        bursts = [(b.t0, b.t1) for b in sig.bursts] if tdd else [(0.0, ctx.cap.duration)]
        out = receiver_signal(bursts, res["span"][0], res["span"][1], res, sig.snr_db, self.name)
        for b in out.bursts:
            b.kind = "continuous" if not tdd else "burst"
        out.features["owns_band"] = True
        return out


PLUGINS = [LteDecoder()]
