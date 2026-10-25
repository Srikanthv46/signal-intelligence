"""Pipeline plugin: NR cells (SSB + MIB) as capture-level signals (owns its band)."""
from __future__ import annotations

import numpy as np
from spectrum_intel.core.bands_india import band_keys, nr_arfcn
from spectrum_intel.core.parallel import pmap

from spectrum_intel.decode.base import Decoder, receiver_signal, result
from .ssb import NR_BANDS, _SSBSlice, _bands_for, _gscn, decode_ssb_multi


class NrDecoder(Decoder):
    name = "nr_ssb"
    after = ("lte",)                 # reads the receiver ppm the lte plugin measures (ctx.cache["ppm"])
    protocol = "3GPP 5G NR SS/PBCH block"
    capture_level = True

    def can_decode(self, sig, ctx):
        return bool(set(band_keys(sig.f0, sig.f1)) & set(NR_BANDS))

    def capture_signals(self, ctx):
        cand = [s for s in ctx.det.signals if not s.flags and s.snr_db >= 3 and s.bw >= 3e6 and self.can_decode(s, ctx)]
        if not cand:
            return []
        cal = ctx.cache.get("ppm") or ctx.cache.get("sweep_ppm")
        ppm = cal["ppm"] if cal else 0.0
        seen, out, todo = set(), [], []
        cap = ctx.cap
        lo_cap, hi_cap = cap.fc - 0.49 * cap.fs, cap.fc + 0.49 * cap.fs     # AD936x passband ~ +-0.49 fs
        for sig in cand:
            key = next(k for k in NR_BANDS if k in band_keys(sig.f0, sig.f1))
            bname, scs = NR_BANDS[key]
            # the SSB need not sit inside the energy box (it may be at the carrier edge or
            # in a fragment): search the raster over the band's part of the capture
            blo = max(lo_cap, min(b[0] for b in _bands_for(key))) + 120 * scs
            bhi = min(hi_cap, max(b[1] for b in _bands_for(key))) - 120 * scs
            for f_ss in _gscn(blo, bhi):
                if round(f_ss / 10e3) in seen:
                    continue
                seen.add(round(f_ss / 10e3))
                todo.append((sig, f_ss, bname, scs))

        def probe(item):
            sig, f_ss, bname, scs = item
            # +-12 kHz: the receiver TCXO drifts ~1-2 ppm during a sweep (3.5 GHz -> ~7 kHz)
            sl = _SSBSlice(ctx.ch, f_ss * (1 + ppm * 1e-6), scs)
            return sl.best(np.arange(-12e3, 12.01e3, 3e3)), sl, sig, f_ss, bname, scs
        scan = pmap(probe, todo)                    # GSCN hypotheses are independent
        if not scan:
            return []
        # z of a loaded TDD band has a heavy tail: an SSB must stand well clear of the scan
        z_thr = max(40.0, 2.5 * float(np.median([x[0]["z"] for x in scan])))
        for b, sl, sig, f_ss, bname, scs in sorted(scan, key=lambda x: -x[0]["z"]):
            if b["z"] < z_thr:
                break
            if any(abs(o.features["ssb_center_hz"] - f_ss) < 240 * scs for o in out):
                continue
            fine = sl.best(np.arange(b["df"] - 2e3, b["df"] + 2.01e3, 500.0))
            b = fine if fine["z"] >= b["z"] else b
            pci, info = decode_ssb_multi(sl, b)
            if pci is None:
                continue
            if not info.get("crc_ok") and (info["sss_z"] < 4.5 or info["sss_z"] / max(info["sss_z2"], 1e-3) < 1.6):
                continue
            b = dict(b, nid2=info["nid2"])
            o = self._signal(ctx, sig, f_ss, scs, bname, b, pci, info, cal)
            o.features["ssb_center_hz"] = round(f_ss)
            out.append(o)
        return out

    def _signal(self, ctx, sig, f_ss, scs, bname, b, pci, info, cal):
        steps = [f"SSB centre hypotheses on the GSCN raster inside {sig.f0 / 1e6:.2f}-{sig.f1 / 1e6:.2f} MHz "
                 f"({scs / 1e3:g} kHz SCS), receiver error {cal['ppm']:+.2f} ppm applied" if cal else
                 "SSB centre hypotheses on the GSCN raster (ppm unknown)",
                 f"PSS N_ID2 {b['nid2']} z {b['z']:.1f} at {b['pos'] / (256 * scs) * 1e3:.3f} ms",
                 f"SSS (differential, channel-free; summed over SSB instances) N_ID1 {info['nid1']} "
                 f"(z {info['sss_z']:.1f} vs {info['sss_z2']:.1f}) -> PCI {pci}",
                 f"{info['n_ssb_instances']} SSB instance(s) (beam sweep / bursts) combined: SSS scores summed; "
                 f"per-SSB PBCH DMRS -> i_SSB {info.get('ssb_indices')}, DMRS channel estimate, MRC",
                 "PBCH LLRs descrambled per i_SSB and summed per half-frame; polar decode + CRC-24C "
                 "(py3gpp.nrBCHDecode)"]
        fields = {"pci": pci, "nid1": info["nid1"], "nid2": b["nid2"], "co_channel_cells": info.get("other_cells"), "ssb_index": info["ibar_ssb"],
                  "ssb_center_hz": round(f_ss), "gscn_freq_hz": round(f_ss), "nr_arfcn": nr_arfcn(f_ss),
                  "ssb_scs_hz": scs, "pss_z": round(b["z"], 1), "sss_z": round(info["sss_z"], 1),
                  "ssb_time_s": round(b["pos"] / (256 * scs), 6), "crc_ok": info["crc_ok"],
                  "receiver_cfo_residual_hz": round(b["df"])}
        if info["crc_ok"]:
            fields["mib"] = info["mib"]
            steps.append(f"CRC-24C OK (L_max {info['lmax']}) -> MIB decoded")
            status, reason, conf = "success", "", 0.97
        else:
            steps.append("CRC-24C failed for all scrambling/L_max hypotheses")
            status, reason, conf = "partial", "PCI from PSS/SSS; PBCH CRC failed (SNR/clipping/channel)", 0.9
        ssb_half = 120 * scs
        res = result(status, "NR SSB receiver (this package; py3gpp sequences + polar/CRC)", steps, fields, reason,
                     label=f"5G NR {bname} PCI {pci} (SSB {f_ss / 1e6:.2f} MHz)", confidence=conf,
                     modulation=f"CP-OFDM {scs / 1e3:g} kHz SCS (SSB: PSS/SSS BPSK, PBCH QPSK)",
                     evidence=f"NR PSS N_ID2 {b['nid2']} (z {b['z']:.1f}) on the GSCN raster, SSS N_ID1 "
                              f"{info['nid1']}" + (", PBCH CRC-24C passed" if info["crc_ok"] else ""),
                     span=(sig.f0, sig.f1))
        # the cell's carrier: the energy signal holding the SSB, else the SSB itself
        host = next((e for e in ctx.det.signals if e.f0 - 1e6 <= f_ss <= e.f1 + 1e6 and e.bw >= 3e6), sig)
        f0, f1 = min(host.f0, f_ss - ssb_half), max(host.f1, f_ss + ssb_half)
        res["span"] = (f0, f1)
        out = receiver_signal([(x.t0, x.t1) for x in host.bursts], f0, f1, res, host.snr_db, self.name)
        out.features["ssb_span_hz"] = [round(f_ss - ssb_half), round(f_ss + ssb_half)]
        out.features["owns_band"] = True
        for x, y in zip(out.bursts, host.bursts):
            x.kind = y.kind
        return out


PLUGINS = [NrDecoder()]
