"""Adapter for receivers whose detect(x, fs, fc, ctx) returns core.emission.Emission
records (one per burst): EmissionReceiver groups the bursts into pipeline signals
(one per link / channel), summarises them into a decode result and keeps the
Emission objects in signal.features["_emissions"] for later plugins (emitter
association). A plugin subclasses it and implements run_detector(); the hooks
keep(), group_key(), summarise(), finish() and in_range() refine the behaviour."""
from __future__ import annotations

from collections import Counter

import numpy as np

from .base import Decoder, overlaps, receiver_signal, result
from ..core.bands_india import band_keys


def _in_ranges(lo, hi, ranges):
    return any(lo < b and hi > a for a, b in ranges)


def _span(e):
    return e.start_us * 1e-6, e.end_us * 1e-6, e.f_lo_mhz * 1e6, e.f_hi_mhz * 1e6


def _clean(d):
    import numpy as np
    if isinstance(d, dict):
        return {k: _clean(v) for k, v in d.items()}
    if isinstance(d, (list, tuple)):
        return [_clean(v) for v in d]
    if isinstance(d, np.generic):
        return d.item()
    return d


class EmissionReceiver(Decoder):
    capture_level = True
    key = ""
    bands: set = set()
    method = ""
    steps: list = []

    def can_decode(self, sig, ctx):
        lo, hi = ctx.cap.fc - ctx.cap.fs / 2, ctx.cap.fc + ctx.cap.fs / 2
        return bool(self.bands & set(band_keys(lo, hi))) and bool(self.bands & set(band_keys(sig.f0, sig.f1)))

    owns_band = False                # True: everything inside the signal's band belongs to it

    def in_range(self, lo, hi):
        return bool(self.bands & set(band_keys(lo, hi)))

    def run_detector(self, ctx) -> list:
        """The receiver's capture-wide detect(): -> list[core.emission.Emission]."""
        raise NotImplementedError

    def emissions(self, ctx):
        lo, hi = ctx.cap.fc - ctx.cap.fs / 2, ctx.cap.fc + ctx.cap.fs / 2
        if not self.in_range(lo, hi):
            return []
        if self.key not in ctx.cache:                   # one capture-wide run, shared by all hooks
            ctx.cache[self.key] = self.run_detector(ctx)
        return ctx.cache[self.key]

    def keep(self, e):
        return e.confidence >= 0.5 or bool((e.decoded or {}).get("crc_ok"))

    def group_key(self, e):
        return (e.label.split("_sf")[0] if "lora" in e.label else e.label, round(e.center_mhz, 1))

    def capture_signals(self, ctx):
        groups = {}
        for e in self.emissions(ctx):
            if not self.keep(e):
                continue                     # weak energy-only candidates stay with the energy detector
            groups.setdefault(self.group_key(e), []).append(e)
        out = []
        for g in groups.values():
            res = self.summarise(g, None)
            f0, f1 = res["span"]
            snr = [e.snr_db for e in g if e.snr_db is not None]
            s = receiver_signal([(e.start_us * 1e-6, e.end_us * 1e-6) for e in g], f0, f1, res,
                                float(np.median(snr)) if snr else None, self.name)
            s.features["_emissions"] = g     # for emitter association (not serialised)
            if self.owns_band:
                s.features["owns_band"] = True
            self.finish(s, g, ctx)
            out.append(s)
        return out

    def finish(self, sig, em, ctx):
        pass

    def decode(self, sig, ctx):
        em = [e for e in self.emissions(ctx) if overlaps(sig, *_span(e))]
        return self.summarise(em, sig) if em else result(
            "not_possible", self.method, self.steps,
            reason=f"{self.protocol} receiver found no matching burst in this signal's span")

    def summarise(self, em, sig):
        labels = Counter(e.label for e in em)
        dec = [e.decoded for e in em if e.decoded]
        ok = [d for d in dec if d.get("crc_ok") or d.get("crc") or d.get("serial")]
        top = labels.most_common(1)[0][0]
        fields = {"n_bursts": len(em), "labels": dict(labels),
                  "bursts": [_clean({"t_start_s": e.start_us * 1e-6, "t_end_s": e.end_us * 1e-6,
                                     "f_center_hz": e.center_mhz * 1e6, "label": e.label, "snr_db": e.snr_db,
                                     "features": {k: v for k, v in e.features.items()
                                                  if not isinstance(v, (list, dict)) or k in ("sync_syms",)},
                                     "decoded": e.decoded}) for e in em[:100]]}
        status = "success" if ok else "partial"
        conf = max(e.confidence for e in em)
        return result(status, self.method, self.steps, fields,
                      "" if ok else "physical layer identified/measured; no CRC-verified payload",
                      label=self.pretty(top, em), confidence=min(0.97, conf + (0.1 if ok else 0)),
                      modulation=self.modulation(em),
                      evidence=f"{len(em)} burst(s) matched by {self.protocol} receiver: {dict(labels)}" +
                               (f"; {len(ok)} CRC-verified payload(s)" if ok else ""),
                      span=(min(e.f_lo_mhz for e in em) * 1e6, max(e.f_hi_mhz for e in em) * 1e6))

    def pretty(self, top, em):
        return top

    def modulation(self, em):
        return em[0].features.get("modulation", "")
