"""Plugin interface between the pipeline and the receivers.

A receiver package (spectrum_intel/receivers/<name>/) keeps its DSP in modules that
import only spectrum_intel.core, plus a plugin.py that adapts it to the pipeline:

    class MyDecoder(Decoder):
        name = "my_protocol"
        capture_level = True                       # defines its own signals
        def capture_signals(self, ctx) -> list[Signal]: ...
    PLUGINS = [MyDecoder()]

or, for a decoder of energy-detected signals,

        def can_decode(self, sig, ctx) -> bool: ...
        def decode(self, sig, ctx) -> dict: ...    # see result() for the shape

Every hook is optional (the defaults do nothing):

    capture_signals(ctx)              receiver-defined signals; energy detections they explain are absorbed
    arbitrate(signals, ctx)           after all capture_signals: drop/merge other receivers' signals this
                                      plugin explains (e.g. analog VTX sidebands); returns the new list
    frequency_label_check(signals)    {"offset_mhz", "suspect", ...} if the plugin's decodes prove the
                                      capture's centre-frequency label wrong (the pipeline re-analyses)
    postprocess(ctx, signals)         report-level extras, e.g. {"emitters": [...]}
    finalize_report(rep)              edit the report dict after signal ids are final

A plugin may cache capture-wide work in ctx.cache. capture_signals() of different
receiver packages run concurrently (analyze._capture_signals), each package on its own
copy of the cache: a plugin that reads what another package's capture_signals() put in
the cache names that package in `after` (e.g. nr: after = ("lte",) for the receiver
ppm). Every other hook runs sequentially in registry order (decode.registry); for
energy signals the first success/partial decode wins, but every attempt is recorded
in decoded.steps.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from ..core.dsp import Channelizer
from ..io.sigmf import Capture
from ..model import Signal


@dataclass
class DecodeContext:
    cap: Capture
    ch: Channelizer
    det: object                     # detect.pipeline.DetectResult
    cache: dict = field(default_factory=dict)


def result(status: str, method: str, steps: list[str], fields: dict | None = None, reason: str = "",
           label: str | None = None, confidence: float | None = None, modulation: str | None = None,
           evidence: str | None = None, span: tuple[float, float] | None = None) -> dict:
    """status: success | partial | not_possible. label/confidence/modulation/
    evidence override the classifier when the decoder knows better; span
    (f_low, f_high) refines the occupied band."""
    return {"status": status, "method": method, "steps": steps, "fields": fields or {}, "reason": reason,
            "label": label, "confidence": confidence, "modulation": modulation, "evidence": evidence,
            "span": span}


class Decoder:
    name = "base"
    protocol = ""
    capture_level = False          # True: defines its own signals via capture_signals()
    after: tuple = ()              # receiver packages whose capture_signals() cache entries this one reads

    def can_decode(self, sig: Signal, ctx: DecodeContext) -> bool:
        return False

    def decode(self, sig: Signal, ctx: DecodeContext) -> dict:
        raise NotImplementedError

    def capture_signals(self, ctx: DecodeContext) -> list[Signal]:
        """Optional: signals defined by this receiver's own decodes (one per
        channel / carrier), built from a capture-wide run. Energy signals they
        explain are absorbed by the pipeline."""
        return []

    def arbitrate(self, signals: list[Signal], ctx: DecodeContext) -> list[Signal]:
        return signals

    def frequency_label_check(self, signals: list[Signal]) -> dict | None:
        return None

    def postprocess(self, ctx: DecodeContext, signals: list[Signal]) -> dict:
        return {}

    def finalize_report(self, rep: dict) -> None:
        pass


Plugin = Decoder                   # a plugin need not decode anything (e.g. emitter association)


def receiver_signal(bursts: list[tuple[float, float]], f0: float, f1: float, res: dict, snr_db: float | None,
                    decoder: str) -> Signal:
    from ..model import Box
    bx = [Box(t0, t1, f0, f1, snr_db or 0.0, snr_db or 0.0, "receiver") for t0, t1 in sorted(bursts)]
    s = Signal(bursts=bx, f0=f0, f1=f1, snr_db=round(snr_db, 1) if snr_db is not None else 0.0)
    s.decoded = {k: res[k] for k in ("status", "reason", "method", "steps", "fields")}
    s.label, s.confidence = res["label"], res["confidence"]
    s.modulation, s.evidence = res["modulation"] or "", res["evidence"] or ""
    s.detection = {"method": f"protocol receiver ({decoder})",
                   "description": f"defined by {len(bx)} burst(s) the {decoder} receiver synchronised to and "
                                  f"decoded (preamble/sync correlation + header checks), not by energy alone",
                   "n_bursts": len(bx)}
    s.features["decoder"] = decoder
    return s


def overlaps(sig: Signal, t0: float, t1: float, f0: float, f1: float, f_frac: float = 0.3) -> bool:
    fov = max(0.0, min(sig.f1, f1) - max(sig.f0, f0))
    if fov < f_frac * min(sig.bw, f1 - f0):
        return False
    return any(min(b.t1, t1) > max(b.t0, t0) - 50e-6 for b in sig.bursts)
