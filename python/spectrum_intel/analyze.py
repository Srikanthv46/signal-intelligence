"""One snapshot end to end: detect -> features -> decode -> classify -> report."""
from __future__ import annotations

import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

from .classify import rules
from .decode import registry
from .decode.base import DecodeContext, Decoder
from .core import spurs
from .core.bands_india import BANDS
from .detect import artifacts
from .detect.merge import merge_pulse_trains
from .detect.pipeline import DetectResult, detect
from .core import parallel
from .core.dsp import Channelizer
from .features import measure
from .io.sigmf import Capture, load
from .model import ARTIFACT_FLAGS, Signal
from .sweep import SweepContext


@dataclass
class Analysis:
    cap: Capture
    det: DetectResult
    signals: list[Signal]
    quality: dict
    label_check: dict | None                            # centre-frequency label check by a plugin (BLE adverts)
    timings: dict = field(default_factory=dict)
    extras: dict = field(default_factory=dict)          # report-level additions from plugins (e.g. emitters)
    plugins: list = field(default_factory=list)         # the plugins that ran (finalize_report hooks)


def _absorb(receivers: list[Signal], energy: list[Signal], cap_s: float, pad: float = 60e-6):
    """Assign energy detections to the receiver-defined signals that explain them.

    burst energy:      >= 50% of its band inside a receiver signal's span and >= 30%
                       of its on-time covered by that receiver's bursts; or the energy box
                       contains the receiver's whole band (centre inside it) and >= 60 % of
                       its on-time is covered: spectral regrowth / shoulders of that emitter
    continuous energy: in a CSMA band the 20th-percentile level is lifted by
                       aggregate traffic; it is explained when receiver spans cover
                       >= 60% of its band and those receivers are active >= 5% of
                       the capture. It is recorded as channel occupancy."""
    taken: dict[int, list[Signal]] = {}
    for e in energy:
        # an energy box cut at the capture edge is still an emission, but only a receiver that owns its band
        # (everything inside belongs to it) explains it; artifact boxes are never absorbed
        edge_cut = e.flags == ["truncated_at_capture_edge"]
        if e.flags and not edge_cut:
            continue
        best, best_score = None, 0.0
        cont = any(b.kind == "continuous" for b in e.bursts)
        cov_f = 0.0
        for k, r in enumerate(receivers):
            if edge_cut and not r.features.get("owns_band"):
                continue
            fov = max(0.0, min(e.f1, r.f1) - max(e.f0, r.f0)) / max(e.bw, 1.0)
            if fov <= 0:
                continue
            on_r = sum(b.t1 - b.t0 for b in r.bursts) / cap_s
            if r.features.get("owns_band") and fov >= 0.8:
                score = fov                    # licensed carrier: everything inside its band is it
            elif cont:
                cov_f += fov if on_r >= 0.05 else 0.0
                score = fov * (on_r >= 0.05)
            else:
                on = sum(b.t1 - b.t0 for b in e.bursts)
                cov = sum(max(0.0, min(b.t1, x.t1 + pad) - max(b.t0, x.t0 - pad))
                          for b in e.bursts for x in r.bursts)
                score = fov if (fov >= 0.5 and on > 0 and cov / on >= 0.3) else 0.0
                if not score and on > 0 and cov / on >= 0.6 and e.f0 <= r.f0 and e.f1 >= r.f1 \
                        and r.f0 <= e.fc <= r.f1:
                    score = 0.5                # receiver band inside a wider, time-aligned box
            if score > best_score:
                best, best_score = k, score
        owned = best is not None and receivers[best].features.get("owns_band")
        if best is not None and (owned or not cont or min(cov_f, 1.0) >= 0.6):
            taken.setdefault(best, []).append(e)
    _absorb_coincident(receivers, energy, taken)
    return taken


COINC_T_S = 20e-6          # pieces of one burst start and end together
COINC_MIN_N = 3            # a wideband burst split at fading notches leaves >= 3 pieces
COINC_MIN_SPAN_HZ = 8e6    # spanning >= 8 MHz together (a narrow hop / RC packet is one piece)
COINC_FOV = 0.8            # each piece >= 80 % inside the receiver signal's span


def _absorb_coincident(receivers: list[Signal], energy: list[Signal], taken: dict[int, list[Signal]]) -> None:
    """One wideband burst (an 802.11 PPDU the receiver could not decode, an OFDM burst) that the energy detector
    split at its fading notches: >= COINC_MIN_N burst-only pieces inside one receiver signal's span whose bursts
    start and end together, spanning >= COINC_MIN_SPAN_HZ, belong to that receiver's emitter."""
    used = {id(e) for v in taken.values() for e in v}

    def match(a, b):                    # every burst of a has a burst of b with the same start and end
        return all(any(abs(x.t0 - y.t0) <= COINC_T_S and abs(x.t1 - y.t1) <= COINC_T_S for y in b.bursts)
                   for x in a.bursts)
    for k, r in enumerate(receivers):
        inside = [e for e in energy if id(e) not in used and not e.flags and e.bursts
                  and not any(b.kind == "continuous" for b in e.bursts)
                  and max(0.0, min(e.f1, r.f1) - max(e.f0, r.f0)) >= COINC_FOV * max(e.bw, 1.0)]
        inside.sort(key=lambda e: -len(e.bursts))
        groups: list[list[Signal]] = []
        for e in inside:
            g = next((g for g in groups if match(e, g[0])), None)
            if g is None:
                groups.append([e])
            else:
                g.append(e)
        for g in groups:
            if len(g) >= COINC_MIN_N and max(e.f1 for e in g) - min(e.f0 for e in g) >= COINC_MIN_SPAN_HZ:
                taken.setdefault(k, []).extend(g)
                used.update(id(e) for e in g)


def analyze_iq(iq, fs: float, fc: float, t0=None, name: str = "live", plugins: list[Decoder] | None = None,
               gain_db: float | None = None, bandwidth_hz: float | None = None, ppm: dict | None = None) -> Analysis:
    """Analyse an in-memory block (e.g. a live receiver buffer); see analyze(). For the
    production entry point (validation, versioned JSON response) see spectrum_intel.api."""
    return analyze(Capture.from_iq(iq, fs, fc, t0=t0, name=name, gain_db=gain_db, bandwidth_hz=bandwidth_hz),
                   None, plugins, ppm=ppm)


def analyze(path: str | Path | Capture, sweep: SweepContext | None = None,
            plugins: list[Decoder] | None = None, ppm: dict | None = None) -> Analysis:
    """Analyse one snapshot (a SigMF path or a Capture) with the given plugins (default: every receiver found by
    decode.registry). If a plugin's decodes prove the capture's centre-frequency label
    wrong (e.g. CRC-verified BLE adverts, > 0.5 MHz off), analyse again at the corrected
    frequency and record both in label_check.

    ppm: the receiver's frequency error, {"ppm": float, "source": str}, when the caller
    knows it (calibration, GPS-disciplined clock); otherwise the sweep's calibration, if
    any, and otherwise receivers estimate it from the capture where they can (LTE PSS)."""
    plugins = registry.plugins() if plugins is None else plugins
    a = _analyze(path, sweep, plugins, ppm=ppm)
    chk = _label_check(a, plugins)
    if chk and chk["suspect"]:
        b = _analyze(path, sweep, plugins, fc_override=a.cap.fc - chk["offset_mhz"] * 1e6, ppm=ppm)
        b.label_check = {**chk, "status": f"corrected_by_{chk['plugin']}",
                         "fc_metadata_hz": float(a.cap.meta["captures"][0]["core:frequency"]),
                         "fc_used_hz": b.cap.fc}
        b.timings["first_pass_s"] = round(sum(v for v in a.timings.values() if isinstance(v, float)), 2)
        return b
    if chk:
        a.label_check = {**chk, "status": "consistent"}
    return a


def _label_check(a: Analysis, plugins: list[Decoder]) -> dict | None:
    for p in plugins:
        try:
            chk = p.frequency_label_check(a.signals)
        except Exception:                                   # noqa: BLE001
            chk = None
        if chk:
            return {**chk, "plugin": p.name}
    return None


def _analyze(path: str | Path | Capture, sweep: SweepContext | None, plugins: list[Decoder],
             fc_override: float | None = None, ppm: dict | None = None) -> Analysis:
    T = {}
    t = time.time()
    if isinstance(path, Capture):
        from dataclasses import replace
        cap = replace(path)                                 # never retune the caller's object
    else:
        cap = load(path)
    if fc_override:
        cap.fc = fc_override
    quality = artifacts.clipping(cap.iq)
    det = detect(cap, sweep.noise_model if sweep else None)
    T["detect_s"] = round(time.time() - t, 2)

    t = time.time()
    ch = Channelizer(cap.iq, cap.fs, cap.fc)
    energy = det.signals
    for s in energy:
        measure.extract(s, ch, cap.n, cap.duration, det.fine.df)
    energy, merged = merge_pulse_trains(energy)
    for s in merged:
        measure.extract(s, ch, cap.n, cap.duration, det.fine.df)
    T["features_s"] = round(time.time() - t, 2)

    # 1) receiver-defined signals absorb the energy signals they explain
    t = time.time()
    ctx = DecodeContext(cap, ch, det)
    cal = ppm or (sweep.ppm if sweep is not None else None)
    if cal:
        ctx.cache["sweep_ppm"] = cal                    # receiver frequency error (caller / sweep calibration)
    # lines from the receiver's own clocks (the caller's / sweep's ppm only: not a receiver's own estimate, so
    # the result does not depend on which receivers ran)
    lines = spurs.receiver_lines(ch.X, cap.fs, cap.fc, cal["ppm"] if cal else None, det.usable_hz)
    ctx.cache["receiver_lines"] = lines                 # receivers skip these (e.g. voice_fm)
    artifacts.flag_receiver_spurs(energy, lines, cap.fs, cal["ppm"] if cal else None)
    out = _capture_signals(ctx, plugins, T)
    for r in out:
        for b in r.bursts:
            b.t1 = min(b.t1, cap.duration)
    for dec in plugins:
        try:
            out = dec.arbitrate(out, ctx)
        except Exception as e:                              # noqa: BLE001
            T[f"{dec.name}_arbitrate_error"] = f"{type(e).__name__}: {e}"[:200]
    absorbed = _absorb(out, energy, cap.duration)
    for k, r in enumerate(out):
        eaten = absorbed.get(k, [])
        if eaten:
            r.detection["absorbed_energy_detections"] = [
                {"f_low_hz": round(e.f0), "f_high_hz": round(e.f1), "t_start_s": round(e.t0, 6),
                 "t_end_s": round(e.t1, 6), "snr_db": e.snr_db,
                 "kind": "channel occupancy (aggregate traffic)" if any(b.kind == "continuous" for b in e.bursts)
                 else "burst"} for e in eaten]
            r.detection["description"] += (f"; explains {len(eaten)} energy detection(s) "
                                           f"(STFT/CFAR) in the same band and time")
            r.snr_db = r.snr_db or max(e.snr_db for e in eaten)
            for key in ("envelope_cv", "cyclic_line_hz", "ofdm_scs_hz"):
                v = next((e.features.get(key) for e in eaten if e.features.get(key) is not None), None)
                if v is not None:
                    r.features.setdefault(key, v)
            energy = [e for e in energy if e not in eaten]
        r.features.update({key: v for key, v in measure.timing(r, cap.duration).items()
                           if key not in r.features or key in ("n_bursts", "burst_period_s")})
    T["receivers_s"] = round(time.time() - t, 2)

    t = time.time()
    extras: dict = {}
    for dec in plugins:
        try:
            extras.update(dec.postprocess(ctx, out))
        except Exception as e:                              # noqa: BLE001
            T[f"{dec.name}_postprocess_error"] = f"{type(e).__name__}: {e}"[:200]
    T["postprocess_s"] = round(time.time() - t, 2)

    # 2) per-signal decoders on what is left, then rules
    t = time.time()
    for s in energy:
        attempts = []
        if not s.flags:
            for dec in plugins:
                try:
                    if dec.capture_level or not dec.can_decode(s, ctx):
                        continue
                    r = dec.decode(s, ctx)
                except Exception as e:                      # noqa: BLE001
                    attempts.append(f"{dec.name}: error {type(e).__name__}: {e}"[:200])
                    continue
                attempts.append(f"{dec.name}: {r['status']}" + (f" ({r['reason']})" if r["reason"] else ""))
                if r["status"] in ("success", "partial"):
                    s.decoded = {k: r[k] for k in ("status", "reason", "method", "steps", "fields")}
                    if r.get("label"):
                        s.label, s.confidence = r["label"], r["confidence"]
                        s.modulation = r.get("modulation") or s.modulation
                        s.evidence = r.get("evidence") or ""
                    if r.get("span"):
                        s.features["decoder_span_hz"] = [round(v) for v in r["span"]]
                    break
        if s.decoded["status"] == "not_possible" or not s.decoded.get("method"):
            rules.classify(s)
            s.decoded = {"status": "not_possible", "reason": _why_not(s, attempts), "method": None,
                         "steps": attempts, "fields": {}}
        else:
            s.decoded["steps"] = attempts + s.decoded["steps"]
        out.append(s)
    T["decode_classify_s"] = round(time.time() - t, 2)

    for s in out:
        _image_caution(s, cap.fc)
        if quality["clipping"]:
            s.flags.append("adc_clipping_in_capture")
    out.sort(key=lambda s: (s.f0, s.t0))
    return Analysis(cap, det, out, quality, None, T, extras, plugins)


def _image_caution(s: Signal, lo: float) -> None:
    """Single-block LO-harmonic hint (core.spurs): an unidentified signal outside the cellular bands that a 3rd or
    5th LO-harmonic image of a cellular downlink would explain gets a caution. A hint only: label and confidence
    stay (one block cannot prove it; the sweep's harmonic check and a second tuning can)."""
    if s.decoded["status"] != "not_possible" or set(ARTIFACT_FLAGS) & set(s.flags):
        return
    if spurs.in_strong_band(s.f0, s.f1, BANDS):
        return
    src = spurs.harmonic_image_sources(s.f0, s.f1, lo, BANDS)
    if not src:
        return
    s.flags.append("possible_lo_harmonic_image")
    s.features["lo_harmonic_image_sources"] = src
    s.evidence = ((s.evidence + "; ") if s.evidence else "") + "CAUTION possible LO-harmonic image of " + ", ".join(
        f"{c['source_low_hz'] / 1e6:.1f}-{c['source_high_hz'] / 1e6:.1f} MHz ({c['band']}, {c['n']} x LO"
        f"{', inverted' if c['inverted'] else ''})" for c in src) + " (no preselector: confirm with a second tuning)"


def _package(dec: Decoder) -> str:
    m = type(dec).__module__.split(".")
    return m[m.index("receivers") + 1] if "receivers" in m else dec.name


def _capture_signals(ctx: DecodeContext, plugins: list[Decoder], T: dict) -> list[Signal]:
    """Every plugin's capture_signals(), grouped by receiver package. Plugins of one
    package run in order; packages run concurrently, each on its own copy of the
    capture cache, which sees only what the packages named in a plugin's `after`
    published (e.g. nr after lte for the receiver ppm). That is exactly what the
    sequential registry order gave, so results do not depend on scheduling. The
    packages' caches are then merged into ctx.cache in registry order."""
    groups: dict[str, list[Decoder]] = {}
    for dec in plugins:
        groups.setdefault(_package(dec), []).append(dec)
    names = list(groups)
    base = dict(ctx.cache)
    caches: dict[str, dict] = {}
    results: dict[str, list[Signal]] = {}
    timings: dict[str, float] = {}
    errors: dict[str, str] = {}
    futures: dict = {}

    def run(name: str) -> None:
        deps = [d for dec in groups[name] for d in getattr(dec, "after", ()) if d in futures and d != name]
        cache = dict(base)
        for d in deps:
            futures[d].result()
            cache.update({k: v for k, v in caches[d].items() if k not in base})
        gctx = DecodeContext(ctx.cap, ctx.ch, ctx.det, cache)
        got: list[Signal] = []
        for dec in groups[name]:
            t_dec = time.time()
            try:
                got += dec.capture_signals(gctx)
            except Exception as e:                          # noqa: BLE001 — a broken plugin must not kill the run
                errors[dec.name] = f"{type(e).__name__}: {e}"[:200]
            timings[dec.name] = round(time.time() - t_dec, 2)
        caches[name], results[name] = cache, got

    n = min(parallel.threads(), len(names))
    if n <= 1:
        for name in names:
            futures[name] = _Done()
            run(name)
    else:
        # a package waits only on packages earlier in registry order; tasks start in
        # submission (FIFO) order, so whatever it waits on is already running or done
        with ThreadPoolExecutor(max_workers=n) as ex:
            for name in names:
                futures[name] = ex.submit(run, name)
            for f in futures.values():
                f.result()
    out: list[Signal] = []
    for name in names:
        out += results.get(name, [])
        ctx.cache.update({k: v for k, v in caches.get(name, {}).items() if k not in base})
    for dec in plugins:
        if dec.name in errors:
            T[f"{dec.name}_error"] = errors[dec.name]
        T.setdefault("receiver_s", {})[dec.name] = timings.get(dec.name, 0.0)
    return out


class _Done:
    @staticmethod
    def result():
        return None


def _why_not(s: Signal, attempts: list[str]) -> str:
    if "lo_leakage" in s.flags or "iq_image" in s.flags:
        return "receiver artifact, not an emission"
    lab = s.label
    if lab.startswith("CW"):
        return "unmodulated carrier: no information content"
    if s.snr_db < 6:
        return f"insufficient SNR ({s.snr_db} dB) for demodulation"
    if "truncated_at_capture_edge" in s.flags:
        return "signal extends beyond the captured bandwidth (truncated)"
    known = {
        "LTE": "no LTE/NR decoder matched (PSS/SSS not found in this span)",
        "GSM": "GSM FCCH/SCH not present in this 32.5 ms window",
        "UMTS": "WCDMA decoding (scrambling-code search, P-CCPCH) not implemented",
        "DVB-T2": "DVB-T2 L1 decoding not implemented; 32 ms covers < 1 T2 frame",
        "PMR": "32.5 ms is too short for voice; digital PMR (DMR/TETRA) vocoder/encryption not handled",
        "TETRA": "TETRA air interface is typically encrypted (TEA); 32.5 ms covers < 1 multiframe",
        "radar": "pulsed radar carries no decodable payload; PRI/pulse width reported",
        "DME": "DME pulse pairs carry no payload; X/Y channel from pair spacing only",
        "Bluetooth": "classic BT / non-advertising BLE uses hopping + unknown access address",
        "ISM": "proprietary SRD protocol (OOK/FSK) — no decoder for the specific device",
    }
    for k, v in known.items():
        if k in lab:
            return v
    if lab == "unknown":
        return "unknown protocol: no decoder applicable"
    return "no decoder implemented for this signal class" + (f"; tried: {', '.join(attempts)}" if attempts else "")
