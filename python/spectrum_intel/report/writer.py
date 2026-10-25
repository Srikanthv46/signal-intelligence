"""Per-snapshot JSON report (user schema + provenance blocks)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ..core.bands_india import describe


def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, float) and not np.isfinite(o):
        return None
    if isinstance(o, bytes):
        return o.hex()
    return o


def signal_dict(i: int, s) -> dict:
    f = s.features
    bursts = sorted(s.bursts, key=lambda b: b.t0)
    feats = {"symbol_rate": f.get("symbol_rate"), "burst_period_s": f.get("burst_period_s"),
             "duty_cycle": f.get("duty_cycle")}
    feats.update({k: v for k, v in f.items()
                  if k not in feats and k != "classifier_candidates" and not k.startswith("_")})
    feats["bursts"] = [{"t_start_s": round(b.t0, 7), "t_end_s": round(b.t1, 7),
                        "f_low_hz": round(b.f0), "f_high_hz": round(b.f1), "snr_db": b.mean_snr_db}
                       for b in bursts[:300]]
    det = dict(s.detection)
    if f.get("classifier_candidates") is not None:
        det["classifier_candidates"] = f["classifier_candidates"]
    return {
        "id": i,
        "label": s.label,
        "description": _description(s.label, feats.get("decoder"), s.flags),
        "confidence": s.confidence,
        "t_start_s": round(s.t0, 7), "t_end_s": round(s.t1, 7),
        "f_low_hz": round(s.f0), "f_high_hz": round(s.f1), "f_center_hz": round(s.fc), "bandwidth_hz": round(s.bw),
        "snr_db": s.snr_db,
        "modulation": s.modulation,
        "features": feats,
        "evidence": s.evidence,
        "band_allocation_india": describe(s.fc),
        "flags": s.flags,
        "detection": det,
        "decoded": {"status": s.decoded["status"], "reason": s.decoded.get("reason", ""),
                    "method": s.decoded.get("method"), "steps": s.decoded.get("steps", []),
                    "fields": s.decoded.get("fields", {}),
                    "validation": _validation(s.label, feats.get("decoder"), s.decoded.get("fields", {}))},
    }


def _description(label, decoder, flags):
    """What the signal is, in plain words (schema/descriptions.json via api.contract)."""
    from ..api import contract as K
    return K.description(label, decoder, flags)


def _validation(label, decoder, fields):
    """proven / research tag (schema/capabilities.json via api.contract)."""
    from ..api import contract as K
    return K.validation(K.technology(label, decoder)[0], decoder, fields)


WEAK_SNR_DB = 6.0
WEAK_NOTE = ("energy-only burst detections below 6 dB in-band SNR with no decoder or rule match: at this "
             "level the CFAR false-alarm rate is set by noise-floor estimate error (~1 dB), so they are "
             "listed for completeness but not reported as signals")


def is_weak(s: dict) -> bool:
    if s["decoded"]["status"] != "not_possible":
        return False
    if "continuous" in s["detection"].get("kinds", []):
        return False
    if {"lo_leakage", "receiver_spur"} & set(s["flags"]):
        return False
    return (s["snr_db"] or 0) < WEAK_SNR_DB


def split_weak(rep: dict) -> dict:
    keep, weak = [], list(rep.get("weak_detections", []))
    for s in rep["signals"]:
        if is_weak(s):
            weak.append({"t_start_s": s["t_start_s"], "t_end_s": s["t_end_s"], "f_low_hz": s["f_low_hz"],
                         "f_high_hz": s["f_high_hz"], "snr_db": s["snr_db"], "rule_label": s["label"],
                         "flags": s["flags"]})
        else:
            keep.append(s)
    for i, s in enumerate(keep):
        s["id"] = i
    rep["signals"] = keep
    rep["weak_detections"] = weak
    rep["weak_detections_note"] = WEAK_NOTE
    return rep


def build(a) -> dict:
    cap = a.cap
    rep = {
        "file": cap.name,
        "source_path": str(cap.path.resolve()) if cap.path.exists() else None,
        "center_freq_hz": cap.fc,
        "sample_rate_hz": cap.fs,
        "capture_time_utc": cap.t0.isoformat() if cap.t0 else "",
        "frequency_label_check": a.label_check or {"status": "not_checked",
                                                   "note": "no plugin evidence (e.g. CRC-valid BLE adverts)"},
        "duration_s": round(cap.duration, 6),
        "gain_db": cap.gain_db,
        "capture_quality": {**a.quality,
                            "noise_floor": a.det.info.get("noise"),
                            "usable_band_hz": [round(v) for v in a.det.usable_hz]},
        "processing": {"timings_s": a.timings,
                       "pipeline": ["detect: coarse+fine STFT, noise floor (sweep template / opening), "
                                    "per-frame gain g(t), continuous + burst CFAR, merge, valley split, "
                                    "group bursts by frequency span, LO-leak / IQ-image flags",
                                    "features: timing, envelope, cyclic lines, OFDM CP autocorrelation, FSK",
                                    "receivers (plugins, see processing.plugins): capture-wide receiver runs; "
                                    "receiver-defined signals absorb matching energy detections; plugin "
                                    "arbitration, post-processing (e.g. emitter association)",
                                    "per-signal decoders on remaining signals, then India band rules"]},
        "signals": [signal_dict(i, s) for i, s in enumerate(a.signals)],
        **a.extras,
    }
    rep["processing"]["plugins"] = [p.name for p in a.plugins]
    rep = split_weak(_clean(rep))
    for p in a.plugins:
        p.finalize_report(rep)
    return rep


def write(a, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / f"{a.cap.name}_report.json"
    p.write_text(json.dumps(build(a), indent=1))
    return p
