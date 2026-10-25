"""Analysis -> response dict (schema/response.schema.json).

The response is built from the per-snapshot report (report.writer.build), so the two
always agree; it re-shapes it for machine consumers: stable technology codes and
families, grouped frequency/time blocks, identifiers, coverage and warnings.
"""
from __future__ import annotations

from collections import Counter
from datetime import datetime, timedelta, timezone

from . import contract as K

_TIMING_FEATURES = {"bursts", "duty_cycle", "burst_period_s", "n_bursts", "continuous", "decoder"}


def _utc(t0: datetime | None, s: float | None) -> str | None:
    if t0 is None or s is None:
        return None
    t = t0 if t0.tzinfo else t0.replace(tzinfo=timezone.utc)
    return (t.astimezone(timezone.utc) + timedelta(seconds=s)).isoformat().replace("+00:00", "Z")


def _frequency(lo: float, hi: float, fc: float) -> dict:
    c = round(0.5 * (lo + hi))
    return {"center_hz": c, "low_hz": lo, "high_hz": hi, "bandwidth_hz": round(hi - lo), "offset_hz": round(c - fc)}


def _time(t0s: float, t1s: float, start: datetime | None) -> dict:
    return {"start_s": t0s, "end_s": t1s, "duration_s": round(t1s - t0s, 7),
            "start_utc": _utc(start, t0s), "end_utc": _utc(start, t1s)}


def signal(s: dict, fc: float, start: datetime | None, detail: str = "full") -> dict:
    f = s["features"]
    tech, fam = K.technology(s["label"], f.get("decoder"))
    artifact = fam == "artifact" or any(x in s["flags"] for x in K.ARTIFACT_FLAGS)
    if artifact and fam != "artifact":
        tech, fam = next((x for x in K.ARTIFACT_FLAGS if x in s["flags"]), tech), "artifact"
    dec = s["decoded"]
    fields = dec.get("fields") or {}
    t = _time(s["t_start_s"], s["t_end_s"], start)
    t.update(n_bursts=f.get("n_bursts", len(f.get("bursts", []))), duty_cycle=f.get("duty_cycle"),
             burst_period_s=f.get("burst_period_s"), continuous=bool(f.get("continuous")))
    out = {
        "id": s["id"],
        "category": "artifact" if artifact else "emission",
        "technology": tech,
        "family": fam,
        "label": s["label"],
        "description": s.get("description") or K.description(s["label"], f.get("decoder"), s["flags"]),
        "confidence": s["confidence"],
        "inferred": dec["status"] == "not_possible" and s["label"] != "unknown",
        "frequency": _frequency(s["f_low_hz"], s["f_high_hz"], fc),
        "time": t,
        "snr_db": s["snr_db"],
        "modulation": s["modulation"],
        "decode": {"status": dec["status"], "reason": dec.get("reason") or "", "method": dec.get("method"),
                   "fields": fields},
        "identifiers": {k: fields[k] for k in K.identifier_keys() if fields.get(k) not in (None, [], {}, "")},
        "validation": K.validation(tech, f.get("decoder"), fields),
        "evidence": s["evidence"],
        "band_allocation": s["band_allocation_india"],
        "flags": s["flags"],
    }
    if detail == "full":
        out["decode"]["steps"] = dec.get("steps", [])
        out["time"]["bursts"] = [{"start_s": b["t_start_s"], "end_s": b["t_end_s"], "low_hz": b["f_low_hz"],
                                  "high_hz": b["f_high_hz"], "snr_db": b["snr_db"]} for b in f.get("bursts", [])]
        out["detection"] = s["detection"]
        out["features"] = {k: v for k, v in f.items() if k not in _TIMING_FEATURES}
    return out


def weak(w: dict, fc: float, start: datetime | None) -> dict:
    return {"frequency": _frequency(w["f_low_hz"], w["f_high_hz"], fc),
            "time": _time(w["t_start_s"], w["t_end_s"], start),
            "snr_db": w["snr_db"], "rule_label": w["rule_label"], "flags": w["flags"]}


def coverage(packages: dict[str, list[str]], fs: float, bw: float | None, band: tuple[float, float],
             dur: float, timings: dict) -> tuple[list[dict], list[dict]]:
    """One entry per enabled receiver package, plus warnings for in-band receivers the
    block cannot serve (too short, sample rate too low)."""
    req = K.requirements()["receivers"]
    rate = min(fs, bw) if bw else fs
    cov, warn = [], []
    for pkg, plugin_names in packages.items():
        r = req.get(pkg, {"min_sample_rate_hz": 0, "min_duration_s": 0, "recommended_duration_s": 0,
                          "bands_hz": [], "why": ""})
        bands = r.get("bands_hz") or []
        in_band = not bands or any(lo < band[1] and hi > band[0] for lo, hi in bands)
        rate_ok = rate >= r["min_sample_rate_hz"]
        dur_ok = dur >= r["min_duration_s"]
        full = dur >= r["recommended_duration_s"]
        err = [timings[f"{n}_error"] for n in plugin_names if f"{n}_error" in timings]
        secs = [timings.get("receiver_s", {}).get(n) for n in plugin_names]
        note = "" if (rate_ok and full) else r.get("why", "")
        cov.append({"receiver": pkg, "status": "error" if err else "ran", "in_band": in_band,
                    "sample_rate_ok": rate_ok, "duration_ok": dur_ok, "full_duration": full,
                    "min_sample_rate_hz": r["min_sample_rate_hz"], "min_duration_s": r["min_duration_s"],
                    "recommended_duration_s": r["recommended_duration_s"],
                    "time_s": round(sum(x for x in secs if x is not None), 3) if any(x is not None for x in secs)
                    else None,
                    "note": note})
        if in_band and not rate_ok:
            warn.append({"code": "SAMPLE_RATE_TOO_LOW", "receiver": pkg,
                         "message": f"{pkg} needs >= {r['min_sample_rate_hz'] / 1e6:g} MHz of bandwidth; "
                                    f"this block has {rate / 1e6:g} MHz"})
        if in_band and not dur_ok:
            warn.append({"code": "BLOCK_TOO_SHORT_FOR_RECEIVER", "receiver": pkg,
                         "message": f"{pkg} needs >= {r['min_duration_s'] * 1e3:g} ms; this block is "
                                    f"{dur * 1e3:.3g} ms ({r.get('why', '')})"})
    return cov, warn


def build(rep: dict, *, analysis, block, engine: dict, packages: dict[str, list[str]], total_s: float,
          detail: str = "full") -> dict:
    """The response for one analysed block. rep = report.writer.build(analysis)."""
    cap = analysis.cap
    fc_used = cap.fc
    start = block.start_time
    q = rep["capture_quality"]
    band = tuple(q["usable_band_hz"])
    dur = rep["duration_s"]
    T = rep["processing"]["timings_s"]
    signals = [signal(s, fc_used, start, detail) for s in rep["signals"]]
    weaks = [weak(w, fc_used, start) for w in rep.get("weak_detections", [])]
    cov, warnings = coverage(packages, cap.fs, block.bandwidth_hz, band, dur, T)
    req = K.requirements()["pipeline"]
    if cap.n < req["recommended_samples"] and dur < req["recommended_duration_s"]:
        warnings.insert(0, {"code": "SHORT_BLOCK", "receiver": None,
                            "message": f"{cap.n} samples ({dur * 1e3:.3g} ms): below the recommended "
                                       f"{req['recommended_samples']} samples / {req['recommended_duration_s'] * 1e3:g}"
                                       f" ms; expect fewer decodes and more spurious energy detections"})
    if q.get("clipping"):
        warnings.append({"code": "ADC_CLIPPING", "receiver": None,
                         "message": f"{q['clipped_fraction'] * 100:.3g} % of samples at ADC full scale: "
                                    f"intermodulation products may appear as signals; lower the gain"})
    lc = rep.get("frequency_label_check") or {}
    if str(lc.get("status", "")).startswith("corrected"):
        warnings.append({"code": "CENTER_FREQUENCY_CORRECTED", "receiver": lc.get("plugin"),
                         "message": f"decoded content proves the block's centre frequency wrong by "
                                    f"{lc.get('offset_mhz')} MHz; analysed at {fc_used / 1e6:.6f} MHz"})
    errors = [{"code": "RECEIVER_ERROR", "receiver": k[:-len("_error")], "message": v}
              for k, v in T.items() if k.endswith("_error") and isinstance(v, str)]
    fam = Counter(s["family"] for s in signals)
    consumed = {"file", "source_path", "center_freq_hz", "sample_rate_hz", "capture_time_utc", "frequency_label_check",
                "duration_s", "gain_db", "capture_quality", "processing", "signals", "weak_detections",
                "weak_detections_note", "emitters", "devices", "findings"}
    stages = {k[:-2]: round(v * 1e3, 1) for k, v in T.items() if k.endswith("_s") and isinstance(v, (int, float))}
    return {
        "schema_version": K.SCHEMA_VERSION,
        "kind": "spectrum_intel.analysis",
        "request_id": block.request_id,
        "status": "ok",
        "engine": engine,
        "input": {"source_id": block.source_id,
                  "center_frequency_hz": block.center_frequency_hz,
                  "center_frequency_used_hz": fc_used,
                  "sample_rate_hz": cap.fs,
                  "bandwidth_hz": block.bandwidth_hz,
                  "n_samples": cap.n,
                  "duration_s": dur,
                  "start_time_utc": _utc(start, 0.0),
                  "gain_db": block.gain_db,
                  "receiver_ppm": block.receiver_ppm,
                  "analysed_band_hz": list(band)},
        "quality": {"adc_clipping": bool(q.get("clipping")), "clipped_fraction": q.get("clipped_fraction"),
                    "peak_abs": q.get("peak_abs"), "full_scale": q.get("full_scale"), "rms_dbfs": q.get("rms_dbfs"),
                    "noise_floor_dbfs_per_bin_median": analysis.det.info.get("noise_floor_dbfs_per_bin_median"),
                    "frequency_label_check": lc or None},
        "coverage": cov,
        "summary": {"n_signals": len(signals), "n_weak_detections": len(weaks),
                    "n_artifacts": sum(s["category"] == "artifact" for s in signals),
                    "n_unknown": sum(s["technology"] == "unknown" for s in signals),
                    "by_status": dict(Counter(s["decode"]["status"] for s in signals)),
                    "by_family": dict(fam),
                    "technologies": sorted({s["technology"] for s in signals})},
        "signals": signals,
        "weak_detections": weaks,
        "emitters": rep.get("emitters", []),
        "devices": rep.get("devices", []),
        "findings": rep.get("findings", []),
        "extensions": {k: v for k, v in rep.items() if k not in consumed},
        "warnings": warnings,
        "errors": errors,
        "timing_ms": {"total": round(total_s * 1e3, 1), "stages": stages,
                      "receivers": {k: round(v * 1e3, 1) for k, v in T.get("receiver_s", {}).items()}},
    }


def failure(block, engine: dict, status: str, code: str, message: str, total_s: float = 0.0) -> dict:
    """A response for a block that could not be analysed (same shape, empty results)."""
    n = int(getattr(block, "n_samples", 0) or 0)
    fs = float(block.sample_rate_hz or 0) or 1.0
    fc = float(block.center_frequency_hz or 0)
    return {
        "schema_version": K.SCHEMA_VERSION, "kind": "spectrum_intel.analysis", "request_id": block.request_id,
        "status": status, "engine": engine,
        "input": {"source_id": block.source_id, "center_frequency_hz": fc, "sample_rate_hz": fs,
                  "bandwidth_hz": block.bandwidth_hz, "n_samples": n, "duration_s": n / fs,
                  "start_time_utc": _utc(block.start_time, 0.0), "gain_db": block.gain_db,
                  "receiver_ppm": getattr(block, "receiver_ppm", None),
                  "analysed_band_hz": [fc - fs / 2, fc + fs / 2]},
        "quality": {}, "coverage": [],
        "summary": {"n_signals": 0, "n_weak_detections": 0, "n_artifacts": 0, "n_unknown": 0, "by_status": {},
                    "by_family": {}, "technologies": []},
        "signals": [], "weak_detections": [], "emitters": [], "devices": [], "findings": [], "extensions": {},
        "warnings": [], "errors": [{"code": code, "receiver": None, "message": message}],
        "timing_ms": {"total": round(total_s * 1e3, 1)},
    }
