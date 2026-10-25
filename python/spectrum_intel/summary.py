"""Sweep-level pass over all per-snapshot reports.

1. Cross-tuning consistency for narrow continuous signals (tones): adjacent
   tunings overlap by ~30 MHz, so a real emitter at RF f shows up at f in every
   tuning that covers f, while a receiver spur keeps its *baseband* offset and
   moves with the LO. Each tone is checked against the other tunings that cover
   it; reports are patched and their PNGs redrawn.
2. Occupancy map 300 MHz - 6 GHz (category-coloured, strongest SNR per bin).
3. summary.json (per-category counts, decoded identities, data-quality notes)
   and signals.csv (one row per signal per snapshot).
"""
from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

CATS = [  # (category, substrings in label, colour)
    ("artifact", ["artifact", "Receiver spur", "LO leakage", "IQ image", "LO-harmonic"], "#9e9e9e"),
    ("drone", ["DJI", "OcuSync", "DroneID", "FPV", "Herelink", "video/data link"], "#ff1744"),
    ("wifi", ["Wi-Fi"], "#2979ff"),
    ("bluetooth", ["Bluetooth", "BLE"], "#00b8d4"),
    ("lte", ["LTE"], "#ff6d00"),
    ("nr", ["5G NR", "NR-like", "NR "], "#d50000"),
    ("gsm", ["GSM"], "#ffab00"),
    ("umts", ["UMTS"], "#ff9100"),
    ("broadcast", ["DVB", "TV"], "#aa00ff"),
    ("aero", ["ADS-B", "Mode S", "SSR", "DME"], "#64dd17"),
    ("radar", ["radar"], "#00c853"),
    ("pmr", ["PMR", "TETRA", "Radiosonde", "FM"], "#c51162"),
    ("ism_srd", ["ISM", "LoRa", "ExpressLRS", "fsk", "FSK"], "#ffd600"),
    ("satellite", ["Iridium", "MSS"], "#6200ea"),
    ("cw", ["CW carrier"], "#bdbdbd"),
    ("unknown", ["unknown"], "#795548"),          # brown: white was invisible on the map
]


def category(label: str) -> str:
    for c, keys, _ in CATS:
        if any(k in label for k in keys):
            return c
    return "other"


def _load(sweep_out: Path):
    reps = []
    for p in sorted(sweep_out.rglob("*_report.json")):
        r = json.loads(p.read_text())
        r["_path"] = str(p)
        reps.append(r)
    return reps


def spur_check(reps, tol=40e3, usable=0.45):
    """Patch narrow continuous signals with cross-tuning evidence. Returns
    the set of report paths that changed."""
    by_fc = defaultdict(list)
    for r in reps:
        by_fc[r["center_freq_hz"]].append(r)
    changed = set()
    for r in reps:
        fc_a, fs = r["center_freq_hz"], r["sample_rate_hz"]
        for s in r["signals"]:
            if s["bandwidth_hz"] > 60e3 or s["t_end_s"] - s["t_start_s"] < 0.8 * r["duration_s"]:
                continue
            if {"lo_leakage", "iq_image"} & set(s["flags"]) or s["decoded"]["status"] != "not_possible":
                continue
            f = s["f_center_hz"]
            off = f - fc_a
            same_rf, same_bb, covered = 0, 0, 0
            for fc_b, rs in by_fc.items():
                if fc_b == fc_a:
                    continue
                for rb in rs:
                    hits_rf = abs(f - fc_b) < usable * fs
                    f_bb = fc_b + off
                    if hits_rf:
                        covered += 1
                        if any(abs(o["f_center_hz"] - f) < tol for o in rb["signals"]):
                            same_rf += 1
                    if abs(off) < usable * fs and any(abs(o["f_center_hz"] - f_bb) < tol and o["bandwidth_hz"] < 80e3
                                                      for o in rb["signals"]):
                        same_bb += 1
            ev = {"other_tunings_covering_rf": covered, "seen_at_same_rf": same_rf,
                  "seen_at_same_baseband_offset_elsewhere": same_bb}
            s["detection"]["sweep_consistency"] = ev
            if covered and same_rf == 0 and same_bb >= 2:
                s["flags"].append("receiver_spur")
                s["label"] = "Receiver spur (moves with LO)"
                s["confidence"] = round(min(0.95, 0.6 + 0.1 * same_bb), 2)
                s["evidence"] = (f"narrow tone at baseband {off / 1e6:+.3f} MHz; absent at the same RF in "
                                 f"{covered} other tuning(s) covering it, but present at the same baseband offset "
                                 f"in {same_bb} other tuning(s): generated inside the receiver")
                s["decoded"]["reason"] = "receiver artifact, not an emission"
                changed.add(r["_path"])
            elif same_rf >= 1:
                s["evidence"] += f"; confirmed external: also seen at the same RF in {same_rf} other tuning(s)"
                if s["label"].startswith("CW carrier"):
                    s["label"] = "CW carrier (external)-like"
                changed.add(r["_path"])
    return changed


def _ovl(a0, a1, b0, b1):
    return max(0.0, min(a1, b1) - max(a0, b0))


def harmonic_check(reps, usable=0.45, harmonics=(3, 5, 7)):
    """LO-harmonic responses. The AD936x mixer also converts RF near n*LO
    (n odd) to baseband, so a strong emitter at F shows up at LO + (F - n*LO)
    (or its mirror). Such a response holds a fixed *baseband* offset relative
    to n*LO, so it is absent at the same RF in other tunings that cover it.
    A signal is flagged when (a) other tunings covering its RF do not see it
    and (b) a >= 10 dB stronger signal exists at n*LO +- its baseband offset."""
    changed = set()
    for r in reps:
        L, fs = r["center_freq_hz"], r["sample_rate_hz"]
        for s in r["signals"]:
            if s["decoded"]["status"] != "not_possible" or s["bandwidth_hz"] < 100e3:
                continue
            if {"lo_leakage", "iq_image", "receiver_spur"} & set(s["flags"]):
                continue
            f0, f1 = s["f_low_hz"], s["f_high_hz"]
            covering = [b for b in reps if b["center_freq_hz"] != L and
                        abs(f0 - b["center_freq_hz"]) < usable * fs and abs(f1 - b["center_freq_hz"]) < usable * fs]
            seen = sum(any(_ovl(f0, f1, o["f_low_hz"], o["f_high_hz"]) >= 0.5 * (f1 - f0) for o in b["signals"])
                       for b in covering)
            src = None
            for n in harmonics:
                for sign in (1, -1):
                    h = sorted((n * L + sign * (f0 - L), n * L + sign * (f1 - L)))
                    for b in reps:
                        if abs(h[0] - b["center_freq_hz"]) > 0.48 * fs or abs(h[1] - b["center_freq_hz"]) > 0.48 * fs:
                            continue
                        for o in b["signals"]:
                            if {"iq_image", "lo_leakage"} & set(o["flags"]):
                                continue
                            if _ovl(h[0], h[1], o["f_low_hz"], o["f_high_hz"]) >= 0.5 * (h[1] - h[0]) and \
                                    (o["snr_db"] or 0) >= (s["snr_db"] or 0) + 10:
                                if src is None or o["snr_db"] > src[2]["snr_db"]:
                                    src = (n, sign, o, h)
            ev = {"other_tunings_covering_rf": len(covering), "seen_at_same_rf": seen}
            if src:
                n, sign, o, h = src
                ev["harmonic_candidate"] = {"n": n, "mirror": sign < 0, "rf_hz": [round(h[0]), round(h[1])],
                                            "source_label": o["label"], "source_snr_db": o["snr_db"]}
            s["detection"]["sweep_consistency"] = {**s["detection"].get("sweep_consistency", {}), **ev}
            if src and covering and seen == 0:
                n, sign, o, h = src
                s["flags"].append("lo_harmonic_response")
                s["label"] = "LO-harmonic response (receiver artifact)"
                s["confidence"] = 0.85
                s["evidence"] = (f"absent at the same RF in {len(covering)} other tuning(s) covering it; maps to "
                                 f"{n}xLO {'-' if sign < 0 else '+'} offset = {h[0] / 1e6:.1f}-{h[1] / 1e6:.1f} MHz where "
                                 f"'{o['label']}' is {o['snr_db'] - s['snr_db']:.0f} dB stronger: mixer harmonic "
                                 f"conversion of that emitter")
                s["decoded"]["reason"] = "receiver artifact, not an emission at this frequency"
                changed.add(r["_path"])
            elif src:
                n, sign, o, h = src
                s["evidence"] += (f"; CAUTION possible {n}xLO-harmonic response of '{o['label']}' "
                                  f"({h[0] / 1e6:.1f}-{h[1] / 1e6:.1f} MHz, {o['snr_db'] - s['snr_db']:.0f} dB stronger); "
                                  f"seen at the same RF in {seen}/{len(covering)} other tuning(s)")
                s["confidence"] = round(min(s["confidence"], 0.45), 2)
                changed.add(r["_path"])
            elif covering and seen == 0 and s["snr_db"] >= 6:
                s["evidence"] += (f"; NOTE not seen at the same RF in {len(covering)} other tuning(s) covering it "
                                  f"(intermittent, or receiver-generated)")
                s["confidence"] = round(min(s["confidence"], 0.4), 2)
                changed.add(r["_path"])
    return changed


def occupancy_map(reps, path: Path, f_lo=300e6, f_hi=6000e6, df=250e3):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    nb = int((f_hi - f_lo) / df)
    snr = np.full(nb, np.nan)
    cat = np.full(nb, "", object)
    for r in reps:
        for s in r["signals"]:
            if {"receiver_spur", "lo_leakage", "iq_image", "lo_harmonic_response"} & set(s["flags"]):
                continue
            i0, i1 = int((s["f_low_hz"] - f_lo) / df), int(np.ceil((s["f_high_hz"] - f_lo) / df))
            i0, i1 = max(0, i0), min(nb, max(i1, i0 + 1))
            v = s["snr_db"] if s["snr_db"] is not None else 0
            for i in range(i0, i1):
                if np.isnan(snr[i]) or v > snr[i]:
                    snr[i], cat[i] = v, category(s["label"])
    col = {c: k for c, _, k in CATS}
    rows = 6
    span = (f_hi - f_lo) / rows
    fig, axs = plt.subplots(rows, 1, figsize=(22, 3.1 * rows))
    f = f_lo + (np.arange(nb) + 0.5) * df
    from .core.bands_india import BANDS
    for k, ax in enumerate(axs):
        a, b = f_lo + k * span, f_lo + (k + 1) * span
        m = (f >= a) & (f < b)
        for c in set(cat[m]) - {""}:
            mm = m & (cat == c)
            ax.bar(f[mm] / 1e6, np.nan_to_num(snr[mm]) + 3, width=df / 1e6, bottom=-3, color=col.get(c, "#fff"),
                   edgecolor="none")
        for lo_, hi_, key, desc, _ in BANDS:
            if hi_ > a and lo_ < b:
                ax.axvspan(max(lo_, a) / 1e6, min(hi_, b) / 1e6, color="#000", alpha=0.04)
                ax.text((max(lo_, a) + min(hi_, b)) / 2e6, 58, key, ha="center", fontsize=6, rotation=90, va="top",
                        color="#555")
        ax.set_xlim(a / 1e6, b / 1e6)
        ax.set_ylim(-3, 60)
        ax.set_ylabel("SNR dB")
        ax.grid(alpha=0.2)
    axs[-1].set_xlabel("frequency (MHz)")
    hs = [Patch(color=k, label=c) for c, _, k in CATS if c in set(cat)]
    axs[0].legend(handles=hs, ncol=8, fontsize=8, loc="upper left")
    axs[0].set_title("Occupancy 300 MHz - 6 GHz (strongest identified signal per 250 kHz; receiver artifacts removed)")
    fig.tight_layout()
    fig.savefig(path, dpi=90)
    plt.close(fig)


def _redraw(job):
    from .io.sigmf import load
    from .viz import annotate
    path, meta = job
    r = json.loads(Path(path).read_text())
    annotate.draw(r, load(meta).iq, Path(path).with_name(f"{r['file']}_annotated.png"))


MAX_DEVICES = 500
SEVERITY_ORDER = {"critical": 0, "high": 1, "medium": 2, "low": 3, "info": 4}


def _devices_across(reps):
    """Merge the per-snapshot device tracks and findings into one picture of the whole sweep.

    Snapshots are analysed in separate processes, so this is the only place an emitter can be followed
    from one to the next. A device heard in many snapshots is a fixture of the area; one heard in a
    single snapshot was passing through, was briefly active, or is at the edge of range -- which is
    itself worth knowing, so `n_snapshots` is reported rather than used to filter.
    """
    devs: dict = {}
    found: dict = {}
    for r in reps:
        when, where = r.get("capture_time_utc"), round(r["center_freq_hz"] / 1e6, 2)
        for d in r.get("devices", []) or []:
            k = d["mac"]
            e = devs.get(k)
            if e is None:
                e = devs[k] = {"mac": k, "vendor": d.get("vendor"), "address": d.get("address"),
                               "role": d.get("role"), "n_snapshots": 0, "n_frames": 0,
                               "snr_db_max": None, "channels": {}, "seen_at_fc_mhz": set(),
                               "first_seen_utc": when, "last_seen_utc": when}
            e["n_snapshots"] += 1
            e["n_frames"] += d.get("n_frames", 0)
            e["seen_at_fc_mhz"].add(where)
            if d.get("snr_db_max") is not None:
                e["snr_db_max"] = max(e["snr_db_max"], d["snr_db_max"]) if e["snr_db_max"] is not None \
                    else d["snr_db_max"]
            for ch, n in (d.get("channels") or {}).items():
                e["channels"][str(ch)] = e["channels"].get(str(ch), 0) + n
            if when:
                e["first_seen_utc"] = min(e["first_seen_utc"] or when, when)
                e["last_seen_utc"] = max(e["last_seen_utc"] or when, when)
            if d.get("role") == "access point":
                e["role"] = "access point"
            for key in ("ssid", "security", "device_name", "drone", "remote_id", "phy",
                        "associated_bssid", "probed_ssids"):
                if d.get(key) is not None and e.get(key) is None:
                    e[key] = d[key]
        for f in r.get("findings", []) or []:
            k = (f["kind"], f.get("ssid") or "", tuple(sorted(f.get("devices", []))))
            e = found.get(k)
            if e is None:
                e = found[k] = {kk: vv for kk, vv in f.items() if kk != "signal_ids"}
                e["n_snapshots"] = 0
                e["seen_at_fc_mhz"] = set()
            e["n_snapshots"] += 1
            e["seen_at_fc_mhz"].add(where)
    out = []
    for e in devs.values():
        e["seen_at_fc_mhz"] = sorted(e["seen_at_fc_mhz"])
        out.append({k: v for k, v in e.items() if v not in (None, {}, [])})
    out.sort(key=lambda d: (-d["n_snapshots"], -(d.get("snr_db_max") or -99)))
    fl = []
    for e in found.values():
        e["seen_at_fc_mhz"] = sorted(e["seen_at_fc_mhz"])
        fl.append(e)
    fl.sort(key=lambda f: (SEVERITY_ORDER.get(f["severity"], 9), -f["n_snapshots"]))
    aps = [d for d in out if d.get("role") == "access point"]
    summary = {
        "n_devices": len(out),
        "n_access_points": len(aps),
        "n_stations": len(out) - len(aps),
        "n_randomized": sum(1 for d in out if (d.get("address") or "").startswith("randomized")),
        "n_in_one_snapshot_only": sum(1 for d in out if d["n_snapshots"] == 1),
        "n_drones": sum(1 for d in out if d.get("drone")),
        "severity_counts": {k: sum(1 for f in fl if f["severity"] == k)
                            for k in SEVERITY_ORDER if any(f["severity"] == k for f in fl)},
        "note": "devices and findings merged across snapshots by MAC address (receivers.devices)",
    }
    return out[:MAX_DEVICES], summary, fl


def build(sweep_out: Path, png: bool = True):
    sweep_out = Path(sweep_out)
    reps = _load(sweep_out)
    if not reps:
        print("no reports")
        return
    from .report.writer import split_weak
    n_before = sum(len(r["signals"]) for r in reps)
    for r in reps:
        k = len(r["signals"])
        split_weak(r)
        r["_weak_changed"] = len(r["signals"]) != k
    print(f"weak-detection split: {n_before} -> {sum(len(r['signals']) for r in reps)} signals", flush=True)
    changed = spur_check(reps) | harmonic_check(reps) | {r["_path"] for r in reps if r.pop("_weak_changed")}
    for r in reps:
        if r["_path"] in changed:
            p = Path(r.pop("_path"))
            p.write_text(json.dumps(r, indent=1))
            r["_path"] = str(p)
    if png:
        run_log = {Path(x["report"]).name: x["meta"] for x in
                   json.loads((sweep_out / "run_log.json").read_text()) if "report" in x} \
            if (sweep_out / "run_log.json").exists() else {}
        def stale_png(rp):
            png = Path(rp).with_name(Path(rp).name.replace("_report.json", "_annotated.png"))
            return not png.exists() or png.stat().st_mtime < Path(rp).stat().st_mtime
        jobs = [(r["_path"], run_log[Path(r["_path"]).name]) for r in reps
                if (r["_path"] in changed or stale_png(r["_path"])) and Path(r["_path"]).name in run_log]
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(2) as ex:
            list(ex.map(_redraw, jobs, chunksize=4))
        print(f"redrew {len(jobs)} annotated PNG(s)", flush=True)
    # CSV
    with open(sweep_out / "signals.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file", "fc_mhz", "id", "label", "category", "confidence",
                    "f_low_mhz", "f_high_mhz", "bw_khz", "t_start_ms", "t_end_ms", "snr_db", "modulation",
                    "decode_status", "flags"])
        for r in reps:
            for s in r["signals"]:
                w.writerow([r["file"], r["center_freq_hz"] / 1e6, s["id"], s["label"], category(s["label"]),
                            s["confidence"], round(s["f_low_hz"] / 1e6, 4), round(s["f_high_hz"] / 1e6, 4),
                            round(s["bandwidth_hz"] / 1e3, 1), round(s["t_start_s"] * 1e3, 3),
                            round(s["t_end_s"] * 1e3, 3), s["snr_db"], s["modulation"], s["decoded"]["status"],
                            ";".join(s["flags"])])
    # identities
    ident = defaultdict(set)
    cats = Counter()
    status = Counter()
    for r in reps:
        for s in r["signals"]:
            cats[category(s["label"])] += 1
            status[s["decoded"]["status"]] += 1
            f = s["decoded"]["fields"]
            for k in ("ssids", "bssids"):
                for v in f.get(k, []) or []:
                    ident[f"wifi_{k}"].add(v)
            for k in ("pci", "bsic", "icao", "earfcn", "arfcn", "nr_arfcn"):
                if f.get(k) is not None:
                    ident[k].add(f"{f[k]}@{s['f_center_hz'] / 1e6:.2f}MHz")
            for b in f.get("bursts", []) or []:
                d = b.get("decoded") or {}
                if d.get("adv_addr"):
                    ident["ble_adv_addr"].add(d["adv_addr"])
    devices, dev_summary, findings = _devices_across(reps)
    lab = Counter((r.get("frequency_label_check") or {}).get("status", "not_checked") for r in reps)
    clip = sorted({round(r["center_freq_hz"] / 1e6, 2) for r in reps if r["capture_quality"].get("clipping")})
    summ = {
        "sweep_output": str(sweep_out),
        "n_snapshots": len(reps),
        "n_signals": sum(len(r["signals"]) for r in reps),
        "n_weak_detections": sum(len(r.get("weak_detections", [])) for r in reps),
        "signals_by_category": dict(cats.most_common()),
        "decode_status": dict(status),
        "decoded_identities": {k: sorted(v) for k, v in ident.items()},
        "devices": devices,
        "device_summary": dev_summary,
        "findings": findings,
        "data_quality": {
            "frequency_label_check": dict(lab),
            "adc_clipping_at_fc_mhz": clip,
            "receiver_spurs_flagged": sum("receiver_spur" in s["flags"] for r in reps for s in r["signals"]),
        },
    }
    (sweep_out / "summary.json").write_text(json.dumps(summ, indent=1))
    if png:
        occupancy_map(reps, sweep_out / "occupancy_map.png")
    print(json.dumps({k: v for k, v in summ.items() if k != "decoded_identities"}, indent=1))
    return summ
