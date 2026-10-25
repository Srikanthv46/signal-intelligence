#!/usr/bin/env python3
"""
drone_survey.py -- which drones are operating in a capture, and which bursts
belong to which drone.

  capture -> spurs + shared channelizer
          -> detectors (each `detect(x, fs, fc, ctx) -> [Emission]`):
               det_dji          OcuSync video / RC uplink / DroneID (+decode)
               det_digital_vtx  Herelink (LTE PSS/SSS, PCI) + blind unknown-OFDM VTX
               det_analog_vtx   FM analog video (PAL/NTSC line-rate)
               det_rc_link      ELRS/LoRa, GFSK FHSS control links
          -> drone_assoc: physical-layer fingerprints (frame-clock lattice, oscillator
             ppm, simultaneity veto, IDs) group bursts into emitters (one radio,
             possibly hopping); emitters sharing an oscillator are linked into one
             drone system (drone + its RC)
          -> report / JSON / SigMF annotations / PNG

A detector that fails or is missing is reported, never fatal: the survey still
covers everything else.

Usage
  python drone_survey.py cap.sigmf-meta [--out prefix] [--only dji,digital_vtx,...]
"""
from __future__ import annotations
import importlib, json, os, time, traceback
from collections import Counter, defaultdict
import numpy as np

from spectrum_intel.core import common as C
from . import spectrum_survey as SS
from spectrum_intel.receivers.emitters.families import DRONE_FAMILIES, LABEL_ALIASES, _assoc_models  # noqa: F401
from spectrum_intel.receivers.rc_link.label_check import BLE_ADV_MHZ, check_frequency_label  # noqa: F401
DETECTORS = [("dji", "spectrum_intel.receivers.dji.detector"),
             ("digital_vtx", "spectrum_intel.receivers.digital_vtx.detector"),
             ("analog_vtx", "spectrum_intel.receivers.analog_vtx.detector"),
             ("rc_link", "spectrum_intel.receivers.rc_link.detector")]


def survey(meta_path, only=None, log=print, fc_override=None):
    x, fs, fc, meta = C.load_capture(meta_path)
    if fc_override:
        log(f"fc override: metadata {fc/1e6:g} MHz -> {fc_override/1e6:g} MHz")
        fc = fc_override
    log(f"{len(x)} samples @ {fs/1e6:g} MSps, fc {fc/1e6:g} MHz, {len(x)/fs*1e3:.2f} ms")
    spurs = SS.find_spurs(x, fs, fc)
    ch = C.Channelizer(x.astype(np.complex128), fs, [(s["offset_hz"], 30e3) for s in spurs])
    ctx = dict(channelizer=ch, spurs=spurs, meta=meta)
    emissions, runs = [], {}
    for name, mod in DETECTORS:
        if only and name not in only:
            continue
        t0 = time.time()
        try:
            det = importlib.import_module(mod)
            em = det.detect(x, fs, fc, ctx)
            runs[name] = dict(ok=True, n=len(em), seconds=round(time.time() - t0, 1),
                              labels=dict(Counter(f"{e.family}/{e.label}" for e in em)))
            emissions += em
        except Exception as e:                           # one broken detector must not hide the rest
            runs[name] = dict(ok=False, error=f"{type(e).__name__}: {e}", seconds=round(time.time() - t0, 1),
                              trace=traceback.format_exc(limit=3))
        log(f"  {name:12s} {runs[name]}")
    emissions.sort(key=lambda e: (e.start_us, e.center_mhz))
    explain_by_continuous(emissions)
    drone_em = [e for e in emissions if e.family in DRONE_FAMILIES and "explained_by" not in e.features]
    other_em = [e for e in emissions if e not in drone_em]
    label_check = check_frequency_label(other_em)
    if label_check:
        log(f"  frequency label check: {label_check['note']} ({label_check['n_adverts']} BLE adverts)")

    # ---- association: which bursts are one radio, which radios are one drone system
    emitters, assoc_err = [], None
    try:
        from spectrum_intel.receivers.emitters import association as A
        A.measure_fingerprints(x, fs, fc, drone_em, ch=ch)
        emitters = A.associate(drone_em, fc, models=_assoc_models())
    except Exception as e:
        assoc_err = f"{type(e).__name__}: {e}"
        log(f"  association failed: {assoc_err}")
    systems = drone_systems(emitters) if not assoc_err else unassociated_systems(drone_em, assoc_err)
    return dict(x=x, fs=fs, fc=fc, meta=meta, spurs=spurs, runs=runs, emissions=drone_em,
                other=other_em, emitters=emitters, systems=systems, assoc_error=assoc_err,
                label_check=label_check, fc_override=bool(fc_override))


# Families whose detections rest on weak/generic evidence and can be produced by the
# sidebands of a strong continuous emitter. DJI / Herelink need hard signatures (ZC/PSS).
WEAK_FAMILIES = {"digital_vtx", "rc_fhss", "telemetry", "narrowband", "unknown", "lora", "bluetooth"}


def explain_by_continuous(emissions):
    """An analog VTX transmits for the whole capture and its FM spectrum (sidebands,
    audio subcarrier, field-blanking transients) reaches beyond its 99.5% band: the line
    period alone looks like a 16 kHz 'OFDM' symbol to a blind numerology estimator.
    Weak-evidence detections inside its band (+ margin) are attributed to it, not
    reported as separate drones. Trade-off: a real RC burst under the VTX is hidden."""
    cont = [e for e in emissions if e.family == "analog_vtx"]
    for e in emissions:
        if e.family not in WEAK_FAMILIES:
            continue
        for a in cont:
            margin = max(6.0, 0.2 * (a.f_hi_mhz - a.f_lo_mhz))
            if a.f_lo_mhz - margin <= e.center_mhz <= a.f_hi_mhz + margin:
                e.features["explained_by"] = f"analog_vtx/{a.label} @ {a.center_mhz:.2f} MHz"
                break


def unassociated_systems(drone_em, err):
    """Association failed: still report every detection, grouped by type, unassociated."""
    groups = defaultdict(list)
    for e in drone_em:
        groups[(e.family, e.label)].append(e)
    return [dict(system_id=f"D{i + 1}", kind=DRONE_FAMILIES.get(f, f), emitters=[], labels=[f"{f}/{l}"],
                 frequencies_mhz=sorted({round(e.center_mhz, 3) for e in es}), n_bursts=len(es), ids={},
                 min_confidence=None, note=f"UNASSOCIATED (association failed: {err})")
            for i, ((f, l), es) in enumerate(sorted(groups.items()))]


def drone_systems(emitters):
    """Union emitters that drone_assoc linked (same oscillator, same vendor family):
    e.g. a DJI uplink + video + DroneID = one drone system."""
    parent = {e.emitter_id: e.emitter_id for e in emitters}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for e in emitters:
        for l in e.linked:
            if l["other"] in parent:
                parent[find(e.emitter_id)] = find(l["other"])
    groups = defaultdict(list)
    for e in emitters:
        groups[find(e.emitter_id)].append(e)
    out = []
    for i, es in enumerate(sorted(groups.values(), key=lambda g: min(e.emitter_id for e in g))):
        fams = sorted({e.family for e in es})
        ids = {}
        for e in es:
            for k in ("serial", "pci", "mac"):
                v = e.fingerprint.get(k) if isinstance(e.fingerprint, dict) else None
                if v is not None:
                    ids[k] = v
        out.append(dict(system_id=f"D{i + 1}", kind=", ".join(DRONE_FAMILIES.get(f, f) for f in fams),
                        emitters=[e.emitter_id for e in es], labels=[f"{e.family}/{e.label}" for e in es],
                        frequencies_mhz=sorted({f for e in es for f in e.frequencies}),
                        n_bursts=sum(len(e.members) for e in es), ids=ids,
                        min_confidence=min((e.confidence for e in es if e.confidence is not None), default=None)))
    return out


# ================================================================ outputs
def _clean(o):
    if isinstance(o, dict): return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [_clean(v) for v in o]
    if isinstance(o, np.integer): return int(o)
    if isinstance(o, np.floating): return float(o)
    if isinstance(o, np.bool_): return bool(o)
    if isinstance(o, np.ndarray): return _clean(o.tolist())
    return o


def write_outputs(r, out, make_plot=True):
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    m = r["meta"]
    inv = dict(capture=dict(fc_mhz=r["fc"] / 1e6, fs_msps=r["fs"] / 1e6, duration_ms=len(r["x"]) / r["fs"] * 1e3,
                            datetime=m["captures"][0].get("core:datetime")),
               detectors=r["runs"], systems=r["systems"], emitters=[e.to_dict() for e in r["emitters"]],
               emissions=[e.to_dict() for e in r["emissions"]],
               non_drone=dict(Counter(f"{e.family}/{e.label}" for e in r["other"])),
               spurs=r["spurs"], association_error=r["assoc_error"],
               frequency_label_check=r["label_check"], fc_overridden=r["fc_override"])
    json.dump(_clean(inv), open(out + "_drones.json", "w"), indent=1)
    mm = json.loads(json.dumps(m))
    eid = {i: e.emitter_id for e in r["emitters"] for i in e.members}
    mm["annotations"] = sorted([{
        "core:sample_start": int(e.start_us * 1e-6 * r["fs"]),
        "core:sample_count": max(1, int((e.end_us - e.start_us) * 1e-6 * r["fs"])),
        "core:freq_lower_edge": e.f_lo_mhz * 1e6, "core:freq_upper_edge": e.f_hi_mhz * 1e6,
        "core:label": f"{e.family}/{e.label}", "core:comment": f"emitter {eid.get(i, '?')}",
    } for i, e in enumerate(r["emissions"])], key=lambda a: a["core:sample_start"])
    json.dump(mm, open(out + "_drones.sigmf-meta", "w"), indent=2)
    if make_plot:
        plot(r, out + "_drones.png")


def plot(r, path, nfft=1024):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    x, fs, fc = r["x"], r["fs"], r["fc"]
    n = len(x) // nfft * nfft
    S = 10 * np.log10(np.abs(np.fft.fftshift(np.fft.fft(x[:n].reshape(-1, nfft) * np.hanning(nfft), axis=1),
                                             axes=1)) ** 2 + 1e-12).T
    fig, ax = plt.subplots(figsize=(18, 8), dpi=110)
    ax.imshow(S, aspect="auto", origin="lower", cmap="gray",
              extent=[0, n / fs * 1e3, (fc - fs / 2) / 1e6, (fc + fs / 2) / 1e6],
              vmin=np.percentile(S, 5), vmax=np.percentile(S, 99.7))
    cmap = plt.get_cmap("tab10")
    eid = {i: k for k, e in enumerate(r["emitters"]) for i in e.members}
    for i, e in enumerate(r["emissions"]):
        k = eid.get(i)
        c = cmap(k % 10) if k is not None else "w"
        ax.add_patch(Rectangle((e.start_us / 1e3, e.f_lo_mhz), (e.end_us - e.start_us) / 1e3,
                               e.f_hi_mhz - e.f_lo_mhz, fill=False, lw=1.3, ec=c))
    for k, e in enumerate(r["emitters"]):
        ax.plot([], [], color=cmap(k % 10), label=f"{e.emitter_id} {e.family}/{e.label} ({len(e.members)})")
    if r["emitters"]:
        ax.legend(loc="upper right", fontsize=7, framealpha=0.7)
    ax.set_xlabel("time (ms)"); ax.set_ylabel("frequency (MHz)")
    ax.set_title(f"Drone survey: {len(r['systems'])} drone system(s), {len(r['emitters'])} emitter(s), "
                 f"{len(r['emissions'])} bursts  (colour = emitter)")
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


def print_report(r, name=""):
    print(f"\n=== {name}  fc {r['fc']/1e6:g} MHz" + ("  (overridden)" if r["fc_override"] else ""))
    if r["label_check"]:
        print(f"  frequency label check: {r['label_check']['note']}")
    for k, v in r["runs"].items():
        print(f"  detector {k:12s} " + (f"{v['n']} bursts {v['labels']} ({v['seconds']} s)" if v["ok"]
                                        else f"FAILED {v['error']}"))
    if not r["systems"]:
        print("  no drone emissions detected")
    for s in r["systems"]:
        print(f"  {s['system_id']}: {s['kind']} -- {s['n_bursts']} bursts, emitters {s['emitters']}, "
              f"{len(s['frequencies_mhz'])} frequencies {s['frequencies_mhz'][:8]}"
              + (f", ids {s['ids']}" if s["ids"] else "")
              + (f", min confidence {s['min_confidence']}" if s["min_confidence"] is not None else "")
              + (f"  [{s['note']}]" if s.get("note") else ""))
    if r["emitters"]:
        from spectrum_intel.receivers.emitters import association as A
        print("  " + A.report(r["emitters"]).replace("\n", "\n  "))
    expl = [e for e in r["other"] if "explained_by" in e.features]
    if expl:
        cnt = Counter(e.family + "/" + e.label + " <- " + e.features["explained_by"] for e in expl)
        print(f"  attributed to a continuous emitter (not separate drones): {dict(cnt)}")
    rest = [e for e in r["other"] if "explained_by" not in e.features]
    if rest:
        print(f"  non-drone: {dict(Counter(f'{e.family}/{e.label}' for e in rest))}")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Detect drones and group their bursts into emitters")
    ap.add_argument("meta", nargs="+")
    ap.add_argument("--out", default="drone_survey")
    ap.add_argument("--only", help="comma list of detectors: " + ",".join(n for n, _ in DETECTORS))
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--fc", type=float, help="override the capture centre frequency (MHz) when the metadata is wrong")
    a = ap.parse_args()
    only = set(a.only.split(",")) if a.only else None
    for mp in a.meta:
        r = survey(mp, only, fc_override=a.fc * 1e6 if a.fc else None)
        prefix = a.out if len(a.meta) == 1 else os.path.join(a.out, os.path.basename(mp).replace(".sigmf-meta", ""))
        write_outputs(r, prefix, not a.no_plot)
        print_report(r, os.path.basename(mp))
